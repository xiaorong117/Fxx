#include "bubble/model.hpp"
#include "bubble/pardiso_backend.hpp"
#include "bubble/analysis.hpp"
#include "bubble/amgx_backend.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace {

using json = nlohmann::json;

std::string step_name(int step, const char* extension) {
  std::ostringstream out;
  out << "state_" << std::setw(6) << std::setfill('0') << step << extension;
  return out.str();
}

const char* state_extension(const bubble::Config& config){
  return config.io.state_format=="compressed_vtp"?".vtp":".vtk";
}
void write_state(const std::filesystem::path& path,const bubble::Solver& solver,
                 const bubble::Assembly* assembly){
  if(solver.config().io.state_format=="compressed_vtp")bubble::write_vtp(path,solver,assembly);
  else {bubble::write_vtk(path,solver,assembly);bubble::audit_vtk(path,solver.network().pores.size(),solver.network().edges.size());}
}

double free_total(const bubble::Solver& solver) {
  double total = 0.0;
  for (std::size_t i = 0; i < solver.state().size(); ++i) {
    const auto& s = solver.state()[i];
    if (s.active_gas) {
      total += static_cast<double>(bubble::free_gas_moles(
          s.pg, s.sw, solver.network().pores[i].volume, solver.config().physics));
    }
  }
  return total;
}

void audit_state(const bubble::Solver& solver, double previous_free) {
  const auto& config = solver.config();
  for (const auto& s : solver.state()) {
    if (!std::isfinite(s.pl) || config.physics.background_abs_Pa + s.pl <= 0.0 || !std::isfinite(s.sw) ||
        !std::isfinite(s.c))
      throw std::runtime_error("state contains invalid absolute pressure or NaN/Inf");
    if (s.c < -1e-12 * config.scaling.concentration || s.sw <= 0.0 || s.sw > 1.0)
      throw std::runtime_error("state violates C/Sw bounds");
    if (s.active_gas && (!std::isfinite(s.pg) || config.physics.background_abs_Pa + s.pg <= 0.0 ||
                         s.sw <= config.active_set.sw_min ||
                         s.sw >= 1.0 - config.active_set.sg_off))
      throw std::runtime_error("active gas state violates Pg/Sw bounds");
  }
  const double now = free_total(solver);
  const double tolerance = 1e-12 * std::max(previous_free, config.scaling.mole_total);
  if (now > previous_free + tolerance)
    throw std::runtime_error("free gas increased in dissolution-only mode");
}

std::string command_line(int argc, char** argv) {
  std::ostringstream out;
  for (int i = 0; i < argc; ++i) {
    if (i) out << ' ';
    out << argv[i];
  }
  return out.str();
}

std::vector<std::string> command_arguments(int argc, char** argv) {
  std::vector<std::string> result;
  result.reserve(static_cast<std::size_t>(argc));
  for (int i = 0; i < argc; ++i) result.emplace_back(argv[i]);
  return result;
}

struct PvdEntry { double time = 0.0; std::string file; };

struct DtDecision { double next = 0.0, gas_cap = -1.0, concentration_cap = -1.0; std::size_t concentration_node=0; double concentration_C=0.0,concentration_dCdt=0.0,concentration_reference=0.0,concentration_tau=0.0; std::string reason = "initial"; };
DtDecision choose_next_dt(const bubble::Solver& solver, const bubble::StepReport& report,
                          double accepted_dt, int consecutive_growth) {
  DtDecision d; const auto& c=solver.config(); d.next=accepted_dt;
  if(!c.time_step_control.enabled){
    if(report.newton_iterations<=c.time.fast_newton_iterations && consecutive_growth>=1){
      d.next=std::min(c.time.max_dt,accepted_dt*c.time.growth_factor);d.reason="easy_newton_growth";
    } else d.reason="normal_newton_keep";
    return d;
  }
  if(report.newton_iterations<=c.time_step_control.easy_newton_iterations &&
     consecutive_growth<c.time_step_control.max_consecutive_growth){d.next=accepted_dt*c.time_step_control.growth_factor;d.reason="easy_newton_growth";}
  else if(report.newton_iterations<=c.time_step_control.normal_newton_iterations){d.reason="normal_newton_keep";}
  else {d.next=accepted_dt*0.8;d.reason="hard_newton_shrink";}
  std::vector<double> gas_taus;
  double min_tau_c=std::numeric_limits<double>::infinity();
  const double c_reference=c.physics.henry_cp*c.physics.background_abs_Pa;
  double total_free=0.0;
  for(std::size_t i=0;i<solver.state().size();++i) if(solver.state()[i].active_gas)
    total_free += static_cast<double>(bubble::free_gas_moles(solver.state()[i].pg,solver.state()[i].sw,solver.network().pores[i].volume,c.physics));
  for(std::size_t i=0;i<solver.state().size();++i){const auto&s=solver.state()[i];if(!s.active_gas||i>=report.final_assembly.transfer.size())continue;
    const double rate=std::abs(report.final_assembly.transfer[i]);if(!(rate>0.0)||!std::isfinite(rate))continue;
    const double ng=static_cast<double>(bubble::free_gas_moles(s.pg,s.sw,solver.network().pores[i].volume,c.physics));
    const double sg=static_cast<double>(1.0L-s.sw);
    // Bubbles already inside the disappearance-event neighborhood must not
    // throttle the global timestep.  They are advanced by the dedicated
    // event locator and retired by the active set on an accepted event step.
    const double mole_fraction = total_free > 0.0 ? ng / total_free : 0.0;
    if(ng>1000.0*c.active_set.ng_off&&sg>1000.0*c.active_set.sg_off &&
       mole_fraction >= c.time_step_control.gas_timescale_min_mole_fraction)
      gas_taus.push_back(ng/rate);
    const double liquid=solver.network().pores[i].volume*static_cast<double>(s.sw);
    if(c.time_step_control.concentration_change_limit_enabled&&liquid>0.0){
      const double dcdt=rate/liquid;const double tau=c_reference/dcdt;
      if(tau<min_tau_c){min_tau_c=tau;d.concentration_cap=c.time_step_control.concentration_change_fraction*tau;
        d.concentration_node=solver.network().pores[i].id;d.concentration_C=static_cast<double>(s.c);
        d.concentration_dCdt=dcdt;d.concentration_reference=c_reference;d.concentration_tau=tau;}
    }
  }
  if(!gas_taus.empty()) {
    std::sort(gas_taus.begin(), gas_taus.end());
    const double q=std::clamp(c.time_step_control.gas_timescale_percentile,0.0,1.0);
    const std::size_t qi=std::min(gas_taus.size()-1, static_cast<std::size_t>(std::floor(q*static_cast<double>(gas_taus.size()-1))));
    const double tau_g=gas_taus[qi];
    d.gas_cap=c.time_step_control.gas_timescale_fraction*tau_g;
    if(d.gas_cap<d.next){d.next=d.gas_cap;d.reason="gas_timescale_cap";}
  }
  if(std::isfinite(min_tau_c)&&d.concentration_cap<d.next){d.next=d.concentration_cap;d.reason="concentration_timescale_cap";}
  d.next=std::clamp(d.next,c.time.min_dt,c.time.max_dt);return d;
}

std::vector<PvdEntry> read_pvd(const std::filesystem::path& path) {
  std::vector<PvdEntry> result;
  std::ifstream input(path);
  if (!input) return result;
  const std::regex pattern("timestep=\\\"([^\\\"]+)\\\"[^>]*file=\\\"([^\\\"]+)\\\"");
  std::string line;
  while (std::getline(input, line)) {
    std::smatch match;
    if (std::regex_search(line, match, pattern))
      result.push_back({std::stod(match[1].str()), match[2].str()});
  }
  return result;
}

void write_pvd(const std::filesystem::path& path, const std::vector<PvdEntry>& entries,
               double mass_transfer_multiplier) {
  std::ofstream output(path);
  if (!output) throw std::runtime_error("cannot write PVD collection");
  output << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <Collection>\n" << std::setprecision(17);
  for (const auto& entry : entries)
    output << "    <DataSet timestep=\"" << entry.time
           << "\" group=\"mass_transfer_multiplier=" << mass_transfer_multiplier
           << "\" part=\"0\" file=\"" << entry.file << "\"/>\n";
  output << "  </Collection>\n</VTKFile>\n";
}

std::uintmax_t directory_bytes(const std::filesystem::path& path) {
  std::uintmax_t total=0;
  for(const auto& item:std::filesystem::recursive_directory_iterator(path))
    if(item.is_regular_file()) total+=item.file_size();
  return total;
}

std::string gpu_inventory() {
  std::string result;std::array<char,512> buffer{};
  FILE* pipe=popen("nvidia-smi --query-gpu=index,name,memory.total,driver_version,compute_cap --format=csv,noheader 2>&1","r");
  if(!pipe)return "nvidia-smi unavailable";
  while(fgets(buffer.data(),static_cast<int>(buffer.size()),pipe))result+=buffer.data();
  const int status=pclose(pipe);if(status!=0)result="nvidia-smi failed: "+result;return result;
}

std::string checkpoint_time_name(double time) {
  std::ostringstream output;
  output << "checkpoint_t" << std::fixed << std::setprecision(6) << time << ".json";
  return output.str();
}

}  // namespace

int main(int argc, char** argv) {
  const auto wall_started = std::chrono::steady_clock::now();
  const std::string started = bubble::utc_timestamp();
  const std::string invoked_as = command_line(argc, argv);
  const auto invoked_argv = command_arguments(argc, argv);
  std::optional<bubble::Config> loaded_config;
  int accepted_steps_before_failure = 0;
  double time_before_failure = 0.0;
  std::ofstream run_log;
  try {
    std::filesystem::path config_path;
    std::filesystem::path restart_path;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--config" && i + 1 < argc) config_path = argv[++i];
      else if (arg == "--restart" && i + 1 < argc) restart_path = argv[++i];
      else throw std::runtime_error("usage: bubble_solver --config case.json [--restart checkpoint.json]");
    }
    if (config_path.empty()) throw std::runtime_error("--config is required");
    auto config = bubble::Config::load(config_path);
    loaded_config = config;
    std::filesystem::create_directories(config.io.output_directory);
    const auto run_log_path = config.io.output_directory / "run.log";
    run_log.open(run_log_path, restart_path.empty() ? std::ios::trunc : std::ios::app);
    if (!run_log) throw std::runtime_error("cannot open run.log");
    run_log << std::scientific << std::setprecision(17);
    run_log << "started_utc=" << started << '\n'
            << "config=" << std::filesystem::absolute(config_path).string() << '\n';
    const auto network = bubble::Network::load(config);
    const auto checkpoint_provenance = bubble::make_checkpoint_provenance(config);
    bubble::Solver solver(network, config);
    solver.initialize();
    std::optional<bubble::CheckpointData> restart;
    std::string restart_input_hash;
    if (!restart_path.empty()) {
      restart_input_hash = bubble::sha256_file(restart_path);
      restart = bubble::read_checkpoint(restart_path, network.pores.size());
      if (restart->background_abs_Pa != config.physics.background_abs_Pa)
        throw std::runtime_error("checkpoint background pressure mismatch");
      if (restart->mass_transfer_multiplier != config.physics.mass_transfer_multiplier)
        throw std::runtime_error("checkpoint mass-transfer multiplier mismatch");
      if (restart->time > config.time.end || restart->next_dt < config.time.min_dt ||
          restart->next_dt > config.time.max_dt) {
        throw std::runtime_error("checkpoint continuation time/dt is incompatible with config");
      }
      if (restart->continuation_config_sha256 !=
              checkpoint_provenance.continuation_config_sha256 ||
          restart->pore_sha256 != checkpoint_provenance.pore_sha256 ||
          restart->throat_sha256 != checkpoint_provenance.throat_sha256 ||
          restart->connect_sha256 != checkpoint_provenance.connect_sha256) {
        throw std::runtime_error("checkpoint provenance does not match config/network inputs");
      }
      solver.restore(restart->state, restart->time, restart->step_number,
                     restart->cumulative_mole_ledger,
                     restart->cumulative_volume_ledger,
                     restart->cumulative_net_liquid_inflow,
                     restart->cumulative_boundary_component_outflow,
                     restart->cumulative_outlet_advective_component_outflow,
                     restart->cumulative_liquid_inflow,
                     restart->cumulative_liquid_outflow,
                     restart->cumulative_transfer_by_node);
      solver.restore_gross_fluxes(restart->cumulative_component_inflow, restart->cumulative_component_outflow);
      solver.restore_inlet_pressure(restart->inlet_pressure);
      solver.set_adaptive_metadata(restart->last_dt_change_reason,
          restart->last_min_free_gas_moles,restart->last_min_sg);
      solver.restore_gas_disappearance_metadata(
          restart->gas_disappearance_events,
          restart->last_gas_disappearance_node_id,
          restart->last_gas_disappearance_time);
    }
    bubble::ScientificRecorder scientific(solver, !restart_path.empty());
    const auto initial_scientific = scientific.snapshot(solver);
    if (config.termination.enabled &&
        config.termination.primary != "residual_free_gas_fraction" &&
        initial_scientific.main_real_gas_saturation <=
            config.termination.target_gas_saturation)
      throw std::runtime_error(
          "initial main-real gas saturation is already at/below termination target");
    if (config.termination.enabled &&
        (config.termination.primary == "residual_free_gas_fraction" ||
         config.termination.primary == "either") &&
        config.termination.residual_free_gas_fraction > 0.0 &&
        initial_scientific.free_gas_mol <= scientific.initial_free_moles() *
            config.termination.residual_free_gas_fraction)
      throw std::runtime_error(
          "initial free gas fraction is already at/below termination target");
    accepted_steps_before_failure = solver.step_number();
    time_before_failure = solver.time();
    const auto degenerate_audit_path =
        config.io.output_directory / "degenerate_geometry_audit.json";
    bubble::write_degenerate_geometry_audit(degenerate_audit_path, solver);
    const auto component_audit_path =
        config.io.output_directory / "component_boundary_and_makeup_audit.json";
    if (restart_path.empty() || !std::filesystem::exists(component_audit_path))
      bubble::write_component_boundary_and_makeup_audit(component_audit_path, solver);
    const auto physical_parameter_audit_path =
        config.io.output_directory / "physical_parameter_audit.json";
    if (config.case_info.kind == "reference_physics_with_unvalidated_closures" &&
        (restart_path.empty() || !std::filesystem::exists(physical_parameter_audit_path)))
      bubble::write_physical_parameter_audit(physical_parameter_audit_path, solver);
    std::vector<std::filesystem::path> outputs;
    outputs.push_back(degenerate_audit_path);
    outputs.push_back(component_audit_path);
    if (std::filesystem::is_regular_file(physical_parameter_audit_path))
      outputs.push_back(physical_parameter_audit_path);
    outputs.push_back(run_log_path);
    std::filesystem::path input_geometry_audit_path;
    std::filesystem::path virtual_boundary_audit_path;
    if (!config.io.input_geometry_audit_file.empty()) {
      input_geometry_audit_path =
          config.io.output_directory / "input_geometry_unit_audit.json";
      std::filesystem::copy_file(config.io.input_geometry_audit_file,
                                 input_geometry_audit_path,
                                 std::filesystem::copy_options::overwrite_existing);
      outputs.push_back(input_geometry_audit_path);
    }
    if (!config.io.virtual_boundary_geometry_audit_file.empty()) {
      virtual_boundary_audit_path =
          config.io.output_directory / "virtual_boundary_geometry_audit.json";
      std::filesystem::copy_file(config.io.virtual_boundary_geometry_audit_file,
                                 virtual_boundary_audit_path,
                                 std::filesystem::copy_options::overwrite_existing);
      outputs.push_back(virtual_boundary_audit_path);
    }
    run_log << "nodes=" << network.pores.size() << " edges=" << network.edges.size()
            << " virtual_nodes=" << network.virtual_boundary_nodes
            << " boundary_edges=" << network.boundary_edges << '\n';
    const auto pvd_path = config.io.output_directory / "states.pvd";
    std::vector<PvdEntry> pvd_entries = restart_path.empty()
        ? std::vector<PvdEntry>{} : read_pvd(pvd_path);

    const auto diagnostics_path = config.io.output_directory / "diagnostics.tsv";
    const bool write_diagnostics_header = restart_path.empty() ||
        !std::filesystem::exists(diagnostics_path) ||
        std::filesystem::file_size(diagnostics_path) == 0;
    std::ofstream diagnostics(diagnostics_path, restart_path.empty() ? std::ios::trunc : std::ios::app);
    if (!diagnostics) throw std::runtime_error("cannot open diagnostics TSV");
    diagnostics << std::scientific << std::setprecision(17);
    if (write_diagnostics_header) {
      diagnostics << "step\ttime_s\trequested_dt_s\tdt_s\ttimestep_retries"
                     "\tnewton_iterations\tlinear_iterations\tlinear_residual"
                     "\tresidual_scaled_inf\tstate_scaled_inf\tupdate_scaled_inf"
                     "\trelative_update_norm\tupdate_tolerance_limit\tupdate_criterion_ratio"
                     "\tresidual_criterion_passed\tupdate_criterion_passed"
                     "\tRl_scaled_inf\tRd_scaled_inf\tRg_scaled_inf\tRc_scaled_inf"
                     "\tRl_max_node_id\tRd_max_node_id\tRg_max_node_id\tRc_max_node_id"
                     "\tboundary_liquid_outflow_m3_s\tdissolved_mol\tfree_mol"
                     "\tboundary_component_outflow_mol_s\tepsilon_N_equation\tepsilon_V_equation"
                     "\tepsilon_N_adjusted\tepsilon_V_adjusted\tactive_gas_count"
                     "\tthreshold_moles_transferred\tactive_set_mole_ledger"
                     "\tactive_set_liquid_volume_ledger"
                     "\tfree_mol_total\tfree_mol_main_flow_component"
                     "\tfree_mol_boundary_connected\tfree_mol_boundaryless"
                     "\tgas_volume_total\tgas_volume_main_flow_component"
                     "\tgas_volume_boundary_connected\tgas_volume_boundaryless"
                     "\tdissolved_mol_total\tboundary_advective_component_outflow"
                     "\tboundary_diffusive_component_outflow"
                     "\tinlet_liquid_inflow\toutlet_liquid_outflow\tnet_liquid_inflow"
                     "\tthroughflow_rate\tcumulative_net_liquid_inflow"
                     "\tliquid_volume_change\tcumulative_liquid_volume_change"
                     "\twrong_direction_boundary_edges"
                     "\tinlet_reverse_flow_edges\toutlet_backflow_edges"
                     "\tinlet_reverse_flow_m3_s\toutlet_backflow_m3_s"
                     "\tinternal_pressure_min_Pa\tinternal_pressure_max_Pa"
                     "\tmaximum_dissolution_driving_force"
                     "\tactive_main_flow\tactive_boundary_connected"
                     "\tactive_boundaryless\tcumulative_boundary_component_outflow_mol"
                     "\tinlet_advective_component_inflow\toutlet_advective_component_outflow"
                     "\tcumulative_outlet_advective_component_outflow_mol"
                     "\tcumulative_component_inflow_mol\tcumulative_component_outflow_mol"
                     "\tactual_liquid_inflow\tactual_liquid_outflow\ttransfer_rate_mol_s"
                     "\tPl_relative_min_node_id\tPl_relative_max_node_id"
                     "\tPg_relative_min_Pa\tPg_relative_max_Pa\tPg_relative_min_node_id\tPg_relative_max_node_id"
                     "\tnext_dt_s\tdt_change_reason\tmin_free_gas_moles\tmin_Sg"
                     "\ttotal_free_gas_moles\ttotal_dissolved_moles\ttotal_gas_moles"
                     "\tgas_timescale_cap_s\tconcentration_timescale_cap_s"
                     "\tgas_disappearance_event\tevent_node_id\tevent_time_s\tevent_dt_s"
                     "\tevent_gas_moles_before\tevent_gas_moles_after"
                     "\tconcentration_limit_node_id\tconcentration_limit_C_mol_m3"
                     "\tconcentration_limit_dC_dt_mol_m3_s\tconcentration_reference_mol_m3"
                     "\tconcentration_tau_s\tphysical_time_s\toverall_Sg"
                     "\tmain_connected_component_Sg\tgas_disappearance_events"
                     "\tmass_transfer_multiplier\ttermination_reason"
                     "\tinlet_pressure_relative_Pa\tinlet_pressure_absolute_Pa\tinlet_flow_target_m3_s"
                     "\tretired_bubble_count\tbatched_retired_bubble_count"
                     "\tbatched_retired_moles\tbatched_retired_gas_volume_m3"
                     "\tbatch_aggregate_free_gas_fraction\tevent_timestep_localized\n";

      // Timing detail remains in global_history.tsv to keep legacy diagnostics readers stable.
    }
    const auto regularization_path = config.io.output_directory / "initial_regularization.tsv";
    if (restart_path.empty()) {
      std::ofstream regularization(regularization_path);
      regularization << "node_id\tinput_Sw\tregularized_Sw\n";
      for (const auto node : solver.regularized_nodes()) {
        regularization << node + 1 << '\t' << network.pores[node].input_sw << '\t'
                       << solver.state()[node].sw << '\n';
      }
      regularization << "# mole_correction\t" << std::setprecision(17)
                     << solver.initial_regularization_moles() << '\n';
      regularization << "# liquid_volume_correction\t"
                     << solver.initial_regularization_volume() << '\n';
      regularization << "# threshold_moles_transferred\t"
                     << solver.initial_threshold_moles_transferred() << '\n';
      const auto initial_vtk = config.io.output_directory / step_name(solver.step_number(), state_extension(config));
      write_state(initial_vtk, solver, nullptr);
      outputs.push_back(initial_vtk);
      pvd_entries.push_back({solver.time(), initial_vtk.filename().string()});
      write_pvd(pvd_path, pvd_entries, config.physics.mass_transfer_multiplier);
    }

    double dt = restart ? restart->next_dt : config.time.initial_dt;
    int consecutive_growth = restart ? restart->consecutive_fast_steps : 0;
    double previous_free = free_total(solver);
    json convergence_steps = json::array();
    double maximum_residual = 0.0;
    int maximum_residual_step = 0;
    double maximum_update_ratio = 0.0;
    int maximum_update_ratio_step = 0;
    std::array<double, 4> maximum_blocks{{0.0, 0.0, 0.0, 0.0}};
    std::array<int, 4> maximum_block_steps{{0, 0, 0, 0}};
    std::array<std::size_t, 4> maximum_block_nodes{{0, 0, 0, 0}};
    bool all_convergence_criteria_passed = true;
    std::size_t next_vtk_time = 0, next_checkpoint_time = 0;
    while (next_vtk_time < config.io.vtk_output_times.size() &&
           config.io.vtk_output_times[next_vtk_time] <= solver.time())
      ++next_vtk_time;
    while (next_checkpoint_time < config.io.checkpoint_output_times.size() &&
           config.io.checkpoint_output_times[next_checkpoint_time] <= solver.time())
      ++next_checkpoint_time;
    bool first_inactive_written = false, main_95_written = false, main_90_written = false;
    std::string termination_reason = "configured_end_time";
    json scientific_events=json::array();
    const std::array<double,6> loss_targets{{0.05,0.10,0.25,0.50,0.75,0.90}};
    std::array<bool,6> loss_written{{false,false,false,false,false,false}};
    auto previous_scientific=initial_scientific;
    std::size_t next_pv_output=0;
    int stopped_transfer_steps = 0;
    bubble::Assembly last_assembly;
    const double end_time_epsilon = 32.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(config.time.end));
    while (solver.time() + end_time_epsilon < config.time.end) {
      const double remaining = config.time.end - solver.time();
      auto report = solver.advance_with_retry(std::min(dt, remaining));
      if (!report.accepted) {
        bubble::write_checkpoint(config.io.output_directory / "checkpoint.json", solver,
            std::max(config.time.min_dt, std::min(dt, config.time.max_dt)), consecutive_growth,
            checkpoint_provenance);
        throw std::runtime_error("timestep failed at t=" + std::to_string(solver.time()) +
                                 ": " + report.failure);
      }
      last_assembly = report.final_assembly;
      accepted_steps_before_failure = solver.step_number();
      time_before_failure = solver.time();
      audit_state(solver, previous_free);
      previous_free = free_total(solver);
      auto decision=choose_next_dt(solver,report,report.dt,consecutive_growth);
      if(report.gas_disappearance_event && report.event_timestep_localized){
        // The tiny event-localized step is not a new global stability scale.
        // Resume from a conservative fraction of the pre-event requested dt;
        // the ordinary retry controller will reduce it if the coupled system
        // cannot yet accept that size.  This avoids thousands of geometric
        // growth steps after each isolated bubble disappearance.
        decision.next=std::clamp(
            std::max(report.dt,
                     report.requested_dt*config.time_step_control.shrink_factor),
            config.time.min_dt,config.time.max_dt);
        decision.reason="gas_disappearance_event";
      }
      report.next_dt=decision.next;report.gas_timescale_cap=decision.gas_cap;
      report.concentration_timescale_cap=decision.concentration_cap;report.dt_change_reason=decision.reason;
      report.concentration_limit_node_id=decision.concentration_node;
      report.concentration_limit_C=decision.concentration_C;report.concentration_limit_dC_dt=decision.concentration_dCdt;
      report.concentration_reference=decision.concentration_reference;report.concentration_tau=decision.concentration_tau;
      if(config.time_step_control.enabled){
        if(decision.reason=="easy_newton_growth")++consecutive_growth;else consecutive_growth=0;
      }else if(decision.reason=="easy_newton_growth")consecutive_growth=0;
      else if(report.newton_iterations<=config.time.fast_newton_iterations)++consecutive_growth;
      else consecutive_growth=0;
      solver.set_adaptive_metadata(decision.reason,report.min_free_gas_moles,report.min_Sg);
      const auto current_scientific=scientific.snapshot(solver);
      const bool gas_saturation_target_enabled=config.termination.enabled &&
          config.termination.primary!="residual_free_gas_fraction";
      const bool target_reached=gas_saturation_target_enabled &&
          current_scientific.main_real_gas_saturation<=config.termination.target_gas_saturation;
      const bool free_fraction_target_enabled=config.termination.enabled &&
          config.termination.residual_free_gas_fraction>0.0 &&
          config.termination.primary!="main_real_gas_saturation";
      const bool free_fraction_reached=free_fraction_target_enabled &&
          current_scientific.free_gas_mol<=scientific.initial_free_moles()*
              config.termination.residual_free_gas_fraction;
      const char* step_termination_reason=free_fraction_reached?
          "residual_free_gas_fraction_reached":
          (target_reached?"main_real_gas_saturation_reached":"running");
      diagnostics << solver.step_number() << '\t' << solver.time() << '\t'
                  << report.requested_dt << '\t' << report.dt << '\t'
                  << report.timestep_retries << '\t'
                  << report.newton_iterations << '\t' << report.linear_iterations << '\t'
                  << report.linear_residual << '\t' << report.residual_norm << '\t'
                  << report.state_norm << '\t' << report.update_norm << '\t'
                  << report.relative_update_norm << '\t' << report.update_tolerance_limit << '\t'
                  << report.update_criterion_ratio << '\t'
                  << static_cast<int>(report.residual_criterion_passed) << '\t'
                  << static_cast<int>(report.update_criterion_passed);
      for (double value : report.final_assembly.residual_norms) diagnostics << '\t' << value;
      for (std::size_t node_id : report.final_assembly.residual_node_ids)
        diagnostics << '\t' << node_id;
      diagnostics << '\t' << report.balance.boundary_liquid_outflow << '\t'
                  << report.balance.total_dissolved << '\t' << report.balance.total_free << '\t'
                  << report.balance.boundary_molar_outflow << '\t'
                  << report.balance.epsilon_moles_equation << '\t'
                  << report.balance.epsilon_volume_equation << '\t'
                  << report.balance.epsilon_moles_adjusted << '\t'
                  << report.balance.epsilon_volume_adjusted << '\t' << report.active_after << '\t'
                  << report.threshold_moles_transferred << '\t' << report.active_set_mole_ledger << '\t'
                  << report.active_set_liquid_volume_ledger << '\t'
                  << report.balance.total_free << '\t' << report.balance.free_main_flow << '\t'
                  << report.balance.free_boundary_connected << '\t'
                  << report.balance.free_boundaryless << '\t'
                  << report.balance.gas_volume_total << '\t'
                  << report.balance.gas_volume_main_flow << '\t'
                  << report.balance.gas_volume_boundary_connected << '\t'
                  << report.balance.gas_volume_boundaryless << '\t'
                  << report.balance.total_dissolved << '\t'
                  << report.balance.boundary_advective_outflow << '\t'
                  << report.balance.boundary_diffusive_outflow << '\t'
                  << report.balance.inlet_liquid_inflow << '\t'
                  << report.balance.outlet_liquid_outflow << '\t'
                  << report.balance.net_liquid_inflow << '\t'
                  << report.balance.throughflow_rate << '\t'
                  << solver.cumulative_net_liquid_inflow() << '\t'
                  << report.balance.liquid_volume_change << '\t'
                  << report.balance.cumulative_liquid_volume_change << '\t'
                  << report.balance.wrong_direction_boundary_edges << '\t'
                  << report.balance.inlet_reverse_flow_edges << '\t'
                  << report.balance.outlet_backflow_edges << '\t'
                  << report.balance.inlet_reverse_flow << '\t'
                  << report.balance.outlet_backflow << '\t'
                  << report.balance.internal_pressure_min << '\t'
                  << report.balance.internal_pressure_max << '\t'
                  << report.balance.maximum_dissolution_driving_force << '\t'
                  << report.balance.active_main_flow << '\t'
                  << report.balance.active_boundary_connected << '\t'
                  << report.balance.active_boundaryless << '\t'
                  << solver.cumulative_boundary_component_outflow() << '\t'
                  << report.balance.inlet_advective_component_inflow << '\t'
                  << report.balance.outlet_advective_component_outflow << '\t'
                  << solver.cumulative_outlet_advective_component_outflow() << '\t'
                  << solver.cumulative_component_inflow() << '\t'
                  << solver.cumulative_component_outflow() << '\t'
                  << report.balance.liquid_inflow << '\t' << report.balance.liquid_outflow << '\t'
                  << report.balance.transfer_rate << '\t'
                  << report.balance.pl_min_id << '\t' << report.balance.pl_max_id << '\t'
                  << report.balance.pg_min << '\t' << report.balance.pg_max << '\t'
                  << report.balance.pg_min_id << '\t' << report.balance.pg_max_id << '\t'
                  << report.next_dt << '\t' << report.dt_change_reason << '\t'
                  << report.min_free_gas_moles << '\t' << report.min_Sg << '\t'
                  << report.balance.total_free << '\t' << report.balance.total_dissolved << '\t'
                  << report.balance.total_free + report.balance.total_dissolved << '\t'
                  << report.gas_timescale_cap << '\t' << report.concentration_timescale_cap << '\t'
                  << static_cast<int>(report.gas_disappearance_event) << '\t'
                  << report.event_node_id << '\t' << report.event_time << '\t' << report.event_dt << '\t'
                  << report.event_gas_moles_before << '\t' << report.event_gas_moles_after << '\t'
                  << report.concentration_limit_node_id << '\t' << report.concentration_limit_C << '\t'
                  << report.concentration_limit_dC_dt << '\t' << report.concentration_reference << '\t'
                  << report.concentration_tau << '\t' << solver.time() << '\t'
                  << current_scientific.mean_sg << '\t'
                  << current_scientific.main_real_gas_saturation << '\t'
                  << solver.gas_disappearance_events() << '\t'
                  << config.physics.mass_transfer_multiplier << '\t'
                  << step_termination_reason << '\t'
                  << solver.inlet_pressure() << '\t'
                  << config.physics.background_abs_Pa + solver.inlet_pressure() << '\t'
                  << config.boundary.inlet_flow_m3_s << '\t'
                  << report.retired_bubble_count << '\t'
                  << report.batched_retired_bubble_count << '\t'
                  << report.batched_retired_moles << '\t'
                  << report.batched_retired_gas_volume << '\t'
                  << report.batch_aggregate_free_gas_fraction << '\t'
                  << static_cast<int>(report.event_timestep_localized) << '\n';
      diagnostics.flush();
      scientific.append(solver, report);
      const double previous_loss=1.0-previous_scientific.free_gas_mol/scientific.initial_free_moles();
      const double current_loss=1.0-current_scientific.free_gas_mol/scientific.initial_free_moles();
      bool scientific_output_event=false;
      for(std::size_t event=0;event<loss_targets.size();++event){
        if(!loss_written[event] && current_loss>=loss_targets[event]){
          const double fraction=current_loss>previous_loss?
            std::clamp((loss_targets[event]-previous_loss)/(current_loss-previous_loss),0.0,1.0):1.0;
          scientific_events.push_back({{"kind","free_gas_loss_fraction"},{"target",loss_targets[event]},
            {"estimated_crossing_time_s",solver.time()-report.dt+fraction*report.dt},
            {"estimated_crossing_injected_pore_volumes",previous_scientific.injected_pore_volumes+
              fraction*(current_scientific.injected_pore_volumes-previous_scientific.injected_pore_volumes)},
            {"accepted_state_step",solver.step_number()},{"accepted_state_time_s",solver.time()}});
          loss_written[event]=true;scientific_output_event=true;
        }
      }
      while(next_pv_output<config.analysis.output_pore_volumes.size() &&
            current_scientific.injected_pore_volumes>=config.analysis.output_pore_volumes[next_pv_output]){
        scientific_output_event=true;++next_pv_output;
      }
      if(target_reached){
        const double before=previous_scientific.main_real_gas_saturation;
        const double after=current_scientific.main_real_gas_saturation;
        const double fraction=before>after?std::clamp((before-config.termination.target_gas_saturation)/(before-after),0.0,1.0):1.0;
        scientific_events.push_back({{"kind","main_real_gas_saturation"},
          {"target",config.termination.target_gas_saturation},
          {"estimated_crossing_time_s",solver.time()-report.dt+fraction*report.dt},
          {"estimated_crossing_injected_pore_volumes",previous_scientific.injected_pore_volumes+
            fraction*(current_scientific.injected_pore_volumes-previous_scientific.injected_pore_volumes)},
          {"accepted_state_step",solver.step_number()},{"accepted_state_time_s",solver.time()}});
        scientific_output_event=true;
      }
      if(free_fraction_reached){
        const double target=config.termination.residual_free_gas_fraction;
        const double before=previous_scientific.free_gas_mol/scientific.initial_free_moles();
        const double after=current_scientific.free_gas_mol/scientific.initial_free_moles();
        const double fraction=before>after?
            std::clamp((before-target)/(before-after),0.0,1.0):1.0;
        scientific_events.push_back({{"kind","residual_free_gas_fraction"},
          {"target",target},
          {"initial_total_free_gas_moles",scientific.initial_free_moles()},
          {"current_total_free_gas_moles",current_scientific.free_gas_mol},
          {"estimated_crossing_time_s",solver.time()-report.dt+fraction*report.dt},
          {"estimated_crossing_injected_pore_volumes",previous_scientific.injected_pore_volumes+
            fraction*(current_scientific.injected_pore_volumes-previous_scientific.injected_pore_volumes)},
          {"accepted_state_step",solver.step_number()},{"accepted_state_time_s",solver.time()}});
        scientific_output_event=true;
      }
      run_log << "step=" << solver.step_number() << " time_s=" << solver.time()
              << " requested_dt_s=" << report.requested_dt << " accepted_dt_s=" << report.dt
              << " retries=" << report.timestep_retries
              << " newton_iterations=" << report.newton_iterations
              << " residual_scaled_inf=" << report.residual_norm
              << " relative_update_norm=" << report.relative_update_norm
              << " epsilon_N=" << report.balance.epsilon_moles_adjusted
              << " epsilon_V=" << report.balance.epsilon_volume_adjusted
              << " inlet_liquid_inflow=" << report.balance.inlet_liquid_inflow
              << " outlet_liquid_outflow=" << report.balance.outlet_liquid_outflow
              << " net_liquid_inflow=" << report.balance.net_liquid_inflow
              << " throughflow_rate=" << report.balance.throughflow_rate
              << " boundary_advective_outflow=" << report.balance.boundary_advective_outflow
              << " boundary_diffusive_outflow=" << report.balance.boundary_diffusive_outflow
              << " gas_volume_main=" << report.balance.gas_volume_main_flow
              << " retired_bubbles=" << report.retired_bubble_count
              << " batched_retired_bubbles=" << report.batched_retired_bubble_count
              << " event_timestep_localized=" << report.event_timestep_localized << '\n';
      for (std::size_t retry = 0; retry < report.rejected_dts.size(); ++retry)
        run_log << "  rejected_dt_s=" << report.rejected_dts[retry]
                << " reason=" << report.retry_failures[retry] << '\n';
      run_log.flush();

      json step_convergence = {
          {"step", solver.step_number()},
          {"requested_dt_s", report.requested_dt},
          {"accepted_dt_s", report.dt},
          {"timestep_retries", report.timestep_retries},
          {"rejected_dt_s", report.rejected_dts},
          {"retry_failures", report.retry_failures},
          {"residual_scaled_inf", report.residual_norm},
          {"state_scaled_inf", report.state_norm},
          {"update_scaled_inf", report.update_norm},
          {"relative_update_norm", report.relative_update_norm},
          {"update_tolerance_limit", report.update_tolerance_limit},
          {"update_criterion_ratio", report.update_criterion_ratio},
          {"residual_criterion_passed", report.residual_criterion_passed},
          {"update_criterion_passed", report.update_criterion_passed},
          {"next_dt_s",report.next_dt},{"dt_change_reason",report.dt_change_reason},
          {"gas_timescale_cap_s",report.gas_timescale_cap},
          {"concentration_timescale_cap_s",report.concentration_timescale_cap},
          {"gas_disappearance_event",report.gas_disappearance_event},
          {"event_node_id",report.event_node_id},{"event_time_s",report.event_time},
          {"event_dt_s",report.event_dt},{"event_gas_moles_before",report.event_gas_moles_before},
          {"event_gas_moles_after",report.event_gas_moles_after},
          {"retired_bubble_count",report.retired_bubble_count},
          {"batched_retired_bubble_count",report.batched_retired_bubble_count},
          {"batched_retired_moles",report.batched_retired_moles},
          {"batched_retired_gas_volume_m3",report.batched_retired_gas_volume},
          {"batch_aggregate_free_gas_fraction",report.batch_aggregate_free_gas_fraction},
          {"event_timestep_localized",report.event_timestep_localized}};
      static constexpr std::array<const char*, 4> block_names{{"Rl", "Rd", "Rg", "Rc"}};
      for (std::size_t block = 0; block < block_names.size(); ++block) {
        step_convergence[std::string(block_names[block]) + "_scaled_inf"] =
            report.final_assembly.residual_norms[block];
        step_convergence[std::string(block_names[block]) + "_max_node_id"] =
            report.final_assembly.residual_node_ids[block];
        if (report.final_assembly.residual_norms[block] > maximum_blocks[block]) {
          maximum_blocks[block] = report.final_assembly.residual_norms[block];
          maximum_block_steps[block] = solver.step_number();
          maximum_block_nodes[block] = report.final_assembly.residual_node_ids[block];
        }
      }
      convergence_steps.push_back(std::move(step_convergence));
      if (report.residual_norm > maximum_residual) {
        maximum_residual = report.residual_norm;
        maximum_residual_step = solver.step_number();
      }
      if (report.update_criterion_ratio > maximum_update_ratio) {
        maximum_update_ratio = report.update_criterion_ratio;
        maximum_update_ratio_step = solver.step_number();
      }
      all_convergence_criteria_passed = all_convergence_criteria_passed &&
          report.residual_criterion_passed && report.update_criterion_passed;

      bool scheduled_vtk = false;
      while (next_vtk_time < config.io.vtk_output_times.size() &&
             solver.time() >= config.io.vtk_output_times[next_vtk_time]) {
        scheduled_vtk = true;
        ++next_vtk_time;
      }
      const double main_fraction = solver.initial_main_gas_volume() > 0.0
          ? report.balance.gas_volume_main_flow / solver.initial_main_gas_volume() : 1.0;
      const bool first_inactive_event = !first_inactive_written &&
          report.active_after < report.active_before;
      const bool main_95_event = !main_95_written && main_fraction <= 0.95;
      const bool main_90_event = !main_90_written && main_fraction <= 0.90;
      if (first_inactive_event) first_inactive_written = true;
      if (main_95_event) main_95_written = true;
      if (main_90_event) main_90_written = true;
      if (solver.step_number() % config.io.vtk_interval_steps == 0 || scheduled_vtk ||
          first_inactive_event || main_95_event || main_90_event || scientific_output_event ||
          solver.time() >= config.time.end) {
        const auto vtk = config.io.output_directory / step_name(solver.step_number(), state_extension(config));
        write_state(vtk, solver, &report.final_assembly);
        outputs.push_back(vtk);
        if (pvd_entries.empty() || pvd_entries.back().file != vtk.filename().string())
          pvd_entries.push_back({solver.time(), vtk.filename().string()});
        write_pvd(pvd_path, pvd_entries, config.physics.mass_transfer_multiplier);
      }
      dt = report.next_dt;
      bool scheduled_checkpoint = false;
      while (next_checkpoint_time < config.io.checkpoint_output_times.size() &&
             solver.time() >= config.io.checkpoint_output_times[next_checkpoint_time]) {
        scheduled_checkpoint = true;
        ++next_checkpoint_time;
      }
      if (solver.step_number() % config.io.checkpoint_interval_steps == 0 ||
          scheduled_checkpoint || scientific_output_event ||
          solver.time() >= config.time.end) {
        const auto checkpoint = config.io.output_directory / "checkpoint.json";
        bubble::write_checkpoint(checkpoint, solver, dt, consecutive_growth,
                                 checkpoint_provenance);
        (void)bubble::read_checkpoint(checkpoint, network.pores.size());
        if (scheduled_checkpoint) {
          const auto retained = config.io.output_directory /
              checkpoint_time_name(solver.time());
          std::filesystem::copy_file(checkpoint, retained,
              std::filesystem::copy_options::overwrite_existing);
          outputs.push_back(retained);
        }
      }
      if (config.time.stop_main_gas_volume_fraction > 0.0 &&
          main_fraction <= config.time.stop_main_gas_volume_fraction) {
        termination_reason = "main_flow_gas_volume_fraction_reached";
        break;
      }
      if(free_fraction_reached){termination_reason="residual_free_gas_fraction_reached";break;}
      if(target_reached){termination_reason="main_real_gas_saturation_reached";break;}
      if (std::filesystem::exists(config.io.output_directory / "STOP")) {
        termination_reason = "safe_stop_requested"; break;
      }
      if(config.termination.enabled && config.termination.maximum_wall_seconds>0 &&
         std::chrono::duration<double>(std::chrono::steady_clock::now()-wall_started).count()>=
             config.termination.maximum_wall_seconds){termination_reason="wall_clock_limit";break;}
      if(config.termination.enabled && directory_bytes(config.io.output_directory)>=
             config.termination.maximum_output_bytes){termination_reason="output_disk_limit";break;}
      if (config.time.equilibrium_tolerance > 0.0 &&
          report.balance.maximum_dissolution_driving_force <= config.time.equilibrium_tolerance &&
          report.balance.transfer_rate <= 1e-16) ++stopped_transfer_steps;
      else stopped_transfer_steps = 0;
      if (stopped_transfer_steps >= 5) {
        termination_reason = "no_positive_transfer_for_five_accepted_steps";
        break;
      }
      previous_scientific=current_scientific;
    }
    if(config.termination.enabled && solver.time()>=config.termination.maximum_physical_time &&
       termination_reason=="configured_end_time") termination_reason="right_censored";
    if (pvd_entries.empty() || pvd_entries.back().time != solver.time()) {
      const auto vtk = config.io.output_directory / step_name(solver.step_number(), state_extension(config));
      write_state(vtk, solver, &last_assembly);
      outputs.push_back(vtk);
      pvd_entries.push_back({solver.time(), vtk.filename().string()});
      write_pvd(pvd_path, pvd_entries, config.physics.mass_transfer_multiplier);
    }
    const auto final_checkpoint = config.io.output_directory / "checkpoint.json";
    bubble::write_checkpoint(final_checkpoint, solver, dt, consecutive_growth,
                             checkpoint_provenance);
    (void)bubble::read_checkpoint(final_checkpoint, network.pores.size());
    const auto scientific_events_path=config.io.output_directory/"termination_events.json";
    const auto final_scientific=scientific.snapshot(solver);
    {std::ofstream events(scientific_events_path);events<<std::setw(2)<<json{{"events",scientific_events},
      {"termination_reason",termination_reason},{"actual_final_time_s",solver.time()},
      {"initial_total_free_gas_moles",scientific.initial_free_moles()},
      {"current_total_free_gas_moles",final_scientific.free_gas_mol},
      {"current_dissolved_gas_moles",final_scientific.dissolved_gas_mol},
      {"actual_final_overall_gas_saturation",final_scientific.mean_sg},
      {"actual_final_main_real_gas_saturation",final_scientific.main_real_gas_saturation},
      {"gas_disappearance_events",solver.gas_disappearance_events()},
      {"last_gas_disappearance_event",{{"node_id",solver.last_gas_disappearance_node_id()},
        {"time_s",solver.last_gas_disappearance_time()}}}}<<'\n';}
    diagnostics.close();

    const double elapsed_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_started).count();
    run_log << "elapsed_seconds=" << elapsed_seconds << '\n';
    run_log.close();

    const auto checkpoint_path = config.io.output_directory / "checkpoint.json";
    outputs.push_back(diagnostics_path);
    outputs.push_back(checkpoint_path);
    outputs.push_back(regularization_path);
    outputs.push_back(pvd_path);
    outputs.push_back(config.io.output_directory/"global_history.tsv");
    outputs.push_back(config.io.output_directory/"segment_history.tsv");
    outputs.push_back(config.io.output_directory/"segment_geometry.tsv");
    outputs.push_back(config.io.output_directory/"constriction_geometry_audit.tsv");
    outputs.push_back(config.io.output_directory/"balance_history.tsv");
    outputs.push_back(config.io.output_directory/"gas_disappearance_events.tsv");
    outputs.push_back(config.io.output_directory/"BALANCE_README.md");
    outputs.push_back(scientific_events_path);
    for (int step = 0; step <= solver.step_number(); ++step) {
      const auto vtk = config.io.output_directory / step_name(step, state_extension(config));
      if (std::filesystem::is_regular_file(vtk)) outputs.push_back(vtk);
    }
    json manifest;
    manifest["case_name"] = config.case_info.name;
    manifest["case_kind"] = config.case_info.kind;
    manifest["scientific_parameters_confirmed"] = config.case_info.scientific_parameters_confirmed;
    manifest["git_commit"] = BUBBLE_GIT_COMMIT;
    manifest["build"] = {
        {"build_type", BUBBLE_BUILD_TYPE},
        {"cxx_standard", 17},
        {"sanitizers_enabled", static_cast<bool>(BUBBLE_SANITIZERS_ENABLED)},
        {"amgx_enabled", static_cast<bool>(BUBBLE_AMGX_ENABLED)},
        {"pardiso_enabled", static_cast<bool>(BUBBLE_PARDISO_ENABLED)},
        {"project_source_tree_sha256", BUBBLE_PROJECT_SOURCE_TREE_SHA256}};
    manifest["compiler"] = {BUBBLE_COMPILER_ID, BUBBLE_COMPILER_VERSION};
    const auto amgx=bubble::amgx_runtime_info();
    const auto pardiso=bubble::pardiso_runtime_info(config.linear.pardiso_threads);
    manifest["linear_backend"]={{"selected",config.linear.backend},
      {"relative_tolerance",config.linear.relative_tolerance},
      {"maximum_iterations",config.linear.max_iterations},
      {"pardiso_compiled",pardiso.compiled},{"pardiso_threads",pardiso.threads},
      {"pardiso_library_path",pardiso.library_path},
      {"amgx_compiled",amgx.compiled},{"amgx_api_version",{amgx.api_major,amgx.api_minor}},
      {"amgx_mode",amgx.mode},{"amgx_library_path",amgx.library_path},
      {"amgx_config_file",config.linear.amgx_config_file.string()},
      {"amgx_config_sha256",config.linear.amgx_config_file.empty()?"":bubble::sha256_file(config.linear.amgx_config_file)},
      {"gpu_inventory",gpu_inventory()},
      {"final_jacobian_rows",last_assembly.jacobian.rows()},
      {"final_jacobian_nonzeros",last_assembly.jacobian.nonZeros()},
      {"csr_index_base",0},{"csr_index_type","signed_int32_checked"},
      {"cpu_residual_recomputed",true}};
    manifest["termination"]={{"primary",config.termination.primary},
      {"target_gas_saturation",config.termination.target_gas_saturation},
      {"residual_free_gas_fraction",config.termination.residual_free_gas_fraction},
      {"maximum_physical_time_s",config.termination.maximum_physical_time},
      {"maximum_wall_seconds",config.termination.maximum_wall_seconds},
      {"maximum_output_bytes",config.termination.maximum_output_bytes},
      {"initial_total_free_gas_moles",scientific.initial_free_moles()},
      {"current_total_free_gas_moles",final_scientific.free_gas_mol},
      {"current_dissolved_gas_moles",final_scientific.dissolved_gas_mol},
      {"overall_gas_saturation",final_scientific.mean_sg},
      {"main_connected_component_gas_saturation",final_scientific.main_real_gas_saturation},
      {"gas_disappearance_events",solver.gas_disappearance_events()},
      {"last_gas_disappearance_event",{{"node_id",solver.last_gas_disappearance_node_id()},
        {"time_s",solver.last_gas_disappearance_time()}}},
      {"actual_reason",termination_reason}};
    manifest["mass_transfer"]={{"configured_kL_m_s",config.physics.mass_transfer_coefficient},
      {"mass_transfer_multiplier",config.physics.mass_transfer_multiplier},
      {"effective_kL_m_s",config.physics.mass_transfer_coefficient*
          config.physics.mass_transfer_multiplier},
      {"coupling","same interphase flux enters Rd with negative sign and Rg with positive sign"}};
    manifest["time_step_control"]={{"enabled",config.time_step_control.enabled},
      {"growth_factor",config.time_step_control.growth_factor},
      {"shrink_factor",config.time_step_control.shrink_factor},
      {"easy_newton_iterations",config.time_step_control.easy_newton_iterations},
      {"normal_newton_iterations",config.time_step_control.normal_newton_iterations},
      {"hard_newton_iterations",config.time_step_control.hard_newton_iterations},
      {"gas_timescale_fraction",config.time_step_control.gas_timescale_fraction},
      {"gas_timescale_percentile",config.time_step_control.gas_timescale_percentile},
      {"gas_timescale_min_mole_fraction",config.time_step_control.gas_timescale_min_mole_fraction},
      {"concentration_timescale_fraction",config.time_step_control.concentration_timescale_fraction},
      {"concentration_change_limit_enabled",config.time_step_control.concentration_change_limit_enabled},
      {"concentration_change_fraction",config.time_step_control.concentration_change_fraction},
      {"concentration_reference_definition","Hcp_mol_m3_Pa * background_abs_Pa"},
      {"max_consecutive_growth",config.time_step_control.max_consecutive_growth}};
    manifest["active_set"]={{"Sw_min",config.active_set.sw_min},
      {"Sg_off",config.active_set.sg_off},{"ng_off_mol",config.active_set.ng_off},
      {"transfer_switch_tolerance_mol_m3",config.active_set.transfer_switch_tolerance},
      {"allow_bubble_reactivation",false},
      {"batch_microbubble_retirement",{
        {"enabled",config.active_set.batch_microbubble_retirement_enabled},
        {"maximum_individual_free_gas_fraction",
          config.active_set.batch_individual_free_gas_fraction},
        {"maximum_batch_free_gas_fraction_per_step",
          config.active_set.batch_total_free_gas_fraction_per_step},
        {"mass_handling","remaining positive free gas is transferred to cluster dissolved inventory"},
        {"volume_handling","retired gas volume enters active-set liquid-volume ledger"},
        {"event_log","gas_disappearance_events.tsv"}}}};
    manifest["dependencies"] = {
        {"Eigen", {{"version", BUBBLE_EIGEN_VERSION},
                    {"header_tree_sha256", BUBBLE_EIGEN_TREE_SHA256}}},
        {"FADBAD++", {{"version", "2.1"},
                       {"distribution_tree_sha256", BUBBLE_FADBAD_TREE_SHA256}}},
        {"nlohmann_json", {{"version", {NLOHMANN_JSON_VERSION_MAJOR,
                                         NLOHMANN_JSON_VERSION_MINOR,
                                         NLOHMANN_JSON_VERSION_PATCH}},
                            {"header_tree_sha256", BUBBLE_NLOHMANN_TREE_SHA256}}}};
    manifest["fadbad_headers"] = {
        {"include_directory", BUBBLE_FADBAD_INCLUDE_DIR},
        {"badiff_h_sha256", bubble::sha256_file(std::filesystem::path(BUBBLE_FADBAD_INCLUDE_DIR) / "badiff.h")},
        {"fadiff_h_sha256", bubble::sha256_file(std::filesystem::path(BUBBLE_FADBAD_INCLUDE_DIR) / "fadiff.h")}};
    manifest["input_unit_conversions"] = {
        {"input_basis", config.geometry.input_basis},
        {"input_length_unit", config.geometry.input_length_unit},
        {"voxel_size_um", config.geometry.voxel_size_um},
        {"image_dimensions_voxels", config.geometry.image_dimensions_voxels},
        {"physical_domain_um", config.geometry.physical_domain_um},
        {"input_length_to_m", config.geometry.length_to_m},
        {"input_area_to_m2", config.geometry.area_to_m2},
        {"input_volume_to_m3", config.geometry.volume_to_m3}};
    manifest["inputs"] = {
        {config.io.pore_file.string(), bubble::sha256_file(config.io.pore_file)},
        {config.io.throat_file.string(), bubble::sha256_file(config.io.throat_file)},
        {config.io.connect_audit_file.string(), bubble::sha256_file(config.io.connect_audit_file)},
        {config_path.string(), bubble::sha256_file(config_path)}};
    const std::array<std::pair<const char*, std::filesystem::path>, 7> provenance_files{{
        {"raw_pore_file", config.input_provenance.raw_pore_file},
        {"raw_throat_file", config.input_provenance.raw_throat_file},
        {"base_pore_file", config.input_provenance.base_pore_file},
        {"base_throat_file", config.input_provenance.base_throat_file},
        {"base_connect_file", config.input_provenance.base_connect_file},
        {"input_geometry_audit_file", config.io.input_geometry_audit_file},
        {"virtual_boundary_geometry_audit_file",
         config.io.virtual_boundary_geometry_audit_file}}};
    for (const auto& item : provenance_files) {
      if (!item.second.empty() && std::filesystem::is_regular_file(item.second))
        manifest["input_provenance"][item.first] = {
            {"path", std::filesystem::absolute(item.second).string()},
            {"sha256", bubble::sha256_file(item.second)}};
    }
    manifest["run_command"] = invoked_as;
    manifest["run_argv"] = invoked_argv;
    if (restart) {
      manifest["restart"] = {
          {"path", std::filesystem::absolute(restart_path).string()},
          {"sha256", restart_input_hash},
          {"time_s", restart->time},
          {"step_number", restart->step_number},
          {"next_dt_s", restart->next_dt},
          {"consecutive_fast_steps", restart->consecutive_fast_steps}};
    }
    manifest["started_utc"] = started;
    manifest["ended_utc"] = bubble::utc_timestamp();
    manifest["elapsed_seconds"] = elapsed_seconds;
    manifest["termination_reason"] = termination_reason;
    manifest["pressure_reference"] = {{"representation", "relative_to_background"},
        {"background_abs_Pa", config.physics.background_abs_Pa},
        {"liquid_unknown_scale_Pa", config.scaling.pressure},
        {"gas_unknown_scale_Pa", config.scaling.gas_pressure},
        {"Rc_residual_scale_Pa", config.scaling.pressure}};
    manifest["exit_status"] = 0;
    manifest["network"] = {{"control_volumes", network.pores.size()},
                            {"half_edges", network.edges.size()},
                            {"effective_zero_distance_clusters", network.zero_distance_clusters.size()},
                            {"production_positive_edges", network.production_positive_edges},
                            {"zero_distance_constraints", network.zero_distance_edges.size()},
                            {"massless_junctions", network.massless_junctions.size()},
                            {"virtual_boundary_nodes", network.virtual_boundary_nodes},
                            {"boundary_edges", network.boundary_edges},
                            {"base_node_count", network.pores.size() - network.virtual_boundary_nodes},
                            {"base_edge_count", network.edges.size() - network.boundary_edges},
                            {"throats_with_two_zero_halves", network.throats_with_two_zero_halves},
                            {"ignored_boundary_to_boundary_edges", network.ignored_boundary_edges},
                            {"connected_components", network.connected_components},
                            {"boundaryless_components", network.boundaryless_components},
                            {"isolated_pores", network.isolated_pores}};
    manifest["network"]["nontrivial_zero_distance_clusters"] = std::count_if(
        network.zero_distance_clusters.begin(), network.zero_distance_clusters.end(),
        [](const auto& cluster) { return cluster.members.size() > 1; });
    manifest["degenerate_geometry_formulation"] =
        "raw geometry unchanged; union-find zero-distance elimination plus massless junctions";
    manifest["virtual_boundary"] = {
        {"enabled", config.virtual_boundary.enabled},
        {"boundary_axis", config.virtual_boundary.boundary_axis},
        {"extension_voxels", config.virtual_boundary.extension_voxels},
        {"voxel_size_um", config.virtual_boundary.voxel_size_um},
        {"boundary_throat_length_um", config.virtual_boundary.extension_voxels *
             config.virtual_boundary.voxel_size_um},
        {"inlet_virtual_nodes", std::count_if(network.pores.begin(), network.pores.end(),
             [&](const auto& pore) { return pore.source_type == 2 &&
                 pore.boundary_flag == config.boundary.inlet_flag; })},
        {"outlet_virtual_nodes", std::count_if(network.pores.begin(), network.pores.end(),
             [&](const auto& pore) { return pore.source_type == 2 &&
                 pore.boundary_flag == config.boundary.outlet_flag; })}};
    manifest["boundary_conditions"] = {
        {"pressure_boundary_type", "Dirichlet"},
        {"inlet_boundary_type", config.boundary.inlet_type},
        {"outlet_boundary_type", config.boundary.outlet_type},
        {"configured_inlet_flow_m3_s", config.boundary.inlet_flow_m3_s},
        {"solved_inlet_pressure_relative_Pa", solver.inlet_pressure()},
        {"solved_inlet_pressure_absolute_Pa", config.physics.background_abs_Pa + solver.inlet_pressure()},
        {"concentration_boundary_type",
         config.boundary.species_mode == "copy_adjacent_concentration"
             ? "copy_adjacent_zero_gradient" : "fixed_reservoir_concentration"},
        {"boundary_diffusive_component_flux",
         config.boundary.species_mode == "copy_adjacent_concentration" ? "zero" : "configured"},
        {"boundary_advective_component_flux", "enabled"},
        {"internal_advection", "enabled"}, {"internal_diffusion", "enabled"}};
    if (config.case_info.kind == "reference_physics_with_unvalidated_closures") {
      manifest["physical_parameter_classification"] = {
          {"reference_properties", {"T_K", "molar_mass_kg_mol", "Z",
             "Hcp_mol_m3_Pa", "Dm_m2_s", "mu_l_Pa_s", "rho_l_kg_m3",
             "sigma_N_m"}},
          {"unvalidated_closures", {"contact_angle_rad", "tortuosity",
             "kr_exponent_m", "qin_b", "constant_kL_Sh2", "Qin_effective_radius"}}};
    }
    if (!virtual_boundary_audit_path.empty()) {
      std::ifstream audit_input(virtual_boundary_audit_path);
      manifest["virtual_boundary"]["geometry_audit"] = json::parse(audit_input);
    }
    const auto final_dofs = solver.make_dof_map(solver.state());
    manifest["network"]["final_pressure_gauges"] = std::count_if(
        final_dofs.node.begin(), final_dofs.node.end(),
        [](const auto& node) { return node.pressure_gauge; });
    manifest["accepted_steps"] = solver.step_number();
    manifest["final_time_s"] = solver.time();
    manifest["configuration"] = json::parse(config.source_text);
    manifest["nonlinear_convergence_audit"] = {
        {"residual_tolerance", config.nonlinear.residual_tolerance},
        {"update_tolerance", config.nonlinear.update_tolerance},
        {"acceptance_rule",
         "residual_scaled_inf <= residual_tolerance AND update_scaled_inf <= update_tolerance * (1 + state_scaled_inf)"},
        {"all_steps_passed_residual_and_update_criteria", all_convergence_criteria_passed},
        {"maximum_residual_scaled_inf", maximum_residual},
        {"maximum_residual_step", maximum_residual_step},
        {"maximum_update_criterion_ratio", maximum_update_ratio},
        {"maximum_update_criterion_ratio_step", maximum_update_ratio_step},
        {"steps", convergence_steps}};
    static constexpr std::array<const char*, 4> manifest_block_names{{"Rl", "Rd", "Rg", "Rc"}};
    for (std::size_t block = 0; block < manifest_block_names.size(); ++block) {
      manifest["nonlinear_convergence_audit"]["block_maxima"][manifest_block_names[block]] = {
          {"scaled_inf", maximum_blocks[block]},
          {"step", maximum_block_steps[block]},
          {"node_id", maximum_block_nodes[block]}};
    }
    for (const auto& path : outputs) {
      if (std::filesystem::is_regular_file(path)) {
        const auto absolute = std::filesystem::absolute(path);
        manifest["outputs"][absolute.string()] = bubble::sha256_file(absolute);
      }
    }
    const auto manifest_path = config.io.output_directory / "run_manifest.json";
    std::ofstream manifest_out(manifest_path);
    manifest_out << std::setw(2) << manifest << '\n';
    if (!manifest_out) throw std::runtime_error("failed writing run manifest");
    manifest_out.close();
    std::filesystem::remove(config.io.output_directory / "run_manifest.failed.json");
    std::cout << "accepted_steps=" << solver.step_number() << " final_time_s=" << solver.time()
              << " output=" << config.io.output_directory << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    if (run_log) {
      run_log << "ERROR: " << error.what() << '\n';
      run_log.close();
    }
    if (loaded_config) {
      try {
        const auto& config = *loaded_config;
        std::filesystem::create_directories(config.io.output_directory);
        json manifest;
        manifest["case_name"] = config.case_info.name;
        manifest["case_kind"] = config.case_info.kind;
        manifest["scientific_parameters_confirmed"] =
            config.case_info.scientific_parameters_confirmed;
        manifest["git_commit"] = BUBBLE_GIT_COMMIT;
        manifest["build"] = {
            {"build_type", BUBBLE_BUILD_TYPE},
            {"cxx_standard", 17},
            {"sanitizers_enabled", static_cast<bool>(BUBBLE_SANITIZERS_ENABLED)},
            {"project_source_tree_sha256", BUBBLE_PROJECT_SOURCE_TREE_SHA256}};
        manifest["compiler"] = {BUBBLE_COMPILER_ID, BUBBLE_COMPILER_VERSION};
        manifest["dependencies"] = {
            {"Eigen", {{"version", BUBBLE_EIGEN_VERSION},
                        {"header_tree_sha256", BUBBLE_EIGEN_TREE_SHA256}}},
            {"FADBAD++", {{"version", "2.1"},
                           {"distribution_tree_sha256", BUBBLE_FADBAD_TREE_SHA256}}},
            {"nlohmann_json", {{"version", {NLOHMANN_JSON_VERSION_MAJOR,
                                             NLOHMANN_JSON_VERSION_MINOR,
                                             NLOHMANN_JSON_VERSION_PATCH}},
                                {"header_tree_sha256", BUBBLE_NLOHMANN_TREE_SHA256}}}};
        manifest["run_command"] = invoked_as;
        manifest["run_argv"] = invoked_argv;
        manifest["started_utc"] = started;
        manifest["ended_utc"] = bubble::utc_timestamp();
        manifest["exit_status"] = 1;
        manifest["accepted_steps"] = accepted_steps_before_failure;
        manifest["failure_time_s"] = time_before_failure;
        manifest["error"] = error.what();
        manifest["configuration"] = json::parse(config.source_text);
        const std::array<std::filesystem::path, 4> inputs{{
            config.io.pore_file, config.io.throat_file,
            config.io.connect_audit_file, config.source_path}};
        for (const auto& path : inputs) {
          if (std::filesystem::is_regular_file(path)) {
            manifest["inputs"][std::filesystem::absolute(path).string()] =
                bubble::sha256_file(path);
          }
        }
        const auto path = config.io.output_directory / "run_manifest.failed.json";
        std::ofstream output(path);
        output << std::setw(2) << manifest << '\n';
      } catch (const std::exception& manifest_error) {
        std::cerr << "ERROR writing failure manifest: " << manifest_error.what() << '\n';
      }
    }
    return 1;
  }
}
