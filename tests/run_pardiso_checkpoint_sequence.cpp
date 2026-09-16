#include "bubble/model.hpp"

#include <nlohmann/json.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using json = nlohmann::json;

void atomic_json(const std::filesystem::path& path, const json& value) {
  auto temporary = path;
  temporary += ".tmp";
  { std::ofstream stream(temporary); stream << std::setw(2) << value << '\n'; }
  std::filesystem::rename(temporary, path);
}

int main(int argc, char** argv) {
  if (argc != 8 && argc != 9)
    throw std::runtime_error(
        "usage: base_config checkpoint output steps dt backend threads [inlet_flow_m3_s]");
  const auto output = std::filesystem::path(argv[3]);
  const int count = std::stoi(argv[4]);
  const double requested_dt = std::stod(argv[5]);
  const std::string backend = argv[6];
  const int threads = std::stoi(argv[7]);
  if (backend != "eigen_bicgstab_ilut" && backend != "pardiso")
    throw std::runtime_error("backend must be eigen_bicgstab_ilut or pardiso");
  if (std::filesystem::exists(output))
    throw std::runtime_error("comparison output already exists");
  std::filesystem::create_directories(output);
  std::ifstream input(argv[1]);
  json resolved = json::parse(input);
  resolved["io"]["output_directory"] = output.string();
  resolved["linear"]["backend"] = backend;
  resolved["linear"]["pardiso_threads"] = threads;
  if (argc == 9) {
    resolved["boundary"]["inlet_type"] = "flow";
    resolved["boundary"]["outlet_type"] = "pressure";
    resolved["boundary"]["inlet_flow_m3_s"] = std::stod(argv[8]);
  }
  resolved["time"]["t_end_s"] = 1e6;
  const auto config_path = output / "resolved_config.json";
  { std::ofstream stream(config_path); stream << std::setw(2) << resolved << '\n'; }
  auto config = bubble::Config::load(config_path);
  auto network = bubble::Network::load(config);
  bubble::Solver solver(network, config);
  solver.initialize();
  const auto checkpoint = bubble::read_checkpoint(argv[2], network.pores.size());
  solver.restore(checkpoint.state, checkpoint.time, checkpoint.step_number,
      checkpoint.cumulative_mole_ledger, checkpoint.cumulative_volume_ledger,
      checkpoint.cumulative_net_liquid_inflow,
      checkpoint.cumulative_boundary_component_outflow,
      checkpoint.cumulative_outlet_advective_component_outflow,
      checkpoint.cumulative_liquid_inflow, checkpoint.cumulative_liquid_outflow,
      checkpoint.cumulative_transfer_by_node);
  solver.restore_gross_fluxes(checkpoint.cumulative_component_inflow,
                              checkpoint.cumulative_component_outflow);
  std::ofstream table(output / "steps.tsv");
  table << "accepted_index\tglobal_step\ttime_s\trequested_dt_s\taccepted_dt_s\tretries"
           "\tstep_wall_seconds\tnewton_iterations\tlinear_iterations"
           "\tcpu_linear_relative_residual\tresidual_scaled_inf\trelative_update_norm"
           "\tepsilon_N\tepsilon_V\tassembly_seconds\tcsr_seconds\tsetup_seconds"
           "\tsymbolic_analysis_seconds\tnumeric_factorization_seconds\tsolve_seconds"
           "\tsymbolic_analysis_calls\tnumeric_factorization_calls\tfactor_nonzeros"
           "\tpeak_memory_kib\tinlet_pressure_relative_Pa"
           "\tinlet_liquid_flow_m3_s\toutlet_liquid_flow_m3_s\ttransfer_rate_mol_s"
           "\tfree_gas_mol\tdissolved_gas_mol\tgas_volume_m3\n"
        << std::scientific << std::setprecision(17);
  json steps = json::array();
  const auto run_started = std::chrono::steady_clock::now();
  for (int accepted = 1; accepted <= count; ++accepted) {
    atomic_json(output / "progress.json", {{"stage", "checkpoint_sequence"},
        {"status", "running"}, {"backend", backend}, {"accepted", accepted - 1},
        {"target", count}, {"time_s", solver.time()}});
    const auto step_started = std::chrono::steady_clock::now();
    const auto report = solver.advance_with_retry(requested_dt);
    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - step_started).count();
    if (!report.accepted || !report.residual_criterion_passed ||
        !report.update_criterion_passed)
      throw std::runtime_error(backend + " sequence failed: " + report.failure);
    table << accepted << '\t' << solver.step_number() << '\t' << solver.time()
          << '\t' << requested_dt << '\t' << report.dt << '\t'
          << report.timestep_retries << '\t' << wall << '\t'
          << report.newton_iterations << '\t' << report.linear_iterations << '\t'
          << report.linear_residual << '\t' << report.residual_norm << '\t'
          << report.relative_update_norm << '\t'
          << report.balance.epsilon_moles_adjusted << '\t'
          << report.balance.epsilon_volume_adjusted << '\t'
          << report.assembly_seconds << '\t' << report.csr_seconds << '\t'
          << report.linear_setup_seconds << '\t'
          << report.linear_symbolic_analysis_seconds << '\t'
          << report.linear_numeric_factorization_seconds << '\t'
          << report.linear_solve_seconds << '\t'
          << report.pardiso_symbolic_analysis_calls << '\t'
          << report.pardiso_numeric_factorization_calls << '\t'
          << report.pardiso_factor_nonzeros << '\t'
          << report.pardiso_peak_memory_kib << '\t'
          << solver.inlet_pressure() << '\t'
          << report.balance.inlet_liquid_inflow << '\t'
          << report.balance.outlet_liquid_outflow << '\t'
          << report.balance.transfer_rate << '\t'
          << report.balance.total_free << '\t'
          << report.balance.total_dissolved << '\t'
          << report.balance.gas_volume_total << '\n';
    table.flush();
    steps.push_back({{"accepted_index", accepted}, {"time_s", solver.time()},
        {"step_wall_seconds", wall}, {"linear_residual", report.linear_residual},
        {"epsilon_N", report.balance.epsilon_moles_adjusted},
        {"epsilon_V", report.balance.epsilon_volume_adjusted}});
  }
  const auto provenance = bubble::make_checkpoint_provenance(config);
  bubble::write_checkpoint(output / "checkpoint.json", solver, requested_dt, 0,
                           provenance);
  const double total_wall = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - run_started).count();
  atomic_json(output / "run_manifest.json", {{"status", "passed"},
      {"backend", backend}, {"pardiso_threads", threads},
      {"total_wall_seconds", total_wall}, {"accepted_steps", count},
      {"initial_checkpoint", std::filesystem::absolute(argv[2]).string()},
      {"initial_checkpoint_sha256", bubble::sha256_file(argv[2])},
      {"resolved_config", config_path.string()},
      {"resolved_config_sha256", bubble::sha256_file(config_path)},
      {"final_checkpoint_sha256", bubble::sha256_file(output / "checkpoint.json")},
      {"steps", steps}});
  atomic_json(output / "progress.json", {{"stage", "complete"},
      {"status", "passed"}, {"backend", backend}, {"accepted", count},
      {"time_s", solver.time()}});
}
