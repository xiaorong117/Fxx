#include "bubble/model.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace bubble {
namespace {

using json = nlohmann::json;

template <class T>
T required(const json& object, const char* key, const std::string& section) {
  if (!object.contains(key) || object.at(key).is_null()) {
    throw std::runtime_error("missing/null required config field: " + section + "." + key);
  }
  try {
    return object.at(key).get<T>();
  } catch (const json::exception& error) {
    throw std::runtime_error("invalid config field " + section + "." + key + ": " + error.what());
  }
}

const json& section(const json& root, const char* key) {
  if (!root.contains(key) || !root.at(key).is_object()) {
    throw std::runtime_error(std::string("missing/non-object config section: ") + key);
  }
  return root.at(key);
}

void positive(double value, const char* name) {
  if (!std::isfinite(value) || value <= 0.0) {
    throw std::runtime_error(std::string(name) + " must be finite and > 0");
  }
}

void nonnegative(double value, const char* name) {
  if (!std::isfinite(value) || value < 0.0) {
    throw std::runtime_error(std::string(name) + " must be finite and >= 0");
  }
}

void open_unit_interval(double value, const char* name) {
  if (!std::isfinite(value) || value <= 0.0 || value >= 1.0) {
    throw std::runtime_error(std::string(name) + " must be finite and lie in (0,1)");
  }
}

void nonblank(const std::string& value, const char* name) {
  const bool has_text = std::any_of(value.begin(), value.end(), [](unsigned char c) {
    return !std::isspace(c);
  });
  if (!has_text) throw std::runtime_error(std::string(name) + " must be nonblank");
}

}  // namespace

Config Config::load(const std::filesystem::path& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open config: " + path.string());
  std::ostringstream buffer;
  buffer << input.rdbuf();
  Config c;
  c.source_path = path;
  c.source_text = buffer.str();

  json root;
  try {
    root = json::parse(c.source_text);
  } catch (const json::exception& error) {
    throw std::runtime_error("JSON parse failure: " + std::string(error.what()));
  }

  const auto& jc = section(root, "case");
  c.case_info.name = required<std::string>(jc, "name", "case");
  c.case_info.kind = required<std::string>(jc, "kind", "case");
  c.case_info.scientific_parameters_confirmed =
      required<bool>(jc, "scientific_parameters_confirmed", "case");
  c.case_info.gas_component = required<std::string>(jc, "gas_component", "case");
  c.case_info.liquid_composition = required<std::string>(jc, "liquid_composition", "case");
  c.case_info.notes = required<std::string>(jc, "notes", "case");
  nonblank(c.case_info.name, "case.name");
  nonblank(c.case_info.gas_component, "case.gas_component");
  nonblank(c.case_info.liquid_composition, "case.liquid_composition");
  nonblank(c.case_info.notes, "case.notes");

  const auto& ji = section(root, "io");
  c.io.pore_file = required<std::string>(ji, "pore_file", "io");
  c.io.throat_file = required<std::string>(ji, "throat_file", "io");
  c.io.connect_audit_file = required<std::string>(ji, "connect_audit_file", "io");
  c.io.output_directory = required<std::string>(ji, "output_directory", "io");
  c.io.state_format = ji.value("state_format", std::string("legacy_vtk"));
  if (ji.contains("input_geometry_audit_file"))
    c.io.input_geometry_audit_file =
        required<std::string>(ji, "input_geometry_audit_file", "io");
  if (ji.contains("virtual_boundary_geometry_audit_file"))
    c.io.virtual_boundary_geometry_audit_file =
        required<std::string>(ji, "virtual_boundary_geometry_audit_file", "io");
  c.io.checkpoint_interval_steps = required<int>(ji, "checkpoint_interval_steps", "io");
  c.io.vtk_interval_steps = required<int>(ji, "vtk_interval_steps", "io");
  if (ji.contains("vtk_output_times_s"))
    c.io.vtk_output_times = required<std::vector<double>>(ji, "vtk_output_times_s", "io");
  if (ji.contains("checkpoint_output_times_s"))
    c.io.checkpoint_output_times =
        required<std::vector<double>>(ji, "checkpoint_output_times_s", "io");
  if (c.io.pore_file.empty() || c.io.throat_file.empty() ||
      c.io.connect_audit_file.empty() || c.io.output_directory.empty())
    throw std::runtime_error("io paths must be nonempty");

  const auto& jg = section(root, "geometry");
  c.geometry.input_basis = required<std::string>(jg, "input_basis", "geometry");
  c.geometry.input_length_unit =
      required<std::string>(jg, "input_length_unit", "geometry");
  c.geometry.voxel_size_um = required<double>(jg, "voxel_size_um", "geometry");
  c.geometry.image_dimensions_voxels =
      required<std::array<int, 3>>(jg, "image_dimensions_voxels", "geometry");
  c.geometry.physical_domain_um =
      required<std::array<double, 3>>(jg, "physical_domain_um", "geometry");
  c.geometry.length_to_m = required<double>(jg, "input_length_to_m", "geometry");
  c.geometry.area_to_m2 = required<double>(jg, "input_area_to_m2", "geometry");
  c.geometry.volume_to_m3 = required<double>(jg, "input_volume_to_m3", "geometry");

  const auto& jvirtual = section(root, "virtual_boundary");
  c.virtual_boundary.enabled = required<bool>(jvirtual, "enabled", "virtual_boundary");
  c.virtual_boundary.boundary_axis =
      required<std::string>(jvirtual, "boundary_axis", "virtual_boundary");
  c.virtual_boundary.extension_voxels =
      required<double>(jvirtual, "extension_voxels", "virtual_boundary");
  c.virtual_boundary.voxel_size_um =
      required<double>(jvirtual, "voxel_size_um", "virtual_boundary");

  if (root.contains("input_provenance")) {
    const auto& provenance = section(root, "input_provenance");
    c.input_provenance.raw_pore_file =
        required<std::string>(provenance, "raw_pore_file", "input_provenance");
    c.input_provenance.raw_throat_file =
        required<std::string>(provenance, "raw_throat_file", "input_provenance");
    c.input_provenance.base_pore_file =
        required<std::string>(provenance, "base_pore_file", "input_provenance");
    c.input_provenance.base_throat_file =
        required<std::string>(provenance, "base_throat_file", "input_provenance");
    c.input_provenance.base_connect_file =
        required<std::string>(provenance, "base_connect_file", "input_provenance");
  }

  const auto& jp = section(root, "physics");
  c.physics.temperature = required<double>(jp, "T_K", "physics");
  c.physics.viscosity = required<double>(jp, "mu_l_Pa_s", "physics");
  c.physics.molecular_diffusivity = required<double>(jp, "Dm_m2_s", "physics");
  c.physics.henry_cp = required<double>(jp, "Hcp_mol_m3_Pa", "physics");
  c.physics.mass_transfer_coefficient = required<double>(jp, "kL_m_s", "physics");
  c.physics.mass_transfer_multiplier = jp.value("mass_transfer_multiplier", 1.0);
  c.physics.surface_tension = required<double>(jp, "sigma_N_m", "physics");
  c.physics.contact_angle = required<double>(jp, "contact_angle_rad", "physics");
  c.physics.molar_mass = required<double>(jp, "molar_mass_kg_mol", "physics");
  c.physics.liquid_density = jp.value("rho_l_kg_m3", 997.0);
  if (required<std::string>(jp, "compressibility_model", "physics") != "constant_Z") {
    throw std::runtime_error("only physics.compressibility_model=constant_Z is implemented");
  }
  c.physics.compressibility = required<double>(jp, "Z", "physics");
  c.physics.kr_exponent = required<double>(jp, "kr_exponent_m", "physics");
  c.physics.diffusion_area_exponent =
      required<double>(jp, "diffusion_area_exponent_alpha", "physics");
  c.physics.tortuosity = required<double>(jp, "tortuosity", "physics");
  c.physics.qin_b = required<double>(jp, "qin_b", "physics");
  if (required<std::string>(jp, "mass_transfer_mode", "physics") != "dissolution_only") {
    throw std::runtime_error("only physics.mass_transfer_mode=dissolution_only is implemented");
  }

  const auto& jb = section(root, "boundary");
  c.boundary.inlet_flag = required<int>(jb, "inlet_flag", "boundary");
  c.boundary.outlet_flag = required<int>(jb, "outlet_flag", "boundary");
  c.boundary.inlet_type = jb.value("inlet_type", std::string("pressure"));
  c.boundary.outlet_type = jb.value("outlet_type", std::string("pressure"));
  c.boundary.inlet_flow_m3_s = jb.value("inlet_flow_m3_s", 0.0);
  if (root.contains("pressure_reference")) {
    c.physics.background_abs_Pa = required<double>(section(root, "pressure_reference"),
        "background_abs_Pa", "pressure_reference");
    positive(c.physics.background_abs_Pa, "pressure_reference.background_abs_Pa");
    c.boundary.pl_in = jb.value("Pl_in_rel_Pa", 0.0);
    c.boundary.pl_out = required<double>(jb, "Pl_out_rel_Pa", "boundary");
  } else {
    c.boundary.pl_in = required<double>(jb, "Pl_in_abs_Pa", "boundary");
    c.boundary.pl_out = required<double>(jb, "Pl_out_abs_Pa", "boundary");
  }
  c.boundary.c_in = required<double>(jb, "C_in_mol_m3", "boundary");
  c.boundary.c_backflow = required<double>(jb, "C_backflow_mol_m3", "boundary");
  c.boundary.species_mode = jb.value("species_mode", std::string("fixed_reservoir_concentration"));
  c.boundary.fixed_concentrations_ignored = jb.value(
      "ignored_when_species_mode_is_copy_adjacent_concentration", false);
  if (required<std::string>(jb, "outlet_diffusive_condition", "boundary") != "zero_gradient") {
    throw std::runtime_error("outlet diffusive condition must be zero_gradient");
  }

  const auto& jinitial = section(root, "initial");
  c.initial.concentration = required<double>(jinitial, "C_initial_mol_m3", "initial");
  if (required<std::string>(jinitial, "pressure_initialization", "initial") != "linear_in_x" ||
      required<std::string>(jinitial, "Sw_source", "initial") != "pore_file") {
    throw std::runtime_error("unsupported initial-condition mode");
  }

  const auto& jt = section(root, "time");
  c.time.end = required<double>(jt, "t_end_s", "time");
  c.time.initial_dt = required<double>(jt, "dt_initial_s", "time");
  c.time.min_dt = required<double>(jt, "dt_min_s", "time");
  c.time.max_dt = required<double>(jt, "dt_max_s", "time");
  c.time.cut_factor = required<double>(jt, "cut_factor", "time");
  c.time.growth_factor = required<double>(jt, "growth_factor", "time");
  c.time.fast_newton_iterations = required<int>(jt, "fast_newton_iterations", "time");
  c.time.max_retries = required<int>(jt, "max_retries_per_step", "time");
  if (jt.contains("stop_main_gas_volume_fraction"))
    c.time.stop_main_gas_volume_fraction =
        required<double>(jt, "stop_main_gas_volume_fraction", "time");
  if (jt.contains("equilibrium_tolerance_mol_m3"))
    c.time.equilibrium_tolerance =
        required<double>(jt, "equilibrium_tolerance_mol_m3", "time");
  if (root.contains("time_step_control")) {
    const auto& jcontrol = section(root, "time_step_control");
    c.time_step_control.enabled = jcontrol.value("enabled", true);
    c.time.min_dt = jcontrol.value("dt_min_s", c.time.min_dt);
    c.time.max_dt = jcontrol.value("dt_max_s", c.time.max_dt);
    c.time_step_control.growth_factor = jcontrol.value("growth_factor", 1.25);
    c.time_step_control.shrink_factor = jcontrol.value("shrink_factor", 0.5);
    c.time_step_control.easy_newton_iterations = jcontrol.value("easy_newton_iterations", 3);
    c.time_step_control.normal_newton_iterations = jcontrol.value("normal_newton_iterations", 7);
    c.time_step_control.hard_newton_iterations = jcontrol.value("hard_newton_iterations", 12);
    c.time_step_control.gas_timescale_fraction = jcontrol.value("gas_timescale_fraction", 0.10);
    c.time_step_control.gas_timescale_percentile = jcontrol.value("gas_timescale_percentile", 0.50);
    c.time_step_control.gas_timescale_min_mole_fraction = jcontrol.value("gas_timescale_min_mole_fraction", 1.0e-5);
    c.time_step_control.concentration_timescale_fraction = jcontrol.value("concentration_timescale_fraction", 0.05);
    c.time_step_control.concentration_change_limit_enabled =
        jcontrol.value("concentration_change_limit_enabled", false);
    c.time_step_control.concentration_change_fraction =
        jcontrol.value("concentration_change_fraction",
                       c.time_step_control.concentration_timescale_fraction);
    c.time_step_control.max_consecutive_growth = jcontrol.value("max_consecutive_growth", 3);
  }

  if (root.contains("termination")) {
    const auto& term = section(root, "termination");
    c.termination.enabled = true;
    c.termination.primary = required<std::string>(term, "primary", "termination");
    c.termination.target_gas_saturation =
        required<double>(term, "target_gas_saturation", "termination");
    c.termination.maximum_physical_time =
        required<double>(term, "maximum_physical_time_s", "termination");
    c.termination.maximum_wall_seconds =
        required<double>(term, "maximum_wall_seconds", "termination");
    c.termination.maximum_output_bytes =
        required<std::uint64_t>(term, "maximum_output_bytes", "termination");
    c.termination.residual_free_gas_fraction =
        term.value("residual_free_gas_fraction", 0.0);
    c.time.end = c.termination.maximum_physical_time;
  }
  if (root.contains("analysis")) {
    const auto& analysis = section(root, "analysis");
    c.analysis.x_segments = required<int>(analysis, "x_segments", "analysis");
    if (analysis.contains("vtk_output_pore_volumes"))
      c.analysis.output_pore_volumes = required<std::vector<double>>(
          analysis, "vtk_output_pore_volumes", "analysis");
  }

  const auto& jn = section(root, "nonlinear");
  c.nonlinear.max_iterations = required<int>(jn, "max_iterations", "nonlinear");
  c.nonlinear.residual_tolerance = required<double>(jn, "residual_tolerance", "nonlinear");
  c.nonlinear.update_tolerance = required<double>(jn, "update_tolerance", "nonlinear");
  c.nonlinear.armijo_c1 = required<double>(jn, "armijo_c1", "nonlinear");
  c.nonlinear.line_search_reduction =
      required<double>(jn, "line_search_reduction", "nonlinear");
  c.nonlinear.line_search_min_lambda =
      required<double>(jn, "line_search_min_lambda", "nonlinear");
  c.nonlinear.fraction_to_boundary =
      required<double>(jn, "fraction_to_boundary", "nonlinear");

  const auto& jl = section(root, "linear");
  c.linear.backend = required<std::string>(jl, "backend", "linear");
  c.linear.relative_tolerance = required<double>(jl, "relative_tolerance", "linear");
  c.linear.max_iterations = required<int>(jl, "max_iterations", "linear");
  c.linear.ilut_drop_tolerance = required<double>(jl, "ilut_drop_tolerance", "linear");
  c.linear.ilut_fill_factor = required<int>(jl, "ilut_fill_factor", "linear");
  c.linear.pardiso_threads = jl.value("pardiso_threads", 16);
  if (jl.contains("amgx_config_file"))
    c.linear.amgx_config_file = required<std::string>(jl, "amgx_config_file", "linear");

  const auto& ja = section(root, "active_set");
  c.active_set.sw_min = required<double>(ja, "Sw_min", "active_set");
  c.active_set.sg_off = required<double>(ja, "Sg_off", "active_set");
  c.active_set.ng_off = required<double>(ja, "ng_off_mol", "active_set");
  c.active_set.transfer_switch_tolerance =
      required<double>(ja, "transfer_switch_tolerance_mol_m3", "active_set");
  if (required<bool>(ja, "allow_bubble_reactivation", "active_set")) {
    throw std::runtime_error("bubble reactivation is prohibited");
  }
  if (ja.contains("batch_microbubble_retirement")) {
    const auto& batch = section(ja, "batch_microbubble_retirement");
    c.active_set.batch_microbubble_retirement_enabled =
        required<bool>(batch, "enabled", "active_set.batch_microbubble_retirement");
    c.active_set.batch_individual_free_gas_fraction = required<double>(
        batch, "maximum_individual_free_gas_fraction",
        "active_set.batch_microbubble_retirement");
    c.active_set.batch_total_free_gas_fraction_per_step = required<double>(
        batch, "maximum_batch_free_gas_fraction_per_step",
        "active_set.batch_microbubble_retirement");
  }

  const auto& js = section(root, "scaling");
  c.scaling.pressure = required<double>(js, "P_ref_Pa", "scaling");
  c.scaling.gas_pressure = js.value("Pg_ref_Pa", c.scaling.pressure);
  c.scaling.concentration = required<double>(js, "C_ref_mol_m3", "scaling");
  c.scaling.volume_rate = required<double>(js, "Q_ref_m3_s", "scaling");
  c.scaling.molar_rate = required<double>(js, "Ndot_ref_mol_s", "scaling");
  c.scaling.mole_total = required<double>(js, "N_scale_mol", "scaling");
  c.scaling.volume_total = required<double>(js, "V_scale_m3", "scaling");

  const auto& jac = section(root, "acceptance");
  c.acceptance.mass_balance_hard_limit =
      required<double>(jac, "mass_balance_hard_limit", "acceptance");
  c.acceptance.volume_balance_hard_limit =
      required<double>(jac, "volume_balance_hard_limit", "acceptance");
  if (!required<bool>(jac, "reject_step_on_balance_failure", "acceptance")) {
    throw std::runtime_error("reject_step_on_balance_failure must be true");
  }

  const auto& notes = section(root, "source_notes");
  nonblank(required<std::string>(notes, "physics", "source_notes"),
           "source_notes.physics");
  nonblank(required<std::string>(notes, "boundary", "source_notes"),
           "source_notes.boundary");
  nonblank(required<std::string>(notes, "initial", "source_notes"),
           "source_notes.initial");
  c.validate_before_mesh();
  return c;
}

void Config::validate_before_mesh() const {
  if (case_info.kind == "research" && !case_info.scientific_parameters_confirmed) {
    throw std::runtime_error("research configuration is not scientifically confirmed");
  }
  if (case_info.kind != "research" && case_info.kind != "numerical_smoke_only" &&
      case_info.kind != "deterministic_test" &&
      case_info.kind != "reference_physics_with_unvalidated_closures") {
    throw std::runtime_error("unsupported case.kind");
  }
  positive(geometry.length_to_m, "geometry.input_length_to_m");
  positive(geometry.area_to_m2, "geometry.input_area_to_m2");
  positive(geometry.volume_to_m3, "geometry.input_volume_to_m3");
  if (geometry.input_basis != "physical")
    throw std::runtime_error("geometry.input_basis must be explicit physical");
  const double expected_length_scale = geometry.input_length_unit == "m" ? 1.0 :
      (geometry.input_length_unit == "um" ? 1.0e-6 :
       (geometry.input_length_unit == "nm" ? 1.0e-9 : 0.0));
  if (expected_length_scale == 0.0 ||
      !std::isfinite(geometry.length_to_m) ||
      std::abs(geometry.length_to_m - expected_length_scale) >
          1e-15 * expected_length_scale)
    throw std::runtime_error("geometry input length unit/conversion mismatch");
  if (std::abs(geometry.area_to_m2 - geometry.length_to_m * geometry.length_to_m) >
          1e-12 * geometry.area_to_m2)
    throw std::runtime_error("geometry area conversion must equal length conversion squared");
  if (std::abs(geometry.volume_to_m3 - geometry.length_to_m * geometry.length_to_m *
                  geometry.length_to_m) > 1e-12 * geometry.volume_to_m3)
    throw std::runtime_error("geometry volume conversion must equal length conversion cubed");
  positive(geometry.voxel_size_um, "geometry.voxel_size_um");
  for (std::size_t axis = 0; axis < 3; ++axis) {
    if (geometry.image_dimensions_voxels[axis] <= 0)
      throw std::runtime_error("geometry image dimensions must be positive");
    positive(geometry.physical_domain_um[axis], "geometry.physical_domain_um");
    const double expected_domain = geometry.voxel_size_um *
        static_cast<double>(geometry.image_dimensions_voxels[axis]);
    if (std::abs(geometry.physical_domain_um[axis] - expected_domain) >
        1e-12 * expected_domain)
      throw std::runtime_error("geometry physical domain must equal image dimensions times voxel size");
  }
  if (virtual_boundary.boundary_axis != "x" &&
      virtual_boundary.boundary_axis != "y" &&
      virtual_boundary.boundary_axis != "z")
    throw std::runtime_error("virtual_boundary.boundary_axis must be x, y, or z");
  positive(virtual_boundary.extension_voxels, "virtual_boundary.extension_voxels");
  positive(virtual_boundary.voxel_size_um, "virtual_boundary.voxel_size_um");
  if (std::abs(virtual_boundary.voxel_size_um - geometry.voxel_size_um) >
      1e-12 * geometry.voxel_size_um)
    throw std::runtime_error("virtual boundary and geometry voxel sizes must match");
  if (virtual_boundary.enabled &&
      (io.input_geometry_audit_file.empty() ||
       io.virtual_boundary_geometry_audit_file.empty() ||
       input_provenance.raw_pore_file.empty() ||
       input_provenance.raw_throat_file.empty() ||
       input_provenance.base_pore_file.empty() ||
       input_provenance.base_throat_file.empty() ||
       input_provenance.base_connect_file.empty()))
    throw std::runtime_error("enabled virtual boundary requires audit and input provenance paths");
  positive(physics.temperature, "physics.T_K");
  positive(physics.viscosity, "physics.mu_l_Pa_s");
  positive(physics.molecular_diffusivity, "physics.Dm_m2_s");
  positive(physics.henry_cp, "physics.Hcp_mol_m3_Pa");
  positive(physics.mass_transfer_coefficient, "physics.kL_m_s");
  positive(physics.mass_transfer_multiplier, "physics.mass_transfer_multiplier");
  positive(physics.surface_tension, "physics.sigma_N_m");
  if (!std::isfinite(physics.contact_angle) || physics.contact_angle < 0.0 ||
      physics.contact_angle > pi)
    throw std::runtime_error("physics.contact_angle_rad must lie in [0,pi]");
  positive(physics.molar_mass, "physics.molar_mass_kg_mol");
  positive(physics.liquid_density, "physics.rho_l_kg_m3");
  positive(physics.compressibility, "physics.Z");
  positive(physics.kr_exponent, "physics.kr_exponent_m");
  positive(physics.diffusion_area_exponent, "physics.diffusion_area_exponent_alpha");
  if (!std::isfinite(physics.tortuosity) || physics.tortuosity < 1.0)
    throw std::runtime_error("physics.tortuosity must be finite and >= 1");
  positive(physics.qin_b, "physics.qin_b");
  positive(physics.background_abs_Pa + boundary.pl_in, "boundary.Pl_in_abs_Pa");
  positive(physics.background_abs_Pa + boundary.pl_out, "boundary.Pl_out_abs_Pa");
  nonnegative(boundary.c_in, "boundary.C_in_mol_m3");
  nonnegative(boundary.c_backflow, "boundary.C_backflow_mol_m3");
  if (boundary.species_mode != "fixed_reservoir_concentration" &&
      boundary.species_mode != "copy_adjacent_concentration")
    throw std::runtime_error("unsupported boundary.species_mode");
  if (boundary.species_mode == "copy_adjacent_concentration") {
    if (!virtual_boundary.enabled)
      throw std::runtime_error("copy-adjacent species mode requires virtual boundaries");
    if (!boundary.fixed_concentrations_ignored)
      throw std::runtime_error("copy-adjacent species mode must mark fixed concentrations ignored");
  } else if (boundary.fixed_concentrations_ignored) {
    throw std::runtime_error("fixed species mode must not mark concentrations ignored");
  }
  nonnegative(initial.concentration, "initial.C_initial_mol_m3");
  positive(time.end, "time.t_end_s");
  positive(time.initial_dt, "time.dt_initial_s");
  positive(time.min_dt, "time.dt_min_s");
  positive(time.max_dt, "time.dt_max_s");
  if (!(time.min_dt <= time.initial_dt && time.initial_dt <= time.max_dt))
    throw std::runtime_error("time steps must satisfy dt_min <= dt_initial <= dt_max");
  open_unit_interval(time.cut_factor, "time.cut_factor");
  if (!std::isfinite(time.growth_factor) || time.growth_factor < 1.0)
    throw std::runtime_error("time.growth_factor must be finite and >= 1");
  if (time.fast_newton_iterations < 0 || time.max_retries < 0 ||
      nonlinear.max_iterations <= 0 || linear.max_iterations <= 0)
    throw std::runtime_error("iteration limits must be positive");
  positive(nonlinear.residual_tolerance, "nonlinear.residual_tolerance");
  positive(nonlinear.update_tolerance, "nonlinear.update_tolerance");
  open_unit_interval(nonlinear.armijo_c1, "nonlinear.armijo_c1");
  open_unit_interval(nonlinear.line_search_reduction,
                     "nonlinear.line_search_reduction");
  if (!std::isfinite(nonlinear.line_search_min_lambda) ||
      nonlinear.line_search_min_lambda <= 0.0 ||
      nonlinear.line_search_min_lambda > 1.0)
    throw std::runtime_error("nonlinear.line_search_min_lambda must lie in (0,1]");
  if (!(nonlinear.fraction_to_boundary > 0.0 && nonlinear.fraction_to_boundary < 1.0))
    throw std::runtime_error("fraction_to_boundary must lie in (0,1)");
  if (linear.backend != "eigen_sparselu" &&
      linear.backend != "eigen_bicgstab_ilut" && linear.backend != "amgx" &&
      linear.backend != "pardiso")
    throw std::runtime_error("unsupported linear.backend");
#if !BUBBLE_AMGX_ENABLED
  if (linear.backend == "amgx")
    throw std::runtime_error(
        "linear.backend=amgx requested, but this binary was built with BUBBLE_ENABLE_AMGX=OFF");
#endif
#if !BUBBLE_PARDISO_ENABLED
  if (linear.backend == "pardiso")
    throw std::runtime_error(
        "linear.backend=pardiso requested, but this binary was built with "
        "BUBBLE_ENABLE_PARDISO=OFF");
#endif
  if (linear.backend == "amgx" && linear.amgx_config_file.empty())
    throw std::runtime_error("linear.amgx_config_file is required for AMGX");
  positive(linear.relative_tolerance, "linear.relative_tolerance");
  nonnegative(linear.ilut_drop_tolerance, "linear.ilut_drop_tolerance");
  if (linear.ilut_fill_factor <= 0)
    throw std::runtime_error("linear.ilut_fill_factor must be positive");
  if (linear.pardiso_threads <= 0)
    throw std::runtime_error("linear.pardiso_threads must be positive");
  if (!(active_set.sw_min > 0.0 && active_set.sw_min < 1.0) ||
      !(active_set.sg_off > 0.0 && active_set.sg_off < 1.0) ||
      active_set.sw_min >= 1.0 - active_set.sg_off)
    throw std::runtime_error("invalid saturation bounds");
  positive(active_set.ng_off, "active_set.ng_off_mol");
  nonnegative(active_set.transfer_switch_tolerance,
              "active_set.transfer_switch_tolerance_mol_m3");
  const auto valid_batch_fraction = [](double value) {
    return std::isfinite(value) && value > 0.0 && value <= 1.0;
  };
  if (!valid_batch_fraction(active_set.batch_individual_free_gas_fraction) ||
      !valid_batch_fraction(active_set.batch_total_free_gas_fraction_per_step) ||
      active_set.batch_total_free_gas_fraction_per_step <
          active_set.batch_individual_free_gas_fraction)
    throw std::runtime_error(
        "active_set batch retirement fractions must lie in (0,1] and the "
        "per-step aggregate limit must not be smaller than the individual limit");
  if (active_set.batch_microbubble_retirement_enabled &&
      !time_step_control.enabled)
    throw std::runtime_error(
        "batch microbubble retirement requires time_step_control.enabled=true");
  positive(scaling.pressure, "scaling.P_ref_Pa");
  positive(scaling.concentration, "scaling.C_ref_mol_m3");
  positive(scaling.volume_rate, "scaling.Q_ref_m3_s");
  positive(scaling.molar_rate, "scaling.Ndot_ref_mol_s");
  positive(scaling.mole_total, "scaling.N_scale_mol");
  positive(scaling.volume_total, "scaling.V_scale_m3");
  positive(acceptance.mass_balance_hard_limit, "acceptance.mass_balance_hard_limit");
  positive(acceptance.volume_balance_hard_limit, "acceptance.volume_balance_hard_limit");
  if (io.checkpoint_interval_steps <= 0 || io.vtk_interval_steps <= 0)
    throw std::runtime_error("output intervals must be positive");
  if(io.state_format!="legacy_vtk"&&io.state_format!="compressed_vtp")
    throw std::runtime_error("io.state_format must be legacy_vtk or compressed_vtp");
  const auto validate_times = [&](const std::vector<double>& times, const char* name) {
    double previous = -1.0;
    for (const double value : times) {
      if (!std::isfinite(value) || value < 0.0 || value <= previous)
        throw std::runtime_error(std::string(name) + " must be finite, nonnegative, and increasing");
      previous = value;
    }
  };
  validate_times(io.vtk_output_times, "io.vtk_output_times_s");
  validate_times(io.checkpoint_output_times, "io.checkpoint_output_times_s");
  if (time.stop_main_gas_volume_fraction != 0.0 &&
      !(time.stop_main_gas_volume_fraction > 0.0 &&
        time.stop_main_gas_volume_fraction < 1.0))
    throw std::runtime_error("time.stop_main_gas_volume_fraction must be zero or lie in (0,1)");
  nonnegative(time.equilibrium_tolerance, "time.equilibrium_tolerance_mol_m3");
  if (time_step_control.enabled) {
    if (!(time_step_control.growth_factor >= 1.0) || !std::isfinite(time_step_control.growth_factor))
      throw std::runtime_error("time_step_control.growth_factor must be finite and >= 1");
    open_unit_interval(time_step_control.shrink_factor, "time_step_control.shrink_factor");
    if (time_step_control.easy_newton_iterations < 0 ||
        time_step_control.normal_newton_iterations < time_step_control.easy_newton_iterations ||
        time_step_control.hard_newton_iterations < time_step_control.normal_newton_iterations ||
        time_step_control.max_consecutive_growth < 0)
      throw std::runtime_error("invalid time_step_control Newton thresholds");
    open_unit_interval(time_step_control.gas_timescale_fraction, "time_step_control.gas_timescale_fraction");
    open_unit_interval(time_step_control.concentration_timescale_fraction,
                       "time_step_control.concentration_timescale_fraction");
    open_unit_interval(time_step_control.concentration_change_fraction,
                       "time_step_control.concentration_change_fraction");
  }
  if (termination.enabled) {
    if (termination.primary != "main_real_gas_saturation" &&
        termination.primary != "residual_free_gas_fraction" &&
        termination.primary != "either")
      throw std::runtime_error("termination.primary must be main_real_gas_saturation, residual_free_gas_fraction, or either");
    if (termination.primary != "residual_free_gas_fraction" &&
        !(termination.target_gas_saturation > 0.0 && termination.target_gas_saturation < 1.0))
      throw std::runtime_error("termination.target_gas_saturation must lie in (0,1) when enabled");
    positive(termination.maximum_physical_time,
             "termination.maximum_physical_time_s");
    nonnegative(termination.maximum_wall_seconds,
                "termination.maximum_wall_seconds");
    if (termination.maximum_output_bytes == 0)
      throw std::runtime_error("termination.maximum_output_bytes must be positive");
    if(termination.residual_free_gas_fraction!=0.0 &&
       !(termination.residual_free_gas_fraction>0.0&&termination.residual_free_gas_fraction<1.0))
      throw std::runtime_error("termination.residual_free_gas_fraction must be zero or lie in (0,1)");
    if ((termination.primary == "residual_free_gas_fraction" ||
         termination.primary == "either") &&
        termination.residual_free_gas_fraction == 0.0)
      throw std::runtime_error("termination.residual_free_gas_fraction must be enabled for the selected primary");
    if (!(time_step_control.gas_timescale_percentile > 0.0 && time_step_control.gas_timescale_percentile <= 1.0))
      throw std::runtime_error("time_step_control.gas_timescale_percentile must lie in (0,1]");
    if (!(time_step_control.gas_timescale_min_mole_fraction >= 0.0 &&
          time_step_control.gas_timescale_min_mole_fraction < 1.0))
      throw std::runtime_error("time_step_control.gas_timescale_min_mole_fraction must lie in [0,1)");
  }
  if (analysis.x_segments != 10 && analysis.x_segments != 20 &&
      analysis.x_segments != 50)
    throw std::runtime_error("analysis.x_segments must be 10, 20, or 50");
  if (boundary.inlet_flag == 0 || boundary.outlet_flag == 0 ||
      boundary.inlet_flag == boundary.outlet_flag)
    throw std::runtime_error("inlet/outlet flags must be distinct and nonzero");
  if (boundary.inlet_type != "pressure" && boundary.inlet_type != "flow")
    throw std::runtime_error("boundary.inlet_type must be pressure or flow");
  if (boundary.outlet_type != "pressure")
    throw std::runtime_error("boundary.outlet_type must be pressure");
  if (boundary.inlet_type == "flow") positive(boundary.inlet_flow_m3_s, "boundary.inlet_flow_m3_s");
}

}  // namespace bubble
