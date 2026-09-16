#pragma once

#include <badiff.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace bubble {

constexpr double gas_constant = 8.31446261815324;
constexpr double pi = 3.141592653589793238462643383279502884;

struct CaseConfig {
  std::string name;
  std::string kind;
  bool scientific_parameters_confirmed = false;
  std::string gas_component;
  std::string liquid_composition;
  std::string notes;
};

struct IoConfig {
  std::string state_format = "legacy_vtk";
  std::filesystem::path pore_file;
  std::filesystem::path throat_file;
  std::filesystem::path connect_audit_file;
  std::filesystem::path output_directory;
  std::filesystem::path input_geometry_audit_file;
  std::filesystem::path virtual_boundary_geometry_audit_file;
  std::filesystem::path physical_parameter_audit_file;
  int checkpoint_interval_steps = 1;
  int vtk_interval_steps = 1;
  std::vector<double> vtk_output_times;
  std::vector<double> checkpoint_output_times;
};

struct GeometryConfig {
  std::string input_basis;
  std::string input_length_unit;
  double voxel_size_um = 0.0;
  std::array<int, 3> image_dimensions_voxels{{0, 0, 0}};
  std::array<double, 3> physical_domain_um{{0.0, 0.0, 0.0}};
  double length_to_m = 1.0;
  double area_to_m2 = 1.0;
  double volume_to_m3 = 1.0;
};

struct VirtualBoundaryConfig {
  bool enabled = false;
  std::string boundary_axis;
  double extension_voxels = 0.0;
  double voxel_size_um = 0.0;
};

struct InputProvenanceConfig {
  std::filesystem::path raw_pore_file, raw_throat_file;
  std::filesystem::path base_pore_file, base_throat_file, base_connect_file;
};

struct PhysicsConfig {
  // Solver pressure variables are offsets from this fixed absolute reference.
  double background_abs_Pa = 0.0;
  double temperature = 0.0;
  double viscosity = 0.0;
  double molecular_diffusivity = 0.0;
  double henry_cp = 0.0;
  double mass_transfer_coefficient = 0.0;
  double mass_transfer_multiplier = 1.0;
  double surface_tension = 0.0;
  double contact_angle = 0.0;
  double molar_mass = 0.0;
  double liquid_density = 0.0;
  double compressibility = 1.0;
  double kr_exponent = 3.0;
  double diffusion_area_exponent = 1.0;
  double tortuosity = 1.0;
  double qin_b = 6.83;
};

struct BoundaryConfig {
  std::string inlet_type = "pressure";
  std::string outlet_type = "pressure";
  double inlet_flow_m3_s = 0.0;
  int inlet_flag = 1;
  int outlet_flag = 2;
  double pl_in = 0.0;
  double pl_out = 0.0;
  double c_in = 0.0;
  double c_backflow = 0.0;
  std::string species_mode;
  bool fixed_concentrations_ignored = false;
};

struct InitialConfig {
  double concentration = 0.0;
};

struct TimeConfig {
  double end = 0.0;
  double initial_dt = 0.0;
  double min_dt = 0.0;
  double max_dt = 0.0;
  double cut_factor = 0.5;
  double growth_factor = 1.2;
  int fast_newton_iterations = 5;
  int max_retries = 16;
  double stop_main_gas_volume_fraction = 0.0;
  double equilibrium_tolerance = 0.0;
};

struct TimeStepControlConfig {
  bool enabled = false;
  double growth_factor = 1.25;
  double shrink_factor = 0.5;
  int easy_newton_iterations = 3;
  int normal_newton_iterations = 7;
  int hard_newton_iterations = 12;
  double gas_timescale_fraction = 0.10;
  double gas_timescale_percentile = 0.50;
  double gas_timescale_min_mole_fraction = 1.0e-5;
  double concentration_timescale_fraction = 0.05;
  bool concentration_change_limit_enabled = false;
  double concentration_change_fraction = 0.05;
  int max_consecutive_growth = 3;
};

struct TerminationConfig {
  bool enabled = false;
  std::string primary = "main_real_gas_saturation";
  double target_gas_saturation = 0.0;
  double maximum_physical_time = 0.0;
  double maximum_wall_seconds = 0.0;
  std::uint64_t maximum_output_bytes = 0;
  double residual_free_gas_fraction = 0.0;
};

struct AnalysisConfig {
  int x_segments = 20;
  std::vector<double> output_pore_volumes;
};

struct NonlinearConfig {
  int max_iterations = 30;
  double residual_tolerance = 1e-8;
  double update_tolerance = 1e-9;
  double armijo_c1 = 1e-4;
  double line_search_reduction = 0.5;
  double line_search_min_lambda = 1e-8;
  double fraction_to_boundary = 0.99;
};

struct LinearConfig {
  std::string backend = "eigen_bicgstab_ilut";
  double relative_tolerance = 1e-10;
  int max_iterations = 2000;
  double ilut_drop_tolerance = 1e-4;
  int ilut_fill_factor = 10;
  int pardiso_threads = 16;
  std::filesystem::path amgx_config_file;
};

struct ActiveSetConfig {
  double sw_min = 1e-6;
  double sg_off = 1e-8;
  double ng_off = 1e-30;
  double transfer_switch_tolerance = 1e-12;
  bool batch_microbubble_retirement_enabled = false;
  double batch_individual_free_gas_fraction = 1e-10;
  double batch_total_free_gas_fraction_per_step = 1e-8;
};

struct ScalingConfig {
  double gas_pressure = 0.0; // zero selects the legacy pressure scale
  double pressure = 1.0;
  double concentration = 1.0;
  double volume_rate = 1.0;
  double molar_rate = 1.0;
  double mole_total = 1.0;
  double volume_total = 1.0;
};

struct AcceptanceConfig {
  double mass_balance_hard_limit = 1e-6;
  double volume_balance_hard_limit = 1e-6;
};

struct Config {
  CaseConfig case_info;
  IoConfig io;
  GeometryConfig geometry;
  VirtualBoundaryConfig virtual_boundary;
  InputProvenanceConfig input_provenance;
  PhysicsConfig physics;
  BoundaryConfig boundary;
  InitialConfig initial;
  TimeConfig time;
  TimeStepControlConfig time_step_control;
  TerminationConfig termination;
  AnalysisConfig analysis;
  NonlinearConfig nonlinear;
  LinearConfig linear;
  ActiveSetConfig active_set;
  ScalingConfig scaling;
  AcceptanceConfig acceptance;
  std::filesystem::path source_path;
  std::string source_text;

  static Config load(const std::filesystem::path& path);
  void validate_before_mesh() const;
};

struct Pore {
  std::size_t id = 0;
  int source_type = 0;
  std::size_t source_id = 0;
  double x = 0.0, y = 0.0, z = 0.0;
  double radius = 0.0, area = 0.0, shape_factor = 0.0, volume = 0.0;
  int boundary_flag = 0;
  double input_sw = 1.0;
};

struct Edge {
  std::size_t id = 0, a = 0, b = 0, source_old_throat_id = 0;
  int side = 0;
  double length = 0.0, input_sw = 1.0, radius = 0.0, area = 0.0;
  double shape_factor = 0.0, source_volume = 0.0;
};

struct ZeroDistanceCluster {
  std::size_t id = 0;
  std::vector<std::size_t> members;
  std::vector<std::size_t> zero_distance_edges;
  std::vector<std::size_t> external_positive_edges;
  int boundary_flag = 0;
};

struct Network {
  std::vector<Pore> pores;
  std::vector<Edge> edges;
  std::vector<std::vector<std::size_t>> adjacency;
  std::vector<std::size_t> cluster_of_pore;
  std::vector<ZeroDistanceCluster> zero_distance_clusters;
  std::vector<std::size_t> zero_distance_edges;
  std::vector<std::size_t> massless_junctions;
  std::size_t production_positive_edges = 0;
  std::size_t throats_with_two_zero_halves = 0;
  std::size_t ignored_boundary_edges = 0;
  std::size_t virtual_boundary_nodes = 0;
  std::size_t boundary_edges = 0;
  std::size_t connected_components = 0;
  std::size_t boundaryless_components = 0;
  std::size_t isolated_pores = 0;
  std::vector<std::size_t> component_id;
  std::vector<int> boundary_connection_class;
  std::size_t main_flow_component_id = 0;
  static Network load(const Config& config);
  void validate(const Config& config) const;
  void build_degenerate_geometry(const Config& config);
};

struct NodeState {
  // Keep the nonlinear state in extended precision.  The full network has
  // high-pressure, nearly-liquid gas controls for which a valid Newton
  // correction is smaller than one binary64 ULP.  The sparse Jacobian and
  // linear solve remain binary64; retaining their updates here prevents the
  // nonlinear residual from stalling solely because the state rounded back
  // to the preceding representable double.
  long double pl = 0.0L;
  long double pg = std::numeric_limits<long double>::quiet_NaN();
  long double sw = 1.0L;
  long double c = 0.0L;
  bool active_gas = false;
  bool permanently_inactive = false;
  bool concentration_derived = false;
};

struct DofMap {
  struct NodeDofs {
    int pl = -1, c = -1, pg = -1, sw = -1;
    int rl = -1, rd = -1, rg = -1, rc = -1;
    bool pressure_gauge = false;
  };
  struct ClusterDofs {
    int pl = -1, c = -1, rl = -1, rd = -1;
    bool pressure_gauge = false;
  };
  std::vector<NodeDofs> node;
  std::vector<ClusterDofs> cluster;
  std::vector<double> unknown_scale, residual_scale;
  int inlet_pressure = -1;
  int inlet_flow_row = -1;
  int size = 0;
};

struct EdgeFlux {
  double q = 0.0, advective = 0.0, diffusive = 0.0;
  bool upwind_a = true;
};

EdgeFlux evaluate_edge_flux(const Edge& edge, double pl_a, double c_a,
                            double sw_a, double pl_b, double c_b, double sw_b,
                            bool upwind_a, bool diffuse,
                            const PhysicsConfig& physics);

struct Assembly {
  Eigen::VectorXd residual;
  Eigen::SparseMatrix<double, Eigen::RowMajor> jacobian;
  std::vector<EdgeFlux> edge_fluxes;
  std::vector<double> transfer;
  std::array<double, 4> residual_norms{{0.0, 0.0, 0.0, 0.0}};
  std::array<std::size_t, 4> residual_node_ids{{0, 0, 0, 0}};
};

struct Balance {
  double pg_min = std::numeric_limits<double>::infinity(), pg_max = -std::numeric_limits<double>::infinity();
  std::size_t pl_min_id=0, pl_max_id=0, pg_min_id=0, pg_max_id=0;
  double component_inflow = 0.0, component_outflow = 0.0;
  double liquid_inflow = 0.0, liquid_outflow = 0.0;
  double transfer_rate = 0.0;
  double total_dissolved = 0.0, total_free = 0.0;
  double boundary_molar_outflow = 0.0, boundary_liquid_outflow = 0.0;
  double boundary_advective_outflow = 0.0, boundary_diffusive_outflow = 0.0;
  double inlet_advective_component_inflow = 0.0;
  double outlet_advective_component_outflow = 0.0;
  double inlet_liquid_inflow = 0.0, outlet_liquid_outflow = 0.0;
  double net_liquid_inflow = 0.0, throughflow_rate = 0.0;
  std::size_t wrong_direction_boundary_edges = 0;
  std::size_t inlet_reverse_flow_edges = 0, outlet_backflow_edges = 0;
  double inlet_reverse_flow = 0.0, outlet_backflow = 0.0;
  double liquid_volume_change = 0.0, cumulative_liquid_volume_change = 0.0;
  double free_main_flow = 0.0, free_boundary_connected = 0.0;
  double free_boundaryless = 0.0;
  double gas_volume_total = 0.0, gas_volume_main_flow = 0.0;
  double gas_volume_boundary_connected = 0.0, gas_volume_boundaryless = 0.0;
  std::size_t active_main_flow = 0, active_boundary_connected = 0;
  std::size_t active_boundaryless = 0;
  double maximum_dissolution_driving_force = 0.0;
  double internal_pressure_min = std::numeric_limits<double>::infinity();
  double internal_pressure_max = -std::numeric_limits<double>::infinity();
  double epsilon_moles_equation = 0.0, epsilon_volume_equation = 0.0;
  double epsilon_moles_adjusted = 0.0, epsilon_volume_adjusted = 0.0;
};

struct GasDisappearanceEventRecord {
  std::size_t node_index = 0;
  std::size_t node_id = 0;
  double estimated_crossing_time = 0.0;
  double accepted_state_time = 0.0;
  double crossing_fraction = 1.0;
  double gas_moles_before = 0.0;
  double gas_moles_after_before_retirement = 0.0;
  double gas_moles_transferred_to_dissolved = 0.0;
  double gas_volume_before = 0.0;
  double gas_volume_after_before_retirement = 0.0;
  double individual_free_gas_fraction = 0.0;
  double batch_aggregate_free_gas_fraction = 0.0;
  bool batched = false;
};

struct StepReport {
  bool accepted = false;
  std::string failure;
  int newton_iterations = 0, linear_iterations = 0;
  double linear_residual = 0.0;
  int amgx_status = -1;
  int amgx_upload_all_calls = 0;
  int amgx_replace_coefficients_calls = 0;
  int amgx_solver_setup_calls = 0;
  double assembly_seconds = 0.0, csr_seconds = 0.0;
  double linear_upload_seconds = 0.0, linear_setup_seconds = 0.0;
  double linear_solve_seconds = 0.0, linear_download_seconds = 0.0;
  double linear_symbolic_analysis_seconds = 0.0;
  double linear_numeric_factorization_seconds = 0.0;
  int pardiso_symbolic_analysis_calls = 0;
  int pardiso_numeric_factorization_calls = 0;
  std::int64_t pardiso_factor_nonzeros = 0;
  std::int64_t pardiso_peak_memory_kib = 0;
  double residual_norm = std::numeric_limits<double>::infinity();
  double update_norm = std::numeric_limits<double>::infinity();
  double state_norm = std::numeric_limits<double>::infinity();
  double relative_update_norm = std::numeric_limits<double>::infinity();
  double update_tolerance_limit = std::numeric_limits<double>::infinity();
  double update_criterion_ratio = std::numeric_limits<double>::infinity();
  bool residual_criterion_passed = false;
  bool update_criterion_passed = false;
  double dt = 0.0;
  double requested_dt = 0.0;
  int timestep_retries = 0;
  std::vector<double> rejected_dts;
  std::vector<std::string> retry_failures;
  std::size_t active_before = 0, active_after = 0;
  double threshold_moles_transferred = 0.0;
  double active_set_mole_ledger = 0.0;
  double active_set_liquid_volume_ledger = 0.0;
  double min_free_gas_moles = 0.0;
  double min_Sg = 0.0;
  std::string dt_change_reason = "initial";
  double next_dt = 0.0;
  double gas_timescale_cap = -1.0;
  double concentration_timescale_cap = -1.0;
  bool gas_disappearance_event = false;
  std::size_t event_node_id = 0;
  double event_time = 0.0;
  double event_dt = 0.0;
  double event_gas_moles_before = 0.0;
  double event_gas_moles_after = 0.0;
  bool event_timestep_localized = false;
  std::size_t retired_bubble_count = 0;
  std::size_t batched_retired_bubble_count = 0;
  double batched_retired_moles = 0.0;
  double batched_retired_gas_volume = 0.0;
  double batch_aggregate_free_gas_fraction = 0.0;
  std::vector<GasDisappearanceEventRecord> gas_disappearance_event_records;
  std::size_t concentration_limit_node_id = 0;
  double concentration_limit_C = 0.0;
  double concentration_limit_dC_dt = 0.0;
  double concentration_reference = 0.0;
  double concentration_tau = 0.0;
  Balance balance;
  Assembly final_assembly;
};

struct AmgxSolveMetrics {
  int iterations = 0;
  int amgx_status = -1;
  int matrix_upload_all_calls = 0;
  int matrix_replace_coefficients_calls = 0;
  int solver_setup_calls = 0;
  double csr_seconds = 0.0, upload_seconds = 0.0, setup_seconds = 0.0;
  double solve_seconds = 0.0, download_seconds = 0.0;
  double cpu_relative_residual = std::numeric_limits<double>::infinity();
};

class AmgxLinearSolver;
class PardisoLinearSolver;

inline double scalar_primal(double value) { return value; }
inline long double scalar_primal(long double value) { return value; }
inline double scalar_primal(const fadbad::B<double>& value) {
  fadbad::B<double> copy = value;
  return copy.x();
}

template <class Scalar>
Scalar one_minus_exp_negative(const Scalar& value) {
  using std::abs;
  using std::exp;
  // expm1 is not supplied by FADBAD++ 2.1.  Directly evaluating
  // 1-exp(-x) loses enough digits at x~1e-6 to create an O(1 Pa) error after
  // division by nanometre radii.  This Horner series is stable for both
  // double and reverse AD and has a remainder below double roundoff at 0.1.
  if (abs(scalar_primal(value)) < 0.1) {
    return value *
           (1.0 +
            value *
                (-1.0 / 2.0 +
                 value *
                     (1.0 / 6.0 +
                      value *
                          (-1.0 / 24.0 +
                           value *
                               (1.0 / 120.0 +
                                value *
                                    (-1.0 / 720.0 +
                                     value *
                                         (1.0 / 5040.0 +
                                          value *
                                              (-1.0 / 40320.0 +
                                               value *
                                                   (1.0 / 362880.0 +
                                                    value *
                                                        (-1.0 /
                                                         3628800.0))))))))));
  }
  return 1.0 - exp(-value);
}

template <class Scalar>
Scalar capillary_pressure(const Scalar& sw, double radius,
                          const PhysicsConfig& p) {
  using std::cos;
  return 2.0 * p.surface_tension * cos(p.contact_angle) /
         (radius * one_minus_exp_negative(p.qin_b * sw));
}

template <class Scalar>
Scalar interface_area(const Scalar& sw, double volume) {
  using std::pow;
  return pow(36.0 * pi, 1.0 / 3.0) * pow(volume * (1.0 - sw), 2.0 / 3.0);
}

template <class Scalar>
Scalar free_gas_moles(const Scalar& pg, const Scalar& sw, double volume,
                      const PhysicsConfig& p) {
  return (p.background_abs_Pa + pg) * volume * (1.0 - sw) /
         (p.compressibility * gas_constant * p.temperature);
}

class Solver {
 public:
  Solver(Network network, Config config);
  ~Solver();
  Solver(Solver&&) noexcept;
  Solver& operator=(Solver&&) noexcept;
  Solver(const Solver&) = delete;
  Solver& operator=(const Solver&) = delete;
  const Network& network() const { return network_; }
  const Config& config() const { return config_; }
  const std::vector<NodeState>& state() const { return state_; }
  std::vector<NodeState>& mutable_state_for_test() { return state_; }
  double time() const { return time_; }
  int step_number() const { return step_number_; }
  double initial_regularization_moles() const {
    return initial_regularization_moles_;
  }
  double initial_regularization_volume() const {
    return initial_regularization_volume_;
  }
  double initial_threshold_moles_transferred() const {
    return initial_threshold_moles_transferred_;
  }
  double initial_main_gas_volume() const { return initial_main_gas_volume_; }
  double initial_total_moles() const { return initial_total_moles_; }
  double initial_liquid_volume() const { return initial_liquid_volume_; }
  const std::vector<long double>& initial_free_gas_by_node() const { return initial_free_gas_by_node_; }
  const std::vector<long double>& initial_sw_by_node() const { return initial_sw_by_node_; }
  const std::vector<long double>& cumulative_transfer_by_node() const { return cumulative_transfer_by_node_; }
  const std::vector<std::size_t>& regularized_nodes() const {
    return regularized_nodes_;
  }

  void initialize();
  DofMap make_dof_map(const std::vector<NodeState>& state) const;
  Assembly assemble(const std::vector<NodeState>& trial,
                    const std::vector<NodeState>& old, double dt,
                    const DofMap& map, const std::vector<bool>& upwind_a,
                    const std::vector<bool>& transfer_active,
                    bool build_jacobian) const;
  StepReport attempt_step(double dt);
  StepReport advance_with_retry(double requested_dt);
  Balance compute_balance(const std::vector<NodeState>& old,
                          const std::vector<NodeState>& current,
                          const Assembly& assembly, double dt,
                          double mole_ledger, double volume_ledger) const;
  void restore(std::vector<NodeState> state, double time, int step_number,
               double cumulative_mole_ledger, double cumulative_volume_ledger,
               double cumulative_net_liquid_inflow = 0.0,
               double cumulative_boundary_component_outflow = 0.0,
               double cumulative_outlet_advective_component_outflow = 0.0,
               double cumulative_liquid_inflow = 0.0,
               double cumulative_liquid_outflow = 0.0,
               std::vector<long double> cumulative_transfer_by_node = {});
  double cumulative_mole_ledger() const { return cumulative_mole_ledger_; }
  double cumulative_volume_ledger() const { return cumulative_volume_ledger_; }
  double cumulative_net_liquid_inflow() const {
    return cumulative_net_liquid_inflow_;
  }
  double cumulative_boundary_component_outflow() const {
    return cumulative_boundary_component_outflow_;
  }
  double cumulative_outlet_advective_component_outflow() const {
    return cumulative_outlet_advective_component_outflow_;
  }
  double cumulative_component_inflow() const { return cumulative_component_inflow_; }
  double cumulative_component_outflow() const { return cumulative_component_outflow_; }
  double cumulative_liquid_inflow() const { return cumulative_liquid_inflow_; }
  double cumulative_liquid_outflow() const { return cumulative_liquid_outflow_; }
  double cumulative_interphase_transfer() const { return cumulative_interphase_transfer_; }
  double inlet_pressure() const { return inlet_pressure_; }
  void restore_inlet_pressure(double value) { inlet_pressure_ = value; }
  void set_adaptive_metadata(std::string reason,double min_ng,double min_sg){
    last_dt_change_reason_=std::move(reason);last_min_free_gas_moles_=min_ng;last_min_sg_=min_sg;
  }
  const std::string& last_dt_change_reason() const{return last_dt_change_reason_;}
  double last_min_free_gas_moles()const{return last_min_free_gas_moles_;}
  double last_min_sg()const{return last_min_sg_;}
  std::uint64_t gas_disappearance_events() const { return gas_disappearance_events_; }
  std::size_t last_gas_disappearance_node_id() const { return last_gas_disappearance_node_id_; }
  double last_gas_disappearance_time() const { return last_gas_disappearance_time_; }
  void restore_gas_disappearance_metadata(std::uint64_t count,
      std::size_t node_id, double event_time) {
    gas_disappearance_events_=count;last_gas_disappearance_node_id_=node_id;
    last_gas_disappearance_time_=event_time;
  }
  void restore_gross_fluxes(double inflow, double outflow) {
    cumulative_component_inflow_ = inflow; cumulative_component_outflow_ = outflow;
  }

 private:
  Network network_;
  Config config_;
  std::vector<NodeState> state_;
  double time_ = 0.0;
  int step_number_ = 0;
  double cumulative_mole_ledger_ = 0.0, cumulative_volume_ledger_ = 0.0;
  double cumulative_net_liquid_inflow_ = 0.0;
  double cumulative_boundary_component_outflow_ = 0.0;
  double cumulative_outlet_advective_component_outflow_ = 0.0;
  double cumulative_component_inflow_ = 0.0, cumulative_component_outflow_ = 0.0;
  double cumulative_liquid_inflow_ = 0.0, cumulative_liquid_outflow_ = 0.0;
  double cumulative_interphase_transfer_ = 0.0;
  double initial_regularization_moles_ = 0.0,
         initial_regularization_volume_ = 0.0;
  double initial_threshold_moles_transferred_ = 0.0;
  double initial_total_moles_ = 0.0, internal_volume_ = 0.0;
  double initial_liquid_volume_ = 0.0;
  double initial_main_gas_volume_ = 0.0;
  std::vector<long double> initial_free_gas_by_node_;
  std::vector<long double> initial_sw_by_node_;
  std::vector<long double> cumulative_transfer_by_node_;
  std::vector<std::size_t> regularized_nodes_;
  std::unique_ptr<AmgxLinearSolver> amgx_solver_;
  std::unique_ptr<PardisoLinearSolver> pardiso_solver_;
  bool linear_system_dumped_ = false;
  std::string last_dt_change_reason_ = "initial";
  double last_min_free_gas_moles_ = 0.0, last_min_sg_ = 0.0;
  double inlet_pressure_ = 0.0;
  std::uint64_t gas_disappearance_events_ = 0;
  std::size_t last_gas_disappearance_node_id_ = 0;
  double last_gas_disappearance_time_ = -1.0;
};

struct CheckpointData {
  double cumulative_component_inflow = 0.0, cumulative_component_outflow = 0.0;
  double background_abs_Pa = 0.0;
  std::vector<NodeState> state;
  double time = 0.0;
  int step_number = 0;
  double cumulative_mole_ledger = 0.0, cumulative_volume_ledger = 0.0;
  double cumulative_net_liquid_inflow = 0.0;
  double cumulative_boundary_component_outflow = 0.0;
  double cumulative_outlet_advective_component_outflow = 0.0;
  double cumulative_liquid_inflow = 0.0, cumulative_liquid_outflow = 0.0;
  std::vector<long double> cumulative_transfer_by_node;
  std::string last_dt_change_reason = "initial";
  double last_min_free_gas_moles = 0.0, last_min_sg = 0.0;
  double mass_transfer_multiplier = 1.0;
  bool batch_microbubble_retirement_enabled = false;
  double batch_individual_free_gas_fraction = 1e-10;
  double batch_total_free_gas_fraction_per_step = 1e-8;
  double inlet_pressure = 0.0;
  std::uint64_t gas_disappearance_events = 0;
  std::size_t last_gas_disappearance_node_id = 0;
  double last_gas_disappearance_time = -1.0;
  // Driver state is part of a mathematically reproducible restart.  Restoring
  // only the primary variables changes the subsequent adaptive timestep
  // sequence whenever growth or retry logic is active.
  double next_dt = 0.0;
  int consecutive_fast_steps = 0;
  std::string continuation_config_sha256;
  std::string pore_sha256, throat_sha256, connect_sha256;
};

struct CheckpointProvenance {
  std::string continuation_config_sha256;
  std::string pore_sha256, throat_sha256, connect_sha256;
};

void write_vtk(const std::filesystem::path& path, const Solver& solver,
               const Assembly* assembly);
void write_vtp(const std::filesystem::path& path, const Solver& solver,
               const Assembly* assembly);
void write_degenerate_geometry_audit(const std::filesystem::path& path,
                                     const Solver& solver);
void write_component_boundary_and_makeup_audit(
    const std::filesystem::path& path, const Solver& solver);
void write_physical_parameter_audit(const std::filesystem::path& path,
                                    const Solver& solver);
void audit_vtk(const std::filesystem::path& path, std::size_t expected_points,
               std::size_t expected_edges);
void write_checkpoint(const std::filesystem::path& path, const Solver& solver,
                      double next_dt, int consecutive_fast_steps,
                      const CheckpointProvenance& provenance);
CheckpointData read_checkpoint(const std::filesystem::path& path,
                               std::size_t expected_nodes);
CheckpointProvenance make_checkpoint_provenance(const Config& config);
std::string sha256_file(const std::filesystem::path& path);
std::string utc_timestamp();

}  // namespace bubble
