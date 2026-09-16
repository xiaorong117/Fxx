#include "bubble/model.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using json = nlohmann::json;

void restore(bubble::Solver& solver, const bubble::CheckpointData& checkpoint) {
  solver.restore(checkpoint.state, checkpoint.time, checkpoint.step_number,
      checkpoint.cumulative_mole_ledger, checkpoint.cumulative_volume_ledger,
      checkpoint.cumulative_net_liquid_inflow,
      checkpoint.cumulative_boundary_component_outflow,
      checkpoint.cumulative_outlet_advective_component_outflow,
      checkpoint.cumulative_liquid_inflow, checkpoint.cumulative_liquid_outflow,
      checkpoint.cumulative_transfer_by_node);
  solver.restore_gross_fluxes(checkpoint.cumulative_component_inflow,
                              checkpoint.cumulative_component_outflow);
}

double main_mean_sg(const bubble::Solver& solver) {
  double volume = 0.0, gas = 0.0;
  for (std::size_t i = 0; i < solver.state().size(); ++i) {
    const auto& pore = solver.network().pores[i];
    if (pore.source_type == 2 || pore.volume == 0.0 ||
        solver.network().component_id[i] != solver.network().main_flow_component_id)
      continue;
    volume += pore.volume;
    gas += pore.volume * static_cast<double>(1.0L - solver.state()[i].sw);
  }
  return gas / volume;
}

int main(int argc, char** argv) {
  if (argc != 5 && argc != 6)
    throw std::runtime_error(
        "usage: config checkpoint output.json pardiso_threads [dt]");
  auto eigen_config = bubble::Config::load(argv[1]);
  auto pardiso_config = eigen_config;
  pardiso_config.linear.backend = "pardiso";
  pardiso_config.linear.pardiso_threads = std::stoi(argv[4]);
  pardiso_config.linear.relative_tolerance = 1e-10;
  const auto checkpoint = bubble::read_checkpoint(
      argv[2], bubble::Network::load(eigen_config).pores.size());
  auto network = bubble::Network::load(eigen_config);
  bubble::Solver eigen(network, eigen_config),
      pardiso(std::move(network), pardiso_config);
  eigen.initialize();
  pardiso.initialize();
  restore(eigen, checkpoint);
  restore(pardiso, checkpoint);
  const double dt = argc == 6 ? std::stod(argv[5]) : 0.01;
  const auto eigen_report = eigen.attempt_step(dt);
  const auto pardiso_report = pardiso.attempt_step(dt);
  if (!eigen_report.accepted || !pardiso_report.accepted)
    throw std::runtime_error("backend step failed: Eigen=" + eigen_report.failure +
                             " PARDISO=" + pardiso_report.failure);

  std::array<double, 4> variable_max{{0.0, 0.0, 0.0, 0.0}};
  double maximum_scaled = 0.0;
  std::size_t worst_node = 0;
  for (std::size_t i = 0; i < eigen.state().size(); ++i) {
    const auto& a = eigen.state()[i];
    const auto& b = pardiso.state()[i];
    const std::array<double, 4> differences{{
        std::abs(static_cast<double>(a.pl - b.pl)) / eigen_config.scaling.pressure,
        a.active_gas ? std::abs(static_cast<double>(a.pg - b.pg)) /
                           eigen_config.scaling.gas_pressure : 0.0,
        std::abs(static_cast<double>(a.sw - b.sw)),
        std::abs(static_cast<double>(a.c - b.c)) /
            eigen_config.scaling.concentration}};
    for (std::size_t j = 0; j < differences.size(); ++j)
      variable_max[j] = std::max(variable_max[j], differences[j]);
    const double value = *std::max_element(differences.begin(), differences.end());
    if (value > maximum_scaled) {
      maximum_scaled = value;
      worst_node = eigen.network().pores[i].id;
    }
  }
  const auto relative_difference = [](double a, double b) {
    return std::abs(a - b) / std::max({std::abs(a), std::abs(b), 1e-300});
  };
  const double eigen_sg = main_mean_sg(eigen);
  const double pardiso_sg = main_mean_sg(pardiso);
  json output{
      {"format", "bubble_pardiso_full_checkpoint_comparison_v1"},
      {"checkpoint", std::filesystem::absolute(argv[2]).string()},
      {"dt_s", dt},
      {"pardiso_threads", pardiso_config.linear.pardiso_threads},
      {"dof_count", eigen_report.final_assembly.jacobian.rows()},
      {"jacobian_nonzeros", eigen_report.final_assembly.jacobian.nonZeros()},
      {"max_scaled_state_difference", maximum_scaled},
      {"maximum_scaled_variable_difference", {{"Pl", variable_max[0]},
                                                {"Pg", variable_max[1]},
                                                {"Sw", variable_max[2]},
                                                {"C", variable_max[3]}}},
      {"worst_state_node_id", worst_node},
      {"main_mean_Sg", {{"eigen", eigen_sg}, {"pardiso", pardiso_sg},
                         {"relative_difference", relative_difference(eigen_sg, pardiso_sg)}}},
      {"free_gas_mol", {{"eigen", eigen_report.balance.total_free},
                         {"pardiso", pardiso_report.balance.total_free},
                         {"relative_difference", relative_difference(
                             eigen_report.balance.total_free,
                             pardiso_report.balance.total_free)}}},
      {"dissolved_gas_mol", {{"eigen", eigen_report.balance.total_dissolved},
                              {"pardiso", pardiso_report.balance.total_dissolved},
                              {"relative_difference", relative_difference(
                                  eigen_report.balance.total_dissolved,
                                  pardiso_report.balance.total_dissolved)}}},
      {"eigen", {{"nonlinear_residual", eigen_report.residual_norm},
                  {"linear_residual", eigen_report.linear_residual},
                  {"newton_iterations", eigen_report.newton_iterations},
                  {"linear_iterations", eigen_report.linear_iterations},
                  {"relative_update_norm", eigen_report.relative_update_norm},
                  {"assembly_seconds", eigen_report.assembly_seconds},
                  {"setup_seconds", eigen_report.linear_setup_seconds},
                  {"solve_seconds", eigen_report.linear_solve_seconds},
                  {"epsilon_N", eigen_report.balance.epsilon_moles_adjusted},
                  {"epsilon_V", eigen_report.balance.epsilon_volume_adjusted}}},
      {"pardiso", {{"nonlinear_residual", pardiso_report.residual_norm},
                    {"linear_residual", pardiso_report.linear_residual},
                    {"newton_iterations", pardiso_report.newton_iterations},
                    {"linear_solve_calls", pardiso_report.linear_iterations},
                    {"relative_update_norm", pardiso_report.relative_update_norm},
                    {"assembly_seconds", pardiso_report.assembly_seconds},
                    {"csr_seconds", pardiso_report.csr_seconds},
                    {"symbolic_analysis_seconds", pardiso_report.linear_symbolic_analysis_seconds},
                    {"numeric_factorization_seconds", pardiso_report.linear_numeric_factorization_seconds},
                    {"setup_seconds", pardiso_report.linear_setup_seconds},
                    {"solve_seconds", pardiso_report.linear_solve_seconds},
                    {"symbolic_analysis_calls", pardiso_report.pardiso_symbolic_analysis_calls},
                    {"numeric_factorization_calls", pardiso_report.pardiso_numeric_factorization_calls},
                    {"factor_nonzeros", pardiso_report.pardiso_factor_nonzeros},
                    {"peak_memory_kib", pardiso_report.pardiso_peak_memory_kib},
                    {"epsilon_N", pardiso_report.balance.epsilon_moles_adjusted},
                    {"epsilon_V", pardiso_report.balance.epsilon_volume_adjusted}}}};
  output["acceptance"] = {
      {"state_difference_pass", maximum_scaled <= 1e-7},
      {"mean_Sg_pass", relative_difference(eigen_sg, pardiso_sg) <= 1e-8},
      {"free_moles_pass", relative_difference(eigen_report.balance.total_free,
                                               pardiso_report.balance.total_free) <= 1e-8},
      {"dissolved_moles_pass", relative_difference(
                                    eigen_report.balance.total_dissolved,
                                    pardiso_report.balance.total_dissolved) <= 1e-8},
      {"both_nonlinear_gates_pass",
       eigen_report.residual_criterion_passed && eigen_report.update_criterion_passed &&
           pardiso_report.residual_criterion_passed && pardiso_report.update_criterion_passed},
      {"conservation_pass",
       eigen_report.balance.epsilon_moles_adjusted <=
               eigen_config.acceptance.mass_balance_hard_limit &&
           pardiso_report.balance.epsilon_moles_adjusted <=
               eigen_config.acceptance.mass_balance_hard_limit &&
           eigen_report.balance.epsilon_volume_adjusted <=
               eigen_config.acceptance.volume_balance_hard_limit &&
           pardiso_report.balance.epsilon_volume_adjusted <=
               eigen_config.acceptance.volume_balance_hard_limit}};
  std::ofstream file(argv[3]);
  file << std::setw(2) << output << '\n';
  std::cout << output.dump(2) << '\n';
  return std::all_of(output["acceptance"].begin(), output["acceptance"].end(),
                     [](const auto& item) { return item.template get<bool>(); }) ? 0 : 2;
}
