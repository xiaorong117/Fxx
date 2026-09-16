#include "bubble/pardiso_backend.hpp"

#include <unsupported/Eigen/SparseExtra>
#include <nlohmann/json.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sys/resource.h>

int main(int argc, char** argv) {
  if (argc != 5)
    throw std::runtime_error("usage: matrix rhs threads result.json");
  const int threads = std::stoi(argv[3]);
  nlohmann::json result{{"backend", "pardiso"}, {"threads", threads}};
  try {
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrix;
    Eigen::VectorXd rhs;
    if (!Eigen::loadMarket(matrix, argv[1]) ||
        !Eigen::loadMarketVector(rhs, argv[2]))
      throw std::runtime_error("MatrixMarket load failed");
    bubble::PardisoLinearSolver solver(threads, 1e-10);
    bubble::PardisoSolveMetrics first, reused;
    const Eigen::VectorXd x1 = solver.solve(matrix, rhs, first);
    const Eigen::VectorXd x2 = solver.solve(matrix, rhs, reused);
    const double state_difference = (x1 - x2).lpNorm<Eigen::Infinity>();
    struct rusage usage {};
    getrusage(RUSAGE_SELF, &usage);
    result.update({
        {"status", "SUCCESS"},
        {"first", {{"cpu_recomputed_relative_residual", first.cpu_relative_residual},
                   {"csr_seconds", first.csr_seconds},
                   {"symbolic_analysis_seconds", first.symbolic_analysis_seconds},
                   {"numeric_factorization_seconds", first.numeric_factorization_seconds},
                   {"solve_seconds", first.solve_seconds},
                   {"symbolic_analysis_calls", first.symbolic_analysis_calls},
                   {"factor_nonzeros", first.factor_nonzeros},
                   {"peak_memory_kib", first.peak_memory_kib}}},
        {"reused", {{"cpu_recomputed_relative_residual", reused.cpu_relative_residual},
                    {"csr_seconds", reused.csr_seconds},
                    {"symbolic_analysis_seconds", reused.symbolic_analysis_seconds},
                    {"numeric_factorization_seconds", reused.numeric_factorization_seconds},
                    {"solve_seconds", reused.solve_seconds},
                    {"symbolic_analysis_calls", reused.symbolic_analysis_calls},
                    {"structure_reused", reused.structure_reused}}},
        {"repeat_solution_max_absolute_difference", state_difference},
        {"maximum_resident_set_kib", usage.ru_maxrss}});
  } catch (const std::exception& error) {
    result.update({{"status", "FAILED"}, {"error", error.what()}});
  }
  std::ofstream output(argv[4]);
  output << std::setw(2) << result << '\n';
  std::cout << result.dump(2) << '\n';
  return result["status"] == "SUCCESS" ? 0 : 2;
}
