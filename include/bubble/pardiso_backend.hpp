#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cstdint>
#include <limits>
#include <memory>
#include <string>

namespace bubble {

struct PardisoSolveMetrics {
  int iterative_refinement_steps = 0;
  int symbolic_analysis_calls = 0;
  int numeric_factorization_calls = 0;
  bool structure_reused = false;
  double csr_seconds = 0.0;
  double symbolic_analysis_seconds = 0.0;
  double numeric_factorization_seconds = 0.0;
  double solve_seconds = 0.0;
  double cpu_relative_residual = std::numeric_limits<double>::infinity();
  std::int64_t factor_nonzeros = 0;
  std::int64_t peak_memory_kib = 0;
};

struct PardisoRuntimeInfo {
  bool compiled = false;
  int threads = 0;
  std::string library_path;
};

class PardisoLinearSolver {
 public:
  PardisoLinearSolver(int threads, double relative_tolerance);
  ~PardisoLinearSolver();
  PardisoLinearSolver(const PardisoLinearSolver&) = delete;
  PardisoLinearSolver& operator=(const PardisoLinearSolver&) = delete;

  Eigen::VectorXd solve(
      const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix,
      const Eigen::VectorXd& rhs, PardisoSolveMetrics& metrics);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

PardisoRuntimeInfo pardiso_runtime_info(int configured_threads = 0);

}  // namespace bubble
