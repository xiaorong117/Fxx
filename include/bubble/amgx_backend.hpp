#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <filesystem>
#include <memory>
#include <string>

#include "bubble/model.hpp"

namespace bubble {

struct AmgxRuntimeInfo {
  bool compiled = false;
  int api_major = 0, api_minor = 0;
  std::string mode = "dDDI";
  std::string library_path;
};

class AmgxLinearSolver {
 public:
  explicit AmgxLinearSolver(const std::filesystem::path& config_file);
  ~AmgxLinearSolver();
  AmgxLinearSolver(const AmgxLinearSolver&) = delete;
  AmgxLinearSolver& operator=(const AmgxLinearSolver&) = delete;

  Eigen::VectorXd solve(const Eigen::SparseMatrix<double, Eigen::RowMajor>& matrix,
                        const Eigen::VectorXd& rhs, AmgxSolveMetrics& metrics);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

AmgxRuntimeInfo amgx_runtime_info();

}  // namespace bubble
