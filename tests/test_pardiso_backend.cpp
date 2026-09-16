#include "bubble/pardiso_backend.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
  Eigen::SparseMatrix<double, Eigen::RowMajor> matrix(6, 6);
  std::vector<Eigen::Triplet<double>> entries;
  for (int i = 0; i < 6; ++i) {
    entries.emplace_back(i, i, 5.0 + i);
    if (i > 0) {
      entries.emplace_back(i, i - 1, -1.3);
      entries.emplace_back(i - 1, i, -0.6);
    }
  }
  matrix.setFromTriplets(entries.begin(), entries.end());
  matrix.makeCompressed();
  Eigen::VectorXd expected(6);
  expected << 1.0, -2.0, 3.0, -4.0, 5.0, -6.0;

  bubble::PardisoLinearSolver solver(2, 1e-10);
  bubble::PardisoSolveMetrics first;
  const Eigen::VectorXd rhs = matrix * expected;
  const Eigen::VectorXd first_solution = solver.solve(matrix, rhs, first);
  if (first.symbolic_analysis_calls != 1 || first.numeric_factorization_calls != 1 ||
      first.structure_reused || first.cpu_relative_residual > 1e-12 ||
      (first_solution - expected).norm() > 1e-11)
    throw std::runtime_error("PARDISO first solve or analysis accounting failed");

  for (int row = 0; row < matrix.outerSize(); ++row)
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(matrix, row);
         it; ++it)
      if (it.row() == it.col()) it.valueRef() += 0.25;
  const Eigen::VectorXd second_rhs = matrix * expected;
  bubble::PardisoSolveMetrics second;
  const Eigen::VectorXd second_solution = solver.solve(matrix, second_rhs, second);
  if (second.symbolic_analysis_calls != 0 ||
      second.numeric_factorization_calls != 1 || !second.structure_reused ||
      second.cpu_relative_residual > 1e-12 ||
      (second_solution - expected).norm() > 1e-11)
    throw std::runtime_error("PARDISO fixed-structure reuse failed");

  matrix.coeffRef(0, 5) = 0.1;
  matrix.makeCompressed();
  const Eigen::VectorXd changed_rhs = matrix * expected;
  bubble::PardisoSolveMetrics changed;
  const Eigen::VectorXd changed_solution = solver.solve(matrix, changed_rhs, changed);
  if (changed.symbolic_analysis_calls != 1 || changed.structure_reused ||
      changed.cpu_relative_residual > 1e-12 ||
      (changed_solution - expected).norm() > 1e-11)
    throw std::runtime_error("PARDISO changed-structure reanalysis failed");

  std::cout << "PASS PARDISO nonsymmetric solve, reuse, and graph-change fallback residual="
            << changed.cpu_relative_residual << '\n';
}
