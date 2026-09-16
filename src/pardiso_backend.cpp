#include "bubble/pardiso_backend.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#if BUBBLE_PARDISO_ENABLED
#include <dlfcn.h>
#endif

namespace bubble {
namespace {

#if BUBBLE_PARDISO_ENABLED
using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point begin) {
  return std::chrono::duration<double>(Clock::now() - begin).count();
}
#endif

#if BUBBLE_PARDISO_ENABLED
using MklInt = int;

using PardisoInitFunction = void (*)(void*, const MklInt*, MklInt*);
using PardisoFunction = void (*)(
    void*, const MklInt*, const MklInt*, const MklInt*, const MklInt*,
    const MklInt*, const void*, const MklInt*, const MklInt*, MklInt*,
    const MklInt*, MklInt*, const MklInt*, void*, void*, MklInt*);

template <class Function>
Function load_symbol(void* library, const char* name) {
  dlerror();
  void* symbol = dlsym(library, name);
  if (const char* error = dlerror())
    throw std::runtime_error(std::string("oneMKL symbol ") + name +
                             " is unavailable: " + error);
  return reinterpret_cast<Function>(symbol);
}

const char* error_name(MklInt error) {
  switch (error) {
    case 0: return "PARDISO_SUCCESS";
    case -1: return "PARDISO_INPUT_INCONSISTENT";
    case -2: return "PARDISO_INSUFFICIENT_MEMORY";
    case -3: return "PARDISO_REORDERING_ERROR";
    case -4: return "PARDISO_ZERO_PIVOT_OR_NUMERICAL_FACTORIZATION_ERROR";
    case -5: return "PARDISO_INTERNAL_ERROR";
    case -6: return "PARDISO_PREORDERING_ERROR";
    case -7: return "PARDISO_DIAGONAL_MATRIX_ERROR";
    case -8: return "PARDISO_32BIT_INTEGER_OVERFLOW";
    case -10: return "PARDISO_NO_LICENSE_FILE";
    case -11: return "PARDISO_LICENSE_EXPIRED";
    case -12: return "PARDISO_LICENSE_USER_HOST_MISMATCH";
    default: return "PARDISO_UNKNOWN_ERROR";
  }
}

[[noreturn]] void throw_error(MklInt error, MklInt phase) {
  std::ostringstream message;
  message << "PARDISO phase " << phase << " failed with " << error_name(error)
          << " (" << error << ')';
  throw std::runtime_error(message.str());
}
#endif

}  // namespace

struct PardisoLinearSolver::Impl {
#if BUBBLE_PARDISO_ENABLED
  void* library = nullptr;
  PardisoInitFunction pardiso_init = nullptr;
  PardisoFunction pardiso_call = nullptr;
  std::array<void*, 64> pt{};
  std::array<MklInt, 64> iparm{};
  std::vector<MklInt> row_ptr;
  std::vector<MklInt> columns;
  MklInt n = 0;
  MklInt nnz = 0;
  MklInt mtype = 11;
  MklInt threads = 1;
  double relative_tolerance = 1e-10;
  bool analyzed = false;
#endif
};

PardisoLinearSolver::PardisoLinearSolver(int threads,
                                         double relative_tolerance)
    : impl_(std::make_unique<Impl>()) {
#if BUBBLE_PARDISO_ENABLED
  if (threads <= 0)
    throw std::runtime_error("linear.pardiso_threads must be positive");
  if (!(relative_tolerance > 0.0) || !std::isfinite(relative_tolerance))
    throw std::runtime_error("PARDISO relative tolerance must be finite and positive");
  impl_->threads = threads;
  impl_->relative_tolerance = relative_tolerance;
  const std::string thread_text = std::to_string(threads);
  setenv("MKL_NUM_THREADS", thread_text.c_str(), 1);
  setenv("OMP_NUM_THREADS", thread_text.c_str(), 1);
  setenv("MKL_DYNAMIC", "FALSE", 1);
  setenv("OMP_DYNAMIC", "FALSE", 1);
  impl_->library = dlopen(BUBBLE_PARDISO_LIBRARY_PATH, RTLD_NOW | RTLD_LOCAL);
  if (!impl_->library)
    throw std::runtime_error(std::string("cannot load oneMKL PARDISO runtime: ") +
                             dlerror());
  try {
    impl_->pardiso_init =
        load_symbol<PardisoInitFunction>(impl_->library, "pardisoinit");
    impl_->pardiso_call =
        load_symbol<PardisoFunction>(impl_->library, "pardiso");
  } catch (...) {
    dlclose(impl_->library);
    impl_->library = nullptr;
    throw;
  }
  impl_->pardiso_init(impl_->pt.data(), &impl_->mtype, impl_->iparm.data());
  // C-style, zero-based CSR input. This is the C iparm[34] entry, documented
  // as Fortran iparm(35).
  impl_->iparm[34] = 1;
  // Permit additional iterative refinement for ill-conditioned coupled
  // Jacobians. PARDISO reports the number actually used in iparm[6].
  impl_->iparm[7] = 8;
  // In-core factorization. The frozen audit establishes that memory is ample.
  impl_->iparm[59] = 0;
#else
  (void)threads;
  (void)relative_tolerance;
  throw std::runtime_error(
      "linear.backend=pardiso requested, but this binary was built with "
      "BUBBLE_ENABLE_PARDISO=OFF");
#endif
}

PardisoLinearSolver::~PardisoLinearSolver() {
#if BUBBLE_PARDISO_ENABLED
  if (!impl_) return;
  if (impl_->analyzed) {
    const MklInt phase = -1, maxfct = 1, mnum = 1, nrhs = 1, msglvl = 0;
    MklInt error = 0, permutation = 0;
    double dummy = 0.0;
    impl_->pardiso_call(impl_->pt.data(), &maxfct, &mnum, &impl_->mtype,
                        &phase, &impl_->n, &dummy, impl_->row_ptr.data(),
                        impl_->columns.data(), &permutation, &nrhs,
                        impl_->iparm.data(), &msglvl, &dummy, &dummy, &error);
  }
  if (impl_->library) dlclose(impl_->library);
#endif
}

Eigen::VectorXd PardisoLinearSolver::solve(
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& input,
    const Eigen::VectorXd& rhs, PardisoSolveMetrics& metrics) {
#if BUBBLE_PARDISO_ENABLED
  const auto csr_started = Clock::now();
  Eigen::SparseMatrix<double, Eigen::RowMajor> original = input;
  original.makeCompressed();
  if (original.rows() != original.cols() || original.rows() != rhs.size() ||
      original.rows() > std::numeric_limits<MklInt>::max() ||
      original.nonZeros() > std::numeric_limits<MklInt>::max())
    throw std::runtime_error("PARDISO CSR dimensions exceed signed 32-bit range");
  const MklInt n = static_cast<MklInt>(original.rows());
  const MklInt block_dimension = n % 2 == 0 ? 2 : 1;
  const MklInt block_rows = n / block_dimension;
  std::vector<std::vector<MklInt>> block_columns(
      static_cast<std::size_t>(block_rows));
  for (MklInt row = 0; row < n; ++row)
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(
             original, row); it; ++it)
      block_columns[static_cast<std::size_t>(row / block_dimension)]
          .push_back(static_cast<MklInt>(it.col()) / block_dimension);
  std::int64_t canonical_nonzeros = 0;
  for (auto& row : block_columns) {
    std::sort(row.begin(), row.end());
    row.erase(std::unique(row.begin(), row.end()), row.end());
    canonical_nonzeros += static_cast<std::int64_t>(row.size()) *
                          block_dimension * block_dimension;
  }
  if (canonical_nonzeros > std::numeric_limits<MklInt>::max())
    throw std::runtime_error(
        "PARDISO canonical block CSR exceeds signed 32-bit range");
  const MklInt nnz = static_cast<MklInt>(canonical_nonzeros);
  std::vector<MklInt> row_ptr(static_cast<std::size_t>(n) + 1);
  std::vector<MklInt> columns;
  std::vector<double> values;
  columns.reserve(static_cast<std::size_t>(nnz));
  values.reserve(static_cast<std::size_t>(nnz));
  for (MklInt block_row = 0; block_row < block_rows; ++block_row) {
    for (MklInt local_row = 0; local_row < block_dimension; ++local_row) {
      const MklInt scalar_row = block_row * block_dimension + local_row;
      row_ptr[static_cast<std::size_t>(scalar_row)] =
          static_cast<MklInt>(columns.size());
      for (MklInt block_column :
           block_columns[static_cast<std::size_t>(block_row)]) {
        for (MklInt local_column = 0; local_column < block_dimension;
             ++local_column) {
          columns.push_back(block_column * block_dimension + local_column);
          values.push_back(0.0);
        }
      }
    }
  }
  row_ptr.back() = static_cast<MklInt>(columns.size());
  for (MklInt row = 0; row < n; ++row) {
    MklInt canonical = row_ptr[static_cast<std::size_t>(row)];
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(
             original, row); it; ++it) {
      while (canonical < row_ptr[static_cast<std::size_t>(row + 1)] &&
             columns[static_cast<std::size_t>(canonical)] < it.col())
        ++canonical;
      if (canonical == row_ptr[static_cast<std::size_t>(row + 1)] ||
          columns[static_cast<std::size_t>(canonical)] != it.col())
        throw std::runtime_error("PARDISO canonical CSR lost an input entry");
      values[static_cast<std::size_t>(canonical)] = it.value();
    }
  }
  if (row_ptr.empty() || row_ptr.front() != 0 || row_ptr.back() != nnz)
    throw std::runtime_error("PARDISO CSR row-pointer endpoints are invalid");
  for (MklInt row = 0; row < n; ++row) {
    if (row_ptr[static_cast<std::size_t>(row)] >
        row_ptr[static_cast<std::size_t>(row + 1)])
      throw std::runtime_error("PARDISO CSR row pointers are not monotone");
    MklInt previous = -1;
    for (MklInt k = row_ptr[static_cast<std::size_t>(row)];
         k < row_ptr[static_cast<std::size_t>(row + 1)]; ++k) {
      const MklInt column = columns[static_cast<std::size_t>(k)];
      if (column < 0 || column >= n || column <= previous ||
          !std::isfinite(values[static_cast<std::size_t>(k)]))
        throw std::runtime_error(
            "PARDISO CSR has an out-of-range, duplicate, unsorted, or nonfinite entry");
      previous = column;
    }
  }
  for (Eigen::Index i = 0; i < rhs.size(); ++i)
    if (!std::isfinite(rhs[i]))
      throw std::runtime_error("PARDISO RHS is nonfinite");
  metrics.csr_seconds += seconds(csr_started);

  const bool same_graph = impl_->analyzed && impl_->n == n && impl_->nnz == nnz &&
                          impl_->row_ptr == row_ptr && impl_->columns == columns;
  const MklInt maxfct = 1, mnum = 1, nrhs = 1, msglvl = 0;
  MklInt permutation = 0;
  auto call = [&](MklInt phase, const double* values, double* b, double* x) {
    MklInt error = 0;
    impl_->pardiso_call(impl_->pt.data(), &maxfct, &mnum, &impl_->mtype,
                        &phase, &n, values, row_ptr.data(), columns.data(),
                        &permutation, &nrhs, impl_->iparm.data(), &msglvl, b,
                        x, &error);
    if (error != 0) throw_error(error, phase);
  };

  if (!same_graph) {
    if (impl_->analyzed) {
      const MklInt release_phase = -1;
      MklInt error = 0;
      double dummy = 0.0;
      impl_->pardiso_call(
          impl_->pt.data(), &maxfct, &mnum, &impl_->mtype, &release_phase,
          &impl_->n, &dummy, impl_->row_ptr.data(), impl_->columns.data(),
          &permutation, &nrhs, impl_->iparm.data(), &msglvl, &dummy, &dummy,
          &error);
      impl_->pt.fill(nullptr);
      impl_->iparm.fill(0);
      impl_->pardiso_init(impl_->pt.data(), &impl_->mtype,
                          impl_->iparm.data());
      impl_->iparm[34] = 1;
      impl_->iparm[7] = 8;
      impl_->iparm[59] = 0;
      impl_->analyzed = false;
    }
    const auto analysis_started = Clock::now();
    double dummy = 0.0;
    call(11, values.data(), &dummy, &dummy);
    metrics.symbolic_analysis_seconds += seconds(analysis_started);
    metrics.symbolic_analysis_calls = 1;
    impl_->n = n;
    impl_->nnz = nnz;
    impl_->row_ptr = row_ptr;
    impl_->columns = columns;
    impl_->analyzed = true;
  } else {
    metrics.structure_reused = true;
  }

  const auto factor_started = Clock::now();
  double dummy = 0.0;
  call(22, values.data(), &dummy, &dummy);
  metrics.numeric_factorization_seconds += seconds(factor_started);
  metrics.numeric_factorization_calls = 1;

  Eigen::VectorXd solution = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd rhs_copy = rhs;
  const auto solve_started = Clock::now();
  call(33, values.data(), rhs_copy.data(), solution.data());
  metrics.solve_seconds += seconds(solve_started);
  metrics.iterative_refinement_steps = impl_->iparm[6];
  metrics.factor_nonzeros = impl_->iparm[17];
  metrics.peak_memory_kib = std::max<std::int64_t>(
      impl_->iparm[14], static_cast<std::int64_t>(impl_->iparm[15]) +
                            static_cast<std::int64_t>(impl_->iparm[16]));
  const Eigen::VectorXd residual = original * solution - rhs;
  metrics.cpu_relative_residual =
      residual.norm() / std::max(rhs.norm(), 1e-300);
  // Match Solver::attempt_step's backend-independent true-residual gate.
  // A stricter private PARDISO floor would create backend-specific timestep
  // rejections even when the common Newton linear acceptance gate passes.
  const double residual_limit =
      std::max(10.0 * impl_->relative_tolerance, 1e-8);
  if (!std::isfinite(metrics.cpu_relative_residual) ||
      metrics.cpu_relative_residual > residual_limit) {
    std::ostringstream message;
    message.precision(17);
    message << "PARDISO CPU-recomputed relative residual "
            << metrics.cpu_relative_residual << " exceeds " << residual_limit;
    throw std::runtime_error(message.str());
  }
  return solution;
#else
  (void)input;
  (void)rhs;
  (void)metrics;
  throw std::runtime_error("PARDISO backend unavailable in this binary");
#endif
}

PardisoRuntimeInfo pardiso_runtime_info(int configured_threads) {
  PardisoRuntimeInfo info;
#if BUBBLE_PARDISO_ENABLED
  info.compiled = true;
  info.threads = configured_threads;
  info.library_path = BUBBLE_PARDISO_LIBRARY_PATH;
#else
  (void)configured_threads;
#endif
  return info;
}

}  // namespace bubble
