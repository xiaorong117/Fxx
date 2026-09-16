#include "bubble/amgx_backend.hpp"

#include <chrono>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#if BUBBLE_AMGX_ENABLED
#include <amgx_c.h>
#endif

namespace bubble {
namespace {
using Clock = std::chrono::steady_clock;
#if BUBBLE_AMGX_ENABLED
double seconds(Clock::time_point begin) {
  return std::chrono::duration<double>(Clock::now() - begin).count();
}
#endif
#if BUBBLE_AMGX_ENABLED
const char* return_code_name(AMGX_RC code) {
  switch (code) {
    case AMGX_RC_OK: return "AMGX_RC_OK";
    case AMGX_RC_BAD_PARAMETERS: return "AMGX_RC_BAD_PARAMETERS";
    case AMGX_RC_UNKNOWN: return "AMGX_RC_UNKNOWN";
    case AMGX_RC_NOT_SUPPORTED_TARGET: return "AMGX_RC_NOT_SUPPORTED_TARGET";
    case AMGX_RC_NOT_SUPPORTED_BLOCKSIZE: return "AMGX_RC_NOT_SUPPORTED_BLOCKSIZE";
    case AMGX_RC_CUDA_FAILURE: return "AMGX_RC_CUDA_FAILURE";
    case AMGX_RC_THRUST_FAILURE: return "AMGX_RC_THRUST_FAILURE";
    case AMGX_RC_NO_MEMORY: return "AMGX_RC_NO_MEMORY";
    case AMGX_RC_IO_ERROR: return "AMGX_RC_IO_ERROR";
    case AMGX_RC_BAD_MODE: return "AMGX_RC_BAD_MODE";
    case AMGX_RC_CORE: return "AMGX_RC_CORE";
    case AMGX_RC_PLUGIN: return "AMGX_RC_PLUGIN";
    case AMGX_RC_BAD_CONFIGURATION: return "AMGX_RC_BAD_CONFIGURATION";
    case AMGX_RC_NOT_IMPLEMENTED: return "AMGX_RC_NOT_IMPLEMENTED";
    case AMGX_RC_LICENSE_NOT_FOUND: return "AMGX_RC_LICENSE_NOT_FOUND";
    case AMGX_RC_INTERNAL: return "AMGX_RC_INTERNAL";
  }
  return "AMGX_RC_UNKNOWN_VALUE";
}
void check(AMGX_RC code, const char* operation) {
  if (code != AMGX_RC_OK)
    throw std::runtime_error(std::string(operation) + " failed with " +
                             return_code_name(code) + " (" +
                             std::to_string(static_cast<int>(code)) + ")");
}
const char* status_name(AMGX_SOLVE_STATUS status){
  switch(status){case AMGX_SOLVE_SUCCESS:return "AMGX_SOLVE_SUCCESS";
    case AMGX_SOLVE_FAILED:return "AMGX_SOLVE_FAILED";
    case AMGX_SOLVE_DIVERGED:return "AMGX_SOLVE_DIVERGED";
    case AMGX_SOLVE_NOT_CONVERGED:return "AMGX_SOLVE_NOT_CONVERGED";}
  return "AMGX_SOLVE_STATUS_UNKNOWN";
}
#endif
}  // namespace

struct AmgxLinearSolver::Impl {
#if BUBBLE_AMGX_ENABLED
  AMGX_config_handle config = nullptr;
  AMGX_resources_handle resources = nullptr;
  AMGX_matrix_handle matrix = nullptr;
  AMGX_vector_handle rhs = nullptr, solution = nullptr;
  AMGX_solver_handle solver = nullptr;
  std::vector<int> row_ptr, columns;
  int rows = 0, nnz = 0, block_dim = 1;
  bool graph_uploaded = false;
  bool dump_written = false;
#endif
};

AmgxLinearSolver::AmgxLinearSolver(const std::filesystem::path& config_file)
    : impl_(std::make_unique<Impl>()) {
#if BUBBLE_AMGX_ENABLED
  if (!std::filesystem::is_regular_file(config_file))
    throw std::runtime_error("AMGX config file is missing: " + config_file.string());
  check(AMGX_initialize(), "AMGX_initialize");
  try {
    check(AMGX_config_create_from_file(&impl_->config, config_file.c_str()),
          "AMGX_config_create_from_file");
    check(AMGX_resources_create_simple(&impl_->resources, impl_->config),
          "AMGX_resources_create_simple");
    check(AMGX_matrix_create(&impl_->matrix, impl_->resources, AMGX_mode_dDDI),
          "AMGX_matrix_create");
    check(AMGX_vector_create(&impl_->rhs, impl_->resources, AMGX_mode_dDDI),
          "AMGX_vector_create(rhs)");
    check(AMGX_vector_create(&impl_->solution, impl_->resources, AMGX_mode_dDDI),
          "AMGX_vector_create(solution)");
    check(AMGX_solver_create(&impl_->solver, impl_->resources, AMGX_mode_dDDI,
                             impl_->config), "AMGX_solver_create");
  } catch (...) {
    if (impl_->solver) AMGX_solver_destroy(impl_->solver);
    if (impl_->solution) AMGX_vector_destroy(impl_->solution);
    if (impl_->rhs) AMGX_vector_destroy(impl_->rhs);
    if (impl_->matrix) AMGX_matrix_destroy(impl_->matrix);
    if (impl_->resources) AMGX_resources_destroy(impl_->resources);
    if (impl_->config) AMGX_config_destroy(impl_->config);
    AMGX_finalize();
    throw;
  }
#else
  (void)config_file;
  throw std::runtime_error(
      "linear.backend=amgx requested, but this binary was built with BUBBLE_ENABLE_AMGX=OFF");
#endif
}

AmgxLinearSolver::~AmgxLinearSolver() {
#if BUBBLE_AMGX_ENABLED
  if (!impl_) return;
  if (impl_->solver) AMGX_solver_destroy(impl_->solver);
  if (impl_->solution) AMGX_vector_destroy(impl_->solution);
  if (impl_->rhs) AMGX_vector_destroy(impl_->rhs);
  if (impl_->matrix) AMGX_matrix_destroy(impl_->matrix);
  if (impl_->resources) AMGX_resources_destroy(impl_->resources);
  if (impl_->config) AMGX_config_destroy(impl_->config);
  AMGX_finalize();
#endif
}

Eigen::VectorXd AmgxLinearSolver::solve(
    const Eigen::SparseMatrix<double, Eigen::RowMajor>& input,
    const Eigen::VectorXd& rhs, AmgxSolveMetrics& metrics) {
#if BUBBLE_AMGX_ENABLED
  const auto csr_started = Clock::now();
  Eigen::SparseMatrix<double, Eigen::RowMajor> original = input;
  original.makeCompressed();
  Eigen::SparseMatrix<double, Eigen::RowMajor> matrix = original;
  matrix.makeCompressed();
  if (matrix.rows() != matrix.cols() || matrix.rows() != rhs.size() ||
      matrix.rows() > std::numeric_limits<int>::max() ||
      matrix.nonZeros() > std::numeric_limits<int>::max())
    throw std::runtime_error("AMGX CSR dimensions exceed signed 32-bit range");
  const int n = static_cast<int>(matrix.rows());
  const int nnz = static_cast<int>(matrix.nonZeros());
  Eigen::VectorXd equilibrated_rhs=rhs;
  std::vector<double> row_scale(static_cast<std::size_t>(n),1.0),
                      column_scale(static_cast<std::size_t>(n),1.0);
  for(int sweep=0;sweep<6;++sweep){
    for(int row=0;row<n;++row){double maximum=0.0;
      for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(matrix,row);it;++it)
        maximum=std::max(maximum,std::abs(it.value()));
      if(!(maximum>0.0)||!std::isfinite(maximum))throw std::runtime_error("AMGX matrix has zero/nonfinite row");
      const double factor=1.0/std::sqrt(maximum);row_scale[static_cast<std::size_t>(row)]*=factor;
      equilibrated_rhs[row]*=factor;
      for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(matrix,row);it;++it)it.valueRef()*=factor;
    }
    std::vector<double> maxima(static_cast<std::size_t>(n),0.0);
    for(int row=0;row<n;++row)for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(matrix,row);it;++it)
      maxima[static_cast<std::size_t>(it.col())]=std::max(maxima[static_cast<std::size_t>(it.col())],std::abs(it.value()));
    for(int col=0;col<n;++col){if(!(maxima[static_cast<std::size_t>(col)]>0.0))throw std::runtime_error("AMGX matrix has zero column");
      column_scale[static_cast<std::size_t>(col)]*=1.0/std::sqrt(maxima[static_cast<std::size_t>(col)]);}
    for(int row=0;row<n;++row)for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(matrix,row);it;++it)
      it.valueRef()/=std::sqrt(maxima[static_cast<std::size_t>(it.col())]);
  }
  const int block_dim=n%2==0?2:1;
  const int block_rows=n/block_dim;
  std::vector<int> row_ptr(static_cast<std::size_t>(block_rows)+1),columns;
  std::vector<double> values;
  std::vector<std::map<int,std::array<double,4>>> blocks(static_cast<std::size_t>(block_rows));
  for (int row = 0; row < n; ++row) {
    int previous = -1;
    for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(matrix,row);it;++it){
      const auto column=it.col();const double value=it.value();
      if (column < 0 || column >= n || column <= previous || !std::isfinite(value))
        throw std::runtime_error("AMGX CSR has out-of-range, duplicate, unsorted, or nonfinite entry");
      previous = static_cast<int>(column);
      auto& block=blocks[static_cast<std::size_t>(row/block_dim)][static_cast<int>(column)/block_dim];
      block[static_cast<std::size_t>((row%block_dim)*block_dim+(column%block_dim))]=value;
    }
  }
  for(int row=0;row<block_rows;++row){row_ptr[static_cast<std::size_t>(row)]=static_cast<int>(columns.size());
    for(const auto&item:blocks[static_cast<std::size_t>(row)]){columns.push_back(item.first);
      for(int k=0;k<block_dim*block_dim;++k)values.push_back(item.second[static_cast<std::size_t>(k)]);}}
  row_ptr.back()=static_cast<int>(columns.size());const int block_nnz=static_cast<int>(columns.size());
  if(row_ptr.front()!=0||row_ptr.back()!=block_nnz)throw std::runtime_error("AMGX block CSR endpoints invalid");
  for (int i = 0; i < n; ++i)
    if (!std::isfinite(rhs[i])) throw std::runtime_error("AMGX RHS is nonfinite");
  if(std::getenv("BUBBLE_AMGX_MATRIX_STATS")){
    double min_abs=std::numeric_limits<double>::infinity(),max_abs=0,min_diag=min_abs,max_diag=0;
    int missing_diag=0;for(int row=0;row<n;++row){bool found=false;
      for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(matrix,row);it;++it){
        const double a=std::abs(it.value());if(a>0){min_abs=std::min(min_abs,a);max_abs=std::max(max_abs,a);}
        if(it.col()==row){found=true;min_diag=std::min(min_diag,a);max_diag=std::max(max_diag,a);}}
      if(!found)++missing_diag;}
    std::cerr.precision(17);std::cerr<<"AMGX_MATRIX_STATS rows="<<n<<" nnz="<<nnz
      <<" missing_diagonal="<<missing_diag<<" min_abs="<<min_abs<<" max_abs="<<max_abs
      <<" min_abs_diagonal="<<min_diag<<" max_abs_diagonal="<<max_diag
      <<" rhs_norm="<<rhs.norm()<<" rhs_inf="<<rhs.lpNorm<Eigen::Infinity>()<<'\n';
  }
  metrics.csr_seconds += seconds(csr_started);
  const char* dump_prefix=std::getenv("BUBBLE_AMGX_DUMP_PREFIX");
  const bool dump_this=dump_prefix&&!impl_->dump_written;
  if(dump_this){
    std::ofstream scaling(std::string(dump_prefix)+"_amgx_ruiz_scaling.tsv");
    scaling<<"scalar_index_zero_based\trow_scale\tcolumn_scale\n"<<std::setprecision(17);
    for(int i=0;i<n;++i)scaling<<i<<'\t'<<row_scale[static_cast<std::size_t>(i)]<<'\t'<<column_scale[static_cast<std::size_t>(i)]<<'\n';
    std::ofstream block(std::string(dump_prefix)+"_block_csr_audit.json");
    block<<"{\n  \"scalar_rows\": "<<n<<",\n  \"scalar_nonzeros\": "<<nnz
      <<",\n  \"block_rows\": "<<block_rows<<",\n  \"block_nonzeros\": "<<block_nnz
      <<",\n  \"block_dimension\": "<<block_dim<<",\n  \"runtime_index_base\": 0,\n  \"integer_width_bits\": 32,\n  \"minimum_column\": 0,\n  \"maximum_column\": "<<block_rows-1<<"\n}\n";
  }

  const bool reuse_structure = std::getenv("BUBBLE_AMGX_DISABLE_STRUCTURE_REUSE") == nullptr;
  const bool same_graph = reuse_structure && impl_->graph_uploaded && impl_->rows == block_rows && impl_->nnz == block_nnz &&
                          impl_->block_dim == block_dim &&
                          impl_->row_ptr == row_ptr && impl_->columns == columns;
  const auto upload_started = Clock::now();
  if (same_graph) {
    check(AMGX_matrix_replace_coefficients(impl_->matrix, block_rows, block_nnz, values.data(), nullptr),
          "AMGX_matrix_replace_coefficients");
    metrics.matrix_replace_coefficients_calls=1;
  } else {
    check(AMGX_matrix_upload_all(impl_->matrix, block_rows, block_nnz, block_dim, block_dim, row_ptr.data(),
                                 columns.data(), values.data(), nullptr),
          "AMGX_matrix_upload_all");
    impl_->rows=block_rows; impl_->nnz=block_nnz; impl_->block_dim=block_dim; impl_->row_ptr=std::move(row_ptr);
    impl_->columns=std::move(columns); impl_->graph_uploaded=true;
    metrics.matrix_upload_all_calls=1;
  }
  check(AMGX_vector_upload(impl_->rhs, block_rows, block_dim, equilibrated_rhs.data()), "AMGX_vector_upload");
  check(AMGX_vector_set_zero(impl_->solution, block_rows, block_dim), "AMGX_vector_set_zero");
  metrics.upload_seconds += seconds(upload_started);

  const auto setup_started = Clock::now();
  check(AMGX_solver_setup(impl_->solver, impl_->matrix), "AMGX_solver_setup");
  metrics.solver_setup_calls=1;
  metrics.setup_seconds += seconds(setup_started);
  Eigen::VectorXd answer=Eigen::VectorXd::Zero(n);
  AMGX_SOLVE_STATUS status = AMGX_SOLVE_FAILED;
  for(int refinement=0;refinement<1;++refinement){
    const auto solve_started = Clock::now();
    check(AMGX_solver_solve_with_0_initial_guess(impl_->solver, impl_->rhs,
                                                 impl_->solution), "AMGX_solver_solve");
    metrics.solve_seconds += seconds(solve_started);
    check(AMGX_solver_get_status(impl_->solver, &status), "AMGX_solver_get_status");
    int iterations=0;check(AMGX_solver_get_iterations_number(impl_->solver, &iterations),
                          "AMGX_solver_get_iterations_number");
    metrics.iterations+=iterations;
    metrics.amgx_status = static_cast<int>(status);
    if(dump_this){
      std::ofstream history(std::string(dump_prefix)+"_amgx_residual_history.tsv");
      history<<"iteration\tamgx_reported_residual\n"<<std::setprecision(17);
      for(int iteration=0;iteration<=iterations;++iteration){double value=0;
        if(AMGX_solver_get_iteration_residual(impl_->solver,iteration,0,&value)==AMGX_RC_OK)
          history<<iteration<<'\t'<<value<<'\n';}
      std::ofstream result(std::string(dump_prefix)+"_amgx_result.json");
      result<<"{\n  \"status_code\": "<<static_cast<int>(status)<<",\n  \"status_name\": \""<<status_name(status)
        <<"\",\n  \"iterations\": "<<iterations<<"\n}\n";
      impl_->dump_written=true;
    }
    Eigen::VectorXd correction(n);
    const auto download_started = Clock::now();
    check(AMGX_vector_download(impl_->solution, correction.data()), "AMGX_vector_download");
    metrics.download_seconds += seconds(download_started);
    for(int i=0;i<n;++i)correction[i]*=column_scale[static_cast<std::size_t>(i)];
    answer+=correction;
    const Eigen::VectorXd residual=rhs-original*answer;
    metrics.cpu_relative_residual=residual.norm()/std::max(rhs.norm(),1e-300);
    if(dump_this){
      std::ofstream solution(std::string(dump_prefix)+"_amgx_solution.mtx");
      solution<<"%%MatrixMarket matrix array real general\n"<<n<<" 1\n"<<std::setprecision(17);
      for(int i=0;i<n;++i)solution<<answer[i]<<'\n';
      std::ofstream result(std::string(dump_prefix)+"_amgx_result.json");
      result<<std::setprecision(17)<<"{\n  \"status_code\": "<<static_cast<int>(status)
        <<",\n  \"status_name\": \""<<status_name(status)<<"\",\n  \"iterations\": "<<iterations
        <<",\n  \"cpu_recomputed_relative_residual\": "<<metrics.cpu_relative_residual<<"\n}\n";
    }
    // AMGX's own convergence decision and the independently recomputed CPU
    // residual are separate acceptance gates.  Never turn NOT_CONVERGED into
    // SUCCESS merely because the CPU-side residual happens to pass.
    if(status==AMGX_SOLVE_SUCCESS && metrics.cpu_relative_residual<=1e-8)break;
    Eigen::VectorXd scaled_residual=residual;
    for(int i=0;i<n;++i)scaled_residual[i]*=row_scale[static_cast<std::size_t>(i)];
    const auto correction_upload=Clock::now();
    check(AMGX_vector_upload(impl_->rhs,block_rows,block_dim,scaled_residual.data()),"AMGX_vector_upload(refinement)");
    check(AMGX_vector_set_zero(impl_->solution,block_rows,block_dim),"AMGX_vector_set_zero(refinement)");
    metrics.upload_seconds+=seconds(correction_upload);
  }
  if (!std::isfinite(metrics.cpu_relative_residual))
    throw std::runtime_error("AMGX CPU-recomputed residual is nonfinite");
  if (status != AMGX_SOLVE_SUCCESS || metrics.cpu_relative_residual > 1e-8)
  { std::ostringstream message;message.precision(17);
    message<<"AMGX solver status="<<status_name(status)<<" ("<<static_cast<int>(status)<<')'
      <<" after iterative refinement iterations="<<metrics.iterations
      <<" CPU relative residual="<<metrics.cpu_relative_residual;
    throw std::runtime_error(message.str()); }
  return answer;
#else
  (void)input; (void)rhs; (void)metrics;
  throw std::runtime_error("AMGX backend unavailable in this binary");
#endif
}

AmgxRuntimeInfo amgx_runtime_info() {
  AmgxRuntimeInfo info;
#if BUBBLE_AMGX_ENABLED
  info.compiled = true;
  AMGX_get_api_version(&info.api_major, &info.api_minor);
  info.library_path = BUBBLE_AMGX_LIBRARY_PATH;
#endif
  return info;
}

}  // namespace bubble
