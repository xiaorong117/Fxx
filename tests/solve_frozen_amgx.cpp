#include "bubble/amgx_backend.hpp"
#include <unsupported/Eigen/SparseExtra>
#include <nlohmann/json.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc,char**argv){if(argc!=5)throw std::runtime_error("usage: matrix rhs config result.json");
  nlohmann::json result={{"backend","amgx"},{"config",std::filesystem::absolute(argv[3]).string()}};
  try{Eigen::SparseMatrix<double,Eigen::RowMajor>a;Eigen::VectorXd b;
    if(!Eigen::loadMarket(a,argv[1])||!Eigen::loadMarketVector(b,argv[2]))throw std::runtime_error("MatrixMarket load failed");
    bubble::AmgxLinearSolver solver(argv[3]);bubble::AmgxSolveMetrics m;const auto x=solver.solve(a,b,m);(void)x;
    result.update({{"status","SUCCESS"},{"iterations",m.iterations},{"cpu_recomputed_relative_residual",m.cpu_relative_residual},
      {"csr_seconds",m.csr_seconds},{"upload_seconds",m.upload_seconds},{"setup_seconds",m.setup_seconds},
      {"solve_seconds",m.solve_seconds},{"download_seconds",m.download_seconds}});
  }catch(const std::exception&e){result.update({{"status","FAILED"},{"error",e.what()}});}
  std::ofstream out(argv[4]);out<<std::setw(2)<<result<<'\n';std::cout<<result.dump(2)<<'\n';return result["status"]=="SUCCESS"?0:2;}
