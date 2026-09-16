#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>
#include <nlohmann/json.hpp>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc,char**argv){if(argc!=4)throw std::runtime_error("usage: matrix rhs result.json");
  Eigen::SparseMatrix<double,Eigen::RowMajor>a;Eigen::VectorXd b;
  if(!Eigen::loadMarket(a,argv[1])||!Eigen::loadMarketVector(b,argv[2]))throw std::runtime_error("MatrixMarket load failed");
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>,Eigen::IncompleteLUT<double>> solver;
  solver.setMaxIterations(5000);solver.setTolerance(1e-14);solver.preconditioner().setDroptol(1e-8);solver.preconditioner().setFillfactor(20);
  const auto begin=std::chrono::steady_clock::now();solver.compute(a);const auto setup=std::chrono::steady_clock::now();Eigen::VectorXd x=solver.solve(b);const auto end=std::chrono::steady_clock::now();
  const double residual=(a*x-b).norm()/std::max(b.norm(),1e-300);
  nlohmann::json result={{"backend","eigen_bicgstab_ilut"},{"status",solver.info()==Eigen::Success?"SUCCESS":"FAILED"},
    {"iterations",solver.iterations()},{"reported_error",solver.error()},{"cpu_recomputed_relative_residual",residual},
    {"setup_seconds",std::chrono::duration<double>(setup-begin).count()},{"solve_seconds",std::chrono::duration<double>(end-setup).count()}};
  std::ofstream out(argv[3]);out<<std::setw(2)<<result<<'\n';std::cout<<result.dump(2)<<'\n';return solver.info()==Eigen::Success&&residual<=1e-8?0:2;}
