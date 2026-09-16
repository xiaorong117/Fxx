#include "bubble/amgx_backend.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc,char**argv){
  if(argc!=2)throw std::runtime_error("AMGX config path required");
  Eigen::SparseMatrix<double,Eigen::RowMajor> a(5,5);
  std::vector<Eigen::Triplet<double>> t;
  for(int i=0;i<5;++i){t.emplace_back(i,i,4.0+i);if(i){t.emplace_back(i,i-1,-1.2);t.emplace_back(i-1,i,-0.7);}}
  a.setFromTriplets(t.begin(),t.end());a.makeCompressed();
  Eigen::VectorXd expected(5);expected<<1,-2,3,-4,5;const Eigen::VectorXd rhs=a*expected;
  bubble::AmgxLinearSolver solver(argv[1]);bubble::AmgxSolveMetrics metrics;
  const auto result=solver.solve(a,rhs,metrics);
  const double residual=(a*result-rhs).norm()/rhs.norm();
  if(metrics.amgx_status!=0||metrics.cpu_relative_residual>1e-9||residual>1e-9||
     (result-expected).norm()>1e-8)
    throw std::runtime_error("AMGX small system CPU residual/state check failed");
  std::cout<<"PASS AMGX small nonsymmetric system iterations="<<metrics.iterations
           <<" cpu_relative_residual="<<metrics.cpu_relative_residual<<'\n';
}
