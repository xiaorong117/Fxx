#include "bubble/model.hpp"

#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using json=nlohmann::json;

void restore(bubble::Solver& s,const bubble::CheckpointData& c){
  s.restore(c.state,c.time,c.step_number,c.cumulative_mole_ledger,c.cumulative_volume_ledger,
    c.cumulative_net_liquid_inflow,c.cumulative_boundary_component_outflow,
    c.cumulative_outlet_advective_component_outflow,c.cumulative_liquid_inflow,
    c.cumulative_liquid_outflow,c.cumulative_transfer_by_node);
  s.restore_gross_fluxes(c.cumulative_component_inflow,c.cumulative_component_outflow);
}

double main_mean_sg(const bubble::Solver&s){double v=0,g=0;for(std::size_t i=0;i<s.state().size();++i){
  const auto&p=s.network().pores[i];if(p.source_type==2||p.volume==0||s.network().component_id[i]!=s.network().main_flow_component_id)continue;
  v+=p.volume;g+=p.volume*static_cast<double>(1-s.state()[i].sw);}return g/v;}

int main(int argc,char**argv){
  if(argc!=5&&argc!=6)throw std::runtime_error("usage: compare config checkpoint amgx_config output.json [dt]");
  auto eigen_config=bubble::Config::load(argv[1]);auto amgx_config=eigen_config;
  amgx_config.linear.backend="amgx";amgx_config.linear.relative_tolerance=1e-10;amgx_config.linear.amgx_config_file=argv[3];
  const auto checkpoint=bubble::read_checkpoint(argv[2],bubble::Network::load(eigen_config).pores.size());
  auto network=bubble::Network::load(eigen_config);bubble::Solver eigen(network,eigen_config),amgx(std::move(network),amgx_config);
  eigen.initialize();amgx.initialize();restore(eigen,checkpoint);restore(amgx,checkpoint);
  const double dt=argc==6?std::stod(argv[5]):std::min(0.1,eigen_config.time.max_dt);const auto er=eigen.attempt_step(dt);const auto ar=amgx.attempt_step(dt);
  if(!er.accepted||!ar.accepted)throw std::runtime_error("backend step failed: Eigen="+er.failure+" AMGX="+ar.failure);
  std::array<double,4> variable_max{{0,0,0,0}};double max_scaled=0;std::size_t worst=0;for(std::size_t i=0;i<eigen.state().size();++i){const auto&a=eigen.state()[i];const auto&b=amgx.state()[i];
    const std::array<double,4> differences{{std::abs(double(a.pl-b.pl))/eigen_config.scaling.pressure,
      a.active_gas?std::abs(double(a.pg-b.pg))/eigen_config.scaling.gas_pressure:0.0,
      std::abs(double(a.sw-b.sw)),std::abs(double(a.c-b.c))/eigen_config.scaling.concentration}};
    for(std::size_t j=0;j<4;++j)variable_max[j]=std::max(variable_max[j],differences[j]);
    const double value=*std::max_element(differences.begin(),differences.end());
    if(value>max_scaled){max_scaled=value;worst=eigen.network().pores[i].id;}}
  const double mean_e=main_mean_sg(eigen),mean_a=main_mean_sg(amgx);
  const auto rel=[](double a,double b){return std::abs(a-b)/std::max({std::abs(a),std::abs(b),1e-300});};
  json out={{"checkpoint",std::filesystem::absolute(argv[2]).string()},{"dt_s",dt},
    {"dof_count",er.final_assembly.jacobian.rows()},{"jacobian_nonzeros",er.final_assembly.jacobian.nonZeros()},
    {"csr_max_column_index",er.final_assembly.jacobian.cols()-1},{"max_scaled_state_difference",max_scaled},
    {"maximum_scaled_variable_difference",{{"Pl",variable_max[0]},{"Pg",variable_max[1]},
      {"Sw",variable_max[2]},{"C",variable_max[3]}}},
    {"worst_state_node_id",worst},{"main_mean_Sg",{{"eigen",mean_e},{"amgx",mean_a},{"relative_difference",rel(mean_e,mean_a)}}},
    {"free_gas_mol",{{"eigen",er.balance.total_free},{"amgx",ar.balance.total_free},{"relative_difference",rel(er.balance.total_free,ar.balance.total_free)}}},
    {"dissolved_gas_mol",{{"eigen",er.balance.total_dissolved},{"amgx",ar.balance.total_dissolved},{"relative_difference",rel(er.balance.total_dissolved,ar.balance.total_dissolved)}}},
    {"residual_blocks_eigen",er.final_assembly.residual_norms},{"residual_blocks_amgx",ar.final_assembly.residual_norms},
    {"eigen",{{"nonlinear_residual",er.residual_norm},{"linear_residual",er.linear_residual},{"newton_iterations",er.newton_iterations},
      {"linear_iterations",er.linear_iterations},{"relative_update_norm",er.relative_update_norm},{"assembly_seconds",er.assembly_seconds},{"setup_seconds",er.linear_setup_seconds},{"solve_seconds",er.linear_solve_seconds},
      {"epsilon_N",er.balance.epsilon_moles_adjusted},{"epsilon_V",er.balance.epsilon_volume_adjusted}}},
    {"amgx",{{"nonlinear_residual",ar.residual_norm},{"linear_residual",ar.linear_residual},{"newton_iterations",ar.newton_iterations},
      {"linear_iterations",ar.linear_iterations},{"amgx_status_code",ar.amgx_status},
      {"amgx_status_name",ar.amgx_status==0?"AMGX_SOLVE_SUCCESS":"NOT_SUCCESS"},
      {"matrix_upload_all_calls",ar.amgx_upload_all_calls},
      {"matrix_replace_coefficients_calls",ar.amgx_replace_coefficients_calls},
      {"solver_setup_calls",ar.amgx_solver_setup_calls},
      {"relative_update_norm",ar.relative_update_norm},{"assembly_seconds",ar.assembly_seconds},{"csr_seconds",ar.csr_seconds},
      {"upload_seconds",ar.linear_upload_seconds},{"setup_seconds",ar.linear_setup_seconds},{"solve_seconds",ar.linear_solve_seconds},
      {"download_seconds",ar.linear_download_seconds},{"epsilon_N",ar.balance.epsilon_moles_adjusted},{"epsilon_V",ar.balance.epsilon_volume_adjusted}}}};
  out["acceptance"]={{"state_difference_pass",max_scaled<=1e-7},{"mean_Sg_pass",rel(mean_e,mean_a)<=1e-8},
    {"free_moles_pass",rel(er.balance.total_free,ar.balance.total_free)<=1e-8},
    {"dissolved_moles_pass",rel(er.balance.total_dissolved,ar.balance.total_dissolved)<=1e-8},
    {"conservation_pass",er.balance.epsilon_moles_adjusted<=eigen_config.acceptance.mass_balance_hard_limit&&
       ar.balance.epsilon_moles_adjusted<=eigen_config.acceptance.mass_balance_hard_limit&&
       er.balance.epsilon_volume_adjusted<=eigen_config.acceptance.volume_balance_hard_limit&&
       ar.balance.epsilon_volume_adjusted<=eigen_config.acceptance.volume_balance_hard_limit}};
  std::ofstream file(argv[4]);file<<std::setw(2)<<out<<'\n';std::cout<<out.dump(2)<<'\n';
  return std::all_of(out["acceptance"].begin(),out["acceptance"].end(),[](const auto&x){return x.template get<bool>();})?0:2;
}
