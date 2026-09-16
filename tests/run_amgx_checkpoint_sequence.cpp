#include "bubble/model.hpp"
#include "bubble/amgx_backend.hpp"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <stdexcept>

using json=nlohmann::json;
void atomic(const std::filesystem::path&p,const json&j){auto q=p;q+=".tmp";{std::ofstream o(q);o<<std::setw(2)<<j<<'\n';}std::filesystem::rename(q,p);}
int main(int argc,char**argv){if(argc!=7&&argc!=8)throw std::runtime_error("usage: base_config checkpoint amgx_config output steps dt [amgx|eigen]");
  const auto output=std::filesystem::path(argv[4]);const int count=std::stoi(argv[5]);const double requested_dt=std::stod(argv[6]);
  const std::string backend=argc==8?argv[7]:"amgx";if(backend!="amgx"&&backend!="eigen")throw std::runtime_error("backend must be amgx or eigen");
  if(std::filesystem::exists(output))throw std::runtime_error("comparison output already exists");
  std::filesystem::create_directories(output);
  std::ifstream in(argv[1]);json resolved=json::parse(in);resolved["io"]["output_directory"]=output.string();resolved["linear"]["backend"]="amgx";
  if(backend=="amgx"){resolved["linear"]["backend"]="amgx";resolved["linear"]["amgx_config_file"]=std::filesystem::absolute(argv[3]).string();}
  else {resolved["linear"]["backend"]="eigen_bicgstab_ilut";resolved["linear"].erase("amgx_config_file");}
  resolved["time"]["t_end_s"]=1e6;
  const auto config_path=output/"resolved_config.json";{std::ofstream o(config_path);o<<std::setw(2)<<resolved<<'\n';}
  auto config=bubble::Config::load(config_path);auto network=bubble::Network::load(config);bubble::Solver solver(network,config);solver.initialize();
  const auto cp=bubble::read_checkpoint(argv[2],network.pores.size());solver.restore(cp.state,cp.time,cp.step_number,cp.cumulative_mole_ledger,
    cp.cumulative_volume_ledger,cp.cumulative_net_liquid_inflow,cp.cumulative_boundary_component_outflow,
    cp.cumulative_outlet_advective_component_outflow,cp.cumulative_liquid_inflow,cp.cumulative_liquid_outflow,cp.cumulative_transfer_by_node);
  solver.restore_gross_fluxes(cp.cumulative_component_inflow,cp.cumulative_component_outflow);
  std::ofstream log(output/"run.log");log<<std::scientific<<std::setprecision(17);
  std::ofstream table(output/"steps.tsv");table<<"accepted_index\tglobal_step\ttime_s\trequested_dt_s\taccepted_dt_s\tretries\tstep_wall_seconds\tnewton_iterations\tlinear_iterations\tamgx_status_name\tcpu_linear_relative_residual\tresidual_scaled_inf\trelative_update_norm\tepsilon_N\tepsilon_V\tassembly_seconds\tcsr_seconds\tupload_seconds\tsetup_seconds\tsolve_seconds\tdownload_seconds\n"<<std::scientific<<std::setprecision(17);
  const auto sequence_started=std::chrono::steady_clock::now();
  json steps=json::array();for(int accepted=1;accepted<=count;++accepted){atomic(output/"progress.json",{{"stage","amgx_checkpoint_sequence"},{"status","running"},{"accepted",accepted-1},{"target",count},{"time_s",solver.time()}});
    const auto step_started=std::chrono::steady_clock::now();auto r=solver.advance_with_retry(requested_dt);const double step_wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-step_started).count();if(!r.accepted)throw std::runtime_error(backend+" sequence failed: "+r.failure);
    if((backend=="amgx"&&r.amgx_status!=0)||(backend=="eigen"&&r.amgx_status!=-1)||r.linear_residual>std::max(10*config.linear.relative_tolerance,1e-8)||!r.residual_criterion_passed||!r.update_criterion_passed)
      throw std::runtime_error("backend/Newton acceptance invariant failed");
    const char* status=backend=="amgx"?"AMGX_SOLVE_SUCCESS":"NOT_APPLICABLE";
    table<<accepted<<'\t'<<solver.step_number()<<'\t'<<solver.time()<<'\t'<<requested_dt<<'\t'<<r.dt<<'\t'<<r.timestep_retries<<'\t'<<step_wall<<'\t'<<r.newton_iterations<<'\t'<<r.linear_iterations<<'\t'<<status<<'\t'<<r.linear_residual<<'\t'<<r.residual_norm<<'\t'<<r.relative_update_norm<<'\t'<<r.balance.epsilon_moles_adjusted<<'\t'<<r.balance.epsilon_volume_adjusted<<'\t'<<r.assembly_seconds<<'\t'<<r.csr_seconds<<'\t'<<r.linear_upload_seconds<<'\t'<<r.linear_setup_seconds<<'\t'<<r.linear_solve_seconds<<'\t'<<r.linear_download_seconds<<'\n';table.flush();
    steps.push_back({{"accepted_index",accepted},{"global_step",solver.step_number()},{"time_s",solver.time()},{"accepted_dt_s",r.dt},
      {"retries",r.timestep_retries},{"step_wall_seconds",step_wall},{"amgx_status",status},{"linear_residual",r.linear_residual},
      {"nonlinear_residual",r.residual_norm},{"relative_update_norm",r.relative_update_norm},
      {"epsilon_N",r.balance.epsilon_moles_adjusted},{"epsilon_V",r.balance.epsilon_volume_adjusted}});log<<steps.back().dump()<<'\n';log.flush();}
  const auto provenance=bubble::make_checkpoint_provenance(config);bubble::write_checkpoint(output/"checkpoint.json",solver,requested_dt,0,provenance);
  const double total_wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-sequence_started).count();json manifest={{"status","passed"},{"backend",backend},{"total_five_step_wall_seconds",total_wall},{"amgx_status_required",backend=="amgx"?"AMGX_SOLVE_SUCCESS":"NOT_APPLICABLE"},
    {"accepted_steps",count},{"initial_checkpoint",std::filesystem::absolute(argv[2]).string()},
    {"initial_checkpoint_sha256",bubble::sha256_file(argv[2])},{"resolved_config",config_path.string()},
    {"resolved_config_sha256",bubble::sha256_file(config_path)},{"amgx_config",std::filesystem::absolute(argv[3]).string()},
    {"amgx_config_sha256",bubble::sha256_file(argv[3])},{"final_checkpoint_sha256",bubble::sha256_file(output/"checkpoint.json")},
    {"steps",steps}};atomic(output/"run_manifest.json",manifest);atomic(output/"progress.json",{{"stage","complete"},{"status","passed"},{"accepted",count},{"time_s",solver.time()},{"checkpoint",(output/"checkpoint.json").string()}});
  std::cout<<manifest.dump(2)<<'\n';}
