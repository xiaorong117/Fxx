#include "bubble/analysis.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <stdexcept>

namespace bubble {
namespace {
double quantile(std::vector<double> values, double p) {
  if (values.empty()) return 0.0;
  std::sort(values.begin(), values.end());
  const double x=p*(values.size()-1); const auto lo=static_cast<std::size_t>(x);
  const auto hi=std::min(lo+1,values.size()-1); return values[lo]+(x-lo)*(values[hi]-values[lo]);
}
double mean(const std::vector<double>& values) {
  return values.empty()?0.0:std::accumulate(values.begin(),values.end(),0.0)/values.size();
}
}

bool ScientificRecorder::sample_node(std::size_t i) const {
  const auto& p=network_->pores[i];
  return p.source_type!=2 && p.volume>0.0 &&
      network_->zero_distance_clusters[network_->cluster_of_pore[i]].boundary_flag==0;
}

double ScientificRecorder::normalized_x(std::size_t i) const {
  return std::clamp((network_->pores[i].x-xmin_)/(xmax_-xmin_),0.0,1.0);
}

double ScientificRecorder::adjacent_constriction_p50(std::size_t i) const {
  return i<adjacent_constriction_p50_.size()?adjacent_constriction_p50_[i]:0.0;
}

ScientificRecorder::ScientificRecorder(const Solver& solver, bool restart)
    : network_(&solver.network()), config_(&solver.config()),
      output_(solver.config().io.output_directory), segments_(solver.config().analysis.x_segments) {
  xmin_=std::numeric_limits<double>::infinity();xmax_=-xmin_;
  for(std::size_t i=0;i<network_->pores.size();++i) if(sample_node(i)){
    xmin_=std::min(xmin_,network_->pores[i].x);xmax_=std::max(xmax_,network_->pores[i].x);
  }
  if(!(xmax_>xmin_)) throw std::runtime_error("sample x extent is empty");
  node_segment_.assign(network_->pores.size(),-1);
  adjacent_constriction_p50_.assign(network_->pores.size(),0.0);
  segment_static_.resize(static_cast<std::size_t>(segments_));
  initial_liquid_=solver.initial_liquid_volume();
  initial_total_=solver.initial_total_moles();
  for(std::size_t i=0;i<network_->pores.size();++i) if(sample_node(i)){
    const int s=std::min(segments_-1,static_cast<int>(normalized_x(i)*segments_));
    node_segment_[i]=s;segment_static_[s].nodes.push_back(i);
    segment_static_[s].volume+=network_->pores[i].volume;sample_volume_+=network_->pores[i].volume;
    const double initial_sg=1.0-static_cast<double>(solver.initial_sw_by_node()[i]);
    initial_gas_volume_+=network_->pores[i].volume*initial_sg;
    initial_free_moles_+=static_cast<double>(solver.initial_free_gas_by_node()[i]);
    if(network_->component_id[i]==network_->main_flow_component_id){
      main_volume_+=network_->pores[i].volume;
      initial_main_sg_+=network_->pores[i].volume*initial_sg;
    }
  }
  initial_main_sg_/=main_volume_;
  std::vector<std::vector<double>> ratios(network_->pores.size());
  for(const auto& e:network_->edges){
    if(e.source_old_throat_id==0 || e.radius<=0) continue;
    for(const auto node:{e.a,e.b}) if(sample_node(node))
      ratios[node].push_back(network_->pores[node].radius/e.radius);
  }
  for(std::size_t i=0;i<ratios.size();++i) adjacent_constriction_p50_[i]=quantile(ratios[i],0.5);
  write_geometry();
  const bool write_headers=!restart||!std::filesystem::exists(output_/"global_history.tsv")||
      std::filesystem::file_size(output_/"global_history.tsv")==0;
  const bool write_event_header=!restart||
      !std::filesystem::exists(output_/"gas_disappearance_events.tsv")||
      std::filesystem::file_size(output_/"gas_disappearance_events.tsv")==0;
  if(write_headers){
    std::ofstream readme(output_/"BALANCE_README.md");
    readme<<"# 守恒记录说明\n\n"
      "`balance_history.tsv` 中入口累计量和出口累计量均以正数记录。组分采用 "
      "`N0 + Nin - Nout + Nledger = Nfree + Ndissolved + E_N`；液体采用 "
      "`Vl0 + Vlin - Vlout + Vledger = Vl(t) + E_V`。`component_balance_error_mol` "
      "和 `liquid_balance_error_m3` 即左侧减右侧。\n\n"
      "气泡缩小时液体饱和体积增加；在不可压缩液体模型中，这一空间必须由边界净流入补充。"
      "气体体积还受绝对压力和毛细压力影响，因此体积减少不能替代自由气体摩尔数减少或界面传质积分。\n\n"
      "PASS使用配置中的 `mass_balance_hard_limit` 和 `volume_balance_hard_limit`，相对误差分别以"
      "初始总组分摩尔数和样品真实控制体积为尺度，并保留求解器逐步守恒硬门限。\n";
  }
  const auto mode=restart?std::ios::app:std::ios::trunc;
  global_.open(output_/"global_history.tsv",mode);
  segment_.open(output_/"segment_history.tsv",mode);
  balance_.open(output_/"balance_history.tsv",mode);
  events_.open(output_/"gas_disappearance_events.tsv",mode);
  if(!global_||!segment_||!balance_||!events_) throw std::runtime_error("cannot open scientific history files");
  for(auto* stream:{&global_,&segment_,&balance_,&events_}) *stream<<std::scientific<<std::setprecision(17);
  if(write_headers){
    global_<<"step\ttime_s\tdt_s\tmain_real_pore_volume_m3\tmain_real_gas_saturation"
      "\tmean_Sw\tmean_Sg\tgas_volume_m3\tfree_gas_mol\tdissolved_gas_mol\ttotal_component_mol"
      "\tcumulative_interphase_transfer_mol\tfree_gas_moles_lost_mol\tgas_volume_retained_fraction"
      "\tfree_gas_moles_retained_fraction\tcumulative_inlet_component_mol\tcumulative_outlet_component_mol"
      "\tinstantaneous_transfer_rate_mol_s\tinlet_liquid_flow_m3_s\toutlet_liquid_flow_m3_s"
      "\tthroughflow_m3_s\tcumulative_inlet_liquid_m3\tcumulative_outlet_liquid_m3\tinjected_pore_volumes"
      "\tactive_bubble_count\tresidual_scaled_inf\trelative_update_norm\tnewton_iterations\tlinear_iterations"
      "\tassembly_seconds\tcsr_seconds\tlinear_upload_seconds\tlinear_setup_seconds\tlinear_solve_seconds"
      "\tlinear_download_seconds\tcpu_linear_relative_residual"
      "\tlinear_symbolic_analysis_seconds\tlinear_numeric_factorization_seconds"
      "\tpardiso_symbolic_analysis_calls\tpardiso_numeric_factorization_calls"
      "\tpardiso_factor_nonzeros\tpardiso_peak_memory_kib"
      "\tretired_bubble_count\tbatched_retired_bubble_count"
      "\tbatched_retired_moles\tbatched_retired_gas_volume_m3"
      "\tbatch_aggregate_free_gas_fraction\tevent_timestep_localized\n";
    segment_<<"step\ttime_s\tsegment_id\tsegment_mean_Sw\tsegment_mean_Sg\tsegment_gas_volume_m3"
      "\tsegment_free_gas_mol\tsegment_dissolved_gas_mol\tsegment_cumulative_interphase_transfer_mol"
      "\tsegment_free_gas_moles_lost_mol\tsegment_active_bubble_count\tsegment_gas_liquid_interfacial_area_m2"
      "\tinterfacial_area_per_gas_volume_1_m\tsegment_mean_C_mol_m3\tsegment_liquid_turnover_rate_1_s"
      "\tleft_advective_component_flux_mol_s\tright_advective_component_flux_mol_s"
      "\tleft_diffusive_component_flux_mol_s\tright_diffusive_component_flux_mol_s"
      "\tsegment_disappearance_events_this_step\tsegment_batched_retirements_this_step\n";
    balance_<<"step\ttime_s\tinitial_total_component_mol\tcurrent_free_gas_mol\tcurrent_dissolved_gas_mol"
      "\tcumulative_component_in_mol\tcumulative_component_out_mol\tactive_set_component_ledger_mol"
      "\tcomponent_balance_error_mol\tcomponent_balance_relative_error\tcomponent_balance_pass"
      "\tinitial_liquid_volume_m3\tcurrent_liquid_volume_m3\tcumulative_liquid_in_m3"
      "\tcumulative_liquid_out_m3\tactive_set_liquid_ledger_m3\tliquid_balance_error_m3"
      "\tliquid_balance_relative_error\tliquid_balance_pass\n";
  }
  if(write_event_header)
    events_<<"step\taccepted_state_time_s\tevent_index_in_step\tevent_node_id\tsegment_id"
      "\tevent_mode\testimated_crossing_time_s\taccepted_dt_s\tcrossing_fraction"
      "\tgas_moles_before\tgas_moles_after_before_retirement"
      "\tgas_moles_transferred_to_dissolved\tgas_volume_before_m3"
      "\tgas_volume_after_before_retirement_m3\tindividual_free_gas_fraction"
      "\tbatch_aggregate_free_gas_fraction\n";
}

ScientificSnapshot ScientificRecorder::snapshot(const Solver& solver) const {
  ScientificSnapshot x; x.main_real_pore_volume_m3=main_volume_;x.sample_pore_volume_m3=sample_volume_;
  double liquid=0,main_gas=0;
  for(std::size_t i=0;i<network_->pores.size();++i) if(sample_node(i)){
    const auto& p=network_->pores[i];const auto& s=solver.state()[i];
    liquid+=p.volume*static_cast<double>(s.sw);
    x.dissolved_gas_mol+=p.volume*static_cast<double>(s.sw*s.c);
    if(s.active_gas){x.gas_volume_m3+=p.volume*static_cast<double>(1-s.sw);
      x.free_gas_mol+=static_cast<double>(free_gas_moles(s.pg,s.sw,p.volume,config_->physics));++x.active_bubbles;}
    if(network_->component_id[i]==network_->main_flow_component_id)
      main_gas+=p.volume*static_cast<double>(1-s.sw);
  }
  x.mean_sw=liquid/sample_volume_;x.mean_sg=x.gas_volume_m3/sample_volume_;
  x.main_real_gas_saturation=main_gas/main_volume_;
  x.injected_pore_volumes=solver.cumulative_liquid_inflow()/sample_volume_;return x;
}

void ScientificRecorder::append(const Solver& solver,const StepReport& r){
  const auto g=snapshot(solver);const double total=g.free_gas_mol+g.dissolved_gas_mol;
  global_<<solver.step_number()<<'\t'<<solver.time()<<'\t'<<r.dt<<'\t'<<g.main_real_pore_volume_m3<<'\t'
    <<g.main_real_gas_saturation<<'\t'<<g.mean_sw<<'\t'<<g.mean_sg<<'\t'<<g.gas_volume_m3<<'\t'
    <<g.free_gas_mol<<'\t'<<g.dissolved_gas_mol<<'\t'<<total<<'\t'<<solver.cumulative_interphase_transfer()<<'\t'
    <<initial_free_moles_-g.free_gas_mol<<'\t'<<g.gas_volume_m3/initial_gas_volume_<<'\t'
    <<g.free_gas_mol/initial_free_moles_<<'\t'<<solver.cumulative_component_inflow()<<'\t'
    <<solver.cumulative_component_outflow()<<'\t'<<r.balance.transfer_rate<<'\t'<<r.balance.inlet_liquid_inflow<<'\t'
    <<r.balance.outlet_liquid_outflow<<'\t'<<r.balance.throughflow_rate<<'\t'<<solver.cumulative_liquid_inflow()<<'\t'
    <<solver.cumulative_liquid_outflow()<<'\t'<<g.injected_pore_volumes<<'\t'<<g.active_bubbles<<'\t'
    <<r.residual_norm<<'\t'<<r.relative_update_norm<<'\t'<<r.newton_iterations<<'\t'<<r.linear_iterations<<'\t'
    <<r.assembly_seconds<<'\t'<<r.csr_seconds<<'\t'<<r.linear_upload_seconds<<'\t'<<r.linear_setup_seconds<<'\t'
    <<r.linear_solve_seconds<<'\t'<<r.linear_download_seconds<<'\t'<<r.linear_residual<<'\t'
    <<r.linear_symbolic_analysis_seconds<<'\t'<<r.linear_numeric_factorization_seconds<<'\t'
    <<r.pardiso_symbolic_analysis_calls<<'\t'<<r.pardiso_numeric_factorization_calls<<'\t'
    <<r.pardiso_factor_nonzeros<<'\t'<<r.pardiso_peak_memory_kib<<'\t'
    <<r.retired_bubble_count<<'\t'<<r.batched_retired_bubble_count<<'\t'
    <<r.batched_retired_moles<<'\t'<<r.batched_retired_gas_volume<<'\t'
    <<r.batch_aggregate_free_gas_fraction<<'\t'<<r.event_timestep_localized<<'\n';global_.flush();

  for(std::size_t event_index=0;event_index<r.gas_disappearance_event_records.size();++event_index){
    const auto& event=r.gas_disappearance_event_records[event_index];
    const int segment_id=event.node_index<node_segment_.size()?node_segment_[event.node_index]+1:0;
    events_<<solver.step_number()<<'\t'<<solver.time()<<'\t'<<event_index+1<<'\t'
      <<event.node_id<<'\t'<<segment_id<<'\t'
      <<(event.batched?"batched_microbubble":"exact_localized")<<'\t'
      <<event.estimated_crossing_time<<'\t'<<r.dt<<'\t'<<event.crossing_fraction<<'\t'
      <<event.gas_moles_before<<'\t'<<event.gas_moles_after_before_retirement<<'\t'
      <<event.gas_moles_transferred_to_dissolved<<'\t'<<event.gas_volume_before<<'\t'
      <<event.gas_volume_after_before_retirement<<'\t'<<event.individual_free_gas_fraction<<'\t'
      <<event.batch_aggregate_free_gas_fraction<<'\n';
  }
  events_.flush();

  std::vector<double> iface_q(segments_+1),iface_a(segments_+1),iface_d(segments_+1);
  for(std::size_t k=0;k<network_->edges.size();++k){const auto&e=network_->edges[k];if(e.length==0)continue;
    int sa=node_segment_[e.a],sb=node_segment_[e.b];double sign=network_->pores[e.a].x<=network_->pores[e.b].x?1:-1;
    int lo=sa,hi=sb;if(sa<0)lo=network_->pores[e.a].x<xmin_?-1:segments_;
    if(sb<0)hi=network_->pores[e.b].x<xmin_?-1:segments_;
    if(lo>hi){std::swap(lo,hi);sign=-sign;} for(int f=lo+1;f<=hi && f<=segments_;++f)if(f>=0){
      iface_q[f]+=sign*r.final_assembly.edge_fluxes[k].q;iface_a[f]+=sign*r.final_assembly.edge_fluxes[k].advective;
      iface_d[f]+=sign*r.final_assembly.edge_fluxes[k].diffusive;}}
  for(int sid=0;sid<segments_;++sid){double liq=0,gas=0,free=0,diss=0,transfer=0,initial_free=0,area=0,cweighted=0;std::size_t active=0,event_count=0,batch_count=0;
    for(const auto& event:r.gas_disappearance_event_records)if(event.node_index<node_segment_.size()&&node_segment_[event.node_index]==sid){++event_count;if(event.batched)++batch_count;}
    for(auto i:segment_static_[sid].nodes){const auto&p=network_->pores[i];const auto&s=solver.state()[i];liq+=p.volume*static_cast<double>(s.sw);
      diss+=p.volume*static_cast<double>(s.sw*s.c);cweighted+=p.volume*static_cast<double>(s.sw*s.c);
      transfer+=static_cast<double>(solver.cumulative_transfer_by_node()[i]);initial_free+=static_cast<double>(solver.initial_free_gas_by_node()[i]);
      if(s.active_gas){gas+=p.volume*static_cast<double>(1-s.sw);free+=static_cast<double>(free_gas_moles(s.pg,s.sw,p.volume,config_->physics));
        area+=static_cast<double>(interface_area(s.sw,p.volume));++active;}}
    const double vol=segment_static_[sid].volume;const double turnover=vol>0?0.5*(std::abs(iface_q[sid])+std::abs(iface_q[sid+1]))/vol:0;
    segment_<<solver.step_number()<<'\t'<<solver.time()<<'\t'<<sid+1<<'\t'<<(vol?liq/vol:0)<<'\t'<<(vol?gas/vol:0)<<'\t'
      <<gas<<'\t'<<free<<'\t'<<diss<<'\t'<<transfer<<'\t'<<initial_free-free<<'\t'<<active<<'\t'<<area<<'\t'
      <<(gas?area/gas:0)<<'\t'<<(liq?cweighted/liq:0)<<'\t'<<turnover<<'\t'<<iface_a[sid]<<'\t'<<iface_a[sid+1]
      <<'\t'<<iface_d[sid]<<'\t'<<iface_d[sid+1]<<'\t'<<event_count<<'\t'<<batch_count<<'\n';}segment_.flush();
  const double en=initial_total_+solver.cumulative_component_inflow()-solver.cumulative_component_outflow()+solver.cumulative_mole_ledger()-total;
  const double ev=initial_liquid_+solver.cumulative_liquid_inflow()-solver.cumulative_liquid_outflow()+solver.cumulative_volume_ledger()-g.mean_sw*sample_volume_;
  const double ern=std::abs(en)/std::max(initial_total_,config_->scaling.mole_total),erv=std::abs(ev)/std::max(sample_volume_,config_->scaling.volume_total);
  balance_<<solver.step_number()<<'\t'<<solver.time()<<'\t'<<initial_total_<<'\t'<<g.free_gas_mol<<'\t'<<g.dissolved_gas_mol<<'\t'
    <<solver.cumulative_component_inflow()<<'\t'<<solver.cumulative_component_outflow()<<'\t'<<solver.cumulative_mole_ledger()<<'\t'
    <<en<<'\t'<<ern<<'\t'<<(ern<=config_->acceptance.mass_balance_hard_limit)<<'\t'<<initial_liquid_<<'\t'
    <<g.mean_sw*sample_volume_<<'\t'<<solver.cumulative_liquid_inflow()<<'\t'<<solver.cumulative_liquid_outflow()<<'\t'
    <<solver.cumulative_volume_ledger()<<'\t'<<ev<<'\t'<<erv<<'\t'<<(erv<=config_->acceptance.volume_balance_hard_limit)<<'\n';balance_.flush();
}

void ScientificRecorder::write_geometry() const {
  std::ofstream out(output_/"segment_geometry.tsv");out<<std::scientific<<std::setprecision(17);
  out<<"segment_id\tx_start_norm\tx_end_norm\treal_pore_count\tinternal_throat_count\tpore_volume_m3"
    "\tpore_radius_mean_m\tpore_radius_p50_m\tpore_radius_p10_m\tpore_radius_p90_m\tthroat_radius_mean_m"
    "\tthroat_radius_p50_m\tthroat_radius_p10_m\tthroat_radius_p90_m\tmean_coordination_number"
    "\tcross_segment_throat_count\thydraulic_conductance_sum_m3_Pa_s\tconstriction_ratio_p50\tconstriction_ratio_p10\tconstriction_ratio_p90\n";
  struct Throat{std::vector<std::size_t> pores;double radius=0,length=0;};std::map<std::size_t,Throat> throats;
  for(const auto&e:network_->edges)if(e.source_old_throat_id){auto&t=throats[e.source_old_throat_id];t.radius=e.radius;t.length+=e.length;
    for(auto i:{e.a,e.b})if(network_->pores[i].source_type==0)t.pores.push_back(i);}
  std::vector<std::vector<double>> tr(segments_),cr(segments_);std::vector<int> tn(segments_),cross(segments_);std::vector<double> conduct(segments_);
  std::ofstream defects(output_/"constriction_geometry_audit.tsv");defects<<"source_old_throat_id\tpore_id_1\tpore_id_2\tthroat_radius_m\tminimum_pore_radius_m\n";
  for(const auto&item:throats){auto p=item.second.pores;std::sort(p.begin(),p.end());p.erase(std::unique(p.begin(),p.end()),p.end());if(p.size()!=2)continue;
    const auto&t=item.second;const double mid=0.5*(network_->pores[p[0]].x+network_->pores[p[1]].x);const int s=std::min(segments_-1,std::max(0,static_cast<int>((mid-xmin_)/(xmax_-xmin_)*segments_)));
    ++tn[s];if(node_segment_[p[0]]!=node_segment_[p[1]])++cross[s];tr[s].push_back(t.radius);const double minp=std::min(network_->pores[p[0]].radius,network_->pores[p[1]].radius);cr[s].push_back(minp/t.radius);
    if(t.length>0)conduct[s]+=pi*std::pow(t.radius,4)/(8*config_->physics.viscosity*t.length);
    if(t.radius>minp)defects<<item.first<<'\t'<<network_->pores[p[0]].id<<'\t'<<network_->pores[p[1]].id<<'\t'<<t.radius<<'\t'<<minp<<'\n';}
  for(int s=0;s<segments_;++s){std::vector<double> pr,coord;for(auto i:segment_static_[s].nodes){pr.push_back(network_->pores[i].radius);coord.push_back(network_->adjacency[i].size());}
    out<<s+1<<'\t'<<double(s)/segments_<<'\t'<<double(s+1)/segments_<<'\t'<<segment_static_[s].nodes.size()<<'\t'<<tn[s]<<'\t'<<segment_static_[s].volume<<'\t'
      <<mean(pr)<<'\t'<<quantile(pr,.5)<<'\t'<<quantile(pr,.1)<<'\t'<<quantile(pr,.9)<<'\t'<<mean(tr[s])<<'\t'<<quantile(tr[s],.5)<<'\t'
      <<quantile(tr[s],.1)<<'\t'<<quantile(tr[s],.9)<<'\t'<<mean(coord)<<'\t'<<cross[s]<<'\t'<<conduct[s]<<'\t'<<quantile(cr[s],.5)<<'\t'<<quantile(cr[s],.1)<<'\t'<<quantile(cr[s],.9)<<'\n';}
}

}  // namespace bubble
