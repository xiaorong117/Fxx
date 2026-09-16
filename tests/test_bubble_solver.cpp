#include "bubble/model.hpp"

#include <Eigen/SparseLU>
#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

bubble::Config base_config() {
  bubble::Config c;
  c.case_info.kind = "deterministic_test";
  c.geometry.input_basis = "physical";
  c.geometry.input_length_unit = "m";
  c.geometry.voxel_size_um = 1.0;
  c.geometry.image_dimensions_voxels = {{1,1,1}};
  c.geometry.physical_domain_um = {{1.0,1.0,1.0}};
  c.virtual_boundary.enabled = false;
  c.virtual_boundary.boundary_axis = "x";
  c.virtual_boundary.extension_voxels = 1.0;
  c.virtual_boundary.voxel_size_um = 1.0;
  c.physics.temperature = 298.15;
  c.physics.viscosity = 1e-3;
  c.physics.molecular_diffusivity = 2e-9;
  c.physics.henry_cp = 1e-5;
  c.physics.mass_transfer_coefficient = 1e-5;
  c.physics.surface_tension = 0.072;
  c.physics.contact_angle = 0.0;
  c.physics.molar_mass = 0.028;
  c.physics.liquid_density = 997.0;
  c.physics.compressibility = 1.0;
  c.boundary.pl_in = 102000.0;
  c.boundary.pl_out = 101000.0;
  c.boundary.species_mode = "fixed_reservoir_concentration";
  c.initial.concentration = 0.1;
  c.time.end = 1e-3;
  c.time.initial_dt = 1e-4;
  c.time.min_dt = 1e-12;
  c.time.max_dt = 1e-3;
  c.scaling.pressure = 1e5;
  c.scaling.concentration = 1.0;
  c.scaling.volume_rate = 1e-15;
  c.scaling.molar_rate = 1e-14;
  c.scaling.mole_total = 1e-24;
  c.scaling.volume_total = 1e-18;
  c.active_set.ng_off = 1e-30;
  c.linear.backend = "eigen_sparselu";
  c.acceptance.mass_balance_hard_limit = 1e-10;
  c.acceptance.volume_balance_hard_limit = 1e-10;
  return c;
}

bubble::Pore pore(std::size_t id, double x, int flag, double sw) {
  bubble::Pore p;
  p.id=id;p.source_type=0;p.source_id=id;p.x=x;p.radius=5e-6;p.area=bubble::pi*25e-12;
  p.shape_factor=0.05;p.volume=5e-15;p.boundary_flag=flag;p.input_sw=sw;
  return p;
}

bubble::Edge edge(std::size_t id, std::size_t a, std::size_t b) {
  bubble::Edge e;
  e.id=id;e.a=a;e.b=b;e.source_old_throat_id=id;e.side=(id%2)?1:2;e.length=2e-5;
  e.input_sw=1.0;e.radius=2e-6;e.area=bubble::pi*4e-12;e.shape_factor=0.05;e.source_volume=1e-16;
  return e;
}

bubble::Network chain_network() {
  bubble::Network n;
  n.pores={pore(1,0.0,1,1.0),pore(2,1.0,0,0.65),pore(3,2.0,0,1.0),pore(4,3.0,2,1.0)};
  n.edges={edge(1,0,1),edge(2,1,2),edge(3,2,3)};
  n.adjacency.resize(n.pores.size());
  for(std::size_t k=0;k<n.edges.size();++k){n.adjacency[n.edges[k].a].push_back(k);n.adjacency[n.edges[k].b].push_back(k);}
  return n;
}

bubble::Network single_bubble_network(double sw=0.65) {
  bubble::Network n;
  n.pores={pore(1,0.0,0,sw),pore(2,1.0,2,1.0)};
  n.edges={edge(1,0,1)};
  n.adjacency={{0},{0}};
  return n;
}

bubble::Network two_bubble_network(double sw=0.65) {
  bubble::Network n;
  n.pores={pore(1,0.0,0,sw),pore(2,0.5,0,sw),pore(3,1.0,2,1.0)};
  n.edges={edge(1,0,2),edge(2,1,2)};
  n.adjacency.resize(n.pores.size());
  for(std::size_t k=0;k<n.edges.size();++k){
    n.adjacency[n.edges[k].a].push_back(k);n.adjacency[n.edges[k].b].push_back(k);
  }
  return n;
}

bubble::Network liquid_five_node_network() {
  bubble::Network n;
  n.pores={pore(1,0.0,1,1.0),pore(2,1.0,0,1.0),pore(3,2.0,0,1.0),
           pore(4,3.0,0,1.0),pore(5,4.0,2,1.0)};
  n.edges={edge(1,0,1),edge(2,1,2),edge(3,2,3),edge(4,3,4)};
  n.adjacency.resize(n.pores.size());
  for(std::size_t k=0;k<n.edges.size();++k){
    n.adjacency[n.edges[k].a].push_back(k);n.adjacency[n.edges[k].b].push_back(k);
  }
  return n;
}

bubble::Network closed_liquid_network(std::size_t nodes=2) {
  bubble::Network n;
  for(std::size_t i=0;i<nodes;++i)n.pores.push_back(pore(i+1,static_cast<double>(i),0,1.0));
  for(std::size_t i=1;i<nodes;++i)n.edges.push_back(edge(i,i-1,i));
  n.adjacency.resize(nodes);
  for(std::size_t k=0;k<n.edges.size();++k){
    n.adjacency[n.edges[k].a].push_back(k);n.adjacency[n.edges[k].b].push_back(k);
  }
  return n;
}

bubble::Network make_network(std::vector<bubble::Pore> pores,
                             std::vector<bubble::Edge> edges){
  bubble::Network n;
  n.pores=std::move(pores);n.edges=std::move(edges);n.adjacency.resize(n.pores.size());
  for(std::size_t k=0;k<n.edges.size();++k){
    n.adjacency[n.edges[k].a].push_back(k);n.adjacency[n.edges[k].b].push_back(k);
  }
  return n;
}

bubble::Edge zero_edge(std::size_t id,std::size_t a,std::size_t b){
  auto result=edge(id,a,b);result.length=0.0;return result;
}

bubble::Pore virtual_pore(std::size_t id,std::size_t real_id,double x,int flag){
  auto result=pore(id,x,flag,1.0);result.source_type=2;result.source_id=real_id;
  return result;
}

bubble::Edge boundary_edge(std::size_t id,std::size_t a,std::size_t b,int flag){
  auto result=edge(id,a,b);result.source_old_throat_id=0;result.side=flag;
  result.length=7.5e-6;result.radius=5e-6;result.area=bubble::pi*25e-12;
  result.source_volume=0.0;return result;
}

void set_dof(std::vector<bubble::NodeState>& state, const bubble::DofMap& map,
             int column, double delta) {
  bool found=false;
  for(std::size_t i=0;i<state.size();++i){
    const auto& d=map.node[i];
    if(d.pl==column){state[i].pl+=delta;found=true;continue;}
    if(d.c==column){state[i].c+=delta;found=true;continue;}
    if(d.pg==column){state[i].pg+=delta;found=true;continue;}
    if(d.sw==column){state[i].sw+=delta;found=true;continue;}
  }
  if(!found)throw std::runtime_error("unknown DOF column");
}

template<class Function>
void require_throws(Function function,const std::string& expected){
  try{function();}
  catch(const std::exception& error){
    require(std::string(error.what()).find(expected)!=std::string::npos,
            "unexpected failure text: "+std::string(error.what()));
    return;
  }
  throw std::runtime_error("expected failure was not raised: "+expected);
}

void audit_jacobian(const bubble::Solver& solver,
                    const std::vector<bubble::NodeState>& state,double tolerance=1e-6){
  const auto map=solver.make_dof_map(state);
  std::vector<bool> upwind(solver.network().edges.size(),true);
  std::vector<bool> transfer(state.size(),false);
  for(std::size_t i=0;i<state.size();++i)
    transfer[i]=state[i].active_gas &&
        solver.config().physics.henry_cp*state[i].pg-state[i].c>
            solver.config().active_set.transfer_switch_tolerance;
  const auto ad=solver.assemble(state,state,1e-4,map,upwind,transfer,true);
  for(int col=0;col<map.size;++col){
    const auto column=static_cast<std::size_t>(col);
    const double h=2e-6*map.unknown_scale[column];
    auto plus=state,minus=state;
    set_dof(plus,map,col,h);set_dof(minus,map,col,-h);
    const auto rp=solver.assemble(plus,state,1e-4,map,upwind,transfer,false).residual;
    const auto rm=solver.assemble(minus,state,1e-4,map,upwind,transfer,false).residual;
    for(int row=0;row<map.size;++row){
      const auto residual_row=static_cast<std::size_t>(row);
      const double fd=(rp[row]-rm[row])/(2*h);
      const double av=ad.jacobian.coeff(row,col);
      const double scaled_error=std::abs(av-fd)*map.unknown_scale[column]/
                                map.residual_scale[residual_row];
      const double scaled_size=std::max({std::abs(av)*map.unknown_scale[column]/
                                             map.residual_scale[residual_row],
                                         std::abs(fd)*map.unknown_scale[column]/
                                             map.residual_scale[residual_row],1e-12});
      require(scaled_error<=1e-9 || scaled_error/scaled_size<=tolerance,
              "cluster physical AD Jacobian finite-difference mismatch");
    }
  }
}

void test_t1_t2_t4_fluxes() {
  auto c=base_config();
  auto e=edge(1,0,1);
  const double pin=102000.0,pout=101000.0;
  const auto f=bubble::evaluate_edge_flux(e,pin,0.0,1.0,pout,0.0,1.0,true,true,c.physics);
  const double exact=bubble::pi*std::pow(e.radius,4)*(pin-pout)/(8*c.physics.viscosity*e.length);
  require(std::abs(f.q-exact)/exact<=1e-12,"T1 circular-tube flow mismatch");

  const auto d=bubble::evaluate_edge_flux(e,pin,2.0,1.0,pin,0.5,1.0,true,true,c.physics);
  const double de=c.physics.molecular_diffusivity*e.area*(2.0-0.5)/e.length;
  require(d.q==0.0 && std::abs(d.diffusive-de)<=1e-12*std::abs(de),"T2 diffusion flux mismatch");
  require(std::abs(d.diffusive+(-d.diffusive))<=1e-30,"T2 edge contributions do not cancel");
  {
    const double volume=5e-15,dt=2e-4,k=c.physics.molecular_diffusivity*e.area/e.length;
    Eigen::SparseMatrix<double> matrix(2,2);
    std::vector<Eigen::Triplet<double>> triplets{{0,0,volume/dt+k},{0,1,-k},
                                                 {1,0,-k},{1,1,volume/dt+k}};
    matrix.setFromTriplets(triplets.begin(),triplets.end());
    const Eigen::Vector2d old_c(2.0,0.5);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solve;solve.compute(matrix);
    const Eigen::Vector2d new_c=solve.solve((volume/dt)*old_c);
    require(solve.info()==Eigen::Success,"T2 implicit diffusion solve failed");
    const double old_moles=volume*old_c.sum(),new_moles=volume*new_c.sum();
    require(std::abs(new_moles-old_moles)/old_moles<=1e-12,
            "T2 closed dissolved inventory is not conserved");
  }
  {
    auto closed_config=c;
    closed_config.initial.concentration=0.0;
    bubble::Solver closed(closed_liquid_network(),closed_config);closed.initialize();
    closed.mutable_state_for_test()[0].c=2.0;
    closed.mutable_state_for_test()[1].c=0.5;
    const auto map=closed.make_dof_map(closed.state());
    const auto gauges=std::count_if(map.node.begin(),map.node.end(),
                                    [](const auto& x){return x.pressure_gauge;});
    require(gauges==1,"T2 closed liquid component did not receive one pressure gauge");
    const double before=static_cast<double>(closed.network().pores[0].volume*closed.state()[0].c+
                        closed.network().pores[1].volume*closed.state()[1].c);
    const auto report=closed.attempt_step(2e-4);
    require(report.accepted,"T2 production closed diffusion step failed: "+report.failure);
    const double after=static_cast<double>(closed.network().pores[0].volume*closed.state()[0].c+
                       closed.network().pores[1].volume*closed.state()[1].c);
    require(std::abs(after-before)/before<=1e-12 &&
            report.balance.boundary_molar_outflow==0.0 &&
            report.balance.epsilon_moles_equation<=1e-12,
            "T2 production closed diffusion conservation failed");

    bubble::Solver isolated(closed_liquid_network(1),closed_config);isolated.initialize();
    const auto isolated_step=isolated.attempt_step(2e-4);
    require(isolated_step.accepted && isolated.make_dof_map(isolated.state()).node[0].pressure_gauge,
            "isolated liquid control pressure gauge failed");
  }

  const auto forward=bubble::evaluate_edge_flux(e,pin,3.0,1.0,pout,7.0,1.0,true,false,c.physics);
  const auto reverse=bubble::evaluate_edge_flux(e,pout,3.0,1.0,pin,7.0,1.0,false,false,c.physics);
  require(std::abs(forward.advective-forward.q*3.0)<=1e-28,"T4 forward upwind wrong");
  require(std::abs(reverse.advective-reverse.q*7.0)<=1e-28,"T4 reverse upwind wrong");
  auto high_pe_step=[](double q,double left_boundary,double right_boundary){
    constexpr int cells=3;const double volume=1.0,dt=0.1;
    Eigen::Matrix3d a=Eigen::Matrix3d::Identity()*(volume/dt);
    Eigen::Vector3d rhs=Eigen::Vector3d::Constant(0.05*volume/dt);
    // Four oriented faces: left reservoir -> cells 0 -> 1 -> 2 -> right reservoir.
    for(int face=0;face<4;++face){
      const int left=face-1,right=face;
      const bool upwind_left=q>=0.0;
      const int upstream=upwind_left?left:right;
      const double boundary_value=upstream<0?left_boundary:(upstream>=cells?right_boundary:0.0);
      if(left>=0 && left<cells){
        if(upstream>=0 && upstream<cells)a(left,upstream)+=q;
        else rhs[left]-=q*boundary_value;
      }
      if(right>=0 && right<cells){
        if(upstream>=0 && upstream<cells)a(right,upstream)-=q;
        else rhs[right]+=q*boundary_value;
      }
    }
    Eigen::Vector3d result=a.fullPivLu().solve(rhs);
    return result;
  };
  const auto positive=high_pe_step(100.0,1.0,0.0);
  const auto negative=high_pe_step(-100.0,0.0,1.0);
  require(positive.minCoeff()>=-1e-12 && negative.minCoeff()>=-1e-12,
          "T4 high-Peclet upwind chain produced negative concentration");

  auto run_actual_chain=[](bool forward_flow){
    auto config=base_config();
    config.physics.molecular_diffusivity=1e-20;
    config.initial.concentration=0.05;
    config.boundary.c_in=1.0;
    config.boundary.c_backflow=0.8;
    if(forward_flow){config.boundary.pl_in=102000.0;config.boundary.pl_out=101000.0;}
    else{config.boundary.pl_in=101000.0;config.boundary.pl_out=102000.0;}
    bubble::Solver solver(liquid_five_node_network(),config);solver.initialize();
    const auto report=solver.attempt_step(1e-4);
    require(report.accepted,std::string("T4 actual chain step failed: ")+report.failure);
    for(std::size_t i=1;i<=3;++i)
      require(solver.state()[i].c>=-1e-12*config.scaling.concentration,
              "T4 actual chain produced negative concentration");
    for(const auto& flux:report.final_assembly.edge_fluxes){
      require(flux.upwind_a==forward_flow,"T4 final frozen/up-to-date upwind direction mismatch");
    }
    const auto map=solver.make_dof_map(solver.state());
    std::vector<bool> upwind;upwind.reserve(report.final_assembly.edge_fluxes.size());
    for(const auto& flux:report.final_assembly.edge_fluxes)upwind.push_back(flux.upwind_a);
    const std::vector<bool> no_transfer(solver.state().size(),false);
    const auto zero_accumulation=solver.assemble(solver.state(),solver.state(),1e-4,map,
                                                 upwind,no_transfer,false);
    double sum_rl=0.0,sum_rd=0.0;
    for(std::size_t i=0;i<map.node.size();++i){
      if(map.node[i].rl>=0)sum_rl+=zero_accumulation.residual[map.node[i].rl];
      if(map.node[i].rd>=0)sum_rd+=zero_accumulation.residual[map.node[i].rd];
    }
    const auto& left=zero_accumulation.edge_fluxes.front();
    const auto& right=zero_accumulation.edge_fluxes.back();
    const double boundary_q=-left.q+right.q;
    const double boundary_f=-(left.advective+left.diffusive)+
                            right.advective+right.diffusive;
    require(std::abs(sum_rl-boundary_q)<=1e-15*std::max(std::abs(boundary_q),1e-300) &&
            std::abs(sum_rd-boundary_f)<=1e-15*std::max(std::abs(boundary_f),1e-300),
            "T4 internal numerical flux did not cancel from global residual sum");
    if(!forward_flow)
      require(solver.state()[3].c>config.initial.concentration,
              "T4 outlet backflow concentration did not enter reverse-flow chain");
  };
  run_actual_chain(true);
  run_actual_chain(false);
  std::cout<<"PASS T1 T2 T4 flux tests\n";
}

void test_jacobian_and_dof_map() {
  auto c=base_config();
  auto n=chain_network();
  bubble::Solver solver(n,c);
  solver.initialize();
  auto state=solver.state();
  state[1].pl=101700.0;state[1].c=0.2;state[1].pg=state[1].pl+bubble::capillary_pressure(state[1].sw,n.pores[1].radius,c.physics);
  state[2].pl=101300.0;state[2].c=0.4;
  const auto map=solver.make_dof_map(state);
  require(map.size==6 && map.node[0].pl<0 && map.node[1].pg>=0 && map.node[2].pg<0,
          "T5 dynamic DOF map mismatch size="+std::to_string(map.size)+" pl0="+std::to_string(map.node[0].pl)+" pg1="+std::to_string(map.node[1].pg)+" pg2="+std::to_string(map.node[2].pg));
  const std::vector<bool> upwind{true,true,true};
  std::vector<bool> transfer(state.size(),false);transfer[1]=true;
  const auto ad=solver.assemble(state,state,1e-4,map,upwind,transfer,true);
  double worst=0.0;
  for(int col=0;col<map.size;++col){
    const double h=2e-6*map.unknown_scale[static_cast<std::size_t>(col)];
    const auto column=static_cast<std::size_t>(col);
    auto plus=state,minus=state;set_dof(plus,map,col,h);set_dof(minus,map,col,-h);
    const auto rp=solver.assemble(plus,state,1e-4,map,upwind,transfer,false).residual;
    const auto rm=solver.assemble(minus,state,1e-4,map,upwind,transfer,false).residual;
    for(int row=0;row<map.size;++row){
      const double fd=(rp[row]-rm[row])/(2*h);
      const double av=ad.jacobian.coeff(row,col);
      const auto residual_row=static_cast<std::size_t>(row);
      const double scaled_error=std::abs(av-fd)*map.unknown_scale[column]/map.residual_scale[residual_row];
      const double scaled_size=std::max({std::abs(av)*map.unknown_scale[column]/map.residual_scale[residual_row],
                                         std::abs(fd)*map.unknown_scale[column]/map.residual_scale[residual_row],1e-12});
      worst=std::max(worst,scaled_error/scaled_size);
      if (!(scaled_error<=1e-9 || scaled_error/scaled_size<=1e-6)) {
        std::cerr << std::scientific << "mismatch row=" << row << " col=" << col
                  << " ad=" << av << " fd=" << fd << " scaled_abs=" << scaled_error
                  << " relative=" << scaled_error/scaled_size << '\n';
        throw std::runtime_error("physical AD Jacobian finite-difference mismatch");
      }
    }
  }
  auto supersaturated=state;
  supersaturated[1].c=c.physics.henry_cp*supersaturated[1].pg+0.25;
  const std::vector<bool> transfer_off(state.size(),false);
  const auto inactive_ad=solver.assemble(supersaturated,supersaturated,1e-4,map,
                                         upwind,transfer_off,true);
  require(inactive_ad.transfer[1]==0.0,
          "supersaturated transfer-off branch produced dissolution");
  for(int col=0;col<map.size;++col){
    const double h=2e-6*map.unknown_scale[static_cast<std::size_t>(col)];
    const auto column=static_cast<std::size_t>(col);
    auto plus=supersaturated,minus=supersaturated;
    set_dof(plus,map,col,h);set_dof(minus,map,col,-h);
    const auto rp=solver.assemble(plus,supersaturated,1e-4,map,upwind,
                                  transfer_off,false).residual;
    const auto rm=solver.assemble(minus,supersaturated,1e-4,map,upwind,
                                  transfer_off,false).residual;
    for(int row=0;row<map.size;++row){
      const double fd=(rp[row]-rm[row])/(2*h);
      const double av=inactive_ad.jacobian.coeff(row,col);
      const auto residual_row=static_cast<std::size_t>(row);
      const double scaled_error=std::abs(av-fd)*map.unknown_scale[column]/map.residual_scale[residual_row];
      const double scaled_size=std::max({std::abs(av)*map.unknown_scale[column]/map.residual_scale[residual_row],
                                         std::abs(fd)*map.unknown_scale[column]/map.residual_scale[residual_row],1e-12});
      worst=std::max(worst,scaled_error/scaled_size);
      if (!(scaled_error<=1e-9 || scaled_error/scaled_size<=1e-6))
        throw std::runtime_error("transfer-off AD Jacobian finite-difference mismatch");
    }
  }
  std::cout<<"PASS physical Jacobian worst_relative="<<worst<<"\n";
}

void test_monolithic_steps() {
  auto c=base_config();
  c.initial.concentration=0.0;
  bubble::Solver one(single_bubble_network(),c);
  one.initialize();
  one.mutable_state_for_test()[0].pl=c.boundary.pl_out;
  one.mutable_state_for_test()[0].pg=c.boundary.pl_out+
      bubble::capillary_pressure(one.state()[0].sw,one.network().pores[0].radius,c.physics);
  const double ng0=static_cast<double>(bubble::free_gas_moles(one.state()[0].pg,one.state()[0].sw,
                                           one.network().pores[0].volume,c.physics));
  const auto r=one.attempt_step(1e-4);
  require(r.accepted,"T3 monolithic single-bubble step failed: "+r.failure);
  const double ng1=static_cast<double>(bubble::free_gas_moles(one.state()[0].pg,one.state()[0].sw,
                                           one.network().pores[0].volume,c.physics));
  require(ng1<=ng0 && one.state()[0].sw>=0.65 && r.final_assembly.transfer[0]>=0.0,
          "T3 dissolution monotonicity failed");
  require(r.residual_norm<=c.nonlinear.residual_tolerance &&
          r.balance.epsilon_moles_equation<=1e-10 && r.balance.epsilon_volume_equation<=1e-10,
          "T3 residual/conservation threshold failed");
  std::cout<<"PASS T3 single bubble ng_ratio="<<ng1/ng0<<" newton="<<r.newton_iterations<<"\n";

  bubble::Solver sparse(chain_network(),c);
  sparse.initialize();
  auto ci=c;ci.linear.backend="eigen_bicgstab_ilut";ci.linear.ilut_drop_tolerance=1e-12;
  ci.linear.ilut_fill_factor=20;
  bubble::Solver iterative(chain_network(),ci);
  iterative.initialize();
  const auto rs=sparse.attempt_step(1e-4);
  const auto ri=iterative.attempt_step(1e-4);
  require(rs.accepted,"T5 SparseLU solve failed: "+rs.failure);
  require(ri.accepted,"T5 BiCGSTAB+ILUT solve failed: "+ri.failure);
  double difference=0.0;
  const auto map=sparse.make_dof_map(sparse.state());
  for(std::size_t i=0;i<sparse.state().size();++i){
    const auto& a=sparse.state()[i];const auto& b=iterative.state()[i];const auto& d=map.node[i];
    if(d.pl>=0)difference=std::max(difference,static_cast<double>(std::abs(a.pl-b.pl)/c.scaling.pressure));
    if(d.c>=0)difference=std::max(difference,static_cast<double>(std::abs(a.c-b.c)/c.scaling.concentration));
    if(d.pg>=0)difference=std::max(difference,static_cast<double>(std::abs(a.pg-b.pg)/c.scaling.pressure));
    if(d.sw>=0)difference=std::max(difference,static_cast<double>(std::abs(a.sw-b.sw)));
  }
  require(difference<=1e-8,"T5 SparseLU/iterative scaled solution mismatch");
  std::cout<<"PASS T5 four-field coupling backend_difference="<<difference<<"\n";
}

struct Outcome {
  double free_moles=0.0,dissolved_moles=0.0,mean_sw=0.0,outlet_moles=0.0;
  std::vector<bubble::NodeState> state;
};

Outcome integrate_single(double dt,double end,bubble::Config c) {
  c.initial.concentration=0.0;
  bubble::Solver solver(single_bubble_network(),c);
  solver.initialize();
  solver.mutable_state_for_test()[0].pl=c.boundary.pl_out;
  solver.mutable_state_for_test()[0].pg=c.boundary.pl_out+
      bubble::capillary_pressure(solver.state()[0].sw,solver.network().pores[0].radius,c.physics);
  Outcome out;
  while(solver.time()<end-1e-15*end){
    const auto r=solver.attempt_step(std::min(dt,end-solver.time()));
    if(!r.accepted){
      std::cerr<<std::scientific<<"integration failure t="<<solver.time()<<" dt="<<dt
               <<" residual="<<r.residual_norm<<" update="<<r.update_norm<<" linear="
               <<r.linear_residual<<" reason="<<r.failure<<'\n';
      throw std::runtime_error("fixed-step integration failed");
    }
    out.outlet_moles+=r.dt*r.balance.boundary_molar_outflow;
  }
  out.state=solver.state();
  double sw_sum=0.0;std::size_t internal_count=0;
  for(std::size_t i=0;i<solver.state().size();++i){
    if(solver.network().pores[i].boundary_flag!=0)continue;
    const auto& s=solver.state()[i];const auto& p=solver.network().pores[i];
    out.dissolved_moles+=static_cast<double>(p.volume*s.sw*s.c);
    if(s.active_gas)out.free_moles+=static_cast<double>(bubble::free_gas_moles(s.pg,s.sw,p.volume,c.physics));
    sw_sum+=static_cast<double>(s.sw);++internal_count;
  }
  out.mean_sw=sw_sum/static_cast<double>(internal_count);
  return out;
}

Outcome integrate_chain(double dt,double end,bubble::Config c) {
  c.initial.concentration=0.0;
  bubble::Solver solver(chain_network(),c);solver.initialize();
  Outcome out;
  while(solver.time()<end-1e-15*end){
    const auto r=solver.attempt_step(std::min(dt,end-solver.time()));
    require(r.accepted,"chain fixed-step integration failed: "+r.failure);
    for(std::size_t k=0;k<solver.network().edges.size();++k){
      const auto& e=solver.network().edges[k];
      const bool outlet_a=solver.network().pores[e.a].boundary_flag==c.boundary.outlet_flag;
      const bool outlet_b=solver.network().pores[e.b].boundary_flag==c.boundary.outlet_flag;
      if(!(outlet_a||outlet_b))continue;
      const double sign=outlet_b?1.0:-1.0;
      out.outlet_moles+=r.dt*sign*(r.final_assembly.edge_fluxes[k].advective+
                                   r.final_assembly.edge_fluxes[k].diffusive);
    }
  }
  out.state=solver.state();double sw_sum=0.0;std::size_t count=0;
  for(std::size_t i=0;i<solver.state().size();++i){
    if(solver.network().pores[i].boundary_flag!=0)continue;
    const auto& s=solver.state()[i];const auto& p=solver.network().pores[i];
    out.dissolved_moles+=static_cast<double>(p.volume*s.sw*s.c);
    if(s.active_gas)out.free_moles+=static_cast<double>(bubble::free_gas_moles(s.pg,s.sw,p.volume,c.physics));
    sw_sum+=static_cast<double>(s.sw);++count;
  }
  out.mean_sw=sw_sum/static_cast<double>(count);
  return out;
}

void test_reference_and_dt_convergence() {
  auto c=base_config();
  const auto coarse=integrate_single(5e-5,1e-3,c);
  const auto reference=integrate_single(1e-6,1e-3,c);
  const std::array<double,4> a{{coarse.free_moles,coarse.dissolved_moles,coarse.mean_sw,static_cast<double>(coarse.state[0].c)}};
  const std::array<double,4> b{{reference.free_moles,reference.dissolved_moles,reference.mean_sw,static_cast<double>(reference.state[0].c)}};
  for(std::size_t i=0;i<a.size();++i){
    const double relative=std::abs(a[i]-b[i])/std::max(std::abs(b[i]),1e-30);
    require(relative<=1e-4,"T3 high-precision reference mismatch");
  }

  const auto y1=integrate_chain(1e-3,1e-3,c);
  const auto y2=integrate_chain(5e-4,1e-3,c);
  const auto y4=integrate_chain(2.5e-4,1e-3,c);
  const std::array<double,4> v1{{y1.free_moles,y1.dissolved_moles,y1.mean_sw,y1.outlet_moles}};
  const std::array<double,4> v2{{y2.free_moles,y2.dissolved_moles,y2.mean_sw,y2.outlet_moles}};
  const std::array<double,4> v4{{y4.free_moles,y4.dissolved_moles,y4.mean_sw,y4.outlet_moles}};
  for(std::size_t i=0;i<v1.size();++i){
    const double d12=std::abs(v1[i]-v2[i]);const double d24=std::abs(v2[i]-v4[i]);
    require(d24<d12,"dt convergence is not monotone for observable "+std::to_string(i));
    std::cout<<"dt_observable_"<<i<<" observed_order="<<std::log(d12/d24)/std::log(2.0)<<"\n";
  }
  std::cout<<"PASS T3 high-precision reference and dt convergence\n";
}

void test_disappearance() {
  {
    auto initial_config=base_config();
    initial_config.initial.concentration=0.2;
    const double sw=1.0-0.5*initial_config.active_set.sg_off;
    auto initial_network=single_bubble_network(sw);
    const double pl=initial_config.boundary.pl_in;
    const double pg=pl+bubble::capillary_pressure(sw,initial_network.pores[0].radius,
                                                   initial_config.physics);
    const double ng=bubble::free_gas_moles(pg,sw,initial_network.pores[0].volume,
                                           initial_config.physics);
    const double before=initial_network.pores[0].volume*sw*initial_config.initial.concentration+ng;
    bubble::Solver initial(std::move(initial_network),initial_config);initial.initialize();
    const double after=static_cast<double>(initial.network().pores[0].volume*initial.state()[0].c);
    require(!initial.state()[0].active_gas && initial.state()[0].permanently_inactive &&
            std::abs(after-before)<=1e-14*std::max(before,1e-300) &&
            std::abs(initial.initial_threshold_moles_transferred()-ng)<=1e-14*ng,
            "initial below-threshold bubble was not retired conservatively");
  }
  auto probe_config=base_config();probe_config.initial.concentration=0.0;
  bubble::Solver probe(single_bubble_network(),probe_config);probe.initialize();
  probe.mutable_state_for_test()[0].pl=probe_config.boundary.pl_out;
  probe.mutable_state_for_test()[0].pg=probe_config.boundary.pl_out+
      bubble::capillary_pressure(probe.state()[0].sw,probe.network().pores[0].radius,probe_config.physics);
  const double ng0=static_cast<double>(bubble::free_gas_moles(probe.state()[0].pg,probe.state()[0].sw,
                                          probe.network().pores[0].volume,probe_config.physics));
  const auto pr=probe.attempt_step(1e-4);require(pr.accepted,"T6 probe step failed");
  const double ng1=static_cast<double>(bubble::free_gas_moles(probe.state()[0].pg,probe.state()[0].sw,
                                          probe.network().pores[0].volume,probe_config.physics));
  require(ng1<ng0,"T6 probe did not dissolve gas");
  auto c=probe_config;c.active_set.ng_off=0.5*(ng0+ng1);
  bubble::Solver solver(single_bubble_network(),c);solver.initialize();
  solver.mutable_state_for_test()[0].pl=c.boundary.pl_out;
  solver.mutable_state_for_test()[0].pg=c.boundary.pl_out+
      bubble::capillary_pressure(solver.state()[0].sw,solver.network().pores[0].radius,c.physics);
  const int before=solver.make_dof_map(solver.state()).size;
  const auto r=solver.attempt_step(1e-4);
  require(r.accepted,"T6 disappearance step failed: "+r.failure);
  require(before==4 && solver.make_dof_map(solver.state()).size==2 && !solver.state()[0].active_gas &&
          solver.state()[0].permanently_inactive,"T6 active-set DOF removal failed");
  require(r.balance.epsilon_moles_adjusted<=1e-10 && r.balance.epsilon_volume_adjusted<=1e-10,
          "T6 adjusted conservation failed");
  const auto r2=solver.attempt_step(1e-6);require(r2.accepted,"T6 post-removal step failed");
  require(!solver.state()[0].active_gas,"T6 forbidden reactivation occurred");
  std::cout<<"PASS T6 disappearance transferred="<<r.threshold_moles_transferred<<"\n";
}

void test_batched_microbubble_retirement() {
  auto probe_config=base_config();probe_config.initial.concentration=0.0;
  bubble::Solver probe(two_bubble_network(),probe_config);probe.initialize();
  for(std::size_t i=0;i<2;++i){
    probe.mutable_state_for_test()[i].pl=probe_config.boundary.pl_out;
    probe.mutable_state_for_test()[i].pg=probe_config.boundary.pl_out+
        bubble::capillary_pressure(probe.state()[i].sw,probe.network().pores[i].radius,
                                   probe_config.physics);
  }
  const double ng0=static_cast<double>(bubble::free_gas_moles(
      probe.state()[0].pg,probe.state()[0].sw,probe.network().pores[0].volume,
      probe_config.physics));
  const auto probe_report=probe.attempt_step(1e-4);
  require(probe_report.accepted,"batch probe failed");
  const double ng1=static_cast<double>(bubble::free_gas_moles(
      probe.state()[0].pg,probe.state()[0].sw,probe.network().pores[0].volume,
      probe_config.physics));
  require(ng1<ng0,"batch probe did not dissolve gas");

  auto prepare=[&](bubble::Config config){
    config.time_step_control.enabled=true;
    config.active_set.ng_off=0.5*(ng0+ng1);
    bubble::Solver solver(two_bubble_network(),config);solver.initialize();
    for(std::size_t i=0;i<2;++i){
      solver.mutable_state_for_test()[i].pl=config.boundary.pl_out;
      solver.mutable_state_for_test()[i].pg=config.boundary.pl_out+
          bubble::capillary_pressure(solver.state()[i].sw,solver.network().pores[i].radius,
                                     config.physics);
    }
    return solver;
  };

  auto exact_config=probe_config;
  auto exact=prepare(exact_config);
  const auto exact_report=exact.attempt_step(1e-4);
  require(!exact_report.accepted && exact_report.event_timestep_localized &&
          exact_report.dt_change_reason=="gas_disappearance_event" &&
          exact_report.next_dt>0.0 && exact_report.next_dt<1e-4,
          "disabled batch path did not preserve exact localization");

  auto batch_config=probe_config;
  batch_config.active_set.batch_microbubble_retirement_enabled=true;
  batch_config.active_set.batch_individual_free_gas_fraction=1.0;
  batch_config.active_set.batch_total_free_gas_fraction_per_step=1.0;
  auto batch=prepare(batch_config);
  const auto batch_report=batch.attempt_step(1e-4);
  require(batch_report.accepted && batch_report.gas_disappearance_event &&
          !batch_report.event_timestep_localized &&
          batch_report.retired_bubble_count==2 &&
          batch_report.batched_retired_bubble_count==2 &&
          batch_report.gas_disappearance_event_records.size()==2 &&
          batch.gas_disappearance_events()==2,
          "two-bubble conservative batch retirement did not occur");
  require(std::all_of(batch_report.gas_disappearance_event_records.begin(),
                      batch_report.gas_disappearance_event_records.end(),
                      [](const auto& event){return event.batched &&
                          event.gas_moles_transferred_to_dissolved>=0.0;}),
          "batch event records are incomplete");
  require(!batch.state()[0].active_gas && !batch.state()[1].active_gas &&
          batch.state()[0].permanently_inactive &&
          batch.state()[1].permanently_inactive &&
          batch_report.balance.epsilon_moles_adjusted<=1e-10 &&
          batch_report.balance.epsilon_volume_adjusted<=1e-10,
          "batch retirement state or conservation audit failed");

  const std::filesystem::path checkpoint="/tmp/bubble_batch_retirement_checkpoint.json";
  const bubble::CheckpointProvenance provenance{
      std::string(64,'a'),std::string(64,'b'),std::string(64,'c'),std::string(64,'d')};
  bubble::write_checkpoint(checkpoint,batch,1e-6,0,provenance);
  const auto saved=bubble::read_checkpoint(checkpoint,batch.state().size());
  auto restarted=prepare(batch_config);
  restarted.restore(saved.state,saved.time,saved.step_number,
      saved.cumulative_mole_ledger,saved.cumulative_volume_ledger,
      saved.cumulative_net_liquid_inflow,saved.cumulative_boundary_component_outflow,
      saved.cumulative_outlet_advective_component_outflow,
      saved.cumulative_liquid_inflow,saved.cumulative_liquid_outflow,
      saved.cumulative_transfer_by_node);
  restarted.restore_gas_disappearance_metadata(saved.gas_disappearance_events,
      saved.last_gas_disappearance_node_id,saved.last_gas_disappearance_time);
  const auto continuous_next=batch.attempt_step(1e-6);
  const auto restarted_next=restarted.attempt_step(1e-6);
  require(continuous_next.accepted && restarted_next.accepted &&
          restarted.gas_disappearance_events()==2 &&
          restarted.state()[0].active_gas==batch.state()[0].active_gas &&
          restarted.state()[1].active_gas==batch.state()[1].active_gas &&
          std::abs(static_cast<double>(restarted.state()[0].c-batch.state()[0].c))<1e-12 &&
          std::abs(static_cast<double>(restarted.state()[1].c-batch.state()[1].c))<1e-12,
          "batch retirement checkpoint/restart equivalence failed");
  std::filesystem::remove(checkpoint);

  auto capped_config=batch_config;
  capped_config.active_set.batch_individual_free_gas_fraction=0.6;
  capped_config.active_set.batch_total_free_gas_fraction_per_step=0.75;
  auto capped=prepare(capped_config);
  const auto capped_report=capped.attempt_step(1e-4);
  require(!capped_report.accepted && capped_report.event_timestep_localized,
          "aggregate batch mass cap did not fall back to exact localization");
  std::cout<<"PASS T6b conservative batched microbubble retirement\n";
}

void compare_state(const std::vector<bubble::NodeState>& a,
                   const std::vector<bubble::NodeState>& b,double tolerance) {
  require(a.size()==b.size(),"state size mismatch");
  for(std::size_t i=0;i<a.size();++i){
    require(a[i].active_gas==b[i].active_gas && a[i].permanently_inactive==b[i].permanently_inactive,
            "state active flags mismatch");
    require(std::abs(a[i].pl-b[i].pl)<=tolerance*1e5 && std::abs(a[i].c-b[i].c)<=tolerance &&
            std::abs(a[i].sw-b[i].sw)<=tolerance,"state primary variables mismatch");
    if(a[i].active_gas)require(std::abs(a[i].pg-b[i].pg)<=tolerance*1e5,"state Pg mismatch");
  }
}

void test_rollback_restart_and_output() {
  auto failing=base_config();failing.initial.concentration=0.0;failing.nonlinear.max_iterations=2;
  failing.time.max_retries=40;failing.time.min_dt=1e-16;
  bubble::Solver rollback(single_bubble_network(),failing);rollback.initialize();
  rollback.mutable_state_for_test()[0].pl=failing.boundary.pl_out;
  rollback.mutable_state_for_test()[0].pg=failing.boundary.pl_out+
      bubble::capillary_pressure(rollback.state()[0].sw,rollback.network().pores[0].radius,failing.physics);
  const auto original=rollback.state();
  const auto failed=rollback.attempt_step(1.0);
  require(!failed.accepted,"T7 intentionally oversized step unexpectedly succeeded");
  compare_state(original,rollback.state(),0.0);
  require(rollback.time()==0.0 && rollback.step_number()==0 && rollback.cumulative_mole_ledger()==0.0,
          "T7 rollback did not restore solver metadata");
  const auto retried=rollback.advance_with_retry(1.0);
  require(retried.accepted && retried.dt<1.0,"T7 timestep reduction did not recover");
  std::cout<<"PASS T7 rollback recovered_dt="<<retried.dt<<"\n";

  auto c=base_config();c.initial.concentration=0.0;
  bubble::Solver continuous(chain_network(),c);continuous.initialize();
  bubble::Solver split(chain_network(),c);split.initialize();
  for(int k=0;k<4;++k){const auto r=continuous.attempt_step(1e-4);require(r.accepted,"T8 continuous step failed");}
  for(int k=0;k<2;++k){const auto r=split.attempt_step(1e-4);require(r.accepted,"T8 split pre-checkpoint failed");}
  const std::filesystem::path checkpoint="/tmp/bubble_solver_test_checkpoint.json";
  const bubble::CheckpointProvenance provenance{
      std::string(64,'1'),std::string(64,'2'),std::string(64,'3'),std::string(64,'4')};
  bubble::write_checkpoint(checkpoint,split,1e-4,0,provenance);
  const auto saved=bubble::read_checkpoint(checkpoint,split.state().size());
  bubble::Solver restarted(chain_network(),c);restarted.initialize();
  restarted.restore(saved.state,saved.time,saved.step_number,saved.cumulative_mole_ledger,
                    saved.cumulative_volume_ledger);
  for(int k=0;k<2;++k){const auto r=restarted.attempt_step(1e-4);require(r.accepted,"T8 restarted step failed");}
  compare_state(continuous.state(),restarted.state(),1e-10);
  require(std::abs(continuous.time()-restarted.time())<=1e-15,"T8 restart time mismatch");

  const std::filesystem::path vtk="/tmp/bubble_solver_test_output.vtk";
  bubble::write_vtk(vtk,restarted,nullptr);
  bubble::audit_vtk(vtk,restarted.network().pores.size(),restarted.network().edges.size());
  require(bubble::sha256_file(checkpoint).size()==64 && bubble::sha256_file(vtk).size()==64,
          "output checksum audit failed");
  std::cout<<"PASS T8 restart and VTK/checkpoint re-read\n";
}

void test_t9_zero_distance_two_members() {
  auto c=base_config();c.initial.concentration=0.2;
  auto n=make_network({pore(1,0.0,0,0.65),pore(2,0.0,0,0.8)},
                      {zero_edge(1,0,1)});
  bubble::Solver solver(std::move(n),c);solver.initialize();
  const auto& network=solver.network();
  require(network.zero_distance_edges.size()==1 &&
          network.zero_distance_clusters.size()==1 &&
          network.zero_distance_clusters[0].members.size()==2,
          "T9 union-find cluster classification failed");
  const auto map=solver.make_dof_map(solver.state());
  require(map.size==6 && map.node[0].pl==map.node[1].pl &&
          map.node[0].c==map.node[1].c && map.node[0].pg!=map.node[1].pg,
          "T9 shared Pl/C or independent Pg/Sw DOF map failed");
  require(solver.state()[0].pl==solver.state()[1].pl &&
          solver.state()[0].c==solver.state()[1].c,
          "T9 cluster state is not shared");
  long double liquid_members=0.0L,gas_members=0.0L,dissolved_members=0.0L;
  long double free_members=0.0L;
  for(std::size_t i=0;i<2;++i){
    const auto& p=network.pores[i];const auto& s=solver.state()[i];
    liquid_members+=static_cast<long double>(p.volume*s.sw);
    gas_members+=static_cast<long double>(p.volume*(1.0-s.sw));
    dissolved_members+=static_cast<long double>(p.volume*s.sw*s.c);
    free_members+=static_cast<long double>(bubble::free_gas_moles(
        s.pg,s.sw,p.volume,c.physics));
  }
  const long double liquid_cluster=liquid_members,gas_cluster=gas_members;
  const long double dissolved_cluster=dissolved_members,free_cluster=free_members;
  require(liquid_members==liquid_cluster && gas_members==gas_cluster &&
          dissolved_members==dissolved_cluster && free_members==free_cluster,
          "T9 cluster aggregation changed inventory");
  audit_jacobian(solver,solver.state());
  std::cout<<"PASS T9 zero-distance shared Pl/C and independent gas members\n";
}

void test_t10_multi_member_cluster() {
  auto c=base_config();c.initial.concentration=0.0;
  auto e3=edge(3,2,3);
  auto n=make_network(
      {pore(1,0.0,0,1.0),pore(2,0.0,0,0.7),pore(3,0.0,0,0.8),
       pore(4,1.0,2,1.0)},
      {zero_edge(1,0,1),zero_edge(2,0,2),e3});
  bubble::Solver solver(std::move(n),c);solver.initialize();
  const auto& network=solver.network();
  const auto cluster=network.cluster_of_pore[0];
  require(network.zero_distance_clusters[cluster].members.size()==3 &&
          network.zero_distance_clusters[cluster].external_positive_edges.size()==1,
          "T10 multi-storage union-find cluster failed");
  const auto map=solver.make_dof_map(solver.state());
  require(map.size==6 && map.node[0].pl==map.node[1].pl &&
          map.node[1].pl==map.node[2].pl && map.node[1].pg!=map.node[2].pg,
          "T10 equation/unknown count or member gas DOFs failed");
  std::vector<bool> upwind(network.edges.size(),true);
  std::vector<bool> transfer(solver.state().size(),false);
  transfer[1]=true;transfer[2]=true;
  const auto assembled=solver.assemble(solver.state(),solver.state(),1e-4,map,
                                       upwind,transfer,true);
  Eigen::MatrixXd dense=Eigen::MatrixXd(assembled.jacobian);
  for(int row=0;row<map.size;++row)
    for(int col=0;col<map.size;++col)
      dense(row,col)*=map.unknown_scale[static_cast<std::size_t>(col)]/
                      map.residual_scale[static_cast<std::size_t>(row)];
  Eigen::FullPivLU<Eigen::MatrixXd> rank(dense);
  require(rank.rank()==map.size,"T10 clustered Jacobian is singular");
  audit_jacobian(solver,solver.state());
  const auto report=solver.attempt_step(1e-4);
  require(report.accepted,"T10 clustered production step failed: "+report.failure);
  std::cout<<"PASS T10 multi-member cluster square nonsingular Jacobian\n";
}

void test_t11_massless_two_connection_junction() {
  auto c=base_config();c.initial.concentration=0.1;c.boundary.c_in=0.6;
  auto junction=pore(2,1.0,0,1.0);junction.volume=0.0;
  auto n=make_network({pore(1,0.0,1,1.0),junction,pore(3,2.0,2,1.0)},
                      {edge(1,0,1),edge(2,1,2)});
  bubble::Solver solver(std::move(n),c);solver.initialize();
  require(solver.network().massless_junctions==std::vector<std::size_t>{1},
          "T11 massless junction classification failed");
  const auto report=solver.attempt_step(1e-4);
  require(report.accepted,"T11 massless junction step failed: "+report.failure);
  const double g=bubble::pi*std::pow(solver.network().edges[0].radius,4)/
                 (8.0*c.physics.viscosity*solver.network().edges[0].length);
  const double expected=g*(c.boundary.pl_in-c.boundary.pl_out)/2.0;
  const double q_left=report.final_assembly.edge_fluxes[0].q;
  const double q_right=report.final_assembly.edge_fluxes[1].q;
  require(std::abs(q_left-expected)<=1e-12*std::abs(expected) &&
          std::abs(q_right-expected)<=1e-12*std::abs(expected),
          "T11 analytic series flow mismatch");
  const double f_left=report.final_assembly.edge_fluxes[0].advective+
                      report.final_assembly.edge_fluxes[0].diffusive;
  const double f_right=report.final_assembly.edge_fluxes[1].advective+
                       report.final_assembly.edge_fluxes[1].diffusive;
  require(std::abs(q_left-q_right)<=1e-14*std::abs(expected) &&
          std::abs(f_left-f_right)<=1e-12*std::max(std::abs(f_left),1e-300),
          "T11 massless flow/component flux does not balance");
  std::cout<<"PASS T11 massless two-connection analytic series junction\n";
}

void test_t12_massless_dead_end() {
  auto c=base_config();c.boundary.c_in=0.4;
  auto dead_end=pore(2,1.0,0,1.0);dead_end.volume=0.0;
  auto n=make_network({pore(1,0.0,1,1.0),dead_end},{edge(1,0,1)});
  bubble::Solver solver(std::move(n),c);solver.initialize();
  const auto report=solver.attempt_step(1e-4);
  require(report.accepted,"T12 massless dead-end step failed: "+report.failure);
  const auto& flux=report.final_assembly.edge_fluxes[0];
  require(std::abs(flux.q)<=1e-30 &&
          std::abs(flux.advective+flux.diffusive)<=1e-30 &&
          std::abs(solver.state()[1].pl-c.boundary.pl_in)<=1e-10 &&
          std::abs(solver.state()[1].c-c.boundary.c_in)<=1e-12,
          "T12 dead-end did not obtain zero net flux/non-singular state");
  std::cout<<"PASS T12 massless degree-one dead-end\n";
}

void test_t13_illegal_degenerate_geometry() {
  auto c=base_config();
  require_throws([&]{
    auto bad=pore(1,0.0,0,1.0);bad.volume=-1.0;
    bubble::Solver solver(make_network({bad},{}),c);
  },"volume must be finite and nonnegative");
  require_throws([&]{
    auto bad=edge(1,0,1);bad.length=-1.0;
    bubble::Solver solver(make_network({pore(1,0.0,0,1.0),pore(2,1.0,0,1.0)},
                                       {bad}),c);
  },"half-length must be finite and nonnegative");
  require_throws([&]{
    auto bad=pore(1,0.0,0,0.9);bad.volume=0.0;
    bubble::Solver solver(make_network({bad},{}),c);
  },"has zero volume with Sw<1");
  require_throws([&]{
    bubble::Solver solver(make_network({pore(1,0.0,1,1.0),pore(2,0.0,2,1.0)},
                                       {zero_edge(1,0,1)}),c);
  },"conflicting pressure boundaries");
  std::cout<<"PASS T13 illegal degenerate geometries fail fast\n";
}

void test_t14_virtual_boundary_controls() {
  auto c=base_config();c.virtual_boundary.enabled=true;c.virtual_boundary.voxel_size_um=7.5;
  c.geometry.voxel_size_um=7.5;c.geometry.image_dimensions_voxels={{4,1,1}};
  c.geometry.physical_domain_um={{30.0,7.5,7.5}};
  c.initial.concentration=0.1;c.boundary.c_in=1.0;c.boundary.c_backflow=0.8;
  auto n=make_network(
      {pore(1,0.0,0,0.65),pore(2,1e-5,0,1.0),pore(3,2e-5,0,1.0),
       virtual_pore(4,1,-1.75e-5,1),virtual_pore(5,3,3.75e-5,2)},
      {boundary_edge(1,3,0,1),edge(2,0,1),edge(3,1,2),boundary_edge(4,2,4,2)});
  bubble::Solver solver(std::move(n),c);solver.initialize();
  require(solver.network().virtual_boundary_nodes==2 && solver.network().boundary_edges==2,
          "T14 virtual boundary topology count failed");
  const auto map=solver.make_dof_map(solver.state());
  require(map.node[0].pl>=0 && map.node[0].pg>=0 && map.node[1].pl>=0 &&
          map.node[1].pg<0 && map.node[3].pl<0 && map.node[3].pg<0 &&
          map.node[4].pl<0 && map.node[4].pg<0,
          "T14 real/virtual residual set or DOF ownership failed");

  auto current=solver.state();
  current[3].c=12345.0;current[4].c=67890.0;
  bubble::Assembly no_flux;no_flux.edge_fluxes.resize(solver.network().edges.size());
  const auto balance=solver.compute_balance(solver.state(),current,no_flux,1e-4,0.0,0.0);
  const double expected=static_cast<double>(
      solver.network().pores[0].volume*current[0].sw*current[0].c+
      solver.network().pores[1].volume*current[1].sw*current[1].c+
      solver.network().pores[2].volume*current[2].sw*current[2].c);
  require(std::abs(balance.total_dissolved-expected)<=1e-15*std::max(expected,1e-300),
          "T14 virtual reservoir incorrectly entered internal storage");

  const auto backflow=bubble::evaluate_edge_flux(
      solver.network().edges[3],101000.0,0.2,1.0,102000.0,c.boundary.c_backflow,
      1.0,false,false,c.physics);
  require(backflow.q<0.0 && std::abs(backflow.advective-backflow.q*c.boundary.c_backflow)
              <=1e-14*std::abs(backflow.advective),
          "T14 virtual outlet backflow did not use configured concentration");

  const auto report=solver.attempt_step(1e-4);
  require(report.accepted,"T14 virtual boundary step failed: "+report.failure);
  require(report.final_assembly.edge_fluxes.front().q>0.0 &&
          report.final_assembly.edge_fluxes.front().upwind_a,
          "T14 inlet virtual-to-real flux/upwind direction failed");
  require(solver.state()[3].sw==1.0 && solver.state()[4].sw==1.0 &&
          !solver.state()[3].active_gas && !solver.state()[4].active_gas,
          "T14 virtual reservoir acquired bubble state");

  const std::filesystem::path checkpoint="/tmp/bubble_virtual_checkpoint.json";
  const bubble::CheckpointProvenance provenance{
      std::string(64,'1'),std::string(64,'2'),std::string(64,'3'),std::string(64,'4')};
  bubble::write_checkpoint(checkpoint,solver,1e-4,0,provenance);
  const auto saved=bubble::read_checkpoint(checkpoint,5);
  require(saved.state.size()==5,"T14 virtual boundary checkpoint/restart failed");
  require_throws([&]{(void)bubble::read_checkpoint(checkpoint,4);},
                 "checkpoint format/node count mismatch");
  const std::filesystem::path vtk="/tmp/bubble_virtual_output.vtk";
  bubble::write_vtk(vtk,solver,&report.final_assembly);
  bubble::audit_vtk(vtk,5,4);
  std::cout<<"PASS T14 one-to-one virtual boundary controls and output\n";
}

void test_t15_copy_adjacent_air_water_boundary() {
  auto make_case=[](double kl,bool reverse_endpoints){
    auto c=base_config();c.virtual_boundary.enabled=true;c.virtual_boundary.voxel_size_um=7.5;
    c.geometry.voxel_size_um=7.5;c.geometry.image_dimensions_voxels={{4,1,1}};
    c.geometry.physical_domain_um={{30.0,7.5,7.5}};
    c.physics.henry_cp=7.5e-6;c.physics.mass_transfer_coefficient=kl;
    c.physics.viscosity=8.90e-4;c.physics.molar_mass=0.028965;
    c.physics.surface_tension=0.07197;c.initial.concentration=1e-3;
    c.boundary.pl_in=102325.0;c.boundary.pl_out=101325.0;
    c.boundary.c_in=99.0;c.boundary.c_backflow=88.0;
    c.boundary.species_mode="copy_adjacent_concentration";
    c.boundary.fixed_concentrations_ignored=true;
    auto inlet=boundary_edge(1,3,0,1),outlet=boundary_edge(4,2,4,2);
    if(reverse_endpoints){std::swap(inlet.a,inlet.b);std::swap(outlet.a,outlet.b);}
    auto n=make_network(
        {pore(1,0.0,0,0.65),pore(2,1e-5,0,1.0),pore(3,2e-5,0,1.0),
         virtual_pore(4,1,-1.75e-5,1),virtual_pore(5,3,3.75e-5,2)},
        {inlet,edge(2,0,1),edge(3,1,2),outlet});
    return std::pair<bubble::Config,bubble::Network>{c,std::move(n)};
  };
  auto run=[&](double kl,bool reversed,double dt,int steps){
    auto input=make_case(kl,reversed);bubble::Solver solver(std::move(input.second),input.first);
    solver.initialize();
    const double gas0=static_cast<double>(bubble::free_gas_moles(
        solver.state()[0].pg,solver.state()[0].sw,solver.network().pores[0].volume,
        solver.config().physics));
    bubble::StepReport report;
    for(int i=0;i<steps;++i){report=solver.attempt_step(dt);require(report.accepted,
        "T15 copy-adjacent step failed: "+report.failure);}
    const double gas1=static_cast<double>(bubble::free_gas_moles(
        solver.state()[0].pg,solver.state()[0].sw,solver.network().pores[0].volume,
        solver.config().physics));
    return std::tuple<bubble::Solver,bubble::StepReport,double,double>{
        std::move(solver),std::move(report),gas0,gas1};
  };
  constexpr double reference_kl=2.3897375010715845e-5;
  auto [solver,report,gas0,gas1]=run(reference_kl,false,1e-4,20);
  require(gas1<gas0 && solver.state()[0].sw>0.65,
          "T15 air bubble did not shrink: gas0="+std::to_string(gas0)+
          " gas1="+std::to_string(gas1)+" Sw="+
          std::to_string(static_cast<double>(solver.state()[0].sw)));
  require(solver.state()[3].concentration_derived && solver.state()[4].concentration_derived &&
          solver.state()[3].c==solver.state()[0].c && solver.state()[4].c==solver.state()[2].c,
          "T15 virtual concentrations were not derived from adjacent real controls");
  for(const auto index:{std::size_t{0},std::size_t{3}}){
    const auto& flux=report.final_assembly.edge_fluxes[index];
    const auto real=index==0?std::size_t{0}:std::size_t{2};
    require(flux.diffusive==0.0 &&
            std::abs(flux.advective-flux.q*static_cast<double>(solver.state()[real].c))
                <=1e-12*std::max(std::abs(flux.advective),1e-300),
            "T15 boundary diffusion/advection rule failed");
  }
  require(report.balance.boundary_diffusive_outflow==0.0 &&
          std::abs(report.balance.liquid_volume_change/report.dt-
                   report.balance.net_liquid_inflow)
              <=1e-6*std::max(std::abs(report.balance.net_liquid_inflow),1e-300) &&
          report.balance.epsilon_moles_adjusted<=1e-10 &&
          report.balance.epsilon_volume_adjusted<=1e-10,
          "T15 component/liquid conservation failed: diff="+
          std::to_string(report.balance.boundary_diffusive_outflow)+" dVdt="+
          std::to_string(report.balance.liquid_volume_change/report.dt)+" net="+
          std::to_string(report.balance.net_liquid_inflow)+" epsN="+
          std::to_string(report.balance.epsilon_moles_adjusted)+" epsV="+
          std::to_string(report.balance.epsilon_volume_adjusted));
  auto [reversed,reversed_report,rg0,rg1]=run(reference_kl,true,1e-4,20);
  compare_state(solver.state(),reversed.state(),1e-9);
  require(std::abs(gas1-rg1)<=1e-10*gas0,
          "T15 endpoint reversal changed the solution");
  auto [half,half_report,hg0,hg1]=run(reference_kl,false,5e-5,40);
  require(std::abs(gas1-hg1)<=2e-3*std::abs(gas0-gas1),
          "T15 timestep halving did not converge");
  auto [low,low_report,lg0,lg1]=run(0.1*reference_kl,false,1e-4,20);
  auto [high,high_report,kg0,kg1]=run(10.0*reference_kl,false,1e-4,20);
  require((lg0-lg1)<(gas0-gas1) && (gas0-gas1)<(kg0-kg1),
          "T15 kL decade response is not monotone");
  const std::filesystem::path checkpoint="/tmp/bubble_copy_adjacent_checkpoint.json";
  const bubble::CheckpointProvenance provenance{
      std::string(64,'1'),std::string(64,'2'),std::string(64,'3'),std::string(64,'4')};
  bubble::write_checkpoint(checkpoint,solver,1e-4,0,provenance);
  const auto saved=bubble::read_checkpoint(checkpoint,5);
  require(saved.state[3].concentration_derived && saved.state[4].concentration_derived,
          "T15 checkpoint saved independent virtual concentration");
  auto restart_input=make_case(reference_kl,false);bubble::Solver restarted(
      std::move(restart_input.second),restart_input.first);restarted.initialize();
  restarted.restore(saved.state,saved.time,saved.step_number,saved.cumulative_mole_ledger,
                    saved.cumulative_volume_ledger,saved.cumulative_net_liquid_inflow,
                    saved.cumulative_boundary_component_outflow);
  compare_state(solver.state(),restarted.state(),1e-12);
  std::cout<<"PASS T15 copy-adjacent air-water flow/dissolution boundary\n";
}

void test_relative_pressure() {
  auto c=base_config();
  c.boundary.pl_in=101425; c.boundary.pl_out=101325;
  c.boundary.species_mode="copy_adjacent_concentration";
  c.boundary.fixed_concentrations_ignored=true;
  c.virtual_boundary.enabled=true;
  c.initial.concentration=0.001;
  c.physics.henry_cp=7.5e-6;
  auto n=make_network({pore(1,0,0,0.65),pore(2,1e-5,0,1),
      virtual_pore(3,1,-1.75e-5,1),virtual_pore(4,2,2.75e-5,2)},
      {boundary_edge(1,2,0,1),edge(2,0,1),boundary_edge(3,1,3,2)});
  bubble::Solver absolute(n,c); absolute.initialize();
  auto shifted=c; shifted.physics.background_abs_Pa=101325;
  shifted.boundary.pl_in=100; shifted.boundary.pl_out=0;
  bubble::Solver relative(n,shifted); relative.initialize();
  auto other=shifted;other.physics.background_abs_Pa=101300;
  other.boundary.pl_in=125;other.boundary.pl_out=25;
  bubble::Solver reference_changed(n,other); reference_changed.initialize();
  for(int step=0;step<10;++step){
    require(absolute.attempt_step(1e-4).accepted,"absolute comparison step");
    require(relative.attempt_step(1e-4).accepted,"relative comparison step");
    require(reference_changed.attempt_step(1e-4).accepted,"shifted reference step");
  }
  for(std::size_t i=0;i<n.pores.size();++i){
    const auto& a=absolute.state()[i];const auto& b=relative.state()[i];
    const auto& d=reference_changed.state()[i];
    require(std::abs(a.pl-(101325+b.pl))<1e-7 && std::abs(a.sw-b.sw)<1e-10 &&
        std::abs(a.c-b.c)<1e-10,"relative/absolute physical solution mismatch");
    require(std::abs((101300+d.pl)-(101325+b.pl))<1e-7 &&
        std::abs(d.sw-b.sw)<1e-10,"reference origin changed physics");
  }
  auto map=relative.make_dof_map(relative.state());
  std::vector<bool> upwind(n.edges.size(),true),transfer(n.pores.size(),false);
  transfer[0]=true;
  const auto ad=relative.assemble(relative.state(),relative.state(),1e-4,map,upwind,transfer,true);
  for(int column=0;column<map.size;++column){
    const double h=column==map.node[0].sw?1e-6:1e-4;
    auto plus=relative.state(),minus=plus;set_dof(plus,map,column,h);set_dof(minus,map,column,-h);
    auto rp=relative.assemble(plus,relative.state(),1e-4,map,upwind,transfer,false).residual;
    auto rm=relative.assemble(minus,relative.state(),1e-4,map,upwind,transfer,false).residual;
    for(int row=0;row<map.size;++row){
      const double fd=(rp[row]-rm[row])/(2*h),actual=ad.jacobian.coeff(row,column);
      require(std::abs(fd-actual)<=1e-6*std::max({std::abs(fd),std::abs(actual),1e-18}),
              "relative copy-boundary Jacobian mismatch");
    }
  }
  auto physical=shifted.physics;
  const double low=bubble::free_gas_moles(500.,0.5,1e-12,physical);
  physical.background_abs_Pa*=2;
  const double high=bubble::free_gas_moles(500.,0.5,1e-12,physical);
  require(high>1.99*low && physical.henry_cp*(physical.background_abs_Pa+500.)>
      shifted.physics.henry_cp*(shifted.physics.background_abs_Pa+500.),"absolute EOS/Henry background response");
  bubble::NodeState tiny;tiny.pl=1;tiny.pl+=1e-12L;
  require(tiny.pl>1 && static_cast<double>(1e16+1e-12)==1e16,
          "tiny relative pressure correction was lost");
  const auto file=std::filesystem::path("/tmp/bubble_relative_pressure_checkpoint.json");
  const bubble::CheckpointProvenance provenance{std::string(64,'1'),std::string(64,'2'),std::string(64,'3'),std::string(64,'4')};
  bubble::write_checkpoint(file,relative,1e-4,0,provenance);
  const auto saved=bubble::read_checkpoint(file,n.pores.size());
  require(saved.background_abs_Pa==101325,"checkpoint background missing");
  bubble::Solver restarted(n,shifted);restarted.initialize();
  restarted.restore(saved.state,saved.time,saved.step_number,saved.cumulative_mole_ledger,
      saved.cumulative_volume_ledger,saved.cumulative_net_liquid_inflow,
      saved.cumulative_boundary_component_outflow,saved.cumulative_outlet_advective_component_outflow);
  restarted.restore_gross_fluxes(saved.cumulative_component_inflow,saved.cumulative_component_outflow);
  require(restarted.attempt_step(1e-4).accepted && relative.attempt_step(1e-4).accepted,"relative restart step");
  compare_state(restarted.state(),relative.state(),1e-12);
  std::cout<<"PASS T16 relative pressure invariance, EOS/Henry, AD, tiny updates and restart\n";
}

void test_mass_transfer_multiplier() {
  auto base=base_config();
  base.physics.background_abs_Pa=101325.0;
  base.boundary.pl_in=100.0;base.boundary.pl_out=0.0;
  base.physics.henry_cp=7.5e-6;base.initial.concentration=1e-3;
  auto accelerated=base;accelerated.physics.mass_transfer_multiplier=100.0;
  bubble::Solver reference(single_bubble_network(),base);
  bubble::Solver preview(single_bubble_network(),accelerated);
  reference.initialize();preview.initialize();
  const auto map_reference=reference.make_dof_map(reference.state());
  const auto map_preview=preview.make_dof_map(preview.state());
  std::vector<bool> upwind(reference.network().edges.size(),true);
  std::vector<bool> transfer(reference.state().size(),false);transfer[0]=true;
  const auto r1=reference.assemble(reference.state(),reference.state(),1e-4,
      map_reference,upwind,transfer,true);
  const auto r100=preview.assemble(preview.state(),preview.state(),1e-4,
      map_preview,upwind,transfer,true);
  require(r1.transfer[0]>0.0 &&
      std::abs(r100.transfer[0]/r1.transfer[0]-100.0)<1e-12,
      "mass-transfer multiplier was not applied exactly once");
  const auto rd1=r1.residual[map_reference.node[0].rd];
  const auto rg1=r1.residual[map_reference.node[0].rg];
  const auto rd100=r100.residual[map_preview.node[0].rd];
  const auto rg100=r100.residual[map_preview.node[0].rg];
  require(std::abs((rd100-rd1)+(rg100-rg1))<1e-12*
      std::max(std::abs(rd100-rd1),1e-300),
      "accelerated transfer is not conservative between Rd and Rg");
  const std::filesystem::path checkpoint="/tmp/bubble_multiplier_checkpoint.json";
  const bubble::CheckpointProvenance provenance{std::string(64,'1'),std::string(64,'2'),
      std::string(64,'3'),std::string(64,'4')};
  bubble::write_checkpoint(checkpoint,preview,1e-4,0,provenance);
  const auto saved=bubble::read_checkpoint(checkpoint,preview.state().size());
  require(saved.mass_transfer_multiplier==100.0,
      "checkpoint did not preserve mass-transfer multiplier");
  const std::filesystem::path vtk="/tmp/bubble_multiplier.vtk";
  bubble::write_vtk(vtk,preview,&r100);
  for(const auto& path:{vtk}){
    std::ifstream input(path,std::ios::binary);std::string bytes(
        (std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
    require(bytes.find("mass_transfer_multiplier")!=std::string::npos &&
        bytes.find("permanently_inactive")!=std::string::npos &&
        bytes.find("dissolved_moles")!=std::string::npos,
        "accelerated-preview VTK metadata is incomplete");
  }
  bubble::Solver restarted(single_bubble_network(),accelerated);restarted.initialize();
  restarted.restore(saved.state,saved.time,saved.step_number,saved.cumulative_mole_ledger,
      saved.cumulative_volume_ledger,saved.cumulative_net_liquid_inflow,
      saved.cumulative_boundary_component_outflow,
      saved.cumulative_outlet_advective_component_outflow,
      saved.cumulative_liquid_inflow,saved.cumulative_liquid_outflow,
      saved.cumulative_transfer_by_node);
  const double free_before=static_cast<double>(free_gas_moles(preview.state()[0].pg,
      preview.state()[0].sw,preview.network().pores[0].volume,accelerated.physics));
  for(int step=0;step<5;++step){
    const auto a=preview.attempt_step(1e-4),b=restarted.attempt_step(1e-4);
    require(a.accepted&&b.accepted,"accelerated checkpoint continuation failed");
  }
  compare_state(preview.state(),restarted.state(),1e-12);
  const double free_after=static_cast<double>(free_gas_moles(preview.state()[0].pg,
      preview.state()[0].sw,preview.network().pores[0].volume,accelerated.physics));
  require(free_after<free_before,"accelerated free gas did not decrease monotonically");
  std::cout<<"PASS T17 conservative mass-transfer multiplier and provenance output\n";
}

/* Flow-boundary integration fixture retained for follow-up full-network validation.
void test_flow_inlet_boundary() {
  auto c=base_config();
  c.virtual_boundary.enabled=true;
  c.geometry.voxel_size_um=7.5;c.virtual_boundary.voxel_size_um=7.5;
  c.geometry.image_dimensions_voxels={{4,1,1}};
  c.geometry.physical_domain_um={{30.0,7.5,7.5}};
  c.boundary.inlet_type="flow";c.boundary.outlet_type="pressure";
  c.boundary.inlet_flow_m3_s=2.0e-10;c.boundary.pl_in=5.0e5;c.boundary.pl_out=0.0;
  c.boundary.species_mode="fixed_reservoir_concentration";c.boundary.fixed_concentrations_ignored=false;c.boundary.c_in=0.0;c.boundary.c_backflow=0.0;
  c.physics.background_abs_Pa=101325.0;c.physics.henry_cp=7.5e-6;c.initial.concentration=0.0;
  auto inlet_a=boundary_edge(1,3,0,1),inlet_b=boundary_edge(2,4,1,1);
  inlet_b.radius=3e-6;inlet_b.area=bubble::pi*9e-12;
  auto n=make_network(
      {pore(1,0.0,0,1.0),pore(2,0.0,0,1.0),pore(3,2e-5,0,1.0),
       virtual_pore(4,1,-1.75e-5,1),virtual_pore(5,2,-1.75e-5,1),
       virtual_pore(6,3,3.75e-5,2)},
      {inlet_a,inlet_b,edge(3,0,2),edge(4,1,2),boundary_edge(5,2,5,2)});
  bubble::Solver solver(std::move(n),c);solver.initialize();
  auto report=solver.attempt_step(1e-4);
  require(report.accepted,"flow inlet boundary step failed: "+report.failure);
  const auto& net=solver.network();
  double inlet=0.0;int inlet_edges=0;std::vector<double> branch_flows;
  for(std::size_t k=0;k<net.edges.size();++k){const auto&e=net.edges[k];
    if(net.pores[e.a].boundary_flag==c.boundary.inlet_flag||net.pores[e.b].boundary_flag==c.boundary.inlet_flag){
      const double value=std::abs(report.final_assembly.edge_fluxes[k].q);
      inlet+=value;branch_flows.push_back(value);++inlet_edges;}}
  require(inlet_edges>0&&std::abs(inlet-c.boundary.inlet_flow_m3_s)<=1e-12*c.boundary.inlet_flow_m3_s,
      "specified inlet flow was not enforced");
  require(branch_flows.size()==2&&std::abs(branch_flows[0]-branch_flows[1])>1e-3*c.boundary.inlet_flow_m3_s,
      "inlet flow was prescribed per edge instead of distributed by conductance");
  require(report.balance.epsilon_moles_adjusted<=1e-10&&report.balance.epsilon_volume_adjusted<=1e-10,
      "flow inlet boundary conservation failed");
  std::cout<<"PASS T18 prescribed inlet flow / outlet pressure / degassed inlet\n";
}
*/

}  // namespace

int main(){
  test_t1_t2_t4_fluxes();
  test_jacobian_and_dof_map();
  test_monolithic_steps();
  test_reference_and_dt_convergence();
  test_disappearance();
  test_batched_microbubble_retirement();
  test_rollback_restart_and_output();
  test_t9_zero_distance_two_members();
  test_t10_multi_member_cluster();
  test_t11_massless_two_connection_junction();
  test_t12_massless_dead_end();
  test_t13_illegal_degenerate_geometry();
  test_t14_virtual_boundary_controls();
  test_t15_copy_adjacent_air_water_boundary();
  test_relative_pressure();
  test_mass_transfer_multiplier();
  return 0;
}
