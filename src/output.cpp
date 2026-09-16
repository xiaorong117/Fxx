#include "bubble/model.hpp"

#include <nlohmann/json.hpp>

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <zlib.h>

namespace bubble {
namespace {

using json = nlohmann::json;

std::string extended_string(long double value) {
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<long double>::max_digits10) << value;
  return out.str();
}

void require_output(std::ofstream& out, const std::filesystem::path& path) {
  if (!out) throw std::runtime_error("failed writing output: " + path.string());
}

template <class Function>
void point_scalar(std::ofstream& out, const char* name, std::size_t n, Function value) {
  out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
  for (std::size_t i = 0; i < n; ++i) out << value(i) << '\n';
}

template <class Function>
void cell_scalar(std::ofstream& out, const char* name, std::size_t n, Function value) {
  out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
  for (std::size_t i = 0; i < n; ++i) out << value(i) << '\n';
}

std::uint32_t rotate_right(std::uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

std::string sha256_bytes(const std::string& input) {
  static constexpr std::array<std::uint32_t, 64> k{{
      0x428a2f98U,0x71374491U,0xb5c0fbcfU,0xe9b5dba5U,0x3956c25bU,0x59f111f1U,0x923f82a4U,0xab1c5ed5U,
      0xd807aa98U,0x12835b01U,0x243185beU,0x550c7dc3U,0x72be5d74U,0x80deb1feU,0x9bdc06a7U,0xc19bf174U,
      0xe49b69c1U,0xefbe4786U,0x0fc19dc6U,0x240ca1ccU,0x2de92c6fU,0x4a7484aaU,0x5cb0a9dcU,0x76f988daU,
      0x983e5152U,0xa831c66dU,0xb00327c8U,0xbf597fc7U,0xc6e00bf3U,0xd5a79147U,0x06ca6351U,0x14292967U,
      0x27b70a85U,0x2e1b2138U,0x4d2c6dfcU,0x53380d13U,0x650a7354U,0x766a0abbU,0x81c2c92eU,0x92722c85U,
      0xa2bfe8a1U,0xa81a664bU,0xc24b8b70U,0xc76c51a3U,0xd192e819U,0xd6990624U,0xf40e3585U,0x106aa070U,
      0x19a4c116U,0x1e376c08U,0x2748774cU,0x34b0bcb5U,0x391c0cb3U,0x4ed8aa4aU,0x5b9cca4fU,0x682e6ff3U,
      0x748f82eeU,0x78a5636fU,0x84c87814U,0x8cc70208U,0x90befffaU,0xa4506cebU,0xbef9a3f7U,0xc67178f2U}};
  std::vector<std::uint8_t> bytes(input.begin(), input.end());
  const std::uint64_t bit_length = static_cast<std::uint64_t>(bytes.size()) * 8U;
  bytes.push_back(0x80U);
  while (bytes.size() % 64 != 56) bytes.push_back(0U);
  for (int shift = 56; shift >= 0; shift -= 8)
    bytes.push_back(static_cast<std::uint8_t>((bit_length >> shift) & 0xffU));
  std::array<std::uint32_t, 8> h{{0x6a09e667U,0xbb67ae85U,0x3c6ef372U,0xa54ff53aU,
                                  0x510e527fU,0x9b05688cU,0x1f83d9abU,0x5be0cd19U}};
  for (std::size_t block = 0; block < bytes.size(); block += 64) {
    std::array<std::uint32_t, 64> w{};
    for (std::size_t i = 0; i < 16; ++i) {
      const std::size_t p = block + 4 * i;
      w[i] = (static_cast<std::uint32_t>(bytes[p]) << 24) |
             (static_cast<std::uint32_t>(bytes[p + 1]) << 16) |
             (static_cast<std::uint32_t>(bytes[p + 2]) << 8) | bytes[p + 3];
    }
    for (std::size_t i = 16; i < 64; ++i) {
      const auto s0 = rotate_right(w[i-15],7) ^ rotate_right(w[i-15],18) ^ (w[i-15] >> 3);
      const auto s1 = rotate_right(w[i-2],17) ^ rotate_right(w[i-2],19) ^ (w[i-2] >> 10);
      w[i] = w[i-16] + s0 + w[i-7] + s1;
    }
    auto a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];
    for (std::size_t i = 0; i < 64; ++i) {
      const auto s1=rotate_right(e,6)^rotate_right(e,11)^rotate_right(e,25);
      const auto ch=(e&f)^((~e)&g);
      const auto t1=hh+s1+ch+k[i]+w[i];
      const auto s0=rotate_right(a,2)^rotate_right(a,13)^rotate_right(a,22);
      const auto maj=(a&b)^(a&c)^(b&c);
      const auto t2=s0+maj;
      hh=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=hh;
  }
  std::ostringstream out;
  out << std::hex << std::setfill('0');
  for (auto value : h) out << std::setw(8) << value;
  return out.str();
}

}  // namespace

void write_vtk(const std::filesystem::path& path, const Solver& solver,
               const Assembly* assembly) {
  std::ofstream out(path);
  require_output(out, path);
  out << std::scientific << std::setprecision(17);
  const auto& network = solver.network();
  const auto& state = solver.state();
  const auto dofs = solver.make_dof_map(state);
  double xmin=std::numeric_limits<double>::infinity(),xmax=-xmin;
  for(const auto&p:network.pores) if(p.source_type!=2 && p.volume>0){xmin=std::min(xmin,p.x);xmax=std::max(xmax,p.x);}
  std::vector<double> constriction(state.size(),0.0);
  for(std::size_t i=0;i<state.size();++i){std::vector<double> ratios;
    if(network.pores[i].source_type!=2) for(auto k:network.adjacency[i]){const auto&e=network.edges[k];
      if(e.source_old_throat_id && e.radius>0)ratios.push_back(network.pores[i].radius/e.radius);}
    if(!ratios.empty()){std::sort(ratios.begin(),ratios.end());constriction[i]=ratios[ratios.size()/2];}}
  const auto displayed_c = [&](std::size_t i) -> long double {
    if (network.pores[i].source_type != 2 ||
        solver.config().boundary.species_mode != "copy_adjacent_concentration")
      return state[i].c;
    const auto edge_index = network.adjacency[i].front();
    const auto& edge = network.edges[edge_index];
    const auto neighbor = edge.a == i ? edge.b : edge.a;
    return state[neighbor].c;
  };
  out << "# vtk DataFile Version 3.0\nbubble dissolution solver\nASCII\nDATASET POLYDATA\n";
  out << "POINTS " << network.pores.size() << " double\n";
  for (const auto& p : network.pores) out << p.x << ' ' << p.y << ' ' << p.z << '\n';
  out << "LINES " << network.edges.size() << ' ' << 3 * network.edges.size() << "\n";
  for (const auto& e : network.edges) out << "2 " << e.a << ' ' << e.b << '\n';
  out << "POINT_DATA " << state.size() << '\n';
  const auto displayed_pl = [&](std::size_t i) -> long double {
    return solver.config().boundary.inlet_type == "flow" &&
        network.pores[i].source_type == 2 &&
        network.pores[i].boundary_flag == solver.config().boundary.inlet_flag
        ? solver.inlet_pressure() : state[i].pl;
  };
  point_scalar(out, "Pl", state.size(), [&](std::size_t i){ return displayed_pl(i); });
  point_scalar(out, "Pl_relative_Pa", state.size(), [&](std::size_t i){ return displayed_pl(i); });
  point_scalar(out, "Pg_relative_Pa", state.size(), [&](std::size_t i){ return state[i].active_gas ? state[i].pg : 0.0L; });
  point_scalar(out, "Pl_absolute_Pa", state.size(), [&](std::size_t i){ return solver.config().physics.background_abs_Pa + displayed_pl(i); });
  point_scalar(out, "Pg_absolute_Pa", state.size(), [&](std::size_t i){ return state[i].active_gas ? solver.config().physics.background_abs_Pa + state[i].pg : 0.0L; });
  point_scalar(out, "Pg", state.size(), [&](std::size_t i){ return state[i].active_gas ? state[i].pg : -1.0; });
  point_scalar(out, "Sw", state.size(), [&](std::size_t i){ return state[i].sw; });
  point_scalar(out, "Sg", state.size(), [&](std::size_t i){ return 1.0-state[i].sw; });
  point_scalar(out, "C", state.size(), [&](std::size_t i){ return displayed_c(i); });
  point_scalar(out, "n_g", state.size(), [&](std::size_t i){
    return state[i].active_gas ? free_gas_moles(state[i].pg,state[i].sw,network.pores[i].volume,solver.config().physics) : 0.0; });
  point_scalar(out, "ntr", state.size(), [&](std::size_t i){ return assembly && i < assembly->transfer.size() && state[i].active_gas ? assembly->transfer[i] : 0.0; });
  point_scalar(out, "active_gas", state.size(), [&](std::size_t i){ return state[i].active_gas ? 1.0 : 0.0; });
  point_scalar(out, "permanently_inactive", state.size(), [&](std::size_t i){ return state[i].permanently_inactive ? 1.0 : 0.0; });
  point_scalar(out, "dissolved_moles", state.size(), [&](std::size_t i){
    return network.pores[i].source_type == 2 ? 0.0L
        : network.pores[i].volume * state[i].sw * displayed_c(i); });
  point_scalar(out, "mass_transfer_multiplier", state.size(), [&](std::size_t){
    return solver.config().physics.mass_transfer_multiplier; });
  point_scalar(out, "batch_microbubble_retirement_enabled", state.size(), [&](std::size_t){
    return solver.config().active_set.batch_microbubble_retirement_enabled ? 1.0 : 0.0; });
  point_scalar(out, "batch_individual_free_gas_fraction", state.size(), [&](std::size_t){
    return solver.config().active_set.batch_individual_free_gas_fraction; });
  point_scalar(out, "batch_total_free_gas_fraction_per_step", state.size(), [&](std::size_t){
    return solver.config().active_set.batch_total_free_gas_fraction_per_step; });
  point_scalar(out, "node_id", state.size(), [&](std::size_t i){ return static_cast<double>(network.pores[i].id); });
  point_scalar(out, "source_type", state.size(), [&](std::size_t i){ return static_cast<double>(network.pores[i].source_type); });
  point_scalar(out, "source_id", state.size(), [&](std::size_t i){ return static_cast<double>(network.pores[i].source_id); });
  point_scalar(out, "boundary_flag", state.size(), [&](std::size_t i){ return static_cast<double>(network.pores[i].boundary_flag); });
  point_scalar(out, "is_virtual_boundary", state.size(), [&](std::size_t i){
    return network.pores[i].source_type == 2 ? 1.0 : 0.0; });
  point_scalar(out, "component_id", state.size(), [&](std::size_t i){
    return i < network.component_id.size() ? static_cast<double>(network.component_id[i]) : 0.0; });
  point_scalar(out, "boundary_connection_class", state.size(), [&](std::size_t i){
    return i < network.boundary_connection_class.size()
        ? static_cast<double>(network.boundary_connection_class[i]) : 0.0; });
  point_scalar(out, "gas_moles", state.size(), [&](std::size_t i){
    return state[i].active_gas ? free_gas_moles(state[i].pg, state[i].sw,
        network.pores[i].volume, solver.config().physics) : 0.0; });
  point_scalar(out, "gas_volume", state.size(), [&](std::size_t i){
    return state[i].active_gas ? network.pores[i].volume * (1.0L - state[i].sw) : 0.0L; });
  point_scalar(out, "bubble_equivalent_radius", state.size(), [&](std::size_t i){
    if (!state[i].active_gas) return 0.0L;
    using std::pow;
    const long double volume = network.pores[i].volume * (1.0L - state[i].sw);
    return pow(3.0L * volume / (4.0L * static_cast<long double>(pi)), 1.0L / 3.0L);
  });
  point_scalar(out, "pressure_gauge", state.size(), [&](std::size_t i){ return dofs.node[i].pressure_gauge ? 1.0 : 0.0; });
  point_scalar(out, "zero_distance_cluster", state.size(), [&](std::size_t i){
    return static_cast<double>(network.cluster_of_pore[i] + 1); });
  point_scalar(out, "massless_junction", state.size(), [&](std::size_t i){
    return network.pores[i].volume == 0.0 ? 1.0 : 0.0; });
  point_scalar(out, "normalized_x", state.size(), [&](std::size_t i){
    return network.pores[i].source_type==2 ? -1.0 : std::clamp((network.pores[i].x-xmin)/(xmax-xmin),0.0,1.0); });
  point_scalar(out, "segment_id", state.size(), [&](std::size_t i){
    if(network.pores[i].source_type==2 || network.pores[i].volume==0.0)return 0.0;
    return 1.0+std::min(solver.config().analysis.x_segments-1,
      static_cast<int>(std::clamp((network.pores[i].x-xmin)/(xmax-xmin),0.0,1.0)*solver.config().analysis.x_segments)); });
  point_scalar(out, "local_Sg", state.size(), [&](std::size_t i){ return 1.0-state[i].sw; });
  point_scalar(out, "initial_Sg", state.size(), [&](std::size_t i){ return 1.0-solver.initial_sw_by_node()[i]; });
  point_scalar(out, "free_gas_moles_lost", state.size(), [&](std::size_t i){
    const long double now=state[i].active_gas?free_gas_moles(state[i].pg,state[i].sw,network.pores[i].volume,solver.config().physics):0.0L;
    return solver.initial_free_gas_by_node()[i]-now; });
  point_scalar(out, "cumulative_interphase_transfer", state.size(), [&](std::size_t i){ return solver.cumulative_transfer_by_node()[i]; });
  point_scalar(out, "adjacent_constriction_ratio_p50", state.size(), [&](std::size_t i){ return constriction[i]; });
  out << "CELL_DATA " << network.edges.size() << '\n';
  cell_scalar(out, "Q", network.edges.size(), [&](std::size_t i){ return assembly ? assembly->edge_fluxes[i].q : 0.0; });
  cell_scalar(out, "liquid_flux_m3_s", network.edges.size(), [&](std::size_t i){ return assembly ? assembly->edge_fluxes[i].q : 0.0; });
  cell_scalar(out, "advective_component_flux_mol_s", network.edges.size(), [&](std::size_t i){ return assembly ? assembly->edge_fluxes[i].advective : 0.0; });
  cell_scalar(out, "diffusive_component_flux_mol_s", network.edges.size(), [&](std::size_t i){ return assembly ? assembly->edge_fluxes[i].diffusive : 0.0; });
  cell_scalar(out, "liquid_flux", network.edges.size(), [&](std::size_t i){ return assembly ? assembly->edge_fluxes[i].q : 0.0; });
  cell_scalar(out, "F_adv", network.edges.size(), [&](std::size_t i){ return assembly ? assembly->edge_fluxes[i].advective : 0.0; });
  cell_scalar(out, "advective_component_flux", network.edges.size(), [&](std::size_t i){ return assembly ? assembly->edge_fluxes[i].advective : 0.0; });
  cell_scalar(out, "F_diff", network.edges.size(), [&](std::size_t i){ return assembly ? assembly->edge_fluxes[i].diffusive : 0.0; });
  cell_scalar(out, "diffusive_component_flux", network.edges.size(), [&](std::size_t i){ return assembly ? assembly->edge_fluxes[i].diffusive : 0.0; });
  cell_scalar(out, "source_old_throat_id", network.edges.size(), [&](std::size_t i){ return static_cast<double>(network.edges[i].source_old_throat_id); });
  cell_scalar(out, "is_boundary_edge", network.edges.size(), [&](std::size_t i){
    const auto& edge = network.edges[i];
    return network.pores[edge.a].source_type == 2 ||
        network.pores[edge.b].source_type == 2 ? 1.0 : 0.0; });
  cell_scalar(out, "edge_length", network.edges.size(), [&](std::size_t i){
    return network.edges[i].length; });
  cell_scalar(out, "edge_length_voxels", network.edges.size(), [&](std::size_t i){
    return network.edges[i].length /
        (solver.config().virtual_boundary.voxel_size_um * 1.0e-6); });
  cell_scalar(out, "side", network.edges.size(), [&](std::size_t i){ return static_cast<double>(network.edges[i].side); });
  cell_scalar(out, "zero_distance_constraint", network.edges.size(), [&](std::size_t i){
    return network.edges[i].length == 0.0 ? 1.0 : 0.0; });
  require_output(out, path);
}

void write_vtp(const std::filesystem::path& path, const Solver& solver,
               const Assembly* assembly) {
  struct Array { std::string name,type; int components=1; std::string bytes,payload; std::uint64_t offset=0; };
  std::vector<Array> points_data,cells_data;const auto& n=solver.network();const auto&s=solver.state();
  auto add_double=[&](std::vector<Array>& arrays,const char*name,int components,auto value,std::size_t count){
    std::vector<double> data(count*static_cast<std::size_t>(components));
    for(std::size_t i=0;i<count;++i)for(int c=0;c<components;++c)data[i*components+c]=static_cast<double>(value(i,c));
    Array a;a.name=name;a.type="Float64";a.components=components;
    a.bytes.assign(reinterpret_cast<const char*>(data.data()),data.size()*sizeof(double));arrays.push_back(std::move(a));};
  double xmin=std::numeric_limits<double>::infinity(),xmax=-xmin;for(const auto&p:n.pores)if(p.source_type!=2&&p.volume>0){xmin=std::min(xmin,p.x);xmax=std::max(xmax,p.x);}
  std::vector<double> constriction(s.size());for(std::size_t i=0;i<s.size();++i){std::vector<double> r;
    if(n.pores[i].source_type!=2){for(auto k:n.adjacency[i])if(n.edges[k].source_old_throat_id&&n.edges[k].radius>0)r.push_back(n.pores[i].radius/n.edges[k].radius);}
    if(!r.empty()){std::sort(r.begin(),r.end());constriction[i]=r[r.size()/2];}}
  const auto dc=[&](std::size_t i){if(n.pores[i].source_type!=2||solver.config().boundary.species_mode!="copy_adjacent_concentration")return s[i].c;
    const auto&e=n.edges[n.adjacency[i].front()];return s[e.a==i?e.b:e.a].c;};
  const auto dp=[&](std::size_t i){return solver.config().boundary.inlet_type=="flow"&&n.pores[i].source_type==2&&n.pores[i].boundary_flag==solver.config().boundary.inlet_flag?static_cast<long double>(solver.inlet_pressure()):s[i].pl;};
  add_double(points_data,"Pl_relative_Pa",1,[&](auto i,int){return dp(i);},s.size());
  add_double(points_data,"Pg_relative_Pa",1,[&](auto i,int){return s[i].active_gas?s[i].pg:0;},s.size());
  add_double(points_data,"Pl_absolute_Pa",1,[&](auto i,int){return solver.config().physics.background_abs_Pa+dp(i);},s.size());
  add_double(points_data,"Pg_absolute_Pa",1,[&](auto i,int){return s[i].active_gas?solver.config().physics.background_abs_Pa+s[i].pg:0;},s.size());
  add_double(points_data,"C",1,[&](auto i,int){return dc(i);},s.size());
  add_double(points_data,"Sw",1,[&](auto i,int){return s[i].sw;},s.size());
  add_double(points_data,"Sg",1,[&](auto i,int){return 1-s[i].sw;},s.size());
  add_double(points_data,"active_gas",1,[&](auto i,int){return s[i].active_gas;},s.size());
  add_double(points_data,"permanently_inactive",1,[&](auto i,int){return s[i].permanently_inactive;},s.size());
  add_double(points_data,"dissolved_moles",1,[&](auto i,int){return n.pores[i].source_type==2?0:n.pores[i].volume*s[i].sw*dc(i);},s.size());
  add_double(points_data,"mass_transfer_multiplier",1,[&](auto,int){return solver.config().physics.mass_transfer_multiplier;},s.size());
  add_double(points_data,"batch_microbubble_retirement_enabled",1,[&](auto,int){return solver.config().active_set.batch_microbubble_retirement_enabled;},s.size());
  add_double(points_data,"batch_individual_free_gas_fraction",1,[&](auto,int){return solver.config().active_set.batch_individual_free_gas_fraction;},s.size());
  add_double(points_data,"batch_total_free_gas_fraction_per_step",1,[&](auto,int){return solver.config().active_set.batch_total_free_gas_fraction_per_step;},s.size());
  add_double(points_data,"node_id",1,[&](auto i,int){return n.pores[i].id;},s.size());
  add_double(points_data,"source_type",1,[&](auto i,int){return n.pores[i].source_type;},s.size());
  add_double(points_data,"boundary_flag",1,[&](auto i,int){return n.pores[i].boundary_flag;},s.size());
  add_double(points_data,"is_virtual_boundary",1,[&](auto i,int){return n.pores[i].source_type==2;},s.size());
  add_double(points_data,"component_id",1,[&](auto i,int){return n.component_id[i];},s.size());
  add_double(points_data,"normalized_x",1,[&](auto i,int){return n.pores[i].source_type==2?-1:std::clamp((n.pores[i].x-xmin)/(xmax-xmin),0.0,1.0);},s.size());
  add_double(points_data,"segment_id",1,[&](auto i,int){return n.pores[i].source_type==2||n.pores[i].volume==0?0:1+std::min(solver.config().analysis.x_segments-1,static_cast<int>(std::clamp((n.pores[i].x-xmin)/(xmax-xmin),0.0,1.0)*solver.config().analysis.x_segments));},s.size());
  add_double(points_data,"local_Sg",1,[&](auto i,int){return 1-s[i].sw;},s.size());
  add_double(points_data,"initial_Sg",1,[&](auto i,int){return 1-solver.initial_sw_by_node()[i];},s.size());
  add_double(points_data,"gas_moles",1,[&](auto i,int){return s[i].active_gas?free_gas_moles(s[i].pg,s[i].sw,n.pores[i].volume,solver.config().physics):0;},s.size());
  add_double(points_data,"gas_volume",1,[&](auto i,int){return s[i].active_gas?n.pores[i].volume*(1-s[i].sw):0;},s.size());
  add_double(points_data,"bubble_equivalent_radius",1,[&](auto i,int){return s[i].active_gas?std::cbrt(3.0*static_cast<double>(n.pores[i].volume*(1-s[i].sw))/(4*pi)):0;},s.size());
  add_double(points_data,"free_gas_moles_lost",1,[&](auto i,int){const auto now=s[i].active_gas?free_gas_moles(s[i].pg,s[i].sw,n.pores[i].volume,solver.config().physics):0;return solver.initial_free_gas_by_node()[i]-now;},s.size());
  add_double(points_data,"cumulative_interphase_transfer",1,[&](auto i,int){return solver.cumulative_transfer_by_node()[i];},s.size());
  add_double(points_data,"adjacent_constriction_ratio_p50",1,[&](auto i,int){return constriction[i];},s.size());
  add_double(cells_data,"liquid_flux_m3_s",1,[&](auto i,int){return assembly?assembly->edge_fluxes[i].q:0;},n.edges.size());
  add_double(cells_data,"advective_component_flux_mol_s",1,[&](auto i,int){return assembly?assembly->edge_fluxes[i].advective:0;},n.edges.size());
  add_double(cells_data,"diffusive_component_flux_mol_s",1,[&](auto i,int){return assembly?assembly->edge_fluxes[i].diffusive:0;},n.edges.size());
  add_double(cells_data,"is_boundary_edge",1,[&](auto i,int){return n.pores[n.edges[i].a].source_type==2||n.pores[n.edges[i].b].source_type==2;},n.edges.size());
  add_double(cells_data,"edge_length_m",1,[&](auto i,int){return n.edges[i].length;},n.edges.size());
  Array coordinates;coordinates.name="Points";coordinates.type="Float64";coordinates.components=3;{std::vector<double>d;d.reserve(n.pores.size()*3);for(const auto&p:n.pores){d.push_back(p.x);d.push_back(p.y);d.push_back(p.z);}coordinates.bytes.assign(reinterpret_cast<const char*>(d.data()),d.size()*8);}
  Array connectivity,offsets;connectivity.name="connectivity";connectivity.type="Int64";offsets.name="offsets";offsets.type="Int64";{std::vector<std::int64_t>c,o;c.reserve(n.edges.size()*2);o.reserve(n.edges.size());for(std::size_t i=0;i<n.edges.size();++i){c.push_back(n.edges[i].a);c.push_back(n.edges[i].b);o.push_back(2*(i+1));}connectivity.bytes.assign(reinterpret_cast<const char*>(c.data()),c.size()*8);offsets.bytes.assign(reinterpret_cast<const char*>(o.data()),o.size()*8);}
  std::vector<Array*> order;for(auto&a:points_data)order.push_back(&a);for(auto&a:cells_data)order.push_back(&a);order.push_back(&coordinates);order.push_back(&connectivity);order.push_back(&offsets);
  std::uint64_t offset=0;for(auto*a:order){uLongf bound=compressBound(a->bytes.size());std::string compressed(bound,'\0');if(compress2(reinterpret_cast<Bytef*>(compressed.data()),&bound,reinterpret_cast<const Bytef*>(a->bytes.data()),a->bytes.size(),Z_BEST_SPEED)!=Z_OK)throw std::runtime_error("zlib compression failed");compressed.resize(bound);
    const std::uint64_t header[4]={1,a->bytes.size(),a->bytes.size(),compressed.size()};a->payload.assign(reinterpret_cast<const char*>(header),sizeof(header));a->payload+=compressed;a->offset=offset;offset+=a->payload.size();}
  std::ofstream out(path,std::ios::binary);if(!out)throw std::runtime_error("cannot write VTP");out<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">\n<PolyData><Piece NumberOfPoints=\""<<n.pores.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\""<<n.edges.size()<<"\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n<PointData>\n";
  const auto tag=[&](const Array&a){out<<"<DataArray type=\""<<a.type<<"\" Name=\""<<a.name<<"\" NumberOfComponents=\""<<a.components<<"\" format=\"appended\" offset=\""<<a.offset<<"\"/>\n";};
  for(const auto&a:points_data)tag(a);
  out<<"</PointData><CellData>\n";
  for(const auto&a:cells_data)tag(a);
  out<<"</CellData><Points>\n";tag(coordinates);out<<"</Points><Lines>\n";
  tag(connectivity);tag(offsets);out<<"</Lines></Piece></PolyData>\n<AppendedData encoding=\"raw\">_";
  for(auto*a:order)out.write(a->payload.data(),a->payload.size());
  out<<"</AppendedData></VTKFile>\n";if(!out)throw std::runtime_error("failed writing compressed VTP");
}

void write_degenerate_geometry_audit(const std::filesystem::path& path,
                                     const Solver& solver) {
  const auto& network = solver.network();
  const auto& state = solver.state();
  const auto dofs = solver.make_dof_map(state);
  json root;
  root["format"] = "bubble_degenerate_geometry_audit_v1";
  root["formulation"] =
      "raw geometry unchanged; union-find zero-distance elimination plus massless junctions";
  root["raw_node_count"] = network.pores.size();
  root["raw_edge_count"] = network.edges.size();
  root["effective_cluster_count"] = network.zero_distance_clusters.size();
  root["production_positive_edge_count"] = network.production_positive_edges;
  root["zero_distance_edge_count"] = network.zero_distance_edges.size();
  root["massless_junction_count"] = network.massless_junctions.size();
  root["throats_with_two_zero_halves"] = network.throats_with_two_zero_halves;
  root["massless_junction_node_ids"] = json::array();
  root["zero_distance_edges"] = json::array();

  for (const auto node : network.massless_junctions)
    root["massless_junction_node_ids"].push_back(network.pores[node].id);
  for (const auto edge_index : network.zero_distance_edges) {
    const auto& edge = network.edges[edge_index];
    root["zero_distance_edges"].push_back({
        {"edge_id", edge.id}, {"pore_id_1", network.pores[edge.a].id},
        {"pore_id_2", network.pores[edge.b].id},
        {"source_old_throat_id", edge.source_old_throat_id}, {"side", edge.side}});
  }
  root["nontrivial_clusters"] = json::array();
  for (const auto& cluster : network.zero_distance_clusters) {
    if (cluster.members.size() == 1 && cluster.zero_distance_edges.empty()) continue;
    json item;
    item["cluster_id"] = cluster.id + 1;
    item["boundary_flag"] = cluster.boundary_flag;
    item["shared_Pl_dof"] = dofs.cluster[cluster.id].pl;
    item["shared_C_dof"] = dofs.cluster[cluster.id].c;
    item["member_node_ids"] = json::array();
    item["zero_distance_edge_ids"] = json::array();
    item["external_positive_edge_ids"] = json::array();
    for (const auto node : cluster.members)
      item["member_node_ids"].push_back(network.pores[node].id);
    for (const auto edge : cluster.zero_distance_edges)
      item["zero_distance_edge_ids"].push_back(network.edges[edge].id);
    for (const auto edge : cluster.external_positive_edges)
      item["external_positive_edge_ids"].push_back(network.edges[edge].id);
    root["nontrivial_clusters"].push_back(std::move(item));
  }

  long double raw_storage = 0.0L, clustered_storage = 0.0L;
  long double raw_liquid = 0.0L, clustered_liquid = 0.0L;
  long double raw_gas_volume = 0.0L, clustered_gas_volume = 0.0L;
  long double raw_moles = 0.0L, clustered_moles = 0.0L;
  for (std::size_t i = 0; i < network.pores.size(); ++i) {
    const auto& pore = network.pores[i];
    if (network.zero_distance_clusters[network.cluster_of_pore[i]].boundary_flag != 0)
      continue;
    raw_storage += static_cast<long double>(pore.volume);
    const auto& value = state[i];
    raw_liquid += static_cast<long double>(pore.volume * value.sw);
    raw_gas_volume += static_cast<long double>(pore.volume * (1.0 - value.sw));
    raw_moles += static_cast<long double>(pore.volume * value.sw * value.c);
    if (value.active_gas)
      raw_moles += static_cast<long double>(free_gas_moles(
          value.pg, value.sw, pore.volume, solver.config().physics));
  }
  for (const auto& cluster : network.zero_distance_clusters) {
    for (const auto node : cluster.members) {
      const auto& pore = network.pores[node];
      if (cluster.boundary_flag != 0) continue;
      clustered_storage += static_cast<long double>(pore.volume);
      const auto& value = state[node];
      clustered_liquid += static_cast<long double>(pore.volume * value.sw);
      clustered_gas_volume += static_cast<long double>(pore.volume * (1.0 - value.sw));
      clustered_moles += static_cast<long double>(pore.volume * value.sw * value.c);
      if (value.active_gas)
        clustered_moles += static_cast<long double>(free_gas_moles(
            value.pg, value.sw, pore.volume, solver.config().physics));
    }
  }
  const auto relative = [](long double a, long double b) {
    return static_cast<double>(std::abs(a - b) /
                               std::max(std::abs(a), static_cast<long double>(1e-300)));
  };
  root["aggregation_invariants"] = {
      {"raw_storage_volume_m3", static_cast<double>(raw_storage)},
      {"clustered_storage_volume_m3", static_cast<double>(clustered_storage)},
      {"storage_relative_error", relative(raw_storage, clustered_storage)},
      {"raw_initial_liquid_volume_m3", static_cast<double>(raw_liquid)},
      {"clustered_initial_liquid_volume_m3", static_cast<double>(clustered_liquid)},
      {"liquid_volume_relative_error", relative(raw_liquid, clustered_liquid)},
      {"raw_initial_gas_volume_m3", static_cast<double>(raw_gas_volume)},
      {"clustered_initial_gas_volume_m3", static_cast<double>(clustered_gas_volume)},
      {"gas_volume_relative_error", relative(raw_gas_volume, clustered_gas_volume)},
      {"raw_initial_total_moles", static_cast<double>(raw_moles)},
      {"clustered_initial_total_moles", static_cast<double>(clustered_moles)},
      {"total_moles_relative_error", relative(raw_moles, clustered_moles)}};
  if (relative(raw_storage, clustered_storage) > 1e-15 ||
      relative(raw_liquid, clustered_liquid) > 1e-15 ||
      relative(raw_gas_volume, clustered_gas_volume) > 1e-15 ||
      relative(raw_moles, clustered_moles) > 1e-15)
    throw std::runtime_error("zero-distance cluster aggregation changed an invariant");

  std::ofstream out(path);
  out << std::setprecision(17) << std::setw(2) << root << '\n';
  require_output(out, path);
}

void write_component_boundary_and_makeup_audit(const std::filesystem::path& path,
                                               const Solver& solver) {
  struct ComponentStats {
    int connection_class = 0;
    std::size_t real_nodes = 0, virtual_nodes = 0, active_bubbles = 0;
    long double free_moles = 0.0L, gas_volume = 0.0L;
  };
  const auto& network = solver.network();
  const auto& state = solver.state();
  std::map<std::size_t, ComponentStats> components;
  for (std::size_t i = 0; i < network.pores.size(); ++i) {
    const std::size_t component = network.component_id[i];
    auto& stats = components[component];
    stats.connection_class = network.boundary_connection_class[i];
    if (network.pores[i].source_type == 2) {
      ++stats.virtual_nodes;
      continue;
    }
    ++stats.real_nodes;
    if (!state[i].active_gas) continue;
    ++stats.active_bubbles;
    stats.gas_volume += static_cast<long double>(network.pores[i].volume) *
                        (1.0L - state[i].sw);
    stats.free_moles += free_gas_moles(state[i].pg, state[i].sw,
                                      network.pores[i].volume,
                                      solver.config().physics);
  }
  static constexpr std::array<const char*, 5> names{{
      "boundaryless", "inlet_only", "outlet_only", "inlet_outlet",
      "isolated_real_virtual"}};
  json root;
  root["format"] = "bubble_component_boundary_and_makeup_audit_v1";
  root["boundary_connection_class_codes"] = {
      {"0", names[0]}, {"1", names[1]}, {"2", names[2]},
      {"3", names[3]}, {"4", names[4]}};
  root["main_flow_component_id"] = network.main_flow_component_id;
  for (std::size_t code = 0; code < names.size(); ++code) {
    std::size_t count = 0, active = 0;
    long double free_moles = 0.0L, gas_volume = 0.0L;
    for (const auto& item : components) {
      if (item.second.connection_class != static_cast<int>(code)) continue;
      ++count;
      active += item.second.active_bubbles;
      free_moles += item.second.free_moles;
      gas_volume += item.second.gas_volume;
    }
    root["classes"][names[code]] = {
        {"component_count", count}, {"active_bubble_count", active},
        {"initial_free_gas_mol", static_cast<double>(free_moles)},
        {"initial_gas_volume_m3", static_cast<double>(gas_volume)},
        {"can_receive_external_liquid", code != 0}};
  }
  const auto main = components.at(network.main_flow_component_id);
  root["main_flow_component"] = {
      {"component_id", network.main_flow_component_id},
      {"real_node_count", main.real_nodes},
      {"virtual_node_count", main.virtual_nodes},
      {"active_bubble_count", main.active_bubbles},
      {"initial_free_gas_mol", static_cast<double>(main.free_moles)},
      {"initial_gas_volume_m3", static_cast<double>(main.gas_volume)}};
  root["components"] = json::array();
  for (const auto& item : components) {
    root["components"].push_back({
        {"component_id", item.first},
        {"boundary_connection_class", item.second.connection_class},
        {"real_node_count", item.second.real_nodes},
        {"virtual_node_count", item.second.virtual_nodes},
        {"active_bubble_count", item.second.active_bubbles},
        {"initial_free_gas_mol", static_cast<double>(item.second.free_moles)},
        {"initial_gas_volume_m3", static_cast<double>(item.second.gas_volume)}});
  }
  std::ofstream out(path);
  out << std::setprecision(17) << std::setw(2) << root << '\n';
  require_output(out, path);
}

void write_physical_parameter_audit(const std::filesystem::path& path,
                                    const Solver& solver) {
  const auto& network = solver.network();
  const auto& state = solver.state();
  std::vector<double> radii;
  std::vector<std::pair<double, double>> weighted;
  json invalid_driving_force = json::array();
  for (std::size_t i = 0; i < state.size(); ++i) {
    if (!state[i].active_gas) continue;
    const double gas_volume = static_cast<double>(
        network.pores[i].volume * (1.0L - state[i].sw));
    const double radius = std::cbrt(3.0 * gas_volume / (4.0 * pi));
    const double moles = static_cast<double>(free_gas_moles(
        state[i].pg, state[i].sw, network.pores[i].volume,
        solver.config().physics));
    radii.push_back(radius);
    weighted.emplace_back(radius, moles);
    if (!(solver.config().physics.henry_cp * (solver.config().physics.background_abs_Pa + state[i].pg) - state[i].c > 0.0))
      invalid_driving_force.push_back(network.pores[i].id);
  }
  if (radii.empty()) throw std::runtime_error("physical audit found no active bubbles");
  std::sort(radii.begin(), radii.end());
  std::sort(weighted.begin(), weighted.end());
  const auto quantile = [&](double probability) {
    const double position = probability * static_cast<double>(radii.size() - 1);
    const auto lower = static_cast<std::size_t>(std::floor(position));
    const auto upper = static_cast<std::size_t>(std::ceil(position));
    const double fraction = position - static_cast<double>(lower);
    return radii[lower] + fraction * (radii[upper] - radii[lower]);
  };
  long double total_weight = 0.0L;
  for (const auto& item : weighted) total_weight += item.second;
  long double cumulative = 0.0L;
  double weighted_median = weighted.back().first;
  for (const auto& item : weighted) {
    cumulative += item.second;
    if (cumulative >= 0.5L * total_weight) {
      weighted_median = item.first;
      break;
    }
  }
  const double reference_kl = solver.config().physics.molecular_diffusivity /
                              weighted_median;
  json root;
  root["format"] = "bubble_air_water_physical_parameter_audit_v1";
  root["reference_properties"] = {
      {"T_K", solver.config().physics.temperature},
      {"effective_dry_air_molar_mass_kg_mol", solver.config().physics.molar_mass},
      {"compressibility_Z", solver.config().physics.compressibility},
      {"Hcp_mol_m3_Pa", solver.config().physics.henry_cp},
      {"Dm_m2_s", solver.config().physics.molecular_diffusivity},
      {"pure_water_viscosity_Pa_s", solver.config().physics.viscosity},
      {"pure_water_density_kg_m3", solver.config().physics.liquid_density},
      {"surface_tension_N_m", solver.config().physics.surface_tension}};
  root["unvalidated_closure_parameters"] = {
      {"contact_angle_rad", solver.config().physics.contact_angle},
      {"tortuosity", solver.config().physics.tortuosity},
      {"relative_permeability_exponent", solver.config().physics.kr_exponent},
      {"qin_b", solver.config().physics.qin_b},
      {"mass_transfer_correlation", "stationary sphere Sh=2; kL=D/Rb_reference"}};
  root["initial_active_bubble_count"] = radii.size();
  root["initial_dissolution_driving_force_failed_node_ids"] = invalid_driving_force;
  root["all_active_bubbles_initially_dissolving"] = invalid_driving_force.empty();
  root["bubble_equivalent_radius_m"] = {
      {"minimum", radii.front()}, {"p10", quantile(0.1)},
      {"median", quantile(0.5)}, {"p90", quantile(0.9)},
      {"maximum", radii.back()},
      {"free_gas_mole_weighted_median", weighted_median}};
  root["kL_reference_m_s"] = reference_kl;
  root["configured_kL_m_s"] = solver.config().physics.mass_transfer_coefficient;
  root["mass_transfer_multiplier"] = solver.config().physics.mass_transfer_multiplier;
  root["effective_kL_m_s"] = solver.config().physics.mass_transfer_multiplier *
                              solver.config().physics.mass_transfer_coefficient;
  root["configured_to_reference_ratio"] =
      solver.config().physics.mass_transfer_coefficient / reference_kl;
  root["reference_Sherwood_number"] = reference_kl * 2.0 * weighted_median /
                                      solver.config().physics.molecular_diffusivity;
  root["background_abs_Pa"] = solver.config().physics.background_abs_Pa;
  root["inlet_pressure_relative_Pa"] = solver.inlet_pressure();
  root["pressure_representation"] = "relative_to_background";
  root["initial_dissolved_moles"] = 0.0;
  long double dissolved = 0.0L, liquid = 0.0L;
  for (std::size_t i=0; i<state.size(); ++i) {
    if(network.pores[i].source_type==2) continue;
    liquid += network.pores[i].volume * state[i].sw;
    dissolved += network.pores[i].volume * state[i].sw * state[i].c;
  }
  root["initial_dissolved_moles"] = static_cast<double>(dissolved);
  root["initial_liquid_volume_m3"] = static_cast<double>(liquid);
  root["initial_free_gas_moles"] = static_cast<double>(total_weight);
  std::ofstream out(path);
  out << std::setprecision(17) << std::setw(2) << root << '\n';
  require_output(out, path);
}

void audit_vtk(const std::filesystem::path& path, std::size_t expected_points,
               std::size_t expected_edges) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot re-read VTK: " + path.string());
  std::ostringstream buffer;
  buffer << input.rdbuf();
  const auto text = buffer.str();
  const std::array<const char*, 21> required{{"POINTS ","LINES ","POINT_DATA ","CELL_DATA ",
      "SCALARS Pl ","SCALARS Pg ","SCALARS Sw ","SCALARS Sg ","SCALARS C ",
      "SCALARS n_g ","SCALARS ntr ","SCALARS active_gas ","SCALARS node_id ",
      "SCALARS source_type ","SCALARS source_id ","SCALARS boundary_flag ",
      "SCALARS is_virtual_boundary ","SCALARS pressure_gauge ",
      "SCALARS zero_distance_cluster ","SCALARS massless_junction ","SCALARS Q "}};
  for (const auto* field : required) {
    if (text.find(field) == std::string::npos) throw std::runtime_error(std::string("missing VTK field: ") + field);
  }
  if (text.find("SCALARS F_adv ") == std::string::npos || text.find("SCALARS F_diff ") == std::string::npos ||
      text.find("SCALARS source_old_throat_id ") == std::string::npos ||
      text.find("SCALARS is_boundary_edge ") == std::string::npos ||
      text.find("SCALARS edge_length ") == std::string::npos ||
      text.find("SCALARS edge_length_voxels ") == std::string::npos ||
      text.find("SCALARS component_id ") == std::string::npos ||
      text.find("SCALARS boundary_connection_class ") == std::string::npos ||
      text.find("SCALARS gas_moles ") == std::string::npos ||
      text.find("SCALARS gas_volume ") == std::string::npos ||
      text.find("SCALARS bubble_equivalent_radius ") == std::string::npos ||
      text.find("SCALARS liquid_flux ") == std::string::npos ||
      text.find("SCALARS advective_component_flux ") == std::string::npos ||
      text.find("SCALARS diffusive_component_flux ") == std::string::npos ||
      text.find("SCALARS side ") == std::string::npos ||
      text.find("SCALARS zero_distance_constraint ") == std::string::npos)
    throw std::runtime_error("missing required VTK edge field");
  if (text.find("POINTS " + std::to_string(expected_points) + " ") == std::string::npos ||
      text.find("LINES " + std::to_string(expected_edges) + " ") == std::string::npos)
    throw std::runtime_error("VTK point/edge count mismatch");
}

void write_checkpoint(const std::filesystem::path& path, const Solver& solver,
                      double next_dt, int consecutive_fast_steps,
                      const CheckpointProvenance& provenance) {
  if (!std::isfinite(next_dt) || next_dt <= 0.0 || consecutive_fast_steps < 0) {
    throw std::runtime_error("invalid adaptive state for checkpoint");
  }
  if (provenance.continuation_config_sha256.size() != 64 ||
      provenance.pore_sha256.size() != 64 || provenance.throat_sha256.size() != 64 ||
      provenance.connect_sha256.size() != 64) {
    throw std::runtime_error("invalid provenance for checkpoint");
  }
  json root;
  root["format"] = "bubble_checkpoint_v6";
  root["pressure_representation"] = "relative_to_background";
  root["background_abs_Pa"] = solver.config().physics.background_abs_Pa;
  root["inlet_pressure_relative_Pa"] = solver.inlet_pressure();
  root["time_s"] = solver.time();
  root["step_number"] = solver.step_number();
  root["cumulative_component_inflow_mol"] = solver.cumulative_component_inflow();
  root["cumulative_component_outflow_mol"] = solver.cumulative_component_outflow();
  root["cumulative_liquid_inflow_m3"] = solver.cumulative_liquid_inflow();
  root["cumulative_liquid_outflow_m3"] = solver.cumulative_liquid_outflow();
  root["cumulative_interphase_transfer_mol"] =
      solver.cumulative_interphase_transfer();
  root["cumulative_interphase_transfer_by_node_mol"] = json::array();
  for (const auto value : solver.cumulative_transfer_by_node())
    root["cumulative_interphase_transfer_by_node_mol"].push_back(extended_string(value));
  root["cumulative_mole_ledger"] = solver.cumulative_mole_ledger();
  root["cumulative_volume_ledger"] = solver.cumulative_volume_ledger();
  root["cumulative_net_liquid_inflow_m3"] = solver.cumulative_net_liquid_inflow();
  root["cumulative_boundary_component_outflow_mol"] =
      solver.cumulative_boundary_component_outflow();
  root["cumulative_outlet_advective_component_outflow_mol"] =
      solver.cumulative_outlet_advective_component_outflow();
  root["next_dt_s"] = next_dt;
  root["consecutive_fast_steps"] = consecutive_fast_steps;
  root["last_dt_change_reason"] = solver.last_dt_change_reason();
  root["last_min_free_gas_moles"] = solver.last_min_free_gas_moles();
  root["last_min_Sg"] = solver.last_min_sg();
  root["mass_transfer_multiplier"] =
      solver.config().physics.mass_transfer_multiplier;
  root["batch_microbubble_retirement_enabled"] =
      solver.config().active_set.batch_microbubble_retirement_enabled;
  root["batch_individual_free_gas_fraction"] =
      solver.config().active_set.batch_individual_free_gas_fraction;
  root["batch_total_free_gas_fraction_per_step"] =
      solver.config().active_set.batch_total_free_gas_fraction_per_step;
  root["gas_disappearance_events"] = solver.gas_disappearance_events();
  root["last_gas_disappearance_node_id"] =
      solver.last_gas_disappearance_node_id();
  root["last_gas_disappearance_time_s"] =
      solver.last_gas_disappearance_time();
  root["provenance"] = {
      {"continuation_config_sha256", provenance.continuation_config_sha256},
      {"pore_sha256", provenance.pore_sha256},
      {"throat_sha256", provenance.throat_sha256},
      {"connect_sha256", provenance.connect_sha256}};
  root["state"] = json::array();
  root["state_extended_precision"] = json::array();
  for (std::size_t i = 0; i < solver.state().size(); ++i) {
    const auto& s = solver.state()[i];
    const bool derived_c = solver.network().pores[i].source_type == 2 &&
        solver.config().boundary.species_mode == "copy_adjacent_concentration";
    json row = {static_cast<double>(s.pl),
                s.active_gas ? static_cast<double>(s.pg) : -1.0,
                static_cast<double>(s.sw), nullptr,
                s.active_gas, s.permanently_inactive};
    if (!derived_c) row[3] = static_cast<double>(s.c);
    root["state"].push_back(std::move(row));
    root["state_extended_precision"].push_back({
        extended_string(s.pl), s.active_gas ? extended_string(s.pg) : "nan",
        extended_string(s.sw), derived_c ? "derived_from_adjacent_real" : extended_string(s.c)});
  }
  std::ofstream out(path);
  out << std::setprecision(17) << root.dump(1) << '\n';
  require_output(out, path);
}

CheckpointData read_checkpoint(const std::filesystem::path& path,
                               std::size_t expected_nodes) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open checkpoint: " + path.string());
  json root = json::parse(input);
  const std::string format = root.at("format").get<std::string>();
  if ((format != "bubble_checkpoint_v5" && format != "bubble_checkpoint_v6") ||
      root.at("state").size() != expected_nodes)
    throw std::runtime_error("checkpoint format/node count mismatch");
  CheckpointData data;
  if (root.at("pressure_representation") != "relative_to_background")
    throw std::runtime_error("checkpoint pressure representation mismatch");
  data.background_abs_Pa = root.at("background_abs_Pa").get<double>();
  data.inlet_pressure = root.value("inlet_pressure_relative_Pa", 0.0);
  data.time = root.at("time_s").get<double>();
  data.step_number = root.at("step_number").get<int>();
  data.cumulative_component_inflow = root.at("cumulative_component_inflow_mol").get<double>();
  data.cumulative_component_outflow = root.at("cumulative_component_outflow_mol").get<double>();
  data.cumulative_mole_ledger = root.at("cumulative_mole_ledger").get<double>();
  data.cumulative_volume_ledger = root.at("cumulative_volume_ledger").get<double>();
  if (format == "bubble_checkpoint_v5" || format == "bubble_checkpoint_v6") {
    data.cumulative_net_liquid_inflow =
        root.at("cumulative_net_liquid_inflow_m3").get<double>();
    data.cumulative_boundary_component_outflow =
        root.at("cumulative_boundary_component_outflow_mol").get<double>();
    data.cumulative_outlet_advective_component_outflow =
        root.at("cumulative_outlet_advective_component_outflow_mol").get<double>();
  }
  if (format == "bubble_checkpoint_v6") {
    data.cumulative_liquid_inflow = root.at("cumulative_liquid_inflow_m3").get<double>();
    data.cumulative_liquid_outflow = root.at("cumulative_liquid_outflow_m3").get<double>();
    const auto& values = root.at("cumulative_interphase_transfer_by_node_mol");
    if (values.size() != expected_nodes)
      throw std::runtime_error("checkpoint cumulative transfer node count mismatch");
    data.cumulative_transfer_by_node.reserve(expected_nodes);
    for (const auto& value : values)
      data.cumulative_transfer_by_node.push_back(std::stold(value.get<std::string>()));
  }
  data.next_dt = root.at("next_dt_s").get<double>();
  data.consecutive_fast_steps = root.at("consecutive_fast_steps").get<int>();
  if(root.contains("last_dt_change_reason"))data.last_dt_change_reason=root.at("last_dt_change_reason").get<std::string>();
  if(root.contains("last_min_free_gas_moles"))data.last_min_free_gas_moles=root.at("last_min_free_gas_moles").get<double>();
  if(root.contains("last_min_Sg"))data.last_min_sg=root.at("last_min_Sg").get<double>();
  data.mass_transfer_multiplier = root.value("mass_transfer_multiplier", 1.0);
  data.batch_microbubble_retirement_enabled =
      root.value("batch_microbubble_retirement_enabled", false);
  data.batch_individual_free_gas_fraction =
      root.value("batch_individual_free_gas_fraction", 1e-10);
  data.batch_total_free_gas_fraction_per_step =
      root.value("batch_total_free_gas_fraction_per_step", 1e-8);
  data.gas_disappearance_events = root.value("gas_disappearance_events", 0ULL);
  data.last_gas_disappearance_node_id =
      root.value("last_gas_disappearance_node_id", std::size_t{0});
  data.last_gas_disappearance_time =
      root.value("last_gas_disappearance_time_s", -1.0);
  const auto& provenance = root.at("provenance");
  data.continuation_config_sha256 =
      provenance.at("continuation_config_sha256").get<std::string>();
  data.pore_sha256 = provenance.at("pore_sha256").get<std::string>();
  data.throat_sha256 = provenance.at("throat_sha256").get<std::string>();
  data.connect_sha256 = provenance.at("connect_sha256").get<std::string>();
  if (!std::isfinite(data.time) || data.time < 0.0 || data.step_number < 0 ||
      !std::isfinite(data.cumulative_mole_ledger) ||
      !std::isfinite(data.cumulative_volume_ledger) ||
      !std::isfinite(data.cumulative_net_liquid_inflow) ||
      !std::isfinite(data.cumulative_boundary_component_outflow) ||
      !std::isfinite(data.cumulative_outlet_advective_component_outflow) ||
      !std::isfinite(data.cumulative_component_inflow) ||
      !std::isfinite(data.cumulative_component_outflow) ||
      !std::isfinite(data.cumulative_liquid_inflow) ||
      !std::isfinite(data.cumulative_liquid_outflow) ||
      !std::isfinite(data.next_dt) || data.next_dt <= 0.0 ||
      !std::isfinite(data.mass_transfer_multiplier) ||
      data.mass_transfer_multiplier <= 0.0 ||
      !std::isfinite(data.batch_individual_free_gas_fraction) ||
      !std::isfinite(data.batch_total_free_gas_fraction_per_step) ||
      data.batch_individual_free_gas_fraction <= 0.0 ||
      data.batch_individual_free_gas_fraction > 1.0 ||
      data.batch_total_free_gas_fraction_per_step <
          data.batch_individual_free_gas_fraction ||
      data.batch_total_free_gas_fraction_per_step > 1.0 ||
      !std::isfinite(data.last_gas_disappearance_time) ||
      data.consecutive_fast_steps < 0 || data.continuation_config_sha256.size() != 64 ||
      data.pore_sha256.size() != 64 || data.throat_sha256.size() != 64 ||
      data.connect_sha256.size() != 64) {
    throw std::runtime_error("checkpoint contains invalid continuation metadata");
  }
  data.state.reserve(expected_nodes);
  const bool has_extended = root.contains("state_extended_precision");
  if (has_extended && root.at("state_extended_precision").size() != expected_nodes)
    throw std::runtime_error("checkpoint extended-precision state count mismatch");
  std::size_t state_index = 0;
  for (const auto& row : root.at("state")) {
    if (!row.is_array() || row.size() != 6) throw std::runtime_error("malformed checkpoint state row");
    NodeState s;
    s.pl=row[0].get<double>();s.pg=row[1].get<double>();s.sw=row[2].get<double>();
    s.c=row[3].is_null()?std::numeric_limits<long double>::quiet_NaN():row[3].get<double>();
    s.active_gas=row[4].get<bool>();s.permanently_inactive=row[5].get<bool>();
    if (has_extended) {
      const auto& precise = root.at("state_extended_precision").at(state_index);
      if (!precise.is_array() || precise.size() != 4)
        throw std::runtime_error("malformed checkpoint extended-precision row");
      s.pl = std::stold(precise[0].get<std::string>());
      if (s.active_gas) s.pg = std::stold(precise[1].get<std::string>());
      s.sw = std::stold(precise[2].get<std::string>());
      const std::string precise_c = precise[3].get<std::string>();
      if (precise_c != "derived_from_adjacent_real") s.c = std::stold(precise_c);
    }
    const bool derived_c = row[3].is_null() && !s.active_gas && has_extended &&
        root.at("state_extended_precision").at(state_index)[3] == "derived_from_adjacent_real";
    if (!std::isfinite(s.pl) || !std::isfinite(s.sw) ||
        (!std::isfinite(s.c) && !derived_c) ||
        (s.active_gas && !std::isfinite(s.pg))) {
      throw std::runtime_error("checkpoint contains non-finite active state");
    }
    if (!s.active_gas) s.pg=std::numeric_limits<double>::quiet_NaN();
    s.concentration_derived = derived_c;
    data.state.push_back(s);
    ++state_index;
  }
  return data;
}

CheckpointProvenance make_checkpoint_provenance(const Config& config) {
  json continuation = json::parse(config.source_text);
  // Extending the requested end time or relocating output does not alter a
  // previously accepted trajectory.  All physics, initial state, timestep,
  // solver, acceptance and input-path settings remain signature-bound.
  continuation.at("time").erase("t_end_s");
  continuation.at("time").erase("stop_main_gas_volume_fraction");
  continuation.at("time").erase("equilibrium_tolerance_mol_m3");
  continuation.at("time").erase("dt_initial_s");
  continuation.at("time").erase("dt_min_s");
  continuation.at("time").erase("dt_max_s");
  continuation.at("time").erase("cut_factor");
  continuation.at("time").erase("growth_factor");
  continuation.at("time").erase("fast_newton_iterations");
  continuation.at("time").erase("max_retries_per_step");
  continuation.erase("time_step_control");
  auto& io = continuation.at("io");
  io.erase("output_directory");
  io.erase("checkpoint_interval_steps");
  io.erase("vtk_interval_steps");
  io.erase("vtk_output_times_s");
  io.erase("checkpoint_output_times_s");
  CheckpointProvenance result;
  result.continuation_config_sha256 = sha256_bytes(continuation.dump());
  result.pore_sha256 = sha256_file(config.io.pore_file);
  result.throat_sha256 = sha256_file(config.io.throat_file);
  result.connect_sha256 = sha256_file(config.io.connect_audit_file);
  return result;
}

std::string sha256_file(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot hash file: " + path.string());
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return sha256_bytes(buffer.str());
}

std::string utc_timestamp() {
  const auto now = std::chrono::system_clock::now();
  const std::time_t value = std::chrono::system_clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  gmtime_s(&tm, &value);
#else
  gmtime_r(&value, &tm);
#endif
  std::ostringstream out;
  out << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
  return out.str();
}

}  // namespace bubble
