#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "bubble/model.hpp"
#include "bubble/amgx_backend.hpp"
#include "bubble/pardiso_backend.hpp"

namespace bubble {
namespace {

using AD = fadbad::B<double>;

bool internal_node(const Network& network, std::size_t node) {
  return network.zero_distance_clusters[network.cluster_of_pore[node]]
             .boundary_flag == 0;
}

template <class Scalar>
struct FluxValues {
  Scalar q;
  Scalar f_adv;
  Scalar f_diff;
};

template <class Scalar>
FluxValues<Scalar> flux_values(const Edge& e, const Scalar& pla,
                               const Scalar& ca, const Scalar& swa,
                               const Scalar& plb, const Scalar& cb,
                               const Scalar& swb, bool upwind_a, bool diffuse,
                               const PhysicsConfig& p) {
  using std::pow;
  const Scalar sw_edge = 0.5 * (swa + swb);
  const double g0 = pi * std::pow(e.radius, 4) / (8.0 * p.viscosity * e.length);
  const Scalar conductance = g0 * pow(sw_edge, p.kr_exponent);
  const Scalar q = conductance * (pla - plb);
  const Scalar cup = upwind_a ? ca : cb;
  const Scalar adv = q * cup;
  const Scalar wet_area = e.area * pow(sw_edge, p.diffusion_area_exponent);
  const Scalar diff = diffuse ? (p.molecular_diffusivity / p.tortuosity) *
                                    wet_area * (ca - cb) / e.length
                              : Scalar(0.0);
  return {q, adv, diff};
}

double inf_norm_scaled(const Eigen::VectorXd& residual, const DofMap& map) {
  double result = 0.0;
  for (int i = 0; i < residual.size(); ++i) {
    if (!std::isfinite(residual[i])) return std::numeric_limits<double>::infinity();
    result = std::max(
        result, std::abs(residual[i] /
                         map.residual_scale[static_cast<std::size_t>(i)]));
  }
  return result;
}

double state_inf_scaled(const std::vector<NodeState>& state,
                        const DofMap& map) {
  double result = 0.0;
  for (std::size_t i = 0; i < state.size(); ++i) {
    const auto& d = map.node[i];
    if (d.pl >= 0)
      result = std::max(
          result, static_cast<double>(std::abs(
                      state[i].pl /
                      map.unknown_scale[static_cast<std::size_t>(d.pl)])));
    if (d.c >= 0)
      result = std::max(
          result,
          static_cast<double>(std::abs(
              state[i].c / map.unknown_scale[static_cast<std::size_t>(d.c)])));
    if (d.pg >= 0)
      result = std::max(
          result, static_cast<double>(std::abs(
                      state[i].pg /
                      map.unknown_scale[static_cast<std::size_t>(d.pg)])));
    if (d.sw >= 0)
      result = std::max(result, static_cast<double>(std::abs(state[i].sw)));
  }
  return result;
}

std::size_t active_count(const std::vector<NodeState>& state) {
  return static_cast<std::size_t>(
      std::count_if(state.begin(), state.end(),
                    [](const NodeState& s) { return s.active_gas; }));
}

bool all_finite(const Eigen::VectorXd& x) {
  for (int i = 0; i < x.size(); ++i) {
    if (!std::isfinite(x[i])) return false;
  }
  return true;
}

void dump_matrix_market(const std::filesystem::path& path,
                        const Eigen::SparseMatrix<double,Eigen::RowMajor>& input){
  Eigen::SparseMatrix<double,Eigen::RowMajor> matrix=input;matrix.makeCompressed();
  std::ofstream out(path);if(!out)throw std::runtime_error("cannot dump matrix: "+path.string());
  out<<"%%MatrixMarket matrix coordinate real general\n"<<matrix.rows()<<' '<<matrix.cols()<<' '<<matrix.nonZeros()<<'\n'<<std::setprecision(17);
  for(int row=0;row<matrix.outerSize();++row)
    for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(matrix,row);it;++it)
      out<<row+1<<' '<<it.col()+1<<' '<<it.value()<<'\n';
}

void dump_vector_market(const std::filesystem::path& path,const Eigen::VectorXd& values){
  std::ofstream out(path);if(!out)throw std::runtime_error("cannot dump vector: "+path.string());
  out<<"%%MatrixMarket matrix array real general\n"<<values.size()<<" 1\n"<<std::setprecision(17);
  for(int i=0;i<values.size();++i)out<<values[i]<<'\n';
}

}  // namespace

EdgeFlux evaluate_edge_flux(const Edge& edge, double pl_a, double c_a,
                            double sw_a, double pl_b, double c_b, double sw_b,
                            bool upwind_a, bool diffuse,
                            const PhysicsConfig& physics) {
  const auto f = flux_values(edge, pl_a, c_a, sw_a, pl_b, c_b, sw_b, upwind_a,
                             diffuse, physics);
  return {f.q, f.f_adv, f.f_diff, upwind_a};
}

Solver::Solver(Network network, Config config)
    : network_(std::move(network)), config_(std::move(config)) {
  network_.validate(config_);
  network_.build_degenerate_geometry(config_);
  if (config_.linear.backend == "amgx")
    amgx_solver_ = std::make_unique<AmgxLinearSolver>(config_.linear.amgx_config_file);
  if (config_.linear.backend == "pardiso")
    pardiso_solver_ = std::make_unique<PardisoLinearSolver>(
        config_.linear.pardiso_threads, config_.linear.relative_tolerance);
}

Solver::~Solver() = default;
Solver::Solver(Solver&&) noexcept = default;
Solver& Solver::operator=(Solver&&) noexcept = default;

void Solver::initialize() {
  state_.assign(network_.pores.size(), {});
  regularized_nodes_.clear();
  initial_regularization_moles_ = 0.0;
  initial_regularization_volume_ = 0.0;
  initial_threshold_moles_transferred_ = 0.0;
  cumulative_mole_ledger_ = 0.0;
  cumulative_volume_ledger_ = 0.0;
  cumulative_net_liquid_inflow_ = 0.0;
  cumulative_boundary_component_outflow_ = 0.0;
  cumulative_outlet_advective_component_outflow_ = 0.0;
  cumulative_component_inflow_ = 0.0;
  cumulative_component_outflow_ = 0.0;
  cumulative_liquid_inflow_ = 0.0;
  cumulative_liquid_outflow_ = 0.0;
  cumulative_interphase_transfer_ = 0.0;
  inlet_pressure_ = config_.boundary.pl_in;
  initial_free_gas_by_node_.assign(network_.pores.size(), 0.0L);
  initial_sw_by_node_.assign(network_.pores.size(), 1.0L);
  cumulative_transfer_by_node_.assign(network_.pores.size(), 0.0L);
  time_ = 0.0;
  step_number_ = 0;

  double xmin = network_.pores.front().x;
  double xmax = xmin;
  for (const auto& p : network_.pores) {
    xmin = std::min(xmin, p.x);
    xmax = std::max(xmax, p.x);
  }
  const double span = xmax - xmin;
  // The zero-distance quotient has one liquid pressure and concentration per
  // cluster.  Use its first (lowest input-ID) member as the deterministic
  // linear-pressure anchor; a boundary member instead imposes the reservoir
  // value on the whole cluster.
  for (const auto& cluster : network_.zero_distance_clusters) {
    const auto representative = cluster.members.front();
    double shared_pl = 0.0;
    double shared_c = config_.initial.concentration;
    if (cluster.boundary_flag == config_.boundary.inlet_flag) {
      shared_pl = config_.boundary.pl_in;
      shared_c = config_.boundary.c_in;
    } else if (cluster.boundary_flag == config_.boundary.outlet_flag) {
      shared_pl = config_.boundary.pl_out;
      shared_c = config_.boundary.c_backflow;
    } else {
      int adjacent_virtual_flag = 0;
      for (const auto node : cluster.members) {
        for (const auto edge_index : network_.adjacency[node]) {
          const auto& edge = network_.edges[edge_index];
          const auto neighbor = edge.a == node ? edge.b : edge.a;
          if (network_.pores[neighbor].source_type != 2) continue;
          const int flag = network_.pores[neighbor].boundary_flag;
          if (adjacent_virtual_flag != 0 && adjacent_virtual_flag != flag)
            throw std::runtime_error(
                "real cluster touches conflicting virtual boundaries");
          adjacent_virtual_flag = flag;
        }
      }
      if (adjacent_virtual_flag == config_.boundary.inlet_flag) {
        shared_pl = config_.boundary.pl_in;
      } else if (adjacent_virtual_flag == config_.boundary.outlet_flag) {
        shared_pl = config_.boundary.pl_out;
      } else {
        const double xi =
            span > 0.0 ? (network_.pores[representative].x - xmin) / span : 0.5;
        shared_pl = config_.boundary.pl_in +
                    xi * (config_.boundary.pl_out - config_.boundary.pl_in);
      }
    }
    for (const auto node : cluster.members) {
      state_[node].pl = shared_pl;
      state_[node].c = shared_c;
      state_[node].sw = 1.0;
    }
    if (cluster.boundary_flag != 0) continue;

    double cluster_dissolved = 0.0;
    double cluster_liquid_volume = 0.0;
    double cluster_retired_moles = 0.0;
    for (const auto node : cluster.members) {
      const auto& p = network_.pores[node];
      auto& s = state_[node];
      double sw = p.input_sw;
      if (sw <= config_.active_set.sw_min) {
        // Keep a small but meaningful interior distance from the active lower
        // bound. One binary64 ULP is smaller than the configured minimum line
        // search step on the physical-micrometre network and can make an
        // otherwise valid first Newton direction impossible to test.
        const double corrected =
            std::min(100.0 * config_.active_set.sw_min,
                     0.5 * (config_.active_set.sw_min + 1.0 -
                            config_.active_set.sg_off));
        initial_regularization_volume_ += p.volume * (corrected - sw);
        initial_regularization_moles_ += p.volume * (corrected - sw) * shared_c;
        sw = corrected;
        regularized_nodes_.push_back(node);
      }
      if (p.volume == 0.0) {
        s.sw = 1.0;
        continue;
      }

      const double dissolved_before_retirement = p.volume * sw * shared_c;
      if (1.0 - sw <= config_.active_set.sg_off) {
        if (sw < 1.0) {
          const double pg =
              shared_pl + capillary_pressure(sw, p.radius, config_.physics);
          if (!std::isfinite(pg) || config_.physics.background_abs_Pa + pg <= 0.0)
            throw std::runtime_error(
                "initial absolute bubble pressure is non-positive at pore " +
                std::to_string(p.id));
          cluster_retired_moles +=
              free_gas_moles(pg, sw, p.volume, config_.physics);
          initial_regularization_volume_ += p.volume * (1.0 - sw);
          s.permanently_inactive = true;
          if (std::find(regularized_nodes_.begin(), regularized_nodes_.end(),
                        node) == regularized_nodes_.end())
            regularized_nodes_.push_back(node);
        }
        s.sw = 1.0;
      } else {
        s.sw = sw;
        s.pg = shared_pl + capillary_pressure(s.sw, p.radius, config_.physics);
        if (!std::isfinite(s.pg) || config_.physics.background_abs_Pa + s.pg <= 0.0)
          throw std::runtime_error(
              "initial absolute bubble pressure is non-positive at pore " +
              std::to_string(p.id));
        const double ng = static_cast<double>(
            free_gas_moles(s.pg, s.sw, p.volume, config_.physics));
        if (ng > config_.active_set.ng_off) {
          s.active_gas = true;
        } else {
          cluster_retired_moles += ng;
          initial_regularization_volume_ += p.volume * (1.0 - sw);
          s.sw = 1.0;
          s.pg = std::numeric_limits<double>::quiet_NaN();
          s.permanently_inactive = true;
          if (std::find(regularized_nodes_.begin(), regularized_nodes_.end(),
                        node) == regularized_nodes_.end())
            regularized_nodes_.push_back(node);
        }
      }
      cluster_dissolved += dissolved_before_retirement;
      cluster_liquid_volume += static_cast<double>(p.volume * s.sw);
    }
    initial_threshold_moles_transferred_ += cluster_retired_moles;
    if (cluster_liquid_volume > 0.0) {
      shared_c =
          (cluster_dissolved + cluster_retired_moles) / cluster_liquid_volume;
      for (const auto node : cluster.members) state_[node].c = shared_c;
    }
  }

  initial_total_moles_ = 0.0;
  internal_volume_ = 0.0;
  initial_liquid_volume_ = 0.0;
  initial_main_gas_volume_ = 0.0;
  for (std::size_t i = 0; i < state_.size(); ++i) {
    if (!internal_node(network_, i)) continue;
    const auto& p = network_.pores[i];
    const auto& s = state_[i];
    internal_volume_ += p.volume;
    initial_liquid_volume_ += static_cast<double>(p.volume * s.sw);
    initial_total_moles_ += static_cast<double>(p.volume * s.sw * s.c);
    if (s.active_gas)
      initial_total_moles_ += static_cast<double>(
          free_gas_moles(s.pg, s.sw, p.volume, config_.physics));
    if (s.active_gas)
      initial_free_gas_by_node_[i] =
          free_gas_moles(s.pg, s.sw, p.volume, config_.physics);
    if (s.active_gas && i < network_.component_id.size() &&
        network_.component_id[i] == network_.main_flow_component_id)
      initial_main_gas_volume_ += static_cast<double>(p.volume * (1.0 - s.sw));
  }
  if (config_.boundary.species_mode == "copy_adjacent_concentration") {
    for (std::size_t i = 0; i < state_.size(); ++i) {
      if (network_.pores[i].source_type != 2) continue;
      const auto edge_index = network_.adjacency[i].front();
      const auto& edge = network_.edges[edge_index];
      const auto neighbor = edge.a == i ? edge.b : edge.a;
      state_[i].c = state_[neighbor].c;
      state_[i].concentration_derived = true;
    }
  }
  for (std::size_t i = 0; i < state_.size(); ++i)
    initial_sw_by_node_[i] = state_[i].sw;
}

DofMap Solver::make_dof_map(const std::vector<NodeState>& state) const {
  if (state.size() != network_.pores.size())
    throw std::runtime_error("state/network size mismatch");
  DofMap map;
  map.node.resize(state.size());
  map.cluster.resize(network_.zero_distance_clusters.size());
  for (const auto& cluster : network_.zero_distance_clusters) {
    auto& shared = map.cluster[cluster.id];
    if (cluster.boundary_flag == 0) {
      shared.pl = shared.rl = map.size++;
      map.unknown_scale.push_back(config_.scaling.pressure);
      map.residual_scale.push_back(config_.scaling.volume_rate);
      shared.c = shared.rd = map.size++;
      map.unknown_scale.push_back(config_.scaling.concentration);
      map.residual_scale.push_back(config_.scaling.molar_rate);
    }
    for (const auto node : cluster.members) {
      auto& d = map.node[node];
      d.pl = shared.pl;
      d.rl = shared.rl;
      d.c = shared.c;
      d.rd = shared.rd;
      if (!state[node].active_gas) continue;
      d.pg = d.rg = map.size++;
      map.unknown_scale.push_back(config_.scaling.gas_pressure > 0.0 ?
          config_.scaling.gas_pressure : config_.scaling.pressure);
      map.residual_scale.push_back(config_.scaling.molar_rate);
      d.sw = d.rc = map.size++;
      map.unknown_scale.push_back(1.0);
      map.residual_scale.push_back(config_.scaling.pressure);
    }
  }
  // A closed, entirely liquid component has one redundant volume equation:
  // pressure is determined only up to an additive constant.  Replace one Rl
  // row by a pressure gauge during assembly.  Components connected to a
  // reservoir, or containing gas whose EOS fixes absolute pressure, need no
  // gauge.
  std::vector<bool> visited(network_.zero_distance_clusters.size(), false);
  for (std::size_t root = 0; root < network_.zero_distance_clusters.size();
       ++root) {
    if (visited[root]) continue;
    std::vector<std::size_t> stack{root};
    visited[root] = true;
    bool has_pressure_boundary = false;
    bool has_active_gas = false;
    std::size_t gauge = network_.zero_distance_clusters.size();
    while (!stack.empty()) {
      const auto cluster_id = stack.back();
      stack.pop_back();
      const auto& cluster = network_.zero_distance_clusters[cluster_id];
      if (cluster.boundary_flag == 0) {
        if (gauge == network_.zero_distance_clusters.size()) gauge = cluster_id;
        for (const auto node : cluster.members)
          has_active_gas = has_active_gas || state[node].active_gas;
      } else {
        has_pressure_boundary = true;
      }
      for (const auto edge_index : cluster.external_positive_edges) {
        const auto& edge = network_.edges[edge_index];
        const auto ca = network_.cluster_of_pore[edge.a];
        const auto cb = network_.cluster_of_pore[edge.b];
        const auto neighbor = ca == cluster_id ? cb : ca;
        if (!visited[neighbor]) {
          visited[neighbor] = true;
          stack.push_back(neighbor);
        }
      }
    }
    if (!has_pressure_boundary && !has_active_gas &&
        gauge != network_.zero_distance_clusters.size()) {
      map.cluster[gauge].pressure_gauge = true;
      map.node[network_.zero_distance_clusters[gauge].members.front()]
          .pressure_gauge = true;
    }
  }
  if (config_.boundary.inlet_type == "flow") {
    map.inlet_pressure = map.inlet_flow_row = map.size++;
    map.unknown_scale.push_back(config_.scaling.pressure);
    map.residual_scale.push_back(config_.scaling.volume_rate);
  }
  return map;
}

Assembly Solver::assemble(const std::vector<NodeState>& trial,
                          const std::vector<NodeState>& old, double dt,
                          const DofMap& map, const std::vector<bool>& upwind_a,
                          const std::vector<bool>& transfer_active,
                          bool build_jacobian) const {
  if (trial.size() != network_.pores.size() || old.size() != trial.size() ||
      upwind_a.size() != network_.edges.size() ||
      transfer_active.size() != trial.size()) {
    throw std::runtime_error("assemble input size mismatch");
  }
  Assembly out;
  out.residual = Eigen::VectorXd::Zero(map.size);
  std::vector<long double> residual_extended(static_cast<std::size_t>(map.size),
                                             0.0L);
  out.edge_fluxes.resize(network_.edges.size());
  out.transfer.assign(trial.size(), 0.0);
  std::vector<Eigen::Triplet<double>> entries;
  if (build_jacobian) entries.reserve(static_cast<std::size_t>(map.size) * 12);

  for (std::size_t i = 0; i < trial.size(); ++i) {
    const auto cluster_id = network_.cluster_of_pore[i];
    const auto& shared = map.cluster[cluster_id];
    if (shared.rl < 0 && !trial[i].active_gas) continue;
    const auto& p = network_.pores[i];
    const auto& now = trial[i];
    const auto& before = old[i];
    const auto& d = map.node[i];
    const bool gas = now.active_gas;

    AD pl = static_cast<double>(now.pl), c = static_cast<double>(now.c);
    AD pg = gas ? static_cast<double>(now.pg) : 0.0;
    AD sw = gas ? static_cast<double>(now.sw) : 1.0;
    const long double sw_old = before.active_gas ? before.sw : 1.0L;
    const double c_old = static_cast<double>(before.c);
    AD ntr = 0.0;
    long double ntr_value = 0.0L;
    if (gas && transfer_active[i]) {
      ntr = config_.physics.mass_transfer_multiplier *
            config_.physics.mass_transfer_coefficient *
            interface_area(sw, p.volume) * (config_.physics.henry_cp * (config_.physics.background_abs_Pa + pg) - c);
      const long double gas_volume = static_cast<long double>(p.volume) *
                                     (1.0L - static_cast<long double>(now.sw));
      const long double area =
          std::pow(36.0L * static_cast<long double>(pi), 1.0L / 3.0L) *
          std::pow(gas_volume, 2.0L / 3.0L);
      ntr_value =
          static_cast<long double>(config_.physics.mass_transfer_multiplier) *
          static_cast<long double>(config_.physics.mass_transfer_coefficient) *
          area *
          (static_cast<long double>(config_.physics.henry_cp) *
               (config_.physics.background_abs_Pa + now.pg) -
           static_cast<long double>(now.c));
    }
    const bool gauge_representative =
        shared.pressure_gauge &&
        i == network_.zero_distance_clusters[cluster_id].members.front();
    AD rl = 0.0;
    if (shared.rl >= 0) {
      if (gauge_representative) {
        rl = config_.scaling.volume_rate * (pl - before.pl) /
             config_.scaling.pressure;
      } else if (!shared.pressure_gauge) {
        rl = p.volume * (sw - sw_old) / dt;
      }
    }
    // Algebraically equivalent increment forms avoid subtracting two nearly
    // equal O(V*C) inventories after a very small timestep.
    AD rd = p.volume * (sw * (c - c_old) + c_old * (sw - sw_old)) / dt - ntr;
    AD rg = 0.0, rc = 0.0;
    if (gas) {
      const double coefficient =
          p.volume / (config_.physics.compressibility * gas_constant *
                      config_.physics.temperature);
      rg = coefficient *
               ((1.0 - sw) * (pg - before.pg) - (config_.physics.background_abs_Pa + before.pg) * (sw - sw_old)) /
               dt +
           ntr;
      rc = pg - pl - capillary_pressure(sw, p.radius, config_.physics);
    }
    if (shared.rl >= 0) {
      const long double rl_value =
          gauge_representative
              ? static_cast<long double>(config_.scaling.volume_rate) *
                    (now.pl - before.pl) / config_.scaling.pressure
              : (shared.pressure_gauge
                     ? 0.0L
                     : static_cast<long double>(p.volume) *
                           (now.sw - static_cast<long double>(sw_old)) / dt);
      const long double rd_value =
          static_cast<long double>(p.volume) *
              (now.sw * (now.c - before.c) +
               before.c * (now.sw - static_cast<long double>(sw_old))) /
              dt -
          ntr_value;
      residual_extended[static_cast<std::size_t>(shared.rl)] += rl_value;
      residual_extended[static_cast<std::size_t>(shared.rd)] += rd_value;
    }
    out.transfer[i] = ntr.x();
    if (gas) {
      const long double coefficient =
          static_cast<long double>(p.volume) /
          (static_cast<long double>(config_.physics.compressibility) *
           static_cast<long double>(gas_constant) *
           static_cast<long double>(config_.physics.temperature));
      const long double rg_value =
          coefficient *
              ((1.0L - static_cast<long double>(now.sw)) *
                   (static_cast<long double>(now.pg) -
                    static_cast<long double>(before.pg)) -
               (config_.physics.background_abs_Pa + before.pg) *
                   (static_cast<long double>(now.sw) -
                    static_cast<long double>(sw_old))) /
              static_cast<long double>(dt) +
          ntr_value;
      residual_extended[static_cast<std::size_t>(d.rg)] += rg_value;
      const long double rc_value =
          now.pg - now.pl -
          capillary_pressure(now.sw, p.radius, config_.physics);
      residual_extended[static_cast<std::size_t>(d.rc)] += rc_value;
    }
    if (build_jacobian) {
      std::vector<std::pair<int, AD*>> equations;
      if (shared.rl >= 0) {
        equations.emplace_back(shared.rl, &rl);
        equations.emplace_back(shared.rd, &rd);
      }
      if (gas) {
        equations.emplace_back(d.rg, &rg);
        equations.emplace_back(d.rc, &rc);
      }
      const auto outputs = static_cast<unsigned int>(equations.size());
      for (unsigned int r = 0; r < outputs; ++r)
        equations[r].second->diff(r, outputs);
      // FADBAD++ B<T> propagates a shared reverse graph when its final
      // non-leaf reference is released. ntr is shared by Rd and Rg, so drop
      // that retained intermediate after seeding all output components.
      ntr = 0.0;
      const std::array<std::pair<int, AD*>, 4> vars{
          {{d.pl, &pl}, {d.c, &c}, {d.pg, &pg}, {d.sw, &sw}}};
      for (unsigned int r = 0; r < outputs; ++r) {
        for (const auto& variable : vars) {
          if (variable.first < 0) continue;
          const double derivative = variable.second->d(r);
          if (derivative != 0.0)
            entries.emplace_back(equations[r].first, variable.first,
                                 derivative);
        }
      }
    }
  }

  for (std::size_t k = 0; k < network_.edges.size(); ++k) {
    const auto& e = network_.edges[k];
    if (e.length == 0.0) continue;
    const auto ca_cluster = network_.cluster_of_pore[e.a];
    const auto cb_cluster = network_.cluster_of_pore[e.b];
    const int flag_a =
        network_.zero_distance_clusters[ca_cluster].boundary_flag;
    const int flag_b =
        network_.zero_distance_clusters[cb_cluster].boundary_flag;
    const bool ia = flag_a == 0, ib = flag_b == 0;
    if (!ia && !ib) continue;
    const bool outlet_a = flag_a == config_.boundary.outlet_flag;
    const bool outlet_b = flag_b == config_.boundary.outlet_flag;
    const bool boundary_edge = !ia || !ib;
    const bool copy_adjacent =
        config_.boundary.species_mode == "copy_adjacent_concentration";
    const bool diffuse =
        copy_adjacent ? !boundary_edge : !outlet_a && !outlet_b;
    AD inlet_pl_ad = inlet_pressure_;

    AD pla =
        flag_a == config_.boundary.inlet_flag
            ? (config_.boundary.inlet_type=="flow" ? inlet_pl_ad : config_.boundary.pl_in)
            : (flag_a == config_.boundary.outlet_flag ? config_.boundary.pl_out
                                                      : trial[e.a].pl);
    AD ca = flag_a == config_.boundary.inlet_flag
                ? config_.boundary.c_in
                : (flag_a == config_.boundary.outlet_flag
                       ? config_.boundary.c_backflow
                       : trial[e.a].c);
    AD swa = ia && trial[e.a].active_gas ? trial[e.a].sw : 1.0;
    AD plb =
        flag_b == config_.boundary.inlet_flag
            ? (config_.boundary.inlet_type=="flow" ? inlet_pl_ad : config_.boundary.pl_in)
            : (flag_b == config_.boundary.outlet_flag ? config_.boundary.pl_out
                                                      : trial[e.b].pl);
    AD cb = flag_b == config_.boundary.inlet_flag
                ? config_.boundary.c_in
                : (flag_b == config_.boundary.outlet_flag
                       ? config_.boundary.c_backflow
                       : trial[e.b].c);
    AD swb = ib && trial[e.b].active_gas ? trial[e.b].sw : 1.0;
    if (copy_adjacent && boundary_edge) {
      if (!ia && ib) ca = cb;
      if (ia && !ib) cb = ca;
    }
    auto flux = flux_values(e, pla, ca, swa, plb, cb, swb, upwind_a[k], diffuse,
                            config_.physics);
    AD total_f = flux.f_adv + flux.f_diff;
    const double q = flux.q.x();
    const double fa = flux.f_adv.x();
    const double fd = flux.f_diff.x();
    out.edge_fluxes[k] = {q, fa, fd, upwind_a[k]};
    const bool inlet_edge = config_.boundary.inlet_type=="flow" &&
        (flag_a==config_.boundary.inlet_flag || flag_b==config_.boundary.inlet_flag);
    if (inlet_edge && map.inlet_flow_row>=0)
      residual_extended[static_cast<std::size_t>(map.inlet_flow_row)] +=
          (flag_a==config_.boundary.inlet_flag ? q : -q);
    if (ia && !map.cluster[ca_cluster].pressure_gauge) {
      residual_extended[static_cast<std::size_t>(map.cluster[ca_cluster].rl)] +=
          q;
      residual_extended[static_cast<std::size_t>(map.cluster[ca_cluster].rd)] +=
          fa + fd;
    } else if (ia) {
      residual_extended[static_cast<std::size_t>(map.cluster[ca_cluster].rd)] +=
          fa + fd;
    }
    if (ib && !map.cluster[cb_cluster].pressure_gauge) {
      residual_extended[static_cast<std::size_t>(map.cluster[cb_cluster].rl)] -=
          q;
      residual_extended[static_cast<std::size_t>(map.cluster[cb_cluster].rd)] -=
          fa + fd;
    } else if (ib) {
      residual_extended[static_cast<std::size_t>(map.cluster[cb_cluster].rd)] -=
          fa + fd;
    }
    if (build_jacobian) {
      flux.q.diff(0, 2);
      total_f.diff(1, 2);
      // advective and diffusive outputs retain the graph used by total_f.
      // Release them before reading derivatives from the leaf variables.
      flux.f_adv = 0.0;
      flux.f_diff = 0.0;
      const std::array<std::pair<int, AD*>, 6> vars{
          {{ia ? map.node[e.a].pl : (flag_a==config_.boundary.inlet_flag&&config_.boundary.inlet_type=="flow"?map.inlet_pressure:-1), &pla},
           {ia ? map.node[e.a].c : -1, &ca},
           {ia ? map.node[e.a].sw : -1, &swa},
           {ib ? map.node[e.b].pl : (flag_b==config_.boundary.inlet_flag&&config_.boundary.inlet_type=="flow"?map.inlet_pressure:-1), &plb},
           {ib ? map.node[e.b].c : -1, &cb},
           {ib ? map.node[e.b].sw : -1, &swb}}};
      for (const auto& variable : vars) {
        if (variable.first < 0) continue;
        const double dq = variable.second->d(0);
        const double df = variable.second->d(1);
        if(inlet_edge && map.inlet_flow_row>=0 && dq!=0.0)
          entries.emplace_back(map.inlet_flow_row, variable.first,
              flag_a==config_.boundary.inlet_flag ? dq : -dq);
        if (ia && !map.cluster[ca_cluster].pressure_gauge) {
          if (dq != 0.0)
            entries.emplace_back(map.cluster[ca_cluster].rl, variable.first,
                                 dq);
        }
        if (ia) {
          if (df != 0.0)
            entries.emplace_back(map.cluster[ca_cluster].rd, variable.first,
                                 df);
        }
        if (ib && !map.cluster[cb_cluster].pressure_gauge) {
          if (dq != 0.0)
            entries.emplace_back(map.cluster[cb_cluster].rl, variable.first,
                                 -dq);
        }
        if (ib) {
          if (df != 0.0)
            entries.emplace_back(map.cluster[cb_cluster].rd, variable.first,
                                 -df);
        }
      }
    }
  }

  if (config_.boundary.inlet_type == "flow" && map.inlet_flow_row >= 0)
    residual_extended[static_cast<std::size_t>(map.inlet_flow_row)] -=
        config_.boundary.inlet_flow_m3_s;

  for (int row = 0; row < map.size; ++row) {
    out.residual[row] =
        static_cast<double>(residual_extended[static_cast<std::size_t>(row)]);
  }
  out.residual_norms = {{0.0, 0.0, 0.0, 0.0}};
  out.residual_node_ids = {{0, 0, 0, 0}};
  for (const auto& cluster : network_.zero_distance_clusters) {
    const auto& shared = map.cluster[cluster.id];
    if (shared.rl < 0) continue;
    const auto node_id = network_.pores[cluster.members.front()].id;
    const std::array<std::pair<int, double>, 2> shared_rows{
        {{shared.rl, config_.scaling.volume_rate},
         {shared.rd, config_.scaling.molar_rate}}};
    for (std::size_t block = 0; block < shared_rows.size(); ++block) {
      const double value = std::abs(out.residual[shared_rows[block].first] /
                                    shared_rows[block].second);
      if (value > out.residual_norms[block]) {
        out.residual_norms[block] = value;
        out.residual_node_ids[block] = node_id;
      }
    }
  }
  for (std::size_t i = 0; i < map.node.size(); ++i) {
    const auto& d = map.node[i];
    if (d.rg < 0) continue;
    const std::array<std::pair<int, double>, 2> gas_rows{
        {{d.rg, config_.scaling.molar_rate}, {d.rc, config_.scaling.pressure}}};
    for (std::size_t offset = 0; offset < gas_rows.size(); ++offset) {
      const std::size_t block = offset + 2;
      const double value = std::abs(out.residual[gas_rows[offset].first] /
                                    gas_rows[offset].second);
      if (value > out.residual_norms[block]) {
        out.residual_norms[block] = value;
        out.residual_node_ids[block] = network_.pores[i].id;
      }
    }
  }
  if (build_jacobian) {
    out.jacobian.resize(map.size, map.size);
    out.jacobian.setFromTriplets(entries.begin(), entries.end(),
                                 [](double a, double b) { return a + b; });
    out.jacobian.makeCompressed();
  }
  return out;
}

Balance Solver::compute_balance(const std::vector<NodeState>& old,
                                const std::vector<NodeState>& current,
                                const Assembly& assembly, double dt,
                                double mole_ledger,
                                double volume_ledger) const {
  Balance b;
  for (const double transfer : assembly.transfer) b.transfer_rate += transfer;
  double old_total = 0.0, old_liquid_volume = 0.0, new_liquid_volume = 0.0;
  for (std::size_t i = 0; i < network_.pores.size(); ++i) {
    if (!internal_node(network_, i)) continue;
    const auto& p = network_.pores[i];
    old_total += static_cast<double>(p.volume * old[i].sw * old[i].c);
    if(current[i].pl < b.internal_pressure_min) b.pl_min_id=p.id;
    if(current[i].pl > b.internal_pressure_max) b.pl_max_id=p.id;
    if(current[i].active_gas){
      if(current[i].pg < b.pg_min){b.pg_min=static_cast<double>(current[i].pg);b.pg_min_id=p.id;}
      if(current[i].pg > b.pg_max){b.pg_max=static_cast<double>(current[i].pg);b.pg_max_id=p.id;}
    }
    b.total_dissolved +=
        static_cast<double>(p.volume * current[i].sw * current[i].c);
    if (old[i].active_gas)
      old_total += static_cast<double>(
          free_gas_moles(old[i].pg, old[i].sw, p.volume, config_.physics));
    if (current[i].active_gas)
      b.total_free += static_cast<double>(free_gas_moles(
          current[i].pg, current[i].sw, p.volume, config_.physics));
    const double gas_volume =
        current[i].active_gas
            ? static_cast<double>(p.volume * (1.0 - current[i].sw))
            : 0.0;
    const double free_moles =
        current[i].active_gas
            ? static_cast<double>(free_gas_moles(current[i].pg, current[i].sw,
                                                 p.volume, config_.physics))
            : 0.0;
    const int connection_class = i < network_.boundary_connection_class.size()
                                     ? network_.boundary_connection_class[i]
                                     : 0;
    const bool main_flow =
        i < network_.component_id.size() &&
        network_.component_id[i] == network_.main_flow_component_id;
    b.gas_volume_total += gas_volume;
    if (connection_class == 0) {
      b.free_boundaryless += free_moles;
      b.gas_volume_boundaryless += gas_volume;
      b.active_boundaryless += static_cast<std::size_t>(current[i].active_gas);
    } else {
      b.free_boundary_connected += free_moles;
      b.gas_volume_boundary_connected += gas_volume;
      b.active_boundary_connected +=
          static_cast<std::size_t>(current[i].active_gas);
    }
    if (main_flow) {
      b.free_main_flow += free_moles;
      b.gas_volume_main_flow += gas_volume;
      b.active_main_flow += static_cast<std::size_t>(current[i].active_gas);
    }
    if (current[i].active_gas)
      b.maximum_dissolution_driving_force = std::max(
          b.maximum_dissolution_driving_force,
          static_cast<double>(config_.physics.henry_cp * (config_.physics.background_abs_Pa + current[i].pg) -
                              current[i].c));
    b.internal_pressure_min =
        std::min(b.internal_pressure_min, static_cast<double>(current[i].pl));
    b.internal_pressure_max =
        std::max(b.internal_pressure_max, static_cast<double>(current[i].pl));
    old_liquid_volume += static_cast<double>(p.volume * old[i].sw);
    new_liquid_volume += static_cast<double>(p.volume * current[i].sw);
  }
  for (std::size_t k = 0; k < network_.edges.size(); ++k) {
    const auto& e = network_.edges[k];
    if (e.length == 0.0) continue;
    const bool ia = internal_node(network_, e.a);
    const bool ib = internal_node(network_, e.b);
    if (ia == ib) continue;
    const double sign = ia ? 1.0 : -1.0;
    const double liquid_out = sign * assembly.edge_fluxes[k].q;
    const double advective_out = sign * assembly.edge_fluxes[k].advective;
    const double diffusive_out = sign * assembly.edge_fluxes[k].diffusive;
    b.boundary_liquid_outflow += liquid_out;
    b.component_inflow += std::max(0.0, -advective_out - diffusive_out);
    b.component_outflow += std::max(0.0, advective_out + diffusive_out);
    b.liquid_inflow += std::max(0.0, -liquid_out);
    b.liquid_outflow += std::max(0.0, liquid_out);
    b.boundary_advective_outflow += advective_out;
    b.boundary_diffusive_outflow += diffusive_out;
    b.boundary_molar_outflow += advective_out + diffusive_out;
    const auto boundary_node = ia ? e.b : e.a;
    const int flag = network_.pores[boundary_node].boundary_flag;
    if (flag == config_.boundary.inlet_flag) {
      b.inlet_liquid_inflow -= liquid_out;
      b.inlet_advective_component_inflow -= advective_out;
      if (liquid_out > 0.0) {
        ++b.wrong_direction_boundary_edges;
        ++b.inlet_reverse_flow_edges;
        b.inlet_reverse_flow += liquid_out;
      }
    } else if (flag == config_.boundary.outlet_flag) {
      b.outlet_liquid_outflow += liquid_out;
      b.outlet_advective_component_outflow += advective_out;
      if (liquid_out < 0.0) {
        ++b.wrong_direction_boundary_edges;
        ++b.outlet_backflow_edges;
        b.outlet_backflow -= liquid_out;
      }
    }
  }
  b.net_liquid_inflow = b.inlet_liquid_inflow - b.outlet_liquid_outflow;
  if(b.pg_min_id==0){b.pg_min=0.0;b.pg_max=0.0;}
  b.throughflow_rate = 0.5 * (b.inlet_liquid_inflow + b.outlet_liquid_outflow);
  b.liquid_volume_change = new_liquid_volume - old_liquid_volume;
  b.cumulative_liquid_volume_change =
      new_liquid_volume - initial_liquid_volume_;
  const double new_total = b.total_dissolved + b.total_free;
  const double mole_raw = new_total - old_total + dt * b.boundary_molar_outflow;
  const double volume_raw =
      new_liquid_volume - old_liquid_volume + dt * b.boundary_liquid_outflow;
  const double mole_den =
      std::max(initial_total_moles_, config_.scaling.mole_total);
  const double volume_den =
      std::max(internal_volume_, config_.scaling.volume_total);
  b.epsilon_moles_equation = std::abs(mole_raw) / mole_den;
  b.epsilon_volume_equation = std::abs(volume_raw) / volume_den;
  b.epsilon_moles_adjusted = std::abs(mole_raw - mole_ledger) / mole_den;
  b.epsilon_volume_adjusted = std::abs(volume_raw - volume_ledger) / volume_den;
  return b;
}

StepReport Solver::attempt_step(double dt) {
  StepReport report;
  const double accepted_inlet_pressure = inlet_pressure_;
  double trial_inlet_pressure = inlet_pressure_;
  bool commit_inlet_pressure = false;
  struct InletPressureRollback {
    double& value; double accepted; bool& commit;
    ~InletPressureRollback(){ if(!commit) value=accepted; }
  } inlet_rollback{inlet_pressure_,accepted_inlet_pressure,commit_inlet_pressure};
  report.dt = dt;
  report.active_before = active_count(state_);
  if (!(dt > 0.0) || state_.empty()) {
    report.failure = "invalid dt or uninitialized state";
    return report;
  }
  const auto old = state_;
  auto trial = old;
  const DofMap map = make_dof_map(trial);
  double last_update = 0.0;
  bool converged = false;
  Assembly final;

  for (int iteration = 0; iteration <= config_.nonlinear.max_iterations;
       ++iteration) {
    std::vector<bool> upwind(network_.edges.size(), true);
    for (std::size_t k = 0; k < network_.edges.size(); ++k) {
      const auto& e = network_.edges[k];
      if (e.length == 0.0) continue;
      const bool ia = internal_node(network_, e.a);
      const bool ib = internal_node(network_, e.b);
      if (!ia && !ib) continue;
      const double pla = static_cast<double>(trial[e.a].pl);
      const double plb = static_cast<double>(trial[e.b].pl);
      const double swa = ia && trial[e.a].active_gas
                             ? static_cast<double>(trial[e.a].sw)
                             : 1.0;
      const double swb = ib && trial[e.b].active_gas
                             ? static_cast<double>(trial[e.b].sw)
                             : 1.0;
      const auto q = flux_values(e, pla, 0.0, swa, plb, 0.0, swb, true, false,
                                 config_.physics)
                         .q;
      upwind[k] = q >= 0.0;
    }
    std::vector<bool> transfer_active(trial.size(), false);
    for (std::size_t i = 0; i < trial.size(); ++i) {
      transfer_active[i] = trial[i].active_gas &&
                           config_.physics.henry_cp * (config_.physics.background_abs_Pa + trial[i].pg) - trial[i].c >
                               config_.active_set.transfer_switch_tolerance;
    }
    const auto assembly_started = std::chrono::steady_clock::now();
    inlet_pressure_ = trial_inlet_pressure;
    Assembly assembled =
        assemble(trial, old, dt, map, upwind, transfer_active, true);
    report.assembly_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - assembly_started).count();
    const double norm = inf_norm_scaled(assembled.residual, map);
    report.residual_norm = norm;
    report.newton_iterations = iteration;
    if (!std::isfinite(norm)) {
      report.failure = "non-finite nonlinear residual";
      return report;
    }
    const double state_norm = state_inf_scaled(trial, map);
    const double update_limit =
        config_.nonlinear.update_tolerance * (1.0 + state_norm);
    report.update_norm = last_update;
    report.state_norm = state_norm;
    report.relative_update_norm = last_update / (1.0 + state_norm);
    report.update_tolerance_limit = update_limit;
    report.update_criterion_ratio = last_update / update_limit;
    report.residual_criterion_passed =
        norm <= config_.nonlinear.residual_tolerance;
    report.update_criterion_passed = last_update <= update_limit;
    if (report.residual_criterion_passed && report.update_criterion_passed) {
      final = std::move(assembled);
      converged = true;
      break;
    }
    if (iteration == config_.nonlinear.max_iterations) break;

    Eigen::SparseMatrix<double, Eigen::RowMajor> scaled = assembled.jacobian;
    for (int row = 0; row < scaled.outerSize(); ++row) {
      for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(
               scaled, row);
           it; ++it) {
        it.valueRef() *= map.unknown_scale[static_cast<std::size_t>(it.col())] /
                         map.residual_scale[static_cast<std::size_t>(row)];
      }
    }
    Eigen::VectorXd rhs(map.size);
    for (int i = 0; i < map.size; ++i)
      rhs[i] = -assembled.residual[i] /
               map.residual_scale[static_cast<std::size_t>(i)];
    if(config_.linear.backend=="amgx"&&!linear_system_dumped_){
      if(const char* requested=std::getenv("BUBBLE_AMGX_DUMP_PREFIX")){
        const std::filesystem::path prefix(requested);
        if(!prefix.parent_path().empty())std::filesystem::create_directories(prefix.parent_path());
        dump_matrix_market(prefix.string()+"_jacobian_unscaled.mtx",assembled.jacobian);
        dump_matrix_market(prefix.string()+"_matrix_scaled.mtx",scaled);
        dump_vector_market(prefix.string()+"_rhs_scaled.mtx",rhs);
        dump_vector_market(prefix.string()+"_initial_solution.mtx",Eigen::VectorXd::Zero(map.size));
        std::ofstream scales(prefix.string()+"_physical_scaling.tsv");
        scales<<"dof_index_zero_based\tunknown_scale\tresidual_scale\n"<<std::setprecision(17);
        for(int i=0;i<map.size;++i)scales<<i<<'\t'<<map.unknown_scale[static_cast<std::size_t>(i)]<<'\t'<<map.residual_scale[static_cast<std::size_t>(i)]<<'\n';
        std::ofstream unknowns(prefix.string()+"_unknown_map.tsv");
        unknowns<<"node_id\tcluster_id\tpl\tc\tpg\tsw\trl\trd\trg\trc\tpressure_gauge\n";
        for(std::size_t i=0;i<map.node.size();++i){const auto&d=map.node[i];unknowns<<network_.pores[i].id<<'\t'<<network_.cluster_of_pore[i]+1<<'\t'
          <<d.pl<<'\t'<<d.c<<'\t'<<d.pg<<'\t'<<d.sw<<'\t'<<d.rl<<'\t'<<d.rd<<'\t'<<d.rg<<'\t'<<d.rc<<'\t'<<d.pressure_gauge<<'\n';}
        std::ofstream metadata(prefix.string()+"_metadata.json");metadata<<std::setprecision(17)
          <<"{\n  \"time_s\": "<<time_<<",\n  \"step_number_before\": "<<step_number_
          <<",\n  \"trial_dt_s\": "<<dt<<",\n  \"newton_iteration_zero_based\": "<<iteration
          <<",\n  \"scalar_rows\": "<<scaled.rows()<<",\n  \"scalar_nonzeros\": "<<scaled.nonZeros()
          <<",\n  \"linear_backend\": \"amgx\",\n  \"amgx_mode\": \"dDDI\",\n  \"block_dimension\": 2,\n  \"matrix_market_index_base\": 1,\n  \"runtime_csr_index_base\": 0,\n  \"amgx_config_file\": \""<<config_.linear.amgx_config_file.string()<<"\"\n}\n";
        std::filesystem::copy_file(config_.linear.amgx_config_file,prefix.string()+"_amgx_config.json",
          std::filesystem::copy_options::overwrite_existing);
      }
      linear_system_dumped_=true;
    }
    Eigen::VectorXd delta_scaled;
    bool linear_ok = false;
    if (config_.linear.backend == "eigen_sparselu") {
      Eigen::SparseLU<Eigen::SparseMatrix<double>> linear;
      Eigen::SparseMatrix<double> column_major = scaled;
      const auto setup_started = std::chrono::steady_clock::now();
      linear.compute(column_major);
      report.linear_setup_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - setup_started).count();
      if (linear.info() == Eigen::Success) {
        const auto solve_started = std::chrono::steady_clock::now();
        delta_scaled = linear.solve(rhs);
        report.linear_solve_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - solve_started).count();
        linear_ok = linear.info() == Eigen::Success;
      }
      report.linear_iterations += 1;
    } else if (config_.linear.backend == "eigen_bicgstab_ilut") {
      Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>,
                      Eigen::IncompleteLUT<double>>
          linear;
      linear.setMaxIterations(config_.linear.max_iterations);
      linear.setTolerance(config_.linear.relative_tolerance);
      linear.preconditioner().setDroptol(config_.linear.ilut_drop_tolerance);
      linear.preconditioner().setFillfactor(config_.linear.ilut_fill_factor);
      const auto setup_started = std::chrono::steady_clock::now();
      linear.compute(scaled);
      report.linear_setup_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - setup_started).count();
      if (linear.info() == Eigen::Success) {
        const auto solve_started = std::chrono::steady_clock::now();
        delta_scaled = linear.solve(rhs);
        report.linear_solve_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - solve_started).count();
        linear_ok = linear.info() == Eigen::Success;
      }
      if (linear.iterations() >
          static_cast<Eigen::Index>(std::numeric_limits<int>::max())) {
        report.failure = "linear iteration count overflow";
        return report;
      }
      report.linear_iterations += static_cast<int>(linear.iterations());
    } else if (config_.linear.backend == "pardiso") {
      try {
        PardisoSolveMetrics metrics;
        delta_scaled = pardiso_solver_->solve(scaled, rhs, metrics);
        linear_ok = true;
        report.linear_iterations += 1;
        report.csr_seconds += metrics.csr_seconds;
        report.linear_symbolic_analysis_seconds +=
            metrics.symbolic_analysis_seconds;
        report.linear_numeric_factorization_seconds +=
            metrics.numeric_factorization_seconds;
        report.linear_setup_seconds += metrics.symbolic_analysis_seconds +
                                       metrics.numeric_factorization_seconds;
        report.linear_solve_seconds += metrics.solve_seconds;
        report.pardiso_symbolic_analysis_calls +=
            metrics.symbolic_analysis_calls;
        report.pardiso_numeric_factorization_calls +=
            metrics.numeric_factorization_calls;
        report.pardiso_factor_nonzeros =
            std::max(report.pardiso_factor_nonzeros,
                     metrics.factor_nonzeros);
        report.pardiso_peak_memory_kib =
            std::max(report.pardiso_peak_memory_kib,
                     metrics.peak_memory_kib);
        report.linear_residual = std::max(report.linear_residual,
                                          metrics.cpu_relative_residual);
      } catch (const std::exception& error) {
        report.failure = error.what();
        return report;
      }
    } else if (config_.linear.backend == "amgx") {
      try {
        AmgxSolveMetrics metrics;
        delta_scaled = amgx_solver_->solve(scaled, rhs, metrics);
        linear_ok = true;
        report.linear_iterations += metrics.iterations;
        report.amgx_status = metrics.amgx_status;
        report.amgx_upload_all_calls += metrics.matrix_upload_all_calls;
        report.amgx_replace_coefficients_calls += metrics.matrix_replace_coefficients_calls;
        report.amgx_solver_setup_calls += metrics.solver_setup_calls;
        report.csr_seconds += metrics.csr_seconds;
        report.linear_upload_seconds += metrics.upload_seconds;
        report.linear_setup_seconds += metrics.setup_seconds;
        report.linear_solve_seconds += metrics.solve_seconds;
        report.linear_download_seconds += metrics.download_seconds;
        report.linear_residual = std::max(report.linear_residual,
                                          metrics.cpu_relative_residual);
      } catch (const std::exception& error) {
        report.failure = error.what();
        return report;
      }
    } else {
      report.failure = "unknown linear backend: " + config_.linear.backend;
      return report;
    }
    if (!linear_ok || delta_scaled.size() != map.size ||
        !all_finite(delta_scaled)) {
      report.failure = "linear solve failed";
      return report;
    }
    const double linear_residual =
        (scaled * delta_scaled - rhs).norm() / std::max(rhs.norm(), 1e-300);
    report.linear_residual = std::max(report.linear_residual, linear_residual);
    if (!std::isfinite(linear_residual) ||
        linear_residual >
            std::max(10.0 * config_.linear.relative_tolerance, 1e-8)) {
      std::ostringstream message;message.precision(17);
      message<<"linear residual exceeds tolerance: "<<linear_residual;
      report.failure = message.str();
      return report;
    }

    Eigen::VectorXd delta(map.size);
    for (int i = 0; i < map.size; ++i)
      delta[i] =
          delta_scaled[i] * map.unknown_scale[static_cast<std::size_t>(i)];
    double lambda = 1.0;
    std::string limiting_variable = "none";
    std::size_t limiting_node_id = 0;
    auto restrict_lambda = [&](double candidate, const char* variable,
                               std::size_t node) {
      if (candidate < lambda) {
        lambda = candidate;
        limiting_variable = variable;
        limiting_node_id = network_.pores[node].id;
      }
    };
    const double fraction = config_.nonlinear.fraction_to_boundary;
    for (std::size_t i = 0; i < trial.size(); ++i) {
      const auto& d = map.node[i];
      if (d.pl >= 0 && delta[d.pl] < 0.0)
        restrict_lambda(
            static_cast<double>(fraction * (config_.physics.background_abs_Pa + trial[i].pl) / (-delta[d.pl])), "Pl",
            i);
      if (d.c >= 0 && delta[d.c] < 0.0)
        restrict_lambda(
            static_cast<double>(fraction * trial[i].c / (-delta[d.c])), "C", i);
      if (d.pg >= 0 && delta[d.pg] < 0.0)
        restrict_lambda(
            static_cast<double>(fraction * (config_.physics.background_abs_Pa + trial[i].pg) / (-delta[d.pg])), "Pg",
            i);
      if (d.sw >= 0) {
        if (delta[d.sw] < 0.0)
          restrict_lambda(
              static_cast<double>(fraction *
                                  (trial[i].sw - config_.active_set.sw_min) /
                                  (-delta[d.sw])),
              "Sw_min", i);
        if (delta[d.sw] > 0.0) {
          // Near a disappearing bubble, the event state itself may lie just
          // beyond the ordinary fraction-to-boundary guard.  Allow this
          // Newton candidate to reach the event surface; the post-convergence
          // active-set and conservation gates then retire the bubble.  Away
          // from that surface retain the original strict safety margin.
          const double sg = 1.0 - static_cast<double>(trial[i].sw);
          const bool event_direction = sg <= 1000.0 * config_.active_set.sg_off &&
              trial[i].sw + delta[d.sw] > 1.0 - config_.active_set.sg_off;
          const double guard = event_direction ? 1.0 :
              1.0 - config_.active_set.sg_off;
          const double event_fraction = event_direction ? 1.0 : fraction;
          restrict_lambda(
              static_cast<double>(event_fraction * (guard - trial[i].sw) /
                                  delta[d.sw]),
              "Sw_max", i);
        }
      }
    }
    lambda = std::min(1.0, std::max(0.0, lambda));
    bool line_ok = false;
    std::vector<NodeState> candidate;
    double candidate_inlet_pressure = trial_inlet_pressure;
    double last_candidate_norm = std::numeric_limits<double>::infinity();
    double last_tested_lambda = lambda;
    double first_candidate_norm = std::numeric_limits<double>::infinity();
    double best_candidate_norm = std::numeric_limits<double>::infinity();
    std::array<double, 4> first_candidate_blocks{{0.0, 0.0, 0.0, 0.0}};
    std::array<double, 4> best_candidate_blocks{{0.0, 0.0, 0.0, 0.0}};
    int first_nonfinite_row = -1;
    bool first_candidate = true;
    while (lambda >= config_.nonlinear.line_search_min_lambda) {
      last_tested_lambda = lambda;
      candidate = trial;
      candidate_inlet_pressure = trial_inlet_pressure +
          (config_.boundary.inlet_type == "flow" ? lambda * delta[map.inlet_pressure] : 0.0);
      if (config_.boundary.inlet_type == "flow") inlet_pressure_ = candidate_inlet_pressure;
      for (std::size_t i = 0; i < candidate.size(); ++i) {
        const auto& d = map.node[i];
        if (d.pl >= 0) candidate[i].pl += lambda * delta[d.pl];
        if (d.c >= 0) candidate[i].c += lambda * delta[d.c];
        if (d.pg >= 0) candidate[i].pg += lambda * delta[d.pg];
        if (d.sw >= 0) candidate[i].sw += lambda * delta[d.sw];
      }
      Assembly candidate_residual =
          assemble(candidate, old, dt, map, upwind, transfer_active, false);
      inlet_pressure_ = trial_inlet_pressure;
      const double candidate_norm =
          inf_norm_scaled(candidate_residual.residual, map);
      if (!std::isfinite(candidate_norm) && first_nonfinite_row < 0) {
        for (int row = 0; row < candidate_residual.residual.size(); ++row) {
          if (!std::isfinite(candidate_residual.residual[row])) {
            first_nonfinite_row = row;
            break;
          }
        }
      }
      last_candidate_norm = candidate_norm;
      if (first_candidate) {
        first_candidate_norm = candidate_norm;
        first_candidate_blocks = candidate_residual.residual_norms;
        first_candidate = false;
      }
      if (candidate_norm < best_candidate_norm) {
        best_candidate_norm = candidate_norm;
        best_candidate_blocks = candidate_residual.residual_norms;
      }
      const double armijo_target =
          std::max(config_.nonlinear.residual_tolerance,
                   (1.0 - config_.nonlinear.armijo_c1 * lambda) * norm);
      if (std::isfinite(candidate_norm) && candidate_norm <= armijo_target) {
        line_ok = true;
        break;
      }
      lambda *= config_.nonlinear.line_search_reduction;
    }
    if (!line_ok) {
      std::ostringstream message;
      message << "Armijo line search failed at Newton iteration " << iteration
              << ": residual=" << norm
              << ", last_candidate_residual=" << last_candidate_norm
              << ", last_tested_lambda=" << last_tested_lambda
              << ", limiting_variable=" << limiting_variable
              << ", limiting_node_id=" << limiting_node_id
              << ", first_nonfinite_row=" << first_nonfinite_row
              << ", first_candidate_residual=" << first_candidate_norm
              << ", best_candidate_residual=" << best_candidate_norm
              << ", update_norm=" << delta_scaled.lpNorm<Eigen::Infinity>()
              << ", max_linear_residual=" << report.linear_residual
              << ", blocks=[" << assembled.residual_norms[0] << ','
              << assembled.residual_norms[1] << ','
              << assembled.residual_norms[2] << ','
              << assembled.residual_norms[3] << ']' << ", first_blocks=["
              << first_candidate_blocks[0] << ',' << first_candidate_blocks[1]
              << ',' << first_candidate_blocks[2] << ','
              << first_candidate_blocks[3] << ']' << ", best_blocks=["
              << best_candidate_blocks[0] << ',' << best_candidate_blocks[1]
              << ',' << best_candidate_blocks[2] << ','
              << best_candidate_blocks[3] << ']';
      int worst_row = 0;
      double worst_scaled = 0.0;
      for (int row = 0; row < map.size; ++row) {
        const double value =
            std::abs(assembled.residual[row] /
                     map.residual_scale[static_cast<std::size_t>(row)]);
        if (value > worst_scaled) {
          worst_scaled = value;
          worst_row = row;
        }
      }
      for (std::size_t node = 0; node < map.node.size(); ++node) {
        const auto& d = map.node[node];
        if (d.rc == worst_row) {
          message << ", worst=Rc(node=" << node + 1 << ", Pl=" << trial[node].pl
                  << ", Pg=" << trial[node].pg << ", Sw=" << trial[node].sw
                  << ", C=" << trial[node].c << ", dPl=" << delta[d.pl]
                  << ", dPg=" << delta[d.pg] << ", dSw=" << delta[d.sw] << ')';
          break;
        }
        if (d.rg == worst_row) {
          const auto& pore = network_.pores[node];
          const long double coefficient =
              static_cast<long double>(pore.volume) /
              (static_cast<long double>(config_.physics.compressibility) *
               static_cast<long double>(gas_constant) *
               static_cast<long double>(config_.physics.temperature));
          const long double accumulation =
              coefficient *
              ((1.0L - static_cast<long double>(trial[node].sw)) *
                   (static_cast<long double>(trial[node].pg) -
                    static_cast<long double>(old[node].pg)) -
               (config_.physics.background_abs_Pa + old[node].pg) *
                   (static_cast<long double>(trial[node].sw) -
                    static_cast<long double>(old[node].sw))) /
              static_cast<long double>(dt);
          message << ", worst=Rg(node=" << node + 1 << ", Pl=" << trial[node].pl
                  << ", Pg=" << trial[node].pg << ", Sw=" << trial[node].sw
                  << ", C=" << trial[node].c
                  << ", accumulation=" << static_cast<double>(accumulation)
                  << ", transfer=" << assembled.transfer[node]
                  << ", dPl=" << delta[d.pl] << ", dPg=" << delta[d.pg]
                  << ", dSw=" << delta[d.sw] << ')';
          break;
        }
      }
      report.failure = message.str();
      return report;
    }
    trial = std::move(candidate);
    trial_inlet_pressure = candidate_inlet_pressure;
    last_update = 0.0;
    for (int i = 0; i < delta_scaled.size(); ++i)
      last_update = std::max(last_update, std::abs(lambda * delta_scaled[i]));
    report.update_norm = last_update;
  }

  if (!converged) {
    report.failure = "Newton iteration limit reached";
    return report;
  }
  const Balance equation_balance =
      compute_balance(old, trial, final, dt, 0.0, 0.0);
  if (equation_balance.epsilon_moles_equation >
          config_.acceptance.mass_balance_hard_limit ||
      equation_balance.epsilon_volume_equation >
          config_.acceptance.volume_balance_hard_limit) {
    report.failure = "converged step rejected by conservation gate";
    report.balance = equation_balance;
    return report;
  }

  report.min_free_gas_moles=std::numeric_limits<double>::infinity();
  report.min_Sg=std::numeric_limits<double>::infinity();
  struct CrossingInfo {
    bool present = false;
    bool individually_eligible = false;
    bool batched = false;
    double fraction = 1.0;
    double old_ng = 0.0;
    double new_ng = 0.0;
    double old_sg = 0.0;
    double new_sg = 0.0;
    double individual_fraction = 0.0;
  };
  std::vector<CrossingInfo> crossings(trial.size());
  double total_old_free_gas = 0.0;
  for (std::size_t i = 0; i < old.size(); ++i) {
    if (!old[i].active_gas) continue;
    const double value = static_cast<double>(free_gas_moles(
        old[i].pg, old[i].sw, network_.pores[i].volume, config_.physics));
    if (!std::isfinite(value) || value < 0.0) {
      report.failure = "old state contains invalid free gas moles";
      return report;
    }
    total_old_free_gas += value;
  }
  double individually_eligible_crossing_moles = 0.0;
  for(std::size_t i=0;i<trial.size();++i){
    if(!old[i].active_gas)continue;
    const auto& pore=network_.pores[i];
    const double old_ng=static_cast<double>(free_gas_moles(old[i].pg,old[i].sw,pore.volume,config_.physics));
    const double new_ng=trial[i].active_gas?static_cast<double>(free_gas_moles(trial[i].pg,trial[i].sw,pore.volume,config_.physics)):0.0;
    const double old_sg=static_cast<double>(1.0L-old[i].sw),new_sg=static_cast<double>(1.0L-trial[i].sw);
    if(!std::isfinite(new_ng)||new_ng<0.0){report.failure="predicted negative free gas moles";report.dt_change_reason="gas_disappearance_event";return report;}
    report.min_free_gas_moles=std::min(report.min_free_gas_moles,new_ng);report.min_Sg=std::min(report.min_Sg,new_sg);
    if(config_.time_step_control.enabled&&(new_ng<config_.active_set.ng_off||new_sg<config_.active_set.sg_off)){
      double fraction=1.0;
      if(new_ng<config_.active_set.ng_off&&old_ng>new_ng)
        fraction=std::min(fraction,(old_ng-config_.active_set.ng_off)/(old_ng-new_ng));
      if(new_sg<config_.active_set.sg_off&&old_sg>new_sg)
        fraction=std::min(fraction,(old_sg-config_.active_set.sg_off)/(old_sg-new_sg));
      auto& crossing = crossings[i];
      crossing.present = true;
      crossing.fraction = std::clamp(fraction, 0.0, 1.0);
      crossing.old_ng = old_ng;
      crossing.new_ng = new_ng;
      crossing.old_sg = old_sg;
      crossing.new_sg = new_sg;
      crossing.individual_fraction = total_old_free_gas > 0.0
          ? old_ng / total_old_free_gas : 0.0;
      crossing.individually_eligible =
          config_.active_set.batch_microbubble_retirement_enabled &&
          crossing.individual_fraction <=
              config_.active_set.batch_individual_free_gas_fraction;
      if (crossing.individually_eligible)
        individually_eligible_crossing_moles += old_ng;
    }
  }
  const double eligible_aggregate_fraction = total_old_free_gas > 0.0
      ? individually_eligible_crossing_moles / total_old_free_gas : 0.0;
  const bool aggregate_batch_allowed =
      config_.active_set.batch_microbubble_retirement_enabled &&
      eligible_aggregate_fraction <=
          config_.active_set.batch_total_free_gas_fraction_per_step;
  report.batch_aggregate_free_gas_fraction =
      aggregate_batch_allowed ? eligible_aggregate_fraction : 0.0;
  const CrossingInfo* earliest_localized = nullptr;
  std::size_t earliest_localized_node = 0;
  for (std::size_t i = 0; i < crossings.size(); ++i) {
    auto& crossing = crossings[i];
    if (!crossing.present) continue;
    crossing.batched = crossing.individually_eligible && aggregate_batch_allowed;
    if (!crossing.batched && crossing.fraction < 0.999 &&
        (!earliest_localized || crossing.fraction < earliest_localized->fraction)) {
      earliest_localized = &crossing;
      earliest_localized_node = i;
    }
  }
  if (earliest_localized &&
      dt > config_.time.min_dt * (1.0 + 1e-12)) {
    report.gas_disappearance_event = true;
    report.event_timestep_localized = true;
    report.event_node_id = network_.pores[earliest_localized_node].id;
    // Locate the accepted event state infinitesimally on the inactive side
    // of the threshold. Stopping short approaches the threshold
    // asymptotically and can collapse to dt_min without switching active set.
    report.event_dt = std::clamp(1.0001 * dt * earliest_localized->fraction,
                                 config_.time.min_dt, dt);
    report.event_time = time_ + dt * earliest_localized->fraction;
    report.event_gas_moles_before = earliest_localized->old_ng;
    report.event_gas_moles_after = earliest_localized->new_ng;
    report.failure = "gas disappearance event requires smaller dt";
    report.dt_change_reason = "gas_disappearance_event";
    report.next_dt = report.event_dt;
    return report;
  }
  if(!std::isfinite(report.min_free_gas_moles)){report.min_free_gas_moles=0.0;report.min_Sg=0.0;}

  for (const auto& cluster : network_.zero_distance_clusters) {
    if (cluster.boundary_flag != 0) continue;
    std::vector<std::size_t> retiring;
    double transferred = 0.0;
    double dissolved = 0.0;
    double liquid_after = 0.0;
    for (const auto node : cluster.members) {
      const auto& p = network_.pores[node];
      const auto& s = trial[node];
      dissolved += static_cast<double>(p.volume * s.sw * s.c);
      if (s.active_gas) {
        const double ng = static_cast<double>(
            free_gas_moles(s.pg, s.sw, p.volume, config_.physics));
        if (1.0 - s.sw < config_.active_set.sg_off ||
            ng < config_.active_set.ng_off) {
          const double old_ng = static_cast<double>(free_gas_moles(
              old[node].pg, old[node].sw, p.volume, config_.physics));
          const bool batched = crossings[node].present && crossings[node].batched;
          GasDisappearanceEventRecord event;
          event.node_index = node;
          event.node_id = p.id;
          event.estimated_crossing_time = crossings[node].present
              ? time_ + dt * crossings[node].fraction : time_ + dt;
          event.accepted_state_time = time_ + dt;
          event.crossing_fraction = crossings[node].present
              ? crossings[node].fraction : 1.0;
          event.gas_moles_before = old_ng;
          event.gas_moles_after_before_retirement = ng;
          event.gas_moles_transferred_to_dissolved = ng;
          event.gas_volume_before =
              p.volume * static_cast<double>(1.0L - old[node].sw);
          event.gas_volume_after_before_retirement =
              p.volume * static_cast<double>(1.0L - s.sw);
          event.individual_free_gas_fraction = crossings[node].present
              ? crossings[node].individual_fraction
              : (total_old_free_gas > 0.0 ? old_ng / total_old_free_gas : 0.0);
          event.batch_aggregate_free_gas_fraction = batched
              ? eligible_aggregate_fraction : 0.0;
          event.batched = batched;
          report.gas_disappearance_event_records.push_back(event);
          report.gas_disappearance_event = true;
          ++report.retired_bubble_count;
          if (batched) {
            ++report.batched_retired_bubble_count;
            report.batched_retired_moles += ng;
            report.batched_retired_gas_volume +=
                event.gas_volume_after_before_retirement;
          } else {
            report.event_timestep_localized = true;
          }
          if (report.event_node_id == 0) {
            report.event_node_id = p.id;
            report.event_time = event.estimated_crossing_time;
            report.event_dt = dt;
            report.event_gas_moles_before = old_ng;
            report.event_gas_moles_after = ng;
          }
          retiring.push_back(node);
          transferred += ng;
          report.active_set_liquid_volume_ledger +=
              static_cast<double>(p.volume * (1.0 - s.sw));
          liquid_after += p.volume;
          continue;
        }
      }
      liquid_after += static_cast<double>(p.volume * s.sw);
    }
    if (retiring.empty()) continue;
    report.threshold_moles_transferred += transferred;
    for (const auto node : retiring) {
      trial[node].sw = 1.0;
      trial[node].pg = std::numeric_limits<double>::quiet_NaN();
      trial[node].active_gas = false;
      trial[node].permanently_inactive = true;
    }
    if (!(liquid_after > 0.0)) {
      report.failure = "retired gas cluster has no liquid storage";
      return report;
    }
    const double shared_c = (dissolved + transferred) / liquid_after;
    for (const auto node : cluster.members) trial[node].c = shared_c;
  }
  report.active_set_mole_ledger = 0.0;
  Balance adjusted =
      compute_balance(old, trial, final, dt, report.active_set_mole_ledger,
                      report.active_set_liquid_volume_ledger);
  adjusted.epsilon_moles_equation = equation_balance.epsilon_moles_equation;
  adjusted.epsilon_volume_equation = equation_balance.epsilon_volume_equation;
  report.balance = adjusted;
  if (adjusted.epsilon_moles_adjusted >
          config_.acceptance.mass_balance_hard_limit ||
      adjusted.epsilon_volume_adjusted >
          config_.acceptance.volume_balance_hard_limit) {
    report.failure = "active-set-adjusted conservation gate failed";
    return report;
  }

  state_ = std::move(trial);
  inlet_pressure_ = trial_inlet_pressure;
  if (config_.boundary.species_mode == "copy_adjacent_concentration") {
    for (std::size_t i = 0; i < state_.size(); ++i) {
      if (network_.pores[i].source_type != 2) continue;
      const auto edge_index = network_.adjacency[i].front();
      const auto& edge = network_.edges[edge_index];
      const auto neighbor = edge.a == i ? edge.b : edge.a;
      state_[i].c = state_[neighbor].c;
      state_[i].concentration_derived = true;
    }
  }
  time_ += dt;
  ++step_number_;
  cumulative_component_inflow_ += dt * report.balance.component_inflow;
  cumulative_component_outflow_ += dt * report.balance.component_outflow;
  cumulative_liquid_inflow_ += dt * report.balance.liquid_inflow;
  cumulative_liquid_outflow_ += dt * report.balance.liquid_outflow;
  for (std::size_t i = 0; i < final.transfer.size(); ++i) {
    const long double amount = static_cast<long double>(dt) * final.transfer[i];
    cumulative_transfer_by_node_[i] += amount;
    cumulative_interphase_transfer_ += static_cast<double>(amount);
  }
  cumulative_mole_ledger_ += report.active_set_mole_ledger;
  cumulative_volume_ledger_ += report.active_set_liquid_volume_ledger;
  cumulative_net_liquid_inflow_ += dt * report.balance.net_liquid_inflow;
  cumulative_boundary_component_outflow_ +=
      dt * report.balance.boundary_molar_outflow;
  cumulative_outlet_advective_component_outflow_ +=
      dt * report.balance.outlet_advective_component_outflow;
  report.active_after = active_count(state_);
  report.final_assembly = std::move(final);
  if (report.gas_disappearance_event && report.active_after < report.active_before) {
    const std::size_t retired = report.retired_bubble_count > 0
        ? report.retired_bubble_count : report.active_before - report.active_after;
    gas_disappearance_events_ += retired;
    if (!report.gas_disappearance_event_records.empty()) {
      const auto& last_event = report.gas_disappearance_event_records.back();
      last_gas_disappearance_node_id_ = last_event.node_id;
      last_gas_disappearance_time_ = last_event.estimated_crossing_time;
    } else {
      last_gas_disappearance_node_id_ = report.event_node_id;
      last_gas_disappearance_time_ = report.event_time;
    }
  }
  report.accepted = true;
  commit_inlet_pressure = true;
  return report;
}

StepReport Solver::advance_with_retry(double requested_dt) {
  double dt = std::min(requested_dt, config_.time.max_dt);
  StepReport last;
  std::ostringstream failures;
  std::vector<double> rejected_dts;
  std::vector<std::string> retry_failures;
  std::string last_retry_reason;
  for (int retry = 0; retry <= config_.time.max_retries; ++retry) {
    if (dt < config_.time.min_dt) {
      last.failure = "dt fell below dt_min";
      last.dt_change_reason = "dt_min_failure";
      last.dt = dt;
      last.requested_dt = requested_dt;
      last.timestep_retries = retry;
      last.rejected_dts = rejected_dts;
      last.retry_failures = retry_failures;
      return last;
    }
    last = attempt_step(dt);
    last.requested_dt = requested_dt;
    last.timestep_retries = retry;
    last.rejected_dts = rejected_dts;
    last.retry_failures = retry_failures;
    if (last.accepted){if(!last_retry_reason.empty())last.dt_change_reason=last_retry_reason;return last;}
    rejected_dts.push_back(dt);
    retry_failures.push_back(last.failure);
    if (retry != 0) failures << " | ";
    failures << "dt=" << dt << " [" << last.failure << ']';
    if(last.dt_change_reason=="gas_disappearance_event")last_retry_reason="gas_disappearance_event";
    else if (last.failure.find("linear") != std::string::npos)last_retry_reason="linear_solver_retry";
    else last_retry_reason="nonlinear_retry";
    if(dt<=config_.time.min_dt*(1.0+1e-12)){
      last.failure="dt_min failure after "+last.failure;last.dt_change_reason="dt_min_failure";
      last.rejected_dts=rejected_dts;last.retry_failures=retry_failures;return last;
    }
    const double reduced=dt*(config_.time_step_control.enabled ? config_.time_step_control.shrink_factor
                                                               : config_.time.cut_factor);
    dt=last_retry_reason=="gas_disappearance_event"&&last.next_dt>0.0?
       std::max(config_.time.min_dt,std::min(last.next_dt,reduced)):std::max(config_.time.min_dt,reduced);
  }
  last.failure = "maximum timestep retries exceeded: " + failures.str();
  last.requested_dt = requested_dt;
  last.timestep_retries = static_cast<int>(rejected_dts.size());
  last.rejected_dts = std::move(rejected_dts);
  last.retry_failures = std::move(retry_failures);
  return last;
}

void Solver::restore(std::vector<NodeState> state, double time, int step_number,
                     double cumulative_mole_ledger,
                     double cumulative_volume_ledger,
                     double cumulative_net_liquid_inflow,
                     double cumulative_boundary_component_outflow,
                     double cumulative_outlet_advective_component_outflow,
                     double cumulative_liquid_inflow,
                     double cumulative_liquid_outflow,
                     std::vector<long double> cumulative_transfer_by_node) {
  if (state.size() != network_.pores.size())
    throw std::runtime_error("checkpoint node count mismatch");
  if (!std::isfinite(time) || time < 0.0 || step_number < 0 ||
      !std::isfinite(cumulative_mole_ledger) ||
      !std::isfinite(cumulative_volume_ledger))
    throw std::runtime_error("checkpoint has invalid solver metadata");
  for (std::size_t i = 0; i < state.size(); ++i) {
    if (network_.pores[i].source_type == 2 && !std::isfinite(state[i].c)) {
      const auto edge_index = network_.adjacency[i].front();
      const auto& edge = network_.edges[edge_index];
      const auto neighbor = edge.a == i ? edge.b : edge.a;
      state[i].c = state[neighbor].c;
      state[i].concentration_derived = true;
    }
    const auto& s = state[i];
    if (!std::isfinite(s.pl) || config_.physics.background_abs_Pa + s.pl <= 0.0 || !std::isfinite(s.c) ||
        s.c < -1e-12 * config_.scaling.concentration || !std::isfinite(s.sw) ||
        s.sw <= 0.0 || s.sw > 1.0)
      throw std::runtime_error(
          "checkpoint state violates Pl/C/Sw bounds at node " +
          std::to_string(i + 1));
    if (!internal_node(network_, i) && s.active_gas)
      throw std::runtime_error("checkpoint has active gas on reservoir");
    if (s.active_gas && !state_[i].active_gas)
      throw std::runtime_error(
          "checkpoint attempts forbidden activation of an initially liquid "
          "node");
    if (state_[i].permanently_inactive && !s.permanently_inactive)
      throw std::runtime_error(
          "checkpoint clears an irreversible inactive marker");
    if (s.permanently_inactive && s.active_gas)
      throw std::runtime_error(
          "checkpoint attempts forbidden bubble reactivation");
    if (s.active_gas && (!std::isfinite(s.pg) || config_.physics.background_abs_Pa + s.pg <= 0.0 ||
                         s.sw <= config_.active_set.sw_min ||
                         s.sw >= 1.0 - config_.active_set.sg_off))
      throw std::runtime_error(
          "checkpoint active state violates Pg/Sw bounds at node " +
          std::to_string(i + 1));
    if ((!internal_node(network_, i) || !s.active_gas) &&
        std::abs(s.sw - 1.0) > 1e-14)
      throw std::runtime_error(
          "checkpoint inactive/reservoir node must have Sw=1 at node " +
          std::to_string(i + 1));
  }
  for (const auto& cluster : network_.zero_distance_clusters) {
    const auto representative = cluster.members.front();
    for (const auto node : cluster.members) {
      const long double pl_scale =
          std::max(std::abs(state[representative].pl),
                   static_cast<long double>(config_.scaling.pressure));
      const long double c_scale =
          std::max(std::abs(state[representative].c),
                   static_cast<long double>(config_.scaling.concentration));
      if (std::abs(state[node].pl - state[representative].pl) >
              1e-13 * pl_scale ||
          std::abs(state[node].c - state[representative].c) > 1e-13 * c_scale)
        throw std::runtime_error(
            "checkpoint violates shared Pl/C in zero-distance cluster " +
            std::to_string(cluster.id + 1));
    }
  }
  state_ = std::move(state);
  time_ = time;
  step_number_ = step_number;
  cumulative_mole_ledger_ = cumulative_mole_ledger;
  cumulative_volume_ledger_ = cumulative_volume_ledger;
  cumulative_net_liquid_inflow_ = cumulative_net_liquid_inflow;
  cumulative_boundary_component_outflow_ =
      cumulative_boundary_component_outflow;
  cumulative_outlet_advective_component_outflow_ =
      cumulative_outlet_advective_component_outflow;
  cumulative_liquid_inflow_ = cumulative_liquid_inflow;
  cumulative_liquid_outflow_ = cumulative_liquid_outflow;
  if (cumulative_transfer_by_node.empty())
    cumulative_transfer_by_node.assign(state_.size(), 0.0L);
  if (cumulative_transfer_by_node.size() != state_.size())
    throw std::runtime_error("checkpoint cumulative transfer node count mismatch");
  cumulative_transfer_by_node_ = std::move(cumulative_transfer_by_node);
  cumulative_interphase_transfer_ = 0.0;
  for (const auto value : cumulative_transfer_by_node_)
    cumulative_interphase_transfer_ += static_cast<double>(value);
}

}  // namespace bubble
