#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

#include "bubble/model.hpp"

namespace bubble {
namespace {

std::vector<std::vector<double>> load_numeric(const std::filesystem::path& path,
                                              std::size_t columns) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open geometry file: " + path.string());
  std::vector<std::vector<double>> rows;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#') continue;
    std::istringstream stream(line);
    std::vector<double> row;
    double value = 0.0;
    while (stream >> value) row.push_back(value);
    if (!stream.eof() || row.size() != columns) {
      throw std::runtime_error(
          path.string() + ":" + std::to_string(line_number) +
          ": expected exactly " + std::to_string(columns) +
          " numeric columns, got " + std::to_string(row.size()));
    }
    rows.push_back(std::move(row));
  }
  return rows;
}

std::size_t exact_id(double value, const std::string& field) {
  if (!std::isfinite(value) || value < 1.0 || std::floor(value) != value ||
      value > static_cast<double>(std::numeric_limits<std::size_t>::max())) {
    throw std::runtime_error(field + " must be a positive integer");
  }
  return static_cast<std::size_t>(value);
}

std::size_t exact_nonnegative_id(double value, const std::string& field) {
  if (!std::isfinite(value) || value < 0.0 || std::floor(value) != value ||
      value > static_cast<double>(std::numeric_limits<std::size_t>::max())) {
    throw std::runtime_error(field + " must be a nonnegative integer");
  }
  return static_cast<std::size_t>(value);
}

int exact_int(double value, const std::string& field) {
  if (!std::isfinite(value) || std::floor(value) != value ||
      value < static_cast<double>(std::numeric_limits<int>::min()) ||
      value > static_cast<double>(std::numeric_limits<int>::max())) {
    throw std::runtime_error(field + " must be an integer");
  }
  return static_cast<int>(value);
}

class UnionFind {
 public:
  explicit UnionFind(std::size_t size) : parent_(size), rank_(size, 0) {
    std::iota(parent_.begin(), parent_.end(), std::size_t{0});
  }
  std::size_t find(std::size_t value) {
    if (parent_[value] != value) parent_[value] = find(parent_[value]);
    return parent_[value];
  }
  void unite(std::size_t a, std::size_t b) {
    a = find(a);
    b = find(b);
    if (a == b) return;
    if (rank_[a] < rank_[b]) std::swap(a, b);
    parent_[b] = a;
    if (rank_[a] == rank_[b]) ++rank_[a];
  }

 private:
  std::vector<std::size_t> parent_;
  std::vector<unsigned int> rank_;
};

void finite_positive(double value, const std::string& field) {
  if (!std::isfinite(value) || value <= 0.0) {
    throw std::runtime_error(field + " must be finite and positive");
  }
}

void audit_connect(const std::filesystem::path& path,
                   const std::vector<std::vector<std::size_t>>& adjacency) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("cannot open connect audit file: " +
                             path.string());
  std::string line;
  std::size_t expected_id = 1;
  while (std::getline(input, line)) {
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#') continue;
    std::istringstream stream(line);
    std::size_t id = 0, degree = 0;
    if (!(stream >> id >> degree) || id != expected_id ||
        id > adjacency.size() || degree != adjacency[id - 1].size()) {
      throw std::runtime_error("connect audit degree/ID mismatch at row " +
                               std::to_string(expected_id));
    }
    ++expected_id;
  }
  if (expected_id != adjacency.size() + 1) {
    throw std::runtime_error("connect audit row count mismatch");
  }
}

}  // namespace

Network Network::load(const Config& config) {
  Network n;
  const auto pore_rows = load_numeric(config.io.pore_file, 12);
  const auto edge_rows = load_numeric(config.io.throat_file, 11);
  n.pores.reserve(pore_rows.size());
  for (std::size_t k = 0; k < pore_rows.size(); ++k) {
    const auto& r = pore_rows[k];
    Pore p;
    p.id = exact_id(r[0], "pore id");
    if (p.id != k + 1)
      throw std::runtime_error("pore IDs must be continuous from 1");
    p.source_type = exact_int(r[1], "source_type");
    p.source_id = exact_id(r[2], "source_id");
    p.x = r[3] * config.geometry.length_to_m;
    p.y = r[4] * config.geometry.length_to_m;
    p.z = r[5] * config.geometry.length_to_m;
    p.radius = r[6] * config.geometry.length_to_m;
    p.area = r[7] * config.geometry.area_to_m2;
    p.shape_factor = r[8];
    p.volume = r[9] * config.geometry.volume_to_m3;
    p.boundary_flag = exact_int(r[10], "boundary_flag");
    p.input_sw = r[11];
    n.pores.push_back(p);
  }
  n.adjacency.resize(n.pores.size());
  n.edges.reserve(edge_rows.size());
  for (std::size_t k = 0; k < edge_rows.size(); ++k) {
    const auto& r = edge_rows[k];
    Edge e;
    e.id = exact_id(r[0], "edge id");
    if (e.id != k + 1)
      throw std::runtime_error("edge IDs must be continuous from 1");
    const auto a = exact_id(r[1], "edge endpoint a");
    const auto b = exact_id(r[2], "edge endpoint b");
    if (a > n.pores.size() || b > n.pores.size() || a == b)
      throw std::runtime_error("invalid edge endpoint at edge " +
                               std::to_string(e.id));
    e.a = a - 1;
    e.b = b - 1;
    e.source_old_throat_id = exact_nonnegative_id(r[3], "source_old_throat_id");
    e.side = exact_int(r[4], "side");
    e.length = r[5] * config.geometry.length_to_m;
    e.input_sw = r[6];
    e.radius = r[7] * config.geometry.length_to_m;
    e.area = r[8] * config.geometry.area_to_m2;
    e.shape_factor = r[9];
    e.source_volume = r[10] * config.geometry.volume_to_m3;
    n.edges.push_back(e);
    n.adjacency[e.a].push_back(k);
    n.adjacency[e.b].push_back(k);
  }
  n.validate(config);
  n.build_degenerate_geometry(config);
  std::vector<bool> visited(n.pores.size(), false);
  n.component_id.assign(n.pores.size(), 0);
  n.boundary_connection_class.assign(n.pores.size(), 0);
  std::size_t largest_main_real_nodes = 0;
  for (std::size_t root = 0; root < n.pores.size(); ++root) {
    if (n.adjacency[root].empty()) ++n.isolated_pores;
    if (visited[root]) continue;
    ++n.connected_components;
    bool has_inlet = false, has_outlet = false;
    std::vector<std::size_t> stack{root};
    std::vector<std::size_t> members;
    visited[root] = true;
    while (!stack.empty()) {
      const auto node = stack.back();
      stack.pop_back();
      members.push_back(node);
      has_inlet = has_inlet ||
                  n.pores[node].boundary_flag == config.boundary.inlet_flag;
      has_outlet = has_outlet ||
                   n.pores[node].boundary_flag == config.boundary.outlet_flag;
      for (const auto edge_index : n.adjacency[node]) {
        const auto& edge = n.edges[edge_index];
        const auto neighbor = edge.a == node ? edge.b : edge.a;
        if (!visited[neighbor]) {
          visited[neighbor] = true;
          stack.push_back(neighbor);
        }
      }
    }
    if (!has_inlet && !has_outlet) ++n.boundaryless_components;
    std::size_t real_nodes = 0, virtual_nodes = 0;
    for (const auto node : members) {
      real_nodes += static_cast<std::size_t>(n.pores[node].source_type != 2);
      virtual_nodes += static_cast<std::size_t>(n.pores[node].source_type == 2);
    }
    int connection_class =
        has_inlet && has_outlet ? 3 : (has_inlet ? 1 : (has_outlet ? 2 : 0));
    if (members.size() == 2 && real_nodes == 1 && virtual_nodes == 1)
      connection_class = 4;
    for (const auto node : members) {
      n.component_id[node] = n.connected_components;
      n.boundary_connection_class[node] = connection_class;
    }
    if (connection_class == 3 && real_nodes > largest_main_real_nodes) {
      largest_main_real_nodes = real_nodes;
      n.main_flow_component_id = n.connected_components;
    }
  }
  audit_connect(config.io.connect_audit_file, n.adjacency);
  return n;
}

void Network::validate(const Config& config) const {
  if (pores.empty()) throw std::runtime_error("network has no pores");
  for (std::size_t k = 0; k < pores.size(); ++k) {
    const auto& p = pores[k];
    if (p.id != k + 1) throw std::runtime_error("non-continuous pore IDs");
    if (p.source_type != 0 && p.source_type != 1 && p.source_type != 2)
      throw std::runtime_error("source_type must be 0, 1, or 2");
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
      throw std::runtime_error("non-finite pore coordinate");
    const std::string label = "pore " + std::to_string(p.id) + " ";
    finite_positive(p.radius, label + "radius");
    finite_positive(p.area, label + "area");
    finite_positive(p.shape_factor, label + "shape_factor");
    if (!std::isfinite(p.input_sw) || p.input_sw < 0.0 || p.input_sw > 1.0)
      throw std::runtime_error("pore Sw must lie in [0,1]");
    if (!std::isfinite(p.volume) || p.volume < 0.0)
      throw std::runtime_error(label + "volume must be finite and nonnegative");
    if (p.volume == 0.0 && p.input_sw < 1.0)
      throw std::runtime_error(label + "has zero volume with Sw<1");
    if (p.boundary_flag != 0 && p.boundary_flag != config.boundary.inlet_flag &&
        p.boundary_flag != config.boundary.outlet_flag)
      throw std::runtime_error("unsupported boundary flag " +
                               std::to_string(p.boundary_flag));
    if (config.virtual_boundary.enabled) {
      if (p.source_type == 2 && p.boundary_flag == 0)
        throw std::runtime_error(
            "virtual boundary node must carry a boundary flag");
      if (p.source_type != 2 && p.boundary_flag != 0)
        throw std::runtime_error(
            "real control volume must not carry a boundary flag");
      if (p.source_type == 2 && p.input_sw != 1.0)
        throw std::runtime_error("virtual boundary node must have Sw=1");
      if (p.source_type == 2 && adjacency[k].size() != 1)
        throw std::runtime_error(
            "virtual boundary node coordination must equal one");
    }
  }
  if (adjacency.size() != pores.size())
    throw std::runtime_error("adjacency size mismatch");
  for (std::size_t k = 0; k < edges.size(); ++k) {
    const auto& e = edges[k];
    if (e.id != k + 1 || e.a >= pores.size() || e.b >= pores.size() ||
        e.a == e.b)
      throw std::runtime_error("invalid edge ID/endpoints");
    if (e.side != 1 && e.side != 2)
      throw std::runtime_error("edge side must be 1 or 2");
    const std::string label = "edge " + std::to_string(e.id) + " ";
    if (!std::isfinite(e.length) || e.length < 0.0)
      throw std::runtime_error(label +
                               "half-length must be finite and nonnegative");
    finite_positive(e.radius, label + "radius");
    finite_positive(e.area, label + "area");
    finite_positive(e.shape_factor, label + "shape_factor");
    if (!std::isfinite(e.input_sw) || e.input_sw < 0.0 || e.input_sw > 1.0)
      throw std::runtime_error("edge Sw must lie in [0,1]");
    if (!std::isfinite(e.source_volume) || e.source_volume < 0.0)
      throw std::runtime_error(
          label + "source volume metadata must be finite and nonnegative");
    if (e.source_volume == 0.0 && e.input_sw < 1.0)
      throw std::runtime_error(label +
                               "has zero source volume metadata with Sw<1");
    const bool virtual_a = pores[e.a].source_type == 2;
    const bool virtual_b = pores[e.b].source_type == 2;
    if (virtual_a || virtual_b) {
      if (!config.virtual_boundary.enabled || virtual_a == virtual_b)
        throw std::runtime_error("invalid virtual boundary edge endpoints");
      if (!(e.length > 0.0) || e.source_old_throat_id != 0)
        throw std::runtime_error(
            "virtual boundary edge must be positive and have no source throat");
      const auto virtual_node = virtual_a ? e.a : e.b;
      const auto real_node = virtual_a ? e.b : e.a;
      if (pores[virtual_node].source_id != pores[real_node].id)
        throw std::runtime_error(
            "virtual boundary source ID does not match its real pore");
    } else if (e.source_old_throat_id == 0) {
      throw std::runtime_error("non-boundary edge lacks a source throat ID");
    }
  }
}

void Network::build_degenerate_geometry(const Config& config) {
  cluster_of_pore.clear();
  zero_distance_clusters.clear();
  zero_distance_edges.clear();
  massless_junctions.clear();
  production_positive_edges = 0;
  throats_with_two_zero_halves = 0;
  ignored_boundary_edges = 0;
  virtual_boundary_nodes = 0;
  boundary_edges = 0;

  UnionFind united(pores.size());
  std::map<std::size_t, std::size_t> zero_halves_per_throat;
  for (std::size_t k = 0; k < edges.size(); ++k) {
    const auto& edge = edges[k];
    if (edge.length == 0.0) {
      united.unite(edge.a, edge.b);
      zero_distance_edges.push_back(k);
      ++zero_halves_per_throat[edge.source_old_throat_id];
    }
  }
  for (const auto& item : zero_halves_per_throat) {
    if (item.second == 2) ++throats_with_two_zero_halves;
    if (item.second > 2)
      throw std::runtime_error(
          "source throat has more than two zero half-edges");
  }

  cluster_of_pore.resize(pores.size());
  std::map<std::size_t, std::size_t> root_to_cluster;
  for (std::size_t node = 0; node < pores.size(); ++node) {
    const auto root = united.find(node);
    auto position = root_to_cluster.find(root);
    if (position == root_to_cluster.end()) {
      const std::size_t id = zero_distance_clusters.size();
      position = root_to_cluster.emplace(root, id).first;
      ZeroDistanceCluster cluster;
      cluster.id = id;
      zero_distance_clusters.push_back(std::move(cluster));
    }
    const auto cluster = position->second;
    cluster_of_pore[node] = cluster;
    zero_distance_clusters[cluster].members.push_back(node);
    if (pores[node].volume == 0.0 && pores[node].source_type != 2)
      massless_junctions.push_back(node);
    if (pores[node].source_type == 2) ++virtual_boundary_nodes;
  }

  for (std::size_t k = 0; k < edges.size(); ++k) {
    const auto& edge = edges[k];
    const auto ca = cluster_of_pore[edge.a];
    const auto cb = cluster_of_pore[edge.b];
    if (edge.length == 0.0) {
      zero_distance_clusters[ca].zero_distance_edges.push_back(k);
    } else if (ca != cb) {
      zero_distance_clusters[ca].external_positive_edges.push_back(k);
      zero_distance_clusters[cb].external_positive_edges.push_back(k);
    }
  }

  for (auto& cluster : zero_distance_clusters) {
    bool has_inlet = false;
    bool has_outlet = false;
    int first_boundary = 0;
    for (const auto node : cluster.members) {
      const int flag = pores[node].boundary_flag;
      if (flag == 0) continue;
      if (first_boundary == 0) first_boundary = flag;
      has_inlet = has_inlet || flag == config.boundary.inlet_flag;
      has_outlet = has_outlet || flag == config.boundary.outlet_flag;
      if (flag != config.boundary.inlet_flag &&
          flag != config.boundary.outlet_flag)
        throw std::runtime_error(
            "zero-distance cluster contains unsupported boundary flag");
    }
    if (has_inlet && has_outlet &&
        config.boundary.pl_in != config.boundary.pl_out) {
      throw std::runtime_error(
          "zero-distance cluster contains conflicting pressure boundaries");
    }
    // With equal boundary pressures, inlet semantics are deterministic and
    // provide the shared concentration value for the degenerate cluster.
    cluster.boundary_flag =
        has_inlet ? config.boundary.inlet_flag : first_boundary;
  }
  for (const auto& edge : edges) {
    if (edge.length == 0.0) continue;
    const auto ca = cluster_of_pore[edge.a];
    const auto cb = cluster_of_pore[edge.b];
    if (ca == cb) continue;
    const bool boundary_a = zero_distance_clusters[ca].boundary_flag != 0;
    const bool boundary_b = zero_distance_clusters[cb].boundary_flag != 0;
    if (pores[edge.a].source_type == 2 || pores[edge.b].source_type == 2)
      ++boundary_edges;
    if (boundary_a && boundary_b) {
      ++ignored_boundary_edges;
    } else {
      ++production_positive_edges;
    }
  }
}

}  // namespace bubble
