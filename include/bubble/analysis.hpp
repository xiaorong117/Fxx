#pragma once

#include <filesystem>
#include <fstream>
#include <vector>

#include "bubble/model.hpp"

namespace bubble {

struct ScientificSnapshot {
  double main_real_pore_volume_m3 = 0.0;
  double main_real_gas_saturation = 0.0;
  double sample_pore_volume_m3 = 0.0;
  double mean_sw = 0.0, mean_sg = 0.0;
  double gas_volume_m3 = 0.0, free_gas_mol = 0.0;
  double dissolved_gas_mol = 0.0;
  std::size_t active_bubbles = 0;
  double injected_pore_volumes = 0.0;
};

class ScientificRecorder {
 public:
  ScientificRecorder(const Solver& solver, bool restart);
  ScientificSnapshot snapshot(const Solver& solver) const;
  void append(const Solver& solver, const StepReport& report);
  double initial_main_gas_saturation() const { return initial_main_sg_; }
  double initial_gas_volume() const { return initial_gas_volume_; }
  double initial_free_moles() const { return initial_free_moles_; }
  const std::vector<int>& node_segments() const { return node_segment_; }
  double normalized_x(std::size_t node) const;
  double adjacent_constriction_p50(std::size_t node) const;

 private:
  struct SegmentStatic {
    std::vector<std::size_t> nodes;
    double volume = 0.0;
  };
  const Network* network_ = nullptr;
  const Config* config_ = nullptr;
  std::filesystem::path output_;
  int segments_ = 20;
  double xmin_ = 0.0, xmax_ = 1.0;
  double sample_volume_ = 0.0, main_volume_ = 0.0;
  double initial_liquid_ = 0.0, initial_total_ = 0.0;
  double initial_gas_volume_ = 0.0, initial_free_moles_ = 0.0;
  double initial_main_sg_ = 0.0;
  std::vector<int> node_segment_;
  std::vector<double> adjacent_constriction_p50_;
  std::vector<SegmentStatic> segment_static_;
  std::ofstream global_, segment_, balance_, events_;

  bool sample_node(std::size_t node) const;
  void write_geometry() const;
};

}  // namespace bubble
