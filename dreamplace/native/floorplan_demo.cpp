// BSD-3-Clause. See LICENSE. Standalone pattern experiment, not a placement
// dependency.
#include "dreamplace/native/floorplan.hpp"
#include <iomanip>
#include <iostream>
using namespace DreamPlace::native;
int main() {
  auto c = [](double x, double y) {
    return FloorplanCircle{"", x, y, 8, 64, 0, -1};
  };
  std::vector<std::pair<std::string, std::vector<FloorplanCircle>>> patterns;
  std::vector<FloorplanCircle> mesh;
  for (int y = 0; y < 3; ++y)
    for (int x = 0; x < 3; ++x)
      mesh.push_back(c(30 + x * 25 + (y - 1) * .7, 30 + y * 25 + (x - 1) * .6));
  patterns.push_back({"mesh", mesh});
  patterns.push_back({"triangle",
                      {c(60.7, 25), c(47, 50.4), c(72.5, 49.7), c(34.9, 75.5),
                       c(60.3, 74.6), c(85.5, 75.2)}});
  patterns.push_back({"tree",
                      {c(70.5, 20), c(42.7, 50.5), c(98, 49.6), c(33.4, 80.1),
                       c(58.9, 79.5), c(83.7, 80.4), c(109.1, 80)}});
  std::cout << std::setprecision(15) << "[";
  bool first = true;
  for (const auto &[name, v] : patterns) {
    auto r = snap_floorplan(v, 0, 0, 145, 120);
    if (!first)
      std::cout << ",";
    first = false;
    std::cout << "{\"name\":\"" << name
              << "\",\"feasible\":" << (r.feasible ? "true" : "false")
              << ",\"seconds\":" << r.seconds << ",\"rejected\":" << r.rejected
              << ",\"maximum_overlap\":" << r.maximum_overlap
              << ",\"maximum_displacement\":" << r.maximum_displacement;
    for (auto key : {"before", "after"}) {
      const auto &points = std::string(key) == "before" ? v : r.circles;
      std::cout << ",\"" << key << "\":[";
      for (size_t i = 0; i < points.size(); ++i) {
        if (i)
          std::cout << ",";
        std::cout << "[" << points[i].x << "," << points[i].y << ","
                  << points[i].radius << "]";
      }
      std::cout << "]";
    }
    std::cout << ",\"relations\":[";
    for (size_t i = 0; i < r.relations.size(); ++i) {
      if (i)
        std::cout << ",";
      const auto &q = r.relations[i];
      std::cout << "{\"kind\":\"" << q.kind << "\",\"axis\":" << q.axis
                << ",\"value\":" << q.value << ",\"nodes\":[";
      for (size_t j = 0; j < q.nodes.size(); ++j) {
        if (j)
          std::cout << ",";
        std::cout << q.nodes[j];
      }
      std::cout << "]}";
    }
    std::cout << "]}";
  }
  std::cout << "]\n";
}
