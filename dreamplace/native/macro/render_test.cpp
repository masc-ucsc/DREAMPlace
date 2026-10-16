// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/render.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace DreamPlace::native;
namespace m = DreamPlace::native::macro;
int main() {
  try {
    Design d;
    d.layout = {0, 0, 40, 20, 1, 1};
    d.num_movable = 2;
    d.node_x = {1, 25};
    d.node_y = {1, 11};
    d.node_size_x = {8, 8};
    d.node_size_y = {8, 8};
    d.build_maps();
    std::vector<m::NodeInfo> info = {{"a.m", "a", "ram", true},
                                     {"b.m", "b", "ram", true}};
    m::Topology t;
    for (int i = 0; i < 2; ++i) {
      m::Group g;
      g.id = info[i].parent;
      g.x = d.node_x[i];
      g.y = d.node_y[i];
      g.width = g.height = 8;
      g.members = {{info[i].name}};
      t.groups.push_back(g);
    }
    t.order = {{"r", "a", "b", Axis::x, {false, 1}}};
    t.anchors = {{"a", "a", Side::left, 2, {false, 1}}};
    const char *dir = std::getenv("TEST_UNDECLARED_OUTPUTS_DIR");
    if (!dir)
      dir = std::getenv("TEST_TMPDIR");
    auto path = std::filesystem::path(dir ? dir : "/tmp") / "h2d.png";
    m::write_png(path, d, info, t, 160);
    std::ifstream input(path, std::ios::binary);
    std::string bytes{std::istreambuf_iterator<char>(input), {}};
    if (bytes.size() < 100 || bytes.substr(1, 3) != "PNG" ||
        bytes.substr(bytes.size() - 8, 4) != "IEND")
      throw std::runtime_error("invalid H2D PNG");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
