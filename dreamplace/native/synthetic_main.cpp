// BSD-3-Clause. See LICENSE.
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include "dreamplace/native/synthetic.hpp"

namespace syn = DreamPlace::native::synthetic;
int main(int argc, char** argv) {
  try {
    std::string which = "all", saved;
    bool animate = false;
    bool utilization_set = false;
    std::filesystem::path out = "results/synthetic";
    syn::Spec mesh, torus, tree;
    mesh.kind = syn::Kind::mesh;
    torus.kind = syn::Kind::torus;
    torus.width = 12;
    torus.height = 8;
    tree.kind = syn::Kind::tree;
    DreamPlace::native::Options options;
    options.max_iterations = 1500;
    options.bins_x = options.bins_y = 16;
    auto integer = [](const std::string& text) {
      size_t end = 0;
      int n = std::stoi(text, &end);
      if (end != text.size()) throw std::invalid_argument("expected integer, got " + text);
      return n;
    };
    auto real = [](const std::string& text) {
      size_t end = 0;
      double n = std::stod(text, &end);
      if (end != text.size()) throw std::invalid_argument("expected number, got " + text);
      return n;
    };
    for (int i = 1; i < argc; ++i) {
      std::string arg = argv[i];
      if (arg == "--help") {
        std::cout
            << "synthetic_place [--case all|mesh|mesh-grid|mesh-grid-io|torus|tree] [--out DIR] [--load SAVED.hhds]\n"
               "  [--mesh-width N] [--mesh-height N] [--torus-width N] [--torus-height N]\n"
               "  [--tree-depth N] [--seed N] [--iterations N] [--bins N] [--threads N]\n"
               "  [--utilization FRACTION] [--animate] [--frame-step N]\n"
               "  [--stop-overflow F] [--target-density F] [--density-weight F] [--max-step F]\n";
        return 0;
      }
      if (arg == "--animate") {
        animate = true;
        continue;
      }
      if (i + 1 == argc) throw std::invalid_argument("missing value for " + arg);
      std::string value = argv[++i];
      if (arg == "--case")
        which = value;
      else if (arg == "--out")
        out = value;
      else if (arg == "--load")
        saved = value;
      else if (arg == "--mesh-width")
        mesh.width = integer(value);
      else if (arg == "--mesh-height")
        mesh.height = integer(value);
      else if (arg == "--torus-width")
        torus.width = integer(value);
      else if (arg == "--torus-height")
        torus.height = integer(value);
      else if (arg == "--tree-depth")
        tree.tree_depth = integer(value);
      else if (arg == "--iterations")
        options.max_iterations = integer(value);
      else if (arg == "--bins")
        options.bins_x = options.bins_y = integer(value);
      else if (arg == "--threads")
        options.num_threads = integer(value);
      else if (arg == "--frame-step")
        options.observer_interval = integer(value);
      else if (arg == "--stop-overflow")
        options.stop_overflow = real(value);
      else if (arg == "--target-density")
        options.target_density = real(value);
      else if (arg == "--density-weight")
        options.initial_density_weight = real(value);
      else if (arg == "--max-step")
        options.max_step = real(value);
      else if (arg == "--utilization") {
        mesh.utilization = torus.utilization = tree.utilization = real(value);
        utilization_set = true;
      } else if (arg == "--seed") {
        int seed = integer(value);
        if (seed < 0) throw std::invalid_argument("seed must be nonnegative");
        mesh.seed = torus.seed = tree.seed = seed;
      } else
        throw std::invalid_argument("unknown option " + arg);
    }
    if (out.is_relative())
      if (const char* workspace = std::getenv("BUILD_WORKSPACE_DIRECTORY"))
        out = std::filesystem::path(workspace) / out;
    if (which != "all" && which != "mesh" && which != "mesh-grid" && which != "mesh-grid-io" && which != "torus" &&
        which != "tree")
      throw std::invalid_argument("case must be all, mesh, mesh-grid, mesh-grid-io, torus or tree");
    auto run = [&](syn::Graph g) {
      auto directory = std::filesystem::absolute(out / g.name);
      auto name = g.name;
      auto result = syn::run(std::move(g), directory, options, animate);
      std::cout << name << ": HPWL " << result.initial_hpwl << " -> " << result.hpwl << " um, legal=" << result.legal
                << ", GP overflow=" << result.overflow << ", converged=" << result.converged << "\n  "
                << directory.string() << "/placement.png\n";
      if (animate) std::cout << "  " << directory.string() << "/animation.html\n";
    };
    if (!saved.empty()) {
      if (utilization_set)
        throw std::invalid_argument("--utilization sizes generated fixtures; --load keeps the saved floorplan");
      std::filesystem::path source = saved;
      if (source.is_relative())
        if (const char* workspace = std::getenv("BUILD_WORKSPACE_DIRECTORY"))
          source = std::filesystem::path(workspace) / source;
      run(syn::load(source));
    } else {
      if (which == "all" || which == "mesh") run(syn::generate(mesh));
      if (which == "all" || which == "mesh-grid") {
        mesh.initialization = syn::Initialization::mesh_grid;
        run(syn::generate(mesh));
      }
      if (which == "all" || which == "mesh-grid-io") {
        mesh.initialization = syn::Initialization::mesh_grid;
        mesh.side_io = true;
        run(syn::generate(mesh));
      }
      if (which == "all" || which == "torus") run(syn::generate(torus));
      if (which == "all" || which == "tree") run(syn::generate(tree));
    }
  } catch (const std::exception& e) {
    std::cerr << "synthetic_place: " << e.what() << '\n';
    return 1;
  }
}
