// BSD-3-Clause. See LICENSE.
#include "thirdparty/flute/flute.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace {

// Deterministic nets spanning the LUT degrees (<= 9) and the merge path (> 9).
std::string route_signature() {
  std::string out;
  unsigned state = 12345;
  const auto next = [&state] {
    state = state * 1103515245u + 12345u;
    return int((state >> 8) % 1000);
  };
  for (int net = 0; net < 400; ++net) {
    const int degree = 2 + net % 24;
    std::vector<flute::DType> x(degree), y(degree);
    for (int i = 0; i < degree; ++i) {
      x[i] = next();
      y[i] = next();
    }
    const auto wl = flute::flute_wl(degree, x.data(), y.data(), 3);
    auto tree = flute::flute(degree, x.data(), y.data(), ACCURACY);
    out += std::to_string(wl) + ':' + std::to_string(tree.length);
    for (int b = 0; b < 2 * tree.deg - 2; ++b) {
      const auto &branch = tree.branch[b];
      out += ' ' + std::to_string(branch.x) + ',' + std::to_string(branch.y) +
             ',' + std::to_string(branch.n);
    }
    out += '\n';
    std::free(tree.branch);
  }
  return out;
}

std::string slurp(const std::string &path) {
  std::ifstream input(path, std::ios::binary);
  return std::string((std::istreambuf_iterator<char>(input)),
                     std::istreambuf_iterator<char>());
}

} // namespace

int main() {
  try {
    const std::string powv = "thirdparty/flute/lut.ICCAD2015/POWV9.dat";
    const std::string post = "thirdparty/flute/lut.ICCAD2015/POST9.dat";
    const auto rejects = [](const std::string &a, const std::string &b) {
      try {
        flute::readLUT(a.c_str(), b.c_str());
      } catch (const std::runtime_error &) {
        return;
      }
      throw std::runtime_error("invalid FLUTE table was accepted");
    };
    rejects("/nonexistent/POWV9.dat", post);
    rejects(powv, "/nonexistent/POST9.dat");
    const auto scratch = std::filesystem::path(std::getenv("TEST_TMPDIR"));
    const auto empty = (scratch / "empty.dat").string();
    std::ofstream(empty).close();
    rejects(empty, post);
    rejects(powv, empty);
    const auto truncated = [&](const std::string &source, const char *name) {
      std::ifstream input(source, std::ios::binary);
      std::string prefix(64, '\0');
      input.read(prefix.data(), prefix.size());
      const auto path = (scratch / name).string();
      std::ofstream output(path, std::ios::binary);
      output.write(prefix.data(), input.gcount());
      output.close();
      return path;
    };
    rejects(truncated(powv, "powv.dat"), post);
    rejects(powv, truncated(post, "post.dat"));
    // The buffer entry point applies the same validation.
    const std::string powv_data = slurp(powv), post_data = slurp(post);
    const auto rejects_buffer = [](const std::string &a, const std::string &b) {
      try {
        flute::readLUT(a.data(), a.size(), b.data(), b.size());
      } catch (const std::runtime_error &) {
        return;
      }
      throw std::runtime_error("invalid FLUTE buffer was accepted");
    };
    rejects_buffer(std::string(), post_data);
    rejects_buffer(powv_data, std::string());
    rejects_buffer(powv_data.substr(0, powv_data.size() / 2), post_data);
    rejects_buffer(powv_data, post_data.substr(0, post_data.size() / 2));
    for (const size_t size : {size_t{0}, size_t{64}}) {
      for (bool null_powv : {false, true}) {
        bool rejected = false;
        try {
          flute::readLUT(null_powv ? nullptr : powv_data.data(),
                         null_powv ? size : powv_data.size(),
                         null_powv ? post_data.data() : nullptr,
                         null_powv ? post_data.size() : size);
        } catch (const std::runtime_error &) {
          rejected = true;
        }
        if (!rejected)
          throw std::runtime_error("null FLUTE buffer was accepted");
      }
    }

    // Tables load once per process, so the file-loaded reference runs in a
    // child. A complete read after all failures must still yield valid tables.
    const auto reference = (scratch / "file_routes.txt").string();
    const pid_t child = fork();
    if (child < 0)
      throw std::runtime_error("fork failed");
    if (child == 0) {
      try {
        flute::readLUT(powv.c_str(), post.c_str());
        std::ofstream(reference, std::ios::binary) << route_signature();
      } catch (...) {
        _exit(1);
      }
      _exit(0);
    }
    int status = 0;
    if (waitpid(child, &status, 0) != child || !WIFEXITED(status) ||
        WEXITSTATUS(status) != 0)
      throw std::runtime_error("file-loaded FLUTE reference failed");
    flute::readLUT_embedded();
    const std::string embedded = route_signature();
    if (embedded.empty() || embedded != slurp(reference))
      throw std::runtime_error("embedded and file-loaded FLUTE trees differ");
    // Four corners plus the center: the optimal tree is a middle bar with two
    // vertical stubs on each side, length 30.
    flute::DType x[] = {0, 10, 0, 10, 5}, y[] = {0, 0, 10, 10, 5};
    auto tree = flute::flute(5, x, y, 3);
    const bool ok = tree.length == 30 && flute::wirelength(tree) == 30;
    std::free(tree.branch);
    if (!ok)
      throw std::runtime_error("unexpected FLUTE wirelength");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
