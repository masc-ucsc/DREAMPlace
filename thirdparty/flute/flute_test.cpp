// BSD-3-Clause. See LICENSE.
#include "thirdparty/flute/flute.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

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
    // A complete read after all failures must still initialize valid tables.
    flute::readLUT("thirdparty/flute/lut.ICCAD2015/POWV9.dat",
                   "thirdparty/flute/lut.ICCAD2015/POST9.dat");
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
