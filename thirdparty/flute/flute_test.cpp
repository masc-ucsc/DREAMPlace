// BSD-3-Clause. See LICENSE.
#include "thirdparty/flute/flute.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

int main() {
  try {
    // A missing table must be reported to the caller (it used to exit the whole
    // process) and must not leave the library marked as loaded.
    bool rejected = false;
    try {
      flute::readLUT("/nonexistent/POWV9.dat", "/nonexistent/POST9.dat");
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    if (!rejected)
      throw std::runtime_error("missing FLUTE table was accepted");
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
