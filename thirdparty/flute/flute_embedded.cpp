// BSD-3-Clause. See LICENSE.
#include <cstddef>

#include "flute.hpp"

namespace flute {
namespace lut_embedded {
// Defined in the flute_embed_gen output (flute_lut_data.cpp).
extern const char powv9[];
extern const std::size_t powv9_size;
extern const char post9[];
extern const std::size_t post9_size;
} // namespace lut_embedded

void readLUT_embedded() {
  readLUT(lut_embedded::powv9, lut_embedded::powv9_size, lut_embedded::post9,
          lut_embedded::post9_size);
}

} // namespace flute
