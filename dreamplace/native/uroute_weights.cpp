// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/uroute_weights.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <map>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace DreamPlace::native {
struct UrouteWeights::Data {
  struct Entry {
    std::vector<int> shape;
    std::span<const float> values;
  };
  void *mapping = MAP_FAILED;
  size_t size = 0;
  std::map<std::string, Entry> tensors;
  ~Data() {
    if (mapping != MAP_FAILED)
      munmap(mapping, size);
  }
};
namespace {
constexpr size_t limit = 256 * 1024 * 1024;
size_t align64(size_t n) { return (n + 63) / 64 * 64; }
void require(bool ok, const char *message) {
  if (!ok)
    throw std::invalid_argument(message);
}
} // namespace
UrouteWeights::UrouteWeights(const std::filesystem::path &path) {
  if constexpr (std::endian::native != std::endian::little)
    throw std::runtime_error(
        "URoute mapped f32 weights require a little-endian host");
  auto d = std::make_shared<Data>();
  int fd = open(path.c_str(), O_RDONLY);
  if (fd < 0)
    throw std::runtime_error("cannot open URoute weights: " + path.string());
  struct stat st{};
  if (fstat(fd, &st) != 0 || !S_ISREG(st.st_mode) || st.st_size < 16 ||
      st.st_size > static_cast<off_t>(limit)) {
    close(fd);
    throw std::invalid_argument("invalid URoute weights file size/type");
  }
  d->size = st.st_size;
  d->mapping = mmap(nullptr, d->size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (d->mapping == MAP_FAILED)
    throw std::runtime_error("cannot map URoute weights");
  const auto *bytes = static_cast<const char *>(d->mapping);
  require(std::memcmp(bytes, "UROUTEW1", 8) == 0,
          "unsupported URoute weights format");
  size_t pos = 8;
  auto read32 = [&] {
    require(pos + 4 <= d->size, "truncated URoute tensor table");
    std::uint32_t v;
    std::memcpy(&v, bytes + pos, 4);
    pos += 4;
    return v;
  };
  require(read32() == uroute_arch_version,
          "unsupported URoute architecture version");
  const auto count = read32();
  require(count >= 1 && count <= 1024, "invalid URoute tensor count");
  struct Descriptor {
    std::string name;
    std::vector<int> shape;
    size_t offset, count;
  };
  std::vector<Descriptor> table;
  for (std::uint32_t i = 0; i < count; ++i) {
    const auto len = read32();
    require(len >= 1 && len <= 128 && pos + len <= d->size,
            "invalid URoute tensor name");
    Descriptor desc{std::string(bytes + pos, len), {}, 0, 1};
    require(std::all_of(desc.name.begin(), desc.name.end(),
                        [](unsigned char c) { return c >= 32 && c < 127; }),
            "non-ASCII URoute tensor name");
    pos += (len + 3) / 4 * 4;
    const auto rank = read32();
    require(rank >= 1 && rank <= 4, "invalid URoute tensor rank");
    for (std::uint32_t j = 0; j < rank; ++j) {
      const auto dim = read32();
      require(dim >= 1 && dim <= 4096, "invalid URoute tensor dimension");
      desc.shape.push_back(dim);
      require(desc.count <= limit / 4 / dim, "oversized URoute tensor");
      desc.count *= dim;
    }
    const auto low = read32(), high = read32();
    const std::uint64_t offset = std::uint64_t(high) << 32 | low;
    require(offset <= d->size, "invalid URoute tensor offset");
    desc.offset = offset;
    table.push_back(std::move(desc));
  }
  size_t end = align64(pos);
  for (const auto &desc : table) {
    require(desc.offset == end && desc.count * 4 <= d->size - end,
            "invalid URoute tensor extent/order");
    const auto *data = reinterpret_cast<const float *>(bytes + desc.offset);
    std::span<const float> values(data, desc.count);
    require(std::all_of(values.begin(), values.end(),
                        [](float v) { return std::isfinite(v); }),
            "nonfinite URoute weights");
    require(
        d->tensors.emplace(desc.name, Data::Entry{desc.shape, values}).second,
        "duplicate URoute tensor");
    end = align64(end + desc.count * 4);
  }
  require(end == d->size, "trailing URoute weights data");
  data_ = d;
  // Freeze the complete architecture, not merely the shapes requested by one
  // execution. Unused or missing tensors are errors, never silently ignored.
  std::map<std::string, std::vector<int>> expected;
  auto conv = [&](const std::string &name, int in, int out, int kernel) {
    expected[name + ".weight"] = {out, in, kernel, kernel};
    expected[name + ".bias"] = {out};
  };
  conv("q_embed", 4, 64, 3);
  conv("s_embed", 5, 64, 3);
  for (int level = 0; level < 4; ++level) {
    const auto index = "." + std::to_string(level);
    for (const auto &prefix : {std::string("query"), std::string("support"),
                               std::string("decoder")}) {
      const int in = prefix == "decoder" ? 128 : 64;
      conv(prefix + index + ".conv1", in, 64, 3);
      conv(prefix + index + ".conv2", 64, 64, 3);
      if (in != 64)
        conv(prefix + index + ".skip", in, 64, 1);
    }
    conv("matches" + index + ".z", 128, 64, 3);
    conv("matches" + index + ".gate", 1, 1, 7);
    conv("matches" + index + ".out", 64, 64, 1);
  }
  conv("output", 64, 1, 1);
  require(expected.size() == d->tensors.size(),
          "URoute weight names do not match architecture v1");
  for (const auto &[name, shape] : expected)
    tensor(name, shape);
}
std::span<const float> UrouteWeights::tensor(const std::string &name,
                                             std::span<const int> shape) const {
  const auto found = data_->tensors.find(name);
  if (found == data_->tensors.end() ||
      !std::equal(shape.begin(), shape.end(), found->second.shape.begin(),
                  found->second.shape.end()))
    throw std::invalid_argument("URoute tensor name/shape mismatch: " + name);
  return found->second.values;
}
} // namespace DreamPlace::native
