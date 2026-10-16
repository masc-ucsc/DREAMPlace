// BSD-3-Clause. See LICENSE.
#pragma once
#include "dreamplace/native/macro/h2d.hpp"
namespace DreamPlace::native::macro {
// Output/cache only. No CLI topology-input or agent-editable JSON path.
std::string report_json(const Result &, const Design &,
                        std::span<const NodeInfo>);
Topology cached_topology(std::string_view, const Design &,
                         std::span<const NodeInfo>, const Parameters &);
} // namespace DreamPlace::native::macro
