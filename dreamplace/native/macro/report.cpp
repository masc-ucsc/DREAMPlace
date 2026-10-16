// BSD-3-Clause. See LICENSE.
#include "dreamplace/native/macro/report.hpp"
#include <cmath>
#include <limits>
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#include <stdexcept>
namespace DreamPlace::native::macro {
namespace {
using Writer = rapidjson::Writer<rapidjson::StringBuffer>;
void text(Writer &w, std::string_view s) { w.String(s.data(), s.size()); }
void key(Writer &w, const char *s) { w.Key(s); }
void coordinate(Writer &w, double v) {
  if (!std::isfinite(v) || std::abs(v) > 1e12 ||
      std::abs(v - std::round(v)) > 1e-6)
    throw std::invalid_argument("topology cache requires integer DBU geometry");
  w.Int64(std::llround(v));
}
void strength(Writer &w, Strength s) {
  key(w, "hard");
  w.Bool(s.hard);
  key(w, "weight");
  w.Double(s.weight);
}
void list(Writer &w, const std::vector<std::string> &v) {
  w.StartArray();
  for (const auto &s : v)
    text(w, s);
  w.EndArray();
}
using Value = rapidjson::Value;
const Value &field(const Value &v, const char *k) {
  if (!v.IsObject() || !v.HasMember(k))
    throw std::invalid_argument(std::string("topology cache missing field: ") +
                                k);
  return v[k];
}
std::string string(const Value &v) {
  if (!v.IsString())
    throw std::invalid_argument("topology cache expected string");
  return {v.GetString(), v.GetStringLength()};
}
double number(const Value &v) {
  if (!v.IsNumber() || !std::isfinite(v.GetDouble()))
    throw std::invalid_argument("topology cache expected finite number");
  return v.GetDouble();
}
int integer(const Value &v) {
  if (!v.IsInt())
    throw std::invalid_argument("topology cache expected integer");
  return v.GetInt();
}
bool boolean(const Value &v) {
  if (!v.IsBool())
    throw std::invalid_argument("topology cache expected boolean");
  return v.GetBool();
}
const Value &array(const Value &v) {
  if (!v.IsArray())
    throw std::invalid_argument("topology cache expected array");
  return v;
}
Strength read_strength(const Value &v) {
  return {boolean(field(v, "hard")), number(field(v, "weight"))};
}
Orientation orientation(const Value &v) {
  int o = integer(v);
  if (o < 0 || o > 3)
    throw std::invalid_argument("topology cache invalid orientation");
  return static_cast<Orientation>(o);
}
} // namespace
std::string report_json(const Result &r, const Design &d,
                        std::span<const NodeInfo> info) {
  validate(r.topology, d, info);
  rapidjson::StringBuffer buffer;
  Writer w(buffer);
  w.StartObject();
  key(w, "schema_version");
  w.Int(1);
  key(w, "kind");
  w.String("h2d_topology");
  key(w, "design");
  text(w, r.topology.design);
  key(w, "hier_delim");
  w.String(".");
  key(w, "macro_fingerprint");
  text(w, r.topology.macro_fingerprint);
  key(w, "netlist_hash");
  text(w, r.topology.netlist_hash);
  key(w, "provider");
  text(w, r.topology.provider);
  key(w, "parameters_key");
  text(w, parameter_key(r.resolved));
  key(w, "params");
  w.StartObject();
#define H2D_INT(n, d, h)                                                       \
  key(w, #n);                                                                  \
  w.Int(r.resolved.n);
#define H2D_REAL(n, d, h)                                                      \
  key(w, #n);                                                                  \
  w.Double(r.resolved.n);
#define H2D_BOOL(n, d, h)                                                      \
  key(w, #n);                                                                  \
  w.Bool(r.resolved.n);
#define H2D_STRING(n, d, h)                                                    \
  key(w, #n);                                                                  \
  text(w, r.resolved.n);
#include "dreamplace/native/macro/params.def"
#undef H2D_INT
#undef H2D_REAL
#undef H2D_BOOL
#undef H2D_STRING
  w.EndObject();
  key(w, "units");
  w.StartObject();
  key(w, "pitch_x");
  w.Double(d.layout.site_width);
  key(w, "pitch_y");
  w.Double(d.layout.row_height);
  w.EndObject();
  key(w, "legal");
  w.Bool(r.legal);
  key(w, "hpwl");
  w.Double(r.hpwl);
  key(w, "rudy_overflow");
  w.Double(r.rudy_overflow);
  key(w, "merit");
  w.Double(r.merit);
  key(w, "macros");
  w.StartArray();
  for (int i = 0; i < int(info.size()); ++i)
    if (info[i].macro) {
      w.StartObject();
      key(w, "name");
      text(w, info[i].name);
      key(w, "master");
      text(w, info[i].master);
      key(w, "width");
      coordinate(w, d.node_size_x[i]);
      key(w, "height");
      coordinate(w, d.node_size_y[i]);
      w.EndObject();
    }
  w.EndArray();
  key(w, "groups");
  w.StartArray();
  for (const auto &g : r.topology.groups) {
    w.StartObject();
    key(w, "id");
    text(w, g.id);
    key(w, "cols");
    w.Int(g.cols);
    key(w, "rows");
    w.Int(g.rows);
    key(w, "x");
    coordinate(w, g.x);
    key(w, "y");
    coordinate(w, g.y);
    key(w, "width");
    coordinate(w, g.width);
    key(w, "height");
    coordinate(w, g.height);
    key(w, "halo_x");
    coordinate(w, g.halo_x);
    key(w, "halo_y");
    coordinate(w, g.halo_y);
    key(w, "orientation");
    w.Int(int(g.orientation));
    key(w, "allowed_orientations");
    w.Uint(g.allowed_orientations);
    key(w, "fixed");
    w.Bool(g.fixed);
    key(w, "freeze_topology");
    w.Bool(g.freeze_topology);
    key(w, "members");
    w.StartArray();
    for (const auto &m : g.members) {
      w.StartObject();
      key(w, "name");
      text(w, m.name);
      key(w, "dx");
      coordinate(w, m.dx);
      key(w, "dy");
      coordinate(w, m.dy);
      key(w, "key_parent");
      text(w, m.key_parent);
      key(w, "key_signature");
      text(w, m.key_signature);
      key(w, "key_ordinal");
      w.Int(m.key_ordinal);
      key(w, "orientation");
      w.Int(int(m.orientation));
      w.EndObject();
    }
    w.EndArray();
    w.EndObject();
  }
  w.EndArray();
  key(w, "anchors");
  w.StartArray();
  for (const auto &a : r.topology.anchors) {
    w.StartObject();
    key(w, "id");
    text(w, a.id);
    key(w, "group");
    text(w, a.group);
    key(w, "side");
    w.Int(int(a.side));
    key(w, "tolerance");
    w.Double(a.tolerance);
    strength(w, a.strength);
    w.EndObject();
  }
  w.EndArray();
  key(w, "order");
  w.StartArray();
  for (const auto &e : r.topology.order) {
    w.StartObject();
    key(w, "id");
    text(w, e.id);
    key(w, "from");
    text(w, e.from);
    key(w, "to");
    text(w, e.to);
    key(w, "axis");
    w.Int(int(e.axis));
    strength(w, e.strength);
    w.EndObject();
  }
  w.EndArray();
  key(w, "fixed");
  w.StartArray();
  for (const auto &f : r.topology.fixed) {
    w.StartObject();
    key(w, "name");
    text(w, f.name);
    key(w, "x");
    coordinate(w, f.x);
    key(w, "y");
    coordinate(w, f.y);
    key(w, "orientation");
    w.Int(int(f.orientation));
    strength(w, f.strength);
    w.EndObject();
  }
  w.EndArray();
  key(w, "blockages");
  w.StartArray();
  for (const auto &b : r.topology.blockages) {
    w.StartObject();
    key(w, "id");
    text(w, b.id);
    key(w, "rect");
    w.StartArray();
    for (double v : {b.rect.xl, b.rect.yl, b.rect.xh, b.rect.yh})
      coordinate(w, v);
    w.EndArray();
    strength(w, b.strength);
    w.EndObject();
  }
  w.EndArray();
  key(w, "dataflow");
  w.StartArray();
  for (const auto &f : r.topology.dataflow) {
    w.StartObject();
    key(w, "from");
    text(w, f.from);
    key(w, "to");
    text(w, f.to);
    key(w, "weight");
    w.Double(f.weight);
    w.EndObject();
  }
  w.EndArray();
  key(w, "realized");
  w.StartObject();
  key(w, "plus");
  list(w, r.topology.realized.plus);
  key(w, "minus");
  list(w, r.topology.realized.minus);
  w.EndObject();
  key(w, "realized_graph");
  w.StartArray();
  for (const auto &e : r.topology.realized_graph) {
    w.StartObject();
    key(w, "from");
    text(w, e.from);
    key(w, "to");
    text(w, e.to);
    key(w, "id");
    text(w, e.id);
    key(w, "axis");
    w.Int(int(e.axis));
    w.EndObject();
  }
  w.EndArray();
  key(w, "violations");
  w.StartArray();
  for (const auto &v : r.violations) {
    w.StartObject();
    key(w, "id");
    text(w, v.id);
    key(w, "kind");
    text(w, v.kind);
    key(w, "amount");
    w.Double(v.amount);
    w.EndObject();
  }
  w.EndArray();
  key(w, "diagnostics");
  list(w, r.topology.diagnostics);
  key(w, "iterations");
  w.StartArray();
  for (const auto &i : r.iterations) {
    w.StartObject();
    key(w, "phases");
    w.StartObject();
    for (const auto &[name, seconds] : i.phases) {
      w.Key(name.c_str());
      w.Double(seconds);
    }
    w.EndObject();
    key(w, "gp_overflow");
    w.Double(i.gp_overflow);
    key(w, "gp_converged");
    w.Bool(i.gp_converged);
    key(w, "screened");
    w.Int(i.screened);
    key(w, "relaxed");
    list(w, i.relaxed);
    key(w, "slack");
    w.StartArray();
    for (const auto &s : i.slack) {
      w.StartObject();
      key(w, "id");
      text(w, s.id);
      key(w, "amount");
      w.Int64(s.amount);
      key(w, "weight");
      w.Double(s.weight);
      w.EndObject();
    }
    w.EndArray();
    key(w, "round");
    w.Int(i.round);
    key(w, "gp_iterations");
    w.Int(i.gp_iterations);
    key(w, "lp_solves");
    w.Int(i.lp_solves);
    key(w, "evaluations");
    w.Int(i.evaluations);
    key(w, "accepted_moves");
    w.Int(i.accepted_moves);
    key(w, "hpwl");
    w.Double(i.hpwl);
    key(w, "rudy_overflow");
    w.Double(i.rudy_overflow);
    key(w, "merit");
    w.Double(i.merit);
    key(w, "anchor_deviation");
    w.Double(i.anchor_deviation);
    key(w, "order_violation");
    w.Double(i.order_violation);
    key(w, "seconds");
    w.Double(i.seconds);
    key(w, "legal");
    w.Bool(i.legal);
    key(w, "selected");
    w.Bool(i.selected);
    w.EndObject();
  }
  w.EndArray();
  w.EndObject();
  if (!w.IsComplete())
    throw std::runtime_error("incomplete H2D report");
  return {buffer.GetString(), buffer.GetSize()};
}
Topology cached_topology(std::string_view json, const Design &d,
                         std::span<const NodeInfo> info, const Parameters &p) {
  rapidjson::Document v;
  v.Parse(json.data(), json.size());
  if (v.HasParseError() || !v.IsObject())
    throw std::invalid_argument("invalid topology cache JSON");
  if (integer(field(v, "schema_version")) != 1 ||
      string(field(v, "kind")) != "h2d_topology" ||
      string(field(v, "parameters_key")) != parameter_key(p) ||
      !boolean(field(v, "legal")))
    throw std::invalid_argument("incompatible topology cache header");
  Topology t;
  t.design = string(field(v, "design"));
  t.macro_fingerprint = string(field(v, "macro_fingerprint"));
  t.netlist_hash = string(field(v, "netlist_hash"));
  if (t.macro_fingerprint != fingerprint(d, info))
    throw std::invalid_argument("macro fingerprint changed");
  t.provider = string(field(v, "provider"));
  for (const auto &x : array(field(v, "groups")).GetArray()) {
    Group g;
    g.id = string(field(x, "id"));
    g.cols = integer(field(x, "cols"));
    g.rows = integer(field(x, "rows"));
    g.x = number(field(x, "x"));
    g.y = number(field(x, "y"));
    g.width = number(field(x, "width"));
    g.height = number(field(x, "height"));
    g.halo_x = number(field(x, "halo_x"));
    g.halo_y = number(field(x, "halo_y"));
    g.orientation = orientation(field(x, "orientation"));
    int mask = integer(field(x, "allowed_orientations"));
    if (mask < 1 || mask > 15)
      throw std::invalid_argument("bad cached orientation mask");
    g.allowed_orientations = mask;
    g.fixed = boolean(field(x, "fixed"));
    g.freeze_topology = boolean(field(x, "freeze_topology"));
    for (const auto &m : array(field(x, "members")).GetArray()) {
      Member member{string(field(m, "name")), number(field(m, "dx")),
                    number(field(m, "dy")),
                    orientation(field(m, "orientation"))};
      if (m.HasMember("key_signature")) {
        member.key_parent = string(field(m, "key_parent"));
        member.key_signature = string(field(m, "key_signature"));
        member.key_ordinal = integer(field(m, "key_ordinal"));
      }
      g.members.push_back(std::move(member));
    }
    t.groups.push_back(std::move(g));
  }
  for (const auto &x : array(field(v, "anchors")).GetArray()) {
    int side = integer(field(x, "side"));
    if (side < 0 || side > 3)
      throw std::invalid_argument("bad cached side");
    t.anchors.push_back({string(field(x, "id")), string(field(x, "group")),
                         static_cast<Side>(side), number(field(x, "tolerance")),
                         read_strength(x)});
  }
  for (const auto &x : array(field(v, "order")).GetArray()) {
    int axis = integer(field(x, "axis"));
    if (axis < 0 || axis > 1)
      throw std::invalid_argument("bad cached axis");
    t.order.push_back({string(field(x, "id")), string(field(x, "from")),
                       string(field(x, "to")), static_cast<Axis>(axis),
                       read_strength(x)});
  }
  for (const auto &x : array(field(v, "fixed")).GetArray())
    t.fixed.push_back({string(field(x, "name")), number(field(x, "x")),
                       number(field(x, "y")),
                       orientation(field(x, "orientation")), read_strength(x)});
  for (const auto &x : array(field(v, "blockages")).GetArray()) {
    const auto &rect = array(field(x, "rect"));
    if (rect.Size() != 4)
      throw std::invalid_argument("bad blockage rectangle");
    t.blockages.push_back(
        {string(field(x, "id")),
         {number(rect[0]), number(rect[1]), number(rect[2]), number(rect[3])},
         read_strength(x)});
  }
  for (const auto &x : array(field(v, "dataflow")).GetArray())
    t.dataflow.push_back({string(field(x, "from")), string(field(x, "to")),
                          number(field(x, "weight"))});
  const auto &realized = field(v, "realized");
  for (const auto &x : array(field(realized, "plus")).GetArray())
    t.realized.plus.push_back(string(x));
  for (const auto &x : array(field(realized, "minus")).GetArray())
    t.realized.minus.push_back(string(x));
  for (const auto &x : array(field(v, "diagnostics")).GetArray())
    t.diagnostics.push_back(string(x));
  if (v.HasMember("realized_graph"))
    for (const auto &x : array(v["realized_graph"]).GetArray()) {
      int axis = integer(field(x, "axis"));
      if (axis < 0 || axis > 1)
        throw std::invalid_argument("bad cached graph axis");
      t.realized_graph.push_back(
          {string(field(x, "from")), string(field(x, "to")),
           string(field(x, "id")), static_cast<Axis>(axis)});
    }
  validate(t, d, info);
  return t;
}
} // namespace DreamPlace::native::macro
