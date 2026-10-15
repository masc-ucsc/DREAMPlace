// BSD-3-Clause. See LICENSE. Optional native frames and local browser playback.
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "dreamplace/native/synthetic.hpp"

namespace DreamPlace::native::synthetic {
namespace {
// Pairwise rectangle intersections, not bin-density overflow or union area.
// Sort by X to avoid checking distant cells. Fixed-fixed intersections do not
// violate movable-cell legality; zero-area IO terminals contribute nothing.
std::pair<size_t, double> overlaps(const Design& d) {
  std::vector<int> order;
  for (size_t i = 0; i < d.node_x.size(); ++i)
    if (d.node_size_x[i] > 0 && d.node_size_y[i] > 0) order.push_back(i);
  std::sort(order.begin(), order.end(), [&](int a, int b) { return d.node_x[a] < d.node_x[b]; });
  size_t pairs = 0;
  double area = 0;
  for (size_t i = 0; i < order.size(); ++i) {
    int a = order[i];
    double xh = d.node_x[a] + d.node_size_x[a];
    for (size_t j = i + 1; j < order.size() && d.node_x[order[j]] < xh; ++j) {
      int b = order[j];
      if (a >= d.num_movable && b >= d.num_movable) continue;
      double w = std::min(xh, d.node_x[b] + d.node_size_x[b]) - d.node_x[b];
      double h =
          std::min(d.node_y[a] + d.node_size_y[a], d.node_y[b] + d.node_size_y[b]) - std::max(d.node_y[a], d.node_y[b]);
      if (w > 1e-6 && h > 1e-6) {
        ++pairs;
        area += w * h;
      }
    }
  }
  return {pairs, area};
}
}  // namespace
void write_animation(const HhdsDesign& initial, const std::vector<AnimationFrame>& frames,
                     const std::filesystem::path& directory, const std::string& title, Initialization initialization) {
  if (frames.empty() || frames.front().stage != PlacementStage::initial ||
      frames.back().stage != PlacementStage::finished)
    throw std::invalid_argument("animation needs a complete placement trace");
  auto folder = directory / "animation_frames";
  std::filesystem::create_directories(folder);
  auto output = [](const std::filesystem::path& path) {
    std::ofstream f(path);
    f.exceptions(std::ios::badbit | std::ios::failbit);
    return f;
  };
  auto csv = output(directory / "animation.csv");
  auto concat = output(directory / "animation.ffconcat");
  csv << "frame,stage,gp_iteration,hpwl_um,legal,duration_seconds,overlap_pairs,pairwise_overlap_area_um2\n"
      << std::setprecision(17);
  concat << "ffconcat version 1.0\n";
  std::ostringstream data;
  data << std::setprecision(17) << '[';
  auto current = initial;
  std::string last;
  for (size_t i = 0; i < frames.size(); ++i) {
    const auto& frame = frames[i];
    if (frame.x.size() != initial.design.node_x.size() || frame.y.size() != frame.x.size())
      throw std::invalid_argument("animation snapshot has nonphysical nodes");
    current.design.node_x = frame.x;
    current.design.node_y = frame.y;
    std::ostringstream name;
    name << "animation_frames/frame_" << std::setfill('0') << std::setw(5) << i << ".png";
    last = name.str();
    auto [overlap_pairs, overlap_area] = overlaps(current.design);
    write_frame_png(initial, current, directory / last, title, initialization, frame, overlap_pairs);
    double duration = frame.stage == PlacementStage::global_placement ? 0.12 : 1.0;
    if (frame.stage == PlacementStage::finished) duration = 2.0;
    double wirelength = hpwl(current.design);
    bool legal = is_legal(current.design);
    csv << i << ',' << stage_name(frame.stage) << ',' << frame.iteration << ',' << wirelength << ',' << legal << ','
        << duration << ',' << overlap_pairs << ',' << overlap_area << '\n';
    concat << "file '" << last << "'\nduration " << duration << '\n';
    if (i) data << ',';
    data << "{\"src\":\"" << last << "\",\"phase\":\"" << stage_name(frame.stage) << "\",\"step\":" << frame.iteration
         << ",\"hpwl\":" << wirelength << ",\"legal\":" << (legal ? "true" : "false") << ",\"ms\":" << duration * 1000
         << ",\"pairs\":" << overlap_pairs << ",\"pairwise_overlap_area\":" << overlap_area << '}';
  }
  // The concat demuxer needs a repeated last image to honor its duration.
  concat << "file '" << last << "'\n";
  data << ']';
  auto html = output(directory / "animation.html");
  html << R"HTML(<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>DREAMPlace iteration playback</title>
<style>
body{font:16px system-ui,sans-serif;background:#f4f6f9;color:#1e2d41;max-width:1200px;margin:24px auto;padding:0 16px}
h1{font-size:24px;margin-bottom:8px}p{line-height:1.5}img{display:block;width:100%;height:auto;border:1px solid #ccd5e0}
.controls{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin:16px 0}
button,select{font:inherit;padding:8px 12px}input{flex:1;min-width:200px}#status{font-variant-numeric:tabular-nums}
</style>
<h1>DREAMPlace iteration playback</h1>
<p>The initial placement stays on the left. Actual sampled iterates appear on the right, followed by legalization
and detailed placement. Playback slows at phase boundaries. Positions are not interpolated.</p>
<img id="frame" alt="Initial placement and current DREAMPlace iteration" src="animation_frames/frame_00000.png">
<div class="controls"><button id="play" type="button">Pause</button>
<button id="previous" type="button">Previous</button><button id="next" type="button">Next</button>
<input id="seek" type="range" min="0" value="0" step="1" aria-label="Animation frame">
<label>Speed <select id="speed"><option value="0.5">0.5x</option><option value="1" selected>1x</option>
<option value="2">2x</option></select></label></div>
<p id="status" aria-live="polite"></p><p><a href="animation.csv">Per-frame metrics</a> ·
<a href="placement.png">Final placement PNG</a></p>
<script id="frames" type="application/json">)HTML"
       << data.str() << R"HTML(</script>
<script>
'use strict';
const frames=JSON.parse(document.getElementById('frames').textContent);
const img=document.getElementById('frame'),seek=document.getElementById('seek'),play=document.getElementById('play');
const speed=document.getElementById('speed'),status=document.getElementById('status');
let index=0,playing=true,timer;
seek.max=frames.length-1;
function show(){const f=frames[index];img.src=f.src;seek.value=index;
status.textContent=`Frame ${index+1}/${frames.length} · ${f.phase} · GP step ${f.step} · HPWL ${f.hpwl.toFixed(2)} µm · ${f.pairs} overlapping pairs · Legal: ${f.legal?'yes':'no'}`;}
function schedule(){clearTimeout(timer);if(playing)timer=setTimeout(()=>{index=(index+1)%frames.length;show();schedule();},frames[index].ms/Number(speed.value));}
function pause(){playing=false;play.textContent='Play';clearTimeout(timer);}
play.onclick=()=>{playing=!playing;play.textContent=playing?'Pause':'Play';schedule();};
seek.oninput=()=>{pause();index=Number(seek.value);show();};
document.getElementById('previous').onclick=()=>{pause();index=Math.max(0,index-1);show();};
document.getElementById('next').onclick=()=>{pause();index=Math.min(frames.length-1,index+1);show();};
speed.onchange=schedule;
show();schedule();
</script></html>
)HTML";
  csv.close();
  concat.close();
  html.close();
}
}  // namespace DreamPlace::native::synthetic
