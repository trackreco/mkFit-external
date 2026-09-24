// seedfind -- standalone driver for the binnor-based geometric seeder.
//
// Loads the geometry plugin and reads events the way mkFit.cc does, builds a
// SeedLayer per layer role, and times ONLY the seeding: the layer fill and the
// search, separately.  Event I/O and any truth evaluation are outside the
// timers.  Run from the standalone BUILD directory (the geometry plugin and its
// .bin are looked up there):
//
//   test-seedgeom/bin/seedfind --input-file <f.bin> [--geom CMS-phase2]
//        [--num-events N] [--reps R] [--dump quads.txt] [--phi-lin MODE MARG]
//        [--qbin-c CM] [--qbin-d CM]

#include "SeedFinder.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace mkfit;
using namespace mkfit::seeding;

namespace {
  using clk = std::chrono::steady_clock;
  double secs(clk::time_point a, clk::time_point b) { return std::chrono::duration<double>(b - a).count(); }

  // Production LayerOfHits axes: phi 2^8 = 256 N-bins, q up to 2^8 N-bins.
  using AxPhi8 = axis_pow2_u1<float, unsigned short, 16, 8>;
  using AxQ = axis<float, unsigned short, 16, 8>;
  using Layer = SeedLayer<AxPhi8, AxQ>;

  unsigned int n_q_bins(const LayerInfo &li, float qbin) {
    return std::max(1u, (unsigned int)std::ceil((li.zmax() - li.zmin()) / qbin));
  }

  void usage() {
    printf(
        "seedfind --input-file F [--geom G] [--num-events N] [--reps R] [--dump FILE]\n"
        "         [--phi-lin MODE MARG] [--qbin-c CM] [--qbin-d CM] [--layers A B C D]\n");
  }
}  // namespace

int main(int argc, char *argv[]) {
  std::string input, geom = "CMS-phase2", dump;
  int n_events = 10, reps = 1;
  int la = 0, lb = 1, lc = 2, ld = 3;
  float qbin_c = -1, qbin_d = -1;
  SeedParams P;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> const char * {
      if (i + 1 >= argc) {
        usage();
        exit(1);
      }
      return argv[++i];
    };
    if (a == "--input-file")
      input = next();
    else if (a == "--geom")
      geom = next();
    else if (a == "--num-events")
      n_events = atoi(next());
    else if (a == "--reps")
      reps = std::max(1, atoi(next()));
    else if (a == "--dump")
      dump = next();
    else if (a == "--phi-lin") {
      P.phi_lin = atoi(next());
      P.phi_lin_marg = atof(next());
    } else if (a == "--qbin-c")
      qbin_c = atof(next());
    else if (a == "--qbin-d")
      qbin_d = atof(next());
    else if (a == "--layers") {
      la = atoi(next());
      lb = atoi(next());
      lc = atoi(next());
      ld = atoi(next());
    } else {
      usage();
      return 1;
    }
  }
  if (input.empty()) {
    usage();
    return 1;
  }

  Config::geomPlugin = geom;
  execTrackerInfoCreatorPlugin(Config::geomPlugin, Config::TrkInfo, Config::ItrInfo);
  const TrackerInfo &ti = Config::TrkInfo;

  DataFile df;
  const int n_in_file = df.openRead(input, ti.n_layers(), ti.geom_version());
  n_events = std::min(n_events, n_in_file);

  const LayerInfo &lia = ti[la], &lib = ti[lb], &lic = ti[lc], &lid = ti[ld];
  if (qbin_c <= 0)
    qbin_c = lic.q_bin();
  if (qbin_d <= 0)
    qbin_d = lid.q_bin();
  // a is iterated whole and b is queried in phi alone, so one q bin each
  Layer ga(lia.zmin(), lia.zmax(), 1), gb(lib.zmin(), lib.zmax(), 1);
  Layer gc(lic.zmin(), lic.zmax(), n_q_bins(lic, qbin_c)), gd(lid.zmin(), lid.zmax(), n_q_bins(lid, qbin_d));
  printf("[seedfind] layers %d %d %d %d; phi bins %u; q bins c %u (%.2f cm) d %u (%.2f cm); phi_lin %d %.4f\n",
         la, lb, lc, ld, gc.n_phi_bins(), gc.n_q_bins(), qbin_c, gd.n_q_bins(), qbin_d, P.phi_lin,
         P.phi_lin_marg);

  std::vector<std::array<unsigned int, 5>> all_quads;
  std::vector<Quad> quads;
  SeedCounters tot;
  double t_fill = 0, t_find = 0;

  for (int iev = 0; iev < n_events; ++iev) {
    Event ev(iev, ti.n_layers());
    ev.read_in(df);

    const auto t0 = clk::now();
    ga.fill(ev.layerHits_[la]);
    gb.fill(ev.layerHits_[lb]);
    gc.fill(ev.layerHits_[lc]);
    gd.fill(ev.layerHits_[ld]);
    const auto t1 = clk::now();
    t_fill += secs(t0, t1);

    double best = std::numeric_limits<double>::max();
    SeedCounters cnt;
    for (int r = 0; r < reps; ++r) {
      quads.clear();
      cnt = SeedCounters();
      const auto s0 = clk::now();
      find_quads(P, ga, gb, gc, gd, quads, cnt);
      best = std::min(best, secs(s0, clk::now()));
    }
    t_find += best;
    tot.add(cnt);
    for (const auto &q : quads)
      all_quads.push_back({(unsigned int)iev, q[0], q[1], q[2], q[3]});
  }

  const double ne = n_events;
  printf("[seedfind] %d events, reps %d (min taken per event)\n", n_events, reps);
  printf("   doublets   %12.0f /ev\n", tot.doublets / ne);
  printf("   c touched  %12.0f /ev  (%.3f per doublet)\n", tot.c_touched / ne, (double)tot.c_touched / tot.doublets);
  printf("   triplets   %12.0f /ev\n", tot.triplets / ne);
  printf("   d touched  %12.0f /ev  (%.3f per triplet)\n", tot.d_touched / ne,
         (double)tot.d_touched / std::max(1L, tot.triplets));
  printf("   QUADS      %12.0f /ev\n", tot.quads / ne);
  printf("   fill       %12.3f ms/ev\n", 1e3 * t_fill / ne);
  printf("   find       %12.3f ms/ev  (%.1f ns per doublet)\n", 1e3 * t_find / ne, 1e9 * t_find / tot.doublets);

  if (!dump.empty()) {
    std::sort(all_quads.begin(), all_quads.end());
    FILE *f = fopen(dump.c_str(), "w");
    if (!f) {
      printf("cannot open %s\n", dump.c_str());
      return 1;
    }
    fprintf(f, "# event ia ib ic id -- original hit indices within layers %d %d %d %d; seedfind\n", la, lb, lc, ld);
    for (const auto &e : all_quads)
      fprintf(f, "%u %u %u %u %u\n", e[0], e[1], e[2], e[3], e[4]);
    fclose(f);
    printf("[seedfind] wrote %zu quads to %s\n", all_quads.size(), dump.c_str());
  }
  return 0;
}
