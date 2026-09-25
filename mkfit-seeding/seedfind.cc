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
//        [--qbin-c CM] [--qbin-d CM] [--arith ref|fast] [--margins REF.txt]
//
// --margins REF: per event, the symmetric difference between this run's quads
// and the list in REF, each differing quad evaluated in both arithmetics (see
// SeedMargins.h), plus how far the fast arithmetic moves every cut margin over
// all reference quads.

#include "SeedFinder.h"
#include "SeedFinderStaged.h"
#include "SeedFinderBMajor.h"
#include "SeedMargins.h"
#include "SeedStats.h"

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
#include <map>
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
        "         [--phi-lin MODE MARG] [--qbin-c CM] [--qbin-d CM] [--layers A B C D]\n"
        "         [--staged | --fuse | --bmajor] [--block N] [--lbin CM]\n"
        "         [--arith ref|fast|fastk] [--stats] [--margins REF] [--eps-cm E] [--eps-rad E] [--margins-print N]\n");
  }

  // The difference tool.  See SeedMargins.h.
  struct MarginStudy {
    std::map<unsigned int, std::vector<Quad>> ref;  // per event, sorted
    double eps_cm = 1e-4, eps_rad = 1e-6;
    int n_print = 50, n_printed = 0;
    long n_ref = 0, n_new = 0, n_ref_only = 0, n_new_only = 0;
    long n_edge = 0, n_far = 0, n_noflip = 0, n_selfbad = 0, n_self = 0;
    long n_ref_near = 0;           // reference quads with some decisive cut within eps: the baseline
    double worst_edge = 0;         // largest |m_ref| / eps among the edge flips
    std::vector<double> dm[MC_N];  // |m_ref - m_fast| over reference quads, both valid
    struct Worst {
      double dm = -1, pt = 0, cot = 0;
      unsigned int ev = 0;
      Quad q{};
    } worst[MC_N];

    bool load(const std::string &fn) {
      FILE *f = fopen(fn.c_str(), "r");
      if (!f)
        return false;
      char line[256];
      while (fgets(line, sizeof line, f)) {
        if (line[0] == '#')
          continue;
        unsigned int e, a, b, c, d;
        if (sscanf(line, "%u %u %u %u %u", &e, &a, &b, &c, &d) == 5)
          ref[e].push_back({a, b, c, d});
      }
      fclose(f);
      for (auto &kv : ref)
        std::sort(kv.second.begin(), kv.second.end());
      return true;
    }

    double scaled(int i, double m) const { return std::abs(m) / (margin_cut_info(i).rad ? eps_rad : eps_cm); }

    void print_eval(const char *tag, const QuadEval &r, const QuadEval &f) {
      printf("      %-5s pT %7.2f cot %+6.3f :", tag, 0.0114 * r.R, r.cot);
      for (int i = 0; i < MC_N; ++i) {
        if (!margin_cut_info(i).decisive && r.c[i].pass == f.c[i].pass)
          continue;
        const bool flip = r.c[i].valid && f.c[i].valid && r.c[i].pass != f.c[i].pass;
        const bool near = r.c[i].valid && scaled(i, r.c[i].m) < 10;
        if (!flip && !near)
          continue;
        printf(" %s%s[%s ref %+.3g fast %+.3g]", flip ? "*" : "", margin_cut_info(i).name,
               margin_cut_info(i).rad ? "rad" : "cm", r.c[i].m, f.c[i].m);
      }
      printf("\n");
    }

    // A: this run's arithmetic.  F: the one compared against the reference in
    // the precision table (A itself unless A is the reference).
    template <class A, class F, typename L>
    void event(unsigned int iev, const SeedParams &P, const L &ga, const L &gb, const L &gc, const L &gd,
               std::vector<Quad> quads) {
      std::sort(quads.begin(), quads.end());
      const std::vector<Quad> &rq = ref[iev];
      n_ref += rq.size();
      n_new += quads.size();
      std::vector<unsigned int> ia(ga.n()), ib(gb.n()), ic(gc.n()), id(gd.n());
      for (unsigned int k = 0; k < ga.n(); ++k)
        ia[ga.orig_[k]] = k;
      for (unsigned int k = 0; k < gb.n(); ++k)
        ib[gb.orig_[k]] = k;
      for (unsigned int k = 0; k < gc.n(); ++k)
        ic[gc.orig_[k]] = k;
      for (unsigned int k = 0; k < gd.n(); ++k)
        id[gd.orig_[k]] = k;
      auto ev2 = [&](const Quad &q, QuadEval &er, QuadEval &ef) {
        eval_quad<ArithRef>(P, ga, gb, gc, gd, ia[q[0]], ib[q[1]], ic[q[2]], id[q[3]], er);
        eval_quad<F>(P, ga, gb, gc, gd, ia[q[0]], ib[q[1]], ic[q[2]], id[q[3]], ef);
      };
      QuadEval er, ef;
      // self-check: the evaluator in this run's arithmetic accepts every quad this run found
      for (const auto &q : quads) {
        ev2(q, er, ef);
        ++n_self;
        if (!(std::is_same_v<A, ArithRef> ? er : ef).pass())
          ++n_selfbad;
      }
      // precision of the fast arithmetic, and the edge baseline, over the reference quads
      for (const auto &q : rq) {
        ev2(q, er, ef);
        bool near = false;
        for (int i = 0; i < MC_N; ++i) {
          // 1e30 is the curvature form's "no degeneracy test" sentinel
          if (er.c[i].valid && ef.c[i].valid && std::abs(er.c[i].m) < 1e29 && std::abs(ef.c[i].m) < 1e29) {
            const double d = std::abs(er.c[i].m - ef.c[i].m);
            dm[i].push_back(d);
            if (d > worst[i].dm)
              worst[i] = {d, 0.0114 * er.R, er.cot, iev, q};
          }
          near |= margin_cut_info(i).decisive && er.c[i].valid && scaled(i, er.c[i].m) <= 1;
        }
        n_ref_near += near;
      }
      // the symmetric difference
      std::vector<Quad> only_ref, only_new;
      std::set_difference(rq.begin(), rq.end(), quads.begin(), quads.end(), std::back_inserter(only_ref));
      std::set_difference(quads.begin(), quads.end(), rq.begin(), rq.end(), std::back_inserter(only_new));
      n_ref_only += only_ref.size();
      n_new_only += only_new.size();
      for (int side = 0; side < 2; ++side) {
        for (const auto &q : side == 0 ? only_ref : only_new) {
          ev2(q, er, ef);
          // side 0: the reference accepts it and this run does not
          const QuadEval &acc = side == 0 ? er : (std::is_same_v<A, ArithRef> ? er : ef);
          const QuadEval &rej = side == 0 ? (std::is_same_v<A, ArithRef> ? er : ef) : er;
          double closest = 1e30;
          bool flip = false;
          for (int i = 0; i < MC_N; ++i) {
            if (!margin_cut_info(i).decisive || !acc.c[i].valid)
              continue;
            if (acc.c[i].pass && !(rej.c[i].valid && rej.c[i].pass)) {
              flip = true;
              closest = std::min(closest, scaled(i, er.c[i].valid ? er.c[i].m : ef.c[i].m));
            }
          }
          const char *cls;
          if (!acc.pass() || rej.pass() || !flip) {
            ++n_noflip;
            cls = "NOFLIP";
          } else if (closest <= 1) {
            ++n_edge;
            worst_edge = std::max(worst_edge, closest);
            cls = "edge";
          } else {
            ++n_far;
            cls = "FAR";
          }
          if (n_printed < n_print || cls[0] != 'e') {
            ++n_printed;
            printf("   %-6s %s ev %u quad %u %u %u %u  closest flip %.3g eps\n", cls,
                   side == 0 ? "ref-only" : "new-only", iev, q[0], q[1], q[2], q[3], flip ? closest : -1.0);
            print_eval("", er, ef);
          }
        }
      }
    }

    void report() const {
      printf("[margins] eps %.3g cm, %.3g rad\n", eps_cm, eps_rad);
      printf("   reference quads %ld, this run %ld; ref-only %ld, new-only %ld\n", n_ref, n_new, n_ref_only,
             n_new_only);
      printf("   differing quads: edge %ld (worst %.3g eps), FAR %ld, NOFLIP %ld\n", n_edge, worst_edge, n_far,
             n_noflip);
      printf("   self-check: %ld of %ld found quads rejected by this run's own evaluator\n", n_selfbad, n_self);
      printf("   baseline: %ld of %ld reference quads (%.3f %%) have a decisive cut within eps\n", n_ref_near, n_ref,
             100.0 * n_ref_near / std::max(1L, n_ref));
      printf("   |margin_ref - margin_fast| over reference quads, per cut:\n");
      printf("     %-11s %4s %8s %11s %11s %11s   %s\n", "cut", "unit", "n", "median", "p99", "max",
             "worst quad: pT, cot, ev, indices");
      for (int i = 0; i < MC_N; ++i) {
        std::vector<double> v = dm[i];
        if (v.empty())
          continue;
        std::sort(v.begin(), v.end());
        const Worst &w = worst[i];
        printf("     %-11s %4s %8zu %11.3g %11.3g %11.3g", margin_cut_info(i).name, margin_cut_info(i).rad ? "rad" : "cm",
               v.size(), v[v.size() / 2], v[(v.size() * 99) / 100], v.back());
        if (v.back() > 0)
          printf("   %7.2f %+6.3f  %u  %u %u %u %u", w.pt, w.cot, w.ev, w.q[0], w.q[1], w.q[2], w.q[3]);
        printf("\n");
      }
    }
  };
}  // namespace

int main(int argc, char *argv[]) {
  std::string input, geom = "CMS-phase2", dump;
  int n_events = 10, reps = 1;
  bool staged = false, fuse = false, bmajor = false;
  int arith = 0;  // 0 ref, 1 fast, 2 fastk
  std::string margins_ref;
  MarginStudy MS;
  bool stats = false;
  SeedStats SS;
  unsigned int block = 64;
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
    else if (a == "--staged")
      staged = true;
    else if (a == "--bmajor")
      staged = bmajor = true;
    else if (a == "--fuse")
      staged = fuse = true;
    else if (a == "--block")
      block = std::max(1, atoi(next()));
    else if (a == "--arith") {
      const std::string v = next();
      if (v == "ref")
        arith = 0;
      else if (v == "fast")
        arith = 1;
      else if (v == "fastk")
        arith = 2;
      else {
        usage();
        return 1;
      }
    } else if (a == "--stats")
      stats = true;
    else if (a == "--margins")
      margins_ref = next();
    else if (a == "--eps-cm")
      MS.eps_cm = atof(next());
    else if (a == "--eps-rad")
      MS.eps_rad = atof(next());
    else if (a == "--margins-print")
      MS.n_print = atoi(next());
    else if (a == "--dump")
      dump = next();
    else if (a == "--phi-lin") {
      P.phi_lin = atoi(next());
      P.phi_lin_marg = atof(next());
    } else if (a == "--lbin")
      P.lbin = atof(next());
    else if (a == "--qbin-c")
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
  if (arith && !staged) {
    printf("[seedfind] --arith fast/fastk needs --staged, --fuse or --bmajor\n");
    return 1;
  }
  if (!margins_ref.empty() && !MS.load(margins_ref)) {
    printf("[seedfind] cannot read %s\n", margins_ref.c_str());
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
  StagedWork work;
  BMajorWork bwork;
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
      if (bmajor && arith == 1)
        find_quads_bmajor<ArithFast>(P, ga, gb, gc, gd, quads, cnt, work, bwork, block);
      else if (bmajor && arith == 2)
        find_quads_bmajor<ArithFastK>(P, ga, gb, gc, gd, quads, cnt, work, bwork, block);
      else if (bmajor)
        find_quads_bmajor<ArithRef>(P, ga, gb, gc, gd, quads, cnt, work, bwork, block);
      else if (staged && arith == 1)
        find_quads_staged<ArithFast>(P, ga, gb, gc, gd, quads, cnt, work, block, fuse);
      else if (staged && arith == 2)
        find_quads_staged<ArithFastK>(P, ga, gb, gc, gd, quads, cnt, work, block, fuse);
      else if (staged)
        find_quads_staged<ArithRef>(P, ga, gb, gc, gd, quads, cnt, work, block, fuse);
      else
        find_quads(P, ga, gb, gc, gd, quads, cnt);
      best = std::min(best, secs(s0, clk::now()));
    }
    t_find += best;
    tot.add(cnt);
    if (stats)
      seed_stats(P, ga, gb, gc, gd, SS);
    if (!margins_ref.empty()) {
      if (arith == 1)
        MS.event<ArithFast, ArithFast>(iev, P, ga, gb, gc, gd, quads);
      else if (arith == 2)
        MS.event<ArithFastK, ArithFastK>(iev, P, ga, gb, gc, gd, quads);
      else
        MS.event<ArithRef, ArithFastK>(iev, P, ga, gb, gc, gd, quads);
    }
    for (const auto &q : quads)
      all_quads.push_back({(unsigned int)iev, q[0], q[1], q[2], q[3]});
  }

  const double ne = n_events;
  printf("[seedfind] %d events, reps %d (min taken per event), arith %s, %s\n", n_events, reps, arith == 2 ? "fastk" : arith == 1 ? "fast" : "ref",
         staged ? ((bmajor ? "b-major, block " : fuse ? "fused, block " : "staged, block ") + std::to_string(block)).c_str()
                : "scalar");
  printf("   doublets   %12.0f /ev\n", tot.doublets / ne);
  printf("   c touched  %12.0f /ev  (%.3f per doublet)\n", tot.c_touched / ne, (double)tot.c_touched / tot.doublets);
  printf("   triplets   %12.0f /ev\n", tot.triplets / ne);
  printf("   d touched  %12.0f /ev  (%.3f per triplet)\n", tot.d_touched / ne,
         (double)tot.d_touched / std::max(1L, tot.triplets));
  printf("   QUADS      %12.0f /ev\n", tot.quads / ne);
  printf("   fill       %12.3f ms/ev\n", 1e3 * t_fill / ne);
  printf("   find       %12.3f ms/ev  (%.1f ns per doublet)\n", 1e3 * t_find / ne, 1e9 * t_find / tot.doublets);
  if (staged) {
    // the counters are those of the LAST repetition of each event
    const char *nm[7] = {"1 doublets", "2 windows", bmajor ? "3 c-lists" : "3 c-cands", bmajor ? "4 c-search" : "4 triplets", "5 helix", "6 d-cands", "7 quads"};
    double ts = 0;
    for (double t : tot.t_stage)
      ts += t;
    for (int i = 0; i < 7; ++i)
      printf("     stage %-11s %9.3f ms/ev  %5.1f %%\n", nm[i], 1e3 * tot.t_stage[i] / ne,
             100 * tot.t_stage[i] / ts);
  }

  if (!margins_ref.empty())
    MS.report();
  if (stats)
    SS.report();

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
