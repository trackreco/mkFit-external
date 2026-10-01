// seedsurf -- driver for the surface finder (SeedSurf.h): quad seeds on any
// four pixel layers, barrel or disc, with truth.  Run from the standalone
// BUILD directory, as seedfind:
//
//   test-seedgeom/bin/seedsurf --input-file F.bin [--num-events N]
//        [--first-event N]   start after N events
//        --pattern A B C D [--pattern ...]   mkFit layer numbers; a pattern with
//                                            +z discs (16-27) also runs mirrored (38-49)
//        [--pt-min GEV] [--d0-max CM] [--zv CM] [--marg-b RAD]
//        [--win-c PHI Q] [--win-d PHI Q] [--no-cut] [--win-scale F]
//        [--pattern-win A B C D PHI_C Q_C PHI_D Q_D]   a pattern with its own windows
//        [--pattern-dwin A B C D PHI_C Q_C APHI_D BPHI_D AQ_D BQ_D]   ... with d windows a + b / pT_est
//        [--pattern-sref A B C D S]   for a listed pattern: b scaled by each candidate's c-d path length / S (cm)
//        [--pattern-dwin-eta A B C D LO HI APHI_D BPHI_D AQ_D BQ_D]   for a listed pattern: d windows in an |eta| slice
//                      of the triplet's helix (outside every slice the pattern's own d windows apply)
//        [--own DELTA] [--own-skip K] [--own-skip-ot K]   phase-space ownership (SurfOwnership in SeedSurf.h), margin in cm
//        [--dedup N]   over all patterns, keep a quad only if it shares < N hits with every better kept quad;
//                      better = fewer outer-tracker layers in the pattern, then smaller
//                      (dq_c/q_c)^2 + (dphi_d/w_phi_d)^2 + (dq_d/w_q_d)^2
//        [--bind CM]   labels bound to geometry: needs SimHitStates in the sample
//        [--chain H] [--chain-holes-ot K] [--chain-hole-always] [--chain-any] [--chain-start-holes K] [--chain-lead-only] [--chain-fast] [--chain-batch] [--chain-fast-check] [--chain-phases]   the feed-forward chain (SurfChain) in
//                      place of the pattern list; the patterns then give window tables and the denominator;
//                      --chain-batch runs the batched float finder (SeedSurfBatch.h) on the same configuration;
//                      --chain-batch-d N its stage d prediction: 0 direct from hit c, 1 one-point cubic, 2 two-point Hermite
//        [--truth OUT.txt] [--resid OUT.txt] [--dump quads.txt] [--quad-dump OUT.txt] [--eta-max E]
//        batch finder fake rejection (SurfChainBatch, README "Fakes in the finder"):
//        [--fk-score S]   a quad's (dq_c/w)^2 + (dphi_c/w)^2 + (dphi_d/w)^2 + (dq_d/w)^2 below S
//        [--shape-win L BINW N LO_0 HI_0 ... LO_N-1 HI_N-1]   barrel pixel layer L: the kept band of the
//                      cluster length along z per |cot theta| bin of width BINW (windows-D121/shape.txt)
//        [--fk-shape]     apply the shape bands
//        [--fk-ot2 F] [--ot2-win APHI BPHI AQ BQ]   a quad with d on OT1-P needs an OT2-P hit within
//                      F x (a + b / pT) in phi [rad] and z [cm] for the helix through b, c, d
//        [--attach-ot1 F] [--ot1-win APHI BPHI AQ BQ]   after the cleaning, for each quad with a pixel d
//                      in OT1's acceptance, the best OT1-P hit (SurfChainBatch::next_hit) if it is within
//                      F x (a + b / pT) in phi and z, as information; with its truth and time
//        [--margins REF.txt] [--margins-print N]   chain only: per event, the symmetric difference of the
//                      chain's quads (before cleaning) against REF, a --dump of a run with the same pattern
//                      options and events; each differing quad's cuts evaluated in double, with margins
//
// --truth: per pattern and for the union of all patterns, by |eta| of the sim
//   track: findable (a hit labelled with the track in each of the pattern's
//   layers, pT > pT_min, production within D0_max of the beam spot in xy),
//   found (a quad whose four hits carry its label); and the quads by the
//   consistent rule (true: four equal labels; fake: max count + unlabelled < 4;
//   undecidable otherwise), binned in the quad's own |eta| (a-d line), and in
//   pT_est.
// --resid: for every findable track and every combination of its labelled hits
//   in the pattern's layers, what each cut sees (SurfEval), so the windows can
//   be set from true quads.  Each row also carries the definite crossings of
//   other layers on its a-d line, before a (lead) and between a and d (inner),
//   so a fit can be restricted to what the chain builds (inner = 0), and the
//   azimuth of the d hit and of the track, and the d hit's radius.

#include "SeedSurf.h"
#include "SeedSurfBatch.h"
#include "RecoTracker/MkFitCore/interface/MkSeeder.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/radix_sort.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace mkfit;
using namespace mkfit::seeding;

namespace {
  using clk = std::chrono::steady_clock;
  double secs(clk::time_point a, clk::time_point b) { return std::chrono::duration<double>(b - a).count(); }

  const char *lname(int l, char *buf) {
    if (l <= 3)
      snprintf(buf, 8, "B%d", l + 1);
    else {
      int m = l >= 38 ? l - 22 : l;
      if (m >= 16 && m <= 23)
        snprintf(buf, 8, "F%d%s", m - 15, l >= 38 ? "-" : "");
      else if (m >= 24 && m <= 27)
        snprintf(buf, 8, "E%d%s", m - 23, l >= 38 ? "-" : "");
      else
        snprintf(buf, 8, "L%d", l);
    }
    return buf;
  }

  struct Pattern {
    std::array<int, 4> l;
    std::string name;
    float win[4] = {-1, -1, -1, -1};  // phi_c q_c phi_d q_d; < 0: the global ones
    float bwin[2] = {0, 0};           // b_phi_d, b_q_d: d windows a + b / pT_est
    float sref = 0;                   // > 0: b scaled by the candidate's c-d path length / sref
    std::vector<SurfParams::EtaWin> eta_win;  // d windows per |eta| slice of the triplet
    bool dynamic = false;             // made by the chain for a layer combination no pattern lists
  };

  constexpr int kNb = 20;
  constexpr double kEtaW = 0.2;
  constexpr int kNpe = 11;
  constexpr double kPeEdge[kNpe + 1] = {0, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0, 1e9};

  struct Stats {
    double den[kNb] = {}, found[kNb] = {}, dup[kNb] = {}, multi[kNb] = {};
    double q_all[kNb] = {}, q_true[kNb] = {}, q_fake[kNb] = {}, q_undec[kNb] = {};
    double pe_all[kNpe] = {}, pe_true[kNpe] = {}, pe_fake[kNpe] = {};
    // the union only: findable, found and extra true quads by the sim track's pT (kPeEdge), in four |eta| regions
    double sp_den[4][kNpe] = {}, sp_found[4][kNpe] = {}, sp_dup[4][kNpe] = {};
    // the union only: kept quads, true and fake, by their own pT_est, in four regions of their own |eta|
    double qp_all[4][kNpe] = {}, qp_true[4][kNpe] = {}, qp_fake[4][kNpe] = {};
    double own_rejected = 0, quads = 0, doublets = 0, triplets = 0, c_touched = 0, d_touched = 0, t_find = 0;
  };

  int ebin(double aeta) {
    const int b = (int)(aeta / kEtaW);
    return b < 0 ? 0 : (b >= kNb ? -1 : b);
  }
  int pbin(double pte) {
    int b = 0;
    while (b < kNpe - 1 && pte >= kPeEdge[b + 1])
      ++b;
    return b;
  }
}  // namespace

int main(int argc, char *argv[]) {
  std::string input, geom = "CMS-phase2", truth_out, resid_out, dump_out, qdump_out;
  int n_events = 10, first_event = 0;
  double eta_max = 4.0;
  double win_scale = 1.0;  // multiplies every c and d window
  double own_delta = -1;   // phase-space ownership margin [cm], < 0 off
  double bind_cm = -1;     // > 0: labels bound to the sim hit position
  int own_skip = 0;        // ownership: skipped definite layers allowed
  int own_skip_ot = -1;    // ... for patterns with an outer-tracker layer; < 0: own_skip
  int own_debug = 0;       // print this many true quads the ownership rejects
  int chain_holes = -1;    // >= 0: the feed-forward chain (SurfChain) instead of the patterns, this many holes
  int chain_hole_always = 0, chain_holes_ot = 0, chain_any = 0, chain_start_holes = -1, chain_lead_only = 0;
  int chain_fast = 0;  // --chain-fast: the float kernels (K2) in the chain
  int chain_batch = 0; // --chain-batch: the batched float finder (SurfChainBatch)
  int chain_batch_d = 0; // --chain-batch-d: its stage d prediction (0 direct, 1 one-point cubic, 2 two-point Hermite)
  int chain_phases = 0;  // --chain-phases: time the chain's phases
  float fk_score = 0, fk_ot2 = 0;
  float ot2_win[4] = {-8.1e-4f, 5.92e-3f, 0.3084f, 0.0909f};  // q97 of true quads, events 0-39
  float ot2_phimin = 1.31e-3f;  // --ot2-phimin: the floor of its phi term, the q97 above 3 GeV
  float ot1_win[4] = {2.16e-3f, 9.99e-3f, 0.2318f, 0.3566f};  // q97 of true quads, events 0-39
  float attach_ot1 = 0;
  bool fk_shape = false;
  SurfChainBatch::ShapeTab shape_tab[4];
  int dedup_n = 0;         // > 0: over all patterns, drop a quad sharing >= N hits with a better kept one
  std::string margins_ref;  // --margins: the reference quad list
  int margins_print = 10;
  SurfParams P;
  std::vector<Pattern> pats;

  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> const char * {
      if (i + 1 >= argc) {
        printf("missing value for %s\n", a.c_str());
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
    else if (a == "--first-event")
      first_event = atoi(next());
    else if (a == "--pattern") {
      Pattern p;
      for (int k = 0; k < 4; ++k)
        p.l[k] = atoi(next());
      pats.push_back(p);
    } else if (a == "--pattern-win") {
      // a pattern with its own windows: A B C D phi_c q_c phi_d q_d
      Pattern p;
      for (int k = 0; k < 4; ++k)
        p.l[k] = atoi(next());
      for (int k = 0; k < 4; ++k)
        p.win[k] = atof(next());
      pats.push_back(p);
    } else if (a == "--pattern-dwin") {
      // A B C D phi_c q_c a_phi_d b_phi_d a_q_d b_q_d: d windows a + b / pT_est
      Pattern p;
      for (int k = 0; k < 4; ++k)
        p.l[k] = atoi(next());
      p.win[0] = atof(next());
      p.win[1] = atof(next());
      p.win[2] = atof(next());
      p.bwin[0] = atof(next());
      p.win[3] = atof(next());
      p.bwin[1] = atof(next());
      pats.push_back(p);
    } else if (a == "--pattern-sref") {
      // A B C D S: for a pattern already listed, scale its d windows' b terms by s_cd / S
      std::array<int, 4> l;
      for (int k = 0; k < 4; ++k)
        l[k] = atoi(next());
      const float sr = atof(next());
      bool found = false;
      for (auto &p : pats)
        if (p.l == l)
          p.sref = sr, found = true;
      if (!found) {
        printf("--pattern-sref: pattern %d %d %d %d not listed before it\n", l[0], l[1], l[2], l[3]);
        return 1;
      }
    } else if (a == "--pattern-dwin-eta") {
      // A B C D ETA_LO ETA_HI APHI_D BPHI_D AQ_D BQ_D: for a listed pattern, d windows in an |eta| slice
      std::array<int, 4> l;
      for (int k = 0; k < 4; ++k)
        l[k] = atoi(next());
      SurfParams::EtaWin w;
      w.lo = atof(next()), w.hi = atof(next());
      w.aphi = atof(next()), w.bphi = atof(next()), w.aq = atof(next()), w.bq = atof(next());
      bool found = false;
      for (auto &p : pats)
        if (p.l == l)
          p.eta_win.push_back(w), found = true;
      if (!found) {
        printf("--pattern-dwin-eta: pattern %d %d %d %d not listed before it\n", l[0], l[1], l[2], l[3]);
        return 1;
      }
    } else if (a == "--win-scale")
      win_scale = atof(next());
    else if (a == "--pt-min")
      P.pt_min = atof(next());
    else if (a == "--d0-max")
      P.d0_max = atof(next());
    else if (a == "--zv")
      P.zv = atof(next());
    else if (a == "--marg-b")
      P.marg_b = atof(next());
    else if (a == "--win-c") {
      P.phi_c = atof(next());
      P.q_c = atof(next());
    } else if (a == "--win-d") {
      P.phi_d = atof(next());
      P.q_d = atof(next());
    } else if (a == "--no-cut")
      P.no_cut = 1;
    else if (a == "--truth")
      truth_out = next();
    else if (a == "--resid")
      resid_out = next();
    else if (a == "--dump")
      dump_out = next();
    else if (a == "--quad-dump")
      qdump_out = next();
    else if (a == "--own")
      own_delta = atof(next());
    else if (a == "--chain")
      chain_holes = atoi(next());
    else if (a == "--chain-start-holes")
      chain_start_holes = atoi(next());
    else if (a == "--chain-lead-only")
      chain_lead_only = 1;
    else if (a == "--chain-fast")
      chain_fast = 1;
    else if (a == "--chain-batch")
      chain_batch = 1;
    else if (a == "--chain-batch-d")
      chain_batch_d = atoi(next());
    else if (a == "--chain-fast-check")
      g_surf_fast_check.on = true;
    else if (a == "--chain-phases")
      chain_phases = 1;
    else if (a == "--fk-score")
      fk_score = atof(next());
    else if (a == "--fk-shape")
      fk_shape = true;
    else if (a == "--fk-ot2")
      fk_ot2 = atof(next());
    else if (a == "--ot2-phimin")
      ot2_phimin = atof(next());
    else if (a == "--ot2-win")
      for (int k = 0; k < 4; ++k)
        ot2_win[k] = atof(next());
    else if (a == "--attach-ot1")
      attach_ot1 = atof(next());
    else if (a == "--ot1-win")
      for (int k = 0; k < 4; ++k)
        ot1_win[k] = atof(next());
    else if (a == "--shape-win") {
      const int l = atoi(next());
      const float bw = atof(next());
      const int nb = atoi(next());
      if (l < 0 || l > 3 || bw <= 0 || nb < 1) {
        printf("--shape-win: bad layer %d, bin width %g or bin count %d\n", l, bw, nb);
        return 1;
      }
      SurfChainBatch::ShapeTab &T = shape_tab[l];
      T.inv_bw = 1.0f / bw;
      T.lo.resize(nb), T.hi.resize(nb);
      for (int b = 0; b < nb; ++b)
        T.lo[b] = atoi(next()), T.hi[b] = atoi(next());
    }
    else if (a == "--chain-any")
      chain_any = 1;
    else if (a == "--chain-hole-always")
      chain_hole_always = 1;
    else if (a == "--chain-holes-ot")
      chain_holes_ot = atoi(next());
    else if (a == "--own-debug")
      own_debug = atoi(next());
    else if (a == "--own-skip")
      own_skip = atoi(next());
    else if (a == "--own-skip-ot")
      own_skip_ot = atoi(next());
    else if (a == "--dedup")
      dedup_n = atoi(next());
    else if (a == "--margins")
      margins_ref = next();
    else if (a == "--margins-print")
      margins_print = atoi(next());
    else if (a == "--bind") {
      // truth BOUND to geometry: a hit keeps its label only if its sim hit (SimHitStates) is within
      // CM of it (pixels; outer tracker: 2.6 cm, the longest strip half-length)
      bind_cm = atof(next());
      Config::readSimHitStates = true;
    } else if (a == "--eta-max")
      eta_max = atof(next());
    else {
      printf("unknown option %s; see the header of seedsurf.cc\n", a.c_str());
      return 1;
    }
  }
  if (input.empty() || pats.empty()) {
    printf("need --input-file and at least one --pattern\n");
    return 1;
  }
  // mirror patterns with +z discs
  {
    std::vector<Pattern> all;
    for (auto p : pats) {
      all.push_back(p);  // the mirror keeps the windows
      bool has = false;
      Pattern m = p;
      for (int &l : m.l)
        if (l >= 16 && l <= 27)
          l += 22, has = true;
      if (has)
        all.push_back(m);
    }
    pats = all;
    char b[8];
    for (auto &p : pats)
      for (int k = 0; k < 4; ++k)
        p.name += std::string(k ? " " : "") + lname(p.l[k], b);
  }

  Config::geomPlugin = geom;
  execTrackerInfoCreatorPlugin(Config::geomPlugin, Config::TrkInfo, Config::ItrInfo);
  const TrackerInfo &ti = Config::TrkInfo;
  DataFile df;
  const int n_in_file = df.openRead(input, ti.n_layers(), ti.geom_version());
  if (first_event > 0) {
    df.skipNEvents(first_event);
    printf("[seedsurf] skipped the first %d events\n", first_event);
  }
  n_events = std::min(n_events, n_in_file - first_event);

  SurfOwnership OWN;
  {
    std::set<int> extra;
    for (const auto &p : pats)
      for (int l : p.l)
        extra.insert(l);
    OWN.setup(ti, extra);
  }
  OWN.delta = own_delta;
  OWN.max_skip = own_skip;
  OWN.max_skip_ot = own_skip_ot;
  // the seeder: its layers, its two finders (+z, -z) and the cleaning; the configuration is built here
  MkSeeder seeder;
  SeedEventOfHits &layers = seeder.hits();
  for (const auto &p : pats)
    for (int l : p.l)
      layers.add_layer(l, ti[l], ti[l].is_barrel() ? 2.0 : 1.0);

  // the chain: every layer it can visit, window tables from the listed patterns
  SurfChain CH[2];
  if (chain_holes >= 0) {
    std::set<int> have;
    for (const auto &p : pats)
      for (int l : p.l)
        have.insert(l);
    for (int l : {0, 1, 2, 3})
      have.insert(l);
    for (int l = 16; l <= 27; ++l)
      have.insert(l), have.insert(l + 22);
    for (int l : have)
      layers.add_layer(l, ti[l], ti[l].is_barrel() ? 2.0 : 1.0);
    {
      SurfOwnership O2;
      O2.setup(ti, have);
      OWN.env = O2.env;  // the chain's crossing test sees every layer
    }
    if (OWN.delta < 0)
      OWN.delta = 0.2;
    for (int sd = 0; sd < 2; ++sd) {
      SurfChain &C = CH[sd];
      C.P = P;
      C.P.phi_c *= win_scale, C.P.q_c *= win_scale, C.P.phi_d *= win_scale, C.P.q_d *= win_scale;
      C.max_holes = chain_holes;
      C.hole_always = chain_hole_always;
      C.max_holes_ot = chain_holes_ot;
      C.known_only = !chain_any;
      C.start_holes = chain_start_holes;
      C.lead_only = chain_lead_only;
      C.fast = chain_fast;
      C.phases = chain_phases;
      C.setup(OWN, sd == 0 ? 1 : -1, have);
    }
    printf("[seedsurf] CHAIN: max holes %d (into OT: %d)%s, crossing margin %.3f cm; start pairs %zu / %zu (+z / -z)\n",
           chain_holes, chain_holes_ot, chain_hole_always ? ", hole also after a found hit" : "", OWN.delta,
           CH[0].starts.size(), CH[1].starts.size());
  }
  if (own_delta >= 0 && chain_holes < 0)
    printf("[seedsurf] phase-space ownership on, delta %.3f cm, max skip %d (patterns with an OT layer: %d)\n", own_delta,
           own_skip, own_skip_ot >= 0 ? own_skip_ot : own_skip);
  if (dedup_n > 0)
    printf("[seedsurf] cleaning over all patterns: a quad sharing >= %d hits with a better one is dropped\n", dedup_n);
  printf("[seedsurf] pt_min %.3f GeV, d0_max %.3f cm, zv %.1f cm, marg_b %.4f; c: phi %.4f q %.4f; d: phi %.4f q %.4f%s\n",
         P.pt_min, P.d0_max, P.zv, P.marg_b, P.phi_c, P.q_c, P.phi_d, P.q_d, P.no_cut ? "; NO CUTS" : "");
  auto params_of = [&](const Pattern &p) {
    SurfParams Q = P;
    if (p.win[0] >= 0)
      Q.phi_c = p.win[0], Q.q_c = p.win[1], Q.phi_d = p.win[2], Q.q_d = p.win[3];
    Q.b_phi_d = p.bwin[0], Q.b_q_d = p.bwin[1];
    Q.s_ref = p.sref;
    for (auto w : p.eta_win) {
      w.aphi *= win_scale, w.bphi *= win_scale, w.aq *= win_scale, w.bq *= win_scale;
      Q.add_eta_win(w);
    }
    Q.phi_c *= win_scale, Q.q_c *= win_scale, Q.phi_d *= win_scale, Q.q_d *= win_scale;
    Q.b_phi_d *= win_scale, Q.b_q_d *= win_scale;
    return Q;
  };
  if (chain_holes >= 0)
    for (const auto &p : pats) {
      const SurfParams Q = params_of(p);
      for (SurfChain &C : CH) {
        auto &wc = C.win_c[{p.l[0], p.l[1], p.l[2]}];
        wc.first = std::max(wc.first, (float)Q.phi_c), wc.second = std::max(wc.second, (float)Q.q_c);
        C.win_d[{p.l[0], p.l[1], p.l[2], p.l[3]}] = {(float)Q.phi_d, (float)Q.b_phi_d, (float)Q.q_d, (float)Q.b_q_d, Q.s_ref,
                                                     std::vector<SurfParams::EtaWin>(Q.eta_win, Q.eta_win + Q.n_eta_win),
                                                     (float)Q.q_c};
      }
    }
  for (const auto &p : pats) {
    const SurfParams Q = params_of(p);
    printf("[seedsurf] pattern %s (%d %d %d %d); c: phi %.4f q %.4f; d: phi %.4f + %.4f/pT q %.4f + %.4f/pT%s\n", p.name.c_str(),
           p.l[0], p.l[1], p.l[2], p.l[3], Q.phi_c, Q.q_c, Q.phi_d, Q.b_phi_d, Q.q_d, Q.b_q_d,
           Q.s_ref > 0 ? (" (b x s_cd / " + std::to_string(Q.s_ref) + " cm)").c_str() : "");
  }

  FILE *fr = nullptr;
  if (!resid_out.empty()) {
    fr = fopen(resid_out.c_str(), "w");
    fprintf(fr, "# seedsurf --resid; pt_min %.3f d0_max %.3f\n", P.pt_min, P.d0_max);
    fprintf(fr, "# R pat ev label eta pt ncomb  z0 dphi_b w_b  dphi_c wd0_c dq_c  dphi_d dq_d pte  hit-sigmas: sphi_d sq_d"
                "  skips: lead inner  phi: d-hit track  r_d  s_cd\n");
  }
  // --resid: the chain's crossing test on each combination's a-d line (margin 0.2 cm, as the chain uses)
  SurfOwnership RO;
  if (fr) {
    std::set<int> have;
    for (const auto &p : pats)
      for (int l : p.l)
        have.insert(l);
    RO.setup(ti, have);
    RO.delta = 0.2;
  }
  FILE *fd = dump_out.empty() ? nullptr : fopen(dump_out.c_str(), "w");
  FILE *fq = qdump_out.empty() ? nullptr : fopen(qdump_out.c_str(), "w");
  if (fq)
    fprintf(fq, "# seedsurf --quad-dump: per event 'E ev n_findable'; per quad kept after the cleaning\n"
                "# Q ip l0 l1 l2 l3  tru fake findable lab  eta pte  d0_abc d0_acd  rc_q rc_phi rd_phi rd_q"
                "  rows0 cols0 rows1 cols1 rows2 cols2 rows3 cols3  ot1 dphi dz n3 same  otz maxc nun\n"
                "# rc_*, rd_*: the c and d residuals over their windows (signed); ot1: -1 not predicted, 0 no hit in the"
                " fetch, 1 best OT1-P hit given, 2 the quad's d is OT1-P; dphi, dz: its residuals to the b-c-d helix;"
                " n3: hits within score 9; same: its label is the quad's; otz: the helix's z at OT1 (mid of the two edges);"
                " maxc: most hits from one sim track, nun: unlabelled hits. The checked layer is OT2-P (6) for a quad whose d is"
                " OT1-P, else OT1-P (4)\n");

  // --margins: the reference list, per event: (pattern index, hits)
  std::map<int, std::set<std::pair<int, Quad>>> mref;
  struct MarginRow {
    int iev, ip, side;  // side: +1 only in this run, -1 only in the reference
    Quad q;
    double rel;         // the smallest |margin| / window over the cuts; NaN: a stage not evaluable
    const char *cut;
  };
  std::vector<MarginRow> mrows;
  long m_same = 0;
  if (!margins_ref.empty()) {
    if (chain_holes < 0) {
      printf("--margins needs --chain\n");
      return 1;
    }
    FILE *fm = fopen(margins_ref.c_str(), "r");
    if (!fm) {
      printf("cannot open %s\n", margins_ref.c_str());
      return 1;
    }
    int e, ip;
    Quad q;
    while (fscanf(fm, "%d %d %u %u %u %u", &e, &ip, &q[0], &q[1], &q[2], &q[3]) == 6)
      mref[e].insert({ip, q});
    fclose(fm);
  }

  SurfChainBatch *const CB[2] = {&seeder.finder(0), &seeder.finder(1)};
  if (chain_holes >= 0 && chain_batch) {
    seeder.setup(CH[0], CH[1]);
    for (int sd = 0; sd < 2; ++sd) {
      SurfChainBatch &B = *CB[sd];
      B.d_mode = chain_batch_d;
      if (g_surf_fast_check.on)
        B.d_check = [](const SeedCand &c, const std::vector<const SurfLayer *> &lay, bool disc, float u, bool ok,
                       float px, float py, float qp, float wq, float wp) {
          surf_check_d(c, lay, disc, u, ok, px, py, qp, wq, wp, g_surf_fast_check);
        };
      B.fk_score = fk_score, B.fk_shape = fk_shape, B.fk_ot2 = fk_ot2;
      B.ot2_aphi = ot2_win[0], B.ot2_bphi = ot2_win[1], B.ot2_aq = ot2_win[2], B.ot2_bq = ot2_win[3];
      B.ot2_phimin = ot2_phimin;
      for (int l = 0; l < 4; ++l)
        B.shape_[l] = shape_tab[l];
    }
  }
  if ((fk_score > 0 || fk_shape || fk_ot2 > 0) && !chain_batch) {
    printf("--fk-* need --chain-batch\n");
    return 1;
  }
  if (fk_shape)
    for (int l = 0; l < 4; ++l)
      if (!shape_tab[l].ok())
        printf("[seedsurf] --fk-shape: no --shape-win for layer %d, no shape cut there\n", l);

  // --quad-dump: the next outer layer a quad's helix is checked against, OT1-P (4) or OT2-P (6);
  // --fk-ot2 and --attach-ot1 need them too
  if (!qdump_out.empty() || fk_ot2 > 0 || attach_ot1 > 0)
    for (int l : {4, 6})
      layers.add_layer(l, ti[l], 2.0);
  // the batch finder is float: no double r, phi cache per hit
  if (chain_holes >= 0 && chain_batch)
    layers.set_with_double(false);

  std::vector<Stats> S(pats.size());
  Stats SU;  // the union
  Stats SC;  // the chain's own counters and time
  // --attach-ot1: kept quads with a pixel d; the helix reaches OT1-P; an OT1-P hit attached; of the true
  // quads that reach it: attached, and the attached hit carries the quad's label
  long at_n = 0, at_reach = 0, at_att = 0, at_true = 0, at_true_att = 0, at_true_same = 0;
  double t_attach = 0;
  double t_fill = 0, t_clean = 0;
  long n_unbound = 0;  // labels dropped by --bind
  int npat = pats.size();

  for (int iev = 0; iev < n_events; ++iev) {
    Event ev(iev, ti.n_layers());
    ev.read_in(df);
    const auto t0 = clk::now();
    seeder.fill(ev.layerHits_);
    t_fill += secs(t0, clk::now());

    const MCHitInfoVec &mc = ev.simHitsInfo_;
    const bool bind = bind_cm > 0 && ev.simHitStates_.size() == mc.size();
    if (bind_cm > 0 && !bind && iev == 0)
      printf("[seedsurf] --bind: no SimHitStates in this sample, labels NOT bound\n");
    auto label_l = [&](const Hit &h, int layer) -> int {
      const int mid = h.mcHitID();
      if (mid < 0 || mid >= (int)mc.size())
        return -1;
      if (bind) {
        const SimHitState &sh = ev.simHitStates_[mid];
        if (!sh.is_valid())
          return -1;
        const double dx = sh.pos[0] - h.x(), dy = sh.pos[1] - h.y(), dz = sh.pos[2] - h.z();
        const double tol = SurfOwnership::is_pix(layer) ? bind_cm : 2.6;
        if (dx * dx + dy * dy + dz * dz > tol * tol) {
          ++n_unbound;
          return -1;
        }
      }
      return mc[mid].mcTrackID();
    };
    auto keep = [&](int l) {
      if (l < 0 || l >= (int)ev.simTracks_.size())
        return false;
      const Track &t = ev.simTracks_[l];
      return t.pT() >= P.pt_min && std::hypot(t.x() - ev.beamSpot_.x, t.y() - ev.beamSpot_.y) <= P.d0_max &&
             std::abs(t.momEta()) < eta_max;
    };
    // labelled hits per layer, per track
    std::map<int, std::unordered_map<int, std::vector<int>>> lab_hits;
    for (auto &kv : layers.layer_map()) {
      const HitVec &hv = ev.layerHits_[kv.first];
      auto &m = lab_hits[kv.first];
      for (int k = 0; k < (int)hv.size(); ++k) {
        const int l = label_l(hv[k], kv.first);
        if (keep(l))
          m[l].push_back(k);
      }
    }
    std::unordered_set<int> fb_any;
    // every quad of every pattern, for the union after the cleaning
    struct Cand {
      int ip, lab, b, pb;
      Quad q;
      double score;
      bool tru, fake;
    };
    std::vector<Cand> cands;
    std::vector<std::vector<Quad>> chain_q;
    std::vector<std::vector<float>> chain_sc;  // the batch finder's cleaning score, parallel to chain_q
    if (chain_holes >= 0) {
      std::vector<std::pair<std::array<int, 4>, Quad>> cq;
      std::vector<float> csc;
      SeedCounters cnt;
      const auto s0 = clk::now();
      const auto &LM = layers.layer_map();
      if (chain_batch)
        seeder.find(cq, cnt, &csc);
      else {
        CH[0].run(LM, cq, cnt);
        CH[1].run(LM, cq, cnt);
      }
      SC.t_find += secs(s0, clk::now());
      SC.quads += cnt.quads, SC.doublets += cnt.doublets, SC.triplets += cnt.triplets;
      SC.c_touched += cnt.c_touched, SC.d_touched += cnt.d_touched;
      // each quad to the pattern of its layers; a combination no pattern lists becomes a dynamic one
      for (size_t ie = 0; ie < cq.size(); ++ie) {
        const auto &e = cq[ie];
        int ip = -1;
        for (int i = 0; i < (int)pats.size() && ip < 0; ++i)
          if (pats[i].l == e.first)
            ip = i;
        if (ip < 0) {
          Pattern np;
          np.l = e.first;
          np.dynamic = true;
          char b[8];
          for (int k = 0; k < 4; ++k)
            np.name += std::string(k ? " " : "") + lname(np.l[k], b);
          np.name += " *";
          pats.push_back(np);
          S.emplace_back();
          ip = pats.size() - 1;
        }
        if ((int)chain_q.size() <= ip)
          chain_q.resize(pats.size()), chain_sc.resize(pats.size());
        chain_q[ip].push_back(e.second);
        if (!csc.empty())
          chain_sc[ip].push_back(csc[ie]);
      }
      chain_q.resize(pats.size());
      chain_sc.resize(pats.size());
      npat = pats.size();

      if (!margins_ref.empty()) {
        std::set<std::pair<int, Quad>> cur;
        for (int ip = 0; ip < (int)chain_q.size(); ++ip)
          for (const auto &q : chain_q[ip])
            cur.insert({ip, q});
        const auto &ref = mref[iev];
        // every cut of one quad in double, with the windows the chain uses for its layers
        auto evaluate = [&](int ip, const Quad &q, int side) {
          const auto &l = pats[ip].l;
          const SurfLayer *Ls[4];
          surf::P3 h[4];
          for (int k = 0; k < 4; ++k) {
            Ls[k] = layers.layer(l[k]);
            const Hit &hh = ev.layerHits_[l[k]][q[k]];
            h[k] = {hh.x(), hh.y(), hh.z()};
          }
          SurfParams Qm = CH[0].P;
          if (auto it = CH[0].win_c.find({l[0], l[1], l[2]}); it != CH[0].win_c.end())
            Qm.phi_c = it->second.first, Qm.q_c = it->second.second;
          if (auto it = CH[0].win_d.find({l[0], l[1], l[2], l[3]}); it != CH[0].win_d.end())
            Qm.phi_d = it->second.aphi, Qm.b_phi_d = it->second.bphi, Qm.q_d = it->second.aq, Qm.b_q_d = it->second.bq,
            Qm.s_ref = it->second.sref;
          SurfEval e;
          surf_eval(Qm, Ls, h, e);
          const double ra = h[0].r(), rb = h[1].r();
          const double wb = e.w_b + Qm.marg_b, wc = e.wd0_c + Qm.phi_c;
          const double wpd = Qm.wphi_d(e.pte, e.s_cd), wqd = Qm.wq_d(e.pte, e.s_cd);
          const std::pair<const char *, double> m[] = {
              {"b_r", (rb - ra - 0.1) / 1.0},
              {"b_phi", (wb - std::abs(e.dphi_b)) / wb},
              {"b_z0", std::min(e.z0 - (Qm.bs_z - Qm.zv), Qm.bs_z + Qm.zv - e.z0) / Qm.zv},
              {"c_phi", (wc - std::abs(e.dphi_c)) / wc},
              {"c_q", (Qm.q_c - std::abs(e.dq_c)) / Qm.q_c},
              {"d_phi", (wpd - std::abs(e.dphi_d)) / wpd},
              {"d_q", (wqd - std::abs(e.dq_d)) / wqd}};
          // the crossing tests: the chain routes on the a-b line (holes before b, the layer after b) and
          // on the a-c line (the layer after c). The smallest distance of any envelope crossing of either
          // line to a state threshold, over the crossing margin.
          double xm = 1e30;
          for (int ln = 0; ln < 2; ++ln) {
            const SurfChain::LineRZ L(h[0], h[ln ? 2 : 1]);
            for (const auto &en : OWN.env) {
              double x1, x2;
              if (!en.disc)
                x1 = L.z0 + L.cot * en.pos_lo, x2 = L.z0 + L.cot * en.pos_hi;
              else {
                if (std::abs(L.cot) < 1e-9 || (en.pos - L.z0) / L.cot <= 0)
                  continue;
                x1 = (en.pos_lo - L.z0) / L.cot, x2 = (en.pos_hi - L.z0) / L.cot;
              }
              const double xl = std::min(x1, x2), xh = std::max(x1, x2), dl = OWN.delta;
              for (double d : {xl - (en.lo + dl), xh - (en.hi - dl), xh - (en.lo - dl), xl - (en.hi + dl)})
                xm = std::min(xm, std::abs(d));
            }
          }
          MarginRow r{iev, ip, side, q, 1e30, "structural"};
          if (xm / OWN.delta < 1e-2)
            r.rel = xm / OWN.delta, r.cut = "crossing";
          for (const auto &c : m) {
            if (!std::isfinite(c.second)) {
              r.rel = NAN, r.cut = "not evaluable";
              break;
            }
            if (std::abs(c.second) < std::abs(r.rel))
              r.rel = c.second, r.cut = c.first;
          }
          mrows.push_back(r);
        };
        for (const auto &x : ref)
          if (!cur.count(x))
            evaluate(x.first, x.second, -1);
          else
            ++m_same;
        for (const auto &x : cur)
          if (!ref.count(x))
            evaluate(x.first, x.second, +1);
      }
    }
    for (int ip = 0; ip < npat; ++ip) {
      const Pattern &p = pats[ip];
      Stats &st = S[ip];
      const SurfLayer *L[4] = {layers.layer(p.l[0]), layers.layer(p.l[1]), layers.layer(p.l[2]), layers.layer(p.l[3])};
      std::vector<Quad> quads;
      SeedCounters cnt;
      const auto s0 = clk::now();
      const SurfParams Q = params_of(p);
      if (chain_holes >= 0) {
        quads = chain_q[ip];
        cnt.quads = quads.size();
      } else
        find_quads_surf(Q, L, quads, cnt);
      st.t_find += secs(s0, clk::now());
      if (OWN.delta >= 0 && chain_holes < 0) {
        // keep the quads this pattern owns, from the a-d line
        const HitVec *Hq[4] = {&ev.layerHits_[p.l[0]], &ev.layerHits_[p.l[1]], &ev.layerHits_[p.l[2]], &ev.layerHits_[p.l[3]]};
        std::vector<Quad> kept;
        for (const auto &q : quads) {
          const Hit &ha = (*Hq[0])[q[0]], &hd = (*Hq[3])[q[3]];
          const double cot = (hd.z() - ha.z()) / (hd.r() - ha.r()), z0 = ha.z() - cot * ha.r();
          if (OWN.owns(p.l.data(), z0, cot))
            kept.push_back(q);
          else if (own_debug > 0) {
            // a rejected TRUE quad: where the line says the track goes
            int l0 = -1;
            bool tru = true;
            for (int k = 0; k < 4; ++k) {
              const int lab = label_l((*Hq[k])[q[k]], p.l[k]);
              tru = tru && lab >= 0 && (k == 0 || lab == l0);
              if (k == 0)
                l0 = lab;
            }
            if (tru && keep(l0)) {
              --own_debug;
              printf("[own] %s eta %.3f z0 %.2f cot %.3f (sim z %.2f) a(r %.2f z %.2f) d(r %.2f z %.2f):", p.name.c_str(),
                     ev.simTracks_[l0].momEta(), z0, cot, ev.simTracks_[l0].z(), ha.r(), ha.z(), hd.r(), hd.z());
              OWN.print(z0, cot);
            }
          }
        }
        st.own_rejected += quads.size() - kept.size();
        quads.swap(kept);
      }
      st.quads += cnt.quads, st.doublets += cnt.doublets, st.triplets += cnt.triplets;
      st.c_touched += cnt.c_touched, st.d_touched += cnt.d_touched;

      // findable for this pattern
      std::unordered_map<int, int> n_true;
      for (const auto &kv : lab_hits[p.l[0]]) {
        bool all = true;
        for (int k = 1; k < 4; ++k)
          all = all && lab_hits[p.l[k]].count(kv.first);
        if (all)
          n_true[kv.first] = 0;
      }
      const HitVec *H[4] = {&ev.layerHits_[p.l[0]], &ev.layerHits_[p.l[1]], &ev.layerHits_[p.l[2]], &ev.layerHits_[p.l[3]]};
      for (const auto &q : quads) {
        int ls[4];
        for (int k = 0; k < 4; ++k)
          ls[k] = label_l((*H[k])[q[k]], p.l[k]);
        int nun = 0, maxc = 0;
        for (int u = 0; u < 4; ++u) {
          nun += ls[u] < 0;
          if (ls[u] < 0)
            continue;
          int c = 0;
          for (int w = 0; w < 4; ++w)
            c += ls[w] == ls[u];
          maxc = std::max(maxc, c);
        }
        const Hit &ha = (*H[0])[q[0]], &hd = (*H[3])[q[3]];
        const int b = ebin(std::abs(std::asinh((hd.z() - ha.z()) / (hd.r() - ha.r()))));
        const bool tru = maxc == 4, fake = maxc + nun < 4;
        double pte;
        {
          surf::P3 h3[3];
          for (int k = 0; k < 3; ++k) {
            const Hit &h = (*H[k])[q[k]];
            const float r = h.r(), ph = h.phi();
            h3[k] = {r * std::cos(ph), r * std::sin(ph), h.z()};
          }
          pte = surf::Helix(h3[0], h3[1], h3[2]).pt();
        }
        const int pb = pbin(pte);
        st.pe_all[pb] += 1;
        st.pe_true[pb] += tru;
        st.pe_fake[pb] += fake;
        if (b >= 0) {
          st.q_all[b] += 1;
          st.q_true[b] += tru;
          st.q_fake[b] += fake;
          st.q_undec[b] += !tru && !fake;
        }
        if (tru) {
          if (auto it = n_true.find(ls[0]); it != n_true.end())
            ++it->second;
        }
        if (chain_holes >= 0 && chain_sc[ip].size() == quads.size()) {
          // the batch finder's own score, in float
          cands.push_back({ip, tru ? ls[0] : -1, b, pb, q, chain_sc[ip][&q - quads.data()], tru, fake});
        } else {
          // the quad's quality, for the cleaning: its own c and d residuals over the windows, in double
          const SurfLayer *Ls[4] = {L[0], L[1], L[2], L[3]};
          surf::P3 h4[4];
          for (int k = 0; k < 4; ++k) {
            const Hit &h = (*H[k])[q[k]];
            const float r = h.r(), ph = h.phi();
            h4[k] = {r * std::cos(ph), r * std::sin(ph), h.z()};
          }
          SurfEval e;
          surf_eval(Q, Ls, h4, e);
          const double wpd = Q.wphi_d(e.pte), wqd = Q.wq_d(e.pte);
          const double sc = (e.dq_c / Q.q_c) * (e.dq_c / Q.q_c) + (e.dphi_d / wpd) * (e.dphi_d / wpd) + (e.dq_d / wqd) * (e.dq_d / wqd);
          cands.push_back({ip, tru ? ls[0] : -1, b, pb, q, std::isfinite(sc) ? sc : 1e30, tru, fake});
        }
        if (fd)
          fprintf(fd, "%d %d %u %u %u %u\n", iev, ip, q[0], q[1], q[2], q[3]);
      }
      for (const auto &kv : n_true) {
        if (!p.dynamic)
          fb_any.insert(kv.first);
        const int b = ebin(std::abs(ev.simTracks_[kv.first].momEta()));
        if (b < 0)
          continue;
        st.den[b] += 1;
        if (kv.second > 0) {
          st.found[b] += 1;
          st.dup[b] += kv.second - 1;
        }
      }
      if (fr) {
        const SurfLayer *Ls[4] = {L[0], L[1], L[2], L[3]};
        for (const auto &kv : n_true) {
          const auto &v0 = lab_hits[p.l[0]][kv.first], &v1 = lab_hits[p.l[1]][kv.first],
                     &v2 = lab_hits[p.l[2]][kv.first], &v3 = lab_hits[p.l[3]][kv.first];
          const int ncomb = v0.size() * v1.size() * v2.size() * v3.size();
          if (ncomb > 64)
            continue;
          const Track &t = ev.simTracks_[kv.first];
          for (int i0 : v0)
            for (int i1 : v1)
              for (int i2 : v2)
                for (int i3 : v3) {
                  const int ii[4] = {i0, i1, i2, i3};
                  surf::P3 h[4];
                  for (int k = 0; k < 4; ++k) {
                    const Hit &hh = (*H[k])[ii[k]];
                    const float r = hh.r(), ph = hh.phi();
                    h[k] = {r * std::cos(ph), r * std::sin(ph), hh.z()};
                  }
                  SurfEval e;
                  surf_eval(P, Ls, h, e);
                  // the d hit's own sigmas in the target's (phi, q)
                  const Hit &hd = (*H[3])[i3];
                  const double x = hd.x(), y = hd.y(), r2 = x * x + y * y;
                  const double sphi = std::sqrt(std::max(0.0, (double)(y * y * hd.exx() - 2 * x * y * hd.exy() + x * x * hd.eyy()))) / r2;
                  const double sq = L[3]->disc ? std::sqrt(std::max(0.0, (double)(x * x * hd.exx() + 2 * x * y * hd.exy() + y * y * hd.eyy()) / r2))
                                               : std::sqrt((double)hd.ezz());
                  const double cot_ad = (h[3].z - h[0].z) / (std::hypot(h[3].x, h[3].y) - std::hypot(h[0].x, h[0].y));
                  const double z0_ad = h[0].z - cot_ad * std::hypot(h[0].x, h[0].y);
                  int lead, inner;
                  RO.skips(p.l.data(), z0_ad, cot_ad, lead, inner);
                  fprintf(fr, "R %d %d %d %.4f %.3f %d  %.4f %.6g %.6g  %.6g %.6g %.6g  %.6g %.6g %.4g  %.4g %.4g  %d %d  %.4f %.4f  %.3f %.3f\n", ip, iev,
                          kv.first, t.momEta(), t.pT(), ncomb, e.z0 - t.z(), e.dphi_b, e.w_b, e.dphi_c, e.wd0_c, e.dq_c,
                          e.dphi_d, e.dq_d, e.pte, sphi, sq, lead, inner, std::atan2(h[3].y, h[3].x), t.momPhi(),
                          std::hypot(h[3].x, h[3].y), e.s_cd);
                }
        }
      }
    }
    // the union over patterns, after the cleaning
    std::vector<char> keep_c(cands.size(), 1);
    const auto tc0 = clk::now();
    if (dedup_n > 0) {
      // the quads in the order of cands, which breaks the cleaning's ties
      std::vector<std::array<int, 4>> cl_l(cands.size());
      std::vector<Quad> cl_q(cands.size());
      std::vector<float> cl_s(cands.size());
      for (size_t i = 0; i < cands.size(); ++i)
        cl_l[i] = pats[cands[i].ip].l, cl_q[i] = cands[i].q, cl_s[i] = (float)cands[i].score;
      seeder.clean(ev.layerHits_, cl_l, cl_q, cl_s, dedup_n, keep_c);
    }
    t_clean += secs(tc0, clk::now());
    if (attach_ot1 > 0) {
      const SurfLayer &LP = *layers.layer(4);
      const auto a0 = clk::now();
      std::vector<int> att(cands.size(), -1);  // the attached OT1-P hit (index into layerHits_[4]), -1 none
      auto fxyz = [&](int l, unsigned int k, float &x, float &y, float &z) {
        const Hit &h = ev.layerHits_[l][k];
        const float r = h.r(), ph = h.phi();
        x = r * std::cos(ph), y = r * std::sin(ph), z = h.z();
      };
      for (int i = 0; i < (int)cands.size(); ++i) {
        if (!keep_c[i])
          continue;
        const Cand &c = cands[i];
        const auto &ll = pats[c.ip].l;
        if (!SurfOwnership::is_pix(ll[3]))
          continue;
        ++at_n;
        float x[4], y[4], z[4];
        for (int k = 0; k < 4; ++k)
          fxyz(ll[k], c.q[k], x[k], y[k], z[k]);
        surfb::HelixF abc, bcd;
        abc.make(x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2]);
        bcd.make(x[1], y[1], z[1], x[2], y[2], z[2], x[3], y[3], z[3]);
        float zm, dp, dz, sc;
        const float pte = abc.ok ? abc.pt() : 1.0f, ip = 1.0f / std::max(0.9f, pte);
        const float wp = attach_ot1 * (ot1_win[0] + ot1_win[1] * ip), wz = attach_ot1 * (ot1_win[2] + ot1_win[3] * ip);
        const int kn = SurfChainBatch::next_hit(LP, bcd, pte, wp, wz, zm, dp, dz, sc);
        // in OT1's acceptance, as for OT2-P in the finder
        if (kn == -2 || std::abs(zm) >= LP.q_hi - 2) {
          att[i] = -2;
          continue;
        }
        ++at_reach;
        if (kn >= 0 && std::abs(dp) < wp && std::abs(dz) < wz)
          att[i] = LP.orig_[kn], ++at_att;
      }
      t_attach += secs(a0, clk::now());
      for (int i = 0; i < (int)cands.size(); ++i)
        if (keep_c[i] && cands[i].tru && SurfOwnership::is_pix(pats[cands[i].ip].l[3]) && att[i] != -2) {
          ++at_true;
          if (att[i] >= 0) {
            ++at_true_att;
            at_true_same += label_l(ev.layerHits_[4][att[i]], 4) == cands[i].lab;
          }
        }
    }
    Stats &su = SU;
    std::unordered_map<int, int> fd_any;
    std::unordered_map<int, std::set<int>> fd_pat;
    for (int i = 0; i < (int)cands.size(); ++i) {
      if (!keep_c[i])
        continue;
      const Cand &c = cands[i];
      su.quads += 1;
      su.pe_all[c.pb] += 1, su.pe_true[c.pb] += c.tru, su.pe_fake[c.pb] += c.fake;
      if (c.b >= 0) {
        su.q_all[c.b] += 1, su.q_true[c.b] += c.tru, su.q_fake[c.b] += c.fake, su.q_undec[c.b] += !c.tru && !c.fake;
        const double ae = (c.b + 0.5) * kEtaW;
        const int rg = ae < 0.8 ? 0 : ae < 1.6 ? 1 : ae < 2.4 ? 2 : 3;
        su.qp_all[rg][c.pb] += 1, su.qp_true[rg][c.pb] += c.tru, su.qp_fake[rg][c.pb] += c.fake;
      }
      if (c.tru && fb_any.count(c.lab)) {
        fd_any[c.lab] += 1;
        fd_pat[c.lab].insert(c.ip);
      }
    }
    for (int l : fb_any) {
      const int b = ebin(std::abs(ev.simTracks_[l].momEta()));
      if (b < 0)
        continue;
      su.den[b] += 1;
      const double ae = std::abs(ev.simTracks_[l].momEta());
      const int rg = ae < 0.8 ? 0 : ae < 1.6 ? 1 : ae < 2.4 ? 2 : 3, pb = pbin(ev.simTracks_[l].pT());
      su.sp_den[rg][pb] += 1;
      if (auto it = fd_any.find(l); it != fd_any.end()) {
        su.found[b] += 1;
        su.dup[b] += it->second - 1;
        su.multi[b] += fd_pat[l].size() > 1;
        su.sp_found[rg][pb] += 1, su.sp_dup[rg][pb] += it->second - 1;
      }
    }
    // --quad-dump: one row per quad kept after the cleaning, with what a fake-rejection cut could use
    if (fq) {
      fprintf(fq, "E %d %zu\n", iev, fb_any.size());
      // the next outer P layer: OT1-P after a pixel d, OT2-P after an OT1-P d
      const SurfLayer *LP4 = layers.layer(4), *LP6 = layers.layer(6);
      const double bx = ev.beamSpot_.x, by = ev.beamSpot_.y;
      auto p3of = [&](int l, unsigned int k) {
        const Hit &h = ev.layerHits_[l][k];
        const float r = h.r(), ph = h.phi();
        return surf::P3(r * std::cos(ph), r * std::sin(ph), h.z());
      };
      // transverse distance of the helix circle from the beam spot
      auto d0of = [&](const surf::Helix &hx) -> double {
        if (!hx.ok)
          return NAN;
        if (std::abs(hx.k) < 1e-9)
          return std::abs(hx.tx * (hx.c.y - by) - hx.ty * (hx.c.x - bx));
        const double Cx = hx.c.x - hx.ty / hx.k, Cy = hx.c.y + hx.tx / hx.k;
        return std::abs(std::hypot(Cx - bx, Cy - by) - 1 / std::abs(hx.k));
      };
      for (int i = 0; i < (int)cands.size(); ++i) {
        if (!keep_c[i])
          continue;
        const Cand &c = cands[i];
        const auto &ll = pats[c.ip].l;
        surf::P3 h[4];
        for (int k = 0; k < 4; ++k)
          h[k] = p3of(ll[k], c.q[k]);
        const surf::Helix abc(h[0], h[1], h[2]), acd(h[0], h[2], h[3]), bcd(h[1], h[2], h[3]);
        const double eta = std::asinh((h[3].z - h[0].z) / (h[3].r() - h[0].r()));
        const SurfParams Q = params_of(pats[c.ip]);
        const SurfLayer *Ls[4] = {layers.layer(ll[0]), layers.layer(ll[1]), layers.layer(ll[2]), layers.layer(ll[3])};
        SurfEval e;
        surf_eval(Q, Ls, h, e);
        const double wc = e.wd0_c + Q.phi_c, wpd = Q.wphi_d(e.pte), wqd = Q.wq_d(e.pte);
        int ot = -1, n3 = 0, same = 0;
        double bdp = NAN, bdz = NAN, otz = NAN;
        const int elay = ll[3] == 4 ? 6 : 4;
        const SurfLayer *LP = elay == 6 ? LP6 : LP4;
        const HitVec &hp4 = ev.layerHits_[elay];
        if (!SurfOwnership::is_pix(ll[3]) && ll[3] != 4)
          ot = 2;
        else if (LP && bcd.ok) {
          double q0, p0, q1, p1;
          const bool ok0 = bcd.predict(false, LP->qbar_lo, q0, p0), ok1 = bcd.predict(false, LP->qbar_hi, q1, p1);
          if (ok0 || ok1) {
            if (!ok0)
              q0 = q1, p0 = p1;
            if (!ok1)
              q1 = q0, p1 = p0;
            ot = 0;
            otz = 0.5 * (q0 + q1);
            const double pte = std::max(0.5, e.pte), sphi = 0.0005 + 3.2e-3 / pte, sz = 0.075 + 0.0316 / pte;
            const double dm = p0 + 0.5 * surf::wrap(p1 - p0), dh = 0.5 * std::abs(surf::wrap(p1 - p0));
            double best = 1e30;
            LP->for_each_in(surf::phi_bins(*LP, dm, dh + 10 * sphi),
                               surf::q_bins(*LP, std::min(q0, q1) - 10 * sz, std::max(q0, q1) + 10 * sz),
                               [&](unsigned int kk) {
                                 double qp, pp;
                                 if (!bcd.predict(false, LP->qbar(kk), qp, pp))
                                   return;
                                 const double dp = surf::wrap(LP->phi_[kk] - pp), dz = LP->z_[kk] - qp;
                                 const double sc = (dp / sphi) * (dp / sphi) + (dz / sz) * (dz / sz);
                                 n3 += sc < 9;
                                 if (sc < best) {
                                   best = sc, bdp = dp, bdz = dz, ot = 1;
                                   same = c.lab >= 0 && label_l(hp4[LP->orig_[kk]], elay) == c.lab;
                                 }
                               });
          }
        }
        // label composition: the largest number of hits from one sim track, and unlabelled hits
        int ls[4], maxc = 0, nun = 0;
        for (int k = 0; k < 4; ++k)
          ls[k] = label_l(ev.layerHits_[ll[k]][c.q[k]], ll[k]);
        for (int u = 0; u < 4; ++u) {
          nun += ls[u] < 0;
          int cc = 0;
          for (int w = 0; w < 4; ++w)
            cc += ls[u] >= 0 && ls[w] == ls[u];
          maxc = std::max(maxc, cc);
        }
        int sp[8];
        for (int k = 0; k < 4; ++k) {
          const Hit &hh = ev.layerHits_[ll[k]][c.q[k]];
          sp[2 * k] = hh.spanRows(), sp[2 * k + 1] = hh.spanCols();
        }
        fprintf(fq, "Q %d %d %d %d %d  %d %d %d %d  %.4f %.4f  %.5f %.5f  %.4f %.4f %.4f %.4f  %d %d %d %d %d %d %d %d  %d %.6f %.5f %d %d  %.3f %d %d\n",
                c.ip, ll[0], ll[1], ll[2], ll[3], c.tru, c.fake, c.tru && fb_any.count(c.lab), c.lab, eta, e.pte,
                d0of(abc), d0of(acd), e.dq_c / Q.q_c, e.dphi_c / wc, e.dphi_d / wpd, e.dq_d / wqd, sp[0], sp[1], sp[2],
                sp[3], sp[4], sp[5], sp[6], sp[7], ot, bdp, bdz, n3, same, otz, maxc, nun);
      }
    }
  }
  if (fr)
    fclose(fr);
  if (fd)
    fclose(fd);
  if (fq)
    fclose(fq);

  const double ne = n_events;
  printf("[seedsurf] %d events; layer fill %.3f ms/ev\n", n_events, 1e3 * t_fill / ne);
  if (dedup_n > 0)
    printf("[seedsurf] cleaning %.3f ms/ev\n", 1e3 * t_clean / ne);
  if (bind_cm > 0)
    printf("[seedsurf] --bind %.3f cm: %.1f label lookups per event dropped (the sim hit farther than that)\n", bind_cm,
           n_unbound / ne);
  // the chain's own counters, from whichever finder ran
  const double ch_fwd = chain_batch ? CB[0]->n_forwarded + CB[1]->n_forwarded : CH[0].n_forwarded + CH[1].n_forwarded;
  const double ch_drop = chain_batch ? CB[0]->n_dropped + CB[1]->n_dropped : CH[0].n_dropped + CH[1].n_dropped;
  const double ch_tst = chain_batch ? CB[0]->t_start + CB[1]->t_start : CH[0].t_start + CH[1].t_start;
  const double ch_cc = chain_batch ? CB[0]->cyc_c + CB[1]->cyc_c : CH[0].cyc_c + CH[1].cyc_c;
  const double ch_cd = chain_batch ? CB[0]->cyc_d + CB[1]->cyc_d : CH[0].cyc_d + CH[1].cyc_d;
  const double ch_nc = chain_batch ? CB[0]->n_c_cand + CB[1]->n_c_cand : CH[0].n_c_cand + CH[1].n_c_cand;
  const double ch_nd = chain_batch ? CB[0]->n_d_cand + CB[1]->n_d_cand : CH[0].n_d_cand + CH[1].n_d_cand;
  if (chain_holes >= 0)
    printf("[seedsurf] CHAIN%s: %.0f doublets, %.0f triplets, %.1f quads per event, %.3f ms/ev; forwarded %.1f, dropped %.1f per event\n",
           chain_batch ? " (batch)" : "", SC.doublets / ne, SC.triplets / ne, SC.quads / ne, 1e3 * SC.t_find / ne,
           ch_fwd / ne, ch_drop / ne);
  if (chain_batch && (fk_score > 0 || fk_shape || fk_ot2 > 0))
    printf("[seedsurf] FAKE CUTS (batch): score < %g, shape %s, OT2-P x%g; per event: %.1f doublets cut on shape,"
           " %.1f OT1-P-d quads tested on OT2-P of which %.1f cut\n", fk_score, fk_shape ? "on" : "off", fk_ot2,
           (CB[0]->n_fk_shape + CB[1]->n_fk_shape) / (double)n_events, (CB[0]->n_ot2_tested + CB[1]->n_ot2_tested) / (double)n_events,
           (CB[0]->n_fk_ot2 + CB[1]->n_fk_ot2) / (double)n_events);
  if (attach_ot1 > 0)
    printf("[seedsurf] ATTACH OT1-P x%g: %.1f kept quads with a pixel d per event, %.1f in OT1-P acceptance, %.1f get a hit;"
           " of %.1f true ones that reach it %.3f get one, and it is the track's own in %.3f; %.3f ms/ev\n",
           attach_ot1, at_n / (double)n_events, at_reach / (double)n_events, at_att / (double)n_events, at_true / (double)n_events,
           at_true_att / std::max(1.0, (double)at_true), at_true_same / std::max(1.0, (double)at_true_att),
           1e3 * t_attach / n_events);
  if (chain_holes >= 0 && chain_phases) {
    // the phases: start doublets by wall clock; stage c and d by cycles, scaled to the rest
    const double ts = ch_tst / ne, tr = SC.t_find / ne - ts;
    const double cc = ch_cc, cd = ch_cd;
    printf("[seedsurf] CHAIN phases per event: start doublets %.1f ms; forward pass %.1f ms, of it stage c %.0f %% "
           "(%.0f candidates) and stage d %.0f %% (%.0f candidates)\n",
           1e3 * ts, 1e3 * tr, 100 * cc / std::max(1.0, cc + cd) * 1.0, ch_nc / ne, 100 * cd / std::max(1.0, cc + cd),
           ch_nd / ne);
    printf("[seedsurf] CHAIN hits touched per event: stage c %.0f, stage d %.0f\n", SC.c_touched / ne, SC.d_touched / ne);
  }
  if (chain_holes >= 0 && g_surf_fast_check.on) {
    const SurfFastCheck &CK = g_surf_fast_check;
    if (chain_batch)
      printf("[seedsurf] FAST CHECK (batch): %ld fetched hits, float single-point prediction vs double Newton: max |dq| %.3g cm,"
             " max |dphi| %.3g rad; max error / window: q %.3g, phi %.3g; %ld hits where only one of them succeeds\n",
             CK.b_cands, CK.b_max_dq, CK.b_max_dphi, CK.b_max_rel_q, CK.b_max_rel_phi, CK.b_fail_mismatch);
    else
      printf("[seedsurf] FAST CHECK: %ld nodes, closed form vs Newton max |dq| %.3g cm, max |dphi| %.3g rad (%ld Newton"
             " failures where the closed form succeeded); %ld candidates, quadratic vs exact max error / window:"
             " q %.3g, phi %.3g\n", CK.nodes, CK.max_node_dq, CK.max_node_dphi, CK.node_fail_mismatch, CK.cands,
             CK.max_rel_q, CK.max_rel_phi);
  }
  if (!margins_ref.empty()) {
    long n_only[2] = {0, 0};
    long hist[2][6] = {};  // |rel| < 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, larger or not evaluable
    for (const auto &r : mrows) {
      const int s = r.side > 0;
      ++n_only[s];
      const double a = std::abs(r.rel);
      const int b = !std::isfinite(a) ? 5 : a < 1e-6 ? 0 : a < 1e-5 ? 1 : a < 1e-4 ? 2 : a < 1e-3 ? 3 : a < 1e-2 ? 4 : 5;
      ++hist[s][b];
    }
    printf("[seedsurf] MARGINS vs %s: %ld quads in both, %ld only in the reference, %ld only in this run\n",
           margins_ref.c_str(), m_same, n_only[0], n_only[1]);
    printf("   closest cut |margin| / window:   < 1e-6   < 1e-5   < 1e-4   < 1e-3   < 1e-2   larger/NaN\n");
    for (int s = 0; s < 2; ++s)
      printf("   %-30s %8ld %8ld %8ld %8ld %8ld %8ld\n", s ? "only in this run" : "only in the reference", hist[s][0],
             hist[s][1], hist[s][2], hist[s][3], hist[s][4], hist[s][5]);
    std::vector<MarginRow> w = mrows;
    std::sort(w.begin(), w.end(), [](const MarginRow &x, const MarginRow &y) {
      const double ax = std::isfinite(x.rel) ? std::abs(x.rel) : 1e31, ay = std::isfinite(y.rel) ? std::abs(y.rel) : 1e31;
      return ax > ay;
    });
    for (int i = 0; i < (int)w.size() && i < margins_print; ++i)
      printf("   worst %2d: ev %d %s %s hits %u %u %u %u  closest cut %s at %.3g\n", i, w[i].iev,
             pats[w[i].ip].name.c_str(), w[i].side > 0 ? "(new)" : "(ref)", w[i].q[0], w[i].q[1], w[i].q[2], w[i].q[3],
             w[i].cut, w[i].rel);
  }
  printf("   %-16s %10s %10s %10s %9s %9s  %s\n", "pattern", "doublets", "triplets", "quads", "ms/ev", "eff", "fake/dec");
  for (int ip = 0; ip < npat; ++ip) {
    const Stats &s = S[ip];
    double D = 0, F = 0, T = 0, K = 0;
    for (int b = 0; b < kNb; ++b)
      D += s.den[b], F += s.found[b], T += s.q_true[b], K += s.q_fake[b];
    printf("   %-16s %10.0f %10.0f %10.1f %9.3f %9.4f  %.4f   findable %.1f /ev, not owned %.1f /ev\n", pats[ip].name.c_str(), s.doublets / ne,
           s.triplets / ne, s.quads / ne, 1e3 * s.t_find / ne, F / std::max(1.0, D), K / std::max(1.0, T + K), D / ne, s.own_rejected / ne);
  }
  if (!truth_out.empty()) {
    FILE *f = fopen(truth_out.c_str(), "w");
    fprintf(f, "# seedsurf --truth; pt_min %.3f d0_max %.3f zv %.1f marg_b %.4f c %.4f %.4f d %.4f %.4f; %d events\n", P.pt_min,
            P.d0_max, P.zv, P.marg_b, P.phi_c, P.q_c, P.phi_d, P.q_d, n_events);
    for (int ip = 0; ip <= npat; ++ip) {
      const Stats &s = ip < npat ? S[ip] : SU;
      fprintf(f, "P %d %s\n", ip, ip < npat ? pats[ip].name.c_str() : "union");
      if (ip == npat)
        fprintf(f, "# union, after ownership and cleaning: extra_true = kept true quads beyond the first; last column = tracks found by >= 2 patterns\n");
      fprintf(f, "# E lo hi findable found extra_true  quads true fake undec   (quads by their own |eta|)\n");
      for (int b = 0; b < kNb; ++b)
        if (s.den[b] || s.q_all[b])
          fprintf(f, "E %.1f %.1f %.0f %.0f %.0f  %.0f %.0f %.0f %.0f  %.0f\n", b * kEtaW, (b + 1) * kEtaW, s.den[b], s.found[b],
                  s.dup[b], s.q_all[b], s.q_true[b], s.q_fake[b], s.q_undec[b], s.multi[b]);
      fprintf(f, "# T pte_lo pte_hi quads true fake\n");
      for (int b = 0; b < kNpe; ++b)
        if (s.pe_all[b])
          fprintf(f, "T %g %g %.0f %.0f %.0f\n", kPeEdge[b], kPeEdge[b + 1], s.pe_all[b], s.pe_true[b], s.pe_fake[b]);
      if (ip == npat) {
        fprintf(f, "# S region pt_lo pt_hi findable found extra_true   (by the sim track's pT; |eta| regions 0-0.8, 0.8-1.6, 1.6-2.4, 2.4-4)\n");
        for (int rg = 0; rg < 4; ++rg)
          for (int b = 0; b < kNpe; ++b)
            if (s.sp_den[rg][b])
              fprintf(f, "S %d %g %g %.0f %.0f %.0f\n", rg, kPeEdge[b], kPeEdge[b + 1], s.sp_den[rg][b], s.sp_found[rg][b],
                      s.sp_dup[rg][b]);
        fprintf(f, "# Q region pte_lo pte_hi quads true fake   (kept quads by their own pT_est and |eta|, same regions)\n");
        for (int rg = 0; rg < 4; ++rg)
          for (int b = 0; b < kNpe; ++b)
            if (s.qp_all[rg][b])
              fprintf(f, "Q %d %g %g %.0f %.0f %.0f\n", rg, kPeEdge[b], kPeEdge[b + 1], s.qp_all[rg][b], s.qp_true[rg][b],
                      s.qp_fake[rg][b]);
      }
    }
    fclose(f);
    printf("[seedsurf] wrote %s\n", truth_out.c_str());
  }
  return 0;
}
