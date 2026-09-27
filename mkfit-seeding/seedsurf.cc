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
//        [--own DELTA] [--own-skip K] [--own-skip-ot K]   phase-space ownership (SurfOwnership in SeedSurf.h), margin in cm
//        [--dedup N]   over all patterns, keep a quad only if it shares < N hits with every better kept quad;
//                      better = fewer outer-tracker layers in the pattern, then smaller
//                      (dq_c/q_c)^2 + (dphi_d/w_phi_d)^2 + (dq_d/w_q_d)^2
//        [--bind CM]   labels bound to geometry: needs SimHitStates in the sample
//        [--chain H] [--chain-holes-ot K] [--chain-hole-always] [--chain-any] [--chain-start-holes K] [--chain-lead-only]   the feed-forward chain (SurfChain) in
//                      place of the pattern list; the patterns then give window tables and the denominator
//        [--truth OUT.txt] [--resid OUT.txt] [--dump quads.txt] [--eta-max E]
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
//   be set from true quads.

#include "SeedSurf.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"

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
  std::string input, geom = "CMS-phase2", truth_out, resid_out, dump_out;
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
  int dedup_n = 0;         // > 0: over all patterns, drop a quad sharing >= N hits with a better kept one
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
    else if (a == "--own")
      own_delta = atof(next());
    else if (a == "--chain")
      chain_holes = atoi(next());
    else if (a == "--chain-start-holes")
      chain_start_holes = atoi(next());
    else if (a == "--chain-lead-only")
      chain_lead_only = 1;
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
  std::map<int, std::unique_ptr<SurfLayer>> layers;
  for (const auto &p : pats)
    for (int l : p.l)
      if (!layers.count(l))
        layers[l] = std::make_unique<SurfLayer>(l, ti[l], ti[l].is_barrel() ? 2.0 : 1.0);

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
      if (!layers.count(l))
        layers[l] = std::make_unique<SurfLayer>(l, ti[l], ti[l].is_barrel() ? 2.0 : 1.0);
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
        C.win_d[{p.l[0], p.l[1], p.l[2], p.l[3]}] = {(float)Q.phi_d, (float)Q.b_phi_d, (float)Q.q_d, (float)Q.b_q_d};
      }
    }
  for (const auto &p : pats) {
    const SurfParams Q = params_of(p);
    printf("[seedsurf] pattern %s (%d %d %d %d); c: phi %.4f q %.4f; d: phi %.4f + %.4f/pT q %.4f + %.4f/pT\n", p.name.c_str(),
           p.l[0], p.l[1], p.l[2], p.l[3], Q.phi_c, Q.q_c, Q.phi_d, Q.b_phi_d, Q.q_d, Q.b_q_d);
  }

  FILE *fr = nullptr;
  if (!resid_out.empty()) {
    fr = fopen(resid_out.c_str(), "w");
    fprintf(fr, "# seedsurf --resid; pt_min %.3f d0_max %.3f\n", P.pt_min, P.d0_max);
    fprintf(fr, "# R pat ev label eta pt ncomb  z0 dphi_b w_b  dphi_c wd0_c dq_c  dphi_d dq_d pte  hit-sigmas: sphi_d sq_d\n");
  }
  FILE *fd = dump_out.empty() ? nullptr : fopen(dump_out.c_str(), "w");

  std::vector<Stats> S(pats.size());
  Stats SU;  // the union
  Stats SC;  // the chain's own counters and time
  double t_fill = 0;
  long n_unbound = 0;  // labels dropped by --bind
  int npat = pats.size();

  for (int iev = 0; iev < n_events; ++iev) {
    Event ev(iev, ti.n_layers());
    ev.read_in(df);
    const auto t0 = clk::now();
    for (auto &kv : layers)
      kv.second->fill(ev.layerHits_[kv.first]);
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
    for (auto &kv : layers) {
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
    if (chain_holes >= 0) {
      std::vector<std::pair<std::array<int, 4>, Quad>> cq;
      SeedCounters cnt;
      const auto s0 = clk::now();
      std::map<int, const SurfLayer *> LM;
      for (auto &kv : layers)
        LM[kv.first] = kv.second.get();
      CH[0].run(LM, cq, cnt);
      CH[1].run(LM, cq, cnt);
      SC.t_find += secs(s0, clk::now());
      SC.quads += cnt.quads, SC.doublets += cnt.doublets, SC.triplets += cnt.triplets;
      // each quad to the pattern of its layers; a combination no pattern lists becomes a dynamic one
      for (const auto &e : cq) {
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
          chain_q.resize(pats.size());
        chain_q[ip].push_back(e.second);
      }
      chain_q.resize(pats.size());
      npat = pats.size();
    }
    for (int ip = 0; ip < npat; ++ip) {
      const Pattern &p = pats[ip];
      Stats &st = S[ip];
      const SurfLayer *L[4] = {layers[p.l[0]].get(), layers[p.l[1]].get(), layers[p.l[2]].get(), layers[p.l[3]].get()};
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
        {
          // the quad's quality, for the cleaning: its own c and d residuals over the windows
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
                  fprintf(fr, "R %d %d %d %.4f %.3f %d  %.4f %.6g %.6g  %.6g %.6g %.6g  %.6g %.6g %.4g  %.4g %.4g\n", ip, iev,
                          kv.first, t.momEta(), t.pT(), ncomb, e.z0 - t.z(), e.dphi_b, e.w_b, e.dphi_c, e.wd0_c, e.dq_c,
                          e.dphi_d, e.dq_d, e.pte, sphi, sq);
                }
        }
      }
    }
    // the union over patterns, after the cleaning
    std::vector<char> keep_c(cands.size(), 1);
    if (dedup_n > 0) {
      std::vector<int> ord(cands.size());
      for (int i = 0; i < (int)ord.size(); ++i)
        ord[i] = i;
      // tier first (a pattern with an outer-tracker layer after every pure-pixel one: its windows are
      // several times wider, so its scores are not comparable), then the score
      auto tier = [&](int i) {
        int t = 0;
        for (int l : pats[cands[i].ip].l)
          t += !SurfOwnership::is_pix(l);
        return t;
      };
      std::sort(ord.begin(), ord.end(), [&](int x, int y) {
        const int tx = tier(x), ty = tier(y);
        return tx != ty ? tx < ty : cands[x].score < cands[y].score;
      });
      std::unordered_map<long, std::vector<int>> by_hit;  // (layer, hit) -> kept quads using it
      auto hkey = [](int layer, unsigned int k) { return (long)layer << 32 | k; };
      std::unordered_map<int, int> shared;
      for (int i : ord) {
        const Cand &c = cands[i];
        const auto &ll = pats[c.ip].l;
        shared.clear();
        bool drop = false;
        for (int k = 0; k < 4 && !drop; ++k)
          if (auto it = by_hit.find(hkey(ll[k], c.q[k])); it != by_hit.end())
            for (int j : it->second)
              if (++shared[j] >= dedup_n) {
                drop = true;
                break;
              }
        if (drop) {
          keep_c[i] = 0;
          continue;
        }
        for (int k = 0; k < 4; ++k)
          by_hit[hkey(ll[k], c.q[k])].push_back(i);
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
      if (auto it = fd_any.find(l); it != fd_any.end()) {
        su.found[b] += 1;
        su.dup[b] += it->second - 1;
        su.multi[b] += fd_pat[l].size() > 1;
      }
    }
  }
  if (fr)
    fclose(fr);
  if (fd)
    fclose(fd);

  const double ne = n_events;
  printf("[seedsurf] %d events; layer fill %.3f ms/ev\n", n_events, 1e3 * t_fill / ne);
  if (bind_cm > 0)
    printf("[seedsurf] --bind %.3f cm: %.1f label lookups per event dropped (the sim hit farther than that)\n", bind_cm,
           n_unbound / ne);
  if (chain_holes >= 0)
    printf("[seedsurf] CHAIN: %.0f doublets, %.0f triplets, %.1f quads per event, %.3f ms/ev; forwarded %.1f, dropped %.1f per event\n",
           SC.doublets / ne, SC.triplets / ne, SC.quads / ne, 1e3 * SC.t_find / ne,
           (CH[0].n_forwarded + CH[1].n_forwarded) / ne, (CH[0].n_dropped + CH[1].n_dropped) / ne);
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
    }
    fclose(f);
    printf("[seedsurf] wrote %s\n", truth_out.c_str());
  }
  return 0;
}
