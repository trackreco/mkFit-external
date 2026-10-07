// ===========================================================================
// ValProp-Eff.cc -- PER-SIM-TRACK RESULTS: val_eff_* and val_chopres_*.
//
// Split out of ValProp.cc, which had reached 3764 lines holding a dozen
// unrelated instrument families. This is the one that answers "where does the
// gain land": efficiency, fakes, duplicates, track length and momentum
// resolution per sim track, resolved in |eta|, pT and hits-per-layer, with a
// paired per-event error bar and a COMMON-subset restriction.
//
// It is a self-contained unit -- its accumulators, binning and helpers are all
// in the anonymous namespace below, and it shares nothing with the rest of
// ValProp.cc but the includes and the declarations in ValProp.h.
//
// Drivers: test/v2p2-eff-cap.sh (the maxCandsPerSeed scan and the production
// reference), v2p2-eff-score.sh (the score-term study), v2p2-eff-fwd.sh,
// v2p2-eff-chop.sh. Output is <prefix>.root + <prefix>.txt; the deck reads a
// copy of the .root directly, see ~/mic-dev/inlayer-slides/plots/420-eff.html.
// ===========================================================================
#include <array>
#include "RecoTracker/MkFitCore/standalone/RdfTrace/ValProp.h"
#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include "RecoTracker/MkFitCore/standalone/DataFormats/RntConversions.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/cms_common_macros.h"
#include "RecoTracker/MkFitCore/interface/PropagationConfig.h"

#include "RecoTracker/MkFitCore/src/Matrix.h"
#include "RecoTracker/MkFitCore/src/PropagationMPlex.h"
#include "RecoTracker/MkFitCore/src/KalmanUtilsMPlex.h"
#include "RecoTracker/MkFitCore/src/MkFinder.h"
#include "RecoTracker/MkFitCore/src/MkFinderV2p2.h"
#include "RecoTracker/MkFitCore/src/V2p2Score.h"
#include "RecoTracker/MkFitCore/src/MkBins.h"

#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"
#include "RecoTracker/MkFitCore/standalone/TrackExtra.h"
#include "RecoTracker/MkFitCMS/standalone/Shell.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"

#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <array>

#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDSym.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>


namespace mkfit {

  // ==========================================================================
  // val_eff -- per-SIM-TRACK efficiency, resolved in |eta|, pT and hits-per-layer.
  //
  // WHY THIS AND NOT THE STANDARD COUNTERS. Everything measured for the in-layer
  // combinatorial search so far is a scalar summed over a run, so none of it says
  // WHERE the gain lands. The two decisions waiting on that -- whether the search
  // should default ON and at which maxCandsPerSeed, and whether the likelihood
  // score's -ln(rho) occupancy term behaves per region -- are both regional.
  //
  // THE METRIC, and why each piece is safe:
  //  - Association is `TrackExtra::setMCTrackIDInfo`, i.e. `2*mccount >= nCandHits`
  //    over the non-seed hits -- exactly what quality-val's "found tracks" counts.
  //    With val_assoc_mtv(true) it is MTV's rule instead, `4*mccount > 3*nCandHits`
  //    over all hits, seed hits included; the report header says which.
  //    A WRONG extra hit raises the denominator only, so a gain here cannot be
  //    bought by taking more hits. `nH >= 80 %` is NOT used anywhere: it tests raw
  //    reco HITS against sim LAYERS and therefore rewards the thing under test.
  //  - The denominator is SIM tracks, not reco tracks, and it is restricted to sim
  //    tracks that a seed actually points at (`Event::simInfoForCurrentSeed`). That
  //    removes the seeding efficiency, which is common to both configurations and
  //    would otherwise dilute the regional shape without changing the difference.
  //  - Duplicates and fakes come off the same pass, so the three move together and
  //    an efficiency gain paid for in fakes cannot hide.
  //
  // PAIRING. Same events, same seeds, same denominator in every configuration, so
  // the error bar quoted on a difference is the spread of the per-EVENT difference
  // in the numerator, not Poisson on the total. The denominator is checked to be
  // identical across configurations and a mismatch is reported.
  namespace {
    // Axis 0: |eta| of the sim track.  Axis 1: its pT.  Axis 2: hits per layer,
    // i.e. how much overlap the sim track's own hit content offers -- the axis
    // along which the in-layer combinatorial is supposed to pay.
    constexpr int VE_NAX = 3;
    constexpr int VE_NB  = 12;
    const int   ve_nbin[VE_NAX] = {12, 10, 6};
    const char *ve_axname[VE_NAX] = {"|eta|", "pT [GeV]", "sim hits / layer"};

    // THE SELECTION IS CMSSW's MTV CONVENTION: |eta| < 2.5 and pT > 0.9, with the
    // cut on the plotted variable RELEASED for that variable's own plot. So the
    // eta axis carries every selected track of any pT above 0.2, the pT axis
    // carries every selected track of any |eta| below 3, and everything else --
    // the totals, the regions, the hits-per-layer axis, the resolution -- uses
    // both cuts. Quoting one efficiency for "the sample" and another for a bin of
    // its own plot is the convention, not an inconsistency.
    constexpr float VE_ETA_CUT = 2.5f;
    constexpr float VE_PT_CUT  = 0.9f;

    const float ve_pt_edge[11] = {0.2f, 0.3f, 0.5f, 0.7f, 0.9f, 1.2f,
                                  1.6f, 2.5f, 4.0f, 10.0f, 1e9f};
    const char *ve_hpl_lab[6]  = {"= 1.00", "1.0-1.1", "1.1-1.2", "1.2-1.35", "1.35-1.5", "> 1.5"};

    int ve_bin_eta(float ae) { int b = (int)(ae / 0.25f); return (b < 0 || b > 11) ? -1 : b; }
    int ve_bin_pt(float pt) {
      if (pt < ve_pt_edge[0]) return -1;
      for (int b = 0; b < 10; ++b) if (pt < ve_pt_edge[b+1]) return b;
      return 9;
    }
    int ve_bin_hpl(float r) {
      if (r < 1.0001f) return 0;
      if (r <= 1.1f)  return 1;
      if (r <= 1.2f)  return 2;
      if (r <= 1.35f) return 3;
      if (r <= 1.5f)  return 4;
      return 5;
    }
    // Regions, for the per-region read the likelihood score needs. Barrel /
    // transition / endcap by |eta| of the sim track, matching the eta ranges the
    // phase-2 seed partitioner uses closely enough to name them.
    const char *ve_regname[3] = {"barrel   |eta|<0.9", "transition 0.9-1.7", "endcap    >1.7"};
    int ve_region(float ae) { return ae < 0.9f ? 0 : (ae < 1.7f ? 1 : 2); }

    // Uncut reco ranges, as mtvsum.C prints them: all, four |eta| ranges, five pT
    // ranges, and 1.7 < |eta| < 2.7 split at pT 0.9.
    constexpr int VE_NU = 12;
    const char *ve_uname[VE_NU] = {"all", "eta 0-0.9", "eta 0.9-1.7", "eta 1.7-2.7", "eta 2.7-4.5",
                                   "pT 0-0.5", "pT 0.5-0.9", "pT 0.9-2", "pT 2-10", "pT > 10",
                                   "1.7-2.7 pT<0.9", "1.7-2.7 pT>0.9"};
    void ve_ufill(std::array<long, VE_NU> &u, float ae, float pt) {
      ++u[0];
      const float er[5] = {0.f, 0.9f, 1.7f, 2.7f, 4.5f};
      for (int r = 0; r < 4; ++r) if (ae >= er[r] && ae < er[r + 1]) ++u[1 + r];
      const float pr[6] = {0.f, 0.5f, 0.9f, 2.f, 10.f, 1e9f};
      for (int r = 0; r < 5; ++r) if (pt >= pr[r] && pt < pr[r + 1]) ++u[5 + r];
      if (ae >= 1.7f && ae < 2.7f) ++u[pt < 0.9f ? 10 : 11];
    }

    struct VeBins { long b[VE_NAX][VE_NB] = {}; long reg[3] = {}; long tot = 0; };
    struct VeRes { int evt; int lbl; int reg; float r; };
    // Same shape as VeBins but summing a weight, for mean track length. Kept
    // separate rather than templating VeBins, which is counted in longs and is
    // the thing every paired statistic runs on.
    struct VeSums { double b[VE_NAX][VE_NB] = {}; double reg[3] = {}; double tot = 0.0; };

    struct VeCfg {
      std::string name;
      VeBins den, num, dup;        // sim-binned; den is the MTV denominator
      // TRACK LENGTH, over the found sim tracks (so the denominator is num).
      // n_found is every hit on the best-matched reco track, seed hits included;
      // n_match is the subset whose mcTrackID is that sim track, EXCLUDING the
      // seed hits, since setMCTrackIDInfo skips them. n_sim is the sim track's
      // own valid-hit count, which is what both should be read against. The
      // matched one is the safe counter: a wrong extra hit cannot raise it.
      VeSums n_found, n_match, n_sim, n_match2;   // n_match2: sum of squares, for the error on the mean
      VeSums n_sim2, n_seedhits, n_seedhits2;     // sim content and SEED SIZE, same treatment
      VeBins dens;                 // ... of which a seed points at them
      VeBins reco, fake;           // reco-binned (axes 0,1 and region only)
      // Seed quality, per region. The search cannot find what it is not seeded
      // for, and MTV's "central" requirement is in practice imposed by the seeds
      // rather than by us -- so the seeding efficiency and the seed purity are
      // the ceiling every number below sits under, and belong in the same report.
      // How often the shared-hit seed lookup succeeded. If it fails often for one
      // collection its association is computed over more hits than the other's,
      // which is not a like-for-like comparison -- so it is reported, not assumed.
      long n_seed_found = 0, n_seed_missing = 0;
      // What the hit-majority vote makes of each seed. A seed whose vote is -1
      // is NOT counted as seeding its sim track, so these two bound how much the
      // quoted seeding efficiency is understated: a TIE is a seed that does
      // point at a sim track and is being dropped, while no matched hit at all
      // is a seed with no truth to point at.
      long n_vote_ok = 0, n_vote_tie = 0, n_vote_notruth = 0;
      // Seeds here are 100 % pixel hits and 4.1-6.2 of them, i.e. initialStep
      // pixel quadruplets. A selected sim track that does not cross 4 distinct
      // PIXEL layers cannot be seeded by one at all, so counting seeding
      // efficiency against every selected sim track charges this iteration for
      // tracks it is not for. Both denominators are reported.
      long n_pix4[3] = {}, n_pix4_seeded[3] = {}, n_pix4_preseeded[3] = {};
      long n_seed[3] = {}, n_seed_pure[3] = {}, n_seed_on_sel[3] = {};
      double sum_seed_gf[3] = {};
      // What the seeds are MADE OF, which is as close as the .bin gets to naming
      // the seeding algorithm: the file carries the track algorithm and the hits,
      // not the producer. Four pixel hits is a pixel quadruplet either way --
      // Patatrack and the standard chain differ in the fit, not in the hit
      // content -- so this bounds the question rather than settling it.
      long sum_seed_hits[3] = {}, sum_seed_pix[3] = {};
      long n_ev = 0;
      std::vector<VeBins> ev_num, ev_den, ev_fake, ev_reco, ev_dup;
      // Reco tracks WITHOUT the selection cuts, in mtvsum.C's ranges, per event.
      // MTV's fake rate is over every reco track, and the pT > 0.9 cut above
      // hides the low-pT fakes from the totals.
      std::vector<std::array<long, VE_NU>> ev_ureco, ev_ufake;
      // d(pT)/pT of the best-matched reco track of each found sim track, one row
      // per found sim track and keyed by (event, sim label) so the report can
      // restrict every configuration to the tracks they ALL found. That
      // restriction is not a refinement -- unrestricted, a configuration that
      // finds more tracks is measured on a harder population, and the difference
      // that produces is larger than the effect being looked for.
      std::vector<VeRes> res;
    };

    std::vector<VeCfg> g_ve;
    std::string g_ve_ref;
    // A SECOND reference, differenced against in its own set of histograms. The
    // tables keep g_ve_ref; this one exists so a plot can show "against
    // production" underneath while the tables still read against what this work
    // replaced. Deltas are PAIRED (shared denominator, per-event numerators).
    std::string g_ve_ref2 = "cmssw_V1";
    FILE *g_ve_log = nullptr;

    int ve_printf(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
    int ve_printf(const char *fmt, ...) {
      va_list ap;  va_start(ap, fmt);  int n = vprintf(fmt, ap);  va_end(ap);
      if (g_ve_log) { va_start(ap, fmt); vfprintf(g_ve_log, fmt, ap); va_end(ap); }
      return n;
    }

    VeCfg &ve_cfg(const char *name) {
      for (auto &c : g_ve) if (c.name == name) return c;
      g_ve.push_back(VeCfg());  g_ve.back().name = name;  return g_ve.back();
    }

    // ae, pt gate which axes this track enters, per the MTV convention above:
    // the eta axis wants the pT cut only, the pT axis the eta cut only, and
    // everything else both.
    void ve_fill(VeBins &v, int be, int bp, int bh, int reg, float ae, float pt) {
      const bool ok_eta = ae < VE_ETA_CUT, ok_pt = pt > VE_PT_CUT;
      if (be >= 0 && ok_pt) ++v.b[0][be];
      if (bp >= 0 && ok_eta) ++v.b[1][bp];
      if (!(ok_eta && ok_pt)) return;
      if (bh >= 0) ++v.b[2][bh];
      if (reg >= 0) ++v.reg[reg];
      ++v.tot;
    }
    void ve_fill_w(VeSums &v, int be, int bp, int bh, int reg, float ae, float pt, double w) {
      const bool ok_eta = ae < VE_ETA_CUT, ok_pt = pt > VE_PT_CUT;
      if (be >= 0 && ok_pt) v.b[0][be] += w;
      if (bp >= 0 && ok_eta) v.b[1][bp] += w;
      if (!(ok_eta && ok_pt)) return;
      if (bh >= 0) v.b[2][bh] += w;
      if (reg >= 0) v.reg[reg] += w;
      v.tot += w;
    }
    void ve_add(VeBins &a, const VeBins &b) {
      for (int x = 0; x < VE_NAX; ++x) for (int i = 0; i < VE_NB; ++i) a.b[x][i] += b.b[x][i];
      for (int r = 0; r < 3; ++r) a.reg[r] += b.reg[r];
      a.tot += b.tot;
    }
    // Paired: mean and sigma-on-the-sum of the per-event difference num(A) - num(B).
    void ve_paired(const std::vector<long> &a, const std::vector<long> &b, double &sum, double &sig) {
      const size_t n = std::min(a.size(), b.size());
      sum = 0.0;  sig = 0.0;
      if (n < 2) return;
      double s = 0.0, s2 = 0.0;
      for (size_t i = 0; i < n; ++i) { const double d = (double)a[i] - (double)b[i]; s += d; s2 += d*d; }
      sum = s;
      const double mean = s / n;
      const double var = (s2 - n*mean*mean) / (n - 1);
      sig = std::sqrt(std::max(0.0, var) * n);   // sigma on the SUM of the n differences
    }
  }

  namespace {
    // Which seed produced this track, found by SHARED HITS rather than by label.
    //
    // It has to be done this way for cmsswTracks_, and is then done this way for
    // ours too so that the association rule is identical for both. The labels
    // cannot be joined: WriteMemoryFile gives a seed the label seedSimIdx[is] --
    // a SIM track index -- while it gives a rec track the label
    // trk_seedIdx->at(ir), the NTUPLE's seed index, and the written seed vector
    // skips every seed that is not initialStep or hltIter0, so the two indices
    // are not the same namespace and the mapping is not in the file.
    //
    // This matters rather than being pedantry: setMCTrackIDInfo excludes seed
    // hits from the association count, and excluding them makes association
    // HARDER (removing a matched hit takes 2m >= n to 2(m-1) >= n-1). Failing to
    // identify a collection's seeds would therefore flatter that collection.
    int ve_seed_of_track(const TrackVec &seeds,
                         const std::map<std::pair<int,int>, std::vector<int>> &hit2seed,
                         const Track &t) {
      std::map<int, int> shared;
      for (int i = 0; i < t.nTotalHits(); ++i) {
        const HitOnTrack hot = t.getHitOnTrack(i);
        if (hot.index < 0 || hot.layer < 0) continue;
        auto it = hit2seed.find({hot.layer, hot.index});
        if (it == hit2seed.end()) continue;
        for (int si : it->second) ++shared[si];
      }
      int best = -1, bn = 0;
      for (const auto &[si, n] : shared)
        if (n > bn || (n == bn && best >= 0 && seeds[si].nFoundHits() < seeds[best].nFoundHits())) {
          if (n > bn) { bn = n; best = si; }
        }
      return bn > 0 ? best : -1;
    }
  }

  // Signal-only denominator: only sim tracks with Track::ProdType::Signal (in-time, event 0),
  // as CMSSW MTV's efficiency (signalOnly). The reco side, and so the fake rate, is unchanged.
  static bool g_ve_signal_only = false;
  void val_eff_signal_only(bool on) {
    g_ve_signal_only = on;
    printf("val_eff_signal_only: %d (efficiency denominator %s)\n", (int) on,
           on ? "signal sim tracks only" : "every selected sim track");
  }

  void val_eff_reset() { g_ve.clear(); g_ve_ref.clear(); }
  void val_eff_ref(const char *cfg) { g_ve_ref = cfg; }
  void val_eff_ref2(const char *cfg) { g_ve_ref2 = cfg; }

  // tracks is the collection under test. It is candidateTracks_ for a v2p2
  // configuration and cmsswTracks_ for the production reference; everything else
  // -- selection, association rule, seed-hit exclusion -- is identical, which is
  // the only way the two can be put on one axis.
  static void ve_accumulate(const Event *ev, const TrackVec &tracks, VeCfg &C) {
    VeBins e_den, e_num, e_dup, e_reco, e_fake, e_dens;
    std::array<long, VE_NU> e_ureco{}, e_ufake{};

    // ---- seeds: which sim tracks did the search actually get a chance at, and
    // where is each seed, so its hits can be excluded from the association count.
    std::map<int, int> seed_by_label;      // seed label -> index in currentSeedTracks()
    std::set<int> seeded_sim;
    std::map<int, int> n_seed_for_sim;     // sim label -> how many seeds point at it
    // The same vote over the file's RAW seed collection, before seed cleaning.
    // currentSeedTracks() is what the cleaner left -- 7498 of 36042 on this
    // sample -- so without this the quoted seeding efficiency cannot tell a seed
    // that was never made from one that was made and then thrown away.
    std::set<int> preclean_sim;
    for (const Track &sd : ev->seedTracks_) {
      const int sl = ev->simInfoForTrack(sd).label;
      if (sl >= 0 && sl < (int) ev->simTracks_.size()) preclean_sim.insert(sl);
    }
    const TrackVec *seeds = nullptr;
    try { seeds = &ev->currentSeedTracks(); } catch (...) { seeds = nullptr; }
    if (seeds) {
      for (int i = 0; i < (int) seeds->size(); ++i) {
        seed_by_label.emplace((*seeds)[i].label(), i);
        const auto sifh = ev->simInfoForCurrentSeed(i);
        const int sl = sifh.label;
        if (sl >= 0) ++C.n_vote_ok;
        else if (sifh.n_match > 0) ++C.n_vote_tie;      // a tie -- Event::simInfoForTrack returns -1
        else ++C.n_vote_notruth;
        if (sl >= 0 && sl < (int) ev->simTracks_.size()) {
          seeded_sim.insert(sl);
          ++n_seed_for_sim[sl];
        }
        const int sr = ve_region(std::abs((*seeds)[i].momEta()));
        if (sr >= 0 && sr < 3) {
          ++C.n_seed[sr];
          C.sum_seed_gf[sr] += sifh.good_frac();
          if (sifh.good_frac() > 0.999f) ++C.n_seed_pure[sr];
          const Track &sd = (*seeds)[i];
          for (int h = 0; h < sd.nTotalHits(); ++h) {
            const HitOnTrack hot = sd.getHitOnTrack(h);
            if (hot.index < 0 || hot.layer < 0) continue;
            ++C.sum_seed_hits[sr];
            if (Config::TrkInfo[hot.layer].is_pixel()) ++C.sum_seed_pix[sr];
          }
        }
      }
    }

    // ---- reco side: association exactly as quality-val defines it.
    std::map<int, int> n_assoc;            // sim label -> number of reco tracks on it
    struct VeBest { float pt = 0.0f; int n_match = -1; int n_found = 0; int n_seedhits = 0; };
    std::map<int, VeBest> best;            // sim label -> its best-matched reco track
    // (layer, index) -> the seeds holding that hit, for the shared-hit seed
    // lookup. Built once per event; seeds are a few hits each, so it is small.
    std::map<std::pair<int,int>, std::vector<int>> hit2seed;
    if (seeds)
      for (int i = 0; i < (int) seeds->size(); ++i)
        for (int h = 0; h < (*seeds)[i].nTotalHits(); ++h) {
          const HitOnTrack hot = (*seeds)[i].getHitOnTrack(h);
          if (hot.index >= 0 && hot.layer >= 0) hit2seed[{hot.layer, hot.index}].push_back(i);
        }
    for (const Track &c : tracks) {
      TrackExtra extra(c.label());
      const int si = seeds ? ve_seed_of_track(*seeds, hit2seed, c) : -1;
      if (si >= 0) { extra.findMatchingSeedHits(c, (*seeds)[si], ev->layerHits_); ++C.n_seed_found; }
      else ++C.n_seed_missing;
      extra.setMCTrackIDInfo(c, ev->layerHits_, ev->simHitsInfo_, ev->simTracks_, false, false);
      const int mc = extra.mcTrackID();

      const float rae = std::abs(c.momEta());
      const int rbe = ve_bin_eta(rae), rbp = ve_bin_pt(c.pT()), rreg = ve_region(rae);
      ve_fill(e_reco, rbe, rbp, -1, rreg, rae, c.pT());
      ve_ufill(e_ureco, rae, c.pT());
      // modifyRefTrackID() codes: >= 0 matched, findable; -2 matched but short;
      // -3/-4 matched to a sim track that is not findable (in the shell: has no
      // seed, i.e. most pileup); -1 and -5..-9 matched to NO sim track; -10
      // duplicate. MTV's fake is "no sim track by hits", so in MTV mode only
      // -1 and -5..-9 are fakes. The default mode keeps every negative code, as
      // all numbers recorded before 2026-09-28 were made that way.
      const bool fake = Config::mtvLikeValidation ? (mc == -1 || (mc <= -5 && mc >= -9))
                                                  : (mc < 0 || mc >= (int) ev->simTracks_.size());
      if (fake) {
        ve_fill(e_fake, rbe, rbp, -1, rreg, rae, c.pT());
        ve_ufill(e_ufake, rae, c.pT());
      } else if (mc < 0)
        ;  // matched to a sim track that is not in the efficiency denominator
      else {
        ++n_assoc[mc];
        // Keep the best-matched reco track per sim track, so a duplicate does not
        // get to vote twice on the resolution.
        auto &b = best[mc];
        if (extra.nHitsMatched() > b.n_match)
          b = {c.pT(), extra.nHitsMatched(), c.nFoundHits(),
               si >= 0 ? (*seeds)[si].nFoundHits() : 0};
      }
    }

    // ---- sim side. The denominator is EVERY selected sim track, not only the
    // seeded ones, which is what MTV means by efficiency. The seeded subset is
    // counted alongside so the seeding ceiling is visible rather than assumed.
    //
    // Selection, CMSSW TrackingParticleSelector's shape: findable, the production
    // vertex central (tip < 3.5 cm, lip < 30 cm), and at least 4 distinct layers,
    // since a seed needs four and a sim track with fewer is not reconstructible
    // by this algorithm at all. The eta and pT cuts are applied per axis inside
    // ve_fill(), so each plot releases the cut on its own variable.
    const int nsim = (int) ev->simTracks_.size();
    for (int L = 0; L < nsim; ++L) {
      const Track &st = ev->simTracks_[L];
      const float ae = std::abs(st.momEta()), pt = st.pT();
      if (!st.isFindable()) continue;
      if (g_ve_signal_only && st.prodType() != Track::ProdType::Signal) continue;
      if (std::hypot(st.x(), st.y()) > 3.5f || std::abs(st.z()) > 30.0f) continue;
      if (ae >= 3.0f || pt < ve_pt_edge[0]) continue;
      const int nlay = st.nUniqueLayers();
      if (nlay < 4) continue;
      int nval = 0;
      for (int i = 0; i < st.nTotalHits(); ++i)
        if (st.getHitOnTrack(i).index >= 0) ++nval;
      const int be = ve_bin_eta(ae), bp = ve_bin_pt(pt), reg = ve_region(ae);
      const int bh = ve_bin_hpl((float) nval / (float) nlay);
      std::set<int> pixlay;
      for (int i = 0; i < st.nTotalHits(); ++i) {
        const HitOnTrack hot = st.getHitOnTrack(i);
        if (hot.index >= 0 && hot.layer >= 0 && Config::TrkInfo[hot.layer].is_pixel())
          pixlay.insert(hot.layer);
      }
      if (reg >= 0 && ae < VE_ETA_CUT && pt > VE_PT_CUT && pixlay.size() >= 4) {
        ++C.n_pix4[reg];
        if (seeded_sim.count(L)) ++C.n_pix4_seeded[reg];
        if (preclean_sim.count(L)) ++C.n_pix4_preseeded[reg];
      }
      ve_fill(e_den, be, bp, bh, reg, ae, pt);
      if (seeded_sim.count(L)) {
        ve_fill(e_dens, be, bp, bh, reg, ae, pt);
        if (reg >= 0 && ae < VE_ETA_CUT && pt > VE_PT_CUT)
          C.n_seed_on_sel[reg] += n_seed_for_sim[L];
      }
      auto it = n_assoc.find(L);
      if (it != n_assoc.end()) {
        ve_fill(e_num, be, bp, bh, reg, ae, pt);
        for (int d = 1; d < it->second; ++d) ve_fill(e_dup, be, bp, bh, reg, ae, pt);
        auto bi = best.find(L);
        if (bi != best.end()) {
          ve_fill_w(C.n_found, be, bp, bh, reg, ae, pt, bi->second.n_found);
          const double nm = std::max(0, bi->second.n_match);
          ve_fill_w(C.n_match, be, bp, bh, reg, ae, pt, nm);
          ve_fill_w(C.n_match2, be, bp, bh, reg, ae, pt, nm*nm);
          ve_fill_w(C.n_sim,   be, bp, bh, reg, ae, pt, nval);
          ve_fill_w(C.n_sim2,  be, bp, bh, reg, ae, pt, (double) nval * nval);
          ve_fill_w(C.n_seedhits,  be, bp, bh, reg, ae, pt, bi->second.n_seedhits);
          ve_fill_w(C.n_seedhits2, be, bp, bh, reg, ae, pt,
                    (double) bi->second.n_seedhits * bi->second.n_seedhits);
          if (pt > 0.0f && ae < VE_ETA_CUT && pt > VE_PT_CUT)
            C.res.push_back({ev->evtID(), L, reg, (bi->second.pt - pt) / pt});
        }
      }
    }

    ve_add(C.den, e_den);  ve_add(C.num, e_num);  ve_add(C.dup, e_dup);
    ve_add(C.dens, e_dens);
    ve_add(C.reco, e_reco); ve_add(C.fake, e_fake);
    C.ev_den.push_back(e_den);  C.ev_num.push_back(e_num);  C.ev_dup.push_back(e_dup);
    C.ev_reco.push_back(e_reco); C.ev_fake.push_back(e_fake);
    C.ev_ureco.push_back(e_ureco); C.ev_ufake.push_back(e_ufake);
    ++C.n_ev;
  }

  void val_eff_event(const Event *ev, const char *cfg) {
    if (ev == nullptr) return;
    ve_accumulate(ev, ev->candidateTracks_, ve_cfg(cfg));
  }

  // The production reference: whatever tracking the job that wrote the ntuple
  // ran, which for these samples is mkFit V1 with prop-to-plane and
  // selectHitIndicesV2. Needs --read-cmssw-tracks AND a .bin converted with
  // --write-rec-tracks; without the latter the section is not in the file at all
  // and this reports an empty collection rather than failing quietly.
  //--------------------------------------------------------------------------
  // val_hitpur -- per mkFit layer, the hits on the reconstructed tracks: how
  // many belong to the track's own particle (TrackExtra's mcTrackID, MTV rule
  // when val_assoc_mtv is on) and how many to another or to none. Only tracks
  // matched to a particle (mcTrackID >= 0) are counted, so "wrong" means a
  // foreign hit on a track that is otherwise that particle's. Seed hits
  // included. For comparing layer plans: where do the wrong hits enter.
  //--------------------------------------------------------------------------
  namespace {
    struct HpCfg {
      std::string name;
      long n_trk = 0;
      long own[64] = {}, wrong[64] = {}, none[64] = {};
    };
    std::vector<HpCfg> g_hp;
  }  // namespace

  //============================================================================
  // val_rank -- the worst end-of-layer selection rank along the ancestry of each
  // seed's final forward-search candidate (TrCandMeta::fwd_rank_*, MKFIT_TRACE,
  // MkFinderV2p2). The layer-step score is additive, so a candidate whose
  // ancestors all ranked <= R at beam width W > R also survives at width R: its
  // competitors at width R are a subset of those at W. So in a run at width 10,
  // the final tracks with fwd_rank_max > 6 are the ones width 6 cannot make.
  // Track-to-meta by TrCandMeta::cand, the export index: run WITHOUT
  // --remove-dup, which erases tracks and shifts it. Truth as val_eff, by the
  // association mode in force (val_assoc_mtv).
  //============================================================================

  namespace {
    struct RkCfg {
      std::string name;
      // [class][rank bin]: class 0 good, 1 fake, 2 other; rank bin 0: <= rank_over,
      // 1..4: rank_over + 1 .. + 4 and above.
      long n[3][5] = {};
      // first layer above rank_over, by group, per class
      long over_grp[3][9] = {};
      // per |eta| region of the reco track: tracks, of them above rank_over
      long reg[3][3][2] = {};
      long n_nometa = 0;
      // per event (in call order), seed index -> (class, hash of the hit list,
      // worst rank), for the per-seed comparison against the first config.
      struct Fin { int cls; size_t hash; int rank; unsigned long long fwd_hash; };
      std::vector<std::map<int, Fin>> per_ev;
    };
    std::vector<RkCfg> g_rk;

    int rk_group(int l) {
      // 0 PixB, 1 TBPS 1st, 2 TBPS 2nd, 3 TB2S 1st, 4 TB2S 2nd, 5 FPix, 6 TEDD 1st, 7 TEDD 2nd
      if (l < 0) return 8;
      if (l >= 38 && l <= 59) l -= 22;
      if (l <= 3) return 0;
      if (l <= 9) return 1 + (l & 1);
      if (l <= 15) return 3 + (l & 1);
      if (l <= 27) return 5;
      return 6 + (l & 1);
    }
  }

  void val_rank_reset() { g_rk.clear(); }

  void val_rank_event(const Event *ev, const char *cfg) {
    if (ev == nullptr) return;
    RkCfg *C = nullptr;
    for (auto &c : g_rk) if (c.name == cfg) C = &c;
    if (!C) { g_rk.push_back(RkCfg()); g_rk.back().name = cfg; C = &g_rk.back(); }
    const int R = Config::V2p2::Diag::rank_over;
    C->per_ev.emplace_back();
    auto &fin = C->per_ev.back();
    std::vector<int> meta_of(ev->candidateTracks_.size(), -1);
    for (int m = 0; m < (int) ev->trCandMetas_.size(); ++m) {
      const int ci = ev->trCandMetas_[m].cand;
      if (ci >= 0 && ci < (int) meta_of.size()) meta_of[ci] = m;
    }
    for (int ci = 0; ci < (int) ev->candidateTracks_.size(); ++ci) {
      if (meta_of[ci] < 0) { ++C->n_nometa; continue; }
      const TrCandMeta &cm = ev->trCandMetas_[meta_of[ci]];
      const Track &c = ev->candidateTracks_[ci];
      TrackExtra extra(c.label());
      extra.setMCTrackIDInfo(c, ev->layerHits_, ev->simHitsInfo_, ev->simTracks_, false, false);
      const int mc = extra.mcTrackID();
      const bool fake = Config::mtvLikeValidation ? (mc == -1 || (mc <= -5 && mc >= -9))
                                                  : (mc < 0 || mc >= (int) ev->simTracks_.size());
      const int cls = fake ? 1 : (mc >= 0 && mc < (int) ev->simTracks_.size()) ? 0 : 2;
      const int rb = cm.fwd_rank_max <= R ? 0 : std::min(4, cm.fwd_rank_max - R);
      ++C->n[cls][rb];
      if (rb > 0) ++C->over_grp[cls][rk_group(cm.fwd_rank_over_layer)];
      const float ae = std::abs(c.momEta());
      const int reg = ae < 0.9f ? 0 : ae < 1.7f ? 1 : 2;
      ++C->reg[cls][reg][0];
      if (rb > 0) ++C->reg[cls][reg][1];
      fin[cm.seed] = {cls, hits_hash(c), cm.fwd_rank_max, cm.fwd_hits_hash};
    }
  }

  void val_rank_report() {
    const int R = Config::V2p2::Diag::rank_over;
    const char *cls_name[3] = {"good", "fake", "other"};
    const char *grp_name[9] = {"PixB", "TBPS1", "TBPS2", "TB2S1", "TB2S2", "FPix", "TEDD1", "TEDD2", "none"};
    printf("\n--- val_rank: worst end-of-layer rank along the final candidate's ancestry (rank_over %d) ---\n", R);
    for (auto &C : g_rk) {
      printf("config %s (tracks without a trace meta: %ld)\n", C.name.c_str(), C.n_nometa);
      printf("  %-6s %9s %9s   rank %d / %d / %d / >=%d\n", "class", "tracks", "<=", R + 1, R + 2, R + 3, R + 4);
      for (int k = 0; k < 3; ++k) {
        long t = 0;
        for (int b = 0; b < 5; ++b) t += C.n[k][b];
        printf("  %-6s %9ld %9ld   %6ld %6ld %6ld %6ld\n", cls_name[k], t, C.n[k][0], C.n[k][1], C.n[k][2],
               C.n[k][3], C.n[k][4]);
      }
      printf("  first layer above rank %d, by group:\n  %-6s", R, "class");
      for (int g = 0; g < 9; ++g) printf(" %6s", grp_name[g]);
      printf("\n");
      for (int k = 0; k < 3; ++k) {
        printf("  %-6s", cls_name[k]);
        for (int g = 0; g < 9; ++g) printf(" %6ld", C.over_grp[k][g]);
        printf("\n");
      }
      printf("  by reco |eta|, tracks / above rank %d:   barrel        transition    endcap\n", R);
      for (int k = 0; k < 3; ++k)
        printf("  %-6s %32s %6ld / %-5ld %6ld / %-5ld %6ld / %-5ld\n", cls_name[k], "", C.reg[k][0][0],
               C.reg[k][0][1], C.reg[k][1][0], C.reg[k][1][1], C.reg[k][2][0], C.reg[k][2][1]);
    }
    // Per seed against the first configuration: the same seeds in every config.
    if (g_rk.size() < 2) return;
    const RkCfg &A = g_rk[0];
    for (size_t ic = 1; ic < g_rk.size(); ++ic) {
      const RkCfg &B = g_rk[ic];
      // [class A][class B][same hits][B rank above]; tf: the same, split by the
      // forward search's final candidate instead of the final track.
      long t[3][3][2][2] = {}, tf[3][3][2][2] = {};
      long only_a = 0, only_b = 0;
      const size_t ne = std::min(A.per_ev.size(), B.per_ev.size());
      for (size_t e = 0; e < ne; ++e) {
        for (auto &[s, fa] : A.per_ev[e]) {
          auto it = B.per_ev[e].find(s);
          if (it == B.per_ev[e].end()) { ++only_a; continue; }
          const auto &fb = it->second;
          ++t[fa.cls][fb.cls][fa.hash == fb.hash][fb.rank > R];
          ++tf[fa.cls][fb.cls][fa.fwd_hash == fb.fwd_hash][fb.rank > R];
        }
        for (auto &[s, fb] : B.per_ev[e])
          if (A.per_ev[e].find(s) == A.per_ev[e].end()) ++only_b;
      }
      printf("per seed, %s -> %s (%zu events; seeds only in %s %ld, only in %s %ld):\n", A.name.c_str(),
             B.name.c_str(), ne, A.name.c_str(), only_a, B.name.c_str(), only_b);
      printf("  %-14s %10s %10s %12s %12s\n", "", "same hits", "", "other hits", "");
      printf("  %-14s %10s %10s %12s %12s\n", "", "rank<=", "rank>", "rank<=", "rank>");
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
          printf("  %-5s -> %-5s  %10ld %10ld %12ld %12ld\n", cls_name[a], cls_name[b], t[a][b][1][0],
                 t[a][b][1][1], t[a][b][0][0], t[a][b][0][1]);
      printf("  split by the FORWARD search's final candidate (same / other forward hits):\n");
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
          printf("  %-5s -> %-5s  %10ld %10ld %12ld %12ld\n", cls_name[a], cls_name[b], tf[a][b][1][0],
                 tf[a][b][1][1], tf[a][b][0][0], tf[a][b][0][1]);
    }
  }

  void val_hitpur_reset() { g_hp.clear(); }

  void val_hitpur_event(const Event *ev, const char *cfg) {
    if (ev == nullptr) return;
    HpCfg *C = nullptr;
    for (auto &c : g_hp) if (c.name == cfg) C = &c;
    if (!C) { g_hp.push_back(HpCfg()); g_hp.back().name = cfg; C = &g_hp.back(); }
    for (const Track &c : ev->candidateTracks_) {
      TrackExtra extra(c.label());
      extra.setMCTrackIDInfo(c, ev->layerHits_, ev->simHitsInfo_, ev->simTracks_, false, false);
      const int mc = extra.mcTrackID();
      if (mc < 0 || mc >= (int) ev->simTracks_.size()) continue;
      ++C->n_trk;
      for (int h = 0; h < c.nTotalHits(); ++h) {
        const HitOnTrack hot = c.getHitOnTrack(h);
        if (hot.index < 0 || hot.layer < 0 || hot.layer >= 64) continue;
        const Hit &hit = ev->layerHits_[hot.layer][hot.index];
        const int id = hit.mcHitID();
        const int st = (id >= 0 && id < (int) ev->simHitsInfo_.size()) ? ev->simHitsInfo_[id].mcTrackID() : -1;
        if (st == mc) ++C->own[hot.layer];
        else if (st < 0) ++C->none[hot.layer];
        else ++C->wrong[hot.layer];
      }
    }
  }

  void val_hitpur_report() {
    printf("\n--- val_hitpur: hits on particle-matched tracks, per mkFit layer: own / foreign / none ---\n");
    printf("%-6s", "layer");
    for (auto &c : g_hp) printf(" | %-30s", c.name.c_str());
    printf("\n");
    for (int l = 0; l < 64; ++l) {
      bool any = false;
      for (auto &c : g_hp) any |= (c.own[l] + c.wrong[l] + c.none[l]) > 0;
      if (!any) continue;
      printf("%-6d", l);
      for (auto &c : g_hp) {
        const long t = c.own[l] + c.wrong[l] + c.none[l];
        printf(" | %8ld %6.2f%% %5.2f%% %5.3f", c.own[l], t ? 100.0 * c.wrong[l] / t : 0.0,
               t ? 100.0 * c.none[l] / t : 0.0, c.n_trk ? (double) t / c.n_trk : 0.0);
      }
      printf("\n");
    }
    printf("  per config: own hits, foreign %%, none %%, hits per matched track; tracks:");
    for (auto &c : g_hp) printf(" %s %ld", c.name.c_str(), c.n_trk);
    printf("\n");
  }

  void val_eff_cmssw_event(const Event *ev, const char *cfg) {
    if (ev == nullptr) return;
    if (ev->cmsswTracks_.empty()) {
      static bool warned = false;
      if (!warned) {
        printf("val_eff_cmssw_event: cmsswTracks_ is EMPTY. Needs --read-cmssw-tracks,\n"
               "  and a .bin written with --write-rec-tracks -- check the header's\n"
               "  'Extra sections' line for CmsswTracks.\n");
        warned = true;
      }
      return;
    }
    ve_accumulate(ev, ev->cmsswTracks_, ve_cfg(cfg));
  }

  namespace {
    // Pull one bin's per-event series out of a config, for the paired statistics.
    std::vector<long> ve_series(const std::vector<VeBins> &ev, int ax, int bin) {
      std::vector<long> v;  v.reserve(ev.size());
      for (const auto &e : ev) v.push_back(ax < 0 ? e.tot : (ax == 3 ? e.reg[bin] : e.b[ax][bin]));
      return v;
    }
    std::string ve_binlabel(int ax, int b) {
      char s[32];
      if (ax == 0) snprintf(s, sizeof(s), "%.2f-%.2f", 0.25*b, 0.25*(b+1));
      else if (ax == 1) {
        if (b == 9) snprintf(s, sizeof(s), "> 10");
        else snprintf(s, sizeof(s), "%.1f-%.1f", ve_pt_edge[b], ve_pt_edge[b+1]);
      }
      else snprintf(s, sizeof(s), "%s", ve_hpl_lab[b]);
      return s;
    }

    // Fill a TH1 with a mean and its ERROR ON THE MEAN from the running sums.
    // Every curve on these plots is a measurement and carries one; drawing a
    // reference line without errors invites it to be read as exact.
    void ve_mean_hist(TH1D *h, int ax, const VeSums &sum, const VeSums &sum2,
                      const VeBins &cnt, int nb) {
      for (int b = 0; b < nb; ++b) {
        const long n = cnt.b[ax][b];
        if (n < 20) continue;
        const double m = sum.b[ax][b] / n;
        const double v = std::max(0.0, sum2.b[ax][b] / n - m * m);
        h->SetBinContent(b + 1, m);
        h->SetBinError(b + 1, std::sqrt(v / n));
      }
    }

    // Mean track length per bin. which = 0 reco hits (seed included), 1 matched
    // hits (seed excluded, the safe counter), 2 the sim track's own hits. The
    // denominator is the FOUND sim tracks, so it is a property of the tracks a
    // configuration reconstructed and moves with the population -- read the
    // matched row, and read it against the sim row in the same column.
    void ve_len_table(int ax, const VeCfg *ref, int which) {
      const int nb = (ax == 3) ? 3 : ve_nbin[ax];
      const char *wn[3] = {"reco hits (incl. seed)", "MATCHED hits (excl. seed)",
                           "sim track's own hits"};
      ve_printf("\n--- mean %s vs %s ---\n", wn[which],
                ax == 3 ? "region" : ve_axname[ax]);
      ve_printf("%-20s %9s", ax == 3 ? "region" : "bin", "found");
      for (const auto &c : g_ve) ve_printf(" | %-10.10s", c.name.c_str());
      ve_printf("\n");
      for (int b = 0; b < nb; ++b) {
        const long dn = (ax == 3) ? ref->num.reg[b] : ref->num.b[ax][b];
        if (dn < 20) continue;
        std::string lab = (ax == 3) ? ve_regname[b] : ve_binlabel(ax, b);
        ve_printf("%-20s %9ld", lab.c_str(), dn);
        for (const auto &c : g_ve) {
          const VeSums &S = (which == 0) ? c.n_found : (which == 1 ? c.n_match : c.n_sim);
          const long n = (ax == 3) ? c.num.reg[b] : c.num.b[ax][b];
          const double v = (ax == 3) ? S.reg[b] : S.b[ax][b];
          ve_printf(" |   %8.3f ", n ? v/n : 0.0);
        }
        ve_printf("\n");
      }
    }

    // One resolved table: efficiency per bin for every configuration, and the
    // paired difference of each against the reference.
    void ve_table(int ax, const VeCfg *ref) {
      const int nb = (ax == 3) ? 3 : ve_nbin[ax];
      ve_printf("\n--- efficiency vs %s ---\n", ax == 3 ? "region" : ve_axname[ax]);
      ve_printf("%-20s %9s", ax == 3 ? "region" : "bin", "sim trks");
      for (const auto &c : g_ve) ve_printf(" | %-10.10s", c.name.c_str());
      ve_printf("\n");
      for (int b = 0; b < nb; ++b) {
        const long den = (ax == 3) ? ref->den.reg[b] : ref->den.b[ax][b];
        if (den < 20) continue;
        std::string lab = (ax == 3) ? ve_regname[b] : ve_binlabel(ax, b);
        ve_printf("%-20s %9ld", lab.c_str(), den);
        for (const auto &c : g_ve) {
          const long d = (ax == 3) ? c.den.reg[b] : c.den.b[ax][b];
          const long n = (ax == 3) ? c.num.reg[b] : c.num.b[ax][b];
          if (&c == ref) ve_printf(" |   %6.2f%% ", d ? 100.0*n/d : 0.0);
          else {
            double sum, sig;
            ve_paired(ve_series(c.ev_num, ax == 3 ? 3 : ax, b),
                      ve_series(ref->ev_num, ax == 3 ? 3 : ax, b), sum, sig);
            const double dp = d ? 100.0*sum/d : 0.0;          // difference in points
            const double sp = d ? 100.0*sig/d : 0.0;
            ve_printf(" | %+6.2f%s%-3.3s", dp,
                      sp > 0 && std::abs(dp) > 3*sp ? "*" : " ",
                      sp > 0 ? (std::abs(dp) > 3*sp ? "sig" : "") : "");
          }
        }
        ve_printf("\n");
      }
      ve_printf("  reference column is the absolute efficiency; the others are the\n"
                "  PAIRED difference in points, * = |delta| > 3 sigma of the per-event spread.\n");
    }

    void ve_ratio_table(int ax, const VeCfg *ref, bool fakes) {
      const int nb = (ax == 3) ? 3 : ve_nbin[ax];
      ve_printf("\n--- %s vs %s ---\n",
                fakes ? "FAKE fraction, per RECO track" : "EXTRA reco tracks per found SIM track",
                ax == 3 ? "region" : ve_axname[ax]);
      ve_printf("%-20s", ax == 3 ? "region" : "bin");
      for (const auto &c : g_ve) ve_printf(" | %-10.10s", c.name.c_str());
      ve_printf("\n");
      for (int b = 0; b < nb; ++b) {
        const long dref = (ax == 3) ? (fakes ? ref->reco.reg[b] : ref->den.reg[b])
                                    : (fakes ? ref->reco.b[ax][b] : ref->den.b[ax][b]);
        if (dref < 20) continue;
        std::string lab = (ax == 3) ? ve_regname[b] : ve_binlabel(ax, b);
        ve_printf("%-20s", lab.c_str());
        for (const auto &c : g_ve) {
          const long d = (ax == 3) ? (fakes ? c.reco.reg[b] : c.den.reg[b])
                                   : (fakes ? c.reco.b[ax][b] : c.den.b[ax][b]);
          const long n = (ax == 3) ? (fakes ? c.fake.reg[b] : c.dup.reg[b])
                                   : (fakes ? c.fake.b[ax][b] : c.dup.b[ax][b]);
          ve_printf(" |   %6.2f%% ", d ? 100.0*n/d : 0.0);
        }
        ve_printf("\n");
      }
    }
  }

  void val_eff_report(const char *prefix) {
    if (g_ve.empty()) { printf("val_eff_report: nothing accumulated.\n"); return; }
    const VeCfg *ref = &g_ve[0];
    if (!g_ve_ref.empty())
      for (const auto &c : g_ve) if (c.name == g_ve_ref) ref = &c;

    const std::string txt = std::string(prefix) + ".txt";
    g_ve_log = fopen(txt.c_str(), "w");

    ve_printf("\n================================================================\n");
    ve_printf("  val_eff -- per-sim-track efficiency, resolved\n");
    ve_printf("================================================================\n");
    ve_printf("DENOMINATOR -- CMSSW MTV's convention. Sim tracks that are findable, have a\n");
    ve_printf("  CENTRAL production vertex (tip < 3.5 cm, lip < 30 cm, as\n");
    ve_printf("  TrackingParticleSelector has it), at least 4 distinct layers, |eta| < 2.5\n");
    ve_printf("  and pT > 0.9 -- with the cut on the PLOTTED variable released for that\n");
    ve_printf("  variable's own plot, so the eta table carries every pT above 0.2 and the pT\n");
    ve_printf("  table every |eta| below 3. Efficiency below is against ALL of them, not\n");
    ve_printf("  only the seeded ones; the seeded subset is in the next block, because the\n");
    ve_printf("  search cannot find what it was not seeded for and MTV's 'central' is in\n");
    ve_printf("  practice imposed by the seeds rather than by us.\n");
    if (Config::mtvLikeValidation) {
      ve_printf("NUMERATOR: >= 1 reco track associated to it by 4*mccount > 3*nCandHits over\n");
      ve_printf("  ALL the track's hits, seed hits included (TrackExtra::setMCTrackIDInfo with\n");
      ve_printf("  Config::mtvLikeValidation, set by val_assoc_mtv()). That is MTV's\n");
      ve_printf("  quickTrackAssociatorByHits: shared / valid reco hits > 0.75, strictly, so\n");
      ve_printf("  3 of 4 FAILS. Fakes use the same rule: a FAKE is a track matched to no sim\n");
      ve_printf("  track at all; one matched to a non-findable (e.g. unseeded pileup) sim\n");
      ve_printf("  track is not a fake. nMinFoundHits = 0, as MTV has no length cut.\n");
    } else {
      ve_printf("NUMERATOR: >= 1 reco track associated to it by 2*mccount >= nCandHits over the\n");
      ve_printf("  non-seed hits (TrackExtra::setMCTrackIDInfo), which is what quality-val's\n");
      ve_printf("  'found tracks' counts.\n");
    }
    ve_printf("  A wrong extra hit raises the denominator of either rule only, so nothing\n");
    ve_printf("  here can be bought by taking more hits.\n");
    ve_printf("  nH >= 80%% is NOT used anywhere.\n");
    ve_printf("Reference configuration: %s\n", ref->name.c_str());

    ve_printf("\n--- the seeding ceiling (identical in every configuration by construction) ---\n");
    ve_printf("%-22s %10s %10s %9s %10s %8s %7s %7s %7s\n", "region", "sim sel", "seeded",
              "seed eff", "all seeds", "per sel", "pure", "hits", "pixel");
    for (int r = 0; r < 3; ++r) {
      const long ds = ref->den.reg[r], ss = ref->dens.reg[r], ns = ref->n_seed[r];
      ve_printf("%-22s %10ld %10ld %8.2f%% %10ld %8.2f %6.1f%% %7.2f %6.1f%%\n",
                ve_regname[r], ds, ss, ds ? 100.0*ss/ds : 0.0, ns,
                ss ? (double) ref->n_seed_on_sel[r]/ss : 0.0,
                ns ? 100.0*ref->n_seed_pure[r]/ns : 0.0,
                ns ? (double) ref->sum_seed_hits[r]/ns : 0.0,
                ref->sum_seed_hits[r] ? 100.0*ref->sum_seed_pix[r]/ref->sum_seed_hits[r] : 0.0);
    }
    ve_printf("  %-22s %10s %10s %9s %10s %9s\n", "... of those that CAN be",
              ">=4 pix lay", "seeded", "seed eff", "pre-clean", "pre eff");
    for (int r = 0; r < 3; ++r)
      ve_printf("  %-22s %10ld %10ld %8.2f%% %10ld %8.2f%%\n", ve_regname[r], ref->n_pix4[r],
                ref->n_pix4_seeded[r],
                ref->n_pix4[r] ? 100.0*ref->n_pix4_seeded[r]/ref->n_pix4[r] : 0.0,
                ref->n_pix4_preseeded[r],
                ref->n_pix4[r] ? 100.0*ref->n_pix4_preseeded[r]/ref->n_pix4[r] : 0.0);
    ve_printf("  'pre-clean' is the same vote over the file's RAW seed collection, before\n");
    ve_printf("  the iteration's seed cleaner ran -- 36042 seeds to 7498 on this sample --\n");
    ve_printf("  so pre eff minus seed eff is what CLEANING removed, and the rest is what\n");
    ve_printf("  seeding never made.\n");
    ve_printf("  The seeds are 100%% pixel hits, i.e. initialStep pixel quadruplets, so a\n");
    ve_printf("  selected sim track that does not cross 4 distinct PIXEL layers cannot be\n");
    ve_printf("  seeded by one and charging this iteration for it is charging it for what\n");
    ve_printf("  later iterations exist to do.\n");
    ve_printf("  'seed eff' is the ceiling on every efficiency below: a track with no seed\n");
    ve_printf("  cannot be found. 'all seeds' is every seed of that region, most of which\n");
    ve_printf("  are on sim tracks OUTSIDE the selection (below 0.9 GeV, mostly), so it is\n");
    ve_printf("  not a duplicate rate; 'per sel' is, being seeds per SELECTED seeded sim\n");
    ve_printf("  track. 'pure' = seeds all of whose valid hits come from one sim track\n");
    {
      const long ok = ref->n_vote_ok, tie = ref->n_vote_tie, nt = ref->n_vote_notruth;
      const long tot = ok + tie + nt;
      if (tot)
        ve_printf("  SEED -> SIM VOTE over %ld seeds: %ld resolved (%.2f%%), %ld TIED (%.2f%%),\n"
                  "  %ld with no truth-matched hit (%.2f%%). A tie returns -1 from\n"
                  "  Event::simInfoForTrack() and is NOT counted as seeding its sim track, so the\n"
                  "  tied fraction is how much 'seed eff' above is UNDERSTATED. The seed's own\n"
                  "  label() is no use for this: WriteMemoryFile sets it to seedSimIdx[is], the\n"
                  "  sim index of the rec track that used the seed, and -1 for every seed that\n"
                  "  produced none.\n",
                  tot, ok, 100.0*ok/tot, tie, 100.0*tie/tot, nt, 100.0*nt/tot);
    }
    ve_printf("  (Event::SimInfoFromHits::good_frac() == 1). 'hits' and 'pixel' say what the\n");
    ve_printf("  seeds are made of; the .bin carries the track algorithm and the hits, not\n");
    ve_printf("  the producer, so this bounds the seeding-algorithm question without\n");
    ve_printf("  settling it -- Patatrack and the standard chain differ in the FIT.\n");
    for (const auto &c : g_ve)
      for (int r = 0; r < 3; ++r)
        if (c.n_seed[r] != ref->n_seed[r])
          ve_printf("  !! %s has %ld seeds in region %d against the reference's %ld\n",
                    c.name.c_str(), c.n_seed[r], r, ref->n_seed[r]);

    ve_printf("\n--- totals ---\n");
    ve_printf("%-22s %7s %9s %9s %8s %8s %9s %9s %8s\n", "configuration", "events",
              "sim sel", "found", "eff", "of seed", "reco trks", "fakes", "dup/sim");
    for (const auto &c : g_ve) {
      ve_printf("%-22s %7ld %9ld %9ld %7.2f%% %7.2f%% %9ld %9ld %7.2f%%\n",
                c.name.c_str(), c.n_ev, c.den.tot, c.num.tot,
                c.den.tot ? 100.0*c.num.tot/c.den.tot : 0.0,
                c.dens.tot ? 100.0*c.num.tot/c.dens.tot : 0.0,
                c.reco.tot, c.fake.tot,
                c.den.tot ? 100.0*c.dup.tot/c.den.tot : 0.0);
      if (c.den.tot != ref->den.tot)
        ve_printf("   !! denominator differs from the reference by %ld -- pairing is NOT exact\n",
                  c.den.tot - ref->den.tot);
    }
    ve_printf("  'eff' is against every selected sim track (MTV); 'of seed' is against the\n");
    ve_printf("  seeded subset, i.e. what the SEARCH alone is responsible for.\n");
    for (const auto &c : g_ve) {
      const long tot = c.n_seed_found + c.n_seed_missing;
      if (tot && c.n_seed_missing * 100 > tot)   // more than 1 % unmatched is worth saying
        ve_printf("  !! %s: the seed could not be identified for %ld of %ld tracks (%.1f%%),\n"
                  "     so their association is computed over more hits than the others' --\n"
                  "     that is NOT like-for-like and flatters this configuration.\n",
                  c.name.c_str(), c.n_seed_missing, tot, 100.0*c.n_seed_missing/tot);
    }
    ve_printf("\n--- paired differences against %s, whole sample ---\n", ref->name.c_str());
    ve_printf("%-22s %14s %14s %14s\n", "configuration", "d found", "d fakes", "d duplicates");
    for (const auto &c : g_ve) {
      if (&c == ref) continue;
      double s1, e1, s2, e2, s3, e3;
      ve_paired(ve_series(c.ev_num, -1, 0),  ve_series(ref->ev_num, -1, 0),  s1, e1);
      ve_paired(ve_series(c.ev_fake, -1, 0), ve_series(ref->ev_fake, -1, 0), s2, e2);
      ve_paired(ve_series(c.ev_dup, -1, 0),  ve_series(ref->ev_dup, -1, 0),  s3, e3);
      ve_printf("%-22s %+8.0f %4.1fs %+8.0f %4.1fs %+8.0f %4.1fs\n", c.name.c_str(),
                s1, e1 > 0 ? s1/e1 : 0.0, s2, e2 > 0 ? s2/e2 : 0.0, s3, e3 > 0 ? s3/e3 : 0.0);
    }
    ve_printf("  's' is the paired significance: the sum of the per-event difference over\n"
              "  the sigma of that sum, so it is the spread of the DIFFERENCE, not Poisson.\n");

    // Fakes over EVERY reco track, no |eta| or pT cut, in mtvsum.C's ranges and
    // by the reco track's own |eta| and pT, as MTV's fake rate is.
    ve_printf("\n--- fakes over ALL reco tracks (no selection cut; reco |eta|, pT) ---\n");
    ve_printf("%-16s", "range");
    for (const auto &c : g_ve) ve_printf(" %22s", c.name.c_str());
    for (const auto &c : g_ve) if (&c != ref) ve_printf(" %17s", ("d fake " + c.name).c_str());
    ve_printf("\n");
    for (int u = 0; u < VE_NU; ++u) {
      ve_printf("%-16s", ve_uname[u]);
      for (const auto &c : g_ve) {
        long r = 0, f = 0;
        for (size_t e = 0; e < c.ev_ureco.size(); ++e) { r += c.ev_ureco[e][u]; f += c.ev_ufake[e][u]; }
        ve_printf(" %6ld/%7ld %6.2f%%", f, r, r ? 100.0 * f / r : 0.0);
      }
      for (const auto &c : g_ve) {
        if (&c == ref) continue;
        std::vector<long> a, b;
        for (auto &x : c.ev_ufake) a.push_back(x[u]);
        for (auto &x : ref->ev_ufake) b.push_back(x[u]);
        double sm, sg;
        ve_paired(a, b, sm, sg);
        ve_printf(" %+9.0f %5.1fs", sm, sg > 0 ? sm / sg : 0.0);
      }
      ve_printf("\n");
    }

    // Resolution, per region. NOT paired -- the population differs between
    // configurations by construction, since a configuration that finds more
    // tracks finds harder ones. Read it as the price of the extra tracks, and
    // read the "common" rows, which restrict every configuration to the sim
    // tracks ALL of them found, as the like-for-like comparison.
    // The COMMON subset: sim tracks every configuration found. Built as the
    // intersection over configurations of the (event, sim label) keys.
    std::map<std::pair<int,int>, int> seen;
    for (const auto &c : g_ve)
      for (const auto &x : c.res) ++seen[{x.evt, x.lbl}];
    const int ncfg = (int) g_ve.size();

    ve_printf("\n--- d(pT)/pT of the best-matched reco track, per region ---\n");
    ve_printf("ALL = every track that configuration found; COMMON = only the sim tracks\n"
              "every configuration found, which is the like-for-like comparison.\n");
    ve_printf("%-22s %-20s %-7s %8s %9s %9s %8s\n", "configuration", "region",
              "sample", "n", "median", "width", "|d|>20%");
    for (const auto &c : g_ve) {
      for (int pass = 0; pass < 2; ++pass) {
        for (int r = 0; r < 4; ++r) {
          std::vector<float> v;
          for (const auto &x : c.res) {
            if (r < 3 && x.reg != r) continue;
            if (pass == 1 && seen[{x.evt, x.lbl}] != ncfg) continue;
            v.push_back(x.r);
          }
          if (v.size() < 50) continue;
          std::sort(v.begin(), v.end());
          const size_t n = v.size();
          const double med = v[n/2];
          const double wid = 0.5 * (v[(size_t)(0.84*n)] - v[(size_t)(0.16*n)]);
          long tail = 0;
          for (float x : v) if (std::abs(x) > 0.2f) ++tail;
          ve_printf("%-22s %-20s %-7s %8zu %+9.5f %9.5f %7.2f%%\n", c.name.c_str(),
                    r == 3 ? "ALL REGIONS" : ve_regname[r], pass ? "COMMON" : "all",
                    n, med, wid, 100.0*tail/n);
        }
      }
    }

    ve_table(3, ref);
    for (int ax = 0; ax < VE_NAX; ++ax) ve_table(ax, ref);
    ve_ratio_table(3, ref, true);
    ve_ratio_table(0, ref, true);
    ve_ratio_table(3, ref, false);
    ve_ratio_table(2, ref, false);

    // TRACK LENGTH. Not paired: the denominator is each configuration's own
    // found tracks, so a configuration that finds more finds shorter ones and
    // the mean moves by composition. The sim row is the same population's truth
    // content and is what the other two should be read against.
    for (int w = 0; w < 3; ++w) ve_len_table(3, ref, w);
    for (int w = 0; w < 2; ++w) { ve_len_table(0, ref, w); ve_len_table(1, ref, w); }
    ve_len_table(0, ref, 2);  ve_len_table(1, ref, 2);

    // ---- the .root file: one efficiency TH1 per configuration per axis, with
    // binomial errors, plus an overlay canvas per axis so show-anrun's TBrowser
    // opens on something readable.
    // ---- the .root file. Two canvases per axis: the efficiency overlay, and the
    // PAIRED DIFFERENCE against the reference. The overlay alone is unreadable
    // once several configurations are in -- they sit within a few points of each
    // other on a 0-1 axis -- and the difference is the quantity with the small
    // error bar, since the denominator is shared and only the numerator moves.
    const VeCfg *ref2 = nullptr;
    for (const auto &c : g_ve) if (c.name == g_ve_ref2) ref2 = &c;

    const std::string rootf = std::string(prefix) + ".root";
    TFile f(rootf.c_str(), "RECREATE");
    static const int kCol[8] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2,
                                kMagenta + 1, kOrange + 7, kCyan + 2, kGray + 2};
    // The x axes carry alphanumeric bin labels ("1.35-1.5", "2.5-4.0"), which
    // need room that ROOT's default bottom margin does not give -- and these are
    // drawn in a deck at slide size, where a clipped label is simply wrong.
    auto ve_pad = [](TCanvas *c) {
      c->SetBottomMargin(0.16);  c->SetLeftMargin(0.11);
      c->SetRightMargin(0.04);   c->SetTopMargin(0.09);
      c->SetGridy(1);
    };
    for (int ax = 0; ax < VE_NAX; ++ax) {
      const int nb = ve_nbin[ax];
      TCanvas *cv = new TCanvas(Form("c_eff_ax%d", ax),
                                Form("efficiency vs %s", ve_axname[ax]), 900, 600);
      ve_pad(cv);
      TCanvas *cd = new TCanvas(Form("c_deff_ax%d", ax),
                                Form("efficiency difference vs %s", ve_axname[ax]), 900, 600);
      ve_pad(cd);
      TLegend *lg = new TLegend(0.60, 0.13, 0.98, 0.13 + 0.05*g_ve.size());
      TLegend *ld = new TLegend(0.60, 0.13, 0.98, 0.13 + 0.05*g_ve.size());
      int ic = 0, id = 0;
      for (const auto &c : g_ve) {
        TH1D *hn = new TH1D(Form("num_ax%d_%d", ax, ic), "", nb, -0.5, nb - 0.5);
        TH1D *hd = new TH1D(Form("den_ax%d_%d", ax, ic), "", nb, -0.5, nb - 0.5);
        for (int b = 0; b < nb; ++b) {
          hn->SetBinContent(b+1, (double) c.num.b[ax][b]);
          hd->SetBinContent(b+1, (double) c.den.b[ax][b]);
        }
        TH1D *he = (TH1D*) hn->Clone(Form("eff_ax%d_%s", ax, c.name.c_str()));
        he->SetTitle(Form("efficiency vs %s;%s;efficiency", ve_axname[ax], ve_axname[ax]));
        he->Divide(hn, hd, 1.0, 1.0, "B");
        for (int b = 0; b < nb; ++b)
          he->GetXaxis()->SetBinLabel(b+1, ve_binlabel(ax, b).c_str());
        he->SetLineColor(kCol[ic % 8]);  he->SetMarkerColor(kCol[ic % 8]);
        he->SetMarkerStyle(20 + (ic % 8));  he->SetLineWidth(2);
        he->SetMinimum(0.0);  he->SetMaximum(1.05);  he->SetStats(0);
        he->Write();
        cv->cd();  he->Draw(ic == 0 ? "E1" : "E1 SAME");
        lg->AddEntry(he, c.name.c_str(), "lp");

        // the paired difference, in points, error = sigma of the per-event sum
        if (&c != ref) {
          TH1D *hdd = new TH1D(Form("d_eff_ax%d_%s", ax, c.name.c_str()),
                               Form("efficiency difference vs %s, paired;%s;points",
                                    ve_axname[ax], ve_axname[ax]), nb, -0.5, nb - 0.5);
          for (int b = 0; b < nb; ++b) {
            const long d = ref->den.b[ax][b];
            hdd->GetXaxis()->SetBinLabel(b+1, ve_binlabel(ax, b).c_str());
            if (d < 20) continue;
            double sum, sig;
            ve_paired(ve_series(c.ev_num, ax, b), ve_series(ref->ev_num, ax, b), sum, sig);
            hdd->SetBinContent(b+1, 100.0*sum/d);
            hdd->SetBinError(b+1, 100.0*sig/d);
          }
          hdd->SetLineColor(kCol[(id+1) % 8]);  hdd->SetMarkerColor(kCol[(id+1) % 8]);
          hdd->SetMarkerStyle(20 + ((id+1) % 8));  hdd->SetLineWidth(2);  hdd->SetStats(0);
          hdd->Write();
          cd->cd();  hdd->Draw(id == 0 ? "E1" : "E1 SAME");
          ld->AddEntry(hdd, c.name.c_str(), "lp");
          ++id;
        }

        // ... and against the SECOND reference, for the panel under the plot.
        if (ref2 && &c != ref2) {
          TH1D *h2 = new TH1D(Form("d2_eff_ax%d_%s", ax, c.name.c_str()),
                              Form("efficiency minus %s, paired;%s;points",
                                   ref2->name.c_str(), ve_axname[ax]), nb, -0.5, nb - 0.5);
          for (int b = 0; b < nb; ++b) {
            const long d = ref2->den.b[ax][b];
            h2->GetXaxis()->SetBinLabel(b+1, ve_binlabel(ax, b).c_str());
            if (d < 20) continue;
            double sum, sig;
            ve_paired(ve_series(c.ev_num, ax, b), ve_series(ref2->ev_num, ax, b), sum, sig);
            h2->SetBinContent(b+1, 100.0*sum/d);
            h2->SetBinError(b+1, 100.0*sig/d);
          }
          h2->SetStats(0);  h2->Write();
        }
        delete hn;  delete hd;
        ++ic;
      }
      // THE SEEDING CEILING, on the same axes. It is identical in every
      // configuration by construction -- the seeds are the same -- so it is one
      // curve off the reference, and it is the line every efficiency below must
      // be read against: a track with no seed cannot be found. Named
      // eff_ax<N>_seeded so it is just another "configuration" to a plotter.
      {
        TH1D *hn = new TH1D(Form("sn_ax%d", ax), "", nb, -0.5, nb - 0.5);
        TH1D *hd = new TH1D(Form("sd_ax%d", ax), "", nb, -0.5, nb - 0.5);
        for (int b = 0; b < nb; ++b) {
          hn->SetBinContent(b+1, (double) ref->dens.b[ax][b]);
          hd->SetBinContent(b+1, (double) ref->den.b[ax][b]);
        }
        TH1D *he = (TH1D*) hn->Clone(Form("eff_ax%d_seeded", ax));
        he->SetTitle(Form("seeding efficiency vs %s;%s;seeded", ve_axname[ax], ve_axname[ax]));
        he->Divide(hn, hd, 1.0, 1.0, "B");
        for (int b = 0; b < nb; ++b)
          he->GetXaxis()->SetBinLabel(b+1, ve_binlabel(ax, b).c_str());
        he->SetLineColor(kGray + 2);  he->SetLineStyle(2);  he->SetLineWidth(2);
        he->SetMarkerStyle(1);  he->SetStats(0);
        he->SetMinimum(0.0);  he->SetMaximum(1.05);
        he->Write();
        cv->cd();  he->Draw("HIST SAME");
        lg->AddEntry(he, "has a seed", "l");
        delete hn;  delete hd;
      }

      cv->cd();  lg->Draw();  cv->Write();
      cd->cd();  ld->Draw();  cd->Write();
    }

    // Mean TRUTH-MATCHED hits per found track. The error is the error on the
    // mean, sqrt((<x^2> - <x>^2)/n), which needs the sum of squares -- without
    // it the spread of track lengths would be mistaken for a measurement error.
    // The sim track's own content goes on the same axes as a dashed reference:
    // these are meaningless read against 100 % and only mean something read
    // against what the track actually left behind.
    for (int ax = 0; ax < 2; ++ax) {
      const int nb = ve_nbin[ax];
      TCanvas *cl = new TCanvas(Form("c_len_ax%d", ax),
                                Form("matched hits per track vs %s", ve_axname[ax]), 900, 600);
      ve_pad(cl);
      TLegend *ll = new TLegend(0.60, 0.13, 0.98, 0.13 + 0.05*(g_ve.size() + 1));
      int ic = 0;
      for (const auto &c : g_ve) {
        TH1D *h = new TH1D(Form("len_ax%d_%s", ax, c.name.c_str()),
                           Form("matched hits per found track vs %s;%s;matched hits",
                                ve_axname[ax], ve_axname[ax]), nb, -0.5, nb - 0.5);
        for (int b = 0; b < nb; ++b)
          h->GetXaxis()->SetBinLabel(b+1, ve_binlabel(ax, b).c_str());
        ve_mean_hist(h, ax, c.n_match, c.n_match2, c.num, nb);
        h->SetLineColor(kCol[ic % 8]);  h->SetMarkerColor(kCol[ic % 8]);
        h->SetMarkerStyle(20 + (ic % 8));  h->SetLineWidth(2);  h->SetStats(0);
        h->SetMinimum(0.0);
        h->Write();
        cl->cd();  h->Draw(ic == 0 ? "E1" : "E1 SAME");
        ll->AddEntry(h, c.name.c_str(), "lp");

        // Delta against the second reference. NOT paired -- the per-event sums
        // for a mean are not kept -- so the error is the two errors in
        // quadrature, which over-states it where the populations overlap, and
        // they overlap heavily here. Read it as an upper bound on the error.
        if (ref2 && &c != ref2) {
          TH1D *h2 = new TH1D(Form("d2_len_ax%d_%s", ax, c.name.c_str()),
                              Form("matched hits minus %s;%s;hits",
                                   ref2->name.c_str(), ve_axname[ax]), nb, -0.5, nb - 0.5);
          TH1D *hr = new TH1D(Form("tmp_ref_ax%d_%d", ax, ic), "", nb, -0.5, nb - 0.5);
          ve_mean_hist(hr, ax, ref2->n_match, ref2->n_match2, ref2->num, nb);
          for (int b = 0; b < nb; ++b) {
            h2->GetXaxis()->SetBinLabel(b+1, ve_binlabel(ax, b).c_str());
            if (h->GetBinContent(b+1) == 0 || hr->GetBinContent(b+1) == 0) continue;
            h2->SetBinContent(b+1, h->GetBinContent(b+1) - hr->GetBinContent(b+1));
            h2->SetBinError(b+1, std::hypot(h->GetBinError(b+1), hr->GetBinError(b+1)));
          }
          h2->SetStats(0);  h2->Write();  delete hr;
        }
        ++ic;
      }
      // Two reference curves, both with errors on the mean: what the sim track
      // left behind, and how much of a found track its SEED already was --
      // matched hits EXCLUDE seed hits, so the seed size says what the search
      // was starting from.
      TH1D *hs = new TH1D(Form("len_ax%d_simref", ax),
                          "sim track's own hits", nb, -0.5, nb - 0.5);
      for (int b = 0; b < nb; ++b) hs->GetXaxis()->SetBinLabel(b+1, ve_binlabel(ax, b).c_str());
      ve_mean_hist(hs, ax, ref->n_sim, ref->n_sim2, ref->num, nb);
      hs->SetLineColor(kGray + 2);  hs->SetLineStyle(2);  hs->SetLineWidth(2);
      hs->SetMarkerStyle(1);  hs->SetStats(0);
      hs->Write();
      cl->cd();  hs->Draw("HIST SAME");
      ll->AddEntry(hs, "sim track has", "l");

      TH1D *hq = new TH1D(Form("len_ax%d_seedref", ax),
                          "hits in the seed", nb, -0.5, nb - 0.5);
      for (int b = 0; b < nb; ++b) hq->GetXaxis()->SetBinLabel(b+1, ve_binlabel(ax, b).c_str());
      ve_mean_hist(hq, ax, ref->n_seedhits, ref->n_seedhits2, ref->num, nb);
      hq->SetLineColor(kGray + 1);  hq->SetLineStyle(3);  hq->SetLineWidth(2);
      hq->SetMarkerStyle(1);  hq->SetStats(0);
      hq->Write();
      cl->cd();  hq->Draw("HIST SAME");
      ll->AddEntry(hq, "seed has", "l");
      ll->Draw();  cl->Write();
    }
    f.Close();

    ve_printf("\nval_eff_report: wrote %s and %s (%zu configurations)\n",
              rootf.c_str(), txt.c_str(), g_ve.size());
    if (g_ve_log) { fclose(g_ve_log); g_ve_log = nullptr; }
  }

  // ==========================================================================
  // val_chopres -- the pT5 pixel-chop recovery, RESOLVED.
  //
  // Same measurement as val_chop_recovery_* and the same reason for preferring
  // it on this sample: the chopped hits were found by the upstream
  // reconstruction, so the comparison is an exact (layer, index) match and NO
  // TRUTH IS INVOLVED. That matters here beyond the usual mc_match caution --
  // the HLT March sample predates the split-cluster arbitration fix, so its
  // rec->sim links are stale and any truth-matched efficiency on it is biased.
  // This metric is immune to that.
  //
  // Resolved in |eta| and pT of the candidate itself, and in chopped hits per
  // chopped LAYER, which is the axis the one-hit-per-layer ceiling lives on.
  namespace {
    struct ChopCfg {
      std::string name;
      VeBins hden, hnum;             // chopped hits, recovered hits
      VeBins tden, tnum;             // tracks with chopped hits, fully recovered
      std::vector<VeBins> ev_hnum, ev_tnum;
      long n_ev = 0;
    };
    std::vector<ChopCfg> g_cr;
    std::string g_cr_ref;

    ChopCfg &cr_cfg(const char *name) {
      for (auto &c : g_cr) if (c.name == name) return c;
      g_cr.push_back(ChopCfg());  g_cr.back().name = name;  return g_cr.back();
    }
    void ve_fill_n(VeBins &v, int be, int bp, int bh, int reg, long n) {  // chopres: no MTV gate
      if (be >= 0) v.b[0][be] += n;
      if (bp >= 0) v.b[1][bp] += n;
      if (bh >= 0) v.b[2][bh] += n;
      if (reg >= 0) v.reg[reg] += n;
      v.tot += n;
    }
  }

  void val_chopres_reset() { g_cr.clear(); g_cr_ref.clear(); }
  void val_chopres_ref(const char *cfg) { g_cr_ref = cfg; }

  void val_chopres_event(const Event *ev, const char *cfg) {
    if (ev == nullptr) return;
    ChopCfg &C = cr_cfg(cfg);
    VeBins e_hnum, e_tnum;
    for (const Track &c : ev->candidateTracks_) {
      auto it = Shell::s_chopped_hits.find(c.label());
      if (it == Shell::s_chopped_hits.end() || it->second.empty()) continue;
      const auto &chopped = it->second;
      std::set<int> clay;
      for (const HitOnTrack &ch : chopped) clay.insert(ch.layer);
      const float ae = std::abs(c.momEta());
      const int be = ve_bin_eta(ae), bp = ve_bin_pt(c.pT()), reg = ve_region(ae);
      const int bh = ve_bin_hpl((float) chopped.size() / (float) std::max<size_t>(1, clay.size()));
      int back = 0;
      for (const HitOnTrack &ch : chopped) {
        for (int i = 0; i < c.nTotalHits(); ++i) {
          const HitOnTrack hot = c.getHitOnTrack(i);
          if (hot.layer == ch.layer && hot.index == ch.index) { ++back; break; }
        }
      }
      ve_fill_n(C.hden, be, bp, bh, reg, (long) chopped.size());
      ve_fill_n(e_hnum, be, bp, bh, reg, back);
      ve_fill(C.tden, be, bp, bh, reg, 0.0f, 1e9f);   // chopres is not MTV-gated
      if (back == (int) chopped.size()) ve_fill(e_tnum, be, bp, bh, reg, 0.0f, 1e9f);
    }
    ve_add(C.hnum, e_hnum);  ve_add(C.tnum, e_tnum);
    C.ev_hnum.push_back(e_hnum);  C.ev_tnum.push_back(e_tnum);
    ++C.n_ev;
    Shell::s_chopped_hits.clear();
  }

  namespace {
    void cr_table(int ax, const ChopCfg *ref, bool track_level) {
      const int nb = (ax == 3) ? 3 : ve_nbin[ax];
      ve_printf("\n--- %s vs %s ---\n",
                track_level ? "tracks FULLY recovered" : "chopped hits recovered",
                ax == 3 ? "region" : (ax == 2 ? "chopped hits / chopped layer" : ve_axname[ax]));
      ve_printf("%-20s %9s", ax == 3 ? "region" : "bin", track_level ? "tracks" : "hits");
      for (const auto &c : g_cr) ve_printf(" | %-10.10s", c.name.c_str());
      ve_printf("\n");
      for (int b = 0; b < nb; ++b) {
        const VeBins &D = track_level ? ref->tden : ref->hden;
        const long den = (ax == 3) ? D.reg[b] : D.b[ax][b];
        if (den < 20) continue;
        std::string lab = (ax == 3) ? ve_regname[b] : ve_binlabel(ax, b);
        ve_printf("%-20s %9ld", lab.c_str(), den);
        for (const auto &c : g_cr) {
          const VeBins &Dc = track_level ? c.tden : c.hden;
          const VeBins &Nc = track_level ? c.tnum : c.hnum;
          const long d = (ax == 3) ? Dc.reg[b] : Dc.b[ax][b];
          const long n = (ax == 3) ? Nc.reg[b] : Nc.b[ax][b];
          if (&c == ref) ve_printf(" |   %6.2f%% ", d ? 100.0*n/d : 0.0);
          else {
            double sum, sig;
            ve_paired(ve_series(track_level ? c.ev_tnum : c.ev_hnum, ax == 3 ? 3 : ax, b),
                      ve_series(track_level ? ref->ev_tnum : ref->ev_hnum, ax == 3 ? 3 : ax, b),
                      sum, sig);
            const double dp = d ? 100.0*sum/d : 0.0;
            const double sp = d ? 100.0*sig/d : 0.0;
            ve_printf(" | %+6.2f%s%-3.3s", dp,
                      sp > 0 && std::abs(dp) > 3*sp ? "*" : " ",
                      sp > 0 && std::abs(dp) > 3*sp ? "sig" : "");
          }
        }
        ve_printf("\n");
      }
    }
  }

  void val_chopres_report(const char *prefix) {
    if (g_cr.empty()) { printf("val_chopres_report: nothing accumulated.\n"); return; }
    const ChopCfg *ref = &g_cr[0];
    if (!g_cr_ref.empty())
      for (const auto &c : g_cr) if (c.name == g_cr_ref) ref = &c;

    const std::string txt = std::string(prefix) + ".txt";
    g_ve_log = fopen(txt.c_str(), "w");

    ve_printf("\n================================================================\n");
    ve_printf("  val_chopres -- pT5 pixel-chop recovery, resolved\n");
    ve_printf("================================================================\n");
    ve_printf("Truth-FREE: exact (layer, index) match against the hits the chop removed.\n");
    ve_printf("Reference configuration: %s\n", ref->name.c_str());
    ve_printf("\n--- totals ---\n");
    ve_printf("%-22s %7s %9s %9s %8s %9s %9s %8s\n", "configuration", "events",
              "hits", "recovered", "frac", "tracks", "full", "frac");
    for (const auto &c : g_cr)
      ve_printf("%-22s %7ld %9ld %9ld %7.2f%% %9ld %9ld %7.2f%%\n",
                c.name.c_str(), c.n_ev, c.hden.tot, c.hnum.tot,
                c.hden.tot ? 100.0*c.hnum.tot/c.hden.tot : 0.0,
                c.tden.tot, c.tnum.tot, c.tden.tot ? 100.0*c.tnum.tot/c.tden.tot : 0.0);
    for (const auto &c : g_cr)
      if (c.hden.tot != ref->hden.tot)
        ve_printf("  !! %s: chopped-hit denominator differs from the reference by %ld\n",
                  c.name.c_str(), c.hden.tot - ref->hden.tot);

    cr_table(3, ref, false);  cr_table(0, ref, false);
    cr_table(1, ref, false);  cr_table(2, ref, false);
    cr_table(3, ref, true);   cr_table(2, ref, true);
    ve_printf("  reference column is absolute; the others are the PAIRED difference in\n"
              "  points, * = |delta| > 3 sigma of the per-event spread.\n");

    const std::string rootf = std::string(prefix) + ".root";
    TFile f(rootf.c_str(), "RECREATE");
    for (int ax = 0; ax < VE_NAX; ++ax) {
      const int nb = ve_nbin[ax];
      TCanvas *cv = new TCanvas(Form("c_chop_ax%d", ax), Form("chop recovery vs axis %d", ax), 900, 600);
      TLegend *lg = new TLegend(0.60, 0.15, 0.98, 0.15 + 0.05*g_cr.size());
      int ic = 0;
      for (const auto &c : g_cr) {
        TH1D *hn = new TH1D(Form("cnum_ax%d_%d", ax, ic), "", nb, -0.5, nb - 0.5);
        TH1D *hd = new TH1D(Form("cden_ax%d_%d", ax, ic), "", nb, -0.5, nb - 0.5);
        for (int b = 0; b < nb; ++b) {
          hn->SetBinContent(b+1, (double) c.hnum.b[ax][b]);
          hd->SetBinContent(b+1, (double) c.hden.b[ax][b]);
        }
        TH1D *he = (TH1D*) hn->Clone(Form("chop_ax%d_%s", ax, c.name.c_str()));
        he->SetTitle(Form("chopped hits recovered vs %s;%s;recovered",
                          ax == 2 ? "chopped hits / layer" : ve_axname[ax],
                          ax == 2 ? "chopped hits / layer" : ve_axname[ax]));
        he->Divide(hn, hd, 1.0, 1.0, "B");
        for (int b = 0; b < nb; ++b) he->GetXaxis()->SetBinLabel(b+1, ve_binlabel(ax, b).c_str());
        he->SetLineColor(1 + ic);  he->SetMarkerColor(1 + ic);  he->SetMarkerStyle(20 + ic);
        he->SetMinimum(0.0);  he->SetMaximum(1.05);
        he->Write();
        cv->cd();  he->Draw(ic == 0 ? "E1" : "E1 SAME");
        lg->AddEntry(he, c.name.c_str(), "lp");
        delete hn;  delete hd;  ++ic;
      }
      lg->Draw();  cv->Write();
    }
    f.Close();
    ve_printf("\nval_chopres_report: wrote %s and %s (%zu configurations)\n",
              rootf.c_str(), txt.c_str(), g_cr.size());
    if (g_ve_log) { fclose(g_ve_log); g_ve_log = nullptr; }
  }


  //--------------------------------------------------------------------------
  // val_trkdump -- one row per (configuration, event, sim track): the best-
  // matched reco track (most matched hits, as val_eff), its pT, its state's
  // position and its hits in HoT order (layer << 20 | index, -1 - layer for a
  // hole). For comparing two configurations track by track on the same sim
  // tracks: same hits and a different pT is the state, different hits the
  // selection.
  //--------------------------------------------------------------------------
  namespace {
    struct TdRow {
      int cfg, ev, sim, nfound;
      float sim_pt, sim_eta, pt, eta, chi2, sx, sy, sz;
      std::vector<int> hits;
      std::vector<float> hr, hz;
      std::vector<int> hown;
    };
    std::vector<std::string> g_td_cfg;
    std::vector<TdRow> g_td;
  }  // namespace

  void val_trkdump_reset() { g_td_cfg.clear(); g_td.clear(); }

  void val_trkdump_event(const Event *ev, const char *cfg) {
    if (ev == nullptr) return;
    int ci = std::find(g_td_cfg.begin(), g_td_cfg.end(), std::string(cfg)) - g_td_cfg.begin();
    if (ci == (int) g_td_cfg.size()) g_td_cfg.push_back(cfg);
    std::map<int, std::pair<int, const Track *>> best;
    for (const Track &c : ev->candidateTracks_) {
      TrackExtra extra(c.label());
      extra.setMCTrackIDInfo(c, ev->layerHits_, ev->simHitsInfo_, ev->simTracks_, false, false);
      const int mc = extra.mcTrackID();
      if (mc < 0 || mc >= (int) ev->simTracks_.size()) continue;
      auto &b = best[mc];
      if (b.second == nullptr || extra.nHitsMatched() > b.first) b = {extra.nHitsMatched(), &c};
    }
    for (auto &[mc, b] : best) {
      const Track &c = *b.second, &st = ev->simTracks_[mc];
      TdRow r{ci, ev->evtID(), mc, c.nFoundHits(), st.pT(), st.momEta(), c.pT(), c.momEta(), c.chi2(),
              c.x(), c.y(), c.z(), {}, {}, {}, {}};
      for (int h = 0; h < c.nTotalHits(); ++h) {
        const HitOnTrack hot = c.getHitOnTrack(h);
        r.hits.push_back(hot.index >= 0 ? (hot.layer << 20 | hot.index) : -1 - hot.layer);
        if (hot.index >= 0) {
          const Hit &hh = ev->layerHits_[hot.layer][hot.index];
          r.hr.push_back(std::hypot(hh.x(), hh.y()));
          r.hz.push_back(hh.z());
          const int id = hh.mcHitID();
          const int st = (id >= 0 && id < (int) ev->simHitsInfo_.size()) ? ev->simHitsInfo_[id].mcTrackID() : -1;
          r.hown.push_back(st == mc ? 1 : (st < 0 ? -1 : 0));
        } else {
          r.hr.push_back(-1);
          r.hz.push_back(0);
          r.hown.push_back(-1);
        }
      }
      g_td.push_back(std::move(r));
    }
  }

  void val_trkdump_write(const char *fn) {
    TFile f(fn, "RECREATE");
    TTree t("trk", "best-matched reco track per sim track");
    TdRow r;
    std::vector<int> *hp = &r.hits;
    std::vector<float> *hrp = &r.hr, *hzp = &r.hz;
    std::vector<int> *hop = &r.hown;
    t.Branch("cfg", &r.cfg);  t.Branch("ev", &r.ev);  t.Branch("sim", &r.sim);
    t.Branch("nfound", &r.nfound);  t.Branch("sim_pt", &r.sim_pt);  t.Branch("sim_eta", &r.sim_eta);
    t.Branch("pt", &r.pt);  t.Branch("eta", &r.eta);  t.Branch("chi2", &r.chi2);
    t.Branch("sx", &r.sx);  t.Branch("sy", &r.sy);  t.Branch("sz", &r.sz);
    t.Branch("hits", &hp);
    t.Branch("hr", &hrp);
    t.Branch("hz", &hzp);
    t.Branch("hown", &hop);
    for (auto &x : g_td) { r = x; t.Fill(); }
    t.Write();
    TNamed n("cfgs", "");
    std::string cs;
    for (auto &c : g_td_cfg) cs += c + " ";
    n.SetTitle(cs.c_str());
    n.Write();
    printf("val_trkdump_write: %zu rows, configs %s-> %s\n", g_td.size(), cs.c_str(), fn);
  }

}  // namespace mkfit
