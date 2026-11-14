#include "RecoTracker/MkFitCMS/standalone/Probe/Probe.h"

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace mkfit {

  //============================================================================
  // Geometry compatibility of a sample -- see the header.
  //============================================================================

  namespace {
    std::vector<double> g_gc_d;
    long g_gc_bad_sid = 0, g_gc_hits = 0, g_gc_evt = 0;
  }

  void probe_geom_check_reset() {
    g_gc_d.clear(); g_gc_bad_sid = 0; g_gc_hits = 0; g_gc_evt = 0;
  }

  void probe_geom_check_event(const Event *ev) {
    if (ev == nullptr) { printf("probe_geom_check_event: null event\n"); return; }
    const TrackerInfo &ti = Config::TrkInfo;
    ++g_gc_evt;
    for (int lay = 0; lay < (int)ev->layerHits_.size(); ++lay) {
      if (lay >= ti.n_layers()) break;
      const LayerInfo &li = ti.layer(lay);
      for (const Hit &h : ev->layerHits_[lay]) {
        ++g_gc_hits;
        const unsigned int sid = h.detIDinLayer();
        if ((int)sid >= li.n_modules()) { ++g_gc_bad_sid; continue; }
        const ModuleInfo &mi = li.module_info(sid);
        const double d = std::fabs((h.x() - mi.pos[0]) * mi.zdir[0] +
                                   (h.y() - mi.pos[1]) * mi.zdir[1] +
                                   (h.z() - mi.pos[2]) * mi.zdir[2]);
        g_gc_d.push_back(d);
      }
    }
  }

  void probe_geom_check_report(const char *sample_name) {
    printf("\n================ geometry compatibility ================\n");
    printf("sample : %s\n", sample_name);
    printf("events : %ld,  hits : %ld\n", g_gc_evt, g_gc_hits);
    printf("module short_id out of range : %ld  (%.4f %%)\n",
           g_gc_bad_sid, g_gc_hits ? 100.0 * g_gc_bad_sid / g_gc_hits : 0.0);
    if (g_gc_d.empty()) { printf("no usable hits\n"); return; }
    std::sort(g_gc_d.begin(), g_gc_d.end());
    auto q = [&](double p) { return g_gc_d[(size_t)(0.01 * p * (g_gc_d.size() - 1))]; };
    printf("|n.(hit - module plane point)|, cm -- a hit is ON its module by construction\n");
    printf("  median %10.4g   p90 %10.4g   p99 %10.4g   max %10.4g\n",
           q(50), q(90), q(99), g_gc_d.back());
    long over = 0; for (double x : g_gc_d) if (x > 0.01) ++over;
    printf("  fraction over 100 um : %.4f %%\n", 100.0 * over / g_gc_d.size());
    printf("\nreference: 2.9e-9 cm median = COMPATIBLE;  0.26 cm median = NOT\n");
    printf("=======================================================\n");
  }

  void probe_sample_info(const Event *ev) {
    if (ev == nullptr) { printf("probe_sample_info: null event\n"); return; }
    auto describe = [](const char *nm, const TrackVec &tv) {
      if (tv.empty()) { printf("  %-18s empty\n", nm); return; }
      std::map<int, int> algo;
      double ptmin = 1e9, ptmax = -1e9;
      long nh = 0, nlab = 0;
      for (const Track &t : tv) {
        algo[t.algoint()]++;
        ptmin = std::min(ptmin, (double)t.pT());
        ptmax = std::max(ptmax, (double)t.pT());
        nh += t.nFoundHits();
        if (t.label() >= 0) ++nlab;
      }
      printf("  %-18s n=%-7zu  pT %.3g .. %.4g  <nhits> %.1f  labelled %.0f%%  algos:",
             nm, tv.size(), ptmin, ptmax, (double)nh / tv.size(), 100.0 * nlab / tv.size());
      for (auto &a : algo) printf(" %d:%d", a.first, a.second);
      printf("\n");
    };
    printf("\n---- sample contents, one event ----\n");
    printf("  %-18s %zu layers\n", "layerHits_", ev->layerHits_.size());
    long nh = 0; for (auto &v : ev->layerHits_) nh += v.size();
    printf("  %-18s %ld hits\n", "total hits", nh);
    printf("  %-18s %zu\n", "simHitsInfo_", ev->simHitsInfo_.size());
    printf("  %-18s %zu   <-- if this equals the number of SIM HITS there is a\n",
           "simTrackStates_", ev->simTrackStates_.size());
    printf("  %-18s     true state at every hit, and the bias measurement needs no fit\n", "");
    {  // how many hits actually carry usable MC info, and do labels line up?
      long with_mc = 0, foreign = 0, checked = 0;
      for (int lay = 0; lay < (int)ev->layerHits_.size(); ++lay)
        for (const Hit &h : ev->layerHits_[lay]) {
          const int mid = h.mcHitID();
          if (mid >= 0 && mid < (int)ev->simHitsInfo_.size()) ++with_mc;
        }
      for (const Track &t : ev->simTracks_) {
        for (int ih = 0; ih < t.nTotalHits(); ++ih) {
          const int idx = t.getHitIdx(ih), lay = t.getHitLyr(ih);
          if (idx < 0 || lay < 0 || lay >= (int)ev->layerHits_.size()) continue;
          if (idx >= (int)ev->layerHits_[lay].size()) continue;
          const Hit &h = ev->layerHits_[lay][idx];
          const int mid = h.mcHitID();
          if (mid < 0 || mid >= (int)ev->simHitsInfo_.size()) continue;
          ++checked;
          if (ev->simHitsInfo_[mid].mcTrackID() != t.label()) ++foreign;
        }
      }
      printf("  %-18s %ld of %ld hits carry MC info\n", "mc coverage", with_mc, nh);
      printf("  %-18s %ld of %ld hits ON a sim track belong to a DIFFERENT particle"
             " (%.2f%%)\n", "foreign hits", foreign, checked,
             checked ? 100.0 * foreign / checked : 0.0);
    }
    describe("simTracks_", ev->simTracks_);
    describe("seedTracks_", ev->seedTracks_);
    describe("cmsswTracks_", ev->cmsswTracks_);
  }

}  // namespace mkfit
