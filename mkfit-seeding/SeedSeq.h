#ifndef mkfit_seeding_SeedSeq_h
#define mkfit_seeding_SeedSeq_h

// The pixel-layer sequence of every sim track, for choosing seeding patterns
// beyond the barrel.  `seedfind --sequences OUT.txt`.
//
// For each sim track with pT > pT_min, production point within D0_max of the
// beam spot in xy and |eta| < 4, one line
//
//   T ev label eta pt nA seqA nB seqB
//
// seqA: the pixel layers (mkFit 0-3, 16-27, 38-49) holding a rec hit whose
//       truth label is the track (rec -> sim, as SeedTruth's findable), in the
//       order the track crosses them;
// seqB: the same from the sim track's own hit list (sim -> rec), which the
//       bestTkIdx arbitration does not gate, but which also carries hits of
//       other particles (delta rays, neighbours) attached to the sim track.
//
// Order: each layer is placed at its hit nearest the production point, and
// layers are sorted by that distance.  At pT > 0.5 GeV no track turns round
// inside the pixel volume (2 R_c > 90 cm against r < 26 cm), so this is the
// crossing order.  A sequence is written as comma-separated mkFit layer
// numbers, "-" if empty.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

namespace mkfit::seeding {

  struct SeedSeq {
    FILE *f = nullptr;
    double eta_max = 4.0;

    static bool is_pixel(int l) { return l <= 3 || (l >= 16 && l <= 27) || (l >= 38 && l <= 49); }
    // pixel only, or every layer (--sequences-all): the outer tracker too
    bool all_layers = false;
    bool recorded(int l) const { return all_layers || is_pixel(l); }

    bool open(const std::string &fn, float pt_min, float d0_max) {
      f = fopen(fn.c_str(), "w");
      if (!f)
        return false;
      fprintf(f, "# seedfind --sequences; pt_min %.3f GeV, d0_max %.3f cm; %s\n", pt_min, d0_max, all_layers ? "all layers" : "pixel layers only");
      fprintf(f, "# T ev label eta pt nA seqA nB seqB   (A: rec hits labelled with the track, B: the sim track's hit list)\n");
      return true;
    }

    void close() {
      if (f)
        fclose(f);
      f = nullptr;
    }

    // per layer, the smallest distance of a hit to the production point; then
    // the layers sorted by it
    static std::string sequence(std::unordered_map<int, double> &dmin, int &n) {
      std::vector<std::pair<double, int>> v;
      for (const auto &kv : dmin)
        v.push_back({kv.second, kv.first});
      std::sort(v.begin(), v.end());
      n = v.size();
      if (v.empty())
        return "-";
      std::string s;
      for (const auto &p : v) {
        if (!s.empty())
          s += ',';
        s += std::to_string(p.second);
      }
      return s;
    }

    void event(int iev, const Event &ev, float pt_min, float d0_max) {
      if (!f)
        return;
      const MCHitInfoVec &mc = ev.simHitsInfo_;
      const int nsim = ev.simTracks_.size();
      auto keep = [&](int l) {
        if (l < 0 || l >= nsim)
          return false;
        const Track &t = ev.simTracks_[l];
        return t.pT() >= pt_min && std::hypot(t.x() - ev.beamSpot_.x, t.y() - ev.beamSpot_.y) <= d0_max &&
               std::abs(t.momEta()) < eta_max;
      };
      auto dist = [&](const Track &t, const Hit &h) {
        return std::sqrt((h.x() - t.x()) * (h.x() - t.x()) + (h.y() - t.y()) * (h.y() - t.y()) +
                         (h.z() - t.z()) * (h.z() - t.z()));
      };
      auto upd = [](std::unordered_map<int, double> &m, int l, double d) {
        auto it = m.find(l);
        if (it == m.end() || d < it->second)
          m[l] = d;
      };
      // A: rec -> sim
      std::unordered_map<int, std::unordered_map<int, double>> A;
      for (int l = 0; l < (int)ev.layerHits_.size(); ++l) {
        if (!recorded(l))
          continue;
        for (const Hit &h : ev.layerHits_[l]) {
          const int mid = h.mcHitID();
          const int lab = (mid >= 0 && mid < (int)mc.size()) ? mc[mid].mcTrackID() : -1;
          if (keep(lab))
            upd(A[lab], l, dist(ev.simTracks_[lab], h));
        }
      }
      for (int lab = 0; lab < nsim; ++lab) {
        if (!keep(lab))
          continue;
        const Track &t = ev.simTracks_[lab];
        // B: sim -> rec
        std::unordered_map<int, double> B;
        for (int i = 0; i < t.nTotalHits(); ++i) {
          const int l = t.getHitLyr(i), k = t.getHitIdx(i);
          if (l < 0 || k < 0 || !recorded(l) || k >= (int)ev.layerHits_[l].size())
            continue;
          upd(B, l, dist(t, ev.layerHits_[l][k]));
        }
        auto ia = A.find(lab);
        std::unordered_map<int, double> none;
        int na = 0, nb = 0;
        const std::string sa = sequence(ia == A.end() ? none : ia->second, na);
        const std::string sb = sequence(B, nb);
        fprintf(f, "T %d %d %.4f %.3f %d %s %d %s\n", iev, lab, t.momEta(), t.pT(), na, sa.c_str(), nb, sb.c_str());
      }
    }
  };

}  // namespace mkfit::seeding

#endif
