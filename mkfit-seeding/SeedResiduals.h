#ifndef mkfit_seeding_SeedResiduals_h
#define mkfit_seeding_SeedResiduals_h

// Residuals of the fixed windows on true quads: `seedfind --residuals OUT.txt`.
//
// For every findable track (SeedTruth's definition), found or not, evaluate
// every cut on the track's OWN four hits with eval_quad<A>(), as SeedMissed
// does, and write one row per track: the truth kinematics, the triplet
// circle's estimate of the curvature and the r-z slope, and the |residual| of
// the three fixed windows.  This is what sizing a window from the track's own
// curvature needs: the residual distribution of true quads against 1/p.
//
// With more than one hit in a layer, the combination kept is the one whose
// first failing cut comes latest in the finder's order; among equals, the one
// with the smallest sum of |residual| / window over c_z, d_phi and d_z.

#include "SeedMissed.h"

#include <cmath>
#include <cstdio>

namespace mkfit::seeding {

  struct SeedResiduals {
    FILE *out = nullptr;
    long n_rows = 0;

    void open(const std::string &fname, const SeedParams &P) {
      out = fopen(fname.c_str(), "w");
      if (!out) {
        fprintf(stderr, "SeedResiduals: cannot open %s\n", fname.c_str());
        return;
      }
      fprintf(out, "# seedfind --residuals: one row per findable track, its own hits\n");
      fprintf(out, "# windows at this run: c_z %.4f cm, d_phi %.5f rad, d_z %.4f cm; phi_lin %d %.4f\n", P.qwin,
              P.phiwin_d, P.qwin_d, P.phi_lin, P.phi_lin_marg);
      fprintf(out, "# stage: first failing cut, index into SeedMissed::kOrder (12 = none fails)\n");
      fprintf(out, "# pt_est = 0.0114 * R [GeV, cm]; residuals are |.|, nan where the cut was not evaluated\n");
      fprintf(out, "# ev label pt p eta found stage R cot c_z d_phi d_z\n");
    }

    template <class A, typename L>
    void event(int iev, const Event &ev, int la, int lb, int lc, int ld, const SeedParams &P, const L &ga,
               const L &gb, const L &gc, const L &gd, const std::vector<Quad> &quads) {
      if (!out)
        return;
      const MCHitInfoVec &mc = ev.simHitsInfo_;
      const HitVec *H[4] = {&ev.layerHits_[la], &ev.layerHits_[lb], &ev.layerHits_[lc], &ev.layerHits_[ld]};
      auto label = [&mc](const Hit &h) -> int {
        const int mid = h.mcHitID();
        return (mid >= 0 && mid < (int)mc.size()) ? mc[mid].mcTrackID() : -1;
      };
      std::unordered_map<int, std::array<std::vector<unsigned int>, 4>> hits;
      for (int i = 0; i < 4; ++i)
        for (unsigned int k = 0; k < H[i]->size(); ++k) {
          const int l = label((*H[i])[k]);
          if (l >= 0)
            hits[l][i].push_back(k);
        }
      std::unordered_map<int, bool> found;
      for (const auto &q : quads) {
        const int l0 = label((*H[0])[q[0]]);
        if (l0 >= 0 && l0 == label((*H[1])[q[1]]) && l0 == label((*H[2])[q[2]]) && l0 == label((*H[3])[q[3]]))
          found[l0] = true;
      }
      std::vector<unsigned int> io[4];
      const L *G[4] = {&ga, &gb, &gc, &gd};
      for (int i = 0; i < 4; ++i) {
        io[i].resize(G[i]->n());
        for (unsigned int k = 0; k < G[i]->n(); ++k)
          io[i][G[i]->orig_[k]] = k;
      }
      auto res = [](const CutVal &cv, double w) { return cv.valid ? w - cv.m : std::nan(""); };
      for (const auto &kv : hits) {
        const int l = kv.first;
        const auto &hl = kv.second;
        if (hl[0].empty() || hl[1].empty() || hl[2].empty() || hl[3].empty() || l >= (int)ev.simTracks_.size())
          continue;
        const Track &t = ev.simTracks_[l];
        if (t.pT() < P.pt_min || std::hypot(t.x() - ev.beamSpot_.x, t.y() - ev.beamSpot_.y) > P.d0_max ||
            std::abs(t.momEta()) > SeedTruth::kEtaMax)
          continue;
        int best_stage = -1;
        double best_sum = 1e30;
        QuadEval best;
        for (unsigned int a : hl[0])
          for (unsigned int b : hl[1])
            for (unsigned int c : hl[2])
              for (unsigned int d : hl[3]) {
                QuadEval e;
                eval_quad<A>(P, ga, gb, gc, gd, io[0][a], io[1][b], io[2][c], io[3][d], e);
                int stage = SeedMissed::kNo;
                for (int s = 0; s < SeedMissed::kNo; ++s) {
                  const CutVal &cv = e.c[SeedMissed::kOrder[s]];
                  if (cv.valid && !cv.pass) {
                    stage = s;
                    break;
                  }
                }
                double sum = 0;
                const double rz = res(e.c[MC_c_z], P.qwin), rp = res(e.c[MC_d_phi], P.phiwin_d),
                             rd = res(e.c[MC_d_z], P.qwin_d);
                sum += std::isnan(rz) ? 1e3 : rz / P.qwin;
                sum += std::isnan(rp) ? 1e3 : rp / P.phiwin_d;
                sum += std::isnan(rd) ? 1e3 : rd / P.qwin_d;
                if (stage > best_stage || (stage == best_stage && sum < best_sum)) {
                  best_stage = stage;
                  best_sum = sum;
                  best = e;
                }
              }
        fprintf(out, "%d %d %.5g %.5g %.4f %d %d %.6g %.6g %.5g %.5g %.5g\n", iev, l, t.pT(), t.p(), t.momEta(),
                found.count(l) ? 1 : 0, best_stage, best.R, best.cot, res(best.c[MC_c_z], P.qwin),
                res(best.c[MC_d_phi], P.phiwin_d), res(best.c[MC_d_z], P.qwin_d));
        ++n_rows;
      }
    }

    void close() {
      if (out) {
        fclose(out);
        printf("   residuals: %ld rows written\n", n_rows);
      }
    }
  };

}  // namespace mkfit::seeding

#endif
