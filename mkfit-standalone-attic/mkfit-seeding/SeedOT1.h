#ifndef mkfit_seeding_SeedOT1_h
#define mkfit_seeding_SeedOT1_h

// Does a quad have a compatible hit in OT1-P?  `seedfind --ot1 OUT.txt`.
//
// A truth study, not a finder stage.  For every quad of the run, the circle
// through its a, c and d hits (signed curvature, tangent at d) and the r-z line
// from a to d along the arc are extended to a target layer, by default mkFit
// layer 4, the P sensors of the first outer-tracker layer.  The prediction is
// made at each candidate hit's OWN radius: the modules are tilted, so a hit's r
// is a measurement, not the layer's.
//
// For each candidate the residuals (dphi, dz) are written together with the
// hit's own error on them, from its full covariance, projected on what a
// residual at fixed r measures:
//
//   res = (J_q - s J_r) . delta,   var = g^T C g,   g = J_q - s J_r
//
// with J_q the gradient of phi or z, J_r that of r, and s the track's d(phi)/dr
// or dz/dr at that radius.  On a sensor facing the IP the r-z term ADDS
// (MKFIT-NOTES, pixel-seeding entry, item 4); the phi term is small but not
// negligible at low pT.
//
// Output, one line per record:
//   Q ev qid cls pte cot eta lab nlab_tgt sim_pt ok ncand
//       cls 0 true, 1 fake (consistent rule), 2 undecidable; lab the label of a
//       true quad, else the label most hits carry (-1 if none); nlab_tgt how
//       many target-layer hits carry lab; ok 1 if the circle reaches the layer
//   H qid dphi dz sphi sz r rel
//       a candidate within the generous fetch window; rel 1 if the hit carries
//       lab, 0 another valid label, -1 unlabelled
//   T qid dphi dz sphi sz r
//       for a true quad, every target-layer hit of its own track, whether
//       fetched or not: the unbiased residual distribution

#include "SeedMath.h"

#include "RecoTracker/MkFitCore/standalone/Event.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

namespace mkfit::seeding {

  struct SeedOT1 {
    FILE *out = nullptr;
    int layer = 4;
    double fetch_phi = 0.02, fetch_z = 0.5;  // generous: rad, cm
    long n_q = 0, n_h = 0, n_t = 0, n_q_true = 0, n_q_true_tgt = 0, n_q_ok = 0;
    long qid = 0;

    struct Pred {
      double k = 0, ux = 0, uy = 0, x0 = 0, y0 = 0, z0 = 0, cots = 0;
    };

    static Pred make_pred(const Hit &a, const Hit &c, const Hit &d) {
      using namespace amath;
      using A = ArithRef;
      const double dx1 = c.x() - a.x(), dy1 = c.y() - a.y(), dx2 = d.x() - c.x(), dy2 = d.y() - c.y();
      const double dx3 = d.x() - a.x(), dy3 = d.y() - a.y();
      const double L1 = std::hypot(dx1, dy1), L2 = std::hypot(dx2, dy2), L3 = std::hypot(dx3, dy3);
      Pred p;
      p.k = 2 * (dx1 * dy2 - dy1 * dx2) / (L1 * L2 * L3);
      const double sa = 0.5 * p.k * L2, ca = std::sqrt(std::max(0.0, 1 - sa * sa));
      const double ex = dx2 / L2, ey = dy2 / L2;
      p.ux = ex * ca - ey * sa;
      p.uy = ex * sa + ey * ca;
      p.x0 = d.x();
      p.y0 = d.y();
      p.z0 = d.z();
      const double s_ad = arc_k<A>(L1, p.k) + arc_k<A>(L2, p.k);
      p.cots = s_ad > 1e-6 ? (d.z() - a.z()) / s_ad : 0;
      return p;
    }

    // (phi, z) of the prediction at radius r; false if the circle does not reach it
    static bool at_r(const Pred &p, double r, double &phi, double &z) {
      using namespace amath;
      using A = ArithRef;
      double px, py, b2;
      if (!curv_cross_r<A>(p.k, p.ux, p.uy, p.x0, p.y0, r, px, py, &b2))
        return false;
      phi = std::atan2(py, px);
      z = p.z0 + p.cots * arc_k<A>(std::hypot(px - p.x0, py - p.y0), p.k);
      return true;
    }

    static double wrap(double d) {
      while (d > M_PI)
        d -= 2 * M_PI;
      while (d < -M_PI)
        d += 2 * M_PI;
      return d;
    }

    // residuals of hit h against the prediction, and the hit's own sigma on each
    static bool residual(const Pred &p, const Hit &h, double &dphi, double &dz, double &sphi, double &sz) {
      const double x = h.x(), y = h.y(), r = std::hypot(x, y);
      double phi0, z0, phi1, z1;
      const double dr = 0.05;
      if (!at_r(p, r, phi0, z0) || !at_r(p, r + dr, phi1, z1))
        return false;
      dphi = wrap(h.phi() - phi0);
      dz = h.z() - z0;
      const double sp = wrap(phi1 - phi0) / dr, sz_ = (z1 - z0) / dr;  // the track's dphi/dr, dz/dr
      const auto &E = h.error();
      const double C[3][3] = {{E.At(0, 0), E.At(0, 1), E.At(0, 2)},
                              {E.At(1, 0), E.At(1, 1), E.At(1, 2)},
                              {E.At(2, 0), E.At(2, 1), E.At(2, 2)}};
      const double Jr[3] = {x / r, y / r, 0};
      const double gp[3] = {-y / (r * r) - sp * Jr[0], x / (r * r) - sp * Jr[1], 0};
      const double gz[3] = {-sz_ * Jr[0], -sz_ * Jr[1], 1};
      double vp = 0, vz = 0;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
          vp += gp[i] * C[i][j] * gp[j];
          vz += gz[i] * C[i][j] * gz[j];
        }
      sphi = std::sqrt(std::max(0.0, vp));
      sz = std::sqrt(std::max(0.0, vz));
      return true;
    }

    void open(const std::string &fname, const SeedParams &P, const LayerInfo &li) {
      out = fopen(fname.c_str(), "w");
      if (!out) {
        fprintf(stderr, "SeedOT1: cannot open %s\n", fname.c_str());
        return;
      }
      fprintf(out, "# seedfind --ot1: target layer %d, r %.3f-%.3f z %.3f-%.3f; fetch |dphi| < %.3f, |dz| < %.3f\n",
              layer, li.rin(), li.rout(), li.zmin(), li.zmax(), fetch_phi, fetch_z);
      fprintf(out, "# windows at this run: c_z %.4f, d_phi %.5f, d_z %.4f; phi_lin %d %.4f; pt_min %.3f\n", P.qwin,
              P.phiwin_d, P.qwin_d, P.phi_lin, P.phi_lin_marg, P.pt_min);
      fprintf(out, "# Q ev qid cls pte cot eta lab nlab_tgt sim_pt ok ncand\n");
      fprintf(out, "# H qid dphi dz sphi sz r rel\n# T qid dphi dz sphi sz r\n");
    }

    template <typename LO>
    void event(int iev, const Event &ev, int la, int lb, int lc, int ld, const LayerInfo &li, const LO &go,
               const std::vector<Quad> &quads) {
      if (!out)
        return;
      const MCHitInfoVec &mc = ev.simHitsInfo_;
      const HitVec *H[4] = {&ev.layerHits_[la], &ev.layerHits_[lb], &ev.layerHits_[lc], &ev.layerHits_[ld]};
      const HitVec &HT = ev.layerHits_[layer];
      auto label = [&mc](const Hit &h) -> int {
        const int mid = h.mcHitID();
        return (mid >= 0 && mid < (int)mc.size()) ? mc[mid].mcTrackID() : -1;
      };
      std::unordered_map<int, std::vector<unsigned int>> tgt;  // label -> target-layer hits
      for (unsigned int k = 0; k < HT.size(); ++k) {
        const int l = label(HT[k]);
        if (l >= 0)
          tgt[l].push_back(k);
      }
      std::vector<unsigned int> cand;
      for (const auto &q : quads) {
        const Hit &ha = (*H[0])[q[0]], &hb = (*H[1])[q[1]], &hc = (*H[2])[q[2]], &hd = (*H[3])[q[3]];
        const int ls[4] = {label(ha), label(hb), label(hc), label(hd)};
        int nun = 0, maxc = 0, lmax = -1;
        for (int u = 0; u < 4; ++u) {
          nun += ls[u] < 0;
          if (ls[u] < 0)
            continue;
          int c = 0;
          for (int w = 0; w < 4; ++w)
            c += ls[w] == ls[u];
          if (c > maxc) {
            maxc = c;
            lmax = ls[u];
          }
        }
        const int cls = maxc == 4 ? 0 : (maxc + nun < 4 ? 1 : 2);
        const Pred p = make_pred(ha, hc, hd);
        const double pte = std::abs(p.k) > 0 ? 0.0114 / std::abs(p.k) : 1e9;
        const double cot = (hd.z() - ha.z()) / (hd.r() - ha.r());
        const auto it = lmax >= 0 ? tgt.find(lmax) : tgt.end();
        const int nlab = it == tgt.end() ? 0 : (int)it->second.size();
        const double sim_pt = (lmax >= 0 && lmax < (int)ev.simTracks_.size()) ? ev.simTracks_[lmax].pT() : -1;

        // fetch: the prediction's span over the layer's radial extent
        double f0, z0, f1, z1;
        const bool ok0 = at_r(p, li.rin(), f0, z0), ok1 = at_r(p, li.rout(), f1, z1);
        const bool ok = ok0 || ok1;
        cand.clear();
        if (ok) {
          if (!ok0)
            f0 = f1, z0 = z1;
          if (!ok1)
            f1 = f0, z1 = z0;
          const double fm = f0 + 0.5 * wrap(f1 - f0), fh = 0.5 * std::abs(wrap(f1 - f0)) + fetch_phi;
          go.for_each_in(go.phi_range(fm - fh, fm + fh),
                         go.q_range(std::min(z0, z1) - fetch_z - 1.0, std::max(z0, z1) + fetch_z + 1.0),
                         [&](unsigned int i) { cand.push_back(go.orig_[i]); });
        }
        struct Row {
          double dphi, dz, sphi, sz, r;
          int rel;
        };
        std::vector<Row> rows;
        for (unsigned int k : cand) {
          Row w;
          if (!residual(p, HT[k], w.dphi, w.dz, w.sphi, w.sz))
            continue;
          if (std::abs(w.dphi) > fetch_phi || std::abs(w.dz) > fetch_z)
            continue;
          w.r = HT[k].r();
          const int l = label(HT[k]);
          w.rel = l < 0 ? -1 : (l == lmax ? 1 : 0);
          rows.push_back(w);
        }
        fprintf(out, "Q %d %ld %d %.5g %.5g %.4f %d %d %.5g %d %zu\n", iev, qid, cls, pte, cot, std::asinh(cot),
                lmax, nlab, sim_pt, ok ? 1 : 0, rows.size());
        for (const Row &w : rows)
          fprintf(out, "H %ld %.4e %.4e %.3e %.3e %.3f %d\n", qid, w.dphi, w.dz, w.sphi, w.sz, w.r, w.rel);
        n_h += rows.size();
        ++n_q;
        n_q_ok += ok;
        if (cls == 0) {
          ++n_q_true;
          n_q_true_tgt += nlab > 0;
          if (it != tgt.end())
            for (unsigned int k : it->second) {
              double dphi, dz, sphi, sz;
              if (!residual(p, HT[k], dphi, dz, sphi, sz))
                continue;
              fprintf(out, "T %ld %.4e %.4e %.3e %.3e %.3f\n", qid, dphi, dz, sphi, sz, HT[k].r());
              ++n_t;
            }
        }
        ++qid;
      }
    }

    void close() {
      if (!out)
        return;
      fclose(out);
      printf("   ot1: layer %d, %ld quads (%ld reach the layer), %.3f candidates per quad in the fetch window\n", layer,
             n_q, n_q_ok, (double)n_h / std::max(1L, n_q));
      printf("   ot1: %ld true quads, %.4f of them with a hit of their own track in the layer; %ld truth rows\n",
             n_q_true, (double)n_q_true_tgt / std::max(1L, n_q_true), n_t);
    }
  };

}  // namespace mkfit::seeding

#endif
