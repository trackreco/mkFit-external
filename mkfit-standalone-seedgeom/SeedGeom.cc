#include "mkfit-standalone-seedgeom/SeedGeom.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"
#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/Track.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <tuple>
#include <vector>
#include <string>

namespace mkfit {

  namespace {

    //--------------------------------------------------------------------------
    // Configuration
    //--------------------------------------------------------------------------
    float g_pt_min = 0.9f;    // GeV -- sets the phi window through R_min
    float g_d0_max = 0.1f;    // cm  -- transverse displacement to tolerate
    float g_qwin = 0.05f;     // cm  -- q tolerance on the third-layer prediction
    int g_la = 0, g_lb = 1, g_lc = 2;
    int g_ld = -1;               // optional 4th layer; < 0 disables
    float g_qwin_d = 0.05f;      // cm,  q window at the 4th layer
    float g_phi_win_d = 0.004f;  // rad, phi window at the 4th layer
    float g_phi_margin = 0.002f; // rad, covers intra-layer radial spread
    int g_nr = 1;                // radial sub-bins per layer (qbar as a bin axis)
    // Layer-c phi window centred on phi extrapolated LINEARLY IN r from the
    // (a,b) pair, instead of on phi_b with the generic pT_min + D0 band.  Linear
    // in r absorbs the curvature (phi ~ phi0 - r/2R); the D0 term goes as D0/r
    // and is not linear, so the tolerance carries D0_max times the second
    // difference of 1/r, plus g_phi_lin_marg for MS and hit resolution.
    int g_phi_lin = 0;  // 0 off, 1 linear window only, 2 linear AND the generic pT_min/D0 band
    float g_phi_lin_marg = 0.003f;  // rad
    int g_nphi = 1024, g_nz = 512;  // grid granularity (nphi must be a power of 2)
    bool g_light = false;           // skip the per-triplet diagnostic vectors
    bool g_do_partb = true;
    bool g_verbose = false;

    constexpr float kBfield = 3.8f;  // T; R[cm] = pT / (0.003 * B)
    inline float r_min_cm() { return g_pt_min / (0.003f * kBfield); }

    constexpr float kPi = 3.14159265358979323846f;
    constexpr float kTwoPi = 2.0f * kPi;

    //--------------------------------------------------------------------------
    // PART A accumulators: truth vertices
    //--------------------------------------------------------------------------
    struct VtxAcc {
      int n_all = 0;
      int n_pt09 = 0;
      int n_pt2 = 0;
      int n_prod[4] = {0, 0, 0, 0};
    };

    int A_n_events = 0;
    long A_n_simtracks = 0;
    std::vector<int> A_nvtx_per_ev, A_nprim_per_ev;
    std::vector<double> A_prim_z, A_prim_r, A_prim_dz;  // dz = nearest-neighbour spacing
    std::vector<double> A_prim_mult, A_prim_mult09;
    std::vector<double> A_sec_r;
    long A_prod_counts[4] = {0, 0, 0, 0};
    long A_n_signal_vtx = 0;
    long A_n_onbeam_points = 0;
    std::vector<int> A_nprim09_per_ev;
    float g_bl_rcut = 0.02f;   // cm, transverse cut defining "on the beamline"
    float g_vtx_gap = 0.005f;  // cm, z gap that separates two vertices
    double A_bs_x = 0, A_bs_y = 0, A_bs_sigz = 0;

    //--------------------------------------------------------------------------
    // PART B accumulators: geometric triplets
    //--------------------------------------------------------------------------
    int B_n_events = 0;
    long B_n_a = 0, B_n_b = 0, B_n_c = 0;   // hits in the three layers
    long B_n_pairs = 0;                     // doublets admitted by the phi window
    long B_n_cells = 0, B_n_ctouch = 0;     // work done at layer c
    double B_t_partb = 0;                   // s, wall clock inside part_b()
    double B_t_search = 0;                  // s, the doublet/triplet/quad loop alone
    long B_n_trip_phi = 0;                  // triplets on the phi window alone (q ignored)
    long B_n_trip = 0;                      // triplets also passing the q prediction
    long B_n_trip_true = 0;                 // ... and all three hits from one sim track
    long B_n_trip_avail = 0;                // sim tracks with a hit in all three layers
    long B_n_trip_found = 0;
    long B_n_cells_d = 0, B_n_dtouch = 0, B_n_quad = 0, B_n_quad_true = 0;
    long B_n_avail_all = 0;
    long B_n_quad_nolink = 0, B_n_quad_fake = 0, B_n_quad_3of4 = 0;
    long B_n_hit_lbl[4] = {0, 0, 0, 0}, B_n_hit_tot[4] = {0, 0, 0, 0};
    // per-event scratch, handed from part_a to the vertex finder
    std::vector<double> A_cur_z;
    std::vector<int> A_cur_n09, A_cur_nall;
    std::vector<double> B_ev_zv;  // z_v votes of this event
    // vertex-finder results
    float g_vf_win = 0.02f;  // cm, half-width of the vote window
    int g_vf_min = 3;        // minimum votes to call a vertex
    float g_vf_match = 0.05f;  // cm, truth match radius
    long V_n_true = 0, V_n_true09 = 0, V_n_found = 0, V_n_match = 0, V_n_match09 = 0, V_n_fake = 0;
    std::vector<double> V_dz, V_eff_mult_num, V_eff_mult_den;
    std::vector<double> B_eff_pt_num, B_eff_pt_den;
    std::vector<double> B_eff_eta_num, B_eff_eta_den;
    // purity resolved in |eta|, taken from the candidate's OWN cot theta so it
    // needs no truth: 5 bins x {all, true, genuine fake, undecidable}
    static constexpr int NEB = 6;
    long B_pur_all[NEB] = {0}, B_pur_true[NEB] = {0}, B_pur_fake[NEB] = {0}, B_pur_nol[NEB] = {0};
    inline int etabin(double e) {
      const double a = std::abs(e);
      return a < 0.4 ? 0 : a < 0.8 ? 1 : a < 1.2 ? 2 : a < 1.6 ? 3 : a < 2.2 ? 4 : 5;
    }
    // per TRUE quad, for splitting the 4th-layer q width into what the hit's
    // CENTROID radius does (|cot| * sigma_r) and what multiple scattering does (1/p)
    std::vector<double> B_q4_dz, B_q4_csr, B_q4_pt, B_q4_sz, B_q4_cot, B_q4_sr, B_q4_sres;
    // with per-sim-hit truth states: the residual split into its two halves, both
    // taken at the TRUE crossing's radius.  e_hit = hit - truth, e_pred = prediction - truth.
    std::vector<double> B_q4_ehit, B_q4_epred, B_q4_tcot, B_q4_tpt, B_q4_ttilt;
    std::vector<double> B_dphi, B_dzd, B_dphi_all, B_dzd_all;                // ... of which at least one true triplet was built

    std::vector<double> B_zres;      // z_c - z_pred, TRUE triplets [cm]
    std::vector<double> B_zres_all;  // ... all triplets
    std::vector<double> B_zv;        // fitted beamline intercept, all triplets
    std::vector<double> B_zv_true;   // ... true triplets
    std::vector<double> B_zvres;     // z_v(fit) - z_v(truth), true triplets, primaries
    std::vector<double> B_zvres_pt;  // matching pT
    std::vector<double> B_ptres;     // (pT_fit - pT_true)/pT_true, true triplets
    std::vector<double> B_slres;     // longitudinal residual of the 3-point test [cm]

    //--------------------------------------------------------------------------
    // A minimal (phi, z) CSR grid -- the binnor in miniature.  Key is
    // iphi * NZ + iz, so a whole phi column is one contiguous range, which is
    // what the doublet stage wants, while a (phi, z) cell is what the
    // third-layer confirmation wants.
    //--------------------------------------------------------------------------
    struct HitGrid {
      int nphi = 1024, nphi_mask = 1023, nz = 512;
      int nr = 1;  // radial sub-bins; qbar (== hit r in the barrel) is the axis
      float rlo = 0, rhi = 0, rfac = 0;
      int ncell = 1024 * 512;
      float zmin = 0, zmax = 0, zfac = 0;
      std::vector<int> start;  // NCELL + 1
      std::vector<int> idx;    // hit indices grouped by cell
      std::vector<float> phi, z, r, invr;
      int n = 0;

      inline int phibin(float p) const {
        int b = (int)std::floor((p + kPi) * (nphi / kTwoPi));
        return b & nphi_mask;
      }
      inline int zbin(float zz) const {
        int b = (int)((zz - zmin) * zfac);
        return b < 0 ? 0 : (b >= nz ? nz - 1 : b);
      }
      inline int rbin(float rr) const {
        int b = (int)((rr - rlo) * rfac);
        return b < 0 ? 0 : (b >= nr ? nr - 1 : b);
      }
      inline int cell(int iphi, int ir, int iz) const { return (iphi * nr + ir) * nz + iz; }
      inline float rsub_mid(int ir) const { return rlo + (ir + 0.5f) / rfac; }
      inline float rsub_half() const { return 0.5f / rfac; }

      void build(const HitVec &hits, int n_r, int n_phi, int n_z) {
        nphi = n_phi;
        nphi_mask = n_phi - 1;
        nz = n_z;
        nr = std::max(1, n_r);
        ncell = nphi * nr * nz;
        n = (int)hits.size();
        phi.resize(n);
        z.resize(n);
        r.resize(n);
        invr.resize(n);
        zmin = 1e9f;
        zmax = -1e9f;
        for (int i = 0; i < n; ++i) {
          phi[i] = hits[i].phi();
          z[i] = hits[i].z();
          r[i] = hits[i].r();
          invr[i] = 1.0f / r[i];
          zmin = std::min(zmin, z[i]);
          zmax = std::max(zmax, z[i]);
        }
        if (zmax <= zmin)
          zmax = zmin + 1.0f;
        zfac = nz / (zmax - zmin + 1e-4f);
        rlo = *std::min_element(r.begin(), r.end());
        rhi = *std::max_element(r.begin(), r.end());
        rfac = nr / (rhi - rlo + 1e-5f);

        start.assign(ncell + 1, 0);
        std::vector<int> cl(n);
        for (int i = 0; i < n; ++i) {
          cl[i] = cell(phibin(phi[i]), rbin(r[i]), zbin(z[i]));
          ++start[cl[i] + 1];
        }
        for (int c = 0; c < ncell; ++c)
          start[c + 1] += start[c];
        idx.resize(n);
        std::vector<int> fill(start.begin(), start.end() - 1);
        for (int i = 0; i < n; ++i)
          idx[fill[cl[i]]++] = i;
      }
    };

    //--------------------------------------------------------------------------
    // helpers
    //--------------------------------------------------------------------------
    // Branch-free: the arguments here are differences of two angles already in
    // (-pi, pi], so at most one wrap is ever needed.
    inline float wrap_pi(float d) {
      d -= kTwoPi * (d > kPi);
      d += kTwoPi * (d < -kPi);
      return d;
    }

    double quant(std::vector<double> &v, double q) {
      if (v.empty())
        return 0;
      size_t k = (size_t)(q * (v.size() - 1));
      std::nth_element(v.begin(), v.begin() + k, v.end());
      return v[k];
    }
    // robust sigma from the inter-quartile range of a Gaussian
    double iqr_sigma(std::vector<double> v) {
      if (v.size() < 4)
        return 0;
      return 0.5 * (quant(v, 0.75) - quant(v, 0.25)) / 0.6745;
    }
    double med(std::vector<double> v) { return quant(v, 0.5); }

    // circle through three 2-D points; returns false if (near-)collinear
    bool circle3(float x1, float y1, float x2, float y2, float x3, float y3, double &cx, double &cy, double &R) {
      double A = x2 - x1, B = y2 - y1, C = x3 - x1, D = y3 - y1;
      double E = A * (x1 + x2) + B * (y1 + y2);
      double F = C * (x1 + x3) + D * (y1 + y3);
      double G = 2.0 * (A * (y3 - y2) - B * (x3 - x2));
      if (std::abs(G) < 1e-12)
        return false;
      cx = (D * E - B * F) / G;
      cy = (A * F - C * E) / G;
      R = std::hypot(x1 - cx, y1 - cy);
      return R > 1e-6 && std::isfinite(R);
    }

    // arc length along a circle of radius R subtending a chord c
    inline double arc(double chord, double R) {
      double h = chord / (2.0 * R);
      if (h > 1.0)
        h = 1.0;
      return 2.0 * R * std::asin(h);
    }

  }  // namespace

  //============================================================================
  // configuration entry points
  //============================================================================
  void sg_set_ptmin(float v) { g_pt_min = v; }
  void sg_set_d0max(float v) { g_d0_max = v; }
  void sg_set_qwin(float v) { g_qwin = v; }
  void sg_set_layers(int la, int lb, int lc) {
    g_la = la;
    g_lb = lb;
    g_lc = lc;
  }
  void sg_set_do_partb(bool on) { g_do_partb = on; }
  void sg_set_verbose(bool on) { g_verbose = on; }
  void sg_set_vtx_gap(float v) { g_vtx_gap = v; }
  void sg_set_bl_rcut(float v) { g_bl_rcut = v; }
  void sg_set_layer4(int ld) { g_ld = ld; }
  void sg_set_qwin_d(float v) { g_qwin_d = v; }
  void sg_set_phiwin_d(float v) { g_phi_win_d = v; }
  void sg_set_phi_margin(float v) { g_phi_margin = v; }
  void sg_set_nr(int n) { g_nr = n; }
  void sg_set_phi_lin(int on, float marg) {
    g_phi_lin = on;
    g_phi_lin_marg = marg;
  }
  void sg_set_light(bool on) { g_light = on; }
  void sg_set_grid(int nphi, int nz) {
    g_nphi = nphi;
    g_nz = nz;
  }
  void sg_set_vf(float win_cm, int min_votes, float match_cm) {
    g_vf_win = win_cm;
    g_vf_min = min_votes;
    g_vf_match = match_cm;
  }

  void sg_reset() {
    A_n_events = 0;
    A_n_simtracks = 0;
    A_nvtx_per_ev.clear();
    A_nprim_per_ev.clear();
    A_prim_z.clear();
    A_prim_r.clear();
    A_prim_dz.clear();
    A_prim_mult.clear();
    A_prim_mult09.clear();
    A_sec_r.clear();
    for (int i = 0; i < 4; ++i)
      A_prod_counts[i] = 0;
    A_n_signal_vtx = 0;
    A_n_onbeam_points = 0;
    A_nprim09_per_ev.clear();

    B_n_events = 0;
    B_n_a = B_n_b = B_n_c = 0;
    B_n_pairs = B_n_cells = B_n_ctouch = 0;
    B_t_partb = B_t_search = 0;
    B_n_trip_phi = B_n_trip = B_n_trip_true = 0;
    B_n_trip_avail = B_n_trip_found = 0;
    B_n_cells_d = B_n_dtouch = B_n_quad = B_n_quad_true = 0;
    B_n_avail_all = 0;
    B_n_quad_nolink = B_n_quad_fake = B_n_quad_3of4 = 0;
    for (int i = 0; i < 4; ++i) {
      B_n_hit_lbl[i] = 0;
      B_n_hit_tot[i] = 0;
    }
    V_n_true = V_n_true09 = V_n_found = V_n_match = V_n_match09 = V_n_fake = 0;
    V_dz.clear();
    V_eff_mult_num.clear();
    V_eff_mult_den.clear();
    B_eff_pt_num.clear();
    B_eff_pt_den.clear();
    B_eff_eta_num.clear();
    B_eff_eta_den.clear();
    for (int i = 0; i < NEB; ++i) {
      B_pur_all[i] = B_pur_true[i] = B_pur_fake[i] = B_pur_nol[i] = 0;
    }
    B_dphi.clear();
    B_dzd.clear();
    B_q4_dz.clear();
    B_q4_csr.clear();
    B_q4_pt.clear();
    B_q4_sz.clear();
    B_q4_cot.clear();
    B_q4_sr.clear();
    B_q4_sres.clear();
    B_q4_ehit.clear();
    B_q4_epred.clear();
    B_q4_tcot.clear();
    B_q4_tpt.clear();
    B_q4_ttilt.clear();
    B_dphi_all.clear();
    B_dzd_all.clear();
    B_zres.clear();
    B_zres_all.clear();
    B_zv.clear();
    B_zv_true.clear();
    B_zvres.clear();
    B_zvres_pt.clear();
    B_ptres.clear();
    B_slres.clear();
    printf("[sg] reset. pt_min=%.2f GeV  d0_max=%.3f cm  layers=(%d,%d,%d)  qwin=%.3f cm\n",
           g_pt_min, g_d0_max, g_la, g_lb, g_lc, g_qwin);
  }

  //============================================================================
  // PART A -- truth vertices by exact group-by on the production point
  //============================================================================
  namespace {
    void part_a(Event *ev) {
      const TrackVec &st = ev->simTracks_;
      A_n_events++;
      A_n_simtracks += (long)st.size();
      A_bs_x = ev->beamSpot_.x;
      A_bs_y = ev->beamSpot_.y;
      A_bs_sigz = ev->beamSpot_.sigmaZ;

      // Group on the exact (x,y,z) bit pattern first -- that part DOES work,
      // the big PU vertices come out as single points of 250-450 tracks.
      std::map<std::tuple<float, float, float>, VtxAcc> vmap;
      for (const auto &t : st) {
        VtxAcc &a = vmap[std::make_tuple(t.x(), t.y(), t.z())];
        a.n_all++;
        float pt = t.pT();
        if (pt > 0.9f)
          a.n_pt09++;
        if (pt > 2.0f)
          a.n_pt2++;
        int pti = (int)t.prodType();
        if (pti >= 0 && pti < 4) {
          a.n_prod[pti]++;
          A_prod_counts[pti]++;
        }
      }

      // A vertex is a z CLUSTER of on-beamline production points, not a single
      // point: alongside each PU vertex sit a few hundred low-multiplicity
      // points at the same z, from decays that happen within the beam width.
      struct P {
        double z;
        VtxAcc a;
      };
      std::vector<P> on;
      for (const auto &kv : vmap) {
        float x = std::get<0>(kv.first), y = std::get<1>(kv.first), z = std::get<2>(kv.first);
        double dr = std::hypot(x - A_bs_x, y - A_bs_y);
        if (dr < g_bl_rcut)
          on.push_back({z, kv.second});
        else
          A_sec_r.push_back(dr);
      }
      std::sort(on.begin(), on.end(), [](const P &a, const P &b) { return a.z < b.z; });

      std::vector<double> zc;
      std::vector<int> nc, nc09, nsig;
      for (size_t i = 0; i < on.size();) {
        size_t j = i;
        double zsum = 0;
        int n = 0, n09 = 0, ns = 0;
        while (j < on.size() && (j == i || on[j].z - on[j - 1].z < g_vtx_gap)) {
          zsum += on[j].z * on[j].a.n_all;
          n += on[j].a.n_all;
          n09 += on[j].a.n_pt09;
          ns += on[j].a.n_prod[(int)Track::ProdType::Signal];
          ++j;
        }
        zc.push_back(zsum / std::max(1, n));
        nc.push_back(n);
        nc09.push_back(n09);
        nsig.push_back(ns);
        i = j;
      }

      int nprim = 0, nprim09 = 0;
      for (size_t k = 0; k < zc.size(); ++k) {
        A_prim_z.push_back(zc[k]);
        A_prim_mult.push_back(nc[k]);
        A_prim_mult09.push_back(nc09[k]);
        ++nprim;
        if (nc09[k] >= 2)
          ++nprim09;
        if (nsig[k] > 0)
          A_n_signal_vtx++;
      }
      for (size_t k = 1; k < zc.size(); ++k)
        A_prim_dz.push_back(zc[k] - zc[k - 1]);

      A_cur_z = zc;
      A_cur_n09.assign(nc09.begin(), nc09.end());
      A_cur_nall.assign(nc.begin(), nc.end());
      A_nvtx_per_ev.push_back((int)vmap.size());
      A_nprim_per_ev.push_back(nprim);
      A_nprim09_per_ev.push_back(nprim09);
      A_n_onbeam_points += (long)on.size();
    }
  }  // namespace

  //============================================================================
  // PART B -- geometric doublets and triplets
  //============================================================================
  namespace {
    // where does the circle (cx,cy,R) cross radius rt?  Returns the crossing
    // closer to (xref,yref) -- the continuation of the track, not the far side
    // of the loop.  Same radical-line construction as propagate_to_r.
    bool circle_cross_r(
        double cx, double cy, double R, double rt, double xref, double yref, double &px, double &py) {
      const double dC = std::hypot(cx, cy);
      if (dC < 1e-9)
        return false;
      const double aa = (rt * rt + dC * dC - R * R) / (2 * dC);
      const double b2 = rt * rt - aa * aa;
      if (b2 < 0)
        return false;
      const double bb = std::sqrt(b2);
      const double ux = cx / dC, uy = cy / dC;  // C-hat
      const double vx = -uy, vy = ux;           // C-hat perp
      const double x1 = aa * ux + bb * vx, y1 = aa * uy + bb * vy;
      const double x2 = aa * ux - bb * vx, y2 = aa * uy - bb * vy;
      if (std::hypot(x1 - xref, y1 - yref) <= std::hypot(x2 - xref, y2 - yref)) {
        px = x1;
        py = y1;
      } else {
        px = x2;
        py = y2;
      }
      return true;
    }

    void record_truth(Event *ev, int label, double zv, double R) {
      if (label < 0 || label >= (int)ev->simTracks_.size())
        return;
      const Track &tr = ev->simTracks_[label];
      if (std::hypot(tr.x() - A_bs_x, tr.y() - A_bs_y) < g_bl_rcut) {
        B_zvres.push_back(zv - tr.z());
        B_zvres_pt.push_back(tr.pT());
      }
      B_ptres.push_back((0.003 * kBfield * R - tr.pT()) / tr.pT());
    }


    //--------------------------------------------------------------------------
    // Vertex finding from the z_v votes alone.  No tracking, no fit -- the
    // intercept of every geometric triplet/quadruplet is one vote, and a
    // vertex is a peak.  Greedy: repeatedly take the densest window, take its
    // weighted mean, remove its votes.
    //--------------------------------------------------------------------------
    void vertex_find() {
      if (A_cur_z.empty())
        return;
      std::vector<double> v = B_ev_zv;
      std::sort(v.begin(), v.end());
      std::vector<char> used(v.size(), 0);
      std::vector<double> found;

      const double w = g_vf_win;
      while (true) {
        // densest window of half-width w over the unused votes
        size_t best_i = 0, best_n = 0, best_j = 0;
        size_t j = 0;
        for (size_t i = 0; i < v.size(); ++i) {
          if (used[i])
            continue;
          if (j < i)
            j = i;
          while (j + 1 < v.size() && v[j + 1] - v[i] <= 2 * w)
            ++j;
          size_t n = 0;
          for (size_t k = i; k <= j; ++k)
            if (!used[k])
              ++n;
          if (n > best_n) {
            best_n = n;
            best_i = i;
            best_j = j;
          }
        }
        if ((int)best_n < g_vf_min)
          break;
        double sum = 0;
        int n = 0;
        for (size_t k = best_i; k <= best_j; ++k)
          if (!used[k]) {
            sum += v[k];
            used[k] = 1;
            ++n;
          }
        found.push_back(sum / n);
      }

      // match to truth
      std::vector<char> tmatched(A_cur_z.size(), 0);
      std::vector<char> fmatched(found.size(), 0);
      for (size_t f = 0; f < found.size(); ++f) {
        int best = -1;
        double bd = 1e9;
        for (size_t t = 0; t < A_cur_z.size(); ++t) {
          if (tmatched[t])
            continue;
          double d = std::abs(found[f] - A_cur_z[t]);
          if (d < bd) {
            bd = d;
            best = (int)t;
          }
        }
        if (best >= 0 && bd < g_vf_match) {
          tmatched[best] = 1;
          fmatched[f] = 1;
          V_dz.push_back(found[f] - A_cur_z[best]);
        }
      }
      V_n_found += (long)found.size();
      for (size_t f = 0; f < found.size(); ++f)
        if (!fmatched[f])
          V_n_fake++;
      for (size_t t = 0; t < A_cur_z.size(); ++t) {
        V_n_true++;
        V_eff_mult_den.push_back(A_cur_n09[t]);
        if (tmatched[t]) {
          V_n_match++;
          V_eff_mult_num.push_back(A_cur_n09[t]);
        }
        if (A_cur_n09[t] >= 2) {
          V_n_true09++;
          if (tmatched[t])
            V_n_match09++;
        }
      }
    }

    void part_b(Event *ev) {
      const int nl = (int)ev->layerHits_.size();
      const bool use_d = (g_ld >= 0);
      if (g_la >= nl || g_lb >= nl || g_lc >= nl || (use_d && g_ld >= nl))
        return;
      const HitVec &ha = ev->layerHits_[g_la];
      const HitVec &hb = ev->layerHits_[g_lb];
      const HitVec &hc = ev->layerHits_[g_lc];
      if (ha.empty() || hb.empty() || hc.empty())
        return;
      const HitVec *hd = use_d ? &ev->layerHits_[g_ld] : nullptr;
      if (use_d && hd->empty())
        return;

      B_n_events++;
      B_n_a += (long)ha.size();
      B_n_b += (long)hb.size();
      B_n_c += (long)hc.size();

      B_ev_zv.clear();
      HitGrid ga, gb, gc, gd;
      ga.build(ha, 1, g_nphi, g_nz);
      gb.build(hb, 1, g_nphi, g_nz);
      gc.build(hc, g_nr, g_nphi, g_nz);
      if (use_d)
        gd.build(*hd, g_nr, g_nphi, g_nz);

      // truth bookkeeping
      const MCHitInfoVec &mc = ev->simHitsInfo_;
      auto mctrk = [&mc](const Hit &h) -> int {
        int mid = h.mcHitID();
        return (mid >= 0 && mid < (int)mc.size()) ? mc[mid].mcTrackID() : -1;
      };
      std::map<int, unsigned> present;
      auto mark = [&](const HitVec &hv, unsigned bit) {
        for (const auto &h : hv) {
          int l = mctrk(h);
          if (l >= 0)
            present[l] |= bit;
        }
      };
      mark(ha, 1u);
      mark(hb, 2u);
      mark(hc, 4u);
      if (use_d)
        mark(*hd, 8u);
      {  // how many hits carry a truth link at all?
        const HitVec *hv[4] = {&ha, &hb, &hc, hd};
        for (int i = 0; i < (use_d ? 4 : 3); ++i) {
          for (const auto &h : *hv[i]) {
            B_n_hit_tot[i]++;
            if (mctrk(h) >= 0)
              B_n_hit_lbl[i]++;
          }
        }
      }
      // FINDABLE denominator.  The phi stencil is built for pT >= pT_min and
      // |D0| <= D0_max; a 0.2 GeV curler inside it is not a miss, it is out of
      // scope.  Quoting efficiency against every sim track with three hits
      // would be quoting against a ceiling the window cannot reach.
      const unsigned want = use_d ? 15u : 7u;
      std::map<int, bool> found_true;
      for (const auto &kv : present) {
        if ((kv.second & want) != want)
          continue;
        const int l = kv.first;
        if (l >= (int)ev->simTracks_.size())
          continue;
        const Track &tr = ev->simTracks_[l];
        B_n_avail_all++;
        if (tr.pT() < g_pt_min)
          continue;
        if (std::hypot(tr.x() - A_bs_x, tr.y() - A_bs_y) > g_d0_max)
          continue;
        B_n_trip_avail++;
        found_true[l] = false;
        B_eff_pt_den.push_back(tr.pT());
        B_eff_eta_den.push_back(tr.momEta());
      }

      // Per-layer radial extent.  The pixel barrel layers are ~1 cm thick and
      // that thickness dominates both windows if it is not accounted for: a
      // prediction made at the layer's MEAN radius is compared against a hit
      // that can sit 0.5 cm inside or outside it, which for |cot theta| ~ 2 is
      // a centimetre of z.
      auto rminmax = [](const HitGrid &g, double &lo, double &hi, double &mean) {
        lo = 1e9;
        hi = -1e9;
        double sum = 0;
        for (int i = 0; i < g.n; ++i) {
          lo = std::min(lo, (double)g.r[i]);
          hi = std::max(hi, (double)g.r[i]);
          sum += g.r[i];
        }
        mean = sum / std::max(1, g.n);
      };
      double ralo, rahi, ram, rblo, rbhi, rbm, rclo, rchi, rcm, rdlo = 0, rdhi = 0, rdm = 0;
      rminmax(ga, ralo, rahi, ram);
      rminmax(gb, rblo, rbhi, rbm);
      rminmax(gc, rclo, rchi, rcm);
      if (use_d)
        rminmax(gd, rdlo, rdhi, rdm);
      const double Rmin = r_min_cm();
      (void)rclo;
      (void)rdlo;

      // phi stencil half-width: curvature term + displacement term, taken at
      // the widest radius pairing the two layers allows.
      // inv_r comes precomputed from the grid, so the window test has no
      // division at all: two FMAs and an add.
      const float inv2R = (float)(1.0 / (2 * Rmin)), d0m = g_d0_max, marg = g_phi_margin;
      const float inv_rbhi = 1.0f / (float)rbhi, inv_rchi = 1.0f / (float)rchi;
      auto wphi_w = [=](float r_in, float inv_in, float r_out, float inv_out) {
        return (r_out - r_in) * inv2R + d0m * (inv_in - inv_out) + marg;
      };

      if (g_verbose)
        printf("[sg]   r: a %.2f (%.2f-%.2f)  b %.2f (%.2f-%.2f)  c %.2f (%.2f-%.2f)%s\n",
               ram, ralo, rahi, rbm, rblo, rbhi, rcm, rclo, rchi,
               use_d ? Form("  d %.2f (%.2f-%.2f)", rdm, rdlo, rdhi) : "");

      const auto t_search0 = std::chrono::steady_clock::now();
      for (int ia = 0; ia < ga.n; ++ia) {
        const float pa = ga.phi[ia], za = ga.z[ia], rra = ga.r[ia];
        const double xa = rra * std::cos(pa), ya = rra * std::sin(pa);
        const float w_ab = wphi_w(rra, ga.invr[ia], (float)rbhi, inv_rbhi);
        const int nb = std::min(gb.nphi, 2 * (int)std::ceil(w_ab * gb.nphi / kTwoPi) + 1);
        const int pb0 = gb.phibin(pa) - nb / 2;

        for (int k = 0; k < nb; ++k) {
          const int c0 = gb.cell((pb0 + k) & gb.nphi_mask, 0, 0);
          for (int s = gb.start[c0]; s < gb.start[c0 + gb.nz]; ++s) {
            const int ib = gb.idx[s];
            const float rrb = gb.r[ib];
            if (rrb <= rra + 0.1f)
              continue;
            const float pbph = gb.phi[ib], zb = gb.z[ib];
            if (std::abs(wrap_pi(pbph - pa)) > wphi_w(rra, ga.invr[ia], rrb, gb.invr[ib]))
              continue;
            B_n_pairs++;

            // q prediction, evaluated at the TARGET HIT's own radius.  The bin
            // range must cover the layer's radial thickness; the per-hit cut
            // must not.
            const double cot = (zb - za) / (rrb - rra);
            // linear-in-r phi extrapolation from (a,b); slope in rad/cm
            const float slope = wrap_pi(pbph - pa) / (rrb - rra);
            const float kab = 1.0f / (rrb - rra);
            const float dinv_ab = gb.invr[ib] - ga.invr[ia];
            // capture SCALARS only: [=] on gb would copy the whole grid per doublet
            const float inv_b = gb.invr[ib], mlin = g_phi_lin_marg;
            auto wlin = [d0m, inv_b, rrb, kab, dinv_ab, mlin](float r_c, float inv_c) {
              return d0m * std::abs(inv_c - inv_b - (r_c - rrb) * kab * dinv_ab) + mlin;
            };
            float w_bc, phic_ctr;
            if (g_phi_lin) {
              // bin range: the centre moves with r across the layer's thickness,
              // and the D0 term is largest at the layer's outer edge
              const float rmid = 0.5f * (float)(rclo + rchi);
              phic_ctr = pbph + slope * (rmid - rrb);
              w_bc = std::max(wlin((float)rclo, 1.0f / (float)rclo), wlin((float)rchi, inv_rchi)) +
                     std::abs(slope) * 0.5f * (float)(rchi - rclo);
            } else {
              phic_ctr = pbph;
              w_bc = wphi_w(rrb, gb.invr[ib], (float)rchi, inv_rchi);
            }
            const int nc2 = std::min(gc.nphi, 2 * (int)std::ceil(w_bc * gc.nphi / kTwoPi) + 1);
            const int pc0 = gc.phibin(phic_ctr) - nc2 / 2;
            const double rsub_h = gc.rsub_half();

            for (int m = 0; m < nc2; ++m) {
              const int iphi = (pc0 + m) & gc.nphi_mask;
              for (int ir = 0; ir < gc.nr; ++ir) {
                // Each radial sub-bin gets its own z window.  The bounding box
                // over the whole layer has z-extent qwin + |cot| * dr_half;
                // split into n_r shells the SCANNED AREA falls as
                // 2 qwin dr + |cot| dr^2 / n_r, so the over-scan converges to 1.
                const double zc_mid = za + cot * (gc.rsub_mid(ir) - rra);
                const double zhalf = g_qwin + std::abs(cot) * rsub_h;
                if (zc_mid + zhalf < gc.zmin || zc_mid - zhalf > gc.zmax)
                  continue;
                const int zlo = gc.zbin(zc_mid - zhalf), zhi = gc.zbin(zc_mid + zhalf);
                const int base = gc.cell(iphi, ir, 0);
                B_n_cells += (zhi - zlo + 1);
                for (int s2 = gc.start[base + zlo]; s2 < gc.start[base + zhi + 1]; ++s2) {
                const int ic = gc.idx[s2];
                B_n_ctouch++;
                const float rrc = gc.r[ic], pcph = gc.phi[ic];
                if (rrc <= rrb + 0.1f)
                  continue;
                if (g_phi_lin) {
                  if (std::abs(wrap_pi(pcph - pbph - slope * (rrc - rrb))) > wlin(rrc, gc.invr[ic]))
                    continue;
                }
                if (g_phi_lin != 1 && std::abs(wrap_pi(pcph - pbph)) > wphi_w(rrb, gb.invr[ib], rrc, gc.invr[ic]))
                  continue;
                B_n_trip_phi++;
                const double zpred = za + cot * (rrc - rra);
                const double dz = gc.z[ic] - zpred;
                if (!g_light)
                  B_zres_all.push_back(dz);
                if (std::abs(dz) > g_qwin)
                  continue;
                B_n_trip++;

                // the one real constraint, done properly: z linear in ARC length
                const double xb = rrb * std::cos(pbph), yb = rrb * std::sin(pbph);
                const double xc = rrc * std::cos(pcph), yc = rrc * std::sin(pcph);
                double cx, cy, R;
                if (!circle3(xa, ya, xb, yb, xc, yc, cx, cy, R))
                  R = 1e6;
                const double s_ab = arc(std::hypot(xb - xa, yb - ya), R);
                const double s_bc = arc(std::hypot(xc - xb, yc - yb), R);
                const double s_ac = s_ab + s_bc;
                const double cot_s = (s_ac > 1e-6) ? (gc.z[ic] - za) / s_ac : 0.0;
                if (!g_light)
                  B_slres.push_back(zb - (za + cot_s * s_ab));
                const double zv = za - cot_s * arc(rra, R);
                if (!g_light)
                  B_zv.push_back(zv);

                const int la = mctrk(ha[ia]), lb = mctrk(hb[ib]), lc = mctrk(hc[ic]);
                const bool true_trip = (la >= 0 && la == lb && la == lc);
                if (true_trip) {
                  B_n_trip_true++;
                  B_zres.push_back(dz);
                }

                //---- optional 4th layer.  The helix is now DETERMINED, so this
                //     is a prediction, not a search.
                if (!use_d) {
                  B_ev_zv.push_back(zv);
                  if (true_trip) {
                    B_zv_true.push_back(zv);
                    if (auto it = found_true.find(la); it != found_true.end())
                        it->second = true;
                    record_truth(ev, la, zv, R);
                  }
                  continue;
                }
                double pxd, pyd;
                if (!circle_cross_r(cx, cy, R, rdm, xc, yc, pxd, pyd))
                  continue;
                const double phid = std::atan2(pyd, pxd);
                const int ndb =
                    std::min(gd.nphi, 2 * (int)std::ceil(g_phi_win_d * gd.nphi / kTwoPi) + 1);
                const int pd0 = gd.phibin((float)phid) - ndb / 2;
                const double rdsub_h = gd.rsub_half();
                for (int q = 0; q < ndb; ++q) {
                  const int iphid = (pd0 + q) & gd.nphi_mask;
                  for (int jr = 0; jr < gd.nr; ++jr) {
                  double pxs, pys;
                  if (!circle_cross_r(cx, cy, R, gd.rsub_mid(jr), xc, yc, pxs, pys))
                    continue;
                  const double zds = gc.z[ic] + cot_s * arc(std::hypot(pxs - xc, pys - yc), R);
                  const double zdh = g_qwin_d + std::abs(cot_s) * rdsub_h;
                  if (zds + zdh < gd.zmin || zds - zdh > gd.zmax)
                    continue;
                  const int zdlo = gd.zbin(zds - zdh), zdhi = gd.zbin(zds + zdh);
                  const int bd = gd.cell(iphid, jr, 0);
                  B_n_cells_d += (zdhi - zdlo + 1);
                  for (int s3 = gd.start[bd + zdlo]; s3 < gd.start[bd + zdhi + 1]; ++s3) {
                    const int id = gd.idx[s3];
                    B_n_dtouch++;
                    const float rrd = gd.r[id];
                    if (rrd <= rrc + 0.1f)
                      continue;
                    double px2, py2;  // re-predict at this hit's own radius
                    if (!circle_cross_r(cx, cy, R, rrd, xc, yc, px2, py2))
                      continue;
                    const double dphi_d = wrap_pi((float)(gd.phi[id] - std::atan2(py2, px2)));
                    const double zd2 = gc.z[ic] + cot_s * arc(std::hypot(px2 - xc, py2 - yc), R);
                    const double dz_d = gd.z[id] - zd2;
                    if (!g_light) {
                      B_dphi_all.push_back(dphi_d);
                      B_dzd_all.push_back(dz_d);
                    }
                    if (std::abs(dphi_d) > g_phi_win_d || std::abs(dz_d) > g_qwin_d)
                      continue;
                    B_n_quad++;
                    B_ev_zv.push_back(zv);
                    const int ld = mctrk((*hd)[id]);
                    {
                      // Impurity has two causes and they need separating: a
                      // GENUINE fake (four valid labels that disagree) and a
                      // LOST TRUTH LINK (bestTkIdx arbitration gave a hit
                      // mcTrackID -1, so a perfectly good quad cannot be
                      // recognised).  Only the first is the algorithm's fault.
                      const int nbad = (la < 0) + (lb < 0) + (lc < 0) + (ld < 0);
                      const int eb = etabin(std::asinh(cot_s));
                      B_pur_all[eb]++;
                      if (nbad > 0) {
                        B_n_quad_nolink++;
                        B_pur_nol[eb]++;
                      } else if (!(la == lb && la == lc && la == ld)) {
                        B_n_quad_fake++;
                        B_pur_fake[eb]++;
                      } else {
                        B_pur_true[eb]++;
                      }
                      // how many of the four agree with the majority label
                      int best = 0;
                      const int ls[4] = {la, lb, lc, ld};
                      for (int u = 0; u < 4; ++u) {
                        if (ls[u] < 0)
                          continue;
                        int cnt = 0;
                        for (int w = 0; w < 4; ++w)
                          if (ls[w] == ls[u])
                            ++cnt;
                        best = std::max(best, cnt);
                      }
                      if (best >= 3)
                        B_n_quad_3of4++;
                    }
                    if (true_trip && ld == la) {
                      B_n_quad_true++;
                      B_dphi.push_back(dphi_d);
                      B_dzd.push_back(dz_d);
                      {
                        // sigma_r of hit d from its covariance, as
                        // HitStructures.cc hit_r_half_extent_of() without hl_fac
                        const Hit &h = (*hd)[id];
                        const double x = h.x(), y = h.y(), r2 = x * x + y * y;
                        const double vr = (x * x * h.exx() + 2 * x * y * h.exy() + y * y * h.eyy()) / r2;
                        B_q4_dz.push_back(dz_d);
                        B_q4_csr.push_back(std::abs(cot_s) * std::sqrt(vr > 0 ? vr : 0.0));
                        B_q4_pt.push_back(0.003 * kBfield * R);
                        B_q4_sz.push_back(std::sqrt(h.ezz()));
                        B_q4_cot.push_back(std::abs(cot_s));
                        B_q4_sr.push_back(std::sqrt(vr > 0 ? vr : 0.0));
                        // the hit's own contribution to the residual z - zpred(r_hit), EXACT:
                        // var(z - cot r) = var_z - 2 cot cov_rz + cot^2 var_r, with the
                        // SIGNED r-z covariance (a tilted sensor moves r and z together)
                        const auto &E = h.error();
                        const double r = std::sqrt(r2);
                        const double crz = (x * E.At(0, 2) + y * E.At(1, 2)) / r;
                        const double vres = h.ezz() - 2 * cot_s * crz + cot_s * cot_s * vr;
                        B_q4_sres.push_back(std::sqrt(vres > 0 ? vres : 0.0));
                        const int mid = h.mcHitID();
                        if (mid >= 0 && mid < (int)ev->simHitStates_.size() && ev->simHitStates_[mid].is_valid()) {
                          const SimHitState &t = ev->simHitStates_[mid];
                          const double rt = std::hypot(t.x(), t.y());
                          // both halves evaluated at the truth radius; their difference is dz_d
                          const double zp_t = zd2 + cot_s * (rt - rrd);
                          B_q4_epred.push_back(zp_t - t.z());
                          B_q4_ehit.push_back((h.z() - t.z()) - cot_s * (rrd - rt));
                          B_q4_tcot.push_back(std::abs(cot_s));
                          B_q4_tpt.push_back(0.003 * kBfield * R);
                          B_q4_ttilt.push_back(std::sqrt(vr > 0 ? vr : 0.0) > 0.2 * std::sqrt(h.ezz()) ? 1 : 0);
                        }
                      }
                      B_zv_true.push_back(zv);
                      if (auto it = found_true.find(la); it != found_true.end())
                        it->second = true;
                      record_truth(ev, la, zv, R);
                    }
                  }
                  }
                }
                }
              }
            }
          }
        }
      }
      B_t_search += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_search0).count();
      for (const auto &kv : found_true)
        if (kv.second) {
          B_n_trip_found++;
          B_eff_pt_num.push_back(ev->simTracks_[kv.first].pT());
          B_eff_eta_num.push_back(ev->simTracks_[kv.first].momEta());
        }
    }
  }  // namespace

  //============================================================================
  // diagnostic: what do the production points actually look like?
  //============================================================================
  void sg_diag(Event *ev) {
    const TrackVec &st = ev->simTracks_;
    printf("\n[sg-diag] %zu sim tracks\n", st.size());
    printf("  beamSpot x %.5f y %.5f z %.5f sigmaZ %.3f wx %.2e wy %.2e\n",
           ev->beamSpot_.x, ev->beamSpot_.y, ev->beamSpot_.z, ev->beamSpot_.sigmaZ,
           ev->beamSpot_.beamWidthX, ev->beamSpot_.beamWidthY);
    std::vector<double> dr;
    dr.reserve(st.size());
    for (const auto &t : st) dr.push_back(std::hypot(t.x(), t.y()));
    std::vector<double> d2 = dr;
    printf("  per-TRACK production radius: p10 %.5f  p50 %.5f  p90 %.3f  p99 %.3f cm\n",
           quant(d2,0.10), quant(d2,0.50), quant(d2,0.90), quant(d2,0.99));
    for (double cut : {1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 1.0}) {
      long n = 0; std::map<std::tuple<float,float,float>,int> g;
      for (const auto &t : st) if (std::hypot(t.x(), t.y()) < cut) { ++n; g[std::make_tuple(t.x(),t.y(),t.z())]++; }
      printf("   dr < %8.5f cm : %7ld tracks in %6zu distinct points\n", cut, n, g.size());
    }
    // group everything, list the biggest
    std::map<std::tuple<float,float,float>,int> g;
    for (const auto &t : st) g[std::make_tuple(t.x(),t.y(),t.z())]++;
    std::vector<std::pair<int,std::tuple<float,float,float>>> v;
    for (auto &kv : g) v.push_back({kv.second, kv.first});
    std::sort(v.begin(), v.end(), [](auto&a, auto&b){ return a.first > b.first; });
    printf("  top 15 production points by multiplicity:\n");
    for (int i = 0; i < 15 && i < (int)v.size(); ++i)
      printf("    n=%5d  x %+.6f  y %+.6f  z %+.5f   r %.6f\n", v[i].first,
             std::get<0>(v[i].second), std::get<1>(v[i].second), std::get<2>(v[i].second),
             std::hypot(std::get<0>(v[i].second), std::get<1>(v[i].second)));
    // Signal tracks only
    std::map<std::tuple<float,float,float>,int> gs;
    for (const auto &t : st) if (t.prodType() == Track::ProdType::Signal) gs[std::make_tuple(t.x(),t.y(),t.z())]++;
    printf("  SIGNAL tracks: %zu distinct production points\n", gs.size());
    int i = 0;
    for (auto &kv : gs) { if (i++ > 9) break;
      printf("    n=%5d  x %+.6f  y %+.6f  z %+.5f   r %.6f\n", kv.second,
             std::get<0>(kv.first), std::get<1>(kv.first), std::get<2>(kv.first),
             std::hypot(std::get<0>(kv.first), std::get<1>(kv.first))); }
  }


  // Radial structure of a layer: is the ~1 cm thickness discrete (ladders) or
  // continuous?  That decides whether a per-hit radial sub-layer index can be
  // made a bin axis and remove the |cot theta| * dr_half over-scan.
  void sg_rdiag(Event *ev) {
    for (int L = 0; L <= 4; ++L) {
      if (L >= (int)ev->layerHits_.size())
        break;
      const HitVec &hv = ev->layerHits_[L];
      if (hv.empty())
        continue;
      std::vector<double> r;
      r.reserve(hv.size());
      double lo = 1e9, hi = -1e9;
      for (const auto &h : hv) {
        double rr = h.r();
        r.push_back(rr);
        lo = std::min(lo, rr);
        hi = std::max(hi, rr);
      }
      const int NB = 60;
      std::vector<int> hh(NB, 0);
      for (double rr : r) {
        int b = (int)((rr - lo) / (hi - lo + 1e-9) * NB);
        hh[std::min(NB - 1, std::max(0, b))]++;
      }
      int mx = *std::max_element(hh.begin(), hh.end());
      printf("\n  layer %d: %zu hits, r %.3f - %.3f cm  (bin %.4f cm)\n", L, hv.size(), lo, hi, (hi - lo) / NB);
      for (int b = 0; b < NB; ++b) {
        int w = mx ? (60 * hh[b]) / mx : 0;
        printf("    %7.3f |%s %d\n", lo + (b + 0.5) * (hi - lo) / NB, std::string(w, '#').c_str(), hh[b]);
      }
    }
  }

  void sg_add_event(Event *ev) {
    if (ev == nullptr)
      return;
    part_a(ev);
    if (g_do_partb) {
      const auto t0 = std::chrono::steady_clock::now();
      part_b(ev);
      B_t_partb += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      vertex_find();
    }
    if (g_verbose)
      printf("[sg] event done: %d sim tracks\n", (int)ev->simTracks_.size());
  }

  //============================================================================
  // report
  //============================================================================
  void sg_report() {
    printf("\n================ SeedGeom report ================\n");
    printf("config: pT_min %.2f GeV (R_min %.1f cm)   D0_max %.3f cm   layers (%d,%d,%d)   qwin %.3f cm\n",
           g_pt_min, r_min_cm(), g_d0_max, g_la, g_lb, g_lc, g_qwin);

    //---- PART A
    printf("\n-- PART A: truth vertices from sim-track production points --\n");
    if (A_n_events == 0) {
      printf("   no events\n");
    } else {
      double nv = 0, np = 0;
      for (int v : A_nvtx_per_ev)
        nv += v;
      for (int v : A_nprim_per_ev)
        np += v;
      printf("   events %d, sim tracks %.0f/ev\n", A_n_events, (double)A_n_simtracks / A_n_events);
      printf("   beamspot (x,y) = (%.4f, %.4f) cm, sigmaZ %.2f cm\n", A_bs_x, A_bs_y, A_bs_sigz);
      double np09 = 0;
      for (int v : A_nprim09_per_ev)
        np09 += v;
      printf("   distinct production points   %8.1f /ev\n", nv / A_n_events);
      printf("   ... on the beamline (dr<%.0f um) %5.1f /ev  (points, pre-clustering)\n",
             1e4 * g_bl_rcut, (double)A_n_onbeam_points / A_n_events);
      printf("   VERTICES after z-clustering (gap %.0f um)  %6.1f /ev\n", 1e4 * g_vtx_gap, np / A_n_events);
      printf("   ... with >= 2 tracks of pT > 0.9      %6.1f /ev\n", np09 / A_n_events);
      printf("   off-beamline points (secondary) %8.1f /ev   <- displaced sample\n", (nv - A_n_onbeam_points) / A_n_events);
      printf("   ProdType  Signal %ld   InTimePU %ld   OutOfTimePU %ld   NotSet %ld  (tracks)\n",
             A_prod_counts[1], A_prod_counts[2], A_prod_counts[3], A_prod_counts[0]);
      printf("   vertices containing a Signal track: %.2f /ev\n", (double)A_n_signal_vtx / A_n_events);
      {
        std::vector<double> z = A_prim_z;
        double m = med(z), s = iqr_sigma(A_prim_z);
        printf("   vertex z: median %+.3f cm, robust sigma %.3f cm\n", m, s);
        printf("   vertex multiplicity: median %.0f all, %.0f with pT>0.9\n",
               med(A_prim_mult), med(A_prim_mult09));
        printf("   nearest-neighbour dz: p10 %.4f  median %.4f  p90 %.4f cm\n",
               quant(A_prim_dz, 0.1), med(A_prim_dz), quant(A_prim_dz, 0.9));
      }
      if (!A_sec_r.empty())
        printf("   secondary production radius: median %.3f  p90 %.3f  p99 %.3f cm\n",
               med(A_sec_r), quant(A_sec_r, 0.9), quant(A_sec_r, 0.99));
    }

    //---- PART B
    if (!g_do_partb || B_n_events == 0) {
      printf("\n-- PART B: not run --\n================================================\n");
      return;
    }
    const double ne = B_n_events;
    printf("\n-- PART B: geometric triplets, layers (%d,%d,%d) --\n", g_la, g_lb, g_lc);
    printf("   events %d;  hits/ev  a %.0f  b %.0f  c %.0f\n", B_n_events, B_n_a / ne, B_n_b / ne, B_n_c / ne);
    printf("   grid: nphi %d  nz %d  nr %d   (cells/layer %d)\n", g_nphi, g_nz, g_nr, g_nphi * g_nz * g_nr);
    printf("   part_b wall clock                 %12.2f ms/ev  (grid build + search + truth)\n", 1e3 * B_t_partb / ne);
    printf("     of which the search loop        %12.2f ms/ev  (%.1f ns per doublet)\n", 1e3 * B_t_search / ne,
           1e9 * B_t_search / std::max(1L, B_n_pairs));
    printf("   doublets (phi window only)        %12.0f /ev\n", B_n_pairs / ne);
    printf("     ... per hit in layer a          %12.2f\n", (double)B_n_pairs / std::max(1L, B_n_a));
    printf("   layer-c cells visited             %12.0f /ev  (%.2f per doublet)\n",
           B_n_cells / ne, (double)B_n_cells / std::max(1L, B_n_pairs));
    printf("   layer-c hits touched              %12.0f /ev  (%.3f per doublet)\n",
           B_n_ctouch / ne, (double)B_n_ctouch / std::max(1L, B_n_pairs));
    printf("   triplets, phi windows only        %12.0f /ev  (%.3f per doublet)\n",
           B_n_trip_phi / ne, (double)B_n_trip_phi / std::max(1L, B_n_pairs));
    printf("   triplets, + q prediction          %12.0f /ev  (%.4f per doublet)\n",
           B_n_trip / ne, (double)B_n_trip / std::max(1L, B_n_pairs));
    printf("     >> q prediction reduction factor  %10.1f x\n",
           (double)B_n_trip_phi / std::max(1L, B_n_trip));
    printf("   of those, all-3-hits-one-simtrack %12.0f /ev  -> purity %.4f\n",
           B_n_trip_true / ne, (double)B_n_trip_true / std::max(1L, B_n_trip));
    if (g_ld >= 0) {
      printf("   -- 4th layer (%d), the helix is DETERMINED so this is a lookup --\n", g_ld);
      printf("   layer-d cells visited             %12.0f /ev  (%.2f per triplet)\n",
             B_n_cells_d / ne, (double)B_n_cells_d / std::max(1L, B_n_trip));
      printf("   layer-d hits touched              %12.0f /ev  (%.3f per triplet)\n",
             B_n_dtouch / ne, (double)B_n_dtouch / std::max(1L, B_n_trip));
      printf("   QUADRUPLETS                       %12.0f /ev  (%.4f per triplet)\n",
             B_n_quad / ne, (double)B_n_quad / std::max(1L, B_n_trip));
      printf("     >> 4th-layer rejection factor     %10.1f x    purity %.4f (true %.0f /ev)\n",
             (double)B_n_trip / std::max(1L, B_n_quad), (double)B_n_quad_true / std::max(1L, B_n_quad),
             B_n_quad_true / ne);
      printf("   quad impurity split: genuine fake (4 valid labels, disagree) %6.0f /ev (%.4f)\n",
             B_n_quad_fake / ne, (double)B_n_quad_fake / std::max(1L, B_n_quad));
      printf("                        >= 1 hit with NO truth link              %6.0f /ev (%.4f)\n",
             B_n_quad_nolink / ne, (double)B_n_quad_nolink / std::max(1L, B_n_quad));
      printf("                        >= 3 of 4 hits share a label             %6.0f /ev (%.4f)\n",
             B_n_quad_3of4 / ne, (double)B_n_quad_3of4 / std::max(1L, B_n_quad));
      {
        printf("   hits carrying a truth link, per layer: ");
        for (int i = 0; i < 4; ++i)
          if (B_n_hit_tot[i])
            printf("%.4f ", (double)B_n_hit_lbl[i] / B_n_hit_tot[i]);
        printf("\n");
        // Do NOT divide the purity by prod(link fraction): most unlabelled
        // hits belong to no sim track at all (noise, sub-threshold secondaries)
        // and could never have been part of a true quad, so that model gives
        // a "corrected purity" above 1.  The decidable ratio is the honest one.
        const long dec = B_n_quad_true + B_n_quad_fake;
        printf("   DECIDABLE purity (true / (true + genuine fake)) = %.4f   [%.4f undecidable]\n",
               (double)B_n_quad_true / std::max(1L, dec),
               (double)B_n_quad_nolink / std::max(1L, B_n_quad));
      }
      // BIAS vs WIDTH.  A wider prediction just costs combinatorics; a biased
      // one is fatal, and uniform B + no material are both bias sources.  The
      // median is the test, in units of its own sigma.
      {
        const double mp = med(B_dphi), sp = iqr_sigma(B_dphi);
        const double mz = med(B_dzd), sz = iqr_sigma(B_dzd);
        printf("   4th-layer residuals, TRUE quads (n %zu):\n", B_dphi.size());
        printf("      dphi  median %+10.6f rad  sigma %.6f  -> bias %+.3f sigma\n",
               mp, sp, sp > 0 ? mp / sp : 0.0);
        printf("      dz    median %+10.6f cm   sigma %.6f  -> bias %+.3f sigma\n",
               mz, sz, sz > 0 ? mz / sz : 0.0);
      }
      // q width vs |cot| sigma_r (centroid-radius term) and vs pT (MS goes as 1/p)
      if (!B_q4_dz.empty()) {
        const double ce[] = {0, 0.002, 0.005, 0.01, 0.02, 0.04, 1e9};
        const double pe[] = {0.9, 1.5, 3.0, 10.0, 1e9};
        constexpr int nc = 6, np = 4;
        printf("   4th-layer q residual of TRUE quads, robust sigma [um] (n), by |cot| sigma_r(hit) and pT:\n");
        printf("   |cot|sig_r [cm]  med[um]");
        for (int j = 0; j < np; ++j)
          printf("   pT %4.1f-%-5.1f", pe[j], pe[j + 1] > 1e8 ? 99.0 : pe[j + 1]);
        printf("       all pT   hit sig_z[um]   pull sigma\n");
        for (int i = 0; i < nc; ++i) {
          std::vector<double> csr_in, sz_in, pull;
          std::vector<double> all;
          std::vector<std::vector<double>> cell(np);
          for (size_t k = 0; k < B_q4_dz.size(); ++k) {
            if (B_q4_csr[k] < ce[i] || B_q4_csr[k] >= ce[i + 1])
              continue;
            csr_in.push_back(B_q4_csr[k]);
            sz_in.push_back(B_q4_sz[k]);
            // pull against the hit's own sigma_z plus the centroid-radius term alone
            pull.push_back(B_q4_dz[k] / std::hypot(B_q4_sz[k], B_q4_csr[k]));
            all.push_back(B_q4_dz[k]);
            for (int j = 0; j < np; ++j)
              if (B_q4_pt[k] >= pe[j] && B_q4_pt[k] < pe[j + 1])
                cell[j].push_back(B_q4_dz[k]);
          }
          printf("   %5.3f-%-6.3f   %7.0f", ce[i], ce[i + 1] > 1e8 ? 9.999 : ce[i + 1],
                 csr_in.empty() ? 0.0 : 1e4 * med(csr_in));
          for (int j = 0; j < np; ++j)
            if (cell[j].size() >= 20)
              printf("   %6.0f (%5zu)", 1e4 * iqr_sigma(cell[j]), cell[j].size());
            else
              printf("   %14s", "-");
          if (all.size() >= 20)
            printf("   %6.0f (%5zu)   %9.0f   %9.3f\n", 1e4 * iqr_sigma(all), all.size(), 1e4 * med(sz_in),
                   iqr_sigma(pull));
          else
            printf("   %14s\n", "-");
        }
      }
      // The along-sensor error of a tilted module moves r AND z together, so it
      // enters the residual at the hit's own r as u (cos tau - cot sin tau), NOT in
      // quadrature: for a module facing the IP the two cancel.  sigma_u, tau from the
      // hit's own sigma_z and sigma_r (the along-sensor term dominates both for P).
      if (!B_q4_dz.empty()) {
        const double ke[] = {0, 0.25, 0.5, 0.75, 1.0, 1.4};
        printf("   same, by |cot theta| and module tilt (tan tau = sig_r/sig_z of the hit):\n");
        printf("   pred = median of the hit's own sigma of (z - cot r), exact from its covariance\n");
        printf("   |cot|        flat: sigma[um] (n)  pred[um]     tilted: sigma[um] (n)  pred[um]  med tau[deg]\n");
        for (int i = 0; i < 5; ++i) {
          std::vector<double> d[2], pr[2], tau_t;
          for (size_t k = 0; k < B_q4_dz.size(); ++k) {
            if (B_q4_cot[k] < ke[i] || B_q4_cot[k] >= ke[i + 1])
              continue;
            const double sz = B_q4_sz[k], sr = B_q4_sr[k], su = std::hypot(sz, sr);
            const double tau = std::atan2(sr, sz);
            const int t = (sr > 0.2 * sz) ? 1 : 0;
            d[t].push_back(B_q4_dz[k]);
            (void)su;
            pr[t].push_back(B_q4_sres[k]);
            if (t)
              tau_t.push_back(tau * 180 / kPi);
          }
          printf("   %4.2f-%-4.2f", ke[i], ke[i + 1]);
          for (int t = 0; t < 2; ++t)
            if (d[t].size() >= 20)
              printf("       %6.0f (%5zu)  %7.0f", 1e4 * iqr_sigma(d[t]), d[t].size(), 1e4 * med(pr[t]));
            else
              printf("       %22s", "-");
          printf("   %8.1f\n", tau_t.empty() ? 0.0 : med(tau_t));
        }
      }
      if (!B_q4_ehit.empty()) {
        const double ke[] = {0, 0.25, 0.5, 0.75, 1.0, 1.4};
        const double pe[] = {0.9, 1.5, 3.0, 10.0, 1e9};
        printf("   TRUTH SPLIT at the true crossing radius (n %zu of %zu true quads have a state), robust sigma [um]:\n",
               B_q4_ehit.size(), B_q4_dz.size());
        printf("   |cot|        flat: e_hit  e_pred (n)          tilted: e_hit  e_pred (n)\n");
        for (int i = 0; i < 5; ++i) {
          std::vector<double> eh[2], ep[2];
          for (size_t k = 0; k < B_q4_ehit.size(); ++k)
            if (B_q4_tcot[k] >= ke[i] && B_q4_tcot[k] < ke[i + 1]) {
              const int t = (int)B_q4_ttilt[k];
              eh[t].push_back(B_q4_ehit[k]);
              ep[t].push_back(B_q4_epred[k]);
            }
          printf("   %4.2f-%-4.2f", ke[i], ke[i + 1]);
          for (int t = 0; t < 2; ++t)
            if (eh[t].size() >= 20)
              printf("       %6.0f %6.0f (%5zu)", 1e4 * iqr_sigma(eh[t]), 1e4 * iqr_sigma(ep[t]), eh[t].size());
            else
              printf("       %20s", "-");
          printf("\n");
        }
        printf("   by pT (all |cot|):   ");
        for (int j = 0; j < 4; ++j) {
          std::vector<double> eh, ep;
          for (size_t k = 0; k < B_q4_ehit.size(); ++k)
            if (B_q4_tpt[k] >= pe[j] && B_q4_tpt[k] < pe[j + 1]) {
              eh.push_back(B_q4_ehit[k]);
              ep.push_back(B_q4_epred[k]);
            }
          if (eh.size() >= 20)
            printf("  pT %.1f-%.1f: e_hit %4.0f e_pred %4.0f (%zu)", pe[j], pe[j + 1] > 1e8 ? 99. : pe[j + 1],
                   1e4 * iqr_sigma(eh), 1e4 * iqr_sigma(ep), eh.size());
        }
        printf("\n");
      }
      printf("                        ALL touched: dphi sigma %.5f rad   dz sigma %.5f cm\n",
             iqr_sigma(B_dphi_all), iqr_sigma(B_dzd_all));
    }
    printf("   sim tracks with a hit in all %d      %10.0f /ev  (any pT)\n",
           g_ld >= 0 ? 4 : 3, B_n_avail_all / ne);
    printf("   ... FINDABLE (pT>%.2f, prod r<%.3f) %10.0f /ev  <- the honest denominator\n",
           g_pt_min, g_d0_max, B_n_trip_avail / ne);
    printf("     of which a true one is built    %12.0f /ev  -> efficiency %.4f\n",
           B_n_trip_found / ne, (double)B_n_trip_found / std::max(1L, B_n_trip_avail));

    {
      const double ptd[] = {0.9, 1.2, 1.6, 2.2, 3.2, 5.0, 1e9};
      const char *ptn[] = {"0.9-1.2", "1.2-1.6", "1.6-2.2", "2.2-3.2", "3.2-5  ", ">5     "};
      const char *etn[] = {"0-0.4", ".4-.8", ".8-1.2", "1.2-1.6", "1.6-2.2", ">2.2"};
      printf("\n   EFFICIENCY, pT x |eta|   (n in parentheses; blank = fewer than 20)\n");
      printf("   %-9s", "pT \\ eta");
      for (int j = 0; j < NEB; ++j)
        printf("%14s", etn[j]);
      printf("%14s\n", "all eta");
      for (int i = 0; i < 6; ++i) {
        printf("   %-9s", ptn[i]);
        for (int j = 0; j <= NEB; ++j) {
          long n = 0, d = 0;
          for (size_t k = 0; k < B_eff_pt_den.size(); ++k)
            if (B_eff_pt_den[k] >= ptd[i] && B_eff_pt_den[k] < ptd[i + 1] &&
                (j == NEB || etabin(B_eff_eta_den[k]) == j))
              ++d;
          for (size_t k = 0; k < B_eff_pt_num.size(); ++k)
            if (B_eff_pt_num[k] >= ptd[i] && B_eff_pt_num[k] < ptd[i + 1] &&
                (j == NEB || etabin(B_eff_eta_num[k]) == j))
              ++n;
          if (d >= 20)
            printf("%9.3f(%4ld)", (double)n / d, d);
          else
            printf("%14s", "-");
        }
        printf("\n");
      }
      printf("   %-9s", "all pT");
      for (int j = 0; j <= NEB; ++j) {
        long n = 0, d = 0;
        for (size_t k = 0; k < B_eff_eta_den.size(); ++k)
          if (j == NEB || etabin(B_eff_eta_den[k]) == j)
            ++d;
        for (size_t k = 0; k < B_eff_eta_num.size(); ++k)
          if (j == NEB || etabin(B_eff_eta_num[k]) == j)
            ++n;
        if (d >= 20)
          printf("%9.3f(%4ld)", (double)n / d, d);
        else
          printf("%14s", "-");
      }
      printf("\n");
      printf("\n   PURITY vs |eta| of the CANDIDATE (from its own cot theta, no truth)\n");
      printf("   %-10s %10s %10s %10s %10s %12s\n", "|eta|", "quads/ev", "true", "fake", "undecid", "decidable");
      for (int j = 0; j < NEB; ++j) {
        if (B_pur_all[j] < 20)
          continue;
        const long dec = B_pur_true[j] + B_pur_fake[j];
        printf("   %-10s %10.0f %10.4f %10.4f %10.4f %12.4f\n", etn[j], B_pur_all[j] / ne,
               (double)B_pur_true[j] / B_pur_all[j], (double)B_pur_fake[j] / B_pur_all[j],
               (double)B_pur_nol[j] / B_pur_all[j], (double)B_pur_true[j] / std::max(1L, dec));
      }
    }
    printf("\n   q-prediction residual z_c - z_pred [cm]\n");
    {
      std::vector<double> az;
      az.reserve(B_zres.size());
      for (double v : B_zres)
        az.push_back(std::abs(v));
      printf("     TRUE triplets : median %+.5f  robust sigma %.5f  p90 |.| %.5f  (n %zu)\n",
             med(B_zres), iqr_sigma(B_zres), quant(az, 0.9), B_zres.size());
    }
    printf("     ALL triplets  : median %+.5f  robust sigma %.5f  (n %zu, pre-q-cut)\n",
           med(B_zres_all), iqr_sigma(B_zres_all), B_zres_all.size());
    printf("   3-point longitudinal residual (arc-length test) [cm]\n");
    printf("     all triplets  : median %+.5f  robust sigma %.5f\n", med(B_slres), iqr_sigma(B_slres));

    printf("\n   beamline intercept z_v [cm]\n");
    printf("     all triplets  : robust sigma %.3f  (n %zu)\n", iqr_sigma(B_zv), B_zv.size());
    printf("     true triplets : robust sigma %.3f  (n %zu)\n", iqr_sigma(B_zv_true), B_zv_true.size());
    printf("     >> z_v RESOLUTION vs truth (true triplets, primary tracks):\n");
    printf("        median %+.5f cm   robust sigma %.5f cm = %.0f um   (n %zu)\n",
           med(B_zvres), iqr_sigma(B_zvres), 1e4 * iqr_sigma(B_zvres), B_zvres.size());
    {
      // resolution in pT bands
      const double edges[] = {0.0, 0.7, 1.0, 2.0, 5.0, 1e9};
      const char *nm[] = {"<0.7", "0.7-1", "1-2", "2-5", ">5"};
      for (int i = 0; i < 5; ++i) {
        std::vector<double> v;
        for (size_t k = 0; k < B_zvres.size(); ++k)
          if (B_zvres_pt[k] >= edges[i] && B_zvres_pt[k] < edges[i + 1])
            v.push_back(B_zvres[k]);
        if (v.size() > 20)
          printf("        pT %-6s  n %7zu   sigma %7.1f um\n", nm[i], v.size(), 1e4 * iqr_sigma(v));
      }
    }
    if (!B_ptres.empty())
      printf("   3-point curvature: (pT_fit-pT_true)/pT_true  median %+.3f  robust sigma %.3f\n",
             med(B_ptres), iqr_sigma(B_ptres));
    printf("\n   -- VERTICES FROM THE VOTES ALONE (window +-%.0f um, >= %d votes, match %.0f um) --\n",
           1e4 * g_vf_win, g_vf_min, 1e4 * g_vf_match);
    printf("      truth vertices    %8.1f /ev     found  %8.1f /ev     fake  %8.1f /ev\n",
           V_n_true / ne, V_n_found / ne, V_n_fake / ne);
    printf("      matched           %8.1f /ev  -> EFFICIENCY  %.4f   (fake fraction %.4f)\n",
           V_n_match / ne, (double)V_n_match / std::max(1L, V_n_true),
           (double)V_n_fake / std::max(1L, V_n_found));
    printf("      truth vtx with >= 2 tracks pT>0.9: %6.1f /ev  -> EFFICIENCY  %.4f\n",
           V_n_true09 / ne, (double)V_n_match09 / std::max(1L, V_n_true09));
    printf("      vertex z resolution: robust sigma %.5f cm = %.0f um  (n %zu)\n",
           iqr_sigma(V_dz), 1e4 * iqr_sigma(V_dz), V_dz.size());
    {
      const double ed[] = {0, 1, 2, 4, 8, 1e9};
      const char *nm[] = {"0", "1", "2-3", "4-7", ">=8"};
      printf("      efficiency vs n(pT>0.9) at the vertex:");
      for (int i = 0; i < 5; ++i) {
        long n = 0, d = 0;
        for (double x : V_eff_mult_den)
          if (x >= ed[i] && x < ed[i + 1])
            ++d;
        for (double x : V_eff_mult_num)
          if (x >= ed[i] && x < ed[i + 1])
            ++n;
        if (d > 20)
          printf("  %s: %.3f (n %ld)", nm[i], (double)n / d, d);
      }
      printf("\n");
    }
    printf("================================================\n");
  }

  void sg_write_root(const char *prefix) {
    TFile f(Form("%s.root", prefix), "recreate");
    auto fill = [](const char *n, const char *t, int nb, double lo, double hi, const std::vector<double> &v) {
      TH1D *h = new TH1D(n, t, nb, lo, hi);
      for (double x : v)
        h->Fill(x);
      h->Write();
    };
    fill("prim_z", "truth primary vertex z;z [cm];vertices", 400, -20, 20, A_prim_z);
    fill("prim_dz", "nearest-neighbour vertex spacing;dz [cm];pairs", 400, 0, 2, A_prim_dz);
    fill("sec_r", "secondary production radius;r [cm];vertices", 400, 0, 40, A_sec_r);
    fill("zv", "triplet beamline intercept, all;z_v [cm];triplets", 800, -20, 20, B_zv);
    fill("zv_true", "triplet beamline intercept, true;z_v [cm];triplets", 800, -20, 20, B_zv_true);
    fill("zvres", "z_v(fit) - z_v(truth);dz [cm];triplets", 400, -0.2, 0.2, B_zvres);
    fill("zres", "z_c - z_pred, true triplets;dz [cm];triplets", 400, -0.3, 0.3, B_zres);
    fill("zres_all", "z_c - z_pred, all;dz [cm];triplets", 400, -0.3, 0.3, B_zres_all);
    fill("slres", "3-point longitudinal residual;dz [cm];triplets", 400, -0.2, 0.2, B_slres);
    fill("ptres", "(pT_fit-pT_true)/pT_true;;triplets", 400, -2, 2, B_ptres);
    fill("vtx_dz", "found vertex - truth vertex;dz [cm];vertices", 400, -0.05, 0.05, V_dz);
    f.Close();
    printf("[sg] wrote %s.root\n", prefix);
  }

}  // namespace mkfit

//==============================================================================
// Short global-scope names for the shell.  ACLiC generates a dictionary for
// this file, so these are callable directly as sg_reset(), sg_ev(s.event()),
// ... with no separate shim translation unit.
//==============================================================================

void sg_reset() { mkfit::sg_reset(); }
void sg_ev(mkfit::Event *ev) { mkfit::sg_add_event(ev); }
void sg_report() { mkfit::sg_report(); }
void sg_write(const char *prefix) { mkfit::sg_write_root(prefix); }
void sg_dg(mkfit::Event *ev) { mkfit::sg_diag(ev); }
void sg_rdg(mkfit::Event *ev) { mkfit::sg_rdiag(ev); }

void sg_ptmin(float v) { mkfit::sg_set_ptmin(v); }
void sg_d0max(float v) { mkfit::sg_set_d0max(v); }
void sg_qwin(float v) { mkfit::sg_set_qwin(v); }
void sg_qwind(float v) { mkfit::sg_set_qwin_d(v); }
void sg_phiwind(float v) { mkfit::sg_set_phiwin_d(v); }
void sg_phimargin(float v) { mkfit::sg_set_phi_margin(v); }
void sg_layers(int a, int b, int c) { mkfit::sg_set_layers(a, b, c); }
void sg_layer4(int d) { mkfit::sg_set_layer4(d); }
void sg_nr(int n) { mkfit::sg_set_nr(n); }
void sg_philin(int b, float m) { mkfit::sg_set_phi_lin(b, m); }
void sg_grid(int nphi, int nz) { mkfit::sg_set_grid(nphi, nz); }
void sg_light(bool b) { mkfit::sg_set_light(b); }
void sg_partb(bool b) { mkfit::sg_set_do_partb(b); }
void sg_verbose(bool b) { mkfit::sg_set_verbose(b); }
void sg_vtxgap(float v) { mkfit::sg_set_vtx_gap(v); }
void sg_blrcut(float v) { mkfit::sg_set_bl_rcut(v); }
void sg_vf(float w, int n, float m) { mkfit::sg_set_vf(w, n, m); }
