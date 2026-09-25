#ifndef mkfit_seeding_SeedLayer_h
#define mkfit_seeding_SeedLayer_h

// A layer of hits binned by the mkFit binnor, with the axis types as template
// parameters.  The registration is LayerOfHits::suckInHits() verbatim; what
// differs is that LayerOfHits fixes its axis types as typedefs, while the seeder
// needs to choose them per layer ROLE: the doublet target is queried in phi
// alone, the third and fourth layers in (phi, q).
//
// Hits are stored as struct-of-arrays in BIN ORDER, i.e. the binnor's m_ranks
// order: q N-bin major, phi fine bin minor.  So for one q N-bin, a range of phi
// N-bins is one contiguous run, and start_[] is a CSR over the N-bins in that
// order.

#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/binnor.h"

#include <cmath>
#include <vector>

namespace mkfit::seeding {

  constexpr float kPi = 3.14159265358979323846f;

  template <typename AxPhi, typename AxQ>
  class SeedLayer {
  public:
    using binnor_t = binnor<unsigned int, AxPhi, AxQ, 18, 14>;

    SeedLayer(float qmin, float qmax, unsigned int nq)
        : ax_phi_(-kPi, kPi), ax_q_(qmin, qmax, nq), binnor_(ax_phi_, ax_q_, true, false) {}

    void fill(const HitVec &hits) {
      n_ = hits.size();
      binnor_.reset_contents();
      binnor_.begin_registration(n_);
      for (unsigned int i = 0; i < n_; ++i)
        binnor_.register_entry_safe(hits[i].phi(), hits[i].z());
      binnor_.finalize_registration();

      phi_.resize(n_);
      z_.resize(n_);
      r_.resize(n_);
      invr_.resize(n_);
      x_.resize(n_);
      y_.resize(n_);
      orig_.resize(n_);
      rlo_ = 1e9;
      rhi_ = -1e9;
      double rsum = 0;
      for (unsigned int i = 0; i < n_; ++i) {
        const unsigned int j = binnor_.m_ranks[i];
        const Hit &h = hits[j];
        phi_[i] = h.phi();
        z_[i] = h.z();
        r_[i] = h.r();
        invr_[i] = 1.0f / r_[i];
        // the finder's own expression, float r times float cos(phi), so a
        // triplet reads the same value it would compute
        x_[i] = r_[i] * std::cos(phi_[i]);
        y_[i] = r_[i] * std::sin(phi_[i]);
        orig_[i] = j;
        // the same double reduction the prototype's rminmax() does
        rlo_ = std::min(rlo_, (double)r_[i]);
        rhi_ = std::max(rhi_, (double)r_[i]);
        rsum += r_[i];
      }
      rmean_ = rsum / std::max(1u, n_);

      const unsigned int nb = n_phi_bins() * n_q_bins();
      start_.assign(nb + 1, 0);
      for (unsigned int k = 0; k < nb; ++k)
        start_[k + 1] = start_[k] + binnor_.m_bins[k].count;
    }

    // [begin, end) of the hits in q N-bin qi and phi N-bins [p1, p2), NO wrap.
    unsigned int run_begin(unsigned int qi, unsigned int p1) const { return start_[qi * n_phi_bins() + p1]; }
    unsigned int run_end(unsigned int qi, unsigned int p2) const { return start_[qi * n_phi_bins() + p2]; }

    // Calls f(i) for every hit in phi N-bins [p.begin, p.end) (half-open, may
    // wrap) and q N-bins [q.begin, q.end).  A range with begin == end on the
    // circle is EMPTY: callers clamp their half-width below pi first.
    template <typename F>
    void for_each_in(typename AxPhi::I_pair p, typename AxQ::I_pair q, F &&f) const {
      const unsigned int np = n_phi_bins();
      for (unsigned int qi = q.begin; qi < q.end; ++qi) {
        if (p.begin < p.end || p.end == 0) {
          const unsigned int e = p.end == 0 ? np : p.end;
          for (unsigned int i = run_begin(qi, p.begin); i < run_end(qi, e); ++i)
            f(i);
        } else if (p.begin > p.end) {
          for (unsigned int i = run_begin(qi, p.begin); i < run_end(qi, np); ++i)
            f(i);
          for (unsigned int i = run_begin(qi, 0); i < run_end(qi, p.end); ++i)
            f(i);
        }
      }
    }

    // As for_each_in(), but hands over whole contiguous runs, f(begin, end).
    template <typename F>
    void for_each_run(typename AxPhi::I_pair p, typename AxQ::I_pair q, F &&f) const {
      const unsigned int np = n_phi_bins();
      for (unsigned int qi = q.begin; qi < q.end; ++qi) {
        if (p.begin < p.end || p.end == 0) {
          f(run_begin(qi, p.begin), run_end(qi, p.end == 0 ? np : p.end));
        } else if (p.begin > p.end) {
          f(run_begin(qi, p.begin), run_end(qi, np));
          f(run_begin(qi, 0), run_end(qi, p.end));
        }
      }
    }

    // lo, hi may lie outside (-pi, pi]; wrap them first, since the axis floors
    // (r - R_min) * fac straight into an unsigned bin index.
    typename AxPhi::I_pair phi_range(float lo, float hi) const {
      return ax_phi_.from_R_minmax_to_N_bins(wrap(lo), wrap(hi));
    }
    static float wrap(float p) {
      if (p > kPi)
        p -= 2 * kPi;
      else if (p < -kPi)
        p += 2 * kPi;
      return p;
    }
    typename AxQ::I_pair q_range(float lo, float hi) const { return ax_q_.from_R_minmax_to_N_bins(lo, hi); }
    typename AxQ::I_pair q_all() const { return typename AxQ::I_pair(0, n_q_bins()); }

    unsigned int n_phi_bins() const { return ax_phi_.size_of_N(); }
    unsigned int n_q_bins() const { return ax_q_.size_of_N(); }
    float q_min() const { return ax_q_.m_R_min; }
    float q_max() const { return ax_q_.m_R_max; }
    unsigned int n() const { return n_; }

    const binnor_t &binnor_ref() const { return binnor_; }

    // Fixed-point copies for the int16 pair tests (branch mkfit-seeding-int16),
    // filled by prep_i16() once per event:
    //   phi16_  phi with 2 pi = 2^16, so a difference wraps by integer overflow
    //   r16_    r in units of 2^-10 cm (int16, so r < 32 cm)
    //   gin16_, gout16_  the two per-hit terms of the phi band, in phi16 units:
    //           w(r_in, r_out) = (r_out - r_in)/2R + D0 (1/r_in - 1/r_out) + marg
    //                          = gout(r_out) + gin(r_in) + marg,  exactly
    static constexpr float kPhi16 = 32768.0f / kPi;
    static constexpr float kR16 = 1024.0f;
    std::vector<unsigned short> phi16_;
    std::vector<short> r16_, gin16_, gout16_;
    void prep_i16(float inv2R, float d0m) {
      phi16_.resize(n_);
      r16_.resize(n_);
      gin16_.resize(n_);
      gout16_.resize(n_);
      for (unsigned int i = 0; i < n_; ++i) {
        phi16_[i] = (unsigned short)(int)std::lrint(phi_[i] * kPhi16);
        r16_[i] = (short)std::lrint(r_[i] * kR16);
        gin16_[i] = (short)std::lrint((d0m * invr_[i] - r_[i] * inv2R) * kPhi16);
        gout16_[i] = (short)std::lrint((r_[i] * inv2R - d0m * invr_[i]) * kPhi16);
      }
    }

    // struct-of-arrays, bin order
    std::vector<float> phi_, z_, r_, invr_, x_, y_;
    std::vector<unsigned int> orig_;
    double rlo_ = 0, rhi_ = 0, rmean_ = 0;

  private:
    AxPhi ax_phi_;
    AxQ ax_q_;
    binnor_t binnor_;
    std::vector<unsigned int> start_;
    unsigned int n_ = 0;
  };

}  // namespace mkfit::seeding

#endif
