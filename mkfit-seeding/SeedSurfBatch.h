#ifndef mkfit_seeding_SeedSurfBatch_h
#define mkfit_seeding_SeedSurfBatch_h

// The feed-forward chain in float, batched, is mkfit::SeedChainFinder in MkFitCore since 2026-10-01
// (src/SeedChainFinder.h, with the description that was here). What stays here is the name it had,
// and --chain-fast-check, which compares its stage d prediction with surf::Helix in double.

#include "SeedSurf.h"
#include "RecoTracker/MkFitCore/src/SeedChainFinder.h"

#include <cmath>
#include <vector>

namespace mkfit::seeding {

  namespace surfb = mkfit::seedchain;
  using SurfChainBatch = SeedChainFinder;

  // --chain-fast-check: the float single-point prediction against surf::Helix::predict in double
  inline void surf_check_d(const SeedCand &c, const std::vector<const SurfLayer *> &lay, bool disc, float u, bool ok, float px,
                      float py, float qp, float wq, float wp, SurfFastCheck &CK) {
    const surf::Helix H(surf_p3(*lay[c.pos[0]], c.k[0]), surf_p3(*lay[c.pos[1]], c.k[1]), surf_p3(*lay[c.pos[2]], c.k[2]));
    double qe, pe;
    const bool oke = H.ok && H.predict(disc, u, qe, pe);
    if (oke != ok) {
      ++CK.b_fail_mismatch;
      return;
    }
    if (!ok)
      return;
    ++CK.b_cands;
    const double dq = std::abs(qp - qe), dp = std::abs(surf::wrap(std::atan2((double)py, (double)px) - pe));
    CK.b_max_dq = std::max(CK.b_max_dq, dq), CK.b_max_dphi = std::max(CK.b_max_dphi, dp);
    CK.b_max_rel_q = std::max(CK.b_max_rel_q, dq / wq), CK.b_max_rel_phi = std::max(CK.b_max_rel_phi, dp / wp);
  }

}  // namespace mkfit::seeding

#endif
