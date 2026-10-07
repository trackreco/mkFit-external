#ifndef RecoTracker_MkFitCore_standalone_DataFormats_ValStructs_h
#define RecoTracker_MkFitCore_standalone_DataFormats_ValStructs_h

// ===========================================================================
// Val* -- persistent records of the synthetic validation of propagation,
// Kalman operations and fitting.
//
// These deliberately hold PLAIN members only (scalars and fixed float arrays),
// so that:
//   - a file written here is readable in a bare root.exe session after
//     gSystem->Load("libMkFitRootDataFormats.so"), with no mkFit build at all;
//   - RDataFrame sees every member as a column with no per-type boilerplate.
//
// No residual, pull or ratio is ever stored: those are RDF Defines over these
// columns, so a residual can be redefined without regenerating the sample.
// ===========================================================================

// ValClosure -- ONE round trip, i.e. one (track, plane configuration) trial.
//
//   O -> A   conditioning leg, result kept as the input state
//   A -> B   the step under test
//   B -> A'  the return, onto the SAME plane A
//
// PLANE A is at the ORIGIN with its normal along the momentum there, so it is
// perfectly conditioned (n.dx/ds = 1) and the return leg contributes no error
// of its own. Everything measured is then a property of the OUTBOUND step.
//
// PLANE B is placed at a chosen TURN ANGLE from A -- not at a chosen radius.
// The propagator solves for a path length s satisfying n.(x(s) - p) = 0; it
// does not propagate "to a radius", and a radius is not a control: it sets both
// how far the step goes and how squarely the track meets the surface, and
// couples them. The two real controls are the turn angle and the incidence
// angle at B, and both are scanned directly. r_b and z_b are reported so the
// configuration stays picturable, not because they steer anything.
//
// The conditioning leg is not optional. A propagate-to-SURFACE has a rank-5
// jacobian, so its output covariance is degenerate along the plane normal; a
// round trip started from a full-rank 6x6 cannot close however correct the
// propagator is. Starting from a state that is itself the output of one
// propagation onto A is also what the real code always feeds in.
//
// With material OFF and a uniform field the propagation is a diffeomorphism,
// so the exact answer for A vs A' is ZERO -- every non-zero value is float
// arithmetic. Units: cm, rad, GeV.
struct ValClosure {
  // ----- identity
  int id = -1;
  int track_id = -1;
  int cfg = -1;  // plane-configuration index, see ValCfgInfo

  // ----- the truth track (duplicated here so no join is ever needed)
  float pt = 0, eta = 0, phi0 = 0;
  float d0 = 0, z0 = 0;  // production vertex offset; 0 = from the origin
  int charge = 0;

  // ----- geometry of this trial (r, z are OUTPUTS here, not controls)
  float dalpha_target = 0;  // requested turn angle A->B  [rad]
  float cosinc_target = 0;  // requested |p^ . n| at B
  float r_a = 0, r_b = 0;   // radius reached at A (=0) and at B  [cm]
  float z_a = 0, z_b = 0;   // z at A (=0) and at B  [cm]
  float pnt_a[3] = {0, 0, 0}, nrm_a[3] = {0, 0, 0};
  float pnt_b[3] = {0, 0, 0}, nrm_b[3] = {0, 0, 0};

  // THE REFERENCE, in double. plane B is constructed to pass through the truth
  // helix's own crossing, so this is the exact answer the propagator should
  // return -- an absolute reference, not a relative one. It must be kept in
  // double: at r ~ 40 cm one float ULP is 3.8e-6 cm, which is the same size as
  // the landing error being measured, so storing it as float would make the
  // measurement report its own rounding.
  double ref_pnt_b[3] = {0, 0, 0};  // exact crossing point at plane B  [cm]
  double ref_s_b = 0;               // exact path length from A to it   [cm]
  float cos_inc_a = 0, cos_inc_b = 0;  // |p_hat . n| at each plane
  float s_ab = 0;                      // truth path length A->B  [cm]
  float dalpha_ab = 0;                 // truth helix angle turned A->B  [rad]

  // ----- propagation status, one per leg (0 = clean)
  int fail_oa = 0, fail_ab = 0, fail_ba = 0;

  // ----- states. par = CCS (x, y, z, 1/pT, phi, theta), err = lower triangle.
  float par_a[6] = {0}, err_a[21] = {0};  // at A, after conditioning
  float par_b[6] = {0}, err_b[21] = {0};  // at B
  float par_c[6] = {0}, err_c[21] = {0};  // back at A'
};

// ValCfgInfo -- what a cfg index means. One entry per configuration, written
// to its own tree so a plot can label itself without a hardcoded legend.
struct ValCfgInfo {
  int cfg = -1;
  float dalpha = 0;   // requested turn angle A->B  [rad]
  float cos_inc = 0;  // requested |p^ . n| at plane B; 1 = plane normal to the track
};

// ValStep -- one propagation step as a REAL track actually presents it.
//
// The synthetic scans choose a turn angle; a detector does not. This records
// what turn angles occur between consecutive hits of reconstructed tracks, so
// that any response measured against turn angle can be weighted by the
// distribution that actually happens. Without it a scan is uniform over a
// variable the experiment is not uniform in.
struct ValStep {
  int   event = -1;
  int   track = -1;      // index in the track collection
  int   step = -1;       // which consecutive-hit step along the track
  int   lay_a = -1, lay_b = -1;
  float pt = 0, eta = 0;
  float r_a = 0, r_b = 0, z_a = 0, z_b = 0;
  float d_perp = 0;      // transverse distance between the two hits [cm]
  float d_3d = 0;        // 3D distance                              [cm]
  float dalpha = 0;      // turn angle of the step                   [rad]
};

// ---------------------------------------------------------------------------
// ValCovXport -- ENSEMBLE covariance transport, one step.
//
// The question: propagate a covariance C_A through the production code path to
// get C_B, and separately propagate 10^4 SAMPLES drawn from C_A and form their
// empirical covariance S. If the transport is right, S == C_B.
//
// Why an ensemble and not a closure test: a round trip cannot see an error made
// identically on both legs, and a marginal-by-marginal comparison cannot see a
// wrong CORRELATION. The ensemble has neither blind spot.
//
// RANK IS THE POINT. Propagating to a surface confines the state to that
// surface, so C_B is rank 5 in 6 dimensions -- position has only 2 free
// directions in the plane. C_B^{-1/2} therefore does not exist, and the
// comparison is made inside C_B's 5-D column space. `lam[]` holds the
// eigenvalues of the whitened empirical covariance there: all 1.000 means the
// transport is right in every direction and correlation, not just on average.
// `null_sig` is the sample spread ALONG C_B's null direction, in units of the
// smallest kept sigma -- it is how a wrong rank would show up.
struct ValCovXport {
  int   cfg = -1;
  int   n_samp = 0;       // samples that propagated without a fail flag
  int   n_fail = 0;
  float dalpha = 0.f;     // turn angle of the step [rad]
  float cos_inc = 0.f;    // |p^ . n| at the target plane
  float pt = 0.f;
  float eta = 0.f;
  float cov_scale = 1.f;  // VARIANCE scale on C_A; the linearity control

  int   rank_cb = 0;      // numerical rank of the transported covariance
  float lam[6] = {0,0,0,0,0,0};        // eigenvalues of C_B^{-1/2} S C_B^{-1/2}
  float ratio[6] = {0,0,0,0,0,0};      // sqrt(S_kk / C_B_kk), marginal, per CCS par
  float mean_pull[6] = {0,0,0,0,0,0};  // (mean_emp - par_B)_k / sqrt(C_B_kk)
  float null_sig = 0.f;   // sample sigma along C_B's null direction / smallest kept sigma
  float lam_min = 0.f, lam_max = 0.f;  // convenience: extremes of lam[0..rank)

  // The SAME comparison, after sliding each sample along its own momentum onto
  // the plane perpendicular to the REFERENCE momentum. This is the curvilinear
  // surface, and it is the surface the transported covariance actually lives
  // on: propagateHelixToPlaneMPlex returns a covariance blind along the
  // momentum, NOT along the target plane's normal, with the surface-crossing
  // term supplied later by jacCurv2Loc. At normal incidence the two surfaces
  // coincide and lam_*_cv == lam_*; off-normal only the _cv pair is a fair test.
  float lam_min_cv = 0.f, lam_max_cv = 0.f;
  int   rank_cv = 0;

  // FULL CHAIN: transport + jacCurv2Loc + the local 2-D projection, as the
  // Kalman update actually sees it. Each sample's landing point is offered to
  // kalmanComputeChi2Plane as a hit with negligible error, against the
  // REFERENCE state and its transported covariance. If that covariance is
  // right, the chi2 are chi2 with 2 d.o.f. -- whose quantiles are closed form,
  // -2*ln(1-p): median 1.3863, p90 4.6052, p99 9.2103, mean 2.
  // This works at ANY incidence, which the bare-transport comparison does not.
  float chi2_med = 0.f, chi2_mean = 0.f, chi2_p90 = 0.f, chi2_p99 = 0.f;
  float chi2_med_ratio = 0.f;   // chi2_med / 1.3863; < 1 means over-covered
  // Direct, decomposition-free check: how far is each propagated sample from
  // the plane it was propagated TO? |n.(x - p)|, cm. A propagation that
  // succeeded must give ~0 here, so this is what says whether a large null_sig
  // is a real defect or an artefact of the whitening.
  float cos_inc_meas = 0.f;  // |p^ . n| ACTUALLY achieved at plane B -- must equal
                             // cos_inc, and is the first thing to check if two
                             // incidences give identical answers.
  float covb_sum = 0.f;      // sum |C_B| elements: a crude fingerprint, so
                             // "did the plane change anything at all" is visible.
  float d_plane_ref = 0.f;   // the mean state's own miss
  float d_plane_med = 0.f, d_plane_p90 = 0.f, d_plane_max = 0.f;

  float cb_spec[6] = {0,0,0,0,0,0};    // eigenvalues of C_B's CORRELATION matrix,
                                       // descending. Dimensionless, sum = 6. The
                                       // rank gap must be visible here, not asserted.
};


// ---------------------------------------------------------------------------
// ValSearchHit -- one scanned hit in the BACKWARD (inward) search, flattened
// out of the MKFIT_TRACE graph so it can be read in a bare root session.
//
// The two questions it exists to answer, on a geometry-verified sample:
//   (1) why is the predicted position, especially in q, so far off?
//   (2) why is the chi2 of MC-matched hits too big?
// Both need the residual AND the error it should be judged against, per hit,
// which is why sigma_q_trk / hit_q_half_len / sigma_phi_trk all travel with it.
struct ValSearchHit {
  int   event = -1;
  int   layer = -1;
  int   hit = -1;
  bool  is_barrel = false;
  bool  is_outward = false;
  bool  mc_match = false;

  // pre-selection residuals, and the two error terms they are cut against
  float dphi = -999.f;
  float dq = -999.f;
  // SIGNED residuals, recomputed here from the predicted point and the hit.
  // The code throws the sign away -- MkFinderV2p2.cc:672 does
  // `ddq = std::abs(q - hit_q)` before the cut and before the trace -- so a
  // systematic OFFSET in q is invisible to both, and shows up only as width.
  // prediction minus hit, so positive = the track is predicted at larger q.
  float dq_s = -999.f;
  float dphi_s = -999.f;
  float theta = -999.f;   // the candidate's theta, to test a theta bias directly
  // |module normal . zhat|. 0 = untilted (normal purely radial), rising to ~1
  // for a strongly tilted module. TBPS rings are quantised 0/40/47/55/60/68/72
  // degrees, so ONE layer holds both flat and tilted modules -- which makes it
  // a control with the subdetector, sensors and material held fixed.
  float mod_tilt = -999.f;
  float hit_q_half_len = -999.f;   // the HIT's own q extent. NOTE the factor is
                                   // hl_fac = 3 for PIXELS and sqrt(3) for STRIPS
                                   // (HitStructures.cc) -- dividing by the wrong
                                   // one mis-states sigma_hit by 1.7x.
  // The hit's own PHI sigma, projected from its covariance. It does NOT exist in
  // LayerOfHits: q is per-hit and covariance-derived there, while phi uses a FLAT
  // HIT_PHI_HALF_EXTENT = 0.0123 rad, which is the recorded q/phi asymmetry. It is
  // needed because core/sigma_trk is biased high -- the measured residual carries
  // the hit error too -- and without it the phi ratio cannot be corrected the way
  // the q one can.
  float hit_phi_sigma = -999.f;
  float sigma_q_trk = -999.f;      // TrLayerSearch::dq_track / 3
  float sigma_phi_trk = -999.f;    // TrLayerSearch::dphi_track / 3
  float cov_xx = -999.f, cov_xy = -999.f, cov_yy = -999.f, cov_zz = -999.f;

  bool  passed_preselect = false;
  bool  passed_pqueue = false;
  int   rank = -1;

  // Kalman, for those that reached it. chi2 is the 2-D module-plane value:
  // expectation 2, median of chi2_2 = 1.386.
  bool  had_kalman = false;
  bool  accepted = false;
  float chi2 = -999.f;
  float residual_x = -999.f;   // across strip / precise / phi direction
  float residual_y = -999.f;   // along strip / coarse
  float residual_z = -999.f;   // off the module plane; must be ~0

  float t_hermite = -999.f;
  float d_plane_h3 = -999.f;

  // ---- Is the Hermite step as good as an exact helix-plane solve?
  //
  // On the search path the propagation takes sPerp (plane::path_from_perp), so it
  // never solves for the plane: the caller hands it propPar from the Hermite's
  // own crossing and sPerp = dalpha/(inv_pt*inv_k). Therefore residual_z ~ 0 is
  // near-circular and proves nothing about the step. These do.
  //
  // The reference is the exact uniform-B helix from the candidate's own state,
  // solved against the module plane in DOUBLE -- and uniform B is the right
  // reference because the mini-propagator uses Config::Bfield (default arg at
  // MkFinderV2p2.cc:413), not the parametric field.
  float d_hermite_exact = -999.f;  // |Hermite point - exact crossing|, cm
  float dalpha_hermite = -999.f;   // turn angle the Hermite actually took, rad
  float dalpha_exact = -999.f;     // turn angle of the exact crossing
  float rel_ds = -999.f;           // (dalpha_h - dalpha_e)/|dalpha_e|: the RELATIVE
                                   // arc-length error, i.e. the error in the sPerp
                                   // the covariance is transported by
  bool  wrong_crossing = false;    // the nearest root is not the smallest-|alpha| one
  int   n_roots = 0;               // plane crossings found in |alpha| < 1.5
  float alpha_small = -999.f;      // the smallest-|alpha| root

  // The layer's own bounding-surface crossings, solved EXACTLY in double from
  // the same state -- i.e. what sp1 / sp2 should have been. If the plane
  // crossing lies outside [alpha_in, alpha_out] the Hermite was extrapolating
  // beyond the two points it interpolates, which is the same thing t_hermite
  // outside [0,1] reports.
  float alpha_in = -999.f, alpha_out = -999.f;

  // provenance, so a bad case can be replayed
  int   step = -1;         // TrCandState::step -- depth of the candidate in THIS
                           // search, i.e. how many layers in. The failure mode
                           // being chased is sequential (a state wrecked at one
                           // layer loses everything after), so step, not layer,
                           // is the axis it lives on.
  int   search_id = -1;    // TrLayerSearch this hit was scanned for -- the unit
                           // the best-hit decision is made over, so it is what
                           // a pickup-efficiency question has to group by
  int   seed = -1;         // index in trSeeds_
  int   global_seed = -1;  // index in seedTracks_
  int   sim = -1;          // index in simTracks_ (TrCandMeta -- observed unset)

  // Truth join done here instead: seed -> its best-matching sim track, then how
  // many hits that sim track actually HAS in this layer. Without it, "no MC hit
  // was scanned" is ambiguous -- the layer plans are deliberately INCLUSIVE (the
  // union over tracks), so a track legitimately skips layers it never crosses,
  // and that is indistinguishable from a window landing in the wrong place.
  int   sim_label = -1;
  int   n_sim_hits_in_layer = -1;
  // Seed PURITY: of the seed's valid hits, the fraction from its best-matching
  // sim track (Event::simInfoForTrack, hit-based). Added 2026-09-21 to separate
  // "the candidate is following the wrong thing" from "the track really kinked":
  // chi2 failures were measured to CLUSTER within a candidate (all-killed
  // candidates 300-25000x more common than binomial), which rules out
  // independent hard scatters and points at a wrong state.
  float seed_good_frac = -1.f;
  int   seed_n_valid = -1, seed_n_match = -1;
  float pt = -999.f, eta = -999.f;
};


// ---------------------------------------------------------------------------
// ValSearchMiss -- ONE ROW PER LAYER-SEARCH, asking the opposite question to
// ValSearchHit: not "of the hits we scanned, which were true", but "where was
// the sim track's own hit, and why did we not even consider it".
//
// Every measurement conditioned on a hit being scanned is blind to the case
// that matters here -- the window opened somewhere the true hit is not, or the
// layer holds no true hit at all because the search is walking a plan layer the
// track never crossed.
//
// Window conventions, from TrLayerSearch: the scan covers
// phi in [phi_center +- phi_delta] and q in [q_min, q_max], where q is z in the
// barrel and r in the endcap. dphi_track / dq_track are 3 sigma of the TRACK
// alone; the q cut additionally allows the hit's own extent.
struct ValSearchMiss {
  int   event = -1;
  int   search_id = -1;
  int   step = -1;
  int   layer = -1;
  bool  is_barrel = false;
  int   sim_label = -1;
  float pt = -999.f, eta = -999.f;

  // The WSR verdict the search itself acted on (TrLayerSearch::wsr): 0 inside,
  // 1 edge, 2 outside. With Config::V2p2::Policy::use_wsr on -- the default -- a
  // WSR_Outside candidate scans NO hits at all, so its row is a record of a
  // search that was declined, not of one that failed. FILTER ON THIS.
  signed char wsr = -1;
  bool  wsr_in_gap = false;

  int   n_sim_in_layer = 0;     // countSimHitsInLayer for the SEARCHED layer
  int   n_scanned = 0;          // hits the search actually visited
  bool  mc_scanned = false;     // any of them was the sim track's
  bool  mc_preselect = false;
  bool  mc_kalman = false;

  // The sim track's own hit in the SEARCHED layer, nearest in phi if several.
  // Residuals are to the window CENTRE, and normalised to the window half-width
  // so that |.| > 1 means "outside the window that was opened".
  bool  has_sim_here = false;
  float sim_dphi = -999.f, sim_dq = -999.f;
  float sim_dphi_norm = -999.f;  // |dphi| / phi_delta
  float sim_dq_norm = -999.f;    // signed position in [q_min, q_max]: 0 = centre,
                                 // +-1 = the edges, beyond = outside
  float sim_d3d = -999.f;        // 3-D distance to the window centre, cm

  // If the searched layer holds no sim hit: where IS the nearest one?
  int   near_layer = -1;
  float near_d3d = -999.f;

  // Why it was not considered. Mutually exclusive, most upstream first:
  //   0 no sim hit in this layer   1 outside the opened window
  //   2 in the window but not scanned (bin range / hit mask)
  //   3 scanned, failed pre-selection
  //   4 pre-selected, evicted before Kalman
  //   5 reached the Kalman, failed the acceptance cut
  //   6 passed the cut, but on no candidate that survived the layer
  //   7 on a candidate that survived the layer
  // The cut is TrKalmanUpdate::chi2_acc < V2p2Params::hit_chi2_cut (chi2 < 30 before
  // 2026-09-26); "survived" is TrKalmanUpdate::accepted. Before 2026-09-30,
  // 6 and 7 were "outranked by a lower chi2" and "lowest chi2", with chi2 < 30.
  int   verdict = -1;
  float mc_chi2 = -999.f;    // best chi2 among the sim track's hits here
  float mc_chi2_acc = -999.f;  // its chi2_acc
  float best_chi2 = -999.f;  // best chi2 among hits that passed the cut
  // Hits that passed the acceptance cut in this search, those of them that are
  // NOT the sim track's, and those of them kept on a surviving candidate. In a
  // verdict-0 search (no sim hit in the layer) every passing hit is a wrong one.
  int   n_pass = 0;
  int   n_pass_wrong = 0;
  int   n_kept_wrong = 0;
  float best_wrong_chi2 = -999.f;      // lowest chi2 among passing hits that are not the sim track's
  float best_wrong_chi2_acc = -999.f;  // lowest chi2_acc among them
};


// ---------------------------------------------------------------------------
// ValCovStep -- ONE KALMAN UPDATE in the search, with the covariance BEFORE and
// AFTER it. This is the instrument for "where does the q covariance contract".
//
// Per layer the covariance should do two things, in opposite directions:
//   propagation + material  GROW it
//   the Kalman update       SHRINK it
// Slide 350 shows the net result is a 11x contraction across |eta| while the
// residual it must contain grows. Only a per-step split can say whether that is
// the update shrinking too hard or the propagation not growing enough -- the
// aggregate cannot distinguish them.
//
// sigma is projected onto the search's own coordinates, not read off a diagonal:
//   barrel  q = z    -> sigma_q^2 = err(2,2)
//   endcap  q = r    -> sigma_q^2 = (x^2 exx + 2xy exy + y^2 eyy) / r^2
//   phi (position)   -> sigma_phi^2 = (y^2 exx - 2xy exy + x^2 eyy) / r^4
// The endcap q form is the true marginal; note LayerOfHits uses the TRACE
// (exx+eyy) for the HIT extent, which is a different, looser quantity (see
// doc/MkFinderV2p2-DesignNotes.md, "Hit extents"). Do not compare the two
// without saying which is which.
struct ValCovStep {
  int   event = -1;
  int   layer = -1;
  int   step = -1;          // depth of the candidate in this search
  bool  is_barrel = false;
  bool  mc_match = false;
  bool  accepted = false;   // this hit won its layer
  float pt = -999.f, eta = -999.f;
  float chi2 = -999.f;

  // TrCandState ids, so the chain can be walked offline: the updated state at
  // one layer is what gets propagated to the next.
  int   state_in = -1, state_out = -1;

  float sig_q_prop = -1.f, sig_q_upd = -1.f;      // cm, before / after the update
  float sig_phi_prop = -1.f, sig_phi_upd = -1.f;  // rad
};

#endif
