//-------------------
// CMS phase2 geometry
//-------------------

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/interface/HitStructures.h"
#include "RecoTracker/MkFitCore/interface/TrackStructures.h"

// missing #include "CMS-phase2-HitSelectionWindows.h"

#include <functional>

using namespace mkfit;

namespace {
// missing #include "CMS-phase2.acc"

  void SetupCoreSteeringParams(IterationConfig &ic) {
    ic.m_region_order[0] = TrackerInfo::Reg_Transition_Pos;
    ic.m_region_order[1] = TrackerInfo::Reg_Transition_Neg;
    ic.m_region_order[2] = TrackerInfo::Reg_Endcap_Pos;
    ic.m_region_order[3] = TrackerInfo::Reg_Endcap_Neg;
    ic.m_region_order[4] = TrackerInfo::Reg_Barrel;

    // NOTE: v2p2 can potentially handle double OT layers, other build methods can not.
    // This means runtime switching between build methods when double layers are active
    // does not work.
    // OBSOLETE COMMENT, kept for the record: "The swap is needed as S layers are
    // before P in layer ordering, probably to be changed." That WAS changed --
    // commit 10781bd48eb (2025-09-16) "Swap PS layer indices so P comes before S"
    // added the `+1 - isStereo` in LayerNumberConverter::convertLayerNumber(), so
    // P (isStereo=1, macro-pixel, 1.5 mm) now gets the LOWER mkFit layer number
    // and S (isStereo=0, strip, 24 mm) the higher. OT_swap_pairs = false is
    // therefore correct now; it no longer needs to compensate for anything.
    //
    // *** DOUBLE LAYERS ARE CURRENTLY OFF. *** With as_single_entry = false,
    // fill_plan_pairs(4,9,...) emits LayerControl(4), LayerControl(5), ... which is
    // identical to fill_plan(4,9). So LayerControl::m_layer_sec is never set,
    // has_second_layer() is always false, and every double-layer path in
    // MkFinderV2p2 (m_rz_limits.m_is_double, BL_s, the m_layer_sec block near
    // MkFinderV2p2.cc:923) is dead code today.
    //
    // Three things to fix BEFORE switching this on -- see RecoTracker/CLAUDE.md,
    // "Double layers: switched off at the config level":
    //
    //  1. TOB 2S (10-15) is never paired, though the geometry says it should be:
    //     is_stereo alternates 0,1,0,1,0,1 across 10-15 exactly as across 4-9, and
    //     the comments below call them "outer 3 double layers" -- but they go
    //     through plain fill_plan(10, 15) while 4-9 go through fill_plan_pairs().
    //     Enabling these flags pairs TBPS and TEC and silently leaves TOB 2S as
    //     six singles, in all three regions that reference it.
    //
    //  2. SetupBackwardSearch() below hardcodes plan-INDEX arithmetic that assumes
    //     the unpaired plan. The transition plan is 4+8+6+6+10 = 34 entries
    //     unpaired but 4+8+3+6+5 = 26 paired, so its index 27 goes out of range --
    //     and SteeringParams::iterator::is_valid() only tests != -1, making
    //     m_layer_plan[27] an unchecked OOB read. Fix: capture m_layer_plan.size()
    //     after each fill_plan* call here, or look the pickup up by layer number.
    //
    //  3. Pairing replaces two propagations -- hence two material applications --
    //     per physical layer with one, so it roughly HALVES the applied multiple
    //     scattering in TBPS and TEC. Direction of the error, carefully:
    //     TrackerInfo::material_radl() is a true bin-local weighted average
    //     (MkFitGeometryESProducer::aggregateMaterialInfo(): rl/weight), not a
    //     per-layer effective value calibrated against this plan. But
    //     applyMaterialEffects() samples exactly ONE bin, at the propagation
    //     destination, with no path-length scaling -- while the 1 cm x 1 cm grid
    //     means a 3.4-6.6 cm thick OT shell spans 3-7 bins. So one sample already
    //     under-integrates the traversal, two samples under-integrate less, and
    //     pairing moves FURTHER from the true integral: less scattering, tighter
    //     covariance, the direction that already hurts. Expect an efficiency
    //     regression from this alone, independent of any bug. The real fix is
    //     path-length scaling / integration, listed in RecoTracker/CLAUDE.md
    //     under "Revisit where material is applied".
    //
    //  Note also: mkFit SPLITS one physical CMS OT layer into two mkFit layers,
    //  using CMSSW's TrackerTopology::isStereo(detid) as the criterion. So "double
    //  layer" here means re-uniting what the layer numbering divided, not
    //  modelling a genuine two-sensor stereo structure -- CMS phase-2 OT pairs are
    //  parallel, no crossing angle, unlike ATLAS ITk. ("stereo" is CMSSW's word;
    //  in phase-2 PS modules it happens to tag the P / macro-pixel sensor.)
    //
    // On OT_swap_pairs: nearly meaningless for this geometry. From the dump the
    // pair members are nested, near-coincident shells -- L4 r(22.14, 28.73) and
    // L5 r(22.39, 28.54), offset 2.5 mm out of a 6.5 cm shell, because TBPS is
    // tilted. The sub-layers are radially INTERLEAVED, so "which one is first" is
    // not well defined at layer granularity; the ordering that matters is the
    // path-length ordering of individual hits, merged across both sub-layers.
    const bool OT_as_single_entry = false;
    const bool OT_swap_pairs = false;

    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Neg];
      sp.reserve_plan(2 + 12);  // BPix + FPix-; BPix3 & 4 are out of acceptance

      sp.fill_plan(0, 1);
      sp.fill_plan(38, 49);     // FPix- all 12

      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Neg];
      sp.reserve_plan(4 + 8 + 12 + 10);     // BPix + FPix- + TOB- +TEC-

      sp.fill_plan( 0,  3);
      sp.fill_plan(38, 45);    // FPix-, first 8 layers

      // TOB, inner 3 double layers, PS
      sp.fill_plan_pairs( 4, 9, OT_as_single_entry, OT_swap_pairs);

      sp.fill_plan_pairs(10, 15, OT_as_single_entry, OT_swap_pairs);  // TOB, outer 3 double layers, 2S

      // TEC, 5 double disks, radially half PS, half 2S
      sp.fill_plan_pairs(50, 59, OT_as_single_entry, OT_swap_pairs);

      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Barrel];
      sp.reserve_plan(4 + 6 + 6);  // BPix + TOB-1 + TOB-2

      sp.fill_plan( 0,  3);       //      [ 0,  3]

      // TOB-1, 6 layers  [ 4,  9] PS
      sp.fill_plan_pairs( 4,  9, OT_as_single_entry, OT_swap_pairs);

      sp.fill_plan_pairs(10, 15, OT_as_single_entry, OT_swap_pairs);  // TOB-2, 6 layers [10,15] 2S

      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Pos];
      sp.reserve_plan(4 + 8 + 12 + 10);  // BPix + FPix+ + TOB+ + TEC+

      sp.fill_plan( 0,  3);
      sp.fill_plan(16, 23);   // FPix-, first 8 layers

      // TOB, inner 3 double layers, PS
      sp.fill_plan_pairs( 4, 9, OT_as_single_entry, OT_swap_pairs);

      sp.fill_plan_pairs(10, 15, OT_as_single_entry, OT_swap_pairs);  // TOB, outer 3 double layers, 2S

      // TEC, 5 double disks, radially half PS, half 2S
      sp.fill_plan_pairs(28, 37, OT_as_single_entry, OT_swap_pairs);

      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Pos];
      sp.reserve_plan(2 + 12);  // BPix + FPix+; BPix3 & 4 are out of acceptance

      sp.fill_plan( 0,  1);
      sp.fill_plan(16, 27);     // FPix- all 12

      sp.set_iterator_limits(2, 0);
    }
  }

  void setup_default_windows(TrackerInfo &ti, IterationConfig &ic) {
    // XXXX To be improved. Also, linear coefs for window functions are NOT set.
    // Loop over layers, setup something based on q-bins / pixel vs not pixel and phase1 settings.
    for (int l = 0; l < ti.n_layers(); ++l) {
      LayerInfo &li = ti.layer_nc(l);
      IterationLayerConfig &ilc = ic.layer(l);

      if (li.is_pixel()) {
        if (li.is_barrel())
          ilc.set_selection_limits(0.01, 0.02, 1.0, 2.0);
        else
          ilc.set_selection_limits(0.01, 0.02, 0.8, 1.6);
      } else {
        if (li.is_barrel())
          ilc.set_selection_limits(0.01, 0.02, 3.0, 5.0);
        else
          ilc.set_selection_limits(0.01, 0.02, 3.0, 5.0);
      }
    }
  }

  void SetupBackwardSearch(IterationConfig &ic) {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_drop_seed_hits = false;
    ic.m_backward_fit_min_hits = 99;
    auto &spv = ic.m_steering_params;
    // XXXX Recheck those limits !!!
    // The bkw-search start layer is set for LST T5 seeds, mostly.
    //
    // These USED TO BE plan INDICES, spelled as arithmetic restating the plan's
    // structure ("4 + 8 + 2*6 + 3" = 27). That is not merely fragile, it is a
    // latent out-of-range read: iterator::is_valid() only tests != -1, and the
    // moment OT_as_single_entry flips, the two transition plans go from 34
    // entries to 26 and index 27 walks off the end of m_layer_plan.
    //
    // Named by LAYER instead, which is a property of the detector and does not
    // move when the plan is rebuilt. plan_index_of_layer() matches either member
    // of a paired entry, so these resolve to the SAME plan entry in both configs
    // (verified: transition 27 -> 22 unpaired -> paired, both being layer 31/53).
    // The layers are unchanged from the indices they replace, so this commit does
    // not move the pickup point.
    spv[TrackerInfo::Reg_Endcap_Neg].set_iterator_limits(2, 0);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(2, 0);
    spv[TrackerInfo::Reg_Barrel].set_iterator_limits(2, 0);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(2, 0);
    spv[TrackerInfo::Reg_Endcap_Pos].set_iterator_limits(2, 0);

    spv[TrackerInfo::Reg_Endcap_Neg    ].set_bkw_search_pickup_at_layer(43); // FPix- TFPX6
    spv[TrackerInfo::Reg_Transition_Neg].set_bkw_search_pickup_at_layer(53); // TEC-  TEDD2
    spv[TrackerInfo::Reg_Barrel        ].set_bkw_search_pickup_at_layer(10); // TB2S  OTLayer4
    spv[TrackerInfo::Reg_Transition_Pos].set_bkw_search_pickup_at_layer(31); // TEC+  TEDD2
    spv[TrackerInfo::Reg_Endcap_Pos    ].set_bkw_search_pickup_at_layer(21); // FPix+ TFPX6
  }

  void SetupIterationParams(IterationParams &ip, unsigned int it = 0) {
    if (it == 0) {
      ip.nlayers_per_seed = 4;
      ip.maxCandsPerSeed = 6;
      ip.maxHolesPerCand = 4;
      ip.maxConsecHoles = 2;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.useHitSelectionV2 = true; // relevant for Std and CE, not for V2
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    }
  }



  void Create_CMS_phase2(TrackerInfo &ti, IterationsInfo &ii, bool verbose) {
    // TrackerInfo needs to be loaded from a bin-file.
    if (ti.n_layers() != 60) {
      fprintf(stderr, "Create_CMS_phase2() FATAL TrackerInfo should have been initialized from a binary file\n"
                       "with the same name as the geometry library and a '.bin' suffix.\n");
      throw std::runtime_error("Create_CMS_phase2 TrackerIngo not initialized");
    }
    ti.print_tracker(1); // 1 - print layers, 2 - print layers and modules

    // In cmssw, this is set in the GeometryESProducer for Phase2.
    // We also have --use-p2p 0|1 in mkFit.exe
    Config::usePropToPlane = true;

    // Likewise set in the GeometryESProducer for Phase2 in cmssw
    // (MkFitGeometryESProducer.cc:563) -- but Config.cc:8 defaults it to FALSE,
    // so standalone was running without it and phase-2 set neither. With it off,
    // applyMaterialEffects() adds multiple scattering only to err(4,4) and
    // err(5,5) and puts NONE into 1/pT -- precisely the element a momentum
    // resolution or chi2 study depends on. Every standalone covariance number
    // taken before 2026-09-10 has that baked in. Also reachable as --use-ptms 0|1.
    Config::usePtMultScat = true;

    PropagationConfig &pconf = ti.prop_config_nc();
    pconf.backward_fit_to_pca = Config::includePCA;
    pconf.finding_requires_propagation_to_hit_pos = true;
    pconf.finding_inter_layer_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    if (Config::usePropToPlane)
      pconf.finding_intra_layer_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    else
      pconf.finding_intra_layer_pflags = PropagationFlags(PF_none);
    pconf.backward_fit_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    pconf.forward_fit_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    pconf.seed_fit_pflags = PropagationFlags(PF_none);
    pconf.pca_prop_pflags = PropagationFlags(PF_none);
    pconf.apply_tracker_info(&ti);

    const bool enable_all_iters_for_seed_cleaning_tests = false;
    // NOTE: if setting the above to true, also set
    //    Config::nItersCMSSW = 10; // or whatever number
    // or use --num-iters-cmssw num

    ii.resize(enable_all_iters_for_seed_cleaning_tests ? 10 : 1);

    ii[0].set_iteration_index_and_track_algorithm(0, (int)TrackBase::TrackAlgorithm::initialStep);
    ii[0].set_num_regions_layers(5, 60); // 16 + 22 + 22

    // Fills TrackerInfo/LayerInfo and default windows of ii[0].m_layer_configs
    setup_default_windows(ti, ii[0]);

    ii[0].m_seed_cleaner_name = "phase1:default";
    // ii[0].m_default_track_scorer_name = "phase1:default";
    ii[0].m_default_track_scorer_name = "phase2:LstIntoPix";

    ii[0].m_seed_partitioner_name = "phase2:1";

    SetupCoreSteeringParams(ii[0]);

    SetupIterationParams(ii[0].m_params, 0);
    ii[0].set_dupl_params(0.24, 0.002, 0.004, 0.008);
    ii[0].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits_pixelseed";

    SetupBackwardSearch(ii[0]);

    // Clone Phase1 seed cleaning
    if (enable_all_iters_for_seed_cleaning_tests) {
      ii[1].set_iteration_index_and_track_algorithm(1, (int)TrackBase::TrackAlgorithm::highPtTripletStep);
      ii[1].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.018, 0.018, 0.036, 0.10, 0.036, 0.10);
      ii[1].m_seed_cleaner_name = "phase1:default";

      ii[2].set_iteration_index_and_track_algorithm(2, (int)TrackBase::TrackAlgorithm::lowPtQuadStep);
      ii[2].set_seed_cleaning_params(0.5, 0.05, 0.05, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10);
      ii[2].m_seed_cleaner_name = "phase1:default";

      ii[3].set_iteration_index_and_track_algorithm(3, (int)TrackBase::TrackAlgorithm::lowPtTripletStep);
      ii[3].set_seed_cleaning_params(0.5, 0.05, 0.05, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10);
      ii[3].m_seed_cleaner_name = "phase1:default";

      ii[4].set_iteration_index_and_track_algorithm(4, (int)TrackBase::TrackAlgorithm::detachedQuadStep);
      ii[4].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10);
      ii[4].m_seed_cleaner_name = "phase1:default";

      ii[5].set_iteration_index_and_track_algorithm(5, (int)TrackBase::TrackAlgorithm::detachedTripletStep);
      ii[5].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10);
      ii[5].m_seed_cleaner_name = "phase1:default";

      ii[6].set_iteration_index_and_track_algorithm(6, (int)TrackBase::TrackAlgorithm::mixedTripletStep);
      ii[6].set_seed_cleaning_params(2.0, 0.05, 0.05, 0.135, 0.135, 0.05, 0.05, 0.135, 0.135);
      ii[6].m_seed_cleaner_name = "phase1:default";

      ii[7].set_iteration_index_and_track_algorithm(7, (int)TrackBase::TrackAlgorithm::pixelLessStep);
      ii[7].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);
      ii[7].m_seed_cleaner_name = ""; // No seed cleaning.

      ii[8].set_iteration_index_and_track_algorithm(8, (int)TrackBase::TrackAlgorithm::tobTecStep);
      ii[8].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);
      ii[8].m_seed_cleaner_name = ""; // No seed cleaning.

      ii[9].set_iteration_index_and_track_algorithm(9, (int)TrackBase::TrackAlgorithm::pixelPairStep);
      ii[9].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);
      ii[9].m_seed_cleaner_name = "phase1:default";
    }

    if (verbose) {
      printf("==========================================================================================\n");
    }

    printf("CMS-phase2 -- Create_TrackerInfo finished\n");

    if (verbose) {
      printf("==========================================================================================\n");
      for (int ii = 0; ii < ti.n_layers(); ++ii)
        ti.layer(ii).print_layer();
      printf("==========================================================================================\n");
    }
  }
}  // namespace

void *TrackerInfoCreator_ptr = (void *)Create_CMS_phase2;
