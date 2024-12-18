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

    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Neg];
      sp.reserve_plan(2 + 12);  // BPix + FPix-; BPix3 & 4 are out of acceptance
      sp.fill_plan(0, 1);
      sp.fill_plan(38, 49);     // FPix- all 12
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Neg];
      sp.reserve_plan(4 + 8 + 6 + 10);  // BPix + FPix- + TOB- +TEC-
      sp.fill_plan( 0,  3);
      sp.fill_plan(38, 45);  // FPix-, first 8 layers
      sp.fill_plan( 4, 15);  // TOB, 6 double layers
      sp.fill_plan(50, 59);  // TEC, 5 double disks
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Barrel];
      sp.reserve_plan(4 + 6 + 6);  // BPix + TOB-1 + TOB-2
      sp.fill_plan( 0,  3);          //                  [ 0,  3]
      sp.fill_plan( 4,  9);          // TOB-1, 6 layers  [ 4,  9]
      sp.fill_plan(10, 15);        // TOB-2, 8 layers  [10, 17]
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Pos];
      sp.reserve_plan(4 + 3 + 6 + 6 + 8 + 18);  // BPix + FPix+ + TIB + TID+ + TOB + TEC+
      sp.fill_plan( 0,  3);
      sp.fill_plan(16, 23);  // FPix-, first 8 layers
      sp.fill_plan( 4, 15);  // TOB, 6 double layers
      sp.fill_plan(28, 37);  // TEC, 5 double disks
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
    spv[TrackerInfo::Reg_Endcap_Neg].set_iterator_limits(2, 0, 3);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(2, 0, 4);
    spv[TrackerInfo::Reg_Barrel].set_iterator_limits(2, 0, 2);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(2, 0, 4);
    spv[TrackerInfo::Reg_Endcap_Pos].set_iterator_limits(2, 0, 3);
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
      ip.useHitSelectionV2 = true;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    }
  }



  void Create_CMS_phase2(TrackerInfo &ti, IterationsInfo &ii, bool verbose) {
    // TrackerInfo needs to be loaded from a bin-file.
    if (ti.n_layers() != 60) {
      fprintf(stderr, "Create_CMS_phase2() FATAL TrackerInfo shold have been initialized from a binary file\n"
                       "with the same name as the geometry library and a '.bin' suffix.\n");
      throw std::runtime_error("Create_CMS_phase2 TrackerIngo not initialized");
    }
    // ti.print_tracker(2); // 1 - print layers, 2 - print layers and modules

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

    ii.resize(enable_all_iters_for_seed_cleaning_tests ? 10 : 1);

    ii[0].set_iteration_index_and_track_algorithm(0, (int)TrackBase::TrackAlgorithm::initialStep);
    ii[0].set_num_regions_layers(5, 60); // 16 + 22 + 22

    // Fills TrackerInfo/LayerInfo and default windows of ii[0].m_layer_configs
    setup_default_windows(ti, ii[0]);

    ii[0].m_seed_cleaner_name = "phase1:default";
    ii[0].m_default_track_scorer_name = "phase1:default";

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
