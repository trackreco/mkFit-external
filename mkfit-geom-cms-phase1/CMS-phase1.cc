//-------------------
// CMS phase1 geometry
//-------------------

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/interface/HitStructures.h"
#include "RecoTracker/MkFitCore/interface/TrackStructures.h"

#include <functional>

using namespace mkfit;

namespace {
#include "CMS-phase1.acc"

  void SetupCoreSteeringParams_PixelQuad(IterationConfig &ic) {
    ic.m_region_order[0] = TrackerInfo::Reg_Transition_Pos;
    ic.m_region_order[1] = TrackerInfo::Reg_Transition_Neg;
    ic.m_region_order[2] = TrackerInfo::Reg_Endcap_Pos;
    ic.m_region_order[3] = TrackerInfo::Reg_Endcap_Neg;
    ic.m_region_order[4] = TrackerInfo::Reg_Barrel;

    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Neg];
      sp.reserve_plan(3 + 3 + 6 + 18);  // BPix + FPix- + TID- + TEC-; BPix4 is out of acceptance
      sp.append_plan(1);
      sp.append_plan(0);
      sp.append_plan(2);
      // sp.fill_plan(0, 2);
      sp.fill_plan(45, 47);
      sp.fill_plan(48, 53);  // TID,  6 disks (3 mono + 3 stereo)
      sp.fill_plan(54, 71);  // TEC, 18 disks (9 mono + 9 stereo)
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Neg];
      sp.reserve_plan(4 + 3 + 6 + 6 + 8 + 18);  // BPix + FPix- + TIB + TID- + TOB + TEC-
      sp.append_plan(1);
      sp.append_plan(0);
      sp.fill_plan(2, 3);
      //sp.fill_plan(0, 3);
      sp.fill_plan(45, 47);
      sp.fill_plan(4, 9);    // TIB,  6 layers (4 mono + 2 stereo)
      sp.fill_plan(48, 53);  // TID,  6 disks  (3 mono + 3 stereo)
      sp.fill_plan(10, 17);  // TOB,  8 layers (6 mono + 2 stereo)
      sp.fill_plan(54, 71);  // TEC, 18 disks  (9 mono + 9 stereo)
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Barrel];
      sp.reserve_plan(4 + 6 + 8);  // BPix + TIB + TOB
      sp.fill_plan(0, 3);          //                                    [ 0,  3]
      sp.fill_plan(4, 9);          // TIB, 6 layers (4 mono + 2 stereo)  [ 4,  9]
      sp.fill_plan(10, 17);        // TOB, 8 layers (6 mono + 2 stereo)  [10, 17]
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Pos];
      sp.reserve_plan(4 + 3 + 6 + 6 + 8 + 18);  // BPix + FPix+ + TIB + TID+ + TOB + TEC+
      sp.append_plan(1);
      sp.append_plan(0);
      sp.fill_plan(2, 3);
      // sp.fill_plan(0, 3);
      sp.fill_plan(18, 20);  //                                     [ 4,  6]
      sp.fill_plan(4, 9);    // TIB,  6 layers (4 mono + 2 stereo)  [ 7, 12]
      sp.fill_plan(21, 26);  // TID,  6 disks  (3 mono + 3 stereo)  [13, 18]
      sp.fill_plan(10, 17);  // TOB,  8 layers (6 mono + 2 stereo)  [19, 26]
      sp.fill_plan(27, 44);  // TEC, 18 disks  (9 mono + 9 stereo)  [27, 44]
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Pos];
      sp.reserve_plan(3 + 3 + 6 + 18);  // BPix + FPix+ + TID+ + TEC+; BPix4 is out of acceptance
      sp.append_plan(1);
      sp.append_plan(0);
      sp.append_plan(2);
      // sp.fill_plan(0, 2);
      sp.fill_plan(18, 20);  //                                     [ 3,  5]
      sp.fill_plan(21, 26);  // TID,  6 disks  (3 mono + 3 stereo)  [ 6, 11]
      sp.fill_plan(27, 44);  // TEC, 18 disks  (9 mono + 9 stereo)  [12, 29]
      sp.set_iterator_limits(2, 0);
    }
  }

  void SetupCoreSteeringParams_Common(IterationConfig &ic) {
    ic.m_region_order[0] = TrackerInfo::Reg_Transition_Pos;
    ic.m_region_order[1] = TrackerInfo::Reg_Transition_Neg;
    ic.m_region_order[2] = TrackerInfo::Reg_Endcap_Pos;
    ic.m_region_order[3] = TrackerInfo::Reg_Endcap_Neg;
    ic.m_region_order[4] = TrackerInfo::Reg_Barrel;

    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Neg];
      sp.reserve_plan(3 + 3 + 6 + 18);  // BPix + FPix- + TID- + TEC-; BPix4 is out of acceptance
      sp.fill_plan(0, 2);
      sp.fill_plan(45, 47);
      sp.fill_plan(48, 53);  // TID,  6 disks (3 mono + 3 stereo)
      sp.fill_plan(54, 71);  // TEC, 18 disks (9 mono + 9 stereo)
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Neg];
      sp.reserve_plan(4 + 3 + 6 + 6 + 8 + 18);  // BPix + FPix- + TIB + TID- + TOB + TEC-
      sp.fill_plan(0, 3);
      sp.fill_plan(45, 47);
      sp.fill_plan(4, 9);    // TIB,  6 layers (4 mono + 2 stereo)
      sp.fill_plan(48, 53);  // TID,  6 disks  (3 mono + 3 stereo)
      sp.fill_plan(10, 17);  // TOB,  8 layers (6 mono + 2 stereo)
      sp.fill_plan(54, 71);  // TEC, 18 disks  (9 mono + 9 stereo)
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Barrel];
      sp.reserve_plan(4 + 6 + 8);  // BPix + TIB + TOB
      sp.fill_plan(0, 3);          //                                    [ 0,  3]
      sp.fill_plan(4, 9);          // TIB, 6 layers (4 mono + 2 stereo)  [ 4,  9]
      sp.fill_plan(10, 17);        // TOB, 8 layers (6 mono + 2 stereo)  [10, 17]
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Transition_Pos];
      sp.reserve_plan(4 + 3 + 6 + 6 + 8 + 18);  // BPix + FPix+ + TIB + TID+ + TOB + TEC+
      sp.fill_plan(0, 3);                       //                                     [ 0,  3]
      sp.fill_plan(18, 20);                     //                                     [ 4,  6]
      sp.fill_plan(4, 9);                       // TIB,  6 layers (4 mono + 2 stereo)  [ 7, 12]
      sp.fill_plan(21, 26);                     // TID,  6 disks  (3 mono + 3 stereo)  [13, 18]
      sp.fill_plan(10, 17);                     // TOB,  8 layers (6 mono + 2 stereo)  [19, 26]
      sp.fill_plan(27, 44);                     // TEC, 18 disks  (9 mono + 9 stereo)  [27, 44]
      sp.set_iterator_limits(2, 0);
    }
    {
      SteeringParams &sp = ic.m_steering_params[TrackerInfo::Reg_Endcap_Pos];
      sp.reserve_plan(3 + 3 + 6 + 18);  // BPix + FPix+ + TID+ + TEC+; BPix4 is out of acceptance
      sp.fill_plan(0, 2);               //                                     [ 0,  2]
      sp.fill_plan(18, 20);             //                                     [ 3,  5]
      sp.fill_plan(21, 26);             // TID,  6 disks  (3 mono + 3 stereo)  [ 6, 11]
      sp.fill_plan(27, 44);             // TEC, 18 disks  (9 mono + 9 stereo)  [12, 29]
      sp.set_iterator_limits(2, 0);
    }
  }

  /*
  ////// Example backward search setup for initialStep iteration (currently 'replaced' by seed duplicate merging)
  void SetupBackwardSearch_Iter0(IterationConfig& ic)
  {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_drop_seed_hits = true;
    ic.m_backward_fit_min_hits   = 8;
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg]    .set_iterator_limits(2, 3, 5);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(2, 3, 7);
    spv[TrackerInfo::Reg_Barrel]        .set_iterator_limits(2, 3, 4);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(2, 3, 7);
    spv[TrackerInfo::Reg_Endcap_Pos]    .set_iterator_limits(2, 3, 5);
  }
*/

  void SetupBackwardSearch_PixelCommon(IterationConfig &ic) {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_drop_seed_hits = false;
    ic.m_backward_fit_min_hits = 99;
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg].set_iterator_limits(2, 0, 3);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(2, 0, 4);
    spv[TrackerInfo::Reg_Barrel].set_iterator_limits(2, 0, 2);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(2, 0, 4);
    spv[TrackerInfo::Reg_Endcap_Pos].set_iterator_limits(2, 0, 3);
  }

  void SetupBackwardSearch_Iter7(IterationConfig &ic) {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_params.maxHolesPerCand = 2;
    ic.m_backward_params.maxConsecHoles = 2;
    // Remove pixel layers from FwdSearch, add them to BkwSearch
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg].set_iterator_limits(8, 6, 19);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(9, 7, 34);
    spv[TrackerInfo::Reg_Barrel].set_iterator_limits(6, 4, 8);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(9, 7, 34);
    spv[TrackerInfo::Reg_Endcap_Pos].set_iterator_limits(8, 6, 19);
  }

  void SetupBackwardSearch_Iter8(IterationConfig &ic) {
    ic.m_backward_params = ic.m_params;
    ic.m_backward_search = true;
    ic.m_backward_params.maxHolesPerCand = 2;
    ic.m_backward_params.maxConsecHoles = 2;
    // Remove pixel/tib/tid layers from FwdSearch, add them to BkwSearch/
    auto &spv = ic.m_steering_params;
    spv[TrackerInfo::Reg_Endcap_Neg].set_iterator_limits(12, 12, 24);
    spv[TrackerInfo::Reg_Transition_Neg].set_iterator_limits(22, 19, 39);
    spv[TrackerInfo::Reg_Barrel].set_iterator_limits(12, 10, 14);
    spv[TrackerInfo::Reg_Transition_Pos].set_iterator_limits(22, 19, 39);
    spv[TrackerInfo::Reg_Endcap_Pos].set_iterator_limits(12, 12, 24);
  }

  void SetupIterationParams(IterationParams &ip, unsigned int it = 0) {
    if (it == 0) {
      ip.nlayers_per_seed = 4;
      ip.maxCandsPerSeed = 5;
      ip.maxHolesPerCand = 4;
      ip.maxConsecHoles = 2;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    } else if (it == 1)  // for triplet steps, nlayers_per_seed=3
    {
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed = 5;
      ip.maxHolesPerCand = 4;
      ip.maxConsecHoles = 2;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    } else if (it == 2) {
      ip.nlayers_per_seed = 4;
      ip.maxCandsPerSeed = 5;
      ip.maxHolesPerCand = 4;
      ip.maxConsecHoles = 2;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    } else if (it == 3)  // for triplet steps, nlayers_per_seed=3
    {
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed = 5;
      ip.maxHolesPerCand = 4;
      ip.maxConsecHoles = 2;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    } else if (it == 4) {
      ip.nlayers_per_seed = 4;
      ip.maxCandsPerSeed = 5;
      ip.maxHolesPerCand = 4;
      ip.maxConsecHoles = 2;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    } else if (it == 5)  // for triplet steps, nlayers_per_seed=3
    {
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed = 5;
      ip.maxHolesPerCand = 4;
      ip.maxConsecHoles = 2;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    } else if (it == 6)  // for triplet steps, nlayers_per_seed=3; for mixeTripletSetp, also maxCandsPerSeed=2
    {
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed = 2;
      ip.maxHolesPerCand = 4;
      ip.maxConsecHoles = 2;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    } else if (it == 7)  // for PixelLess step, maxCandsPerSeed=2 and maxHolesPerCand=maxConsecHoles=0
    {
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed = 2;
      ip.maxHolesPerCand = 0;
      ip.maxConsecHoles = 1;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    } else if (it == 8)  // for TobTec step, maxCandsPerSeed=2 and maxHolesPerCand=maxConsecHoles=0
    {
      ip.nlayers_per_seed = 3;
      ip.maxCandsPerSeed = 2;
      ip.maxHolesPerCand = 0;
      ip.maxConsecHoles = 1;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    } else if (it == 9)  // addign also pixel pair step - algo -> 6
    {
      ip.nlayers_per_seed = 2;
      ip.maxCandsPerSeed = 3;
      ip.maxHolesPerCand = 4;
      ip.maxConsecHoles = 2;
      ip.chi2Cut_min = 15.0;
      ip.chi2CutOverlap = 3.5;
      ip.pTCutOverlap = 0.0;
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    }
  }

  void fill_hit_selection_windows_params(IterationsInfo &ii) {
    // Read meta and data vectors.
#include "CMS-phase1-HitSelectionWindows.h"

    // Build algo to iteration index map.
    std::map<int, int> algo2idx;
    for (int i = 0; i < (int) ii.size(); ++i)
      algo2idx.insert({ii[i].m_track_algorithm, i});

    int data_pos = 0;

    for (int l = 0; l < (int) meta.size(); l += 3) {
      int algo = meta[l];
      auto li = algo2idx.find(algo);
      if (li == algo2idx.end())
        throw std::runtime_error("unknown algorithm");
      auto &iinfo = ii[li->second];

      int layer = meta[l+1];
      if (layer < 0 || layer >= (int) iinfo.m_layer_configs.size())
        throw std::runtime_error("layer out of range");
      auto  &linfo = iinfo.m_layer_configs[layer];

      int fwd = meta[l+2];
      auto &v = fwd ? linfo.m_winpars_fwd : linfo.m_winpars_bkw;
      
      v.reserve(12);
      for (int i = 0; i < 12; ++i)
        v.push_back(data[data_pos++]);
    }
  }

  void Create_CMS_phase1(TrackerInfo &ti, IterationsInfo &ii, bool verbose) {
    // TrackerInfo needs to be loaded from a bin-file.
    if (ti.n_layers() != 72) {
      fprintf(stderr, "Create_CMS_phase1() FATAL TrackerInfo shold have been initialized from a binary file\n"
                       "with the same name as the geometry library and a '.bin' suffix.\n");
      throw std::runtime_error("Create_CMS_phase1 TrackerIngo not initialized");
    }
    // ti.print_tracker(2); // 1 - print layers, 2 - print layers and modules

    PropagationConfig pconf;
    pconf.backward_fit_to_pca = Config::includePCA;
    pconf.finding_requires_propagation_to_hit_pos = true;
    pconf.finding_inter_layer_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    pconf.finding_intra_layer_pflags = PropagationFlags(PF_none);
    // PF_copy_input_state_on_fail is only set for barrel in MkFinder::bkFitFitTracks() itself.
    pconf.backward_fit_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material); // | PF_copy_input_state_on_fail);
    pconf.forward_fit_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    pconf.seed_fit_pflags = PropagationFlags(PF_none);
    pconf.pca_prop_pflags = PropagationFlags(PF_none);
    pconf.set_as_default();

    const int N_iter = 10;

    ii.resize(N_iter);
    ii[0].set_iteration_index_and_track_algorithm(0, (int)TrackBase::TrackAlgorithm::initialStep);
    ii[0].set_num_regions_layers(5, 72);

    // Fills TrackerInfo/LayerInfo and default windows of ii[0].m_layer_configs
    Create_CMS_phase1_AutoGen(ti, ii);
    ii[0].m_seed_cleaner_name = "phase1:default";
    ii[0].m_seed_partitioner_name = "phase1:1";
    ii[0].m_default_track_scorer_name = "phase1:default";
    SetupCoreSteeringParams_PixelQuad(ii[0]);

    // At this point copy out layer/steering stuff for reuse in later iterations.
    IterationConfig def_itconf_pixelquad;
    def_itconf_pixelquad.cloneLayerSteerCore(ii[0]);

    SetupIterationParams(ii[0].m_params, 0);
    ii[0].set_dupl_params(0.24, 0.002, 0.004, 0.008);
    ii[0].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits_pixelseed";
    SetupBackwardSearch_PixelCommon(ii[0]);

    ii[1].set_num_regions_layers(5, 72);
    ii[1].m_layer_configs = ii[0].m_layer_configs;
    ii[1].m_seed_cleaner_name = "phase1:default";
    ii[1].m_seed_partitioner_name = "phase1:1";
    ii[1].m_default_track_scorer_name = "phase1:default";

    SetupCoreSteeringParams_Common(ii[1]);

    // At this point copy out layer/steering stuff for reuse in later iterations.
    IterationConfig def_itconf_common;
    def_itconf_common.cloneLayerSteerCore(ii[1]);

    SetupIterationParams(ii[1].m_params, 1);
    ii[1].set_iteration_index_and_track_algorithm(1, (int)TrackBase::TrackAlgorithm::highPtTripletStep);
    ii[1].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.018, 0.018, 0.036, 0.10, 0.036, 0.10);
    ii[1].set_dupl_params(0.24, 0.03, 0.05, 0.08);
    ii[1].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits_pixelseed";
    SetupBackwardSearch_PixelCommon(ii[1]);

    ii[2].cloneLayerSteerCore(def_itconf_pixelquad);
    SetupIterationParams(ii[2].m_params, 2);
    ii[2].set_iteration_index_and_track_algorithm(2, (int)TrackBase::TrackAlgorithm::lowPtQuadStep);
    ii[2].set_seed_cleaning_params(0.5, 0.05, 0.05, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10);
    ii[2].set_dupl_params(0.5, 0.01, 0.03, 0.05);
    ii[2].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits_pixelseed";
    SetupBackwardSearch_PixelCommon(ii[2]);

    ii[3].cloneLayerSteerCore(def_itconf_common);
    SetupIterationParams(ii[3].m_params, 3);
    ii[3].set_iteration_index_and_track_algorithm(3, (int)TrackBase::TrackAlgorithm::lowPtTripletStep);
    ii[3].set_seed_cleaning_params(0.5, 0.05, 0.05, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10);
    ii[3].set_dupl_params(0.33, 0.018, 0.05, 0.018);
    ii[3].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits_pixelseed";
    SetupBackwardSearch_PixelCommon(ii[3]);

    ii[4].cloneLayerSteerCore(def_itconf_pixelquad);
    SetupIterationParams(ii[4].m_params, 4);
    ii[4].set_iteration_index_and_track_algorithm(4, (int)TrackBase::TrackAlgorithm::detachedQuadStep);
    ii[4].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10);
    ii[4].set_dupl_params(0.24, 0.018, 0.05, 0.05);
    ii[4].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits_pixelseed";
    SetupBackwardSearch_PixelCommon(ii[4]);

    ii[5].cloneLayerSteerCore(def_itconf_common);
    SetupIterationParams(ii[5].m_params, 5);
    ii[5].set_iteration_index_and_track_algorithm(5, (int)TrackBase::TrackAlgorithm::detachedTripletStep);
    ii[5].set_seed_cleaning_params(2.0, 0.018, 0.018, 0.05, 0.05, 0.10, 0.10, 0.10, 0.10);
    ii[5].m_pre_bkfit_filter_name = ""; // no pre-filter
    ii[5].m_post_bkfit_filter_name = "phase1:qfilter_n_layers";
    ii[5].set_dupl_params(0.24, 0.01, 0.01, 0.1);
    ii[5].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits_pixelseed";
    SetupBackwardSearch_PixelCommon(ii[5]);

    ii[6].cloneLayerSteerCore(def_itconf_common);
    SetupIterationParams(ii[6].m_params, 6);
    ii[6].set_iteration_index_and_track_algorithm(6, (int)TrackBase::TrackAlgorithm::mixedTripletStep);
    ii[6].set_seed_cleaning_params(2.0, 0.05, 0.05, 0.135, 0.135, 0.05, 0.05, 0.135, 0.135);
    ii[6].set_dupl_params(0.2, 0.05, 0.05, 0.05);
    ii[6].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits_pixelseed";
    SetupBackwardSearch_PixelCommon(ii[6]);

    ii[7].cloneLayerSteerCore(def_itconf_common);
    SetupIterationParams(ii[7].m_params, 7);
    ii[7].set_iteration_index_and_track_algorithm(7, (int)TrackBase::TrackAlgorithm::pixelLessStep);
    ii[7].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);
    ii[7].m_seed_cleaner_name = ""; // No seed cleaning.
    ii[7].m_requires_seed_hit_sorting = true;
    ii[7].m_pre_bkfit_filter_name = "phase1:qfilter_pixelLessFwd";
    ii[7].m_post_bkfit_filter_name = "phase1:qfilter_pixelLessBkwd";
    ii[7].dc_fracSharedHits = 0.14;
    ii[7].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits";
    SetupBackwardSearch_Iter7(ii[7]);

    ii[8].cloneLayerSteerCore(def_itconf_common);
    SetupIterationParams(ii[8].m_params, 8);
    ii[8].set_iteration_index_and_track_algorithm(8, (int)TrackBase::TrackAlgorithm::tobTecStep);
    ii[8].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);
    ii[8].m_seed_cleaner_name = ""; // No seed cleaning.
    ii[8].m_requires_seed_hit_sorting = true;
    ii[8].m_params.minHitsQF = 4;
    ii[8].m_pre_bkfit_filter_name = "phase1:qfilter_n_hits";
    ii[8].m_post_bkfit_filter_name = ""; // no post-filter
    ii[8].dc_fracSharedHits = 0.25;
    ii[8].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits";
    SetupBackwardSearch_Iter8(ii[8]);

    ii[9].cloneLayerSteerCore(def_itconf_common);
    SetupIterationParams(ii[9].m_params, 9);
    ii[9].set_iteration_index_and_track_algorithm(9, (int)TrackBase::TrackAlgorithm::pixelPairStep);
    ii[9].set_seed_cleaning_params(2.0, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135);
    ii[9].m_params.minHitsQF = 3;
    ii[9].m_pre_bkfit_filter_name = "phase1:qfilter_n_hits_pixseed";
    ii[9].m_post_bkfit_filter_name = ""; // no post-filter
    ii[9].set_dupl_params(0.5, 0.03, 0.05, 0.05);
    ii[9].m_duplicate_cleaner_name = "phase1:clean_duplicates_sharedhits_pixelseed";
    SetupBackwardSearch_PixelCommon(ii[9]);

    fill_hit_selection_windows_params(ii);

    if (verbose) {
      printf("==========================================================================================\n");
    }

    printf("CMS-phase1 -- Create_TrackerInfo finished\n");

    if (verbose) {
      printf("==========================================================================================\n");
      for (int li = 0; li < ti.n_layers(); ++li)
        ti.layer(li).print_layer();
      printf("==========================================================================================\n");
    }
  }
}  // namespace

void *TrackerInfoCreator_ptr = (void *)Create_CMS_phase1;
