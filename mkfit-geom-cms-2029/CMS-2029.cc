//-------------------
// CMS 2029 geometry
//-------------------

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/standalone/ConfigStandalone.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/interface/HitStructures.h"
#include "RecoTracker/MkFitCore/interface/TrackStructures.h"

// missing #include "CMS-2029-HitSelectionWindows.h"

#include <functional>

using namespace mkfit;

namespace {
// missing #include "CMS-2029.acc"

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
    // Loop over layers, setup something based on q-bins / pixel vs not pixel and 2017 settings.
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
      ip.minPtCut = 0.0;
      ip.maxClusterSize = 8;
    }
  }


  //=================
  // partitionSeeds1
  //=================

  [[maybe_unused]] void partitionSeeds1(const TrackerInfo &trk_info,
                                       const TrackVec &in_seeds,
                                       const EventOfHits &eoh,
                                       IterationSeedPartition &part) {
    // Seeds are placed into eta regions and sorted on region + eta.

    // Merge mono and stereo limits for relevant layers / parameters.
    // TrackerInfo could hold joint limits for sub-detectors.
    const auto &L = trk_info;
    const float tecp1_rin = std::min(L[28].rin(), L[29].rin());
    const float tecp1_rout = std::max(L[28].rout(), L[29].rout());
    const float tecp1_zmin = std::min(L[28].zmin(), L[29].zmin());

    const float tecp2_rin = std::min(L[30].rin(), L[31].rin());
    const float tecp2_zmax = std::max(L[30].zmax(), L[31].zmax());

    const float tecn1_rin = std::min(L[50].rin(), L[51].rin());
    const float tecn1_rout = std::max(L[50].rout(), L[51].rout());
    const float tecn1_zmax = std::max(L[50].zmax(), L[51].zmax());

    const float tecn2_rin = std::min(L[52].rin(), L[53].rin());
    const float tecn2_zmin = std::min(L[52].zmin(), L[53].zmin());

    const float tec_z_extra = 0.0f;  // 10.0f;

    const int size = in_seeds.size();

    auto barrel_pos_check = [](const Track &S, float maxR, float rin, float zmax) -> bool {
      bool inside = maxR > rin && S.zAtR(rin) < zmax;
      return inside;
    };

    auto barrel_neg_check = [](const Track &S, float maxR, float rin, float zmin) -> bool {
      bool inside = maxR > rin && S.zAtR(rin) > zmin;
      return inside;
    };

    auto endcap_pos_check = [](const Track &S, float maxR, float rout, float rin, float zmin) -> bool {
      bool inside = maxR > rout ? S.zAtR(rout) > zmin : (maxR > rin && S.zAtR(maxR) > zmin);
      return inside;
    };

    auto endcap_neg_check = [](const Track &S, float maxR, float rout, float rin, float zmax) -> bool {
      bool inside = maxR > rout ? S.zAtR(rout) < zmax : (maxR > rin && S.zAtR(maxR) < zmax);
      return inside;
    };

    for (int i = 0; i < size; ++i) {
      const Track &S = in_seeds[i];

      HitOnTrack hot = S.getLastHitOnTrack();
      float eta = eoh[hot.layer].refHit(hot.index).eta();

      // Region to be defined by propagation / intersection tests
      TrackerInfo::EtaRegion reg;

      const bool z_dir_pos = S.pz() > 0;
      const float maxR = S.maxReachRadius();

      if (z_dir_pos) {
        bool in_tec_as_brl = barrel_pos_check(S, maxR, tecp2_rin, tecp2_zmax);

        if (!in_tec_as_brl) {
          reg = TrackerInfo::Reg_Endcap_Pos;
        } else {
          bool in_tec = endcap_pos_check(S, maxR, tecp1_rout, tecp1_rin, tecp1_zmin - tec_z_extra);

          if (!in_tec) {
            reg = TrackerInfo::Reg_Barrel;
          } else {
            reg = TrackerInfo::Reg_Transition_Pos;
          }
        }
      } else {
        bool in_tec_as_brl = barrel_neg_check(S, maxR, tecn2_rin, tecn2_zmin);

        if (!in_tec_as_brl) {
          reg = TrackerInfo::Reg_Endcap_Neg;
        } else {
          bool in_tec = endcap_neg_check(S, maxR, tecn1_rout, tecn1_rin, tecn1_zmax + tec_z_extra);

          if (!in_tec) {
            reg = TrackerInfo::Reg_Barrel;
          } else {
            reg = TrackerInfo::Reg_Transition_Neg;
          }
        }
      }

      part.m_region[i] = reg;
      if (part.m_phi_eta_foo)
        part.m_phi_eta_foo(eoh[hot.layer].refHit(hot.index).phi(), eta);
    }
  }

  //======================
  // partitionSeeds1debug
  //======================

  [[maybe_unused]] void partitionSeeds1debug(const TrackerInfo &trk_info,
                                             const TrackVec &in_seeds,
                                             const EventOfHits &eoh,
                                             IterationSeedPartition &part) {
    // Seeds are placed into eta regions and sorted on region + eta.

    // Merge mono and stereo limits for relevant layers / parameters.
    // TrackerInfo could hold joint limits for sub-detectors.
    const auto &L = trk_info;
    const float tecp1_rin = std::min(L[28].rin(), L[29].rin());
    const float tecp1_rout = std::max(L[28].rout(), L[29].rout());
    const float tecp1_zmin = std::min(L[28].zmin(), L[29].zmin());

    const float tecp2_rin = std::min(L[30].rin(), L[31].rin());
    const float tecp2_zmax = std::max(L[30].zmax(), L[31].zmax());

    const float tecn1_rin = std::min(L[50].rin(), L[51].rin());
    const float tecn1_rout = std::max(L[50].rout(), L[51].rout());
    const float tecn1_zmax = std::max(L[50].zmax(), L[51].zmax());

    const float tecn2_rin = std::min(L[52].rin(), L[53].rin());
    const float tecn2_zmin = std::min(L[52].zmin(), L[53].zmin());

    const float tec_z_extra = 0.0f;  // 10.0f;

    const int size = in_seeds.size();

    auto barrel_pos_check = [](const Track &S, float maxR, float rin, float zmax, const char *det) -> bool {
      bool inside = maxR > rin && S.zAtR(rin) < zmax;

      printf("  in_%s=%d  maxR=%7.3f, rin=%7.3f -- ", det, inside, maxR, rin);
      if (maxR > rin) {
        printf("maxR > rin:   S.zAtR(rin) < zmax  -- %.3f <? %.3f\n", S.zAtR(rin), zmax);
      } else {
        printf("maxR < rin: no pie.\n");
      }

      return inside;
    };

    auto barrel_neg_check = [](const Track &S, float maxR, float rin, float zmin, const char *det) -> bool {
      bool inside = maxR > rin && S.zAtR(rin) > zmin;

      printf("  in_%s=%d  maxR=%7.3f, rin=%7.3f -- ", det, inside, maxR, rin);
      if (maxR > rin) {
        printf("maxR > rin:   S.zAtR(rin) > zmin  -- %.3f >? %.3f\n", S.zAtR(rin), zmin);
      } else {
        printf("maxR < rin: no pie.\n");
      }

      return inside;
    };

    auto endcap_pos_check = [](const Track &S, float maxR, float rout, float rin, float zmin, const char *det) -> bool {
      bool inside = maxR > rout ? S.zAtR(rout) > zmin : (maxR > rin && S.zAtR(maxR) > zmin);

      printf("  in_%s=%d  maxR=%7.3f, rout=%7.3f, rin=%7.3f -- ", det, inside, maxR, rout, rin);
      if (maxR > rout) {
        printf("maxR > rout:  S.zAtR(rout) > zmin  -- %.3f >? %.3f\n", S.zAtR(rout), zmin);
      } else if (maxR > rin) {
        printf("maxR > rin:   S.zAtR(maxR) > zmin) -- %.3f >? %.3f\n", S.zAtR(maxR), zmin);
      } else {
        printf("maxR < rin: no pie.\n");
      }

      return inside;
    };

    auto endcap_neg_check = [](const Track &S, float maxR, float rout, float rin, float zmax, const char *det) -> bool {
      bool inside = maxR > rout ? S.zAtR(rout) < zmax : (maxR > rin && S.zAtR(maxR) < zmax);

      printf("  in_%s=%d  maxR=%7.3f, rout=%7.3f, rin=%7.3f -- ", det, inside, maxR, rout, rin);
      if (maxR > rout) {
        printf("maxR > rout:  S.zAtR(rout) < zmax  -- %.3f <? %.3f\n", S.zAtR(rout), zmax);
      } else if (maxR > rin) {
        printf("maxR > rin:   S.zAtR(maxR) < zmax  -- %.3f <? %.3f\n", S.zAtR(maxR), zmax);
      } else {
        printf("maxR < rin: no pie.\n");
      }

      return inside;
    };

    for (int i = 0; i < size; ++i) {
      const Track &S = in_seeds[i];

      HitOnTrack hot = S.getLastHitOnTrack();
      float eta = eoh[hot.layer].refHit(hot.index).eta();
      // float  eta = S.momEta();

      // Region to be defined by propagation / intersection tests
      TrackerInfo::EtaRegion reg;

      const bool z_dir_pos = S.pz() > 0;
      const float maxR = S.maxReachRadius();

      printf("partitionSeeds1debug seed index %d, z_dir_pos=%d (pz=%.3f), maxR=%.3f\n", i, z_dir_pos, S.pz(), maxR);

      if (z_dir_pos) {
        bool in_tec_as_brl = barrel_pos_check(S, maxR, tecp2_rin, tecp2_zmax, "TECasBarrelp");

        if (!in_tec_as_brl) {
          reg = TrackerInfo::Reg_Endcap_Pos;
          printf("  --> region = %d, endcap pos\n", reg);
        } else {
          bool in_tec = endcap_pos_check(S, maxR, tecp1_rout, tecp1_rin, tecp1_zmin - tec_z_extra, "TECp");

          if (!in_tec) {
            reg = TrackerInfo::Reg_Barrel;
            printf("  --> region = %d, barrel\n", reg);
          } else {
            reg = TrackerInfo::Reg_Transition_Pos;
            printf("  --> region = %d, transition pos\n", reg);
          }
        }
      } else {
        bool in_tec_as_brl = barrel_neg_check(S, maxR, tecn2_rin, tecn2_zmin, "TECasBarreln");

        if (!in_tec_as_brl) {
          reg = TrackerInfo::Reg_Endcap_Neg;
          printf("  --> region = %d, endcap neg\n", reg);
        } else {
          bool in_tec = endcap_neg_check(S, maxR, tecn1_rout, tecn1_rin, tecn1_zmax + tec_z_extra, "TECn");

          if (!in_tec) {
            reg = TrackerInfo::Reg_Barrel;
            printf("  --> region = %d, barrel\n", reg);
          } else {
            reg = TrackerInfo::Reg_Transition_Neg;
            printf("  --> region = %d, transition neg\n", reg);
          }
        }
      }

      part.m_region[i] = reg;
      if (part.m_phi_eta_foo)
        part.m_phi_eta_foo(eoh[hot.layer].refHit(hot.index).phi(), eta);
    }
  }

  void Create_CMS_2029(TrackerInfo &ti, IterationsInfo &ii, bool verbose) {
    // TrackerInfo needs to be loaded from a bin-file.
    if (ti.n_layers() != 60) {
      fprintf(stderr, "Create_CMS_2029() FATAL TrackerInfo shold have been initialized from a binary file\n"
                       "with the same name as the geometry library and a '.bin' suffix.\n");
      throw std::runtime_error("Create_CMS_2029 TrackerIngo not initialized");
    }
    // ti.print_tracker(2); // 1 - print layers, 2 - print layers and modules

    PropagationConfig pconf;
    pconf.backward_fit_to_pca = Config::includePCA;
    pconf.finding_requires_propagation_to_hit_pos = true;
    pconf.finding_inter_layer_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    pconf.finding_intra_layer_pflags = PropagationFlags(PF_none);
    pconf.backward_fit_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    pconf.forward_fit_pflags = PropagationFlags(PF_use_param_b_field | PF_apply_material);
    pconf.seed_fit_pflags = PropagationFlags(PF_none);
    pconf.pca_prop_pflags = PropagationFlags(PF_none);
    pconf.set_as_default();

    ii.resize(1);
    ii[0].set_iteration_index_and_track_algorithm(0, (int)TrackBase::TrackAlgorithm::initialStep);
    ii[0].set_num_regions_layers(5, 60); // 16 + 22 + 22

    // Fills TrackerInfo/LayerInfo and default windows of ii[0].m_layer_configs
    setup_default_windows(ti, ii[0]);

    ii[0].m_seed_partitioner = partitionSeeds1;

    SetupCoreSteeringParams(ii[0]);

    SetupIterationParams(ii[0].m_params, 0);
    ii[0].set_dupclean_flag();
    ii[0].set_dupl_params(0.24, 0.002, 0.004, 0.008);

    SetupBackwardSearch(ii[0]);

    if (verbose) {
      printf("==========================================================================================\n");
    }

    printf("CMS-2029 -- Create_TrackerInfo finished\n");

    if (verbose) {
      printf("==========================================================================================\n");
      for (int ii = 0; ii < ti.n_layers(); ++ii)
        ti.layer(ii).print_layer();
      printf("==========================================================================================\n");
    }
  }
}  // namespace

void *TrackerInfoCreator_ptr = (void *)Create_CMS_2029;
