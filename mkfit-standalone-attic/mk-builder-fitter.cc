// Code from MkBuilder and MkFitter that was used in original build-test-mplex
// steering functions, test best-hit, standard, and clone-engine -- before
// the advent of multiple iterations and MIMI steering.

namespace mkfit {

    //----------------------------
    // Carryover from MkBuilder.h
    //----------------------------

    void create_seeds_from_sim_tracks();
    void find_seeds();
    void fit_seeds();
    void PrepareSeeds();

    // This one is here mostly to keep code for printing overlap hit truth - DUMP_OVERLAP_RTTS.
    quality_store_tracks(const EventOfCombCandidates &eoccs, TrackVec & tracks);


    //---------------------------
    // Carryover from MkFitter.h
    //---------------------------

    void FitTracks(const int N_proc, const Event * ev, const PropagationFlags pflags);
    void FitTracksSteered(const bool is_barrel[], const int N_proc, const Event * ev, const PropagationFlags pflags);
    void CollectFitValidation(const int hi, const int N_proc, const Event * ev) const;

}


//=============================================================================

//==============================
// Carryover from MkBuilder.cc
//==============================

//------------------------------------------------------------------------------
// Seeding functions: importing, finding and fitting
//------------------------------------------------------------------------------

void MkBuilder::create_seeds_from_sim_tracks()
{
  // Import from simtrack snatching first Config::nlayers_per_seed hits.
  //
  // Reduce number of hits, pick endcap over barrel when both are available.
  //   This is done by assumin endcap hit is after the barrel, no checks.

  // bool debug = true;

  TrackerInfo &trk_info = Config::TrkInfo;
  TrackVec    &sims     = m_event->simTracks_;
  TrackVec    &seeds    = m_event->seedTracks_;

  const int size = sims.size();
  seeds.clear();       // Needed when reading from file and then recreating from sim.
  seeds.reserve(size);

  dprintf("MkBuilder::create_seeds_from_sim_tracks processing %d simtracks.\n", size);

  for (int i = 0; i < size; ++i)
  {
    const Track &src = sims[i];

    dprintf("  [%d] pT=%f eta=%f n_hits=%d lbl=%d\n", i, src.pT(), src.momEta(), src.nFoundHits(), src.label());

    if (src.isNotFindable())
    {
      dprintf("  [%d] not findable.\n", i);
      continue;
    }

    int h_sel = 0, h = 0;
    const HitOnTrack *hots = src.getHitsOnTrackArray();
    HitOnTrack  new_hots[ Config::nlayers_per_seed_max ];

    // Exit condition -- need to check one more hit after Config::nlayers_per_seed
    // good hits are found.
    bool last_hit_check = false;

    while ( ! last_hit_check && h < src.nTotalHits())
    {
      assert (hots[h].index >= 0 && "Expecting input sim tracks (or seeds later) to not have holes");

      if (h_sel == Config::nlayers_per_seed) last_hit_check = true;

      // Check if hit is on a sibling layer given the previous one.
      if (h_sel > 0 && trk_info.are_layers_siblings(new_hots[h_sel - 1].layer, hots[h].layer))
      {
        dprintf("    [%d] Sibling layers %d %d ... overwriting with new one\n", i,
                new_hots[h_sel - 1].layer, hots[h].layer);

        new_hots[h_sel - 1] = hots[h];
      }
      // Drop further hits on the same layer. It seems hard to select the best one (in any way).
      else if (h_sel > 0 && new_hots[h_sel - 1].layer == hots[h].layer)
      {
        dprintf("    [%d] Hits on the same layer %d ... keeping the first one\n", i, hots[h].layer);
      }
      else if ( ! last_hit_check)
      {
        new_hots[h_sel++] = hots[h];
      }

      ++h;
    }

    if (h_sel < Config::nlayers_per_seed)
    {
      printf("MkBuilder::create_seeds_from_sim_tracks simtrack %d only yielded %d hits. Skipping ...\n",
             src.label(), h_sel);
      continue;
    }

    seeds.emplace_back( Track(src.state(), 0, src.label(), h_sel, new_hots) );

    dprintf("  [%d->%d] Seed nh=%d, last_lay=%d, last_idx=%d\n", i, (int) seeds.size() - 1,
            seeds.back().nTotalHits(), seeds.back().getLastHitLyr(), seeds.back().getLastHitIdx());
    // dprintf("  "); for (int i=0; i<dst.nTotalHits();++i) printf(" (%d/%d)", dst.getHitIdx(i), dst.getHitLyr(i)); printf("\n");
  }

  dprintf("MkBuilder::create_seeds_from_sim_tracks finished processing of %d sim tracks - created %d seeds.\n",
          size, (int) seeds.size());
}

void MkBuilder::find_seeds()
{
  fprintf(stderr, "__FILE__::__LINE__ Needs fixing for B/E support, search for XXMT4K\n");
  exit(1);

#ifdef DEBUG
  bool debug(false);
#endif
  typedef tbb::concurrent_vector<TripletIdx> TripletIdxConVec;
  TripletIdxConVec seed_idcs;

  //double time = dtime();
  findSeedsByRoadSearch(seed_idcs,m_event_of_hits.m_layers_of_hits,m_event->layerHits_[1].size(),m_event);
  //time = dtime() - time;

  // use this to initialize tracks
  // XXMT4K  ... configurable input layers ... or hardcode something else for endcap.
  // Could come from TrackerInfo ...
  // But what about transition ... TrackerInfo as well or arbitrary combination of B/E seed layers ????

  const LayerOfHits &loh0 = m_event_of_hits.m_layers_of_hits[0];
  const LayerOfHits &loh1 = m_event_of_hits.m_layers_of_hits[1];
  const LayerOfHits &loh2 = m_event_of_hits.m_layers_of_hits[2];

  // make seed tracks
  TrackVec & seedtracks = m_event->seedTracks_;
  seedtracks.resize(seed_idcs.size());
  for (size_t iseed = 0; iseed < seedtracks.size(); iseed++)
  {
    auto & seedtrack = seedtracks[iseed];
    seedtrack.setLabel(iseed);

    // use to set charge
    const Hit & hit0 = loh0.refHit(seed_idcs[iseed][0]);
    const Hit & hit1 = loh1.refHit(seed_idcs[iseed][1]);
    const Hit & hit2 = loh2.refHit(seed_idcs[iseed][2]);

    seedtrack.setCharge(calculateCharge(hit0,hit1,hit2));

    for (int ihit = 0; ihit < Config::nlayers_per_seed; ihit++)
    {
      // XXMT4K  - ihit to layer[ihit]
      seedtrack.addHitIdx(seed_idcs[iseed][ihit], ihit, 0.0f);
    }

    for (int ihit = Config::nlayers_per_seed; ihit < Config::nLayers; ihit++)
    {
      seedtrack.setHitIdxLyr(ihit, -1, -1);
    }

    dprint("iseed: " << iseed << " mcids: " << hit0.mcTrackID(m_event->simHitsInfo_) << " " <<
	   hit1.mcTrackID(m_event->simHitsInfo_) << " " << hit1.mcTrackID(m_event->simHitsInfo_));
  }
}

} // end namespace mkfit

namespace
{
  void fill_seed_layer_sig(const Track& trk, int n_hits, bool is_brl[])
  {
    const TrackerInfo &trk_info = Config::TrkInfo;

    for (int i = 0; i < n_hits; ++i)
    {
      is_brl[i] = trk_info.layer( trk.getHitLyr(i) ).is_barrel();
    }
  }

  bool are_seed_layer_sigs_equal(const Track& trk, int n_hits, const bool is_brl_ref[])
  {
    const TrackerInfo &trk_info = Config::TrkInfo;

    for (int i = 0; i < n_hits; ++i)
    {
      if(trk_info.layer( trk.getHitLyr(i) ).is_barrel() != is_brl_ref[i]) return false;
    }

    return true;
  }
}

namespace mkfit {

void MkBuilder::fit_seeds()
{
  // Expect seeds to be sorted in eta (in some way) and that Event::seedEtaSeparators_[]
  // array holds starting indices of 5 eta regions.
  // Within each region it vectorizes the fit as long as layer indices of all seeds match.
  // See layer_sig_change label below.
  // Alternatively, we could be using the layer plan (but it might require addition of
  // a new flag in LayerControl (well, should really change those bools to a bitfield).

  // debug = true;

  TrackVec& seedtracks = m_event->seedTracks_;

  dcall(print_seeds(seedtracks));

  TBB_PARALLEL_FOR_EACH(m_regions.begin(), m_regions.end(),
  [&](int reg)
  {
    RegionOfSeedIndices rosi(m_seedEtaSeparators, reg);

    TBB_PARALLEL_FOR(rosi.tbb_blk_rng_vec(),
      [&](const tbb::blocked_range<int>& blk_rng)
    {
      // printf("TBB seeding krappe -- range = %d to %d - extent = %d ==> %d to %d - extent %d\n",
      //        i.begin(), i.end(), i.end() - i.begin(), beg, std::min(end,theEnd), std::min(end,theEnd) - beg);

      // printf("Seed info pos(  x       y       z        r;     eta    phi)   mom(  pt      pz;     eta    phi)\n");

      FITTER( mkfttr );

      RangeOfSeedIndices rng = rosi.seed_rng(blk_rng);

      while (rng.valid())
      {
#ifdef DEBUG
        // MT dump seed so i see if etas are about right
        for (int i = rng.m_beg; i < rng.m_end; ++i)
        {
          auto &t = seedtracks[i];
          auto &dst = t;
          dprintf("Seed %4d lbl=%d pos(%+7.3f %+7.3f %+7.3f; %+7.3f %+6.3f %+6.3f) mom(%+7.3f %+7.3f; %+6.3f %+6.3f)\n",
                  i, t.label(), t.x(), t.y(), t.z(), t.posR(), t.posEta(), t.posPhi(),
                  t.pT(), t.pz(), t.momEta(), t.momPhi());
          dprintf("  Idx/lay for above track:"); for (int i=0; i<dst.nTotalHits();++i) dprintf(" (%d/%d)", dst.getHitIdx(i), dst.getHitLyr(i)); dprintf("\n");
        }
#endif

        // We had seeds sorted in eta_mom ... but they have dZ displacement ...
        // so they can go through "semi random" barrel/disk pattern close to
        // transition region for overall layers 2 and 3 where eta of barrel is
        // larger than transition region.
        // E.g., for 10k tracks in endcap/barrel the break happens ~250 times,
        // often several times witin the same NN range (5 time is not rare with NN=8).
        //
        // Sorting on eta_pos of the last seed hit yields ~50 breaks on the same set.
        // This is being used now (in import_seeds()).
        //
        // In the following we make sure seed range passed to vectorized
        // function has compatible layer signatures (barrel / endcap).

      layer_sig_change:

        bool is_brl[Config::nlayers_per_seed_max];

        fill_seed_layer_sig(seedtracks[rng.m_beg], Config::nlayers_per_seed, is_brl);

        for (int i = rng.m_beg + 1; i < rng.m_end; ++i)
        {
          if ( ! are_seed_layer_sigs_equal(seedtracks[i], Config::nlayers_per_seed, is_brl))
          {
            dprintf("Breaking seed range due to different layer signature at %d (%d, %d)\n", i, rng.m_beg, rng.m_end);

            fit_one_seed_set(seedtracks, rng.m_beg, i, mkfttr.get(), is_brl);

            rng.m_beg = i;
            goto layer_sig_change;
          }
        }

        fit_one_seed_set(seedtracks, rng.m_beg, rng.m_end, mkfttr.get(), is_brl);

        ++rng;
      }
    });
  });
}

inline void MkBuilder::fit_one_seed_set(TrackVec& seedtracks, int itrack, int end,
                                        MkFitter *mkfttr, const bool is_brl[])
{
  // debug=true;

  mkfttr->setNhits(Config::nlayers_per_seed);
  mkfttr->inputTracksAndHits(seedtracks, m_event_of_hits.m_layers_of_hits, itrack, end);

  if (Config::cf_seeding) mkfttr->ConformalFitTracks(false, itrack, end);

  mkfttr->FitTracksSteered(is_brl, end - itrack, m_event, PropagationConfig::get_default().seed_fit_pflags);

  mkfttr->outputFittedTracksAndHitIdx(m_event->seedTracks_, itrack, end, false);
}


//------------------------------------------------------------------------------
// PrepareSeeds
//------------------------------------------------------------------------------

void MkBuilder::PrepareSeeds()
{
  // {
  //   TrackVec  &tv = m_event->seedTracks_;
  //   char pref[80];
  //   for (int i = 0; i < (int) tv.size(); ++i)
  //   {
  //     sprintf(pref, "Pre-cleaning seed silly value check event=%d index=%d:", m_event->evtID(), i);
  //     tv[i].hasSillyValues(true, false, pref);
  //   }
  // }

  if (Config::seedInput == simSeeds)
  {
    if (Config::clean_cms_simtracks_for_seeding)
    {
      m_event->clean_cms_simtracks();

      // printf("\n* Simtracks after cleaning:\n");
      // m_event->print_tracks(m_event->simTracks_, true);
      // printf("\n");
    }
    // create_seeds_from_sim_tracks();

    seed_post_cleaning(m_event->seedTracks_);
  }
  else if (Config::seedInput == cmsswSeeds)
  {
    m_event->relabel_bad_seedtracks();

    // want to make sure we mark which sim tracks are findable based on cmssw seeds BEFORE seed cleaning
    if (Config::sim_val || Config::quality_val)
    {
      prep_simtracks();
    }

    // need to make a map of seed ids to cmssw tk ids BEFORE seeds are sorted
    if (Config::cmssw_val)
    {
      m_event->validation_.makeSeedTkToCMSSWTkMap(*m_event);
    }

    // this is a hack that allows us to associate seeds with cmssw tracks for the text dump plots
    if (Config::dumpForPlots && Config::readCmsswTracks)
    {
      for (size_t itrack = 0; itrack < m_event->cmsswTracks_.size(); itrack++)
      {
        const auto &cmsswtrack = m_event->cmsswTracks_[itrack];
        const auto cmsswlabel = cmsswtrack.label();
        auto &seedtrack = m_event->seedTracks_[cmsswlabel];
        seedtrack.setLabel(cmsswlabel);
      }
    }

    if (Config::seedCleaning == cleanSeedsN2)
    {
      m_event->clean_cms_seedtracks();

      // Select specific cmssw seed for detailed debug.
      // {
      //   Track xx = m_event->seedTracks_[6];
      //   m_event->seedTracks_.clear();
      //   m_event->seedTracks_.push_back(xx);
      // }
    }
    else if (Config::seedCleaning == cleanSeedsPure)
    {
      m_event->use_seeds_from_cmsswtracks();
    }
    else if (Config::seedCleaning == cleanSeedsBadLabel)
    {
      m_event->clean_cms_seedtracks_badlabel();
    }
    else if (Config::seedCleaning != noCleaning)
    {
      std::cerr << "Specified reading cmssw seeds, but an incorrect seed cleaning option! Exiting..." << std::endl;
      exit(1);
    }

    seed_post_cleaning(m_event->seedTracks_);

    // in rare corner cases, seed tracks could be fully cleaned out: skip mapping if so
    if (m_event->seedTracks_.empty()) return;
  }
  else if (Config::seedInput == findSeeds)
  {
    // MIMI - doesnotwork
    // find_seeds();
  }
  else
  {
    std::cerr << "No input seed collection option selected!! Exiting..." << std::endl;
    exit(1);
  }

  // Do not refit cmssw seeds (this if was nested in fit_one_seed_set() until now).
  // Eventually we can add force-refit option.
  if (Config::seedInput != cmsswSeeds)
  {
    // MIMI - doesnotwork
    // fit_seeds();
  }
}


// Overlap hit truth dumper that was part of MkBuilder::quality_store_tracks()

// #define DUMP_OVERLAP_RTTS

void quality_store_tracks(const EventOfCombCandidates &eoccs, TrackVec &tracks) {

#ifdef DUMP_OVERLAP_RTTS

      // SIMTRACK DUMPERS

      static bool first = true;
      if (first) {
        // ./mkFit ... | perl -ne 'if (/^ZZZ_OVERLAP/) { s/^ZZZ_OVERLAP //og; print; }' > ovlp.rtt
        printf("SSS_OVERLAP label/I:prod_type/I:is_findable/I:layer/I:pt/F:eta/F:phi/F\n");

        printf(
            "SSS_TRACK "
            "label/I:prod_type/I:is_findable/I:pt/F:eta/F:phi/F:nhit_sim/I:nlay_sim/I:novlp/I:novlp_pix/I:novlp_strip/"
            "I:novlp_stereo/I\n");

        first = false;
      }

      for (int i = 0; i < (int)m_event->simTracks_.size(); ++i) {
        Track &bb = m_event->simTracks_[i];

        if (bb.prodType() == Track::ProdType::Signal) {
          bb.sortHitsByLayer();

          int no = 0, npix = 0, nstrip = 0, nstereo = 0, prev_lay = -1, last_ovlp = -1;

          for (int hi = 0; hi < bb.nTotalHits(); ++hi) {
            HitOnTrack hot = bb.getHitOnTrack(hi);

            if (hot.layer == prev_lay && hot.layer != last_ovlp) {
              last_ovlp = hot.layer;

              ++no;

              const LayerInfo &li = Config::TrkInfo.layer(hot.layer);

              if (li.is_pixb_lyr() || li.is_pixe_lyr()) {
                ++npix;
              } else {
                ++nstrip;
              }

              if (li.is_stereo())
                ++nstereo;

              printf("SSS_OVERLAP %d %d %d %d %f %f %f\n",
                     bb.label(),
                     (int)bb.prodType(),
                     bb.isFindable(),
                     hot.layer,
                     bb.pT(),
                     bb.posEta(),
                     bb.posPhi());
            }
            prev_lay = hot.layer;
          }

          printf("SSS_TRACK %d %d %d %f %f %f %d %d %d %d %d %d\n",
                 bb.label(),
                 (int)bb.prodType(),
                 bb.isFindable(),
                 bb.pT(),
                 bb.momEta(),
                 bb.momPhi(),
                 bb.nTotalHits(),
                 bb.nUniqueLayers(),
                 no,
                 npix,
                 nstrip,
                 nstereo);
        }
      }

#endif

      int chi2_500_cnt = 0, chi2_nan_cnt = 0;

      for (int i = 0; i < eoccs.m_size; i++) {
        // See MT-RATS comment below.
        assert(!eoccs.m_candidates[i].empty() && "BackwardFitBH requires output tracks to align with seeds.");

        // take the first one!
        if (!eoccs.m_candidates[i].empty()) {
          const TrackCand &bcand = eoccs.m_candidates[i].front();

          if (std::isnan(bcand.chi2()))
            ++chi2_nan_cnt;
          if (bcand.chi2() > 500)
            ++chi2_500_cnt;

#ifdef DUMP_OVERLAP_RTTS
          // DUMP overlap hits
          int no_good = 0;
          int no_bad = 0;
          int no = 0;  // total, esp for tracks that don't have good label
          const HoTNode *hnp = &bcand.refLastHoTNode();
          while (true) {
            if (hnp->m_index_ovlp >= 0) {
              static bool first = true;
              if (first) {
                // ./mkFit ... | perl -ne 'if (/^ZZZ_OVERLAP/) { s/^ZZZ_OVERLAP //og; print; }' > ovlp.rtt
                printf(
                    "ZZZ_OVERLAP label/I:prod_type/I:is_findable/I:layer/I:pt/F:eta/F:phi/F:"
                    "chi2/F:chi2_ovlp/F:module/I:module_ovlp/I:hit_label/I:hit_label_ovlp/I\n");
                first = false;
              }

              auto &LoH = m_event_of_hits.m_layers_of_hits[hnp->m_hot.layer];

              const Hit &h = LoH.refHit(hnp->m_hot.index);
              const MCHitInfo &mchi = m_event->simHitsInfo_[h.mcHitID()];
              const Hit &o = LoH.refHit(hnp->m_index_ovlp);
              const MCHitInfo &mcoi = m_event->simHitsInfo_[o.mcHitID()];

              const TrackBase &bb =
                  (bcand.label() >= 0) ? (const TrackBase &)m_event->simTracks_[bcand.label()] : bcand;

              if (bcand.label() >= 0) {
                if (bcand.label() == mcoi.mcTrackID())
                  ++no_good;
                else
                  ++no_bad;
              }
              ++no;

              // label/I:can_idx/I:layer/I:pt/F:eta/F:phi/F:chi2/F:chi2_ovlp/F:module/I:module_ovlp/I:hit_label/I:hit_label_ovlp/I
              printf("ZZZ_OVERLAP %d %d %d %d %f %f %f %f %f %u %u %d %d\n",
                     bb.label(),
                     (int)bb.prodType(),
                     bb.isFindable(),
                     hnp->m_hot.layer,
                     bb.pT(),
                     bb.posEta(),
                     bb.posPhi(),
                     hnp->m_chi2,
                     hnp->m_chi2_ovlp,
                     h.detIDinLayer(),
                     o.detIDinLayer(),
                     mchi.mcTrackID(),
                     mcoi.mcTrackID());
            }

            if (hnp->m_prev_idx >= 0)
              hnp = &eoccs.m_candidates[i].m_hots[hnp->m_prev_idx];
            else
              break;
          }

          if (bcand.label() >= 0) {
            static bool first = true;
            if (first) {
              // ./mkFit ... | perl -ne 'if (/^ZZZ_TRACK/) { s/^ZZZ_TRACK //og; print; }' > track.rtt
              printf(
                  "ZZZ_TRACK "
                  "label/I:prod_type/I:is_findable/I:pt/F:eta/F:phi/F:nhit_sim/I:nlay_sim/I:nhit_rec/I:nhit_miss_rec/"
                  "I:novlp/I:novlp_good/I:novlp_bad/I\n");
              first = false;
            }

            const Track &bb = m_event->simTracks_[bcand.label()];

            printf("ZZZ_TRACK %d %d %d %f %f %f %d %d %d %d %d %d %d\n",
                   bb.label(),
                   (int)bb.prodType(),
                   bb.isFindable(),
                   bb.pT(),
                   bb.momEta(),
                   bb.momPhi(),
                   bb.nTotalHits(),
                   bb.nUniqueLayers(),
                   bcand.nFoundHits(),
                   bcand.nMissingHits(),
                   no,
                   no_good,
                   no_bad);
          }
            // DUMP END
#endif

          tracks.emplace_back(bcand.exportTrack());

#ifdef DEBUG_BACKWARD_FIT_BH
          printf("CHITRK %d %g %g %g %g %g\n",
                 bcand.nFoundHits(),
                 bcand.chi2(),
                 bcand.chi2() / (bcand.nFoundHits() * 3 - 6),
                 bcand.pT(),
                 bcand.momPhi(),
                 bcand.theta());
#endif
        }
      }

      if (!Config::silent && (chi2_500_cnt > 0 || chi2_nan_cnt > 0)) {
        std::lock_guard<std::mutex> printlock(Event::printmutex);
        printf("MkBuilder::quality_store_tracks bad track chi2 (backward fit?). is-nan=%d, gt-500=%d.\n",
               chi2_nan_cnt,
               chi2_500_cnt);
      }
    }
}


//=============================================================================

//=============================
// Carryover from MkFitter.cc
//=============================

void MkFitter::FitTracks(const int N_proc, const Event * ev, const PropagationFlags pflags)
{
  // Fitting loop.

  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Note, charge is not passed (line propagation).
    // propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
    //                       Err[iP], Par[iP]);

    propagateTracksToHitR(msPar[hi], N_proc, pflags);

    kalmanUpdate(Err[iP], Par[iP], msErr[hi], msPar[hi],
                 Err[iC], Par[iC], N_proc);

    if (Config::fit_val) MkFitter::CollectFitValidation(hi,N_proc,ev);
  }
  // XXXXX What's with chi2?
}

void MkFitter::CollectFitValidation(const int hi, const int N_proc, const Event * ev) const
{
  for (int n = 0; n < N_proc; ++n)
  {
    const float upt = 1.f/Par[iC](n,3,0);
    const FitVal tmpfitval
    {
      Par[iP].constAt(n,2,0),
      std::sqrt(Err[iP].constAt(n,2,2)),
      getPhi(Par[iP].constAt(n,0,0),Par[iP].constAt(n,1,0)),
      std::sqrt(getPhiErr2(Par[iP](n,0,0),Par[iP](n,1,0),Err[iP](n,0,0),Err[iP](n,1,1),Err[iP](n,0,1))),
      upt,
      std::sqrt(Err[iC](n,3,3))*upt*upt,
      Par[iC](n,4,0),
      std::sqrt(Err[iC](n,4,4)),
      getEta(Par[iC](n,5,0)),
      std::sqrt(Err[iC](n,5,5)/std::sin(Par[iC](n,5,0)))
    };

    ev->validation_.collectFitInfo(tmpfitval,Label(n,0,0),hi);
  }
}

void MkFitter::FitTracksSteered(const bool is_barrel[], const int N_proc, const Event * ev, const PropagationFlags pflags)
{
  // Fitting loop.

  dprintf("MkFitter::FitTracksSteered %d %d %d\n", is_barrel[0], is_barrel[1], is_barrel[2]);

  for (int hi = 0; hi < Nhits; ++hi)
  {
    // Note, charge is not passed (line propagation).
    // propagateLineToRMPlex(Err[iC], Par[iC], msErr[hi], msPar[hi],
    //                       Err[iP], Par[iP]);

    if (is_barrel[hi])
    {
      propagateTracksToHitR(msPar[hi], N_proc, pflags);

      kalmanUpdate(Err[iP], Par[iP], msErr[hi], msPar[hi],
                   Err[iC], Par[iC], N_proc);
    }
    else
    {
      propagateTracksToHitZ(msPar[hi], N_proc, pflags);

      kalmanUpdateEndcap(Err[iP], Par[iP], msErr[hi], msPar[hi],
                         Err[iC], Par[iC], N_proc);
    }

    if (Config::fit_val) MkFitter::CollectFitValidation(hi,N_proc,ev);
  }
  // XXXXX What's with chi2?
}
