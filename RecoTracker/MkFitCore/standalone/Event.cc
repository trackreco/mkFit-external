#include "Event.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"

//#define DEBUG
#include "RecoTracker/MkFitCore/src/Debug.h"

#include <cstring>
#include <memory>
#include <set>

#include <format>

namespace {
  std::unique_ptr<mkfit::Validation> dummyValidation(mkfit::Validation::make_validation("dummy", nullptr));
}

namespace mkfit {

  std::mutex Event::printmutex;

  Event::Event(int evtID, int nLayers) : validation_(*dummyValidation), evtID_(evtID) {
    layerHits_.resize(nLayers);
    layerHitMasks_.resize(nLayers);
  }

  Event::Event(Validation &v, int evtID, int nLayers) : validation_(v), evtID_(evtID) {
    layerHits_.resize(nLayers);
    layerHitMasks_.resize(nLayers);
    validation_.resetValidationMaps();  // need to reset maps for every event.
  }

  void Event::reset(int evtID) {
    evtID_ = evtID;

    for (auto &&l : layerHits_) {
      l.clear();
    }
    for (auto &&l : layerHitMasks_) {
      l.clear();
    }

    simHitsInfo_.clear();
    simTrackStates_.clear();
    simHitStates_.clear();
    simTracks_.clear();
    simTracksExtra_.clear();
    seedTracks_.clear();
    seedTracksExtra_.clear();
    candidateTracks_.clear();
    candidateTracksExtra_.clear();
    fitTracks_.clear();
    fitTracksExtra_.clear();
    cmsswTracks_.clear();
    cmsswTracksExtra_.clear();
    beamSpot_ = {};

    validation_.resetValidationMaps();  // need to reset maps for every event.

  #ifdef MKFIT_TRACE
    trCandMetas_.clear();
    trCandStages_.clear();
    trCandStates_.clear();
    trLayerSearches_.clear();
    trHitMatches_.clear();
    trKalmanUpdates_.clear();
    trBkFitUpdates_.clear();
    trSeeds_.clear();
  #endif
  }

  void Event::validate() {
    // special map needed for sim_val_for_cmssw + set the track scores
    if (Config::sim_val_for_cmssw) {
      validation_.makeRecoTkToSeedTkMapsDumbCMSSW(*this);
      validation_.setTrackScoresDumbCMSSW(*this);
    }

    // standard eff/fr/dr validation
    if (Config::sim_val || Config::sim_val_for_cmssw) {
      validation_.setTrackExtras(*this);
      validation_.makeSimTkToRecoTksMaps(*this);
      validation_.makeSeedTkToRecoTkMaps(*this);
      validation_.fillEfficiencyTree(*this);
      validation_.fillFakeRateTree(*this);
    }

    // special cmssw to mkfit validation
    if (Config::cmssw_val) {
      validation_.makeCMSSWTkToSeedTkMap(*this);
      validation_.makeRecoTkToRecoTkMaps(*this);
      validation_.setTrackExtras(*this);
      validation_.makeCMSSWTkToRecoTksMaps(*this);
      validation_.fillCMSSWEfficiencyTree(*this);
      validation_.fillCMSSWFakeRateTree(*this);
    }

    if (Config::fit_val) {  // fit val for z-phi tuning
      validation_.fillFitTree(*this);
    }
  }

  void Event::printStats(const TrackVec &trks, TrackExtraVec &trkextras) {
    int miss(0), found(0), fp_10(0), fp_20(0), hit8(0), h8_10(0), h8_20(0);

    for (auto &&trk : trks) {
      auto &&extra = trkextras[trk.label()];
      extra.setMCTrackIDInfo(trk, layerHits_, simHitsInfo_, simTracks_, false, true);
      if (extra.mcTrackID() < 0) {
        ++miss;
      } else {
        auto &&mctrk = simTracks_[extra.mcTrackID()];
        auto pr = trk.pT() / mctrk.pT();
        found++;
        bool h8 = trk.nFoundHits() >= 8;
        bool pt10 = pr > 0.9 && pr < 1.1;
        bool pt20 = pr > 0.8 && pr < 1.2;
        fp_10 += pt10;
        fp_20 += pt20;
        hit8 += h8;
        h8_10 += h8 && pt10;
        h8_20 += h8 && pt20;
      }
    }
    std::cout << "found tracks=" << found << "  in pT 10%=" << fp_10 << "  in pT 20%=" << fp_20
              << "     no_mc_assoc=" << miss << std::endl
              << "  nH >= 8   =" << hit8 << "  in pT 10%=" << h8_10 << "  in pT 20%=" << h8_20 << std::endl;
  }

  void Event::write_out(DataFile &data_file) {
    FILE *fp = data_file.f_fp;

    static std::mutex writemutex;
    std::lock_guard<std::mutex> writelock(writemutex);

    auto start = ftell(fp);
    int evsize = sizeof(int);
    fwrite(&evsize, sizeof(int), 1, fp);  // this will be overwritten at the end

    evsize += write_tracks(fp, simTracks_);

    if (data_file.hasSimTrackStates()) {
      int nts = simTrackStates_.size();
      fwrite(&nts, sizeof(int), 1, fp);
      fwrite(&simTrackStates_[0], sizeof(TrackState), nts, fp);
      evsize += sizeof(int) + nts * sizeof(TrackState);
    }

    if (data_file.hasSimHitStates()) {
      int nhs = simHitStates_.size();
      fwrite(&nhs, sizeof(int), 1, fp);
      if (nhs > 0)
        fwrite(&simHitStates_[0], sizeof(SimHitState), nhs, fp);
      evsize += sizeof(int) + nhs * sizeof(SimHitState);
    }

    int nl = layerHits_.size();
    fwrite(&nl, sizeof(int), 1, fp);
    evsize += sizeof(int);
    for (int il = 0; il < nl; ++il) {
      int nh = layerHits_[il].size();
      fwrite(&nh, sizeof(int), 1, fp);
      fwrite(&layerHits_[il][0], sizeof(Hit), nh, fp);
      evsize += sizeof(int) + nh * sizeof(Hit);
    }

    if (data_file.hasHitIterMasks()) {
      //sizes are the same as in layerHits_
      for (int il = 0; il < nl; ++il) {
        int nh = layerHitMasks_[il].size();
        assert(nh == (int)layerHits_[il].size());
        fwrite(&layerHitMasks_[il][0], sizeof(uint64_t), nh, fp);
        evsize += nh * sizeof(uint64_t);
      }
    }

    int nm = simHitsInfo_.size();
    fwrite(&nm, sizeof(int), 1, fp);
    fwrite(&simHitsInfo_[0], sizeof(MCHitInfo), nm, fp);
    evsize += sizeof(int) + nm * sizeof(MCHitInfo);

    if (data_file.hasSeeds()) {
      evsize += write_tracks(fp, seedTracks_);
    }

    if (data_file.hasCmsswTracks()) {
      evsize += write_tracks(fp, cmsswTracks_);
    }

    if (data_file.hasBeamSpot()) {
      fwrite(&beamSpot_, sizeof(BeamSpot), 1, fp);
      evsize += sizeof(BeamSpot);
    }

    fseek(fp, start, SEEK_SET);
    fwrite(&evsize, sizeof(int), 1, fp);
    fseek(fp, 0, SEEK_END);

    //layerHitMap_ is recreated afterwards

    /*
  printf("write %i tracks\n",nt);
  for (int it = 0; it<nt; it++) {
    printf("track with pT=%5.3f\n",simTracks_[it].pT());
    for (int ih=0; ih<simTracks_[it].nTotalHits(); ++ih) {
      printf("hit lyr:%2d idx=%i\n", simTracks_[it].getHitLyr(ih), simTracks_[it].getHitIdx(ih));
    }
  }
  printf("write %i layers\n",nl);
  for (int il = 0; il<nl; il++) {
    printf("write %i hits in layer %i\n",layerHits_[il].size(),il);
    for (int ih = 0; ih<layerHits_[il].size(); ih++) {
      printf("hit with r=%5.3f x=%5.3f y=%5.3f z=%5.3f\n",layerHits_[il][ih].r(),layerHits_[il][ih].x(),layerHits_[il][ih].y(),layerHits_[il][ih].z());
    }
  }
  */
  }

  // #define DUMP_SEEDS
  // #define DUMP_SEED_HITS
  // #define DUMP_TRACKS
  // #define DUMP_TRACK_HITS
  // #define DUMP_LAYER_HITS
  // #define DUMP_REC_TRACKS
  // #define DUMP_REC_TRACK_HITS

  void Event::read_in(DataFile &data_file, FILE *in_fp) {
    FILE *fp = in_fp ? in_fp : data_file.f_fp;

    data_file.advancePosToNextEvent(fp);

    int nt = read_tracks(fp, simTracks_);
    Config::nTracks = nt;

    if (data_file.hasSimTrackStates()) {
      int nts;
      fread(&nts, sizeof(int), 1, fp);
      simTrackStates_.resize(nts);
      fread(&simTrackStates_[0], sizeof(TrackState), nts, fp);
    }

    // Genuinely optional on the read side, unlike simTrackStates_ above: without
    // --read-sim-hit-states the section is seeked past, so a file carrying it
    // costs disk and nothing else.
    if (data_file.hasSimHitStates()) {
      int nhs;
      fread(&nhs, sizeof(int), 1, fp);
      if (Config::readSimHitStates) {
        simHitStates_.resize(nhs);
        if (nhs > 0)
          fread(&simHitStates_[0], sizeof(SimHitState), nhs, fp);
      } else {
        simHitStates_.clear();
        std::fseek(fp, (long) nhs * sizeof(SimHitState), SEEK_CUR);
      }
    }

    int nl;
    fread(&nl, sizeof(int), 1, fp);
    layerHits_.resize(nl);
    layerHitMasks_.resize(nl);
    for (int il = 0; il < nl; ++il) {
      int nh;
      fread(&nh, sizeof(int), 1, fp);
      layerHits_[il].resize(nh);
      layerHitMasks_[il].resize(nh, 0);  //init to 0 by default
      fread(&layerHits_[il][0], sizeof(Hit), nh, fp);
    }

    if (data_file.hasHitIterMasks()) {
      for (int il = 0; il < nl; ++il) {
        int nh = layerHits_[il].size();
        fread(&layerHitMasks_[il][0], sizeof(uint64_t), nh, fp);
      }
    }

    int nm;
    fread(&nm, sizeof(int), 1, fp);
    simHitsInfo_.resize(nm);
    fread(&simHitsInfo_[0], sizeof(MCHitInfo), nm, fp);

    if (data_file.hasSeeds()) {
      int ns = read_tracks(fp, seedTracks_, Config::seedInput != cmsswSeeds);
      (void)ns;

#ifdef DUMP_SEEDS
      printf("Read %i seedtracks (neg value means actual reading was skipped)\n", ns);
      for (int it = 0; it < ns; it++) {
        const Track &ss = seedTracks_[it];
        printf("  %-3i q=%+i pT=%7.3f eta=% 7.3f nHits=%i label=%4i algo=%2i\n",
               it,
               ss.charge(),
               ss.pT(),
               ss.momEta(),
               ss.nFoundHits(),
               ss.label(),
               (int)ss.algorithm());
#ifdef DUMP_SEED_HITS
        for (int ih = 0; ih < seedTracks_[it].nTotalHits(); ++ih) {
          int lyr = seedTracks_[it].getHitLyr(ih);
          int idx = seedTracks_[it].getHitIdx(ih);
          if (idx >= 0) {
            const Hit &hit = layerHits_[lyr][idx];
            printf("    hit %2d lyr=%3d idx=%4d pos r=%7.3f z=% 8.3f   mc_hit=%3d mc_trk=%3d\n",
                   ih,
                   lyr,
                   idx,
                   layerHits_[lyr][idx].r(),
                   layerHits_[lyr][idx].z(),
                   hit.mcHitID(),
                   hit.mcTrackID(simHitsInfo_));
          } else
            printf("    hit %2d idx=%i\n", ih, seedTracks_[it].getHitIdx(ih));
        }
#endif
      }
#endif
    }

    int nert = -99999;
    if (data_file.hasCmsswTracks()) {
      nert = read_tracks(fp, cmsswTracks_, !Config::readCmsswTracks);
      (void)nert;
    }

    /*
    // HACK TO ONLY SELECT ONE PROBLEMATIC TRACK.
    // Note that MC matching gets screwed.
    // Works for MC seeding.
    //
    printf("************** SIM SELECTION HACK IN FORCE ********************\n");
    TrackVec x;
    x.push_back(simTracks_[3]);
    simTracks_.swap(x);
    nt = 1;
  */

#ifdef DUMP_TRACKS
    printf("Read %i simtracks\n", nt);
    for (int it = 0; it < nt; it++) {
      const Track &t = simTracks_[it];
      printf("  %-3i q=%+i pT=%7.3f eta=% 7.3f nHits=%2d  label=%4d\n",
             it,
             t.charge(),
             t.pT(),
             t.momEta(),
             t.nFoundHits(),
             t.label());
#ifdef DUMP_TRACK_HITS
      for (int ih = 0; ih < t.nTotalHits(); ++ih) {
        int lyr = t.getHitLyr(ih);
        int idx = t.getHitIdx(ih);
        if (idx >= 0) {
          const Hit &hit = layerHits_[lyr][idx];
          printf("    hit %2d lyr=%2d idx=%3d pos r=%7.3f x=% 8.3f y=% 8.3f z=% 8.3f   mc_hit=%3d mc_trk=%3d\n",
                 ih,
                 lyr,
                 idx,
                 layerHits_[lyr][idx].r(),
                 layerHits_[lyr][idx].x(),
                 layerHits_[lyr][idx].y(),
                 layerHits_[lyr][idx].z(),
                 hit.mcHitID(),
                 hit.mcTrackID(simHitsInfo_));
        } else
          printf("    hit %2d idx=%i\n", ih, t.getHitIdx(ih));
      }
#endif
    }
#endif
#ifdef DUMP_LAYER_HITS
    // clang-format off
    printf("Read %i layers\n", nl);
    int total_hits = 0;
    float min_delta = 100, max_delta = -100;
    for (int il = 0; il < nl; il++) {
      if (layerHits_[il].empty())
        continue;

      const LayerInfo &linfo = Config::TrkInfo[il];
      printf("Read %i hits in layer %i, is_brl=%d is_pix=%d\n",
             (int)layerHits_[il].size(), il, linfo.is_barrel(), linfo.is_pixel());
      total_hits += layerHits_[il].size();
      for (int ih = 0; ih < (int)layerHits_[il].size(); ih++) {
        const Hit &hit = layerHits_[il][ih];
        unsigned int mid = hit.detIDinLayer();
        const ModuleInfo &mi = linfo.module_info(mid);

        if ( ! linfo.is_pixel() && ! linfo.is_barrel()) {
          float delta = mi.half_length - std::sqrt(3)*std::sqrt(hit.exx() + hit.eyy());
          min_delta = std::min(min_delta, delta);
          max_delta = std::max(max_delta, delta);
        }
        // continue;

        printf("  mcHitID=%-5d r=%7g x=%8g y=%8g z=%8g"
               "  sx=%10.4g sy=%10.4e sz=%10.4e sxy=%10.4e, mhl=%10.4e, delta=%10.4e"
               "  chg_pcm=%u (%u - %u)\n",
               hit.mcHitID(), hit.r(), hit.x(), hit.y(), hit.z(),
               std::sqrt(hit.exx()), std::sqrt(hit.eyy()), std::sqrt(hit.ezz()),
               std::sqrt(hit.exx() + hit.eyy()),
               mi.half_length,
               mi.half_length - std::sqrt(3)*std::sqrt(hit.exx() + hit.eyy()),
               hit.chargePerCM(), hit.minChargePerCM(), hit.maxChargePerCM());
      }
    }
    printf("Total hits in all layers = %d; endcap strips: min_delta=%.5f  max_delta=%.5f\n",
           total_hits, min_delta, max_delta);
    // clang-format on
#endif
#ifdef DUMP_REC_TRACKS
    printf("Read %i rectracks\n", nert);
    for (int it = 0; it < nert; it++) {
      const Track &t = cmsswTracks_[it];
      printf("  %-3i with q=%+i pT=%7.3f eta=% 7.3f nHits=%2d  label=%4d algo=%2d\n",
             it,
             t.charge(),
             t.pT(),
             t.momEta(),
             t.nFoundHits(),
             t.label(),
             (int)t.algorithm());
#ifdef DUMP_REC_TRACK_HITS
      for (int ih = 0; ih < t.nTotalHits(); ++ih) {
        int lyr = t.getHitLyr(ih);
        int idx = t.getHitIdx(ih);
        if (idx >= 0) {
          const Hit &hit = layerHits_[lyr][idx];
          printf("    hit %2d lyr=%2d idx=%3d pos r=%7.3f z=% 8.3f   mc_hit=%3d mc_trk=%3d\n",
                 ih,
                 lyr,
                 idx,
                 hit.r(),
                 hit.z(),
                 hit.mcHitID(),
                 hit.mcTrackID(simHitsInfo_));
        } else
          printf("    hit %2d        idx=%i\n", ih, t.getHitIdx(ih));
      }
#endif
    }
#endif

    if (data_file.hasBeamSpot()) {
      fread(&beamSpot_, sizeof(BeamSpot), 1, fp);
    }

    if (!Config::silent)
      printf("Read complete, %d simtracks on file.\n", nt);
  }

  //------------------------------------------------------------------------------

  int Event::write_tracks(FILE *fp, const TrackVec &tracks) {
    // Returns total number of bytes written.

    int n_tracks = tracks.size();
    fwrite(&n_tracks, sizeof(int), 1, fp);

    auto start = ftell(fp);
    int data_size = 2 * sizeof(int) + n_tracks * sizeof(Track);
    fwrite(&data_size, sizeof(int), 1, fp);

    fwrite(tracks.data(), sizeof(Track), n_tracks, fp);

    for (int i = 0; i < n_tracks; ++i) {
      fwrite(tracks[i].beginHitsOnTrack(), sizeof(HitOnTrack), tracks[i].nTotalHits(), fp);
      data_size += tracks[i].nTotalHits() * sizeof(HitOnTrack);
    }

    fseek(fp, start, SEEK_SET);
    fwrite(&data_size, sizeof(int), 1, fp);
    fseek(fp, 0, SEEK_END);

    return data_size;
  }

  int Event::read_tracks(FILE *fp, TrackVec &tracks, bool skip_reading) {
    // Returns number of read tracks (negative if actual reading was skipped).

    int n_tracks, data_size;
    fread(&n_tracks, sizeof(int), 1, fp);
    fread(&data_size, sizeof(int), 1, fp);

    if (skip_reading) {
      fseek(fp, data_size - 2 * sizeof(int), SEEK_CUR);  // -2 because data_size counts itself and n_tracks too
      n_tracks = -n_tracks;
    } else {
      tracks.resize(n_tracks);

      fread(tracks.data(), sizeof(Track), n_tracks, fp);

      for (int i = 0; i < n_tracks; ++i) {
        auto &hot_vec = const_cast<std::vector<HitOnTrack> &>(tracks[i].refHitsOnTrackVector());
        memset((void *)&hot_vec, 0, sizeof(hot_vec));
        hot_vec.resize(tracks[i].nTotalHits());
        fread(hot_vec.data(), sizeof(HitOnTrack), tracks[i].nTotalHits(), fp);
      }
    }

    return n_tracks;
  }

  //------------------------------------------------------------------------------

  void Event::setInputFromCMSSW(std::vector<HitVec> hits, TrackVec seeds) {
    layerHits_ = std::move(hits);
    seedTracks_ = std::move(seeds);
  }

  //------------------------------------------------------------------------------

  int Event::clean_cms_simtracks() {
    // Sim tracks from cmssw have the following issues:
    // - hits are not sorted by layer;
    // - there are tracks with too low number of hits, even 0;
    // - even with enough hits, there can be too few layers (esp. in endcap);
    // - tracks from secondaries can have extremely low pT.
    // Possible further checks:
    // - make sure enough hits exist in seeding layers.
    //
    // What is done:
    // 1. Hits are sorted by layer;
    // 2. Non-findable tracks are marked with Track::Status::not_findable flag.
    //
    // Returns number of passed simtracks.

    dprintf("Event::clean_cms_simtracks processing %lu simtracks.\n", simTracks_.size());

    int n_acc = 0;
    [[maybe_unused]] int i = -1;
    for (Track &t : simTracks_) {
      i++;

      t.sortHitsByLayer();

      const int lyr_cnt = t.nUniqueLayers();

      //const int lasthit = t.getLastFoundHitPos();
      //const float eta = layerHits_[t.getHitLyr(lasthit)][t.getHitIdx(lasthit)].eta();

      if (lyr_cnt < Config::cmsSelMinLayers)  // || Config::TrkInfo.is_transition(eta))
      {
        dprintf("Rejecting simtrack %d, n_hits=%d, n_layers=%d, pT=%f\n", i, t.nFoundHits(), lyr_cnt, t.pT());
        t.setNotFindable();
      } else {
        dprintf("Accepting simtrack %d, n_hits=%d, n_layers=%d, pT=%f\n", i, t.nFoundHits(), lyr_cnt, t.pT());
        ++n_acc;
      }
    }

    return n_acc;
  }

  int Event::clean_cms_seedtracks(TrackVec *seed_ptr) {
    const float etamax_brl = Config::c_etamax_brl;
    const float dpt_common = Config::c_dpt_common;
    const float dzmax_brl = Config::c_dzmax_brl;
    const float drmax_brl = Config::c_drmax_brl;
    const float ptmin_hpt = Config::c_ptmin_hpt;
    const float dzmax_hpt = Config::c_dzmax_hpt;
    const float drmax_hpt = Config::c_drmax_hpt;
    const float dzmax_els = Config::c_dzmax_els;
    const float drmax_els = Config::c_drmax_els;

    const float dzmax2_inv_brl = 1.f / (dzmax_brl * dzmax_brl);
    const float drmax2_inv_brl = 1.f / (drmax_brl * drmax_brl);
    const float dzmax2_inv_hpt = 1.f / (dzmax_hpt * dzmax_hpt);
    const float drmax2_inv_hpt = 1.f / (drmax_hpt * drmax_hpt);
    const float dzmax2_inv_els = 1.f / (dzmax_els * dzmax_els);
    const float drmax2_inv_els = 1.f / (drmax_els * drmax_els);

    TrackVec &seeds = (seed_ptr != nullptr) ? *seed_ptr : seedTracks_;
    const int ns = seeds.size();

    TrackVec cleanSeedTracks;
    cleanSeedTracks.reserve(ns);
    std::vector<bool> writetrack(ns, true);

    const float invR1GeV = 1.f / Config::track1GeVradius;

    std::vector<int> nHits(ns);
    std::vector<int> charge(ns);
    std::vector<float> oldPhi(ns);
    std::vector<float> pos2(ns);
    std::vector<float> eta(ns);
    std::vector<float> ctheta(ns);
    std::vector<float> invptq(ns);
    std::vector<float> pt(ns);
    std::vector<float> x(ns);
    std::vector<float> y(ns);
    std::vector<float> z(ns);

    for (int ts = 0; ts < ns; ts++) {
      const Track &tk = seeds[ts];
      nHits[ts] = tk.nFoundHits();
      charge[ts] = tk.charge();
      oldPhi[ts] = tk.momPhi();
      pos2[ts] = std::pow(tk.x(), 2) + std::pow(tk.y(), 2);
      eta[ts] = tk.momEta();
      ctheta[ts] = 1.f / std::tan(tk.theta());
      invptq[ts] = tk.charge() * tk.invpT();
      pt[ts] = tk.pT();
      x[ts] = tk.x();
      y[ts] = tk.y();
      z[ts] = tk.z();
    }

    for (int ts = 0; ts < ns; ts++) {
      if (not writetrack[ts])
        continue;  //FIXME: this speed up prevents transitive masking; check build cost!

      const float oldPhi1 = oldPhi[ts];
      const float pos2_first = pos2[ts];
      const float Eta1 = eta[ts];
      const float Pt1 = pt[ts];
      const float invptq_first = invptq[ts];

      //#pragma simd /* Vectorization via simd had issues with icc */
      for (int tss = ts + 1; tss < ns; tss++) {
        const float Pt2 = pt[tss];

        ////// Always require charge consistency. If different charge is assigned, do not remove seed-track
        if (charge[tss] != charge[ts])
          continue;

        const float thisDPt = std::abs(Pt2 - Pt1);
        ////// Require pT consistency between seeds. If dpT is large, do not remove seed-track.
        if (thisDPt > dpt_common * (Pt1))
          continue;

        const float Eta2 = eta[tss];
        const float deta2 = std::pow(Eta1 - Eta2, 2);

        const float oldPhi2 = oldPhi[tss];

        const float pos2_second = pos2[tss];
        const float thisDXYSign05 = pos2_second > pos2_first ? -0.5f : 0.5f;

        const float thisDXY = thisDXYSign05 * sqrt(std::pow(x[ts] - x[tss], 2) + std::pow(y[ts] - y[tss], 2));

        const float invptq_second = invptq[tss];

        const float newPhi1 = oldPhi1 - thisDXY * invR1GeV * invptq_first;
        const float newPhi2 = oldPhi2 + thisDXY * invR1GeV * invptq_second;

        const float dphi = cdist(std::abs(newPhi1 - newPhi2));

        const float dr2 = deta2 + dphi * dphi;

        const float thisDZ = z[ts] - z[tss] - thisDXY * (ctheta[ts] + ctheta[tss]);
        const float dz2 = thisDZ * thisDZ;

        ////// Reject tracks within dR-dz elliptical window.
        ////// Adaptive thresholds, based on observation that duplicates are more abundant at large pseudo-rapidity and low track pT
        if (std::abs(Eta1) < etamax_brl) {
          if (dz2 * dzmax2_inv_brl + dr2 * drmax2_inv_brl < 1.0f)
            writetrack[tss] = false;
        } else if (Pt1 > ptmin_hpt) {
          if (dz2 * dzmax2_inv_hpt + dr2 * drmax2_inv_hpt < 1.0f)
            writetrack[tss] = false;
        } else {
          if (dz2 * dzmax2_inv_els + dr2 * drmax2_inv_els < 1.0f)
            writetrack[tss] = false;
        }
      }

      if (writetrack[ts])
        cleanSeedTracks.emplace_back(seeds[ts]);
    }

    seeds.swap(cleanSeedTracks);

#ifdef DEBUG
    {
      const int ns2 = seeds.size();
      printf("Number of CMS seeds before %d --> after %d cleaning\n", ns, ns2);

      for (int it = 0; it < ns2; it++) {
        const Track &ss = seeds[it];
        printf("  %3i q=%+i pT=%7.3f eta=% 7.3f nHits=%i label=% i\n",
               it,
               ss.charge(),
               ss.pT(),
               ss.momEta(),
               ss.nFoundHits(),
               ss.label());
      }
    }
#endif

    return seeds.size();
  }

  int Event::select_tracks_iter(unsigned int n) {
    if (n == 0)
      return 1;

    unsigned int algorithms[] = {4, 22, 23, 5, 24, 7, 8, 9, 10, 6};  //to be stored somewhere common

    //saving seeds by algorithm
    const int ns = seedTracks_.size();

    TrackVec cleanSeedTracks;
    cleanSeedTracks.reserve(ns);

    for (int ts = 0; ts < ns; ts++) {
      const Track &tk = seedTracks_[ts];
      unsigned int algo = (unsigned int)tk.algorithm();
      if (std::find(algorithms, algorithms + n, algo) != algorithms + n)
        cleanSeedTracks.emplace_back(seedTracks_[ts]);
    }
    seedTracks_.swap(cleanSeedTracks);

    //saving tracks by algorithm
    const int nt = cmsswTracks_.size();

    TrackVec cleanTracks;
    cleanTracks.reserve(nt);

    for (int ts = 0; ts < nt; ts++) {
      const Track &tk = cmsswTracks_[ts];
      unsigned int algo = (unsigned int)tk.algorithm();
      if (std::find(algorithms, algorithms + n, algo) != algorithms + n)
        cleanTracks.emplace_back(cmsswTracks_[ts]);
    }
    cmsswTracks_.swap(cleanTracks);
    return cmsswTracks_.size() + seedTracks_.size();
  }

  int Event::clean_cms_seedtracks_badlabel() {
    printf("***\n*** REMOVING SEEDS WITH BAD LABEL. This is a development hack. ***\n***\n");
    TrackVec buf;
    seedTracks_.swap(buf);
    std::copy_if(
        buf.begin(), buf.end(), std::back_inserter(seedTracks_), [](const Track &t) { return t.label() >= 0; });
    return seedTracks_.size();
  }

  int Event::use_seeds_from_cmsswtracks() {
    int ns = seedTracks_.size();

    TrackVec cleanSeedTracks;
    cleanSeedTracks.reserve(ns);

    for (auto &&cmsswtrack : cmsswTracks_) {
      cleanSeedTracks.emplace_back(seedTracks_[cmsswtrack.label()]);
    }

    seedTracks_.swap(cleanSeedTracks);

    return seedTracks_.size();
  }

  void Event::relabel_bad_seedtracks() {
    int newlabel = 0;
    for (auto &&track : seedTracks_) {
      if (track.label() < 0)
        track.setLabel(--newlabel);
    }
  }

  void Event::relabel_cmsswtracks_from_seeds() {
    std::map<int, int> cmsswLabelMap;
    for (size_t iseed = 0; iseed < seedTracks_.size(); iseed++) {
      for (size_t icmssw = 0; icmssw < cmsswTracks_.size(); icmssw++) {
        if (cmsswTracks_[icmssw].label() == static_cast<int>(iseed)) {
          cmsswLabelMap[icmssw] = seedTracks_[iseed].label();
          break;
        }
      }
    }
    for (size_t icmssw = 0; icmssw < cmsswTracks_.size(); icmssw++) {
      cmsswTracks_[icmssw].setLabel(cmsswLabelMap[icmssw]);
    }
  }

  //==============================================================================
  // HitMask handling
  //==============================================================================

  void Event::fill_hitmask_bool_vectors(int track_algo, std::vector<std::vector<bool>> &layer_masks) {
    // Convert from per-hit uint64_t to per layer bool-vectors for given
    // iteration.

    uint64_t iter_mask = uint64_t(1) << track_algo;

    const int n_lay = (int)layerHits_.size();
    layer_masks.resize(n_lay);

    for (int l = 0; l < n_lay; ++l) {
      const int n_hit = (int)layerHits_[l].size();
      layer_masks[l].resize(n_hit);

      for (int i = 0; i < n_hit; ++i) {
        layer_masks[l][i] = layerHitMasks_[l][i] & iter_mask;
      }
    }
  }

  void Event::fill_hitmask_bool_vectors(std::vector<int> &track_algo_vec, std::vector<std::vector<bool>> &layer_masks) {
    // Convert from per-hit uint64_t to per layer bool-vectors for a list of
    // iterations.
    // A hit mask is set if it is set for _all_ listed iterations.

    uint64_t iter_mask = 0;
    for (auto ta : track_algo_vec)
      iter_mask |= uint64_t(1) << ta;

    const int n_lay = (int)layerHits_.size();
    layer_masks.resize(n_lay);

    for (int l = 0; l < n_lay; ++l) {
      const int n_hit = (int)layerHits_[l].size();
      layer_masks[l].resize(n_hit);

      for (int i = 0; i < n_hit; ++i) {
        uint64_t hitmasks = layerHitMasks_[l][i];
        layer_masks[l][i] = ((iter_mask ^ hitmasks) & iter_mask) == 0;
      }
    }
  }

  //==============================================================================
  // Helpers for extracting information about tracks
  //==============================================================================

  Event::SimInfoFromHits Event::simInfoForTrack(const Track &s) const {
    int n_total = s.nTotalHits();
    int n_valid = 0;
    int n_pix_total = 0, n_strip_total = 0;

    struct LabelCount {
      int n_match = 0, n_pix_match = 0, n_strip_match = 0;
    };
    std::map<int, LabelCount> lab_cnt;

    // Loop over hits on track and count matches to MC truth.
    for (int hi = 0; hi < n_total; ++hi) {
      auto hot = s.getHitOnTrack(hi);
      // printf(" %d", hot.index);

      // Ignore invalid hits
      if (hot.index < 0)
        continue;

      ++n_valid;

      bool is_pixel = Config::TrkInfo[hot.layer].is_pixel();
      if (is_pixel)
        ++n_pix_total;
      else
        ++n_strip_total;

      const Hit &h = layerHits_[hot.layer][hot.index];
      int hl = simHitsInfo_[h.mcHitID()].mcTrackID_;
      // printf(" (%d)", hl);
      if (hl >= 0) {
        LabelCount &lc = lab_cnt[hl];
        ++lc.n_match;
        if (is_pixel)
          ++lc.n_pix_match;
        else
          ++lc.n_strip_match;
      }
    }

    // Find the label and hit-counts with maximum number of matches.
    int max_c = -1, max_c_pix = -1, max_c_strip = -1, max_l = -1;
    for (auto &x : lab_cnt) {
      if (x.second.n_match > max_c) {
        max_l = x.first;
        max_c = x.second.n_match;
        max_c_pix = x.second.n_pix_match;
        max_c_strip = x.second.n_strip_match;
      } else if (x.second.n_match == max_c) {
        max_l = -1;
      }
    }
    if (max_c < 0) {
      max_c = max_c_pix = max_c_strip = 0;
      max_l = -1;
    }
    // printf(" ] -> %d %d => %d\n", s.nTotalHits(), max_c, max_l);
    return { max_l, n_total, n_valid, max_c, n_pix_total, max_c_pix, n_strip_total, max_c_strip };
  }

  Event::SimInfoFromHits Event::simInfoForTrack(Track &s, bool relabel) {
    auto sifh = simInfoForTrack(s);
    if (relabel) {
      s.setLabel(sifh.label);
    }
    return sifh;
  }

  //==============================================================================

  int Event::countSimHitsInLayer(int label, int layer) const {
    if (label < 0 || label >= (int) simTracks_.size())
      return -1;
    const Track &s = simTracks_[label];
    int count = 0;
    for (int hi = 0; hi < s.nTotalHits(); ++hi) {
      HitOnTrack hot = s.getHitOnTrack(hi);
      if (hot.layer == layer && hot.index >= 0)
        ++count;
    }
    return count;
  }

  int Event::countPixelHits(const Track &track, bool inner_only) const{
    int npix = 0;
    for (int i = 0; i < track.nTotalHits(); ++i) {
      int lay = track.getHitLyr(i);
      if (lay >= 0) {
        if (Config::TrkInfo[lay].is_pixel())
          ++npix;
        else if (inner_only)
          break;
      }
    }
    return npix;
  }

  int Event::countPixelLayers(const Track &track, bool inner_only) const {
    std::set<int> layers;
    for (int i = 0; i < track.nTotalHits(); ++i) {
      int lay = track.getHitLyr(i);
      if (lay >= 0) {
        if (Config::TrkInfo[lay].is_pixel())
          layers.insert(lay);
        else if (inner_only)
          break;
      }
    }
    return layers.size();
  }

  int Event::firstInnerPixelLayer(const Track &track, int offset) const {
    if (offset < 0) throw std::runtime_error("Invalid offset in Event::firstInnerPixelLayer");

    std::vector<int> pix_layers = getInnerPixelLayers(track);
    int n_pix = pix_layers.size();
    return (offset >= n_pix) ? -1 : pix_layers[offset];
  }

  int Event::lastInnerPixelLayer(const Track &track, int offset) const {
    if (offset < 0) throw std::runtime_error("Invalid offset in Event::lastInnerPixelLayer");

    std::vector<int> pix_layers = getInnerPixelLayers(track);
    int n_pix = pix_layers.size();
    return (offset >= n_pix) ? -1 : pix_layers[n_pix - offset - 1];
  }

  std::vector<int> Event::getInnerPixelLayers(const Track &track) const {
    std::vector<int> pix_layers;
    pix_layers.reserve(6);
    int last_pix = -1;
    for (int i = 0; i < track.nTotalHits(); ++i) {
      int lay = track.getHitLyr(i);
      if (lay >= 0) {
        if (Config::TrkInfo[lay].is_pixel() ) {
          if (lay != last_pix)
            pix_layers.push_back(lay);
          last_pix = lay;
        } else {
          break;
        }
      }
    }
    return pix_layers;
  }

  int Event::countStripHits(const Track &track, bool outer_only) const {
    int nhit = 0;
    for (int i = track.nTotalHits() - 1; i >= 0 ; --i) {
      int lay = track.getHitLyr(i);
      if (lay >= 0) {
        if (Config::TrkInfo[lay].is_strip())
          ++nhit;
        else if (outer_only)
          break;
      }
    }
    return nhit;
  }

  int Event::countStripLayers(const Track &track, bool outer_only) const {
    std::set<int> layers;
    for (int i = track.nTotalHits() - 1; i >= 0 ; --i) {
      int lay = track.getHitLyr(i);
      if (lay >= 0) {
        if (Config::TrkInfo[lay].is_strip())
          layers.insert(lay);
        else if (outer_only)
          break;
      }
    }
    return layers.size();
  }

  int Event::firstInnerStripLayer(const Track &track) const {
    int first_layer = -1;
    for (int i = 0; i < track.nTotalHits(); ++i) {
      int lay = track.getHitLyr(i);
      if (lay >= 0) {
        if (Config::TrkInfo[lay].is_strip()) {
          first_layer = lay;
          break;
        }
      }
    }
    return first_layer;
  }

  //==============================================================================
  // Handling of current seed vectors and MC label reconstruction from hit data
  //==============================================================================

  void Event::setCurrentSeedTracks(const TrackVec &seeds) {
    currentSeedTracks_ = &seeds;

    currentSeedSimFromHits_.resize(currentSeedTracks_->size());
    for (int si = 0; si < (int)currentSeedTracks_->size(); ++si) {
      const Track &s = currentSeed(si);
      // printf("%3d (%d): [", si, s.label());
      currentSeedSimFromHits_[si] = simInfoForTrack(s);
    }
  }

  void Event::resetCurrentSeedTracks() {
    currentSeedTracks_ = nullptr;
    currentSeedSimFromHits_.clear();
  }

  // The following functions are needed for mkFit@HLT phase2 (and maybe beyond)

  void Event::relabelSeedTracksSequentially() {
    const int ns = seedTracks_.size();
    for (int i = 0; i < ns; ++i)
      seedTracks_[i].setLabel(i);
  }

  void Event::filterOutMislabeledHitsInSimTracks() {
    for (Track &track : simTracks_) {
      const int label = track.label();

      track.filterHits([this, label](const HitOnTrack& hot) {
        if (hot.index < 0) return false;  // Keep invalid hits

        const Hit& hit = layerHits_[hot.layer][hot.index];
        const int hit_label = simHitsInfo_[hit.mcHitID()].mcTrackID();

        return hit_label != label;  // Remove mismatched, typically hit_label == -1
      });
    }
  }

  //==============================================================================
  // Print / memUsage
  //==============================================================================

  void Event::print_tracks(const TrackVec &tracks, bool print_hits) const {
    const int nt = tracks.size();
    auto score_func = IterationConfig::get_track_scorer("default");
    //WARNING: Printouts for hits will not make any sense if mkFit is not run with a validation flag such as --quality-val
    printf("Event::print_tracks printing %d tracks %s hits:\n", nt, (print_hits ? "with" : "without"));
    for (int it = 0; it < nt; it++) {
      const Track &t = tracks[it];
      printf("  %i with q=%+i pT=%7.3f eta=% 7.3f nHits=%2d  label=%4d findable=%d score=%7.3f chi2=%7.3f\n",
             it,
             t.charge(),
             t.pT(),
             t.momEta(),
             t.nFoundHits(),
             t.label(),
             t.isFindable(),
             getScoreCand(score_func, t),
             t.chi2());

      if (print_hits) {
        for (int ih = 0; ih < t.nTotalHits(); ++ih) {
          int lyr = t.getHitLyr(ih);
          int idx = t.getHitIdx(ih);
          if (idx >= 0) {
            const Hit &hit = layerHits_[lyr][idx];
            printf("    hit %2d lyr=%2d idx=%3d pos r=%7.3f z=% 8.3f   mc_hit=%3d mc_trk=%3d\n",
                   ih,
                   lyr,
                   idx,
                   layerHits_[lyr][idx].r(),
                   layerHits_[lyr][idx].z(),
                   hit.mcHitID(),
                   hit.mcTrackID(simHitsInfo_));
          } else
            printf("    hit %2d lyr=%2d idx=%3d\n", ih, t.getHitLyr(ih), t.getHitIdx(ih));
        }
      }
    }
  }

  size_t Event::memUsage() const {
    size_t total = 0;

    // Layer hits: Iterate over layers to count hits in each buffer
    for (const auto& hv : layerHits_) {
      total += hv.capacity() * sizeof(Hit); // Sum capacity of hits across all layers
    }

    // Track vectors
    total += simTracks_.size() * sizeof(Track);
    total += seedTracks_.size() * sizeof(Track);
    total += candidateTracks_.size() * sizeof(Track);
    total += fitTracks_.size() * sizeof(Track);

    // Trace data (if MKFIT_TRACE is defined)
    #ifdef MKFIT_TRACE
      total += trCandMetas_.size() * sizeof(TrCandMeta);
      total += trCandStates_.size() * sizeof(TrCandState);
      total += trHitMatches_.size() * sizeof(TrHitMatch);
      total += trKalmanUpdates_.size() * sizeof(TrKalmanUpdate);
      total += trBkFitUpdates_.size() * sizeof(TrBkFitUpdate);
    #endif

    return total;
  }

  void Event::printMemUsage() const {
    printf("Event Memory Breakdown:\n");
    printf("========================\n");

    // Calculate hits size dynamically
    size_t totalHits = 0;
    for (const auto& hv : layerHits_)
      totalHits += hv.size();
    printf("Layer Hits:  %5zu * %3zu = %zu bytes\n", totalHits, sizeof(Hit), totalHits * sizeof(Hit));

    printf("Sim tracks:  %5zu * %3zu = %zu bytes\n", simTracks_.size(), sizeof(Track), simTracks_.size() * sizeof(Track));
    printf("Seed tracks: %5zu * %3zu = %zu bytes\n", seedTracks_.size(), sizeof(Track), seedTracks_.size() * sizeof(Track));
    printf("Cand tracks: %5zu * %3zu = %zu bytes\n", candidateTracks_.size(), sizeof(Track), candidateTracks_.size() * sizeof(Track));
    printf("Fit tracks:  %5zu * %3zu = %zu bytes\n", fitTracks_.size(), sizeof(Track), fitTracks_.size() * sizeof(Track));

    #ifdef MKFIT_TRACE
      printf("Trace cand metas:     %6zu * %3zu = %zu bytes\n", trCandMetas_.size(), sizeof(TrCandMeta), trCandMetas_.size() * sizeof(TrCandMeta));
      printf("Trace cand states:    %6zu * %3zu = %zu bytes\n", trCandStates_.size(), sizeof(TrCandState), trCandStates_.size() * sizeof(TrCandState));
      printf("Trace hit matches:    %6zu * %3zu = %zu bytes\n", trHitMatches_.size(), sizeof(TrHitMatch), trHitMatches_.size() * sizeof(TrHitMatch));
      printf("Trace kalman updates: %6zu * %3zu = %zu bytes\n", trKalmanUpdates_.size(), sizeof(TrKalmanUpdate), trKalmanUpdates_.size() * sizeof(TrKalmanUpdate));
      printf("Trace bkfit updates: %6zu * %3zu = %zu bytes\n", trBkFitUpdates_.size(), sizeof(TrBkFitUpdate), trBkFitUpdates_.size() * sizeof(TrBkFitUpdate));
    #endif

    printf("========================\n");
    printf("Total estimated size (without vector overhead): %zu bytes\n", memUsage());
  }

  //==============================================================================
  // Trace and RDF stuff
  //==============================================================================

#ifdef MKFIT_TRACE

  void Event::build_trace_maps_etc() {
    trChildrenByState_.clear();
    trRootCands_.clear();
    for (auto& c : trCandStates_) {
        if (c.parent_id >= 0)
          trChildrenByState_[c.parent_id].push_back(c.id);
        else
          trRootCands_.push_back(c.id);
    }
    for (auto &s : trLayerSearches_) {
      trCandStates_[s.state_id].search_id = s.id;
    }
    trHitMatchesByState_.clear();
    for (auto& h : trHitMatches_) {
      trHitMatchesByState_[h.state_id].push_back(h.id);
    }
    trKalmanUpdatesByState_.clear();
    for (auto& k : trKalmanUpdates_) {
      trKalmanUpdatesByState_[k.state_id_in].push_back(k.id);
    }
    trBkFitUpdatesByState_.clear();
    for (auto& b : trBkFitUpdates_) {
      trBkFitUpdatesByState_[b.state_id_in].push_back(b.id);
    }

    int msize = trCandMetas_.size();
    trSIFHforSeedByMeta_.clear();
    trSIFHforSeedByMeta_.resize(msize);
    trSIFHforCandByMeta_.clear();
    trSIFHforCandByMeta_.resize(msize);

    for (int i = 0; i < (int)trCandMetas_.size(); ++i) {
      auto& c = trCandMetas_[i];
      trSIFHforSeedByMeta_[i] = simInfoForCurrentSeed(c.seed);
      // c.cand stays -1 for a meta whose candidate never got exported -- it was
      // filtered out, or the driver did not publish its tracks into
      // candidateTracks_ at all. candidateTracks_[] is unchecked, so guard it.
      if (c.cand >= 0 && c.cand < (int)candidateTracks_.size())
        trSIFHforCandByMeta_[i] = simInfoForTrack(candidateTracks_[c.cand]);
      c.sim = trSIFHforSeedByMeta_[i].label;
    }

    // copy currentSeeds out, they might be temporary or modified
    trSeeds_ = *currentSeedTracks_;
  }

#endif

  //==============================================================================
  // DataFile
  //==============================================================================

  int DataFile::openRead(const std::string &fname, int expected_n_layers,
                         const std::string &expected_geom_version) {
    constexpr int min_ver = 7;  // v7/v8 files stay readable; they simply carry
    constexpr int max_ver = 9;  // no geometry stamp, so none is compared.

    f_fp = fopen(fname.c_str(), "r");
    assert(f_fp != 0 && "Opening of input file failed.");

    // Magic and version are the first two ints of the file so that the rest can
    // be read knowing what it is. Read those two alone, decide, then read the
    // header for THAT version.
    int magic = 0, version = 0;
    fread(&magic, sizeof(int), 1, f_fp);
    fread(&version, sizeof(int), 1, f_fp);

    if (magic != 0xBEEF) {
      fprintf(stderr, "Incompatible input file (wrong magick).\n");
      exit(1);
    }
    if (version < min_ver || version > max_ver) {
      fprintf(stderr,
              "Unsupported file version %d. Supported versions are from %d to %d.\n",
              version,
              min_ver,
              max_ver);
      exit(1);
    }

    const size_t hdr_size = DataFileHeader::size_of_version(version);
    f_header = DataFileHeader();  // f_geom_version stays empty unless the file has one
    fseek(f_fp, 0, SEEK_SET);
    fread(&f_header, hdr_size, 1, f_fp);

    f_data_start = hdr_size;
    f_pos = f_data_start;

    if (f_header.f_sizeof_track != sizeof(Track)) {
      fprintf(stderr,
              "sizeof(Track) on file (%d) different from current value (%d).\n",
              f_header.f_sizeof_track,
              (int)sizeof(Track));
      exit(1);
    }
    if (f_header.f_sizeof_hit != sizeof(Hit)) {
      fprintf(stderr,
              "sizeof(Hit) on file (%d) different from current value (%d).\n",
              f_header.f_sizeof_hit,
              (int)sizeof(Hit));
      exit(1);
    }
    if (f_header.f_sizeof_hot != sizeof(HitOnTrack)) {
      fprintf(stderr,
              "sizeof(HitOnTrack) on file (%d) different from current value (%d).\n",
              f_header.f_sizeof_hot,
              (int)sizeof(HitOnTrack));
      exit(1);
    }
    if (f_header.f_n_layers != expected_n_layers) {
      fprintf(stderr,
              "Number of layers on file (%d) is different from current TrackerInfo (%d).\n",
              f_header.f_n_layers,
              expected_n_layers);
      exit(1);
    }

    // The check the SAMPLE PROVENANCE disaster needed and did not have: a stale
    // sample and the current geometry agree on n_layers (both phase-2 geometries
    // have 60), so nothing caught it. The names differ, and that is cheap.
    if (!expected_geom_version.empty() && f_header.f_geom_version[0] != '\0' &&
        expected_geom_version != f_header.f_geom_version) {
      fprintf(stderr,
              "GEOMETRY MISMATCH: sample '%s' was written against geometry '%s', "
              "but the loaded geometry is '%s'.\n"
              "Module short-ids are assigned per geometry, so every module-frame "
              "quantity from this pairing would be meaningless.\n",
              fname.c_str(), f_header.f_geom_version, expected_geom_version.c_str());
      exit(1);
    }

    printf("Opened file '%s', format version %d, n_layers %d, n_events %d\n",
           fname.c_str(),
           f_header.f_format_version,
           f_header.f_n_layers,
           f_header.f_n_events);
    if (f_header.f_extra_sections) {
      if (f_header.f_geom_version[0] != '\0')
        printf("  Geometry: %s\n", f_header.f_geom_version);
      printf("  Extra sections:");
      if (f_header.f_extra_sections & ES_SimTrackStates)
        printf(" SimTrackStates");
      if (f_header.f_extra_sections & ES_Seeds)
        printf(" Seeds");
      if (f_header.f_extra_sections & ES_CmsswTracks)
        printf(" CmsswTracks");
      if (f_header.f_extra_sections & ES_HitIterMasks)
        printf(" HitIterMasks");
      if (f_header.f_extra_sections & ES_BeamSpot)
        printf(" BeamSpot");
      if (f_header.f_extra_sections & ES_SimHitStates)
        printf(" SimHitStates%s", Config::readSimHitStates ? "" : "(skipped)");
      printf("\n");
    }

    if (Config::seedInput == cmsswSeeds && !hasSeeds()) {
      fprintf(stderr, "Reading of CmsswSeeds requested but data not available on file.\n");
      exit(1);
    }

    if (Config::readCmsswTracks && !hasCmsswTracks()) {
      fprintf(stderr, "Reading of CmsswTracks requested but data not available on file.\n");
      exit(1);
    }

    return f_header.f_n_events;
  }

  void DataFile::openWrite(const std::string &fname, int n_layers, int n_ev, int extra_sections,
                           const std::string &geom_version) {
    f_fp = fopen(fname.c_str(), "w");
    f_header.f_n_layers = n_layers;
    f_header.f_n_events = n_ev;
    f_header.f_extra_sections = extra_sections;

    if (!geom_version.empty()) {
      std::strncpy(f_header.f_geom_version, geom_version.c_str(), sizeof(f_header.f_geom_version) - 1);
      f_header.f_geom_version[sizeof(f_header.f_geom_version) - 1] = '\0';
    }
    fwrite(&f_header, sizeof(DataFileHeader), 1, f_fp);

    f_data_start = ftell(f_fp);
    f_pos = f_data_start;
  }

  void DataFile::rewind() {
    std::lock_guard<std::mutex> readlock(f_next_ev_mutex);
    f_pos = f_data_start;
    fseek(f_fp, f_pos, SEEK_SET);
  }

  int DataFile::advancePosToNextEvent(FILE *fp) {
    int evsize;

    std::lock_guard<std::mutex> readlock(f_next_ev_mutex);

    fseek(fp, f_pos, SEEK_SET);
    fread(&evsize, sizeof(int), 1, fp);
    if (Config::loopOverFile) {
      // File ended, rewind back to beginning
      if (feof(fp) != 0) {
        f_pos = f_data_start;
        fseek(fp, f_pos, SEEK_SET);
        fread(&evsize, sizeof(int), 1, fp);
      }
    }

    f_pos += evsize;

    return evsize;
  }

  void DataFile::skipNEvents(int n_to_skip) {
    int evsize;

    std::lock_guard<std::mutex> readlock(f_next_ev_mutex);

    while (n_to_skip-- > 0) {
      fseek(f_fp, f_pos, SEEK_SET);
      fread(&evsize, sizeof(int), 1, f_fp);
      f_pos += evsize;
    }
  }

  void DataFile::close() {
    fclose(f_fp);
    f_fp = 0;
    f_header = DataFileHeader();
  }

  void DataFile::CloseWrite(int n_written) {
    if (f_header.f_n_events != n_written) {
      fseek(f_fp, 0, SEEK_SET);
      f_header.f_n_events = n_written;
      fwrite(&f_header, sizeof(DataFileHeader), 1, f_fp);
    }
    close();
  }

  //==============================================================================
  // Misc debug / printout
  //==============================================================================
  // clang-format off

  void print(std::string pfx, int itrack, const Track &trk, const Event &ev) {
    std::cout << pfx << ": " << itrack << " hits: " << trk.nFoundHits() << " label: " << trk.label()
              << " algo: " << trk.algoint() << "\n";
    print(trk.state());

    for (int i = 0; i < trk.nTotalHits(); ++i) {
      auto hot = trk.getHitOnTrack(i);
      printf("  %2d: lyr %2d idx %5d", i, hot.layer, hot.index);
      if (hot.index >= 0) {
        auto &h = ev.layerHits_[hot.layer][hot.index];
        int hl = ev.simHitsInfo_[h.mcHitID()].mcTrackID_;
        printf("  %4d  x=%8.3f y=%8.3f z=%8.3f r=%8.3f | e_z=%8.3g e_r=%8.3g | pix=%d str=%d brl=%d\n",
                hl, h.x(), h.y(), h.z(), h.r(),
                std::sqrt(h.ezz()),
                std::sqrt(getRadErr2(h.x(), h.y(), h.exx(), h.eyy(), h.exy())),
                Config::TrkInfo[hot.layer].is_pixel(),
                Config::TrkInfo[hot.layer].is_stereo(),
                Config::TrkInfo[hot.layer].is_barrel()
                );
      } else {
        printf("\n");
      }
    }
  }

  void print(std::string pfx, int itrack, const Track &trk, int hit_begin, int hit_end, const Event &ev) {
    std::cout << pfx << ": " << itrack << " hits: " << trk.nFoundHits() << " label: " << trk.label()
              << " algo: " << trk.algoint() << "\n";
    print(trk.state());

    for (int i = hit_begin; i < hit_end; ++i) {
      auto hot = trk.getHitOnTrack(i);
      printf("  %2d: lyr %2d idx %5d", i, hot.layer, hot.index);
      if (hot.index >= 0) {
        auto &h = ev.layerHits_[hot.layer][hot.index];
        int hl = ev.simHitsInfo_[h.mcHitID()].mcTrackID_;
        printf("  %4d  x=%8.3f y=%8.3f z=%8.3f r=%8.3f | e_z=%8.3g e_r=%8.3g | pix=%d str=%d brl=%d\n",
                hl, h.x(), h.y(), h.z(), h.r(),
                std::sqrt(h.ezz()),
                std::sqrt(getRadErr2(h.x(), h.y(), h.exx(), h.eyy(), h.exy())),
                Config::TrkInfo[hot.layer].is_pixel(),
                Config::TrkInfo[hot.layer].is_stereo(),
                Config::TrkInfo[hot.layer].is_barrel()
                );
      } else {
        printf("\n");
      }
    }
  }

  void print(std::string pfx, const TrackVec &tvec, const Event &ev) {
    int nt = tvec.size();
    for (int i = 0; i < nt; ++i) {
      print(pfx, i, tvec[i], ev);
    }
  }

  void print(std::string pfx, const Event::SimInfoFromHits &si) {
    printf("%s: label=%5d n_valid=%2d n_match=%2d n_invalid=%2d good_frac()=%.4f n_pix=%2d n_pix_match=%2d n_strip=%2d n_strip_match=%2d\n",
           pfx.c_str(), si.label, si.n_valid, si.n_match, si.n_invalid(), si.good_frac(),
           si.n_pix, si.n_pix_match, si.n_strip, si.n_strip_match);
  }

#ifdef MKFIT_TRACE

  std::string format(const ::EVec3 &v, int width, int prec, char feg) {
    if (feg == 'f')
      return std::format("({0: {3}.{4}f}, {1: {3}.{4}f}, {2: {3}.{4}f}; {5: {3}.{4}f})", v.fX, v.fY, v.fZ, width, prec, v.R());
    else if(feg == 'e')
      return std::format("({0: {3}.{4}e}, {1: {3}.{4}e}, {2: {3}.{4}e}; {5: {3}.{4}e})", v.fX, v.fY, v.fZ, width, prec, v.R());
    else if(feg == 'g')
      return std::format("({0: {3}.{4}g}, {1: {3}.{4}g}, {2: {3}.{4}g}; {5: {3}.{4}g})", v.fX, v.fY, v.fZ, width, prec, v.R());
    else
      return std::format("({0: {3}.{4}}, {1: {3}.{4}}, {2: {3}.{4}}; {5: {3}.{4}})", v.fX, v.fY, v.fZ, width, prec, v.R());
  }

  void print(std::string prefix, const ::EBiVec3 &s, const std::string postfix) {
    printf("%s: pos=(% 8.3f,% 8.3f,% 8.3f)  mom=(% 8.3f,% 8.3f,% 8.3f)%s",
           prefix.c_str(), s.pos.fX, s.pos.fY, s.pos.fZ, s.mom.fX, s.mom.fY, s.mom.fZ, postfix.c_str());
  }

  void print(std::string pfx, const TrCandMeta &cm, const Event *ev) {
    printf("%s: id=%d event=%d seed=%d global_seed=%d sim=%d cand=%d "
      "stage_id_0=%d stage_id_1=%d stage_id_2=%d \n",
      pfx.c_str(), cm.id, cm.event,
      cm.seed, cm.global_seed, cm.sim, cm.cand,
      cm.stage_ids[0], cm.stage_ids[1], cm.stage_ids[2]);
      if (ev != nullptr) {
        print("seed_sim_info", ev->trSIFHforSeedByMeta_[cm.id]);
        print("cand_sim_info", ev->trSIFHforCandByMeta_[cm.id]);
      }
  }

  void print(std::string pfx, const TrCandStage &cs) {
    printf("%s: id=%d meta_id=%d parent_stage_id=%d stage=%d "
      "root_state_id=%d final_state_id=%d\n",
      pfx.c_str(), cs.id, cs.meta_id, cs.parent_stage_id, cs.stage,
      cs.root_state_id, cs.final_state_id);
  }

  void print(std::string pfx, const TrCandState &cs) {
    printf("%s: id=%d parent_id=%d meta_id=%d stage_id=%d layer=%d step=%d "
        "has_children=%d on_final_path=%d\n",
        pfx.c_str(), cs.id, cs.parent_id, cs.meta_id, cs.stage_id,
        cs.layer, cs.step, cs.has_children, cs.on_final_path);
        print("    kine", cs.kine);
  }

  void print(std::string pfx, const TrLayerSearch &ls) {
    printf("%s: id=%d state_id=%d layer=%d layer_sec=%d %s %s\n",
           pfx.c_str(), ls.id, ls.state_id, ls.layer, ls.layer_sec,
           ls.is_barrel ? "barrel" : "endcap", ls.is_outward ? "outward" : "inward");
    printf("    entry=%s dalpha=% f fail=%d\n"
           "    exit =%s dalpha=% f fail=%d\n",
           format(ls.prop_entry.pos, 10, 8, ' ').c_str(), ls.prop_entry.dalpha, ls.prop_entry.fail_flag,
           format(ls.prop_exit.pos, 10, 8, ' ').c_str(), ls.prop_exit.dalpha, ls.prop_exit.fail_flag);
    printf("    phi_c=% f phi_d=% f dphi_track=% f | q_c=% f q=[% f,% f] dq_track=% f\n",
           ls.phi_center, ls.phi_delta, ls.dphi_track,
           ls.q_center(), ls.q_min, ls.q_max, ls.dq_track);
    printf("    cov xx=% .4g xy=% .4g yy=% .4g zz=% .4g\n",
           ls.cov_xx, ls.cov_xy, ls.cov_yy, ls.cov_zz);
    printf("    bins p=[%u,%u) q=[%u,%u)", ls.p1, ls.p2, ls.q1, ls.q2);
    if (ls.layer_sec >= 0)
      printf(" | sec p=[%u,%u) q=[%u,%u)", ls.p1_sec, ls.p2_sec, ls.q1_sec, ls.q2_sec);
    printf("\n    hits scanned=%d masked=%d presel=%d pqueue=%d\n",
           ls.n_hits_scanned, ls.n_hits_masked, ls.n_hits_presel, ls.n_hits_pqueue);
  }

  void print(std::string pfx, const TrHitMatch &hm) {
    printf("%s: id=%d state_id=%d search_id=%d layer=%d hit=%d mc_match=%d "
          "score=% f dphi=% f dq=% f passed_preselect=%d "
          "rank=%d passed_pqueue=%d kalman_id=%d\n",
          pfx.c_str(), hm.id, hm.state_id, hm.search_id, hm.layer,
          hm.hit, hm.mc_match, hm.score, hm.dphi, hm.dq, hm.passed_preselect,
          hm.sub_rank, hm.passed_pqueue, hm.kalman_id);
    printf("    t_hermite=% f d_plane_h3=% .3e\n", hm.t_hermite, hm.d_plane_h3);
    printf("    residual_xyz=(res_x=% f res_y=% f res_z=% f); ", hm.residual_x, hm.residual_y, hm.residual_z);
    print("kine", hm.kine_on_plane);
  }

  void print(std::string pfx, const TrKalmanUpdate &ku) {
    printf("%s: id=%d hit_match_id=%d state_id_in=%d state_id_out=%d  "
            "chi2=% f chi2_trk=% f accepted=%d\n",
            pfx.c_str(), ku.id, ku.hit_match_id, ku.state_id_in, ku.state_id_out,
            ku.chi2, ku.chi2_trk, ku.accepted);
#ifdef MKFIT_TRACE_KALMAN_DEBUG
    print("    prop ", ku.propagated_state);
    print("    updt ", ku.updated_state);
#endif
  }

#endif

}  // end namespace mkfit
