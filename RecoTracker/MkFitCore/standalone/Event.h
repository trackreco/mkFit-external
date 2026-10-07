#ifndef RecoTracker_MkFitCore_standalone_Event_h
#define RecoTracker_MkFitCore_standalone_Event_h

#include "RecoTracker/MkFitCore/interface/Config.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/interface/BeamSpot.h"
#include "Validation.h"

#ifdef MKFIT_TRACE
#include "RecoTracker/MkFitCore/standalone/DataFormats/RntStructs.h"
#endif

#include <mutex>

namespace mkfit {

  //==============================================================================
  // SimHitState
  //==============================================================================

  // Truth state at a SIM HIT: position and momentum, and nothing else.
  //
  // Indexed by mcHitID, i.e. parallel to Event::simHitsInfo_ -- the same
  // convention the legacy Event::simTrackStates_ uses.
  //
  // Deliberately NOT a TrackState, which is 112 B:
  //  - no covariance. A Geant truth state has none, and the sim TRACK's own
  //    covariance is a documented placeholder (err(i,i) = value^2, a 100 %
  //    relative error, singular at the origin). Writing it would be 84 bytes of
  //    zeros per sim hit.
  //  - no charge. It is a per-TRACK property, reachable from the same index as
  //      simTracks_[ simHitsInfo_[mcHitID].mcTrackID() ].charge()
  //
  // That is 24 B against 112, which is what makes the section affordable: the
  // April PU sample carries 383k-471k sim hits per event, so the full TrackState
  // form would be 49 MB/event and more than double the file, against 10.5
  // MB/event (+30 %) for this.
  //
  // INVALID IS ZERO MOMENTUM. There is no separate valid flag -- a real sim hit
  // never has |p| = 0, so `mom` all-zero means "no truth state for this hit",
  // which is the case for a rec hit whose sim link was not established. Test it
  // with is_valid() rather than by reading mcTrackID, because the two are not
  // equivalent: bestTkIdx() can clear the track link while the sim hit itself is
  // perfectly well defined (see bestTkIdx() in WriteMemoryFile.cc).
  struct SimHitState {
    SVector3 pos;
    SVector3 mom;

    SimHitState() : pos(0.f, 0.f, 0.f), mom(0.f, 0.f, 0.f) {}
    SimHitState(const SVector3 &p, const SVector3 &m) : pos(p), mom(m) {}
    SimHitState(float x, float y, float z, float px, float py, float pz)
        : pos(x, y, z), mom(px, py, pz) {}

    bool is_valid() const { return mom[0] != 0.f || mom[1] != 0.f || mom[2] != 0.f; }

    float x() const { return pos[0]; }
    float y() const { return pos[1]; }
    float z() const { return pos[2]; }
    float px() const { return mom[0]; }
    float py() const { return mom[1]; }
    float pz() const { return mom[2]; }

    float r() const { return std::hypot(pos[0], pos[1]); }
    float pT() const { return std::hypot(mom[0], mom[1]); }
    float p() const { return std::sqrt(mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]); }
    float momPhi() const { return std::atan2(mom[1], mom[0]); }
    float momEta() const {
      const float pt = pT();
      return std::log((p() + mom[2]) / (pt > 0.f ? pt : 1e-9f));
    }
  };

  typedef std::vector<SimHitState> SHSVec;

  struct DataFile;

  class Event {
  public:
    explicit Event(int evtID, int nLayers);
    Event(Validation &v, int evtID, int nLayers);

    void reset(int evtID);
    void validate();
    void printStats(const TrackVec &, TrackExtraVec &);

    int evtID() const { return evtID_; }
    void resetLayerHitMap(bool resetSimHits);

    void write_out(DataFile &data_file);
    void read_in(DataFile &data_file, FILE *in_fp = 0);
    int write_tracks(FILE *fp, const TrackVec &tracks);
    int read_tracks(FILE *fp, TrackVec &tracks, bool skip_reading = false);

    void setInputFromCMSSW(std::vector<HitVec> hits, TrackVec seeds);

    int use_seeds_from_cmsswtracks();  //special mode --> use only seeds which generated cmssw reco track
    int clean_cms_simtracks();
    int clean_cms_seedtracks(
        TrackVec *seed_ptr = nullptr);    //operates on seedTracks_; returns the number of cleaned seeds
    int clean_cms_seedtracks_badlabel();  //operates on seedTracks_, removes those with label == -1;
    void relabel_bad_seedtracks();
    void relabel_cmsswtracks_from_seeds();

    int select_tracks_iter(unsigned int n = 0);  //for cmssw input

    void fill_hitmask_bool_vectors(int track_algo, std::vector<std::vector<bool>> &layer_masks);
    void fill_hitmask_bool_vectors(std::vector<int> &track_algo_vec, std::vector<std::vector<bool>> &layer_masks);

    Validation &validation_;

    // For seed access in deep data dumpers.
    struct SimInfoFromHits {
      int label = -1, n_hits = 0, n_valid = 0, n_match = 0;
      int n_pix = 0, n_pix_match = 0;
      int n_strip = 0, n_strip_match = 0;
      float good_frac() const { return (float)n_match / n_valid; }
      int n_invalid() const { return n_hits - n_valid; }
      int n_pix_bad() const { return n_pix - n_pix_match; }
      int n_strip_bad() const { return n_strip - n_strip_match; }
      bool is_set() const { return label >= 0; }
    };
    SimInfoFromHits simInfoForTrack(const Track &s) const;
    SimInfoFromHits simInfoForTrack(Track &s, bool relabel);

    int countSimHitsInLayer(int label, int layer) const;

    int countPixelHits(const Track &track, bool inner_only) const;
    int countInnerPixelHits(const Track &track) const { return countPixelHits(track, true); }
    int countAllPixelHits(const Track &track) const { return countPixelHits(track, false); }

    int countPixelLayers(const Track &track, bool inner_only) const;
    int countInnerPixelLayers(const Track &track) const { return countPixelLayers(track, true); }
    int countAllPixelLayers(const Track &track) const { return countPixelLayers(track, false); }

    int firstInnerPixelLayer(const Track &track, int offset = 0) const;
    int lastInnerPixelLayer(const Track &track, int offset = 0) const;
    std::vector<int> getInnerPixelLayers(const Track &track) const;

    int countStripHits(const Track &track, bool outer_only) const;
    int countOuterStripHits(const Track &track) const { return countStripHits(track, true); }
    int countAllStripHits(const Track &track) const { return countStripHits(track, false); }

    int countStripLayers(const Track &track, bool outer_only) const;
    int countOuterStripLayers(const Track &track) const { return countStripLayers(track, true); }
    int countAllStripLayers(const Track &track) const { return countStripLayers(track, false); }

    int firstInnerStripLayer(const Track &track) const;

    void setCurrentSeedTracks(const TrackVec &seeds);
    void resetCurrentSeedTracks();
    const Track &currentSeed(int i) const { return (*currentSeedTracks_)[i]; }
    SimInfoFromHits simInfoForCurrentSeed(int i) const { return currentSeedSimFromHits_[i]; }
    const TrackVec& currentSeedTracks() const { return *currentSeedTracks_; }

    void relabelSeedTracksSequentially();
    void filterOutMislabeledHitsInSimTracks();

    void print_tracks(const TrackVec &tracks, bool print_hits) const;

    size_t memUsage() const;
    void printMemUsage() const;

  private:
    int evtID_;

  public:
    BeamSpot beamSpot_;  // XXXX Read/Write of BeamSpot + file-version bump or extra-section to be added.
    std::vector<HitVec> layerHits_;
    std::vector<std::vector<uint64_t>> layerHitMasks_;  //aligned with layerHits_
    MCHitInfoVec simHitsInfo_;

    TrackVec simTracks_, seedTracks_, candidateTracks_, fitTracks_;
    TrackVec cmsswTracks_;
    // validation sets these, so needs to be mutable
    mutable TrackExtraVec simTracksExtra_, seedTracksExtra_, candidateTracksExtra_, fitTracksExtra_;
    mutable TrackExtraVec cmsswTracksExtra_;

    TSVec simTrackStates_;
    // Truth state per SIM HIT, indexed by mcHitID (parallel to simHitsInfo_).
    // Optional on both sides: written only with --write-sim-hit-states, and read
    // only with --read-sim-hit-states (otherwise the section is seeked past, so
    // it costs nothing but disk). See SimHitState in TrackState.h for why it is
    // not a TrackState.
    SHSVec simHitStates_;

    const TrackVec *currentSeedTracks_ = nullptr;
    mutable std::vector<SimInfoFromHits> currentSeedSimFromHits_;

  #ifdef MKFIT_TRACE
    // Not thread safe within event, multiple Events ok.
    mutable std::vector<TrCandMeta> trCandMetas_;
    mutable std::vector<TrCandStage> trCandStages_;
    mutable std::vector<TrCandState> trCandStates_;
    mutable std::vector<TrLayerSearch> trLayerSearches_;
    mutable std::vector<TrHitMatch> trHitMatches_;
    mutable std::vector<TrKalmanUpdate> trKalmanUpdates_;
    mutable std::vector<TrBkFitUpdate> trBkFitUpdates_;

    mutable TrackVec trSeeds_;
    /* *** for multiple iteration tracing ***
      could be made std::vector<TrackVec> trSeedsPerIter_;
      but then would also need:
      struct IterTraceInfo {
        int meta_begin = 0;
        int state_begin = 0;
        int seed_begin = 0;  // Index into trSeedsPerIter_
        int iteration_idx = 0;
        int search_direction = 0;
      };
      mutable std::vector<IterTraceInfo> trIterInfos_;
    */

    mutable SeedVecInsp seedVecInsp_; // describe mixture of input seeds for HLT setup

    TrCandMeta& tr_candmeta(int i) const { return trCandMetas_[i]; }
    TrCandStage& tr_candstage(int i) const { return trCandStages_[i]; }
    TrCandState& tr_candstate(int i) const { return trCandStates_[i]; }
    TrLayerSearch& tr_layersearch(int i) const { return trLayerSearches_[i]; }
    TrHitMatch& tr_hitmatch(int i) const { return trHitMatches_[i]; }
    TrKalmanUpdate& tr_kalmanupdate(int i) const { return trKalmanUpdates_[i]; }
    TrBkFitUpdate& tr_bkfitupdate(int i) const { return trBkFitUpdates_[i]; }

    TrCandMeta& trace_candmeta(TrCandMeta && cm) const {
      int s = trCandMetas_.size();
      auto &t = trCandMetas_.emplace_back(cm);
      t.id = s;
      return t;
    }
    TrCandStage& trace_candstage(TrCandStage && cs) const {
      int s = trCandStages_.size();
      auto &t = trCandStages_.emplace_back(cs);
      t.id = s;
      return t;
    }
    TrCandState& trace_candstate(TrCandState && cs) const {
      int s = trCandStates_.size();
      auto &t = trCandStates_.emplace_back(cs);
      t.id = s;
      return t;
    }
    TrLayerSearch &trace_layersearch(TrLayerSearch && ls) const {
      int s = trLayerSearches_.size();
      auto &t = trLayerSearches_.emplace_back(ls);
      t.id = s;
      return t;
    }
    TrHitMatch& trace_hitmatch(TrHitMatch && hm) const {
      int s = trHitMatches_.size();
      auto &t = trHitMatches_.emplace_back(hm);
      t.id = s;
      return t;
    }
    TrKalmanUpdate& trace_kalmanupdate(TrKalmanUpdate && ku) const {
      int s = trKalmanUpdates_.size();
      auto &t = trKalmanUpdates_.emplace_back(ku);
      t.id = s;
      return t;
    }
    TrBkFitUpdate& trace_bkfitupdate(TrBkFitUpdate && bu) const {
      int s = trBkFitUpdates_.size();
      auto &t = trBkFitUpdates_.emplace_back(bu);
      t.id = s;
      return t;
    }

    int trace_new_cand_meta(int event, int seed_index) const {
      auto &cm = trace_candmeta({ -1, event, seed_index });
      return cm.id;
    }
    // always do stage and initial state together
    // int trace_new_cand_stage(int meta_id, int parent_stage_id, int stage) {
    //   auto &cstg = trace_candstage({ -1, meta_id, parent_stage_id, stage });
    //   return cstg.id;
    // }
    std::pair<int,int>
    trace_new_cand_stage_and_state(int meta_id, int parent_stage_id, int stage, int layer, const EBiVec3 &kine, const TrackState &state) const {
      assert(stage >= 0 && stage <= 2 && "stage expected to be between 0 and 2");
      auto &cstage = trace_candstage({ -1, meta_id, parent_stage_id, stage });
      auto &cstate = trace_candstate({ -1, -1, meta_id, cstage.id, layer, 0, kine, state });
      cstage.root_state_id = cstate.id;
      return { cstage.id, cstate.id };
    }
    int trace_new_cand_state(int parent_state_id, int layer, const EBiVec3 &kine, const TrackState &state) const {
      auto &pcs = trCandStates_[parent_state_id];
      pcs.has_children = true;
      auto &cs = trace_candstate({ -1, parent_state_id, pcs.meta_id, pcs.stage_id, layer, pcs.step + 1, kine, state });
      return cs.id;
    }
    int trace_new_kalman_update(int hit_match_id, int state_id_in, float chi2, float chi2_trk) const {
      auto &ku = trace_kalmanupdate({ -1, hit_match_id, state_id_in, -1, chi2, chi2_trk });
      return ku.id;
    }

    // Aggregators, maps
    void build_trace_maps_etc();

    std::vector<int> trRootCands_;
    std::unordered_map<int, std::vector<int>> trChildrenByState_;
    std::unordered_map<int, std::vector<int>> trHitMatchesByState_;
    std::unordered_map<int, std::vector<int>> trKalmanUpdatesByState_;
    std::unordered_map<int, std::vector<int>> trBkFitUpdatesByState_;

    std::vector<SimInfoFromHits> trSIFHforSeedByMeta_;
    std::vector<SimInfoFromHits> trSIFHforCandByMeta_;
  #endif

    static std::mutex printmutex;
  };

  typedef std::vector<Event> EventVec;

  struct DataFileHeader {
    int f_magic = 0xBEEF;
    int f_format_version = 9;  //v9 adds f_geom_version; v8 added ES_SimHitStates
    int f_sizeof_track = sizeof(Track);
    int f_sizeof_hit = sizeof(Hit);
    int f_sizeof_hot = sizeof(HitOnTrack);
    int f_n_layers = -1;
    int f_n_events = -1;

    int f_extra_sections = 0;

    // ---- everything above is the v7/v8 header, byte for byte ----
    // New fields go BELOW, and openRead() reads the prefix first and the rest
    // only for versions that have it; a blind sizeof-sized fread would swallow
    // event data from an older file.

    // Identity of the geometry this sample was written against, copied from the
    // geometry binary's own stamp by writeMemoryFile. Empty means the file
    // predates the stamp. This is the thing whose absence let two 2024 samples
    // pass every check while sitting 0.26 cm off their own module planes.
    static constexpr int s_geom_version_size = 64;
    char f_geom_version[s_geom_version_size] = {0};

    DataFileHeader() = default;

    // Size of the v7/v8 header: 7 ints through f_n_events, plus f_extra_sections.
    static constexpr size_t s_v8_size = 8 * sizeof(int);

    // Header size ON FILE for a given format version. Magic and version are the
    // first two ints of the file precisely so this can be asked before anything
    // else is read; every future growth adds a case here and nothing else in the
    // reader has to know. Never use sizeof(DataFileHeader) for this -- that is
    // the CURRENT version's size, and using it both over-reads an older file's
    // header and seeks past its first event.
    static constexpr size_t size_of_version(int v) {
      return v >= 9 ? s_v8_size + s_geom_version_size : s_v8_size;
    }
  };

  struct DataFile {
    enum ExtraSection {
      ES_SimTrackStates = 0x1,
      ES_Seeds = 0x2,
      ES_CmsswTracks = 0x4,
      ES_HitIterMasks = 0x8,
      ES_BeamSpot = 0x10,
      ES_SimHitStates = 0x20
    };

    FILE *f_fp = 0;
    // Byte offset of the first event, i.e. the header size AS IT IS ON THIS FILE.
    // Not sizeof(DataFileHeader): that grows with every format version, and a v7
    // or v8 file has a shorter header, so positioning off sizeof() would seek
    // past the first event's size word and read nothing at all.
    long f_data_start = sizeof(DataFileHeader);
    long f_pos = sizeof(DataFileHeader);

    DataFileHeader f_header;

    std::mutex f_next_ev_mutex;

    // ----------------------------------------------------------------

    bool hasSimTrackStates() const { return f_header.f_extra_sections & ES_SimTrackStates; }
    bool hasSeeds() const { return f_header.f_extra_sections & ES_Seeds; }
    bool hasCmsswTracks() const { return f_header.f_extra_sections & ES_CmsswTracks; }
    bool hasHitIterMasks() const { return f_header.f_extra_sections & ES_HitIterMasks; }
    bool hasBeamSpot() const { return f_header.f_extra_sections & ES_BeamSpot; }
    bool hasSimHitStates() const { return f_header.f_extra_sections & ES_SimHitStates; }
    const char* geomVersion() const { return f_header.f_geom_version; }

    int openRead(const std::string &fname, int expected_n_layers,
                 const std::string &expected_geom_version = "");
    void openWrite(const std::string &fname, int n_layers, int n_ev, int extra_sections = 0,
                   const std::string &geom_version = "");

    void rewind();

    int advancePosToNextEvent(FILE *fp);

    void skipNEvents(int n_to_skip);

    void close();
    void CloseWrite(int n_written);  //override nevents in the header and close
  };

  void print(std::string pfx, int itrack, const Track &trk, const Event &ev);

  void print(std::string pfx, int itrack, const Track &trk, int hit_begin, int hit_end, const Event &ev);

  void print(std::string pfx, const TrackVec &tvec, const Event &ev);

  void print(std::string pfx, const Event::SimInfoFromHits &si);

#ifdef MKFIT_TRACE
  std::string format(const ::EVec3 &v, int width=8, int prec=3, char feg=' ');
  void print(std::string prefix, const ::EBiVec3 &s, std::string postfix="\n");
  void print(std::string pfx, const TrCandMeta &cm, const Event *ev);
  void print(std::string pfx, const TrCandStage &cs);
  void print(std::string pfx, const TrCandState &cs);
  void print(std::string pfx, const TrLayerSearch &ls);
  void print(std::string pfx, const TrHitMatch &hm);
  void print(std::string pfx, const TrKalmanUpdate &ku);
#endif

}  // end namespace mkfit
#endif
