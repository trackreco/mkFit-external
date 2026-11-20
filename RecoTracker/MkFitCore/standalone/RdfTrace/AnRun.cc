#include "AnRun.h"


#include "RecoTracker/MkFitCore/standalone/RdfTrace/RdfVectorSource.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TBrowser.h"

#include <format>
#include <cstdarg>
#include <mutex>

using RNode = ROOT::RDF::RNode;

//====================================================================================
// an_printf -- tee of stdout into <prefix>.txt
//====================================================================================

namespace {
  FILE *g_an_log = nullptr;
  // an_printf() is called from RDF Foreach bodies, which run concurrently once
  // implicit MT is on. The lock keeps each line whole in both destinations.
  std::mutex g_an_print_mutex;
}

void an_log_open(const std::string &fname) {
  an_log_close();
  g_an_log = fopen(fname.c_str(), "w");
  if (g_an_log == nullptr)
    printf("*** an_log_open: failed to open '%s' for writing.\n", fname.c_str());
}

void an_log_close() {
  if (g_an_log) {
    fclose(g_an_log);
    g_an_log = nullptr;
  }
}

int an_printf(const char *fmt, ...) {
  std::lock_guard<std::mutex> plock(g_an_print_mutex);
  va_list ap;
  va_start(ap, fmt);
  int n = vprintf(fmt, ap);
  va_end(ap);
  if (g_an_log) {
    va_start(ap, fmt);
    vfprintf(g_an_log, fmt, ap);
    va_end(ap);
    fflush(g_an_log);
  }
  return n;
}

AnRun::~AnRun() {
  for (auto ev_ptr : m_ev_vec) {
    delete ev_ptr;
  }
  an_log_close();
}

void AnRun::SetPrefix(const std::string &p) {
  m_prefix = p;
  an_log_open(m_prefix + ".txt");
  an_printf("### AnRun output prefix '%s' -- canvases to %s.root, this log to %s.txt\n\n",
            m_prefix.c_str(), m_prefix.c_str(), m_prefix.c_str());
}

//====================================================================================
#pragma region Vector Source
//====================================================================================

#define EV CTX.ev

void AnRun::RunOldVecBased() {
  m_rdf_hitmatch = mkfit::RdfVectorSources::MakeTrHitMatchDF(CTX.ev)
  .Define("C", [this](int cid) -> const TrCandState* { return & EV.trCandStates_[cid]; }, {"state_id"})
  .Define("cand_pt", [this](int cid) { return EV.trCandStates_[cid].kine.mom.R(); }, {"state_id"})
  .Define("cand_step", [](const TrCandState* C) { return C->step; }, {"C"})
  .Define("cand_layer", [](const TrCandState* C) { return C->layer; }, {"C"})
  .Define("meta_id", [](const TrCandState* C) { return C->meta_id; }, {"C"} )
  .Define("seed_gf", [this](int mid) { return EV.trSIFHforSeedByMeta_[mid].good_frac(); }, {"meta_id"})
  .Define("cand_gf", [this](int mid) { return EV.trSIFHforCandByMeta_[mid].good_frac(); }, {"meta_id"})
  // chi2 lives on TrKalmanUpdate now, not on TrHitMatch -- join back via kalman_id.
  // kalman_id is -1 for hits that never reached a Kalman update, and -1 is also
  // the low edge the chi2 histograms below already use.
  .Define("kalman_chi2", [this](int kid) {
      return kid >= 0 ? EV.trKalmanUpdates_[kid].chi2 : -1.0f; }, {"kalman_id"})
  .Define("kalman_accepted", [this](int kid) {
      return kid >= 0 && EV.trKalmanUpdates_[kid].accepted; }, {"kalman_id"})
  ;

  // Some ranges, could just go over a subset. Booked here, but read at the very
  // end of this function: the first read runs the event loop, so reading before
  // the canvas groups below are booked would cost an extra pass over the data.
  std::vector<ROOT::RDF::RResultPtr<TStatistic>> stats;
  { auto r = *m_rdf_hitmatch;
    stats = { r.Stats("cand_pt"), r.Stats("cand_step"), r.Stats("cand_layer") };
  }

  { auto &C = NewCanvasGroup(2, 2, "cands", "Candidate properties");
    auto r = *m_rdf_hitmatch;
    C.Add(r.Histo1D("cand_pt"), "s").add_pre(CGrp::logy);
    C.Add(r.Histo1D("cand_step"));
    C.Add(r.Histo1D("cand_layer"));
    C.Add(r.Histo2D({"layer_v_step", "layer vs step", 60, 0, 60, 20, 0, 20}, "cand_layer", "cand_step"), "colz");
  }
  { auto &C = NewCanvasGroup(2, 2, "residuals_chi2", "Residuals and chi2");
    auto r = *m_rdf_hitmatch;
    C.Add(r.Histo1D("residual_x"), "s").add_pre(CGrp::logy);
    C.Add(r.Histo1D("residual_y"), "s").add_pre(CGrp::logy);
    C.Add(r.Histo1D("residual_z"), "s").add_pre(CGrp::logy);
    C.Add(r.Histo1D("kalman_chi2"), "s").add_pre(CGrp::logy);
  }

  { auto &C = NewCanvasGroup(2, 2, "chi2_by_step", "Chi2 by step");
    auto r = *m_rdf_hitmatch;
    C.Add(r.Histo2D({"chi2_vs_layer", "", 60, 0,60, 101, -1, 100}, "cand_layer", "kalman_chi2" ), "s").add_pre(CGrp::logy);
    C.Add(r.Filter("cand_step == 0").Histo1D({"chi2_step0", "chi2 on step 0", 101, -1, 100}, "kalman_chi2"), "s").add_pre(CGrp::logy);
    C.Add(r.Filter("cand_step == 1").Histo1D({"chi2_step1", "chi2 on step 1", 101, -1, 100}, "kalman_chi2"), "s").add_pre(CGrp::logy);
    C.Add(r.Filter("cand_step == 2").Histo1D({"chi2_step2", "chi2 on step 2", 101, -1, 100}, "kalman_chi2"), "s").add_pre(CGrp::logy);
  }

  // ==== large y residuals ====
  { auto &C = NewCanvasGroup(2, 2, "resy.gt.400", "Y Residuals greater than 400");

    auto r = (*m_rdf_hitmatch)
    .Filter("std::abs(residual_y) > 400")
    .Define("parent_layer", [this](const TrCandState* C) { return C->parent_id >=0 ? EV.trCandStates_[C->parent_id].layer : -1; }, {"C"} );

    C.Add(r.Histo1D("cand_layer"));
    C.Add(r.Histo1D("parent_layer"));
    C.Add(r.Histo2D({"parent_layer_v_cand_layer", "parent layer vs cand layer", 61, -1, 60, 61, -1, 60}, "parent_layer", "cand_layer"));
    C.Add(r.Histo1D("residual_x"), "s").add_pre(CGrp::logy);

    f_resdy_400 = r;
  }

  // ==== going into pixels layers ====
  // works for direct backward search only (cand_step == X)... to be tuned further
  // Post-Kalman errors used to come off TrHitMatch::kalman_state, which no
  // longer exists -- they are now on TrKalmanUpdate::updated_state, reached via
  // kalman_id, and only present in MKFIT_TRACE_KALMAN_DEBUG builds.
  // TrKalmanUpdate::propagated_state is available the same way, if the
  // pre-update errors are wanted here too.
#ifdef MKFIT_TRACE_KALMAN_DEBUG
  #define KU_ERR_DEFINE(NAME)                                                  \
    .Define(#NAME, [this](int kid) {                                           \
        return kid >= 0 ? EV.trKalmanUpdates_[kid].updated_state.NAME() : -1.0f; }, \
        {"kalman_id"})
#endif

  auto define_kalman_stuff = [this](RNode r) -> RNode {
#ifdef MKFIT_TRACE_KALMAN_DEBUG
    return r
    KU_ERR_DEFINE(exx)
    KU_ERR_DEFINE(eyy)
    KU_ERR_DEFINE(ezz)
    KU_ERR_DEFINE(epT)
    KU_ERR_DEFINE(etheta)
    KU_ERR_DEFINE(emomPhi);
#else
    return r;
#endif
  };
  auto plot_kalman_stuff = [](RNode r, CanvasGroup &C) -> void {
    C.Add(r.Histo1D("residual_x"));
    C.Add(r.Histo1D("residual_y"));
    C.Add(r.Histo1D("residual_z"));
    C.Add(r.Histo1D("dphi"));
    C.Add(r.Histo1D("dq"));
    C.Add(r.Histo1D("rank"));
    C.Add(r.Histo1D("kalman_accepted"));
    C.Add(r.Histo1D("kalman_chi2"));
#ifdef MKFIT_TRACE_KALMAN_DEBUG
    C.Add(r.Histo1D("exx"));
    C.Add(r.Histo1D("eyy"));
    C.Add(r.Histo1D("ezz"));
    C.Add(r.Histo1D("epT"));
    C.Add(r.Histo1D("etheta"));
    C.Add(r.Histo1D("emomPhi"));
#endif
  };
  { auto &C = NewCanvasGroup(5,3, "SeedIntoLay3", "Seed propagated into layer 3");
    auto r = m_rdf_hitmatch->Filter("cand_step == 0 && seed_gf > 0.9 && layer == 3 && mc_match");
    plot_kalman_stuff(define_kalman_stuff(r), C);
  }
  { auto &C = NewCanvasGroup(5,3, "SeedIntoLay2", "Seed propagated into layer 2");
    auto r = m_rdf_hitmatch->Filter("cand_step == 1 && seed_gf > 0.9 && layer == 2 && mc_match");
    plot_kalman_stuff(define_kalman_stuff(r), C);
  }

  // ==== hit goodness etc via meta ====
  { auto &C = NewCanvasGroup(2,2, "good_frac_by_meta", "good_fraction by meta");
    m_rdf_meta = mkfit::RdfVectorSources::MakeTrCandMetaDF(CTX.ev)
    .Define("seed_gf", [this](int mid) { return EV.trSIFHforSeedByMeta_[mid].good_frac(); }, { "id" })
    .Define("cand_gf", [this](int mid) { return EV.trSIFHforCandByMeta_[mid].good_frac(); }, { "id" })
    .Define("cand_good_pix", [this](int mid) { return EV.trSIFHforCandByMeta_[mid].n_pix_match; }, { "id" })
    .Define("cand_bad_pix", [this](int mid) { return EV.trSIFHforCandByMeta_[mid].n_pix_bad(); }, { "id" })
    ;
    auto r = *m_rdf_meta;
    C.Add(r.Histo1D("seed_gf"));
    C.Add(r.Histo1D("cand_gf"));
    C.Add(r.Histo1D("cand_good_pix"));
    C.Add(r.Histo1D("cand_bad_pix"));
  }

  // Read the ranges booked at the top -- everything on m_rdf_hitmatch is booked
  // by now, so this triggers a single event loop for all of it.
  an_printf("cand_pt: %f -> %f, cand_step: %f -> %f, cand_layer: %f - %f\n",
        stats[0]->GetMin(), stats[0]->GetMax(),
        stats[1]->GetMin(), stats[1]->GetMax(),
        stats[2]->GetMin(), stats[2]->GetMax());
}

#undef EV

#pragma endregion
//==============================================================================
#include "AnRunEvMap.h"
//==============================================================================
#pragma region Event Source
//==============================================================================

/* Multi selection pre-proto
  struct SelectionInfo {
    std::string cut;
    std::string label;
    Color_t color;
  };

  std::vector<SelectionInfo> selections = {
    {"seed_gf > 0.9", "gf > 0.9", kRed},
    {"seed_gf > 0.7", "gf > 0.7", kBlue},
    {"seed_gf > 0.5", "gf > 0.5", kGreen},
  };
*/

//====================================================================================

void AnRun::SetupRdfEvent(std::vector<const mkfit::Event*>& ev_vec) {
  m_ev_vec.swap( ev_vec );
  m_rdf_event = mkfit::RdfSources::MakeEventDF(m_ev_vec)
  .Define("evtID", [](const mkfit::Event* ev) { return ev->evtID(); }, {"event"})
  .Define("meta_id", EV_MAP(trCandMetas_, id))
  ;
}

const mkfit::Event* AnRun::get_event_ptr(int event_id) const {
  if (event_id > 0 && event_id <= (int) m_ev_vec.size())
    return m_ev_vec.at(event_id - 1);
  return nullptr;
}

//====================================================================================

void AnRun::RunBasicSeedCandCheck() {

  auto r = (*m_rdf_event)
  // seeds, directly from trSeeds_ -- to compare to stuff from trCandMetas_
  .Define("seed_pt",  EV_MAP(trSeeds_, pT()))
  .Define("seed_eta", EV_MAP(trSeeds_, momEta()))
  .Define("seed_n_hits", EV_MAP(trSeeds_, nFoundHits()))

  .Define("seed_gf", EV_GATHER(trSIFHforSeedByMeta_, good_frac(), "meta_id"))
  .Define("cand_gf", EV_GATHER(trSIFHforCandByMeta_, good_frac(), "meta_id" ))
  .Define("cand_good_pix", EV_GATHER(trSIFHforCandByMeta_, n_pix_match, "meta_id" ))
  .Define("cand_bad_pix", EV_GATHER(trSIFHforCandByMeta_, n_pix_bad(), "meta_id" ))
  ;

  { auto &C = NewCanvasGroup(4, 3, "seeds", "seed properties");
    C.AddRealH1D(r, "seed_pt", 100, 0, 10, "s").add_pre(CGrp::logy);
    C.AddRealH1D(r, "seed_eta", 80, -4, 4, "s");
    C.AddIntH1D(r, "seed_n_hits", 6, 12, "s");
    C.Add(r.Histo2D({"seed_n_hits_vs_eta", "seed_n_hits_vs_eta;seed_eta;seed_n_hits",
                    100, -4, 4, 7, 5.5, 12.5}, "seed_eta", "seed_n_hits" ), "s colz");
    C.AddRealH1D(r, "seed_gf", 55, -0.05, 1.05, "s");
    C.AddRealH1D(r, "cand_gf", 55, -0.05, 1.05, "s");
    C.AddIntH1D(r, "cand_good_pix", 0, 10, "s");
    C.AddIntH1D(r, "cand_bad_pix", 0, 10, "s");
    C.AddRealIntH2D(r, "seed_pt", "cand_good_pix", 100, 0, 10, 0, 10, "s lego2");
    C.AddRealIntH2D(r, "seed_eta", "cand_good_pix", 60, -3, 3, 0, 10, "s lego2");
  }
}

void AnRun::RunMetaVsSeedDuplicateCheck() {

  auto r = (*m_rdf_event)
  .Define("meta_seed", EV_MAP(trCandMetas_, seed))
  .Define("n_metas", "(int) meta_id.size()")
  .Define("n_seeds", "(int) seed_pt.size()")
  .Define("n_unique_seeds",
        [](const ROOT::RVec<int>& v) -> int {
          std::set<int> s(v.begin(), v.end());
          return (int) s.size();
        },
        {"meta_seed"})
  .Define("has_duplicate_metas", "n_metas > n_unique_seeds")
  .Define("metas_per_seed_ratio", "(float)n_metas / n_seeds")

  .Define("min_seed", "Min(meta_seed)")
  .Define("max_seed", "Max(meta_seed)")
  .Define("seed_range", "max_seed - min_seed + 1")
  ;

  // === Find events with duplicates ===
  auto dup_events = r.Filter("has_duplicate_metas == true");

  // Book both counts before reading either one -- reading triggers the event
  // loop, so the inline GetValue() pair used to cost two loops.
  auto rp_n_dup = dup_events.Count();
  auto rp_n_all = r.Count();

  an_printf("Events with duplicate metas: %llu / %llu\n",
         rp_n_dup.GetValue(), rp_n_all.GetValue());

  r.Foreach([](int evtID, int n_seeds, int n_metas, int seed_range,
                          int n_unique, bool has_dup, float ratio) {
    an_printf("%-10d ns=%-10d nm=%-10d ss_range=%-10d n_uniq_ss=%-10d has_dup=%-10s ratio=%-10.2f\n",
         evtID, n_seeds, n_metas, seed_range, n_unique, has_dup ? "YES" : "NO", ratio);
    }, {"evtID", "n_seeds", "n_metas", "seed_range",
        "n_unique_seeds", "has_duplicate_metas", "metas_per_seed_ratio"});

  // === Print details for duplicate events ===
  dup_events
    .Define("duplicate_seeds",
            [](const ROOT::RVec<int>& v) {
              std::map<int, int> counts;
              for (int x : v) counts[x]++;
              ROOT::RVec<int> dups;
              for (const auto& [val, cnt] : counts) {
                if (cnt > 1) dups.push_back(val);
              }
              return dups;
            },
            {"meta_seed"})
    .Foreach([](const ROOT::RVec<int>& dup_ss,
                const ROOT::RVec<int>& meta_id,
                const ROOT::RVec<int>& meta_seed,
                int evtID) {
      an_printf("Event %d: %lu duplicate seeds\n", evtID, dup_ss.size());
      for (int ss : dup_ss) {
        an_printf("  seed=%d: meta_ids={", ss);
        for (size_t i = 0; i < meta_seed.size(); ++i) {
          if (meta_seed[i] == ss) {
            an_printf("%d ", meta_id[i]);
          }
        }
        an_printf("}\n");
      }
    }, {"duplicate_seeds", "meta_id", "meta_seed", "evtID"});

  { auto &C = NewCanvasGroup(2, 2, "meta_dup", "search for duplicate metas");
    C.Add(r.Histo1D("n_metas"));
    C.Add(r.Histo1D("n_unique_seeds"));
    C.Add(r.Histo2D({"n_metas_vs_n_seeds", "n_metas vs n_seeds", 100, 0, 300, 100, 0, 300},
                    "n_seeds", "n_metas"), "colz");
    C.Add(r.Histo1D("has_duplicate_metas"));
  }
}

//==============================================================================

void AnRun::Run_T5_vs_pT5_AsSeeds_DuplicateCount() {

  auto r = (*m_rdf_event)
  .Define("full_seed_gf", EV_GATHER(trSIFHforSeedByMeta_, good_frac(), "meta_id"))
  .Define("meta_mask", "full_seed_gf > 0.9")
  .Define("selected_metas", EV_COMPRESS(trCandMetas_, id, "meta_mask"))
  .Define("selected_seeds", EV_COMPRESS(trCandMetas_, seed, "meta_mask"))
  .Define("selected_sims", EV_COMPRESS(trCandMetas_, sim, "meta_mask"))

  .Define("seed_pt", EV_GATHER(trSeeds_, pT(), "selected_seeds"))
  .Define("seed_eta", EV_GATHER(trSeeds_, momEta(), "selected_seeds"))
  .Define("seed_n_hits", EV_GATHER(trSeeds_, nFoundHits(), "selected_seeds"))

  .Define("sim_pt", EV_GATHER(simTracks_, pT(), "selected_sims"))

  .Define("sim_n_pix_hits", EV_GATHER_FUNC(simTracks_, countAllPixelHits, "selected_sims"))
  .Define("sim_n_strip_hits", EV_GATHER_FUNC(simTracks_, countAllStripHits, "selected_sims"))

  .Define("seed_n_dups", [](const mkfit::Event* ev, const ROOT::RVec<int> meta_idcs) {
    ROOT::RVec<int> out;
    std::map<int, int> lbl_to_count;
    for (auto idx : meta_idcs) {
      int label = ev->trSIFHforSeedByMeta_[idx].label;
      ++lbl_to_count[label];
    }
    for (auto [label, count] : lbl_to_count) {
      out.push_back(count - 1);
    }
    return out;
  }, {"event", "selected_metas"})

  .Define("seed_n_dups_re_pTs", [](const mkfit::Event* ev, const ROOT::RVec<int> meta_idcs) {
    ROOT::RVec<int> out;
    std::set<int> pT_lbls;
    for (int i = 0; i < ev->seedVecInsp_.n_pTNs; i++) {
      pT_lbls.insert(ev->simInfoForTrack(ev->seedTracks_[i]).label);
    }
    std::map<int, int> lbl_to_dupstate;
    for (auto idx : meta_idcs) {
      int label = ev->trSIFHforSeedByMeta_[idx].label;
      lbl_to_dupstate[label] = (pT_lbls.find(label) == pT_lbls.end()) ? 0 : 1;
    }
    for (auto [label, dupstate] : lbl_to_dupstate) {
      out.push_back(dupstate);
    }
    return out;
  }, {"event", "selected_metas"})
  ;

  r = r
  .Define("delta_pt_o_sim_pt", "(seed_pt - sim_pt)/sim_pt");

  { auto &C = NewCanvasGroup(5, 2, "Tx_as_seed_sim_stuff", "Tx as seeds, duplicates -- with good_hit_frac > 0.9");
    C.AddRealH1D(r, "delta_pt_o_sim_pt", 100, -10, 10, "s").add_pre(CGrp::logy);
    C.AddRealH1D(r, "seed_pt", 100, 0, 30, "s").add_pre(CGrp::logy);
    C.AddRealH1D(r, "sim_pt", 100, 0, 30, "s").add_pre(CGrp::logy);
    C.Add(r.Histo1D("sim_n_pix_hits")).add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_pix_hits", 0, 12, "s");
    C.Add(r.Histo1D("sim_n_strip_hits")).add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_strip_hits", 0, 50, "s");
    C.AddIntH1D(r, "seed_n_dups", 0, 5, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "seed_n_dups_re_pTs", 0, 5, "s").add_pre(CGrp::logy);
  }

  // ==== Print events with high hit counts for manual investigation ====
  auto r_high_hits = r
  .Filter([](const ROOT::RVec<int>& pix_hits, const ROOT::RVec<int>& strip_hits) {
    for (size_t i = 0; i < pix_hits.size(); ++i) {
      if (pix_hits[i] > 30 || strip_hits[i] > 50) return true;
    }
    return false;
  }, {"sim_n_pix_hits", "sim_n_strip_hits"})
  ;

  an_printf("\n=== Events with sim_n_pix_hits > 30 OR sim_n_strip_hits > 50 ===\n");
  an_printf("%-10s %-12s %-15s %-15s\n", "evtID", "sim_idx", "pix_hits", "strip_hits");
  an_printf("%-10s %-12s %-15s %-15s\n", "-----", "-------", "--------", "----------");

  // Booked before the Foreach: Foreach is an instant action that runs the event
  // loop, so booking the count first gets it filled in that same loop instead
  // of costing a second one.
  auto rp_n_high_hits = r_high_hits.Count();

  r_high_hits.Foreach([](int evtID,
                        const ROOT::RVec<int>& sim_indices,
                        const ROOT::RVec<int>& pix_hits,
                        const ROOT::RVec<int>& strip_hits) {
    for (size_t i = 0; i < pix_hits.size(); ++i) {
      if (pix_hits[i] > 30 || strip_hits[i] > 50) {
        an_printf("%-10d %-12d %-15d %-15d\n",
              evtID, sim_indices[i], pix_hits[i], strip_hits[i]);
      }
    }
  }, {"evtID", "selected_sims", "sim_n_pix_hits", "sim_n_strip_hits"});

  an_printf("=== Total events with high-hit tracks: %llu ===\n\n", rp_n_high_hits.GetValue());
}

//==============================================================================

void AnRun::Run_T5s_into_Pix() {

  auto r_top = (*m_rdf_event)
  .Define("full_seed_gf", EV_GATHER(trSIFHforSeedByMeta_, good_frac(), "meta_id"))

  // Watch it, there are some good_frac = 0 cases, sim tracks NOT defined there!
  // .Define("pre_meta_mask", "full_seed_gf >= 1.0f")
  // .Define("pre_meta_mask", "full_seed_gf >= 0.5f && full_seed_gf < 1.0f")
  .Define("pre_meta_mask", "full_seed_gf >= 0.5f")
  .Define("pre_selected_metas", EV_COMPRESS(trCandMetas_, id, "pre_meta_mask"))
  .Define("pre_selected_seeds", EV_COMPRESS(trCandMetas_, seed, "pre_meta_mask"))
  .Define("pre_selected_sims", EV_COMPRESS(trCandMetas_, sim, "pre_meta_mask"))

  .Define("pre_sim_n_pix_hits", EV_GATHER_FUNC(simTracks_, countInnerPixelHits, "pre_selected_sims"))
  .Define("pre_sim_n_pix_all_hits", EV_GATHER_FUNC(simTracks_, countAllPixelHits, "pre_selected_sims"))
  .Define("pre_sim_n_pix_layers", EV_GATHER_FUNC(simTracks_, countInnerPixelLayers, "pre_selected_sims"))
  .Define("pre_sim_n_pix_all_layers", EV_GATHER_FUNC(simTracks_, countAllPixelLayers, "pre_selected_sims"))

  .Define("pre_sim_n_strip_hits", EV_GATHER_FUNC(simTracks_, countOuterStripHits, "pre_selected_sims"))
  .Define("pre_sim_n_strip_all_hits", EV_GATHER_FUNC(simTracks_, countAllStripHits, "pre_selected_sims"))
  .Define("pre_sim_n_strip_layers", EV_GATHER_FUNC(simTracks_, countOuterStripLayers, "pre_selected_sims"))
  .Define("pre_sim_n_strip_all_layers", EV_GATHER_FUNC(simTracks_, countAllStripLayers, "pre_selected_sims"))

  .Define("pre_sim_pT", EV_GATHER(simTracks_, pT(), "pre_selected_sims"))
  .Define("pre_sim_eta", EV_GATHER(simTracks_, momEta(), "pre_selected_sims"))
  .Define("pre_seed_first_layer", EV_GATHER(trSeeds_, getHitLyr(0), "pre_selected_seeds"))
  .Define("pre_sim_last_inner_pixel_layer", EV_GATHER_FUNC(simTracks_, lastInnerPixelLayer, "pre_selected_sims"))
  ;

  // "Final" filter for barrel, 4 sim pixel layers, first seed hit in barrel, last sim pixel in barrel
  r_top = r_top

  // Barrel T5 into outer pixel barrel with at least 4 sim pixel layers
  .Define("T5_into_pix_barrel_mask", "    pre_sim_n_pix_layers >= 4"
                                     "&& (pre_seed_first_layer == 4 || pre_seed_first_layer == 5)"
                                     "&&  pre_sim_last_inner_pixel_layer == 3")

  // More relaxed, exploration of strip barrel into pixel endcap
  .Define("T5_into_pix_ec_mask", "    pre_sim_n_pix_layers > 3"
                                 "&& (pre_seed_first_layer == 4 || pre_seed_first_layer == 5)"
                                 "&&  pre_sim_last_inner_pixel_layer > 3")

  // .Define("selection_mask", "T5_into_pix_barrel_mask")
  // .Define("selection_mask", "!T5_into_pix_barrel_mask")
  .Define("selection_mask", "T5_into_pix_ec_mask")
  ;

  r_top = r_top
  .Define("selected_sims", "pre_selected_sims[selection_mask]")
  .Define("selected_metas", "pre_selected_metas[selection_mask]")

  .Define("sim_n_pix_hits", "pre_sim_n_pix_hits[selection_mask]")
  .Define("sim_n_pix_all_hits", "pre_sim_n_pix_all_hits[selection_mask]")
  .Define("sim_n_pix_layers", "pre_sim_n_pix_layers[selection_mask]")
  .Define("sim_n_pix_all_layers", "pre_sim_n_pix_all_layers[selection_mask]")

  .Define("sim_n_strip_hits", "pre_sim_n_strip_hits[selection_mask]")
  .Define("sim_n_strip_all_hits", "pre_sim_n_strip_all_hits[selection_mask]")
  .Define("sim_n_strip_layers", "pre_sim_n_strip_layers[selection_mask]")
  .Define("sim_n_strip_all_layers", "pre_sim_n_strip_all_layers[selection_mask]")

  .Define("sim_pT", "pre_sim_pT[selection_mask]")
  .Define("sim_eta", "pre_sim_eta[selection_mask]")
  .Define("seed_first_layer", "pre_seed_first_layer[selection_mask]")
  ;

  { auto &C = NewCanvasGroup(4, 2, "t5intoPix_full_n_pixel_hits", "Tracing Tx from barrel layers 4/5 into pixels -- full pixel N_hits");
    auto r = r_top;
    C.AddIntH1D(r, "pre_sim_n_pix_hits", 0, 20, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "pre_sim_n_pix_layers", 0, 12, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "pre_sim_n_pix_all_hits", 0, 20, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "pre_sim_n_pix_all_layers", 0, 12, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_pix_hits", 0, 20, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_pix_layers", 0, 12, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_pix_all_hits", 0, 20, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_pix_all_layers", 0, 12, "s").add_pre(CGrp::logy);
  }
  { auto &C = NewCanvasGroup(4, 2, "t5intoPix_full_n_strip_hits", "Tracing Tx from barrel layers 4/5 into pixels -- full strip N_hits");
    auto r = r_top;
    C.AddIntH1D(r, "pre_sim_n_strip_hits", 0, 50, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "pre_sim_n_strip_layers", 0, 24, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "pre_sim_n_strip_all_hits", 0, 50, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "pre_sim_n_strip_all_layers", 0, 24, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_strip_hits", 0, 50, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_strip_layers", 0, 24, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_strip_all_hits", 0, 50, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_strip_all_layers", 0, 24, "s").add_pre(CGrp::logy);
  }

  { auto &C = NewCanvasGroup(4, 2, "t5intoPix", "Tracing Tx from barrel layers 4/5 into pixels -- selections");
    auto r = r_top;
    C.AddIntH1D(r, "sim_n_pix_hits", 0, 20, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_pix_layers", 0, 12, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_strip_hits", 0, 50, "s").add_pre(CGrp::logy);
    C.AddIntH1D(r, "sim_n_strip_layers", 0, 24, "s").add_pre(CGrp::logy);

    C.Add(r.Histo1D("sim_pT"), "s").add_pre(CGrp::logy);
    C.Add(r.Histo1D("sim_eta"), "s");
    C.Add(r.Histo1D("pre_sim_pT"), "s").add_pre(CGrp::logy);
    C.Add(r.Histo1D("pre_sim_eta"), "s");
  }

  // ===== HitMatch plots
  // By specific layers (could be a set)
  // for (auto layer : { 2, 3 })
  // {
  //   // First hit-match pass, define hm_indices and prepare index columns for sub-selections.
  //   auto r_layer = define_hitmatch_indices_by_layers(r_top, { layer });
  //   define_and_plot_hitmatch_selections(r_layer, std::format("L{}", layer), std::format("Layer {}", layer) );
  // }
  // By pixel layer sequence on the sim track.
  for (auto layer_offset : { 0, 1, 2 })
  {
    // First hit-match pass, define hm_indices and prepare index columns for sub-selections.
    auto r_layer = define_hitmatch_indices_by_sim_pixel_layer(r_top, layer_offset, false);
    define_and_plot_hitmatch_selections(r_layer, std::format("PTL{}", layer_offset), std::format("Pixel Tail Layer {}", layer_offset) );
  }

  // ======= Candidate Track Pixel Hit Analysis =======
  {
    // Gather candidate track info for selected metas
    auto r = r_top
    .Define("cand_n_pix", EV_GATHER(trSIFHforCandByMeta_, n_pix, "selected_metas"))
    .Define("cand_n_pix_match", EV_GATHER(trSIFHforCandByMeta_, n_pix_match, "selected_metas"))
    .Define("cand_n_pix_bad", EV_GATHER(trSIFHforCandByMeta_, n_pix_bad(), "selected_metas"))
    .Define("cand_n_strip", EV_GATHER(trSIFHforCandByMeta_, n_strip, "selected_metas"))
    .Define("cand_n_strip_match", EV_GATHER(trSIFHforCandByMeta_, n_strip_match, "selected_metas"))
    .Define("cand_n_strip_bad", EV_GATHER(trSIFHforCandByMeta_, n_strip_bad(), "selected_metas"))
    .Define("cand_good_frac", EV_GATHER(trSIFHforCandByMeta_, good_frac(), "selected_metas"))
    ;

    { auto &C = NewCanvasGroup(3, 3, "cand_pix_hits", "Candidate Track Pixel Hits (Selected Metas)");
      // Pixel hit distributions
      C.AddIntH1D(r, "cand_n_pix", 0, 12, "s").add_pre(CGrp::logy);
      C.AddIntH1D(r, "cand_n_pix_match", 0, 12, "s").add_pre(CGrp::logy);
      C.AddIntH1D(r, "cand_n_pix_bad", 0, 12, "s").add_pre(CGrp::logy);

      // Strip hit distributions
      C.AddIntH1D(r, "cand_n_strip", 0, 30, "s").add_pre(CGrp::logy);
      C.AddIntH1D(r, "cand_n_strip_match", 0, 30, "s").add_pre(CGrp::logy);
      C.AddIntH1D(r, "cand_n_strip_bad", 0, 30, "s").add_pre(CGrp::logy);

      // Good fraction (matching / total)
      C.Add(r.Histo1D({"cand_good_frac", "Candidate Good Hit Fraction;Fraction;Events", 55, -0.05, 1.05}, "cand_good_frac"), "s");
    }

    // ======= 2D Correlations: Pixel Hits vs Selection Variables =======
    { auto &C = NewCanvasGroup(3, 2, "cand_pix_corr", "Candidate Pixel Hits vs Sim/Seed Properties");

      // Pixel hits vs sim pT
      C.Add(r.Histo2D({"cand_n_pix_vs_sim_pT", "Pixel Hits vs Sim pT;Sim pT (GeV);N Pixel Hits",
                            60, 0, 30, 15, -0.5, 14.5}, "sim_pT", "cand_n_pix"), "colz");
      // Pixel hits vs sim eta
      C.Add(r.Histo2D({"cand_n_pix_vs_sim_eta", "Pixel Hits vs Sim #eta;Sim #eta;N Pixel Hits",
                            80, -4, 4, 15, -0.5, 14.5}, "sim_eta", "cand_n_pix"), "colz");
      // Matching Pixel hits vs sim eta
      C.Add(r.Histo2D({"cand_n_pix_matched_vs_sim_eta", "Matched Pixel Hits vs Sim #eta;Sim #eta;N Pixel Hits",
                            80, -4, 4, 15, -0.5, 14.5}, "sim_eta", "cand_n_pix_match"), "colz");

      // Pixel match fraction vs sim eta
      C.Add(r.Histo2D({"cand_good_frac_vs_sim_eta", "Match Fraction vs Sim #eta;Sim #eta;Good Fraction",
                            80, -4, 4, 22, -0.1, 1.1}, "sim_eta", "cand_good_frac"), "colz");

      // Bad pixel hits vs seed first layer
      C.Add(r.Histo2D({"cand_n_pix_bad_vs_seed_layer", "Bad Pixels vs Seed First Layer;Seed First Layer;N Bad Pixels",
                            10, 0.5, 10.5, 15, -0.5, 14.5}, "seed_first_layer", "cand_n_pix_bad"), "colz");

      // Good fraction vs sim pixel layers
      C.Add(r.Histo2D({"cand_good_frac_vs_sim_pix_layers", "Match Fraction vs Sim Pixel Layers;Sim Pixel Layers;Good Fraction",
                            12, -0.5, 11.5, 22, -0.1, 1.1}, "sim_n_pix_layers", "cand_good_frac"), "colz");
    }

    // Assign to member for manual processing

    m_T5 = r;

    // ======= Print Summary Statistics =======
    an_printf("\n=== Candidate Track Pixel Hit Summary (Selected Metas) ===\n");

    // Book all three, then read: a GetValue() in between bookings triggers the
    // event loop, so the inline form costs one full loop per Stats.
    auto rp_n_pix = r.Stats("cand_n_pix");
    auto rp_n_pix_match = r.Stats("cand_n_pix_match");
    auto rp_n_pix_bad = r.Stats("cand_n_pix_bad");

    auto stats_n_pix = rp_n_pix.GetValue();
    auto stats_n_pix_match = rp_n_pix_match.GetValue();
    auto stats_n_pix_bad = rp_n_pix_bad.GetValue();

    an_printf("cand_n_pix:       Mean=%.2f, RMS=%.2f, Min=%.0f, Max=%.0f\n",
           stats_n_pix.GetMean(), stats_n_pix.GetRMS(), stats_n_pix.GetMin(), stats_n_pix.GetMax());
    an_printf("cand_n_pix_match: Mean=%.2f, RMS=%.2f, Min=%.0f, Max=%.0f\n",
           stats_n_pix_match.GetMean(), stats_n_pix_match.GetRMS(), stats_n_pix_match.GetMin(), stats_n_pix_match.GetMax());
    an_printf("cand_n_pix_bad:   Mean=%.2f, RMS=%.2f, Min=%.0f, Max=%.0f\n",
           stats_n_pix_bad.GetMean(), stats_n_pix_bad.GetRMS(), stats_n_pix_bad.GetMin(), stats_n_pix_bad.GetMax());

    an_printf("==========================================================\n\n");
  }
}



// Guarded gather of a per-state quantity out of trCandStates_, indexed by a
// state-id column that can contain -1 (no such stage / state). Note that
// gather_with_member() / EV_GATHER do NOT bounds-check, so they must not be
// used with stage2_root_state_id.
template <typename F, typename T>
static auto gather_state_guarded(const mkfit::Event* ev, const ROOT::RVec<int>& sid, F f, T sentinel) {
  ROOT::RVec<T> out(sid.size());
  for (size_t i = 0; i < sid.size(); ++i) {
    out[i] = (sid[i] >= 0 && sid[i] < (int) ev->trCandStates_.size())
             ? f(ev->trCandStates_[sid[i]]) : sentinel;
  }
  return out;
}

void AnRun::Run_Stage2_RootState_QualityCheck() {
  // Stage 2 = BkwSearch (index 2 in stage_ids array)
  auto r = (*m_rdf_event)
  .Define("stage2_idx",
    [](const mkfit::Event* ev) {
      return map_with_member(ev->trCandMetas_, [](auto& m){ return m.stage_ids[2]; });
    }, {"event"})
  .Define("stage2_root_state_id",
    [](const mkfit::Event* ev, const ROOT::RVec<int>& stage2_idx) {
      ROOT::RVec<int> root_ids(stage2_idx.size());
      for (size_t i = 0; i < stage2_idx.size(); ++i) {
        root_ids[i] = (stage2_idx[i] >= 0 && stage2_idx[i] < (int)ev->trCandStages_.size())
                      ? ev->trCandStages_[stage2_idx[i]].root_state_id
                      : -1;
      }
      return root_ids;
    }, {"event", "stage2_idx"})
  .Define("stage2_final_state_id",
    [](const mkfit::Event* ev, const ROOT::RVec<int>& stage2_idx) {
      ROOT::RVec<int> final_ids(stage2_idx.size());
      for (size_t i = 0; i < stage2_idx.size(); ++i) {
        final_ids[i] = (stage2_idx[i] >= 0 && stage2_idx[i] < (int)ev->trCandStages_.size())
                       ? ev->trCandStages_[stage2_idx[i]].final_state_id
                       : -1;
      }
      return final_ids;
    }, {"event", "stage2_idx"})
  // trSIFHforSeedByMeta_ is meta-indexed, so it gathers by "meta_id".
  // trSeeds_ is seed-indexed, so it must gather by the seed index taken off the meta.seed
  .Define("seed_id", EV_GATHER(trCandMetas_, seed, "meta_id"))
  .Define("seed_gf", EV_GATHER(trSIFHforSeedByMeta_, good_frac(), "meta_id"))
  .Define("seed_eta", EV_GATHER(trSeeds_, momEta(), "seed_id"))
  .Define("seed_pt", EV_GATHER(trSeeds_, pT(), "seed_id"))
  .Define("has_nan",
    [](const mkfit::Event* ev, const ROOT::RVec<int>& root_state_id) {
      ROOT::RVec<int> nans(root_state_id.size());
      for (size_t i = 0; i < root_state_id.size(); ++i) {
        nans[i] = (root_state_id[i] >= 0 && root_state_id[i] < (int)ev->trCandStates_.size() &&
                   ev->trCandStates_[root_state_id[i]].state.hasNanNSillyValues()) ? 1 : 0;
      }
      return nans;
    }, {"event", "stage2_root_state_id"})

  // Where is the root state -- layer and position. layer is an int and always
  // finite, so it stays meaningful even if the state parameters are NaN too.
  .Define("s2_layer", [](const mkfit::Event* ev, const ROOT::RVec<int>& sid) {
      return gather_state_guarded(ev, sid, [](const TrCandState& s){ return s.layer; }, -1);
    }, {"event", "stage2_root_state_id"})
  .Define("s2_r", [](const mkfit::Event* ev, const ROOT::RVec<int>& sid) {
      return gather_state_guarded(ev, sid, [](const TrCandState& s){ return s.state.posR(); }, -999.0f);
    }, {"event", "stage2_root_state_id"})
  .Define("s2_z", [](const mkfit::Event* ev, const ROOT::RVec<int>& sid) {
      return gather_state_guarded(ev, sid, [](const TrCandState& s){ return s.state.z(); }, -999.0f);
    }, {"event", "stage2_root_state_id"})

  .Define("valid_s2", "stage2_root_state_id >= 0")

  // Create masks for selections (RVec<bool>)
  .Define("nan_mask", "valid_s2 && has_nan == 1")
  .Define("ok_mask", "valid_s2 && has_nan == 0")

  // Extract masked RVecs for histogramming
  .Define("seed_gf_nan", "seed_gf[nan_mask]")
  .Define("seed_gf_ok", "seed_gf[ok_mask]")
  .Define("seed_eta_nan", "seed_eta[nan_mask]")
  .Define("seed_eta_ok", "seed_eta[ok_mask]")
  .Define("seed_pt_nan", "seed_pt[nan_mask]")
  .Define("seed_pt_ok", "seed_pt[ok_mask]")
  // Valid-only, for the "where did the probe land" plots (drops the -1/-999 sentinels).
  .Define("s2_layer_v", "s2_layer[valid_s2]")
  .Define("s2_r_v", "s2_r[valid_s2]")
  .Define("s2_z_v", "s2_z[valid_s2]")
  // 1 when the stage-2 root state is also its final state, i.e. the stage was
  // created but the candidate never advanced past it. Valid entries only.
  .Define("s2_root_eq_final",
    // NOTE: RVec comparison operators yield RVec<int>, not RVec<bool>.
    [](const ROOT::RVec<int>& root_id, const ROOT::RVec<int>& final_id, const ROOT::RVec<int>& valid) {
      ROOT::RVec<int> same;
      same.reserve(root_id.size());
      for (size_t i = 0; i < root_id.size(); ++i) {
        if (valid[i]) same.push_back(root_id[i] == final_id[i] ? 1 : 0);
      }
      return same;
    }, {"stage2_root_state_id", "stage2_final_state_id", "valid_s2"})
  ;

  m_Quality = r;

  auto &C = NewCanvasGroup(4, 3, "stage2_root_state_quality", "Stage 2 Root State Quality");

  // Count candidates (not events!)
  // Keep the result handle: the summary below reads its bins, so we neither
  // book a second histogram nor trigger a second event loop.
  auto h_nan = r.Histo1D({"has_nan", "NaN Flag;Status;Candidates", 3, -1.5, 1.5}, "has_nan");
  C.Add(h_nan);
  C.Add(r.Histo1D({"seed_gf_all", "Seed GF (All);GF;Candidates", 55, -0.05, 1.05}, "seed_gf"));
  C.Add(r.Histo1D({"seed_gf_nan", "Seed GF (NaN);GF;Candidates", 55, -0.05, 1.05}, "seed_gf_nan"));
  C.Add(r.Histo1D({"seed_eta_all", "Seed #eta (All);#eta;Candidates", 80, -4, 4}, "seed_eta"));
  C.Add(r.Histo1D({"seed_eta_nan", "Seed #eta (NaN);#eta;Candidates", 80, -4, 4}, "seed_eta_nan"));
  C.Add(r.Histo2D({"nan_vs_seedgf", "NaN vs Seed GF;Seed GF;NaN", 55, -0.05, 1.05, 3, -1.5, 1.5}, "seed_gf", "has_nan"), "colz");
  C.Add(r.Histo2D({"nan_vs_seedeta", "NaN vs Seed #eta;Seed #eta;NaN", 80, -4, 4, 3, -1.5, 1.5}, "seed_eta", "has_nan"), "colz");
  C.Add(r.Histo2D({"gf_vs_eta", "Seed GF vs #eta;#eta;GF", 80, -4, 4, 55, -0.05, 1.05}, "seed_eta", "seed_gf"), "colz");

  // Where the probe landed -- layer and position of the stage-2 root state.
  C.Add(r.Histo1D({"s2_layer", "Root state layer;layer;Candidates", 60, -0.5, 59.5}, "s2_layer_v"));
  C.Add(r.Histo1D({"s2_r", "Root state r;r (cm);Candidates", 120, 0, 120}, "s2_r_v"));
  C.Add(r.Histo1D({"s2_z", "Root state z;z (cm);Candidates", 140, -280, 280}, "s2_z_v"));
  C.Add(r.Histo2D({"s2_rz", "Root state r-z;z (cm);r (cm)", 140, -280, 280, 60, 0, 120}, "s2_z_v", "s2_r_v"), "colz");

  // Root state == final state means the stage was set up but never extended.
  auto h_same = r.Histo1D({"s2_root_eq_final", "Root state == final state;same;Candidates", 2, -0.5, 1.5},
                          "s2_root_eq_final");

  // Summary stats.
  // NOTE: bins must be read off an explicitly-modelled histogram. Booking
  // r.Histo1D("has_nan") with no model gives default binning, so the old
  // Integral(1,1)/Integral(2,2) here always reported 0.
  an_printf("\n=== Stage 2 Root State NaN Check ===\n");
  double ok_count = h_nan->GetBinContent(2);   // bin 2 == value 0
  double nan_count = h_nan->GetBinContent(3);  // bin 3 == value 1
  double total = nan_count + ok_count;
  an_printf("Valid S2 candidates: %.0f, NaN: %.0f (%.2f%%), OK: %.0f (%.2f%%)\n",
         total, nan_count, total > 0 ? 100.0*nan_count/total : 0.0,
         ok_count, total > 0 ? 100*ok_count/total : 0.0);
  an_printf("=====================================\n\n");

  // Root vs final state: how many candidates never got extended past the stage root.
  an_printf("=== Stage 2 Root vs Final State ===\n");
  double n_same = h_same->GetBinContent(2);  // bin 2 == value 1 (root == final)
  double n_ext = h_same->GetBinContent(1);   // bin 1 == value 0 (root != final)
  double n_tot = n_same + n_ext;
  an_printf("Valid S2 candidates: %.0f\n", n_tot);
  an_printf("  root == final (never extended): %.0f (%.2f%%)\n",
            n_same, n_tot > 0 ? 100.0*n_same/n_tot : 0.0);
  an_printf("  root != final (extended):       %.0f (%.2f%%)\n",
            n_ext, n_tot > 0 ? 100.0*n_ext/n_tot : 0.0);
  an_printf("===================================\n\n");

  // Print sample NaN cases
  an_printf("\n=== Sample NaN Cases (first 20) ===\n");
  std::atomic<int> nan_counter{0};
  r.Foreach([&](int evtID, const ROOT::RVec<int>& meta_id, const ROOT::RVec<int>& has_nan,
                const ROOT::RVec<float>& seed_gf, const ROOT::RVec<float>& seed_eta) {
    for (size_t i = 0; i < meta_id.size(); ++i) {
      if (nan_counter >= 20) return;
      if (has_nan[i] == 1) {
        an_printf("Event %d: meta=%d seed_gf=%.3f seed_eta=%.3f\n",
               evtID, meta_id[i], seed_gf[i], seed_eta[i]);
        ++nan_counter;
      }
    }
  }, {"evtID", "meta_id", "has_nan", "seed_gf", "seed_eta"});
  an_printf("===============================\n\n");
}

inline ROOT::RVec<float> log10_safe(const ROOT::RVec<float>& cov) {
  ROOT::RVec<float> logcov(cov.size());
  for (size_t i = 0; i < cov.size(); ++i) {
    logcov[i] = (cov[i] > 0) ? std::log10(cov[i]) : -10.0f;
  }
  return logcov;
}

#define DEF_LOG10(name) \
  .Define("log10_" #name, \
    [](const ROOT::RVec<float>& v) -> ROOT::RVec<float> { \
      ROOT::RVec<float> r(v.size()); \
      for(size_t i=0; i<v.size(); ++i) r[i] = v[i]>0 ? std::log10(v[i]) : -10.0f; \
      return r; \
    }, {#name})

void AnRun::Run_Stage2_RootState_Covariance_Check() {
  // Stage 2 = BkwSearch (index 2 in stage_ids array)
  // Extract and plot diagonal elements of covariance matrix
  auto r = (*m_rdf_event)
  .Define("seed_id", EV_GATHER(trCandMetas_, seed, "meta_id"))
  .Define("sim_id", EV_GATHER(trCandMetas_, sim, "meta_id"))
  .Define("stage2_idx", EV_GATHER(trCandMetas_, stage_ids[2], "meta_id"))
  .Define("stage2_root_state_id",
    [](const mkfit::Event* ev, const ROOT::RVec<int>& stage2_idx) {
      ROOT::RVec<int> root_ids(stage2_idx.size());
      for (size_t i = 0; i < stage2_idx.size(); ++i) {
        root_ids[i] = ev->trCandStages_[stage2_idx[i]].root_state_id;
      }
      return root_ids;
    }, {"event", "stage2_idx"})
   .Define("seed_gf", EV_GATHER(trSIFHforSeedByMeta_, good_frac(), "meta_id"))
   .Define("seed_eta", EV_GATHER(trSeeds_, momEta(), "seed_id"))
  .Define("seed_pt", EV_GATHER(trSeeds_, pT(), "seed_id"))
  .Define("valid_s2", "stage2_root_state_id >= 0")
  // Extract 6 diagonal covariance elements
  .Define("cov_00", EV_GATHER(trCandStates_, state.errors(0, 0), "stage2_root_state_id"))
  .Define("cov_11", EV_GATHER(trCandStates_, state.errors(1, 1), "stage2_root_state_id"))
  .Define("cov_22", EV_GATHER(trCandStates_, state.errors(2, 2), "stage2_root_state_id"))
  .Define("cov_33", EV_GATHER(trCandStates_, state.errors(3, 3), "stage2_root_state_id"))
  .Define("cov_44", EV_GATHER(trCandStates_, state.errors(4, 4), "stage2_root_state_id"))
  .Define("cov_55", EV_GATHER(trCandStates_, state.errors(5, 5), "stage2_root_state_id"))

  // Log10 versions (using Transform - cleaner)
  // .Define("log10_cov_00", "log10_safe(cov_00)")
  // .Define("log10_cov_11", "log10_safe(cov_11)")
  // .Define("log10_cov_22", "log10_safe(cov_22)")
  // .Define("log10_cov_33", "log10_safe(cov_33)")
  // .Define("log10_cov_44", "log10_safe(cov_44)")
  // .Define("log10_cov_55", "log10_safe(cov_55)")
  DEF_LOG10(cov_00)
  DEF_LOG10(cov_11)
  DEF_LOG10(cov_22)
  DEF_LOG10(cov_33)
  DEF_LOG10(cov_44)
  DEF_LOG10(cov_55)
  ;

  m_Cov = r;

  // Create 6 canvases, one for each diagonal element
  auto &C00 = NewCanvasGroup(2, 2, "cov_00_xx", "Cov[0][0] - X variance");
  auto &C11 = NewCanvasGroup(2, 2, "cov_11_yy", "Cov[1][1] - Y variance");
  auto &C22 = NewCanvasGroup(2, 2, "cov_22_zz", "Cov[2][2] - Z variance");
  auto &C33 = NewCanvasGroup(2, 2, "cov_33_pxpx", "Cov[3][3] - 1/pT variance");
  auto &C44 = NewCanvasGroup(2, 2, "cov_44_pypy", "Cov[4][4] - phi variance");
  auto &C55 = NewCanvasGroup(2, 2, "cov_55_pzpz", "Cov[5][5] - theta variance");

  // Common histogram definitions for each covariance element
  const char* cov_names[6] = {"log10_cov_00", "log10_cov_11", "log10_cov_22", "log10_cov_33", "log10_cov_44", "log10_cov_55"};
  const char* cov_raw_names[6] = {"cov_00", "cov_11", "cov_22", "cov_33", "cov_44", "cov_55"};
  const char* cov_titles[6] = {"log10 Cx", "log10 Cy", "log10 Cz", "log10 Cinv_pt", "log10 Cphi", "log10 Ctheta"};
  const char* cov_units[6] = {"cm^{2}", "cm^{2}", "cm^{2}", "(GeV/c)^{-2}", "(rad)^{2}", "(rad)^{2}"};

  CanvasGroup* canvases[6] = {&C00, &C11, &C22, &C33, &C44, &C55};

  for (int i = 0; i < 6; ++i) {
    auto& C = *canvases[i];

    // Distribution of covariance element
    C.Add(r.Histo1D({Form("h_%s_all", cov_names[i]), Form("%s (All);%s;Candidates", cov_titles[i], cov_units[i]), 110, -11, 0}, cov_names[i]));

    // Vs seed quality
    C.Add(r.Histo2D({Form("h_%s_vs_gf", cov_names[i]), Form("%s vs Seed GF;GF;%s", cov_titles[i], cov_units[i]), 55, -0.05, 1.05, 110, -11, 0}, "seed_gf", cov_names[i]), "colz");

    // Vs seed eta
    C.Add(r.Histo2D({Form("h_%s_vs_eta", cov_names[i]), Form("%s vs Seed #eta;#eta;%s", cov_titles[i], cov_units[i]), 80, -4, 4, 110, -11, 0}, "seed_eta", cov_names[i]), "colz");

    // Vs seed pT
    C.Add(r.Histo2D({Form("h_%s_vs_pt", cov_names[i]), Form("%s vs Seed p_{T};p_{T} (GeV);%s", cov_titles[i], cov_units[i]), 100, 0, 20, 110, -11, 0}, "seed_pt", cov_names[i]), "colz");
  }

  // Summary statistics
  an_printf("\n=== Stage 2 Root State Covariance Diagonal Check ===\n");

  // Explicit binning: a model-less Histo1D() gets default binning, the same
  // trap that made the NaN summary report zeros.
  auto h_valid = r.Histo1D({"valid_s2_count", "valid_s2;valid;Candidates", 2, -0.5, 1.5}, "valid_s2");

  // Report both scales, because they answer different questions and used to be
  // silently mixed here:
  //  - log10 mean/rms  -- the typical order of magnitude; the right summary for
  //    a quantity spread over decades. This is what the old mean/rms printed,
  //    even though the line was labelled "Cov[i][i]".
  //  - raw mean/min/max -- the actual values, in the units of the covariance.
  //    min/max were previously taken from TH1::GetMinimum()/GetMaximum(), which
  //    return bin *contents*, not the value range, hence the nonsense numbers.
  // Book everything first, then read: calling GetValue() inside the loop would
  // trigger a separate event loop per element.
  std::vector<ROOT::RDF::RResultPtr<TStatistic>> st_raw, st_log;
  for (int i = 0; i < 6; ++i) {
    st_raw.push_back(r.Stats(cov_raw_names[i]));
    st_log.push_back(r.Stats(cov_names[i]));
  }

  // First read -- this is what actually runs the loop, for h_valid and all the
  // Stats and canvas histograms booked above.
  an_printf("Valid S2 candidates: %.0f\n", h_valid->Integral());

  for (int i = 0; i < 6; ++i) {
    an_printf("Cov[%d][%d]: log10 mean=%8.4f rms=%7.4f | raw mean=%.4e min=%.4e max=%.4e\n",
              i, i, st_log[i]->GetMean(), st_log[i]->GetRMS(),
              st_raw[i]->GetMean(), st_raw[i]->GetMin(), st_raw[i]->GetMax());
  }
  an_printf("====================================================\n\n");

  // Print sample extreme cases (largest covariance values)
  // an_printf("\n=== Sample Large Covariance Cases (first 10 per element) ===\n");
  // for (int i = 0; i < 6; ++i) {
  //   an_printf("\n--- Cov[%d][%d] ---\n", i, i);
  //   std::atomic<int> counter{0};
  //   r.Foreach([&](int evtID, const ROOT::RVec<int>& meta_id,
  //                 const ROOT::RVec<float>& seed_pt, const ROOT::RVec<ROOT::RVec<float>>& cov_diag) {
  //     for (size_t j = 0; j < meta_id.size(); ++j) {
  //       if (counter >= 10) return;
  //       if (cov_diag[j][i] > 1.0) {  // Threshold for "large"
  //         an_printf("Event %d: meta=%d pt=%.2f cov[%d][%d]=%.4e\n",
  //                evtID, meta_id[j], seed_pt[j], i, i, cov_diag[j][i]);
  //         ++counter;
  //       }
  //     }
  //   }, {"evtID", "meta_id", "seed_pt", "cov_diag"});
  // }
  // an_printf("====================================================\n\n");
}
#pragma endregion
//==============================================================================
#pragma region HitMatch & Kalman
//==============================================================================

RNode AnRun::define_hitmatch_indices_by_layers(RNode r, const std::unordered_set<int> &layers) {
  // Extract hit-match indices by "selected_metas" and a set of layers -------
  // Output: new column "hm_indices"
  return r
  .Define("hm_indices", [layers](const mkfit::Event* ev, const RVecI& sel_metas) {
    std::unordered_set<int> sel_set(sel_metas.begin(), sel_metas.end());
    RVecI hm_indices;
    hm_indices.reserve(sel_metas.size() * layers.size() * 20);  // rough heuristic
    for (size_t i = 0; i < ev->trHitMatches_.size(); ++i) {
      const auto& hm = ev->trHitMatches_[i];
      if (layers.count(hm.layer) && sel_set.count( ev->trCandStates_[hm.state_id].meta_id ))
        hm_indices.push_back(i);
    }
    return hm_indices;
  }, {"event", "selected_metas"});
}

RNode AnRun::define_hitmatch_indices_by_sim_pixel_layer(RNode r, int layer_offset, bool use_head) {
  // Build a vector from meta_id → target_pixel_layer once
  // Then iterate through hit-matches with O(1) lookup
  return r.Define("hm_indices", [layer_offset, use_head](const mkfit::Event *ev,
                                                         const RVecI& sel_metas) {
    // Step 1: Build meta_id → target_layer vector (~2000 max, very fast)
    std::vector<int> meta_to_layer(ev->trCandMetas_.size(), -1);

    for (int meta_id : sel_metas) {
      const TrCandMeta &meta = ev->trCandMetas_[meta_id];
      const int sim_idx = meta.sim;

      if (sim_idx < 0 || sim_idx >= (int)ev->simTracks_.size()) continue;

      const mkfit::Track &sim_track = ev->simTracks_[sim_idx];
      int target_layer = use_head ? ev->firstInnerPixelLayer(sim_track, layer_offset)
                                  : ev->lastInnerPixelLayer(sim_track, layer_offset);

      meta_to_layer[meta_id] = target_layer;  // -1 if no valid layer found
    }

    // Step 2: Iterate through hit-matches with O(1) lookup
    RVecI hm_indices;
    hm_indices.reserve(sel_metas.size() * 6);

    for (int i = 0; i < (int) ev->trHitMatches_.size(); ++i) {
      const TrHitMatch &hm = ev->trHitMatches_[i];
      const int meta_id = ev->trCandStates_[hm.state_id].meta_id;
      assert (meta_id >= 0 && meta_id < (int) meta_to_layer.size());

      int target_layer = meta_to_layer[meta_id];
      if (target_layer >= 0 && hm.layer == target_layer) {
        hm_indices.push_back(i);
      }
    }

    return hm_indices;
  }, {"event", "selected_metas"});
}

RNode AnRun::define_hitmatch_stuff(RNode r, const std::string &idx_column, const std::string &pref) {
  // Define hit match stuff
  return r
  .Define(pref + "mc_match", EV_GATHER(trHitMatches_, mc_match, idx_column))
  .Define(pref + "passed_preselect", EV_GATHER(trHitMatches_, passed_preselect, idx_column))
  .Define(pref + "passed_pqueue", EV_GATHER(trHitMatches_, passed_pqueue, idx_column))
  .Define(pref + "dphi", EV_GATHER(trHitMatches_, dphi, idx_column))
  .Define(pref + "dq", EV_GATHER(trHitMatches_, dq, idx_column))
  .Define(pref + "rank", EV_GATHER(trHitMatches_, sub_rank, idx_column))
  ;
}

RNode AnRun::define_kalmanupdate_stuff(RNode r, const std::string &idx_column, const std::string &pref) {
  // Define Kalman update stuff
  return r
  .Define(pref + "ku_chi2", EV_GATHER(trKalmanUpdates_, chi2, idx_column))
  .Define(pref + "ku_accepted", EV_GATHER(trKalmanUpdates_, accepted, idx_column))
  // .Define(pref + "ku_chi2_trk", EV_GATHER(trKalmanUpdates_, chi2_trk, idx_column))
  ;
}

void AnRun::plot_hitmatch_layer_stuff(RNode r, CanvasGroup &C, const std::string &pref) {
  // Plot stuff relevant to hit-matching
  C.AddIntH1D(r, pref + "mc_match", 0, 1, "s");
  C.AddIntH1D(r, pref + "passed_preselect", 0, 1, "s");
  C.AddIntH1D(r, pref + "passed_pqueue", 0, 1, "s");
  C.AddIntH1D(r, pref + "rank", -1, 8, "s").add_pre(CGrp::logy);
  C.Add(r.Histo1D(pref + "dphi"), "s").add_pre(CGrp::logy);
  C.Add(r.Histo1D(pref + "dq"), "s").add_pre(CGrp::logy);
}

void AnRun::plot_kalmanupdate_layer_stuff(RNode r, CanvasGroup &C, const std::string &pref) {
  C.Add(r.Histo1D(pref + "ku_chi2"), "s").add_pre(CGrp::logy);
  C.Add(r.Histo1D({(pref + "ku_chi2_zoom_100").c_str(), "chi2 zoom to 100", 101, -1, 100},
                  pref + "ku_chi2"), "s").add_pre(CGrp::logy);
  C.AddIntH1D(r, pref + "ku_accepted", 0, 1, "s");
  // C.Add(r.Histo1D(pref + "ku_chi2_trk"), "s").add_pre(CGrp::logy);
}

namespace {

  struct HmIdxSelection {
    std::string sel_prefix; // Prefix of index sub-selection, e.g., <sel_prefix>_hm_indices. For primaries same as col_prefix.
    std::string col_prefix; // Column prefix for variables exported to top RDF (only local if empty)
    std::string title;      // Title postfix for canvas
    std::string index_cut;  // Cut applied to index column, e.g., hm_indices. Empty for primaries.

    // Primary selection is expected to be defined by the selection index base, without sel_prefix.
    // Also, it is expected that it has non-zero col_prefix and variables there already pre-populated,
    // as it is assumed they must be present for secondary selection filters.

    bool is_primary()        const { return index_cut.empty(); } // if true, columns are in the master rdf, prefixed with name
    bool has_column_prefix() const { return !col_prefix.empty(); }

    std::string selection_name(const std::string &idx_col_base) const {
      // Secondary selections are derived from primary ones done "by hand".
      return is_primary() ? idx_col_base : (sel_prefix + "_" + idx_col_base);
    }
    std::string column_prefix() const {
      // Primary columns are defined in the master rdf with name prefix, others in derived rdf without a prefix.
      return col_prefix.empty() ? "" : (col_prefix + "_");
    }
  };

  const std::vector<HmIdxSelection> kSelections = {
    { "all",                  "all",  "ALL",                              "" },
    { "match",                "",     "MC-Match",                         "all_mc_match" },
    { "match_pass_preselect", "",     "MC-match && presel && NOT pqueue", "all_mc_match && all_passed_preselect && ! all_passed_pqueue"},
    { "match_pass_all",       "good", "MC-match && presel && pqueue",     "all_mc_match && all_passed_preselect && all_passed_pqueue"},
    { "nomatch_pass_all",     "evil", "NO MC-match && presel && pqueue",  " ! all_mc_match && all_passed_preselect && all_passed_pqueue"}
  };

  struct KuExportOpts {
    std::string prefix;      // "good_" or "evil_"
    std::string hm_idx_col;  // Which hit-match indices to use
    std::string title;       // Title postfix for canvas

    std::string index_column(const std::string &idx_col_base) const {
      return prefix + "_" + idx_col_base;
    }
    std::string column_prefix() const {
      return prefix + "_";
    }
  };

  const std::vector<KuExportOpts> kKuExports = {
    { "good", "match_pass_all_hm_indices",   "Matched that passed preselect and pqueue" },
    { "evil", "nomatch_pass_all_hm_indices", "NOT Matched that passed preselect and pqueue" }
  };

  const std::string hm_idx_col_base("hm_indices");
  const std::string ku_idx_col_base("ku_indices");
  // c++-23 is ok with regular strings ... cling is getting there, May 2026.
  constexpr const char *canvas_name_fmt = "{}_hit_match_quality_{}";
  constexpr const char *canvas_title_fmt = "{} Hit Match Quality";
  constexpr const char *canvas_ku_name_fmt = "{}_hit_kalman_update_quality_{}";
  constexpr const char *canvas_ku_title_fmt = "{} Kalman Update Quality";

} // end anon namespace

RNode AnRun::define_and_plot_hitmatch_selections(RNode r_hm, const std::string &short_name, const std::string &long_name) {

  r_hm = define_hitmatch_stuff(r_hm, hm_idx_col_base, "all_");

  // Now define index selection columns in the top rdf.
  for (auto const &sel : kSelections) {
    if (sel.is_primary()) continue;
    r_hm = r_hm.Define(sel.selection_name(hm_idx_col_base), std::format("{} [ {} ]", hm_idx_col_base, sel.index_cut));
  }

  // Second hit-match pass, define derived columns (export those with col_prefix), plot stuff.
  for (auto const &sel : kSelections) {
    RNode r = r_hm;
    // Redefine also for primary, so the "internal" column names are without the prefix.
    r = define_hitmatch_stuff(r_hm, sel.selection_name(hm_idx_col_base));
    if ( ! sel.is_primary() && sel.has_column_prefix()) {
      // If we have column prefix on the secondary, we also export it to layer RDF for later comparisons & such.
      r_hm = define_hitmatch_stuff(r_hm, sel.selection_name(hm_idx_col_base), sel.column_prefix());
    }
    auto &C = NewCanvasGroup(3, 2, std::format(canvas_name_fmt, short_name, sel.sel_prefix), std::format(canvas_title_fmt, long_name), sel.title);
    // We plot stuff without the prefix (it's not even available as we've but it into r_hm afterwards).
    plot_hitmatch_layer_stuff(r, C);
  }

  // First pass of linking kalman-update for hit-matches that passed preselect and pqueue.
  // Here we do it by hand, for now -- and export them all back into top level RDF.
  for (const auto& ku_opt : kKuExports) {
    r_hm = r_hm.Define(ku_opt.index_column(ku_idx_col_base), EV_GATHER(trHitMatches_, kalman_id, ku_opt.hm_idx_col));
  }
  // And the second pass, definitions (with export for all (for now)) + plotting
  for (const auto& ku_opt : kKuExports) {
    // Define Kalman columns and export back to r_hm
    r_hm = define_kalmanupdate_stuff(r_hm, ku_opt.index_column(ku_idx_col_base), ku_opt.column_prefix());

    // Plot
    auto &C = NewCanvasGroup(2, 2, std::format(canvas_ku_name_fmt, short_name, ku_opt.prefix),
                                   std::format(canvas_ku_title_fmt, long_name), ku_opt.title);
    // Always plot with prefix (this is different than for hit-matches -- could change that, let's see how it looks)
    plot_kalmanupdate_layer_stuff(r_hm, C, ku_opt.column_prefix());
  }

  return r_hm;
}

std::vector<AnRun::KuStartingPoint> AnRun::CollectKuStartingPoints(RNode r, const std::string& ku_idx_col, bool accepted_only) {

  std::vector<KuStartingPoint> points;

  r.Foreach([&](mkfit::Event* ev, const ROOT::RVec<int>& ku_idx) {
  for (size_t i = 0; i < ku_idx.size(); ++i) {
    int ku_id = ku_idx[i];
    const auto& ku = ev->trKalmanUpdates_[ku_id];
    bool accepted = ku.accepted;
    if (accepted || ! accepted_only) {
      int state_id = ku.state_id_in;
      const auto& state = ev->trCandStates_[state_id];
      int meta_id = state.meta_id;

      points.emplace_back(ev->evtID(), meta_id, state_id, ku_id, accepted);
    }
  }
  }, {"event", ku_idx_col});

  return points;
}

void AnRun::TraceAndPrintAcceptedPath(const KuStartingPoint& start) {
  const auto* ev = get_event_ptr(start.evtID);

  an_printf("\n");
  an_printf("================================================================================\n");
  an_printf("=== Accepted Kalman Path Trace\n");
  an_printf("=== Event %d | Meta %d | Start State %d | Start KU %d | Accepted=%d\n",
         start.evtID, start.meta_id, start.state_id, start.ku_id, start.accepted);
  an_printf("================================================================================\n");

  int current_ku_id = start.ku_id;
  int current_state_id = start.state_id;
  int step = 0;

  while (current_ku_id >= 0 && current_ku_id < (int)ev->trKalmanUpdates_.size()) {
    const auto& ku = ev->trKalmanUpdates_[current_ku_id];
    const auto& state = ev->trCandStates_[current_state_id];
    const auto& hm = ev->trHitMatches_[ku.hit_match_id];

    an_printf("\n--- Step %d ---\n", step++);

    an_printf("[Candidate State %d]\n", current_state_id);
    mkfit::print("  ", state);

    an_printf("[Hit Match %d]\n", ku.hit_match_id);
    mkfit::print("  ", hm);

    an_printf("[Kalman Update %d]\n", current_ku_id);
    mkfit::print("  ", ku);

    // Check if this KU created an output state
    if (ku.state_id_out >= 0 && ku.state_id_out < (int)ev->trCandStates_.size()) {
      // const auto& next_state = ev->trCandStates_[ku.state_id_out];

      // Find next accepted KU from this state
      int next_ku_id = -1;
      auto ku_it = ev->trKalmanUpdatesByState_.find(ku.state_id_out);
      if (ku_it != ev->trKalmanUpdatesByState_.end()) {
        for (int kid : ku_it->second) {
          if (ev->trKalmanUpdates_[kid].accepted) {
            next_ku_id = kid;
            break;
          }
        }
      }

      if (next_ku_id < 0) {
        an_printf("  → END OF ACCEPTED PATH (no further accepted KUs)\n");
        break;
      }

      current_ku_id = next_ku_id;
      current_state_id = ku.state_id_out;
    } else {
      an_printf("  → END OF PATH (no output state)\n");
      break;
    }
  }

  an_printf("\n================================================================================\n");
  an_printf("=== Path Trace Complete | Total Steps: %d\n", step);
  an_printf("================================================================================\n\n");
}

#pragma endregion
//==============================================================================
#pragma region Canvas stuff
//==============================================================================

void AnRun::DrawCanvasGroups() {
  for (auto &cg : m_canvas_groups) {
    cg->Draw();
  }
}

void AnRun::WriteCanvasGroupsToFile(const std::string &fname) const {
  TFile f(fname.c_str(), "RECREATE");
  for (auto &cg : m_canvas_groups) {
    // cg->m_canvas->Write(cg->m_canvas->GetName());

    // Generate unique name for canvases in the file so they can be drawn again.
    // Note that just using a different key name is not enough -- the original name is saved in the file.
    std::string name(cg->m_canvas->GetName());
    std::string new_name = name + "_";
    cg->m_canvas->SetName(new_name.c_str());
    cg->m_canvas->Write();
    cg->m_canvas->SetName(name.c_str());
  }
  f.Close();
  an_printf("%s Wrote %d canvases into file '%s'\n", __func__, (int) m_canvas_groups.size(), fname.c_str());

  // Grr, shows empty histograms. It's not global canvas names ... what else could it be?
  new TFile(fname.c_str());
  new TBrowser();
}

void AnRun::WriteCanvasGroupsToFile() {
  WriteCanvasGroupsToFile(m_prefix + ".root");
  // Log is complete at this point -- flush and close it so <prefix>.txt is final.
  an_log_close();
}

#pragma endregion

//==============================================================================
