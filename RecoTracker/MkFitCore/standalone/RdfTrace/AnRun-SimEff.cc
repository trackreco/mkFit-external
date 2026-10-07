#include "AnRun.h"
#include "AnRunEvMap.h"

#include "RecoTracker/MkFitCore/standalone/TrackExtra.h"

#include "TH1D.h"

#include <map>

#pragma region Sim-track efficiency
//==============================================================================

// Per-SIM-TRACK efficiency, MTV's convention, as an RDataFrame over the Event.
//
// Nothing here touches the trace graph, and it does not need to: one RDF entry
// is one Event and the Event carries simTracks_, so EV_MAP projects a sim-track
// quantity into a column exactly as it does a seed one. The graph is rooted at
// seeds and searches and has no sim-track node, which is a statement about the
// GRAPH and not about what this dataframe can reach.
//
// The one thing that is not a projection is the reco->sim association, because
// it needs TrackExtra and the seed-hit exclusion. That is one Define returning a
// per-sim-track count, computed once per Event.
//
// WHAT THIS DELIBERATELY DOES NOT DO: compare configurations. The error bar that
// matters on a difference is the spread of the per-EVENT difference between two
// configurations, which is a join across separate runs on event id, not a
// projection within one -- and the Event only ever holds the last
// configuration's candidateTracks_. That half lives in mkfit::val_eff_*.

namespace {
  // sim label -> number of reco tracks associated to it by quality-val's rule,
  // 2*mccount >= nCandHits over the non-seed hits.
  ROOT::RVec<int> an_sim_assoc_counts(const mkfit::Event *ev) {
    ROOT::RVec<int> n(ev->simTracks_.size(), 0);
    const mkfit::TrackVec *seeds = nullptr;
    try { seeds = &ev->currentSeedTracks(); } catch (...) { seeds = nullptr; }
    std::map<int, int> seed_by_label;
    if (seeds)
      for (int i = 0; i < (int) seeds->size(); ++i)
        seed_by_label.emplace((*seeds)[i].label(), i);
    for (const mkfit::Track &c : ev->candidateTracks_) {
      mkfit::TrackExtra extra(c.label());
      auto si = seed_by_label.find(c.label());
      if (seeds && si != seed_by_label.end())
        extra.findMatchingSeedHits(c, (*seeds)[si->second], ev->layerHits_);
      extra.setMCTrackIDInfo(c, ev->layerHits_, ev->simHitsInfo_, ev->simTracks_, false, false);
      const int mc = extra.mcTrackID();
      if (mc >= 0 && mc < (int) n.size()) ++n[mc];
    }
    return n;
  }

  // Number of seeds pointing at each sim track -- the ceiling on the efficiency,
  // since a track with no seed cannot be found.
  ROOT::RVec<int> an_sim_seed_counts(const mkfit::Event *ev) {
    ROOT::RVec<int> n(ev->simTracks_.size(), 0);
    const mkfit::TrackVec *seeds = nullptr;
    try { seeds = &ev->currentSeedTracks(); } catch (...) { return n; }
    for (int i = 0; i < (int) seeds->size(); ++i) {
      const int sl = ev->simInfoForCurrentSeed(i).label;
      if (sl >= 0 && sl < (int) n.size()) ++n[sl];
    }
    return n;
  }

  // num/den with binomial errors, kept alive for the canvas and the file.
  TH1D *an_eff_ratio(TH1 *num, TH1 *den, const char *name, const char *title) {
    TH1D *e = (TH1D*) num->Clone(name);
    e->SetTitle(title);
    e->Divide(num, den, 1.0, 1.0, "B");
    e->SetMinimum(0.0);  e->SetMaximum(1.05);  e->SetStats(0);
    e->SetLineWidth(2);  e->SetMarkerStyle(20);
    return e;
  }
}

void AnRun::RunSimTrackEfficiency() {
  auto r = (*m_rdf_event)
    .Define("sim_eta",    EV_MAP(simTracks_, momEta()))
    .Define("sim_pt",     EV_MAP(simTracks_, pT()))
    .Define("sim_nlay",   EV_MAP(simTracks_, nUniqueLayers()))
    .Define("sim_nhit",   EV_MAP(simTracks_, nFoundHits()))
    .Define("sim_vx",     EV_MAP(simTracks_, x()))
    .Define("sim_vy",     EV_MAP(simTracks_, y()))
    .Define("sim_vz",     EV_MAP(simTracks_, z()))
    .Define("sim_findable", EV_MAP(simTracks_, isFindable()))
    .Define("sim_nassoc", an_sim_assoc_counts, {"event"})
    .Define("sim_nseed",  an_sim_seed_counts,  {"event"})

    .Define("sim_aeta",   "abs(sim_eta)")
    // TrackingParticleSelector's "central": tip < 3.5 cm, lip < 30 cm.
    .Define("sim_central", "sqrt(sim_vx*sim_vx + sim_vy*sim_vy) < 3.5f && abs(sim_vz) < 30.0f")
    .Define("sim_base",   "sim_findable && sim_central && sim_nlay >= 4 && sim_aeta < 3.0f && sim_pt > 0.2f")
    .Define("sim_found",  "sim_nassoc > 0")
    .Define("sim_seeded", "sim_nseed > 0")
    // MTV releases the cut on the plotted variable, and only that one.
    .Define("sel_eta",    "sim_base && sim_pt > 0.9f")
    .Define("sel_pt",     "sim_base && sim_aeta < 2.5f")
    .Define("sel_all",    "sim_base && sim_pt > 0.9f && sim_aeta < 2.5f")

    .Define("eta_den",    "sim_aeta[sel_eta]")
    .Define("eta_num",    "sim_aeta[sel_eta && sim_found]")
    .Define("eta_seed",   "sim_aeta[sel_eta && sim_seeded]")
    .Define("pt_den",     "sim_pt[sel_pt]")
    .Define("pt_num",     "sim_pt[sel_pt && sim_found]")
    .Define("pt_seed",    "sim_pt[sel_pt && sim_seeded]")
    // hits per layer -- how much overlap the sim track's own content offers,
    // which is the axis the in-layer combinatorial pays along.
    .Define("sim_hpl",    "(ROOT::RVecF) sim_nhit / (ROOT::RVecF) sim_nlay")
    .Define("hpl_den",    "sim_hpl[sel_all]")
    .Define("hpl_num",    "sim_hpl[sel_all && sim_found]")
    // Counts on the FULL selection, both cuts applied, so the summary below is
    // the same population mkfit::val_eff_* reports and the two cannot disagree.
    .Define("n_sel",      "(int) Sum(sel_all)")
    .Define("n_seeded",   "(int) Sum(sel_all && sim_seeded)")
    .Define("n_found",    "(int) Sum(sel_all && sim_found)")
    ;

  auto h_eta_den  = r.Histo1D({"eta_den",  "sim tracks;|#eta|;tracks",      12, 0.0, 3.0}, "eta_den");
  auto h_eta_num  = r.Histo1D({"eta_num",  "found;|#eta|;tracks",           12, 0.0, 3.0}, "eta_num");
  auto h_eta_seed = r.Histo1D({"eta_seed", "seeded;|#eta|;tracks",          12, 0.0, 3.0}, "eta_seed");
  auto h_pt_den   = r.Histo1D({"pt_den",   "sim tracks;p_{T} [GeV];tracks", 40, 0.0, 10.0}, "pt_den");
  auto h_pt_num   = r.Histo1D({"pt_num",   "found;p_{T} [GeV];tracks",      40, 0.0, 10.0}, "pt_num");
  auto h_pt_seed  = r.Histo1D({"pt_seed",  "seeded;p_{T} [GeV];tracks",     40, 0.0, 10.0}, "pt_seed");
  auto h_hpl_den  = r.Histo1D({"hpl_den",  "sim tracks;sim hits / layer;tracks", 25, 1.0, 2.0}, "hpl_den");
  auto h_hpl_num  = r.Histo1D({"hpl_num",  "found;sim hits / layer;tracks",      25, 1.0, 2.0}, "hpl_num");
  auto n_sel = r.Sum<int>("n_sel"), n_seeded = r.Sum<int>("n_seeded"), n_found = r.Sum<int>("n_found");

  // One trigger of the loop, then the ratios, which cannot be RResultPtrs.
  TH1D *e_eta  = an_eff_ratio(h_eta_num.GetPtr(), h_eta_den.GetPtr(),
                              "eff_vs_eta", "efficiency vs |#eta|;|#eta|;efficiency");
  TH1D *s_eta  = an_eff_ratio(h_eta_seed.GetPtr(), h_eta_den.GetPtr(),
                              "seedeff_vs_eta", "seeding efficiency vs |#eta|;|#eta|;seeded");
  TH1D *e_pt   = an_eff_ratio(h_pt_num.GetPtr(), h_pt_den.GetPtr(),
                              "eff_vs_pt", "efficiency vs p_{T};p_{T} [GeV];efficiency");
  TH1D *s_pt   = an_eff_ratio(h_pt_seed.GetPtr(), h_pt_den.GetPtr(),
                              "seedeff_vs_pt", "seeding efficiency vs p_{T};p_{T} [GeV];seeded");
  TH1D *e_hpl  = an_eff_ratio(h_hpl_num.GetPtr(), h_hpl_den.GetPtr(),
                              "eff_vs_hitsperlayer",
                              "efficiency vs sim hits per layer;sim hits / layer;efficiency");

  { auto &C = NewCanvasGroup(3, 3, "simeff", "per-sim-track efficiency, MTV selection");
    C.AddTH1(e_eta, "E1");
    C.AddTH1(e_pt,  "E1");
    C.AddTH1(e_hpl, "E1");
    C.AddTH1(s_eta, "E1");
    C.AddTH1(s_pt,  "E1");
    C.Add(h_eta_den, "");
    C.Add(h_pt_den,  "").add_pre(CGrp::logy);
    C.Add(h_hpl_den, "");
  }

  an_printf("\n=== per-sim-track efficiency, MTV selection (|eta| < 2.5, pT > 0.9,\n"
            "    central vertex, >= 4 layers; the cut on the plotted variable released) ===\n");
  an_printf("  selected sim tracks : %d\n", n_sel.GetValue());
  an_printf("  ... seeded          : %d   (%.2f%%)\n", n_seeded.GetValue(),
            n_sel.GetValue() ? 100.0*n_seeded.GetValue()/n_sel.GetValue() : 0.0);
  an_printf("  ... found           : %d   (%.2f%%)\n", n_found.GetValue(),
            n_sel.GetValue() ? 100.0*n_found.GetValue()/n_sel.GetValue() : 0.0);
  an_printf("  Both cuts applied here, so this is the same population val_eff_* reports;\n"
            "  the eta and pT PLOTS each release their own cut, which is why their\n"
            "  integrals are larger.\n");
  an_printf("  A configuration COMPARISON does not belong here -- the Event holds only the\n"
            "  last configuration's candidateTracks_, and a paired error bar is a join\n"
            "  across runs on event id. That is mkfit::val_eff_*.\n");
}

#pragma endregion
