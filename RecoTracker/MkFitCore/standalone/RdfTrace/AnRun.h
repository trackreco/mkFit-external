#ifndef RecoTracker_MkFitCore_standalone_RdfTrace_AnRun_h
#define RecoTracker_MkFitCore_standalone_RdfTrace_AnRun_h

#include "RecoTracker/MkFitCore/standalone/RdfTrace/RdfSources.h"
#include "RecoTracker/MkFitCore/standalone/RdfTrace/CanvasGroup.h"

#include "ROOT/RDataFrame.hxx"

#include <optional>

class CanvasGroup;

// Tee-printf for analysis output: always goes to stdout, and additionally into
// the log file when one is open (see AnRun::SetPrefix()). Free function, so it
// is also usable from non-capturing lambdas passed to RDF.
int an_printf(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
void an_log_open(const std::string &fname);
void an_log_close();

struct AnRun {
  using ANode = std::optional<ROOT::RDF::RNode>;
  using RNode = ROOT::RDF::RNode;
  using RVecI = ROOT::RVec<int>;
  using RVecF = ROOT::RVec<float>;

  std::vector<std::unique_ptr<CanvasGroup>> m_canvas_groups;

  mkfit::RdfCtx CTX;

  AnRun(const mkfit::TrackerInfo& ti) : CTX(mkfit::MakeCtx(nullptr, ti))
  {}

  AnRun(const mkfit::Event* ev, const mkfit::TrackerInfo& ti) :
    CTX(mkfit::MakeCtx(ev, ti))
  {}

  ~AnRun();

  // ----- Output naming: <prefix>.root for canvases, <prefix>.txt for printouts.

  std::string m_prefix = "mkfit";

  // Sets the prefix and (re)opens <prefix>.txt as the an_printf() log.
  void SetPrefix(const std::string &p);
  const std::string& Prefix() const { return m_prefix; }

  // -----

  ANode m_rdf_hitmatch;
  ANode m_rdf_meta;
  ANode f_resdy_400;

  void RunOldVecBased();

  // -----

  ANode m_rdf_event;
  std::vector<const mkfit::Event*> m_ev_vec;
  const mkfit::Event* get_event_ptr(int event_id) const;

  // Swaps out ev_vec and owns it.
  void SetupRdfEvent(std::vector<const mkfit::Event*>& ev_vec);

  void RunBasicSeedCandCheck();

  // Per-sim-track efficiency (MTV selection) straight off the Event, plus the
  // seeding efficiency that bounds it. Single configuration by construction --
  // see the comment on the definition.
  void RunSimTrackEfficiency();

  void RunMetaVsSeedDuplicateCheck();

  void Run_T5_vs_pT5_AsSeeds_DuplicateCount();

  void Run_T5s_into_Pix();
  ANode m_T5;

  void Run_Stage2_RootState_QualityCheck();
  ANode m_Quality;

  void Run_Stage2_RootState_Covariance_Check();
  ANode m_Cov;

  // ----- HitMatch helpers

  static RNode define_hitmatch_indices_by_layers(RNode r, const std::unordered_set<int> &layers);

  static RNode define_hitmatch_indices_by_sim_pixel_layer(RNode r, int layer_offset, bool use_head = true);

  static RNode define_hitmatch_stuff(RNode r, const std::string &idx_column, const std::string &pref = "");
  static RNode define_kalmanupdate_stuff(RNode r, const std::string &idx_column, const std::string &pref = "");

  static void plot_hitmatch_layer_stuff(RNode r, CanvasGroup &C, const std::string &pref = "");
  static void plot_kalmanupdate_layer_stuff(RNode r, CanvasGroup &C, const std::string &pref = "");

  RNode define_and_plot_hitmatch_selections(RNode r_hm, const std::string &short_name, const std::string &long_name);

  // --- Kalman update helpers

  struct KuStartingPoint {
    int evtID;
    int meta_id;
    int state_id;
    int ku_id;
    bool accepted;
  };

  std::vector<KuStartingPoint> CollectKuStartingPoints(RNode r, const std::string& ku_idcx_col, bool accepted_only);

  void TraceAndPrintAcceptedPath(const KuStartingPoint& start);

  // ===== CanvasGroup management

  CanvasGroup& NewCanvasGroup(const std::string &n="", const std::string &t="", const std::string &pfx="") {
    m_canvas_groups.emplace_back( std::make_unique<CanvasGroup>(n,t,pfx) );
    return *m_canvas_groups.back();
  }
  CanvasGroup& NewCanvasGroup(int dx=1, int dy=1, const std::string &n="", const std::string &t="", const std::string &pfx="") {
    m_canvas_groups.emplace_back( std::make_unique<CanvasGroup>(dx, dy, n, t, pfx) );
    return *m_canvas_groups.back();
  }

  void DrawCanvasGroups();
  void WriteCanvasGroupsToFile(const std::string &fname) const;
  // Writes <prefix>.root and closes the <prefix>.txt log.
  void WriteCanvasGroupsToFile();
};

#endif
