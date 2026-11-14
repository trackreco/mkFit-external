#ifndef RecoTracker_MkFitCore_standalone_RdfTrace_CanvasGroup_h
#define RecoTracker_MkFitCore_standalone_RdfTrace_CanvasGroup_h

#include "RecoTracker/MkFitCore/standalone/Event.h"

#include "ROOT/RDataFrame.hxx"

class TCanvas;
class TVirtualPad;
class TPaveText;

extern TCanvas *canvas_ptr;
extern int canvas_i, canvas_imax;

void np();
void nc(int dx=1, int dy=1, const char *n=0, const char *t=0, const char *pfx=0);

// CanvasCroup -- to hold a group of histograms and draw them together
struct CanvasGroup
{
  using predraw_mod_func = std::function<void(TH1*, TVirtualPad*)>;
  using postdraw_mod_func = std::function<void(TH1*, TVirtualPad*)>; // return bool to signal modified/update?
  using vec_predraw_mod_func = std::vector<predraw_mod_func>;
  using vec_postdraw_mod_func = std::vector<postdraw_mod_func>;

  struct Entry {
    ROOT::RDF::RResultPtr<TH1> histo;
    // An efficiency is a RATIO of two booked histograms, so it does not exist
    // until after the event loop and can never be an RResultPtr. Entries made by
    // AddTH1() carry the finished object here instead, and Draw() prefers it.
    TH1 *raw = nullptr;
    std::string options;
    std::vector<predraw_mod_func> pre_funcs;
    std::vector<postdraw_mod_func> post_funcs;

    Entry& add_pre(predraw_mod_func func) { pre_funcs.push_back(func); return *this;}
    Entry& add_post(postdraw_mod_func func) { post_funcs.push_back(func); return *this; }
    Entry& add_pre(vec_predraw_mod_func &funcs) { pre_funcs.insert(pre_funcs.end(), funcs.begin(), funcs.end()); return *this; }
    Entry& add_post(vec_postdraw_mod_func funcs) { post_funcs.insert(post_funcs.end(), funcs.begin(), funcs.end()); return *this; }
  };

  std::vector<Entry> m_entries;
  TCanvas *m_canvas = nullptr;
  TVirtualPad *m_title_pad = nullptr;
  TVirtualPad *m_plot_pad = nullptr;
  TPaveText *m_title_text = nullptr;
  int m_n_divs = -1;

  // CanvasGroup(const char *n=0, const char *t=0, const char *pfx=0);
  // CanvasGroup(int dx, int dy, const char *n=0, const char *t=0, const char *pfx=0);

  CanvasGroup(const std::string &n="", const std::string &t="", const std::string &pfx="");
  CanvasGroup(int dx, int dy, const std::string &n="", const std::string &t="", const std::string &pfx="");


  Entry& Add(ROOT::RDF::RResultPtr<TH1> histo, const std::string &opts = "") {
    m_entries.push_back( { histo, nullptr, opts, {}, {} } );
    return m_entries.back();
  }

  // For histograms the dataframe cannot produce -- ratios, efficiencies,
  // differences between two configurations. The group does NOT take ownership;
  // the caller keeps the object alive, which is what writing it to a TFile does.
  Entry& AddTH1(TH1 *histo, const std::string &opts = "") {
    m_entries.push_back( { {}, histo, opts, {}, {} } );
    return m_entries.back();
  }

  Entry& AddRealH1D(ROOT::RDF::RNode &r, const std::string &column, int n_bins, double min, double max, const std::string &opts = "") {
    std::string t = column + ";" + column + ";N";
    return Add(r.Histo1D({column.c_str(), t.c_str(), n_bins, min, max}, column), opts);
  }
  Entry& AddIntH1D(ROOT::RDF::RNode &r, const std::string &column, int min, int max, const std::string &opts = "") {
    return AddRealH1D(r, column, max - min + 1, min - 0.5, max + 0.5, opts);
  }

  Entry& AddRealH2D(ROOT::RDF::RNode &r, const std::string &column_x, const std::string &column_y, int n_bins_x, double min_x, double max_x, int n_bins_y, double min_y, double max_y, const std::string &opts = "") {
    std::string n = column_y + "_VS_" + column_x;
    std::string t = n + ";" + column_x + ";" + column_y + ";N";
    return Add(r.Histo2D({n.c_str(), t.c_str(), n_bins_x, min_x, max_x, n_bins_y, min_y, max_y}, column_x, column_y), opts);
  }
  Entry& AddRealIntH2D(ROOT::RDF::RNode &r, const std::string &column_x, const std::string &column_y, int n_bins_x, double min_x, double max_x, int min_y, int max_y, const std::string &opts = "") {
    return AddRealH2D(r, column_x, column_y, n_bins_x, min_x, max_x, max_y - min_y + 1, min_y - 0.5, max_y + 0.5, opts);
  }
  Entry& AddIntRealH2D(ROOT::RDF::RNode &r, const std::string &column_x, const std::string &column_y, int min_x, int max_x, int n_bins_y, double min_y, double max_y, const std::string &opts = "") {
    return AddRealH2D(r, column_x, column_y, max_x - min_x + 1, min_x - 0.5, max_x + 0.5, n_bins_y, min_y, max_y, opts);
  }

  Entry& Entry() {
    if ( ! m_entries.empty() )
      return m_entries.back();
    throw std::runtime_error("CanvasGroup empty when requesting current entry.");
  }

  CanvasGroup& MakeCanvas(const std::string &n="", const std::string &t="", const std::string &pfx="");
  void Divide(int dx, int dy);
  void DivideSquare(int n = 0);

  void Draw();

  static int s_canvas_width;
  static int s_canvas_height;
  static int s_canvas_counter;
  static bool s_do_titles;

  static const predraw_mod_func stats;
  static const predraw_mod_func logx, logy, logz;

  // Useful?
  // static const postdraw_mod_func post_stats = [](TH1 *h, TVirtualPad *p) { h->SetStats(1); p->Modified(); p->Update(); };
};

typedef CanvasGroup CGrp;

#endif
