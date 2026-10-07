#include "CanvasGroup.h"

#include "TCanvas.h"
#include "TPaveText.h"

TCanvas *canvas_ptr = nullptr;
int canvas_i = -1, canvas_imax = -1;

void np() {
  if (++canvas_i > canvas_imax) canvas_i = 1;
  canvas_ptr->cd(canvas_i);
}

void nc(int dx, int dy, const char *n, const char *t, const char *pfx) {
  static int n_canvases = 0;
  ++n_canvases;
  TString name, title;
  if (n == 0) {
    name.Form("c_%d", n_canvases);
  } else {
    name = n;
  }
  if (t == 0) {
    title.Form("Canvas %d", n_canvases);
  } else {
    if (pfx)
      title.Form("%s -- %s", t, pfx);
    else
      title = t;
  }
  canvas_ptr = new TCanvas(name, title);
  canvas_imax = dx * dy;
  canvas_i = 0;
  if (dx>1 || dy>1) {
    canvas_ptr->Divide(dx, dy);
    np();
  }
}

// void logx(int f=1) { gPad->SetLogx(f); }
// void logy(int f=1) { gPad->SetLogy(f); }
// void logz(int f=1) { gPad->SetLogz(f); }

// ================================================================

int CanvasGroup::s_canvas_width = 1200;
int CanvasGroup::s_canvas_height = 800;
int CanvasGroup::s_canvas_counter = 0;
bool CanvasGroup::s_do_titles = true;

const CanvasGroup::predraw_mod_func CanvasGroup::stats = [](TH1 *h, TVirtualPad *p) { h->SetStats(1); };
const CanvasGroup::predraw_mod_func CanvasGroup::logx = [](TH1 *h, TVirtualPad *p) { p->SetLogx(1); };
const CanvasGroup::predraw_mod_func CanvasGroup::logy = [](TH1 *h, TVirtualPad *p) { p->SetLogy(1); };
const CanvasGroup::predraw_mod_func CanvasGroup::logz = [](TH1 *h, TVirtualPad *p) { p->SetLogz(1); };

CanvasGroup::CanvasGroup(const std::string &n, const std::string &t, const std::string &pfx)
{
  MakeCanvas(n, t, pfx);
}

CanvasGroup::CanvasGroup(int dx, int dy, const std::string &n, const std::string &t, const std::string &pfx)
{
  MakeCanvas(n, t, pfx);
  Divide(dx, dy);
}

CanvasGroup& CanvasGroup::MakeCanvas(const std::string &n, const std::string &t, const std::string &pfx) {
  s_canvas_counter++;
  TString name, title;
  if (n.empty()) {
    name.Form("cg_%d", s_canvas_counter);
  } else {
    name = n;
  }
  if (t.empty()) {
    title.Form("Canvas Group %d", s_canvas_counter);
  } else {
    if ( ! pfx.empty())
      title.Form("%s -- %s", t.c_str(), pfx.c_str());
    else
      title = t;
  }
  m_canvas = new TCanvas(name, title, s_canvas_width, s_canvas_height);

  if (s_do_titles)
  {
    // Title pad
    m_canvas->cd();
    m_title_pad = new TPad(name + "_title", "title pad", 0.0, 0.95, 1.0, 1.0);
    m_title_pad->SetBorderSize(1);
    m_title_pad->SetFillColor(0);
    m_title_pad->Draw();
    // Title text
    m_title_pad->cd();
    m_title_text = new TPaveText(0.0, 0.0, 1.0, 1.0, "NDC");
    m_title_text->AddText(title);
    m_title_text->SetTextAlign(22);  // Center
    m_title_text->SetTextFont(42);
    m_title_text->SetTextSize(0.5);
    m_title_text->SetFillColor(0);
    m_title_text->SetBorderSize(0);
    m_title_text->Draw();

    // Plot pad for histograms (95% of canvas)
    m_canvas->cd();
    m_plot_pad = new TPad(name + "_plot", "plot pad", 0.0, 0.0, 1.0, 0.95);
    m_plot_pad->Draw();
  }
  else
  {
    m_plot_pad = m_canvas;
  }

  return *this;
}

void CanvasGroup::Divide(int dx, int dy) {
  if (dx < 1 || dy < 1) throw std::runtime_error("attempting to divide canvas with non-positive dimensions");
  m_plot_pad->Divide(dx, dy);
  m_n_divs = dx * dy;
}

void CanvasGroup::DivideSquare(int n) {
  if (n <= 0)
    n = m_entries.size();
  int w = 1, h = 1;
  double nsqrt = std::sqrt((double)n);
  if (m_canvas->GetWindowWidth() > m_canvas->GetWindowHeight()) {
    w = std::ceil(nsqrt);
    h = std::floor(nsqrt);
    if (w*h < n) w++;
  } else {
    h = std::ceil(nsqrt);
    w = std::floor(nsqrt);
    if (w*h < n) h++;
  }
  Divide(w, h);
}

void CanvasGroup::Draw() {
  if ( ! m_canvas)
    MakeCanvas();
  if (m_n_divs < 1)
    DivideSquare();

  int N_draw = std::min(m_n_divs, (int) m_entries.size());
  if (N_draw < (int) m_entries.size())
    printf("*** WARNING *** CanvasGroup::Draw: For canvas '%s' only drawing %d/%zu histograms!\n",
           m_canvas->GetName(), N_draw, m_entries.size());

  for (int i = 0; i < N_draw; ++i)
  {
    auto &entry = m_entries[i];
    if ( ! entry.raw && ! entry.histo ) continue;

    TVirtualPad* pad = m_plot_pad->cd(i + 1);
    TH1* h_ptr = entry.raw ? entry.raw : entry.histo.GetPtr();

    for (auto& func : entry.pre_funcs)
      func(h_ptr, pad);

    h_ptr->Draw(entry.options.c_str());

    // ??? This might need to be done AFTER drawin, in another pass
    // Also, will probably need to do pad->Modified(), pad/canvas->Update() after each modification, or at least at the end of all modifications
    // Return from post-func could steer this decision
    for (auto& func : entry.post_funcs) {
      func(h_ptr, pad);
    }
  }
}
