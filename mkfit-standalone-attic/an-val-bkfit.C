// chi2 through the backward fit, hit by hit, on plain T5 seeds.
//
//   LD_LIBRARY_PATH=. root.exe -l -b -q \
//     '../RecoTracker/MkFitCore/standalone/test/an-val-bkfit.C("val-bkfit-t5.root")'
//
// Expectation: each hit is a TWO-dimensional measurement in the module plane,
// so chi2 per hit ~ 2 and the track total ~ 2*n_hits. Step 0 is excluded: the
// seed state already sits on the outermost hit's plane (to 5e-19 cm), so that
// update re-uses the hit that produced the state and its chi2 is not a
// measurement of anything.

#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"

void an_val_bkfit(const char *in_file = "val-bkfit-t5.root",
                  const char *prefix = "200-bkfit") {
  gSystem->Load("libMkFitRootDataFormats.so");
  FILE *log = fopen(Form("%s.txt", prefix), "w");
  auto P = [&](const char *fmt, ...) { va_list a; va_start(a,fmt); vprintf(fmt,a); va_end(a);
    va_start(a,fmt); vfprintf(log,fmt,a); va_end(a); };

  // --- read it all in; it is small, and the track grouping is easier by hand
  TFile f(in_file); TTree *t = (TTree*)f.Get("bkfit");
  ValFitHit *h = nullptr; t->SetBranchAddress("h", &h);
  struct Trk { std::vector<ValFitHit> hits; };
  std::map<std::pair<int,int>, Trk> trks;
  const Long64_t N = t->GetEntries();
  for (Long64_t i = 0; i < N; ++i) { t->GetEntry(i); trks[{h->event, h->label}].hits.push_back(*h); }

  P("### %s -- chi2 through the backward fit, T5 seeds only\n###\n", prefix);
  P("### definition : per-hit chi2 from kalmanOperationPlaneLocal, and the track\n");
  P("###              total. Each hit is a 2-D measurement in the module plane,\n");
  P("###              so chi2/hit has expectation ~2 and the total ~2*n_hits.\n");
  P("### population : %lld hit records, %zu tracks, 20 events, LST T5 seeds,\n",
    (long long)N, trks.size());
  P("###              backward fit only -- no forward search (it is commented\n");
  P("###              out in Shell-LST.cc), so the hits are LST's, not ours.\n");
  P("### step 0     : EXCLUDED. The seed state lies on the outermost hit's plane\n");
  P("###              to 5e-19 cm, so that update re-uses the hit that made it.\n###\n");

  std::vector<double> red;            // chi2/(2*n) per track
  std::map<std::pair<int,int>, double> track_red;
  int n_short = 0;
  for (auto &kv : trks) {
    auto &v = kv.second.hits;
    std::sort(v.begin(), v.end(), [](const ValFitHit&a, const ValFitHit&b){ return a.step < b.step; });
    double c2 = 0; int n = 0;
    for (auto &x : v) if (x.step > 0 && std::isfinite(x.chi2)) { c2 += x.chi2; ++n; }
    if (n < 2) { ++n_short; continue; }
    const double r = c2 / (2.0 * n);
    red.push_back(r); track_red[kv.first] = r;
  }
  std::sort(red.begin(), red.end());
  auto q=[&](std::vector<double>&u,double p){ return u.empty()?0.0:u[(size_t)(0.01*p*(u.size()-1))]; };
  P("tracks with >=2 usable hits: %zu  (%d too short)\n", red.size(), n_short);
  P("chi2 / (2 * n_hits), per track:\n");
  P("  p10 %.3g   p25 %.3g   median %.3g   p75 %.3g   p90 %.3g   p99 %.3g   max %.3g\n",
    q(red,10), q(red,25), q(red,50), q(red,75), q(red,90), q(red,99), red.empty()?0:red.back());
  P("  expectation ~1.0\n\n");

  // --- three categories by the SEED's truth purity, not by chi2 itself.
  // good_frac is the fraction of the seed's valid hits that come from its
  // best-matching sim track (Event::simInfoForTrack). Splitting on chi2 would
  // only slice the symptom by itself; splitting on purity tests a cause.
  const char *cname[3] = {"pure  (gf = 1.0)", "mostly (gf >= 0.8)", "dirty (gf < 0.8)"};
  auto cat_of = [](float gf) { return gf >= 0.9999f ? 0 : (gf >= 0.8f ? 1 : 2); };

  std::vector<double> per_hit[3], red_cat[3], per_hit_mc[3][2];
  long nfail[3] = {0,0,0}, ntot[3] = {0,0,0}, ntrk[3] = {0,0,0}, n_nogf = 0;
  for (auto &kv : trks) {
    const float gf = kv.second.hits.empty() ? -1.f : kv.second.hits[0].good_frac;
    if (gf < 0) { ++n_nogf; continue; }
    const int cat = cat_of(gf);
    auto it = track_red.find(kv.first);
    if (it != track_red.end()) { red_cat[cat].push_back(it->second); ++ntrk[cat]; }
    for (auto &x : kv.second.hits) {
      if (x.step == 0 || !std::isfinite(x.chi2)) continue;
      per_hit[cat].push_back(x.chi2);
      ++ntot[cat]; if (x.fail) ++nfail[cat];
      per_hit_mc[cat][(x.mc_track_id == x.sim_label) ? 1 : 0].push_back(x.chi2);
    }
  }
  P("tracks with no truth match at all: %ld\n\n", n_nogf);
  P("PER TRACK, chi2/(2*n_hits) -- expectation 1.0\n");
  P("  %-20s %7s %9s %9s %9s %9s %9s\n","category","tracks","p25","median","p75","p90","max");
  for (int c = 0; c < 3; ++c) {
    auto &u = red_cat[c]; std::sort(u.begin(), u.end());
    P("  %-20s %7ld %9.3g %9.3g %9.3g %9.3g %9.3g\n", cname[c], ntrk[c],
      q(u,25), q(u,50), q(u,75), q(u,90), u.empty()?0:u.back());
  }
  P("\nPER HIT, chi2 -- expectation ~2 (a 2-D measurement), median of chi2_2 = 1.39\n");
  P("  %-20s %7s %9s %9s %9s %9s %9s\n","category","hits","p25","median","p75","p90","max");
  for (int c = 0; c < 3; ++c) {
    auto &u = per_hit[c]; std::sort(u.begin(), u.end());
    P("  %-20s %7zu %9.3g %9.3g %9.3g %9.3g %9.3g\n", cname[c], u.size(),
      q(u,25), q(u,50), q(u,75), q(u,90), u.empty()?0:u.back());
  }
  P("\nPER HIT, split by whether THAT HIT's mcTrackID matches the track's sim label\n");
  P("  %-20s %8s %9s %9s | %8s %9s %9s\n","category",
    "n match","med","p90","n other","med","p90");
  for (int c = 0; c < 3; ++c) {
    auto &m = per_hit_mc[c][1]; std::sort(m.begin(), m.end());
    auto &o = per_hit_mc[c][0]; std::sort(o.begin(), o.end());
    P("  %-20s %8zu %9.3g %9.3g | %8zu %9.3g %9.3g\n", cname[c],
      m.size(), q(m,50), q(m,90), o.size(), q(o,50), q(o,90));
  }
  P("\n  propagation failures: pure %ld/%ld, mostly %ld/%ld, dirty %ld/%ld\n",
    nfail[0],ntot[0], nfail[1],ntot[1], nfail[2],ntot[2]);
  fclose(log);
}
