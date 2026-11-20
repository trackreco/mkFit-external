// Shell driver of the probes, Probe.h. The work is in libMkFitCMS.so; this
// file only gives cling the declarations to call through, so that no Matriplex
// header has to be parsed by the interpreter. Load it in mkFit's --shell, from
// a build directory, with <SACMS> the directory its Makefile exports as SACMS:
//
//   unset DISPLAY; export LD_LIBRARY_PATH=.
//   echo .q | ./mkFit --geom CMS-phase2 --input-file <any .bin> --num-events 1 \
//     --num-thr 1 --shell \
//     --shell-command 'gROOT->SetBatch(kTRUE)' \
//     --shell-command 'gROOT->ProcessLine(".L <SACMS>/Probe/probe.C")' \
//     --shell-command 'probe_gen("probe-closure.root")'
//
// No geometry is used by the generator; the shell just wants a file to start.
// The probes of an event take s.event() after s.ProcessEventStd(), e.g.
// probe_eff_ev(s.event(), "cfg").

#include "RecoTracker/MkFitCMS/standalone/Probe/Probe.h"

// The probes are called by their library names, without mkfit::. The knobs
// are Shell methods, s.SetMaxCands(6) and so on (Shell.h, section Knobs).
using namespace mkfit;

void probe_gen(const char *out = "probe-closure.root") { probe_gen_closure(out); }

void probe_dalpha_ev(const Event *ev, int which = 0, float pt_min = 0.5) { probe_dalpha_add_event(ev, which, pt_min); }

void probe_info(const Event *ev) { probe_sample_info(ev); }

void probe_covx(const char *f = "probe-covxport.root", int n = 10000) { probe_cov_transport(f, n); }
void probe_trkdump_ev(const Event *ev, const char *cfg) { probe_trkdump_event(ev, cfg); }
void probe_te_reset() { probe_track_eff_reset(); }
void probe_te_ev(const Event *e) { probe_track_eff_event(e); }
void probe_te_report(const char *t) { probe_track_eff_report(t); }

void probe_geom_reset() { probe_geom_check_reset(); }
void probe_geom_ev(const Event *ev) { probe_geom_check_event(ev); }
void probe_geom_report(const char *name = "") { probe_geom_check_report(name); }

void probe_eff_ev(const Event *e, const char *c) { probe_eff_event(e, c); }
void probe_eff_cmssw_ev(const Event *e, const char *c) { probe_eff_cmssw_event(e, c); }

#ifdef MKFIT_TRACE
void probe_rank_ev(const Event *ev, const char *cfg) { probe_rank_event(ev, cfg); }
void probe_kinfo_ev(const Event *e) { probe_kinfo_event(e); }
void probe_sr_reset() { probe_search_reset(); }
void probe_sr_ev(const Event *e, int i) { probe_search_event(e, i); }
void probe_sr_write(const char *f = "probe-search.root") { probe_search_write(f); }
void probe_bt_reset() { probe_bkfit_trace_reset(); }
void probe_bt_ev(const Event *e) { probe_bkfit_trace_event(e); }
void probe_bt_report(const char *p = "probe-bkfit-trace") { probe_bkfit_trace_report(p); }
#endif
