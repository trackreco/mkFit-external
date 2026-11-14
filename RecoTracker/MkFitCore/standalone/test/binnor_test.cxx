// Mini correctness suite for binnor.h -- the axes, their range helpers and the
// binnor's own region query. Companion to binnor_demo.cxx, which measures SPEED
// and asserts nothing.
//
//   build:  c++ -o binnor_test -O2 -std=c++20 -I../../../.. \
//             binnor_test.cxx ../../src/radix_sort.cc
//   run:    ./binnor_test          # exit 0 = all required checks passed
//
// WHY IT EXISTS. The phi bin range in MkBins.cc was hand-rolled without the
// "+1" that makes a half-open [b, e) range cover its upper end, so the top bin
// was never scanned. It cost a full spare phi bin of PHI_BIN_EXTRA_FAC to hide,
// and it was found by scanning a tuning constant rather than by a test. Every
// check below is one line and would have caught it immediately.
//
// THE PROPERTY UNDER TEST IS COVERAGE, and it is the only one that matters for
// a search window: for a range [lo, hi], the returned half-open bin range must
// contain the bin of EVERY real in [lo, hi]. Anything narrower silently drops
// hits; anything wider only costs time.
//
// Checks marked KNOWN are documented defects, reported but not failing the
// suite. When one turns green, promote it to CHECK.
//
// WHERE TO ADD MORE. The harness is the three functions below check()/
// check_known()/section() and nothing else -- no framework, no dependencies.
// A new area is a new t_*() function plus one line in main(). The seed maker
// will want its own section here: it uses the same axes and the same region
// query, so the coverage predicate range_holds() applies to it unchanged, and
// the brute-force comparison in t_binnor_query_roundtrip() is the template for
// "the fast path returns exactly what the slow path does".
//
// NOTE ON THE CIRCLE. A range on S^1 is an ARC. The full circle is not an arc
// and is deliberately NOT expressible: -pi and +pi are the same point, so
// (-pi, +pi) asks for a point and correctly returns one bin. What that makes a
// caller responsible for is the PRECONDITION -- a half-width at or above pi
// wraps and silently returns an arbitrary small arc instead of everything.

#include "RecoTracker/MkFitCore/interface/binnor.h"
#include "RecoTracker/MkFitCore/interface/radix_sort.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <algorithm>
#include <type_traits>
#include <vector>

using namespace mkfit;

//------------------------------------------------------------------ harness

static int g_pass = 0, g_fail = 0, g_known = 0, g_fixed = 0;

static void check(bool ok, const std::string &what, const std::string &detail = "") {
  if (ok) { ++g_pass; return; }
  ++g_fail;
  printf("  FAIL  %s%s%s\n", what.c_str(), detail.empty() ? "" : "  --  ", detail.c_str());
}

// A documented defect: expected to fail. Reports loudly if it starts passing.
static void check_known(bool ok, const std::string &what, const std::string &why) {
  if (ok) { ++g_fixed; printf("  FIXED %s  -- was: %s  (promote to check())\n", what.c_str(), why.c_str()); }
  else    { ++g_known; printf("  known %s  -- %s\n", what.c_str(), why.c_str()); }
}

static void section(const char *name) { printf("\n[%s]\n", name); }

//------------------------------------------------------- coverage predicate

// Does the half-open bin range [b, e) contain bin `t`?  `periodic` selects wrap
// semantics, in which b == e is read as EMPTY -- which is what the scan loops
// `for (pi = p1; pi != p2; pi = mask(pi+1))` actually do.
template <typename I>
static bool range_holds(I b, I e, I t, unsigned n_bins, bool periodic) {
  if (!periodic) return t >= b && t < e;
  if (b == e) return false;               // the ambiguity: empty, not full
  if (b < e)  return t >= b && t < e;
  return t >= b || t < e;                 // wrapped
}

template <typename I>
static unsigned range_size(I b, I e, unsigned n_bins, bool periodic) {
  if (!periodic) return e > b ? unsigned(e - b) : 0u;
  if (b == e) return 0u;
  return b < e ? unsigned(e - b) : unsigned(n_bins - b + e);
}

//--------------------------------------------------------------- the axes

using phi_axis_t = axis_pow2_u1<float, unsigned short, 16, 8>;   // as HitStructures.h:30
using q_axis_t   = axis<float, unsigned short, 16, 8>;

static constexpr unsigned kPhiBins = 1u << 8;

//------------------------------------------------------------------ tests

// 1. The periodic axis's own range helper must cover.  This is the idiom that
//    exists in binnor.h and that NOTHING in the tree calls.
static void t_periodic_helper_covers() {
  section("periodic axis: from_R_minmax_to_N_bins covers its range");
  phi_axis_t ax(-M_PI, M_PI);
  std::mt19937 rng(1234);
  std::uniform_real_distribution<float> uphi(-M_PI, M_PI);
  std::uniform_real_distribution<float> uw(0.0f, 1.0f);   // half-width, rad

  int bad = 0, bad_wrap = 0;
  for (int i = 0; i < 20000; ++i) {
    float c = uphi(rng), w = uw(rng);
    float lo = c - w, hi = c + w;
    auto p = ax.from_R_minmax_to_N_bins(lo, hi);
    // sample the closed range densely and require every sample's bin to be in
    for (int k = 0; k <= 64; ++k) {
      float x = lo + (hi - lo) * k / 64.0f;
      unsigned short t = ax.from_R_to_N_bin_safe(x);
      if (!range_holds<unsigned short>(p.begin, p.end, t, kPhiBins, true)) {
        ++bad;
        if (p.begin > p.end) ++bad_wrap;
        break;
      }
    }
  }
  check(bad == 0, "periodic helper covers", "misses on " + std::to_string(bad) +
        " of 20000 ranges (" + std::to_string(bad_wrap) + " wrapped)");
}

// 2. The bounded axis's helper must cover, and must clamp rather than wrap.
static void t_bounded_helper_covers() {
  section("bounded axis: from_R_minmax_to_N_bins covers and clamps");
  q_axis_t ax(-120.0f, 120.0f, 1u << 8);
  std::mt19937 rng(99);
  std::uniform_real_distribution<float> uq(-140.0f, 140.0f);   // deliberately past the ends

  int bad = 0;
  for (int i = 0; i < 20000; ++i) {
    float a = uq(rng), b = uq(rng);
    float lo = std::min(a, b), hi = std::max(a, b);
    auto p = ax.from_R_minmax_to_N_bins(lo, hi);
    for (int k = 0; k <= 64; ++k) {
      float x = lo + (hi - lo) * k / 64.0f;
      if (x < -120.0f || x > 120.0f) continue;               // outside the axis
      unsigned short t = ax.from_R_to_N_bin_safe(x);
      if (!range_holds<unsigned short>(p.begin, p.end, t, 1u << 8, false)) { ++bad; break; }
    }
  }
  check(bad == 0, "bounded helper covers", std::to_string(bad) + " misses");
}

// 3. THE REGRESSION. Compare the three call-site idioms against coverage.
static void t_callsite_idioms() {
  section("call-site idioms for a phi range");
  phi_axis_t ax(-M_PI, M_PI);
  std::mt19937 rng(7);
  std::uniform_real_distribution<float> uphi(-M_PI, M_PI);
  std::uniform_real_distribution<float> uw(0.001f, 0.20f);

  int bad_mkbins = 0, bad_v1 = 0, bad_helper = 0;
  for (int i = 0; i < 20000; ++i) {
    float c = uphi(rng), w = uw(rng);
    float lo = c - w, hi = c + w;

    // MkBins.cc: phiBinChecked(lo) .. phiBinChecked(hi)   -- NO +1
    unsigned short mb_b = ax.from_R_to_N_bin_safe(lo), mb_e = ax.from_R_to_N_bin_safe(hi);
    // MkFinder.cc:374 (V1): phiMaskApply(phiBin(hi) + 1)  -- mask AFTER add
    unsigned short v1_b = ax.from_R_to_N_bin_safe(lo);
    unsigned short v1_e = (unsigned short)((ax.from_R_to_N_bin(hi) + 1) & ax.c_N_mask);
    auto hp = ax.from_R_minmax_to_N_bins(lo, hi);

    for (int k = 0; k <= 32; ++k) {
      float x = lo + (hi - lo) * k / 32.0f;
      unsigned short t = ax.from_R_to_N_bin_safe(x);
      if (!range_holds<unsigned short>(mb_b, mb_e, t, kPhiBins, true)) { ++bad_mkbins; break; }
    }
    for (int k = 0; k <= 32; ++k) {
      float x = lo + (hi - lo) * k / 32.0f;
      unsigned short t = ax.from_R_to_N_bin_safe(x);
      if (!range_holds<unsigned short>(v1_b, v1_e, t, kPhiBins, true)) { ++bad_v1; break; }
    }
    for (int k = 0; k <= 32; ++k) {
      float x = lo + (hi - lo) * k / 32.0f;
      unsigned short t = ax.from_R_to_N_bin_safe(x);
      if (!range_holds<unsigned short>(hp.begin, hp.end, t, kPhiBins, true)) { ++bad_helper; break; }
    }
  }
  check(bad_v1 == 0,     "V1 idiom covers",              std::to_string(bad_v1) + " misses");
  check(bad_helper == 0, "axis helper covers",           std::to_string(bad_helper) + " misses");
  check_known(bad_mkbins == 0, "legacy MkBins idiom (no +1) covers",
              "drops the top bin on " + std::to_string(bad_mkbins) + " of 20000 ranges"
              " -- REPLACED, kept as the regression it was");
}

// 3b. The CURRENT MkBins form: the axis helper over [span +- cut tolerance],
//     then widened by whole BINS on the index. The contract it must satisfy is
//     stronger than plain coverage of the span -- it must fetch every hit the
//     per-hit cut could ACCEPT, i.e. anything within `tol` of any point of the
//     span. That is the invariant that makes "cut wider than fetch" impossible.
static void t_current_fetch_contract() {
  section("current MkBins fetch: covers everything the cut can accept");
  phi_axis_t ax(-M_PI, M_PI);
  const float bin_w = 2.0f * float(M_PI) / kPhiBins;
  std::mt19937 rng(31337);
  std::uniform_real_distribution<float> uphi(-M_PI, M_PI);
  std::uniform_real_distribution<float> uspan(0.0f, 0.30f);
  std::uniform_real_distribution<float> utol(0.0f, 0.05f);

  for (int extra_bins : {0, 1, 2}) {
    int bad = 0;
    for (int i = 0; i < 20000; ++i) {
      float c = uphi(rng), half = uspan(rng), tol = utol(rng);
      float lo = c - half, hi = c + half;

      auto pr = ax.from_R_minmax_to_N_bins(lo - tol, hi + tol);
      unsigned short b = (unsigned short)((pr.begin - extra_bins) & ax.c_N_mask);
      unsigned short e = (unsigned short)((pr.end   + extra_bins) & ax.c_N_mask);

      // every phi the cut could accept: within tol of some point of [lo, hi]
      for (int k = 0; k <= 64; ++k) {
        float x = (lo - tol) + ((hi + tol) - (lo - tol)) * k / 64.0f;
        unsigned short t = ax.from_R_to_N_bin_safe(x);
        if (!range_holds<unsigned short>(b, e, t, kPhiBins, true)) { ++bad; break; }
      }
    }
    check(bad == 0, "fetch covers the cut, extra_bins = " + std::to_string(extra_bins),
          std::to_string(bad) + " of 20000");
  }

  // And the extender really is in bin units: n extra bins must widen the range
  // by exactly 2n bins, with no rounding slop.
  int bad_width = 0;
  for (int i = 0; i < 5000; ++i) {
    float c = uphi(rng), half = uspan(rng);
    auto pr = ax.from_R_minmax_to_N_bins(c - half, c + half);
    unsigned base = range_size<unsigned short>(pr.begin, pr.end, kPhiBins, true);
    if (base == 0 || base + 4 > kPhiBins) continue;
    unsigned short b = (unsigned short)((pr.begin - 2) & ax.c_N_mask);
    unsigned short e = (unsigned short)((pr.end   + 2) & ax.c_N_mask);
    if (range_size<unsigned short>(b, e, kPhiBins, true) != base + 4) ++bad_width;
  }
  check(bad_width == 0, "bin extender widens by exactly 2n bins",
        std::to_string(bad_width) + " off");
}

// 4. Degenerate and full-circle ranges -- the two ends of the b == e ambiguity.
static void t_degenerate_and_full_circle() {
  section("degenerate and full-circle phi ranges");
  phi_axis_t ax(-M_PI, M_PI);

  // lo == hi must still yield exactly one bin.
  int bad_point = 0;
  for (int i = 0; i < 1000; ++i) {
    float x = -M_PI + 2.0f * M_PI * i / 1000.0f;
    auto p = ax.from_R_minmax_to_N_bins(x, x);
    if (range_size<unsigned short>(p.begin, p.end, kPhiBins, true) != 1) ++bad_point;
  }
  check(bad_point == 0, "point range gives exactly 1 bin", std::to_string(bad_point) + " bad");

  // COINCIDENT ENDPOINTS ARE A POINT, NOT A CIRCLE. On S^1, -pi and +pi are the
  // same point, so from_R_minmax_to_N_bins(-pi, +pi) asks for a range from a
  // point to itself and one bin is the CORRECT answer. The full circle is not
  // expressible here at all, and that is a domain boundary rather than a defect:
  // a range on a circle is an ARC, and the circle is not an arc -- two endpoints
  // carry no third datum saying "go all the way round".
  auto wrapped = ax.from_R_minmax_to_N_bins(-M_PI, M_PI);
  check(range_size<unsigned short>(wrapped.begin, wrapped.end, kPhiBins, true) == 1,
        "coincident endpoints (-pi, +pi) give one bin");

  // THE PRECONDITION THAT FOLLOWS, and it is the caller's. A half-width at or
  // above pi cannot be represented, and the failure is SILENT: the request wraps
  // and comes back as some arbitrary small arc rather than as everything. This
  // check documents the boundary; MkBins clamps before it is reached.
  int wide_undersized = 0;
  for (int i = 0; i < 1000; ++i) {
    float c = -M_PI + 2.0f * M_PI * i / 1000.0f;
    float w = 1.05f * M_PI;                       // deliberately over the limit
    auto p2 = ax.from_R_minmax_to_N_bins(c - w, c + w);
    if (range_size<unsigned short>(p2.begin, p2.end, kPhiBins, true) < kPhiBins / 2)
      ++wide_undersized;
  }
  check(wide_undersized > 0,
        "half-width >= pi is out of domain and fails SILENTLY (caller must clamp)",
        std::to_string(wide_undersized) + " of 1000 came back under half the circle");
}

// 5. The binnor itself: a region query must return exactly the brute-force set.
static void t_binnor_query_roundtrip() {
  section("binnor: region query against brute force");
  phi_axis_t ax_phi(-M_PI, M_PI);
  q_axis_t   ax_q(-120.0f, 120.0f, 1u << 8);
  binnor<unsigned int, phi_axis_t, q_axis_t, 18, 14> b(ax_phi, ax_q);

  std::mt19937 rng(4242);
  std::uniform_real_distribution<float> uphi(-M_PI, M_PI), uq(-119.0f, 119.0f);
  const int NP = 20000;
  std::vector<float> P(NP), Q(NP);

  b.begin_registration(NP);
  for (int i = 0; i < NP; ++i) { P[i] = uphi(rng); Q[i] = uq(rng); b.register_entry(P[i], Q[i]); }
  b.finalize_registration();

  int bad = 0;
  for (int t = 0; t < 200; ++t) {
    float pc = uphi(rng), pw = 0.05f + 0.20f * (t % 5);
    float qc = uq(rng),   qw = 2.0f + 8.0f * (t % 3);

    // brute force over the CLOSED window
    std::vector<int> want;
    for (int i = 0; i < NP; ++i) {
      float d = P[i] - pc;
      while (d >  M_PI) d -= 2 * M_PI;
      while (d < -M_PI) d += 2 * M_PI;
      if (std::abs(d) <= pw && std::abs(Q[i] - qc) <= qw) want.push_back(i);
    }

    // binnor: walk the bin ranges and keep what really falls inside
    auto pr = ax_phi.from_R_minmax_to_N_bins(pc - pw, pc + pw);
    auto qr = ax_q.from_R_minmax_to_N_bins(qc - qw, qc + qw);
    std::vector<int> got;
    for (unsigned short qi = qr.begin; qi != qr.end; ++qi) {
      for (unsigned short pi = pr.begin; pi != pr.end; pi = (unsigned short)((pi + 1) & ax_phi.c_N_mask)) {
        auto cbi = b.get_content(pi, qi);
        for (unsigned j = 0; j < cbi.count; ++j) {
          int idx = b.m_ranks[cbi.first + j];
          float d = P[idx] - pc;
          while (d >  M_PI) d -= 2 * M_PI;
          while (d < -M_PI) d += 2 * M_PI;
          if (std::abs(d) <= pw && std::abs(Q[idx] - qc) <= qw) got.push_back(idx);
        }
      }
    }
    std::sort(want.begin(), want.end());
    std::sort(got.begin(), got.end());
    if (want != got) ++bad;
  }
  check(bad == 0, "binnor query == brute force", std::to_string(bad) + " of 200 queries differ");
}

// 6. radix_sort: the OPCODE port. finalize_registration() uses it, so the query
//    round-trip already covers it end to end -- but a direct check localises a
//    failure instead of making it look like a binnor bug, and it can reach the
//    edge cases a random fill never produces.
template <typename V, typename R>
static void radix_case(const std::vector<V> &vals, const char *what) {
  std::vector<R> ranks;
  radix_sort<V, R> rs;
  rs.sort(vals, ranks);

  bool ok = ranks.size() == vals.size();
  // ranks must be a PERMUTATION ...
  if (ok) {
    std::vector<char> seen(vals.size(), 0);
    for (auto r : ranks) { if (r >= vals.size() || seen[r]) { ok = false; break; } seen[r] = 1; }
  }
  // ... that orders the values non-decreasingly, and is STABLE on ties.
  if (ok)
    for (size_t i = 1; i < ranks.size(); ++i)
      if (vals[ranks[i - 1]] > vals[ranks[i]] ||
          (vals[ranks[i - 1]] == vals[ranks[i]] && ranks[i - 1] > ranks[i])) { ok = false; break; }

  check(ok, std::string("radix_sort ") + what, "n = " + std::to_string(vals.size()));
}

static void t_radix_sort() {
  section("radix_sort (OPCODE port): permutation, order, stability");
  std::mt19937 rng(20260923);
  using u32 = unsigned int;
  using u16 = unsigned short;

  { std::vector<u32> v;                 radix_case<u32, u32>(v, "empty"); }
  { std::vector<u32> v{42};             radix_case<u32, u32>(v, "single"); }
  { std::vector<u32> v(1000, 7u);       radix_case<u32, u32>(v, "all equal (stability)"); }
  { std::vector<u32> v(1000); for (u32 i = 0; i < 1000; ++i) v[i] = i;
    radix_case<u32, u32>(v, "already sorted"); }
  { std::vector<u32> v(1000); for (u32 i = 0; i < 1000; ++i) v[i] = 999 - i;
    radix_case<u32, u32>(v, "reverse sorted"); }
  { std::vector<u32> v(5000); std::uniform_int_distribution<u32> d(0, 0xffffffffu);
    for (auto &x : v) x = d(rng); radix_case<u32, u32>(v, "random u32, full range"); }
  { // values differing only in the TOP byte, so the last radix pass decides
    std::vector<u32> v(256); for (u32 i = 0; i < 256; ++i) v[i] = (255u - i) << 24;
    radix_case<u32, u32>(v, "top-byte-only spread"); }
  { // and only in the LOW byte, so the first pass decides and later ones are no-ops
    std::vector<u32> v(256); for (u32 i = 0; i < 256; ++i) v[i] = 255u - i;
    radix_case<u32, u32>(v, "low-byte-only spread"); }
  { std::vector<u32> v(4000); std::uniform_int_distribution<u32> d(0, 3);
    for (auto &x : v) x = d(rng); radix_case<u32, u32>(v, "many ties (stability)"); }
  // the u16-RANK instantiation, which is what binnor uses for small fills
  { std::vector<u32> v(3000); std::uniform_int_distribution<u32> d(0, 0xffffffffu);
    for (auto &x : v) x = d(rng); radix_case<u32, u16>(v, "u16 ranks, random u32"); }
  { std::vector<u32> v(2000); std::uniform_int_distribution<u32> d(0, 3);
    for (auto &x : v) x = d(rng); radix_case<u32, u16>(v, "u16 ranks, many ties"); }

  // NOTE: radix_sort's static_assert admits V = unsigned short, but radix_sort.cc
  // instantiates only <u32,u32> and <u32,u16>, so a u16-VALUE sort is a link
  // error rather than a compile error. Either the assert is too permissive or an
  // instantiation is missing; not exercised here because it cannot be linked.
}

//-------------------------------------------------------------------- main

int main() {
  printf("binnor_test -- correctness suite for binnor.h\n");
  printf("phi axis: %u bins over 2pi, bin width %.6f rad, half-bin %.6f\n",
         kPhiBins, 2.0 * M_PI / kPhiBins, M_PI / kPhiBins);

  t_periodic_helper_covers();
  t_bounded_helper_covers();
  t_callsite_idioms();
  t_current_fetch_contract();
  t_degenerate_and_full_circle();
  t_binnor_query_roundtrip();
  t_radix_sort();

  printf("\n----------------------------------------------------------\n");
  printf("passed %d, FAILED %d, known-broken %d, newly-fixed %d\n",
         g_pass, g_fail, g_known, g_fixed);
  if (g_fixed) printf("NOTE: a known-broken check now passes -- promote it.\n");
  return g_fail == 0 ? 0 : 1;
}
