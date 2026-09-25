#ifndef mkfit_seeding_SeedMath_h
#define mkfit_seeding_SeedMath_h

// The per-triplet and per-quad arithmetic, templated on an arithmetic policy.
//
//   ArithRef   double, std::hypot / asin / atan2.  Bit-identical to the
//              prototype's expressions (and to detail:: in SeedFinder.h), so
//              it reproduces the reference quad list exactly.
//   ArithFast  float, sqrt for hypot, vdt::fast_asinf / fast_atan2f, and the
//              same circle-CENTRE form as ArithRef.  Kept as the record of why
//              that form cannot go to float: for a stiff track the centre is
//              ~R away and dC^2 - R^2 cancels.  Measured, 50 events: margins
//              move by up to 1.5 mrad at pT ~2 TeV, 3 quads differ.
//   ArithFastK float and vdt as ArithFast, with the circle given by its signed
//              CURVATURE k and the point and tangent at the third hit.  The
//              crossing with |X| = r is the line X.(k P0 + n) = k (r^2 + r0^2)/2
//              + P0.n -- the circle equation multiplied through by k -- whose
//              terms all stay O(r) as k -> 0.
//
// finish_triplets<A>() and the margin evaluator (SeedMargins.h) both call the
// functions below, so the evaluator cannot drift from the finder: a quad the
// finder accepts in arithmetic A passes every cut of eval_quad<A>().
//
// Every cut also reports a signed MARGIN, positive on the passing side, in the
// cut's own unit (cm or rad), so a quad found in one arithmetic and not in the
// other can be traced to the cut that flipped and to how close it was.

#include "SeedLayer.h"

#include "vdt/asin.h"
#include "vdt/atan2.h"

#include <cmath>

namespace mkfit::seeding {

  namespace vdtv {
    // vdt::fast_asinf with the same polynomial and the same selects, but the
    // sign restored with copysign instead of through vdt's ieee754 union.
    // OR-ing the sign bit onto a non-negative result IS copysign, so the
    // values are identical; the union's store and reload is what stops GCC
    // from vectorising a loop that calls it.
    inline float fast_asinf(float x) {
      const float a = std::abs(x);
      const bool big = a > 0.5f;
      const float z0 = big ? 0.5f * (1.0f - a) : a * a;
      const float xx = big ? std::sqrt(z0) : a;
      const float z = ((((4.2163199048E-2f * z0 + 2.4181311049E-2f) * z0 + 4.5470025998E-2f) * z0 +
                        7.4953002686E-2f) * z0 + 1.6666752422E-1f) * z0 * xx + xx;
      const float tmp = 1.57079632679489661923f - (z + z);  // vdt::details::PIO2F
      const float res = a < 1e-4f ? a : (big ? tmp : z);
      return std::copysign(res, x);
    }
  }  // namespace vdtv

  struct ArithRef {
    using real = double;
    static constexpr const char *name = "ref";
    static double hypot(double a, double b) { return std::hypot(a, b); }
    static double sqrt(double a) { return std::sqrt(a); }
    static double asin(double a) { return std::asin(a); }
    static double atan2(double y, double x) { return std::atan2(y, x); }
    static bool finite(double a) { return std::isfinite(a); }
    static constexpr bool curvature_form = false;
  };

  struct ArithFast {
    using real = float;
    static constexpr const char *name = "fast";
    static float hypot(float a, float b) { return std::sqrt(a * a + b * b); }
    static float sqrt(float a) { return std::sqrt(a); }
    static float asin(float a) { return vdtv::fast_asinf(a); }
    static float atan2(float y, float x) { return vdt::fast_atan2f(y, x); }
    // -Ofast implies -ffinite-math-only, so std::isfinite is folded to true;
    // test the exponent bits instead
    static bool finite(float a) {
      unsigned int u;
      __builtin_memcpy(&u, &a, 4);
      return (u & 0x7f800000u) != 0x7f800000u;
    }
    static constexpr bool curvature_form = false;
  };

  struct ArithFastK : ArithFast {
    static constexpr const char *name = "fastk";
    static constexpr bool curvature_form = true;
  };

  namespace amath {
    constexpr float kTwoPi = 2.0f * kPi;
    constexpr float kCoverEps = 1e-4f;  // cm / rad, as the prototype's cover mode

    inline float wrap_pi(float d) {
      d -= kTwoPi * (d > kPi);
      d += kTwoPi * (d < -kPi);
      return d;
    }

    // The coordinates come in as float (the stored hit x, y) in both policies.
    // g_out: the determinant the degeneracy test is made on.
    template <class A>
    inline bool circle3(float x1,
                        float y1,
                        float x2,
                        float y2,
                        float x3,
                        float y3,
                        typename A::real &cx,
                        typename A::real &cy,
                        typename A::real &R,
                        typename A::real *g_out = nullptr) {
      using T = typename A::real;
      T A_ = x2 - x1, B = y2 - y1, C = x3 - x1, D = y3 - y1;
      T E = A_ * (x1 + x2) + B * (y1 + y2);
      T F = C * (x1 + x3) + D * (y1 + y3);
      T G = T(2.0) * (A_ * (y3 - y2) - B * (x3 - x2));
      if (g_out)
        *g_out = G;
      if (std::abs(G) < T(1e-12))
        return false;
      cx = (D * E - B * F) / G;
      cy = (A_ * F - C * E) / G;
      R = A::hypot(x1 - cx, y1 - cy);
      return R > T(1e-6) && A::finite(R);
    }

    template <class A>
    inline typename A::real arc(typename A::real chord, typename A::real R) {
      using T = typename A::real;
      T h = chord / (T(2.0) * R);
      if (h > T(1.0))
        h = T(1.0);
      return T(2.0) * R * A::asin(h);
    }

    // b2_out: rt^2 - aa^2, the quantity whose sign decides the crossing.
    template <class A>
    inline bool circle_cross_r(typename A::real cx,
                               typename A::real cy,
                               typename A::real R,
                               typename A::real rt,
                               typename A::real xref,
                               typename A::real yref,
                               typename A::real &px,
                               typename A::real &py,
                               typename A::real *b2_out = nullptr) {
      using T = typename A::real;
      const T dC = A::hypot(cx, cy);
      if (b2_out)
        *b2_out = -1;
      if (dC < T(1e-9))
        return false;
      const T aa = (rt * rt + dC * dC - R * R) / (2 * dC);
      const T b2 = rt * rt - aa * aa;
      if (b2_out)
        *b2_out = b2;
      if (b2 < 0)
        return false;
      const T bb = A::sqrt(b2);
      const T ux = cx / dC, uy = cy / dC;
      const T vx = -uy, vy = ux;
      const T x1 = aa * ux + bb * vx, y1 = aa * uy + bb * vy;
      const T x2 = aa * ux - bb * vx, y2 = aa * uy - bb * vy;
      if (A::hypot(x1 - xref, y1 - yref) <= A::hypot(x2 - xref, y2 - yref)) {
        px = x1;
        py = y1;
      } else {
        px = x2;
        py = y2;
      }
      return true;
    }

    // Arc length of a chord L on a circle of curvature k: 2 asin(|k| L / 2) / |k|,
    // written as L asin(h)/h so that it is exact at k = 0.
    template <class A>
    inline typename A::real arc_k(typename A::real L, typename A::real k) {
      using T = typename A::real;
      T h = std::abs(k) * L * T(0.5);
      if (h > T(1.0))
        h = T(1.0);
      return h < T(1e-3) ? L * (T(1) + h * h * T(1.0 / 6)) : L * A::asin(h) / h;
    }

    // Crossing of |X| = rt with the circle through (x0, y0), unit tangent
    // (ux, uy) there, signed curvature k (> 0 turning left).  The crossing
    // closest to (x0, y0) is taken, as circle_cross_r() does.
    template <class A>
    inline bool curv_cross_r(typename A::real k,
                             typename A::real ux,
                             typename A::real uy,
                             typename A::real x0,
                             typename A::real y0,
                             typename A::real rt,
                             typename A::real &px,
                             typename A::real &py,
                             typename A::real *b2_out = nullptr) {
      using T = typename A::real;
      const T nx = -uy, ny = ux;
      const T mx = k * x0 + nx, my = k * y0 + ny;  // k * (centre), finite at k = 0
      const T m2 = mx * mx + my * my;
      const T c = T(0.5) * k * (rt * rt + x0 * x0 + y0 * y0) + (x0 * nx + y0 * ny);
      const T im = T(1) / A::sqrt(m2);
      const T cn = c * im;  // distance of the radical line from the origin
      const T b2 = rt * rt - cn * cn;
      if (b2_out)
        *b2_out = b2;
      // branch-free, so that it vectorises: without a crossing, (px, py) is
      // the finite foot point and the return value says so
      const T hh = A::sqrt(std::max(b2, T(0)));
      const T fx = mx * (cn * im), fy = my * (cn * im);  // foot point
      const T wx = -my * im, wy = mx * im;               // along the line
      const T sg = ((fx - x0) * wx + (fy - y0) * wy) > 0 ? T(-1) : T(1);
      px = fx + sg * hh * wx;
      py = fy + sg * hh * wy;
      return b2 >= 0;
    }
  }  // namespace amath

  // A cut value: m is the signed margin (>= 0 passes, except the strict r
  // ordering cuts, where > 0 passes); pass is the finder's own decision.
  struct CutVal {
    double m = 0;
    bool pass = false;
    bool valid = false;  // false: not evaluated, an earlier stage made it undefined
  };

  // The helix of a triplet and the 4th-layer window.
  template <class A>
  struct TripletHelix {
    using T = typename A::real;
    T cx = 0, cy = 0, R = 0, cots = 0, xc = 0, yc = 0;
    T k = 0, ux = 0, uy = 0;  // curvature form
    T phid = 0, phidh = 0, zlo = 0, zhi = 0;
    bool ok = false;
    CutVal c3, reach;  // circle degeneracy, the 4th layer reachable at rdlo or rdhi
  };

  // The curvature-form helix of a triplet and the 4th-layer window, written
  // without branches so that finish_triplets() can run it in a simd loop
  // across triplets.  helix_triplet<A>() calls the same function, so the
  // margin evaluator sees the finder's arithmetic.
  template <class A>
  struct HelixK {
    using T = typename A::real;
    T k, ux, uy, cots, phid, phidh, zlo, zhi, b20, b21;
    bool ok;
  };

  // flatten: inline everything it calls, vdt included, so that the simd
  // loop in finish_triplets() sees no calls
  template <class A>
  [[gnu::always_inline, gnu::flatten]] inline void helix_k(const float qwin_d,
                      const float phiwin_d,
                      const float xa,
                      const float ya,
                      const float za,
                      const float xb,
                      const float yb,
                      const float xc,
                      const float yc,
                      const float zc,
                      const typename A::real rdlo,
                      const typename A::real rdhi,
                      HelixK<A> &o) {
    using namespace amath;
    using T = typename A::real;
    const T dx1 = xb - xa, dy1 = yb - ya, dx2 = xc - xb, dy2 = yc - yb, dx3 = xc - xa, dy3 = yc - ya;
    const T L1 = A::hypot(dx1, dy1), L2 = A::hypot(dx2, dy2), L3 = A::hypot(dx3, dy3);
    const T cr = dx1 * dy2 - dy1 * dx2;
    const T k = T(2) * cr / (L1 * L2 * L3);  // signed Menger curvature
    // tangent at c: the chord b->c turned by half its turn angle
    const T sa = T(0.5) * k * L2, ca = A::sqrt(std::max(T(0), T(1) - sa * sa));
    const T ex = dx2 / L2, ey = dy2 / L2;
    const T ux = ex * ca - ey * sa, uy = ex * sa + ey * ca;
    const T s_ac = arc_k<A>(L1, k) + arc_k<A>(L2, k);
    const T cots = (s_ac > T(1e-6)) ? T(zc - za) / s_ac : T(0.0);
    // locals, not &o.b20: taking a member's address keeps o in memory and
    // turns its stores into scatters, which stops the simd loop vectorising
    T p0x, p0y, p1x, p1y, b20, b21;
    const bool ok0 = curv_cross_r<A>(k, ux, uy, xc, yc, rdlo, p0x, p0y, &b20);
    const bool ok1 = curv_cross_r<A>(k, ux, uy, xc, yc, rdhi, p1x, p1y, &b21);
    o.b20 = b20;
    o.b21 = b21;
    const T a0 = A::atan2(p0y, p0x), a1 = A::atan2(p1y, p1x);
    const T f0 = ok0 ? a0 : a1;
    const T f1 = ok1 ? a1 : f0;
    const T dfh = T(0.5) * wrap_pi((float)(f1 - f0));
    const T arc0 = arc_k<A>(A::hypot(p0x - T(xc), p0y - T(yc)), k);
    const T arc1 = arc_k<A>(A::hypot(p1x - T(xc), p1y - T(yc)), k);
    const T zz0 = ok0 ? T(zc) + cots * arc0 : T(zc) + cots * arc1;
    const T zz1 = ok1 ? T(zc) + cots * arc1 : zz0;
    const T zdh = qwin_d + kCoverEps;
    o.k = k;
    o.ux = ux;
    o.uy = uy;
    o.cots = cots;
    o.ok = ok0 | ok1;
    o.phid = f0 + dfh;
    o.phidh = T(phiwin_d) + std::abs(dfh) + T(kCoverEps);
    o.zlo = std::min(zz0, zz1) - zdh;
    o.zhi = std::max(zz0, zz1) + zdh;
  }

  template <class A>
  inline void helix_triplet(const float qwin_d,
                            const float phiwin_d,
                            const float xa,
                            const float ya,
                            const float za,
                            const float xb,
                            const float yb,
                            const float xc,
                            const float yc,
                            const float zc,
                            const typename A::real rdlo,
                            const typename A::real rdhi,
                            TripletHelix<A> &h) {
    using namespace amath;
    using T = typename A::real;
    if constexpr (A::curvature_form) {
      HelixK<A> o;
      helix_k<A>(qwin_d, phiwin_d, xa, ya, za, xb, yb, xc, yc, zc, rdlo, rdhi, o);
      h.k = o.k;
      h.ux = o.ux;
      h.uy = o.uy;
      h.R = T(1) / std::abs(o.k);
      h.c3 = {1e30, true, true};  // no degeneracy: k = 0 is a straight line
      h.cots = o.cots;
      h.xc = xc;
      h.yc = yc;
      h.ok = o.ok;
      h.reach = {std::max(double(o.b20) / (2 * double(rdlo)), double(o.b21) / (2 * double(rdhi))), h.ok, true};
      if (!h.ok)
        return;
      h.phid = o.phid;
      h.phidh = o.phidh;
      h.zlo = o.zlo;
      h.zhi = o.zhi;
      return;
    }
    T cx = 0, cy = 0, R, G = 0;
    const bool c3ok = circle3<A>(xa, ya, xb, yb, xc, yc, cx, cy, R, &G);
    if (!c3ok)
      R = T(1e6);
    h.c3 = {std::abs(G) - 1e-12, c3ok, true};
    const T s_ab = arc<A>(A::hypot(T(xb) - T(xa), T(yb) - T(ya)), R);
    const T s_bc = arc<A>(A::hypot(T(xc) - T(xb), T(yc) - T(yb)), R);
    const T s_ac = s_ab + s_bc;
    const T cot_s = (s_ac > T(1e-6)) ? T(zc - za) / s_ac : T(0.0);
    T p0x, p0y, p1x, p1y, b20, b21;
    const bool ok0 = circle_cross_r<A>(cx, cy, R, rdlo, xc, yc, p0x, p0y, &b20);
    const bool ok1 = circle_cross_r<A>(cx, cy, R, rdhi, xc, yc, p1x, p1y, &b21);
    h.ok = ok0 || ok1;
    // b2 / 2r ~ r - |aa| near the edge: a length in cm
    h.reach = {std::max(double(b20) / (2 * double(rdlo)), double(b21) / (2 * double(rdhi))), h.ok, true};
    h.cx = cx;
    h.cy = cy;
    h.R = R;
    h.cots = cot_s;
    h.xc = xc;
    h.yc = yc;
    if (!h.ok)
      return;
    const T f0 = ok0 ? A::atan2(p0y, p0x) : A::atan2(p1y, p1x);
    const T f1 = ok1 ? A::atan2(p1y, p1x) : f0;
    const T dfh = T(0.5) * wrap_pi((float)(f1 - f0));
    const T zz0 = ok0 ? T(zc) + cot_s * arc<A>(A::hypot(p0x - T(xc), p0y - T(yc)), R)
                      : T(zc) + cot_s * arc<A>(A::hypot(p1x - T(xc), p1y - T(yc)), R);
    const T zz1 = ok1 ? T(zc) + cot_s * arc<A>(A::hypot(p1x - T(xc), p1y - T(yc)), R) : zz0;
    const T zdh = qwin_d + kCoverEps;
    h.phid = f0 + dfh;
    h.phidh = T(phiwin_d) + std::abs(dfh) + T(kCoverEps);
    h.zlo = std::min(zz0, zz1) - zdh;
    h.zhi = std::max(zz0, zz1) + zdh;
  }

  // The 4th-layer per-hit cuts of one (triplet, d-hit) pair.
  struct QuadCuts {
    CutVal cd_r, d_cross, d_phi, d_z;
    bool pass() const { return cd_r.pass && d_cross.pass && d_phi.pass && d_z.pass; }
  };

  // full = false: stop at the first failing cut, as the finder does.
  template <class A>
  inline bool quad_cuts(const float qwin_d,
                        const float phiwin_d,
                        const TripletHelix<A> &h,
                        const float zc,
                        const float rrc,
                        const float rrd,
                        const float phi_d,
                        const float z_d,
                        QuadCuts *qc) {
    using namespace amath;
    using T = typename A::real;
    const bool r_ok = !(rrd <= rrc + 0.1f);
    if (qc)
      qc->cd_r = {double(rrd - (rrc + 0.1f)), r_ok, true};
    else if (!r_ok)
      return false;
    T px2, py2, b2;
    bool x_ok;
    if constexpr (A::curvature_form)
      x_ok = curv_cross_r<A>(h.k, h.ux, h.uy, h.xc, h.yc, T(rrd), px2, py2, &b2);
    else
      x_ok = circle_cross_r<A>(h.cx, h.cy, h.R, T(rrd), h.xc, h.yc, px2, py2, &b2);
    if (qc)
      qc->d_cross = {double(b2) / (2 * double(rrd)), x_ok, true};
    if (!x_ok)
      return false;
    const T dphi_d = wrap_pi((float)(phi_d - A::atan2(py2, px2)));
    T zd2;
    if constexpr (A::curvature_form)
      zd2 = T(zc) + h.cots * arc_k<A>(A::hypot(px2 - h.xc, py2 - h.yc), h.k);
    else
      zd2 = T(zc) + h.cots * arc<A>(A::hypot(px2 - h.xc, py2 - h.yc), h.R);
    const T dz_d = z_d - zd2;
    const bool p_ok = !(std::abs(dphi_d) > phiwin_d), z_ok = !(std::abs(dz_d) > qwin_d);
    if (qc) {
      qc->d_phi = {double(phiwin_d) - double(std::abs(dphi_d)), p_ok, true};
      qc->d_z = {double(qwin_d) - double(std::abs(dz_d)), z_ok, true};
    }
    return r_ok && p_ok && z_ok;
  }

}  // namespace mkfit::seeding

#endif
