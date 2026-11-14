#include "RecoTracker/MkFitCore/src/MiniPropagators.h"
#include <cstdio>
#include <cmath>

// Independent reference: RK4 integration of dr/ds = phat, dphat/ds = kappa*(phat_y,-phat_x,0).
// kappa = sign(q) * 0.01*sol*B / |p|  [1/cm].  Shares no algebra with propagate_to_r()
// nor with the code's closed-form helix -- it only knows F = q v x B.
namespace {
  struct RkState { double x,y,z, ux,uy,uz; };

  RkState rk_deriv(const RkState &s, double kap) {
    return { s.ux, s.uy, s.uz, kap*s.uy, -kap*s.ux, 0.0 };
  }
  void rk_step(RkState &s, double h, double kap) {
    auto add=[&](const RkState&a,const RkState&d,double f){
      return RkState{a.x+f*d.x,a.y+f*d.y,a.z+f*d.z,a.ux+f*d.ux,a.uy+f*d.uy,a.uz+f*d.uz}; };
    RkState k1=rk_deriv(s,kap);
    RkState k2=rk_deriv(add(s,k1,h/2),kap);
    RkState k3=rk_deriv(add(s,k2,h/2),kap);
    RkState k4=rk_deriv(add(s,k3,h),kap);
    s.x  += h/6*(k1.x +2*k2.x +2*k3.x +k4.x);
    s.y  += h/6*(k1.y +2*k2.y +2*k3.y +k4.y);
    s.z  += h/6*(k1.z +2*k2.z +2*k3.z +k4.z);
    s.ux += h/6*(k1.ux+2*k2.ux+2*k3.ux+k4.ux);
    s.uy += h/6*(k1.uy+2*k2.uy+2*k3.uy+k4.uy);
    s.uz += h/6*(k1.uz+2*k2.uz+2*k3.uz+k4.uz);
  }
  // integrate until hypot(x,y) crosses R; return path length and final state
  bool rk_to_r(RkState s, double kap, double R, double &s_out, RkState &out) {
    const double h = 1e-4;               // cm
    double r0 = std::hypot(s.x,s.y), sl = 0;
    for (int i = 0; i < 20000000; ++i) {
      RkState prev = s; double rp = std::hypot(s.x,s.y);
      rk_step(s, h, kap); sl += h;
      double r = std::hypot(s.x,s.y);
      if ((rp-R)*(r-R) <= 0 && i > 0) {   // bracketed -- bisect this step
        double lo=0, hi=h;
        for (int j=0;j<60;++j){ double m=0.5*(lo+hi); RkState t=prev; rk_step(t,m,kap);
          if ((rp-R)*(std::hypot(t.x,t.y)-R) <= 0) hi=m; else lo=m; }
        RkState t=prev; rk_step(t,0.5*(lo+hi),kap);
        out=t; s_out = sl - h + 0.5*(lo+hi); return true;
      }
      (void)r0;
    }
    return false;
  }
}

void cc_unit() {
  namespace mp = mkfit::mini_propagators;
  const double B = mkfit::Config::Bfield;
  const double kfac = 0.01 * mkfit::Const::sol * B;
  printf("\nUNIT  Bfield=%.3f  0.01*sol*B=%.6f  (R_c[cm] = pT/%.6f)\n", B, kfac, kfac);
  printf("UNIT  %-34s %12s %12s %12s %10s\n", "case", "|dpos| cm", "|dmom| GeV", "dalpha", "fail");

  struct Case { double x,y,z, pT,phi,theta; int q; double R; const char *nm; };
  Case cases[] = {
    {  3.0,  0.0,  0.0,  1.0, 0.30, 1.20, +1,  20.0, "q+ pT=1  outward 3->20" },
    {  3.0,  0.0,  0.0,  1.0, 0.30, 1.20, -1,  20.0, "q- pT=1  outward 3->20" },
    { 25.0, 10.0, 30.0, 10.0, 1.10, 0.80, +1,  60.0, "q+ pT=10 outward 27->60" },
    { 60.0,-20.0,-40.0,  2.0,-2.00, 2.10, -1,  30.0, "q- pT=2  inward 63->30" },
    { 15.0,  4.0,  5.0,  0.8, 0.90, 1.57, +1,  40.0, "q+ pT=0.8 outward, barrel" },
    { 90.0, 30.0, 10.0,  1.5, 2.60, 1.40, -1,  40.0, "q- pT=1.5 inward 95->40" },
  };
  for (auto &c : cases) {
    double px = c.pT*std::cos(c.phi), py = c.pT*std::sin(c.phi), pz = c.pT/std::tan(c.theta);
    double pmag = std::sqrt(px*px+py*py+pz*pz);
    // --- mkfit closed form
    mp::State s0(c.x,c.y,c.z, px,py,pz, 0.f, 0);
    mp::InitialState is(s0, (short)c.q, 1.0f/c.pT, c.theta);
    mp::State got;
    is.propagate_to_r(mp::PA_Exact, c.R, got, true);
    // --- RK4 reference
    double kap = (c.q > 0 ? +1.0 : -1.0) * kfac / pmag;
    RkState r0{c.x,c.y,c.z, px/pmag,py/pmag,pz/pmag}, rf; double sl=0;
    if (!rk_to_r(r0, kap, c.R, sl, rf)) {
      // RK4 never reaches R -- the target is genuinely unreachable, so this tests the
      // clamp instead. Analytic reachable band from the circle geometry:
      double kk = (c.q > 0 ? -1.0 : +1.0) / kfac;      // sign(k) = -sign(q)
      double cx = c.x - kk*py, cy = c.y + kk*px;
      double dc = std::hypot(cx,cy), Rc = std::fabs(kk)*c.pT;
      double rlo = std::fabs(dc-Rc), rhi = dc+Rc;
      double rgot = std::hypot(got.x,got.y);
      double want = (c.R < rlo) ? rlo + mkfit::mini_propagators::kReachMargin
                                : rhi - mkfit::mini_propagators::kReachMargin;
      printf("UNIT  %-34s UNREACHABLE band=[%.3f,%.3f] R=%.1f -> clamped r=%.4f "
             "(expect %.4f, d=%.2e)  fail=%d %s\n",
             c.nm, rlo, rhi, c.R, rgot, want, std::fabs(rgot-want), got.fail_flag,
             (got.fail_flag==2 && std::fabs(rgot-want)<1e-3) ? "OK" : "** CHECK **");
      continue;
    }
    double dpos = std::sqrt(std::pow(got.x-rf.x,2)+std::pow(got.y-rf.y,2)+std::pow(got.z-rf.z,2));
    double dmom = std::sqrt(std::pow(got.px-rf.ux*pmag,2)+std::pow(got.py-rf.uy*pmag,2)+std::pow(got.pz-rf.uz*pmag,2));
    // alpha from RK4: transverse momentum azimuth rotates by +alpha in the code's convention
    double a_rk = std::remainder(std::atan2(rf.uy,rf.ux) - std::atan2(py,px), 2*M_PI);
    printf("UNIT  %-34s %12.3e %12.3e %12.3e %10d   (a_mkfit=%+.6f a_rk4=%+.6f)\n",
           c.nm, dpos, dmom, got.dalpha - a_rk, got.fail_flag, got.dalpha, a_rk);
  }
}

// Test B: the one-point (extrapolating) Hermite against RK4.
// Note it needs nothing but a state, inv_k and a chosen dalpha -- no prior two-point
// call, and the state can be any position+momentum pair you like.
void cc_unit_h3() {
  namespace mp = mkfit::mini_propagators;
  const double B = mkfit::Config::Bfield;
  const double kfac = 0.01 * mkfit::Const::sol * B;
  printf("\nUNIT-H3  one-point Hermite vs RK4, |dpos| in cm\n");
  printf("UNIT-H3  %-26s %10s %10s %10s %10s %10s\n",
         "case (dalpha span)", "t=0", "t=0.25", "t=0.50", "t=0.75", "t=1.00");
  struct Case { double x,y,z,pT,phi,theta; int q; double dal; const char *nm; };
  Case cases[] = {
    { 20.0, 5.0, 10.0, 1.0, 0.30, 1.20, +1, 0.05, "pT=1   dal=0.05" },
    { 20.0, 5.0, 10.0, 1.0, 0.30, 1.20, +1, 0.20, "pT=1   dal=0.20" },
    { 20.0, 5.0, 10.0, 1.0, 0.30, 1.20, -1, 0.20, "pT=1 q- dal=0.20" },
    { 60.0,-20.0, 30.0,10.0, 1.10, 0.80, +1, 0.05, "pT=10  dal=0.05" },
    { 8.0,  2.0,  1.0, 0.5,-0.70, 1.50, +1, 0.40, "pT=0.5 dal=0.40" },
  };
  for (auto &c : cases) {
    double px = c.pT*std::cos(c.phi), py = c.pT*std::sin(c.phi), pz = c.pT/std::tan(c.theta);
    double pmag = std::sqrt(px*px+py*py+pz*pz);
    mp::State s0(c.x,c.y,c.z, px,py,pz, 0.f, 0);
    mp::InitialState is(s0, (short)c.q, 1.0f/c.pT, c.theta);
    // Build the cubic from that single state.
    mp::StatePlex sp; sp.copyIn(0, (const mp::State &) is);
    mp::MPF ik(0.0f), dal(0.0f);
    ik[0] = is.inv_k; dal[0] = (float) c.dal;
    mp::Hermite3D H;
    H.calculate_coeffs(sp, ik, dal);
    printf("UNIT-H3  %-26s", c.nm);
    double kap = (c.q > 0 ? +1.0 : -1.0) * kfac / pmag;
    for (double t : {0.0, 0.25, 0.50, 0.75, 1.00}) {
      mp::MPF tp(0.0f), hx, hy, hz; tp[0] = (float) t;
      H.evaluate(tp, hx, hy, hz);
      double alpha = c.dal * t;
      // SIGNED path length. s_3D = k * |p| * alpha with k = 1/inv_k, and sign(k) =
      // -sign(q): for a positive charge, moving FORWARD makes alpha NEGATIVE. Using
      // |alpha| here compares the cubic at -alpha against RK4 at +s and "fails" for
      // exactly one charge sign -- which is what it did on the first attempt.
      double sl = alpha * pmag / is.inv_k;
      RkState r{c.x,c.y,c.z, px/pmag,py/pmag,pz/pmag};
      int n = 4000; double h = sl/n;
      for (int i=0;i<n;++i) rk_step(r, h, kap);
      double d = std::sqrt(std::pow(hx[0]-r.x,2)+std::pow(hy[0]-r.y,2)+std::pow(hz[0]-r.z,2));
      printf(" %10.2e", d);
    }
    printf("   (s_span=%+.2f cm)\n", c.dal*pmag/is.inv_k);
  }
}
