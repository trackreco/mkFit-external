// geocmp -- compare two TrackerInfo geometry binaries, module by module.
//
// Answers "is this geometry the same tracker as that one", which the geometry
// VERSION STAMP deliberately does not: the stamp records provenance, and two
// different names can carry identical geometry (D110, D121 and D127 do).
//
// Compare through the accessors, NOT with cmp. Until 2026-09-22 ModuleInfo wrote
// two uninitialised pad bytes per module, so two dumps of the SAME geometry
// differed in ~46 kB; any binary written before that fix still carries them.
//
// Build by hand -- test/ is not swept by the Makefile's $(wildcard ${SADIR}/*.cc):
//
//   . ~/root7.env
//   S=../../../../standalone      # the build dir holding libMkFitCore.so
//   g++ -std=c++20 -O2 -I../../../.. -I${ROOTSYS}/include geocmp.cc -o geocmp \
//       -L$S -lMkFitCore -lMkFitCMS -Wl,-rpath,$S \
//       -L${ROOTSYS}/lib -lCore -lMathCore -ltbb -Wl,-rpath,${ROOTSYS}/lib
//   LD_LIBRARY_PATH=$S:${ROOTSYS}/lib ./geocmp A.bin B.bin
//
// Compare two TrackerInfo binaries module by module, through the real accessors.
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace mkfit;

static double q(std::vector<double>& v, double f) {
  if (v.empty()) return 0;
  size_t i = std::min(v.size() - 1, (size_t)(f * v.size()));
  std::nth_element(v.begin(), v.begin() + i, v.end());
  return v[i];
}

int main(int argc, char** argv) {
  if (argc < 3) { printf("usage: geocmp A.bin B.bin\n"); return 1; }
  TrackerInfo a, b;
  a.read_bin_file(argv[1]);
  b.read_bin_file(argv[2]);
  printf("A geom '%s'  layers %d\n", a.geom_version().c_str(), a.n_layers());
  printf("B geom '%s'  layers %d\n", b.geom_version().c_str(), b.n_layers());
  if (a.n_layers() != b.n_layers()) { printf("LAYER COUNT DIFFERS\n"); return 2; }

  long n_mod = 0, n_moved = 0, n_detid = 0, n_rot = 0;
  std::vector<double> dpos, drot;
  double maxpos = 0, maxrot = 0;
  int worst_lay = -1;

  for (int l = 0; l < a.n_layers(); ++l) {
    const LayerInfo& la = a[l];
    const LayerInfo& lb = b[l];
    if (la.n_modules() != lb.n_modules()) {
      printf("layer %2d: module count %d vs %d\n", l, la.n_modules(), lb.n_modules());
      continue;
    }
    for (int m = 0; m < la.n_modules(); ++m) {
      const ModuleInfo& ma = la.module_info(m);
      const ModuleInfo& mb = lb.module_info(m);
      ++n_mod;
      if (ma.detid != mb.detid) ++n_detid;
      double d = std::sqrt((ma.pos[0]-mb.pos[0])*(ma.pos[0]-mb.pos[0]) +
                           (ma.pos[1]-mb.pos[1])*(ma.pos[1]-mb.pos[1]) +
                           (ma.pos[2]-mb.pos[2])*(ma.pos[2]-mb.pos[2]));
      double dz = std::sqrt((ma.zdir[0]-mb.zdir[0])*(ma.zdir[0]-mb.zdir[0]) +
                            (ma.zdir[1]-mb.zdir[1])*(ma.zdir[1]-mb.zdir[1]) +
                            (ma.zdir[2]-mb.zdir[2])*(ma.zdir[2]-mb.zdir[2]));
      if (d > 0) { ++n_moved; dpos.push_back(d); if (d > maxpos) { maxpos = d; worst_lay = l; } }
      if (dz > 0) { ++n_rot; drot.push_back(dz); maxrot = std::max(maxrot, dz); }
    }
  }
  printf("\nmodules compared        : %ld\n", n_mod);
  printf("detid differs           : %ld\n", n_detid);
  printf("position differs        : %ld  (%.2f %%)\n", n_moved, 100.0*n_moved/std::max(1L,n_mod));
  printf("  |dpos| um  p50 %.4f  p90 %.4f  p99 %.4f  max %.4f   (worst layer %d)\n",
         1e4*q(dpos,0.50), 1e4*q(dpos,0.90), 1e4*q(dpos,0.99), 1e4*maxpos, worst_lay);
  printf("normal differs          : %ld  (%.2f %%)\n", n_rot, 100.0*n_rot/std::max(1L,n_mod));
  printf("  |dzdir| mrad p50 %.4f  p90 %.4f  max %.4f\n",
         1e3*q(drot,0.50), 1e3*q(drot,0.90), 1e3*maxrot);

  // ---- material grid, which the binary also carries ----
  printf("\nmaterial grid A: %d x %d bins, range z %.1f r %.1f\n",
         a.mat_nbins_z(), a.mat_nbins_r(), a.mat_range_z(), a.mat_range_r());
  printf("material grid B: %d x %d bins, range z %.1f r %.1f\n",
         b.mat_nbins_z(), b.mat_nbins_r(), b.mat_range_z(), b.mat_range_r());
  if (a.mat_nbins_z() == b.mat_nbins_z() && a.mat_nbins_r() == b.mat_nbins_r()) {
    long nb = 0, nd = 0, nz_a = 0, nz_b = 0;
    std::vector<double> rel;
    double maxrel = 0; int mz = -1, mr = -1;
    for (int iz = 0; iz < a.mat_nbins_z(); ++iz)
      for (int ir = 0; ir < a.mat_nbins_r(); ++ir) {
        double ra = a.material_radl(iz, ir), rb = b.material_radl(iz, ir);
        ++nb;
        if (ra > 0) ++nz_a;
        if (rb > 0) ++nz_b;
        if (ra != rb) {
          ++nd;
          double den = std::max(std::abs(ra), std::abs(rb));
          double r = den > 0 ? std::abs(ra - rb) / den : 0;
          rel.push_back(r);
          if (r > maxrel) { maxrel = r; mz = iz; mr = ir; }
        }
      }
    printf("radl bins             : %ld   non-empty A %ld  B %ld\n", nb, nz_a, nz_b);
    printf("radl bins differing   : %ld  (%.2f %% of all, %.2f %% of non-empty A)\n",
           nd, 100.0*nd/nb, 100.0*nd/std::max(1L,nz_a));
    printf("  relative diff  p50 %.4f  p90 %.4f  p99 %.4f  max %.4f  at (z %d, r %d)\n",
           q(rel,0.50), q(rel,0.90), q(rel,0.99), maxrel, mz, mr);
  }
  return 0;
}
