// ACLiC loader for the pure-geometry seeding study.
//
//   from the standalone BUILD directory, inside an mkFit --shell:
//     .x <mkFit-external checkout>/mkfit-standalone-seedgeom/seedgeom-load.C
//   The study is compiled from the directory this file sits in, so any
//   checkout works -- including the mkfit-seeding worktree next to the build's
//   own mkFit-external.
//   then sg_reset(), sg_ev(s.event()), sg_report(), ...
//
// Nothing of this study lives in the CMSSW source tree; it is compiled on
// demand by ACLiC, which also generates the dictionary that makes the entry
// points callable from cling.
//
// THE FLAGS ARE NOT HARDCODED.  They come from `make echo-aclic`, because a
// stale copy here is silent corruption, not a build error: standalone/Event.h
// gains members under MKFIT_TRACE, so an ACLiC build compiled without the
// identical -D set reads every member at the wrong offset.  ACLIC_OPT matters
// too -- -Ofast implies fast-math, so without it borderline float comparisons
// differ and timings are not comparable.
//
// The split of responsibility: the MAKEFILE owns mkFit's -D and codegen flags;
// ACLiC owns ROOT's -I and -std, taken from the RUNNING ROOT.  Do not let the
// Makefile supply ROOT flags -- invoked without ROOTSYS its root-config falls
// back to system ROOT and would inject the wrong headers and -std.

void seedgeom_load() {
  gSystem->Load("libMkFitCore.so");

  const TString bld = gSystem->pwd();                // .../standalone
  const TString src = gSystem->DirName(bld.Data());  // .../src
  // Where THIS study lives.  __FILE__ is the path .x was given; make it absolute.
  TString here = gSystem->DirName(__FILE__);
  if (!gSystem->IsAbsoluteFileName(here))
    here = bld + "/" + here;
  const TString ext = gSystem->DirName(here.Data());  // the mkFit-external checkout holding it

  // SADIR lives in the build-dir Makefile that configure wrote; passing it on
  // the command line works for build dirs generated before `echo-aclic`
  // existed, which the plain `make echo-aclic` would not.
  TString sadir = gSystem->GetFromPipe("sed -n 's/^export SADIR *:= *//p' Makefile");
  sadir = sadir.Strip(TString::kBoth);
  if (sadir.IsNull())
    sadir = src + "/RecoTracker/MkFitCore/standalone";

  // Run make under the ROOT WE ARE RUNNING.  Without this the sub-shell's
  // root-config resolves to whatever is first on PATH -- system ROOT here --
  // and Makefile.config derives CXX_STD (and potentially more) from it, so the
  // flags would describe a different build than the library we link against.
  const TString rootsys = gROOT->GetRootSys();
  const TString cmd = Form(
      "ROOTSYS=%s PATH=%s/bin:$PATH "
      "make -s --no-print-directory -C objs-Core -f %s/Makefile SADIR=%s SRCDIR=%s BLDDIR=%s echo-aclic",
      rootsys.Data(), rootsys.Data(), sadir.Data(), sadir.Data(), src.Data(), bld.Data());
  const TString out = gSystem->GetFromPipe(cmd);

  TString defs, opt;
  std::unique_ptr<TObjArray> lines(out.Tokenize("\n"));
  for (auto *o : TRangeDynCast<TObjString>(*lines)) {
    const TString &l = o->GetString();
    if (l.BeginsWith("ACLIC_DEFS="))
      defs = l(11, l.Length());
    else if (l.BeginsWith("ACLIC_OPT="))
      opt = l(10, l.Length());
  }
  if (defs.IsNull()) {
    Error("seedgeom_load",
          "could not get flags from `make echo-aclic`; refusing to guess -- a wrong -D set\n"
          "         does not fail to link, it silently misreads every Event member.\n"
          "         Tried: %s",
          cmd.Data());
    return;
  }

  // ext first, so "mkfit-standalone-seedgeom/..." resolves to THIS checkout; the
  // build's own mkFit-external second, for everything the library was built with.
  gSystem->AddIncludePath(
      Form(" -I%s -I%s -I%s/mkFit-external %s ", src.Data(), ext.Data(), bld.Data(), defs.Data()));
  gSystem->SetFlagsOpt(opt);

  // Keep the generated .so / dictionary out of the mkFit-external checkout.
  // ACLiC compiles with cwd = build dir and hands the linker the mkFit libs by
  // the RELATIVE names the shell loaded them under ("./libMkFitCore.so", from
  // LD_LIBRARY_PATH=.), so those names have to resolve there too.
  gSystem->mkdir("test-seedgeom/aclic", kTRUE);
  for (auto l : {"libMkFitCore.so", "libMkFitCMS.so", "libMkFitRootDataFormats.so", "CMS-phase2.so"})
    gSystem->Symlink(Form("%s/%s", bld.Data(), l), Form("test-seedgeom/aclic/%s", l));
  gSystem->SetBuildDir("test-seedgeom/aclic", kTRUE);

  printf("[seedgeom] defs: %s\n", defs.Data());
  printf("[seedgeom] opt : %s\n", opt.Data());
  printf("[seedgeom] src : %s\n", here.Data());
  gROOT->ProcessLine(Form(".L %s/SeedGeom.cc+O", here.Data()));
  printf("[seedgeom] loaded via ACLiC.\n");
}
