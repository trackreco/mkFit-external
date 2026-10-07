#ifndef RecoTracker_MkFitCore_standalone_RdfTrace_RdfVectorSource_h
#define RecoTracker_MkFitCore_standalone_RdfTrace_RdfVectorSource_h

#include "RecoTracker/MkFitCore/standalone/DataFormats/RntStructs.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"
#include "RecoTracker/MkFitCore/standalone/RdfTrace/RdfSources.h"

#include "ROOT/RDataFrame.hxx"

namespace mkfit {

  // ---------------------------------------------------------------------------
  // VectorBackedRDataSource -- one RDF entry per element of a std::vector<T>.
  //
  // Columns are discovered from the ROOT dictionary of T (TClass data members)
  // and read by pointer arithmetic on the member offsets, so any dictionary-ed
  // struct works with no per-type boilerplate. T must be in a LinkDef, see
  // DataFormats/RootDFs_LinkDef.h.
  //
  // STATUS: not used by the current analysis path, which is one-entry-per-Event
  // via RdfSources::MakeEventDF() / EventSource instead. Kept because it is a
  // complete, compiling, runnable example of a custom ROOT RDataSource -- worth
  // more as working code than as notes. Its smoke test is
  // Shell::TestVectorSource() -> AnRun::RunOldVecBased().
  //
  // The template itself lives in RdfVectorSource.cc: it is only ever
  // instantiated by the factories below, so it needs no header exposure.
  // ---------------------------------------------------------------------------

  class RdfVectorSources {
  public:
    static ROOT::RDataFrame MakeTrCandMetaDF(const Event &ev);
    static ROOT::RDataFrame MakeTrCandStateDF(const Event &ev);
    static ROOT::RDataFrame MakeTrHitMatchDF(const Event &ev);

    static ROOT::RDataFrame MakeTrackDF(const TrackVec &tvec);
    static ROOT::RDataFrame MakeSeedDF(const Event &ev);
  };

}  // end namespace mkfit

#endif
