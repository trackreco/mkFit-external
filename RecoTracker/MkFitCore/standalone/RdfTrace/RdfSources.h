// filepath: /foo/matevz/mic-dev/current/src/RecoTracker/MkFitCore/standalone/RntDumper/RntSources.h
#ifndef RecoTracker_MkFitCore_standalone_RdfTrace_RdfSources_h
#define RecoTracker_MkFitCore_standalone_RdfTrace_RdfSources_h

#include "RecoTracker/MkFitCore/standalone/DataFormats/RntStructs.h"
#include "RecoTracker/MkFitCore/standalone/Event.h"
#include "ROOT/RDataFrame.hxx"

#include <functional>

namespace mkfit {

  // LambdaColumnReader -- trivial RDF column reader that defers to a lambda
  // returning the address of the value for the current entry. Shared by the
  // data sources in RdfSources.cc and RdfVectorSource.cc.

  class LambdaColumnReader : public ROOT::Detail::RDF::RColumnReaderBase {
    private:
    std::function<void* ()> lambda_;

    public:
    explicit LambdaColumnReader(std::function<void* ()> lambda) : lambda_(std::move(lambda)) {}

    void *GetImpl(Long64_t gimpl_entry) override {
      return lambda_();
    }
  };

  // RdfSources, static factory for creating RDataFrames from the Event.
  //
  // One RDF entry = one mkfit::Event: the source exposes a single "event"
  // column and everything else is a Define() projection off that pointer.
  //
  // For the older one-entry-per-vector-element sources, see
  // RdfVectorSource.h / mkfit::RdfVectorSources.

  class RdfSources {
    public:
    static ROOT::RDataFrame MakeEventDF(std::vector<const Event*>& events);
  };

  // RdfCtx -- to access event and potentially other things in lambda column readers, without making them capture-heavy.
  // Use as:
  // auto C = mkfit::MakeCtx(s.event());
  // r.Define("kalmanIdx", [C](int id) { return C.ev->kalmanByState_[id]; }, {"id"});
  //
  // Also, can do this, in C++ and TRint, and use EV-> instead of C.ev-> in the lambda bodies:
  // #define EV C.ev

  struct RdfCtx {
    const Event& ev;
    const TrackerInfo& trk_info;
  };

  inline RdfCtx MakeCtx(const Event* ev, const TrackerInfo& ti) { return {*ev, ti}; }

}  // end namespace mkfit

#endif
