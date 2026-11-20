#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ namespace mkfit::Const;
#pragma link C++ defined_in namespace mkfit::Const;

#pragma link C++ namespace mkfit::Config;
#pragma link C++ defined_in namespace mkfit::Config;

#pragma link C++ typedef mkfit::SVector6;
#pragma link C++ typedef mkfit::SMatrixSym66;
#pragma link C++ typedef mkfit::SMatrix66;

#pragma link C++ class mkfit::Shell;
#pragma link C++ class mkfit::Event;
#pragma link C++ class mkfit::MeasurementState + ;
#pragma link C++ class mkfit::Hit + ;
#pragma link C++ class mkfit::Hit::PackedData + ;
#pragma link C++ class mkfit::HitOnTrack + ;
#pragma link C++ class std::vector < mkfit::HitOnTrack> + ;
#pragma link C++ class mkfit::TrackState + ;
#pragma link C++ class mkfit::SimHitState + ;
#pragma link C++ class std::vector < mkfit::SimHitState> + ;
#pragma link C++ typedef mkfit::SHSVec;
#pragma link C++ class mkfit::TrackBase + ;
#pragma link C++ class mkfit::TrackBase::Status + ;
#pragma link C++ enum  mkfit::TrackBase::TrackAlgorithm + ;
#pragma link C++ class mkfit::Track + ;
#pragma link C++ class std::vector < mkfit::Track> + ;
#pragma link C++ typedef mkfit::HitVec;
#pragma link C++ typedef mkfit::HoTVec;
#pragma link C++ typedef mkfit::TrackVec;

#pragma link C++ class mkfit::ModuleShape + ;
#pragma link C++ class mkfit::ModuleInfo + ;
#pragma link C++ class mkfit::LayerInfo + ;
#pragma link C++ class mkfit::TrackerInfo + ;

#pragma link C++ class mkfit::MkBuilder + ;
#pragma link C++ class mkfit::MkJob + ;

// RntDumper

#pragma link C++ class RntDumper;

// RDF Trace etc

#ifdef MKFIT_TRACE
#pragma link C++ class mkfit::RdfSources - ;
#pragma link C++ class mkfit::RdfVectorSources - ;
#pragma link C++ class mkfit::RdfCtx - ;

#pragma link C++ class CanvasGroup - ;
// #pragma link C++ class CanvasGroup::Entry - ;
#pragma link C++ typedef CGrp;

#pragma link C++ class AnRun - ;

#pragma link C++ function np();
#pragma link C++ function nc(int, int, const char *, const char *, const char *);
#endif
