// #pragma link C++ class ROOT::Experimental::REveVector+;

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ typedef EVec3;
#pragma link C++ class EBiVec3 + ;

#pragma link C++ class PropInfo + ;

// RntDumper

#pragma link C++ class HeaderLayer + ;
#pragma link C++ class SimSeedInfo + ;
#pragma link C++ class BinSearch + ;

#pragma link C++ class mkfit::IdxChi2List + ;
#pragma link C++ class HitInfo + ;
#pragma link C++ class HitMatchInfo + ;
#pragma link C++ class std::vector < HitInfo> + ;
#pragma link C++ class std::vector < HitMatchInfo> + ;

#pragma link C++ class CandInfo + ;
#pragma link C++ class std::vector < CandInfo> + ;

#pragma link C++ class FailedPropInfo + ;
#pragma link C++ class std::vector < FailedPropInfo> + ;

// RDF Trace

#pragma link C++ class TrCandMeta + ;
#pragma link C++ class std::vector < TrCandMeta > +;

#pragma link C++ class TrCandStage + ;
#pragma link C++ class std::vector < TrCandStage > +;

#pragma link C++ class TrCandState + ;
#pragma link C++ class std::vector < TrCandState > +;

#pragma link C++ class TrLayerSearch + ;
#pragma link C++ class std::vector < TrLayerSearch > +;

#pragma link C++ class TrHitMatch + ;
#pragma link C++ class std::vector < TrHitMatch > +;

#pragma link C++ class TrKalmanUpdate + ;
#pragma link C++ class std::vector < TrKalmanUpdate > +;

#pragma link C++ class TrBkFitUpdate + ;
#pragma link C++ class std::vector < TrBkFitUpdate > +;

// #pragma link C++ class KalmanInfo + ;

// Probe -- records written by the probes, MkFitCMS/standalone/Probe

#pragma link C++ class ProbeClosure + ;
#pragma link C++ class std::vector < ProbeClosure > +;

#pragma link C++ class ProbeStep + ;
#pragma link C++ class std::vector < ProbeStep > +;


#pragma link C++ class ProbeCfgInfo + ;
#pragma link C++ class std::vector < ProbeCfgInfo > +;

#pragma link C++ class ProbeCovXport + ;
#pragma link C++ class std::vector < ProbeCovXport > +;

#pragma link C++ class ProbeSearchHit + ;
#pragma link C++ class std::vector < ProbeSearchHit > +;

#pragma link C++ class ProbeSearchMiss + ;
#pragma link C++ class std::vector < ProbeSearchMiss > +;

#pragma link C++ class ProbeCovStep + ;
#pragma link C++ class std::vector < ProbeCovStep > +;

#pragma link C++ class SeedVecInsp + ;
