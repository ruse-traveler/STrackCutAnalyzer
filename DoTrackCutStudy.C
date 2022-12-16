// 'DoTrackCutStudy.C'
// Derek Anderson
// 12.15.2022
//
// Runs the 'STrackCutStudy'
// class.

#ifndef DOTRACKCUTSTUDY_C
#define DOTRACKCUTSTUDY_C

// root includes
#include <TROOT.h>
#include <TString.h>
// user includes
#include </sphenix/u/danderson/install/include/strackcutstudy/STrackCutStudy.h>

// load libraries
R__LOAD_LIBRARY(/sphenix/u/danderson/install/lib/libstrackcutstudy.so)



void DoTrackCutStudy() {

  // lower verbosity
  gErrorIgnoreLevel = kWarning;

  // i/o parameters
  const TString sOutFile("trackCutStudy.onlySinglePiM.d16m12y2022.root");
  const TString sInFileEO("input/merge/sPhenixG4_forTrackCutStudy_embedOnly0t199_g4svtxeval.d5m12y2022.root");
  const TString sInFilePU("input/merge/sPhenixG4_forTrackCutStudy_withPileup0t99_g4svtxeval.d5m12y2022.root");
  const TString sInTupleEO("ntp_track");
  const TString sInTuplePU("ntp_gtrack");

  // calculation parameters
  // TODO: rename these to be more accurate
  const Double_t weirdPtFracMin(0.20);
  const Double_t weirdPtFracMax(1.20);

  // run track cut study
  STrackCutStudy *study = new STrackCutStudy();
  study -> SetInputOutputFiles(sInFileEO, sInFilePU, sOutFile);
  study -> SetInputTuples(sInTupleEO, sInTuplePU);
  study -> SetWeirdFractionCuts(weirdPtFracMin, weirdPtFracMax);
  study -> Init();
  study -> Analyze();
  study -> End();

}  // end 'DoTrackCutStudy()'

#endif

// end ------------------------------------------------------------------------
