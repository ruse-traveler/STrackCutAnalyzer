// ----------------------------------------------------------------------------
// 'DoTrackCutStudy.C'
// Derek Anderson
// 12.15.2022
//
// Runs the 'STrackCutStudy'
// class.
// ----------------------------------------------------------------------------

#ifndef DOTRACKCUTSTUDY_C
#define DOTRACKCUTSTUDY_C

// standard c includes
#include <cstdlib>
#include <utility>
// root includes
#include <TROOT.h>
#include <TString.h>
// user includes
#include </sphenix/user/danderson/install/include/strackcutstudy/STrackCutStudy.h>

using namespace std;

// load libraries
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libstrackcutstudy.so)

// global constants
static const Ssiz_t NTxt = 3;



void DoTrackCutStudy() {

  // lower verbosity
  gErrorIgnoreLevel = kWarning;

  // i/o parameters
  const TString sOutFile("trackCutStudy.forMvtxCheck_withMvtxCut_finePtBinsWithNoIntNorm.pt020n5pim.d21m2y2023.root");
  const TString sInFileEO("input/embed_only/final_merge/sPhenixG4_forTrackCutStudy_embedOnly0t1099_g4svtxeval.pt020n5pim.d12m1y2023.root");
  const TString sInFilePU("input/test/sPhenixG4_testWithPileup001_g4svtxEval.d18m12y2022.root");
  const TString sInTupleEO("ntp_track");
  const TString sInTuplePU("ntp_gtrack");
  const TString sInClusterEO("ntp_cluster");

  // study parameters
  const Bool_t   doIntNorm(false);
  const Bool_t   doAvgClusterCalc(true);
  const Double_t normalPtFracMin(0.20);
  const Double_t normalPtFracMax(1.20);

  // cut flags
  const Bool_t doPrimaryCut = true;
  const Bool_t doMVtxCut    = true;
  const Bool_t doVzCut      = true;
  const Bool_t doDcaXyCut   = true;
  const Bool_t doDcaZcut    = true;
  const Bool_t doQualityCut = true;

  // track cuts
  const pair<UInt_t,   UInt_t>   nMVtxRange   = {0,    100};
  const pair<Double_t, Double_t> vzRange      = {-5.,  5.};
  const pair<Double_t, Double_t> dcaXyRange   = {-20., 20.};
  const pair<Double_t, Double_t> dcaZrange    = {-20., 20.};
  const pair<Double_t, Double_t> qualityRange = {0.,   10.};

  // text for plot
  const TString sTxtEO[NTxt] = {"#bf{#it{sPHENIX}} Simulation", "single #pi^{-}, p_{T} #in (0, 20) GeV/c", "#bf{Embedded Only Tracks}"};
  const TString sTxtPU[NTxt] = {"#bf{#it{sPHENIX}} Simulation", "0-20 fm Hijing, 50 kHz pileup", "#bf{With Pileup Tracks}"};

  // run track cut study
  STrackCutStudy *study = new STrackCutStudy();
  study -> SetInputOutputFiles(sInFileEO, sInFilePU, sOutFile);
  study -> SetInputTuples(sInTupleEO, sInTuplePU);
  study -> SetStudyParameters(doIntNorm, doAvgClusterCalc, normalPtFracMin, normalPtFracMax);
  study -> SetCutFlags(doPrimaryCut, doMVtxCut, doVzCut, doDcaXyCut, doDcaZcut, doQualityCut);
  study -> SetTrackCuts(nMVtxRange, vzRange, dcaXyRange, dcaZrange, qualityRange);
  study -> SetPlotText(NTxt, NTxt, sTxtEO, sTxtPU);
  study -> Init();
  study -> Analyze();
  study -> End();

}  // end 'DoTrackCutStudy()'

#endif

// end ------------------------------------------------------------------------
