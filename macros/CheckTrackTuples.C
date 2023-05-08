// ----------------------------------------------------------------------------
// 'CheckTrackTuples.C'
// Derek Anderson
// 05.08.2023
//
// Creates a chain out of the ntp_track and ntp_gtrack
// evaluator tuples and draws some leaves to check the
// output.
// ----------------------------------------------------------------------------

// standard c includes
#include <cassert>
#include <fstream>
#include <iostream>
// root includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;

// i/o parameters
static const TString SInListDef("sPhenixG4_forPtCheck_embedOnly0300s_g4svtxeval.run6n100pt020pim.d8m5y2023.list");
static const TString SGoodListDef("checkingTrackTuples.goodFiles_embedOnly0300s.run6n100pt020pim.d8m5y2023.list");
static const TString SBadListDef("checkingTrackTuples.badFiles_embedOnly0300s.run6n100pt020pim.d8m5y2023.list");
static const TString SOutFileDef("checkingTrackTuples.embedOnly0300s.run6n100pt020pim.d8m5y2023.root");



void CheckTrackTuples(const TString sInList = SInListDef, const TString sGoodList = SGoodListDef, const TString sBadList = SBadListDef, const TString sOutFile = SOutFileDef) {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Checking track study tuples..." << endl;

  // initialize output
  TFile  *fOut    = new TFile(sOutFile.Data(), "recreate");
  TChain *tTrack  = new TChain("ntp_track");
  TChain *tGTrack = new TChain("ntp_gtrack");
  cout << "    Opened output file and declared chains." << endl;

  // open streams
  ifstream files(sInList.Data());
  ofstream good(sGoodList.Data());
  ofstream bad(sBadList.Data());
  if (!files) {
    cerr << "PANIC: input file couldn't be opened!\n" << endl;
    return;
  }
  cout << "    Opened streams.\n"
       << "    Reading in files..."
       << endl;

  // loop over files
  Int_t  trkBytes;
  Int_t  gtrkBytes;
  string sLine;
  while (files) {

    // stream file name
    files >> sLine;
    TString sFileName(sLine);

    // try adding tuples to chain
    trkBytes  = tTrack  -> Add(sFileName.Data(), 0);
    gtrkBytes = tGTrack -> Add(sFileName.Data(), 0);

    // check if good
    const Bool_t isTrkTupleGood  = (trkBytes  > 0);
    const Bool_t isGTrkTupleGood = (gtrkBytes > 0);
    const Bool_t isFileGood      = (isTrkTupleGood && isGTrkTupleGood);
    if (isFileGood) {
      cout << "      Added file '" << sLine << "'..." << endl;
      good << sLine;
      good << endl;
    } else {
      cout << "      Bad file:  '" << sLine << "'..." << endl;
      bad << sLine;
      bad << endl;
    }
  }  // end file loop
  cout << "    Finished reading in files." << endl;

  // check entries
  TCanvas *cTrkPt = new TCanvas("cTrkPt", "pt from ntp_track", 700, 500);
  cTrkPt -> cd();
  tTrack -> Draw("ntp_track.pt");
  fOut   -> cd();
  cTrkPt -> Write();
  cTrkPt -> Close();

  TCanvas *cGTrkPt = new TCanvas("cGTrkPt", "gpt from ntp_gtrack", 700, 500);
  cGTrkPt -> cd();
  tGTrack -> Draw("ntp_gtrack.gpt");
  fOut    -> cd();
  cGTrkPt -> Write();
  cGTrkPt -> Close();
  cout << "    Made plots for checking." << endl; 

  // save chains and close output
  fOut    -> cd();
  tTrack  -> Write();
  tGTrack -> Write();
  fOut    -> Close();

  cout << "  Finished checking tuples!\n" << endl;
  return;

}

// end ------------------------------------------------------------------------
