// ----------------------------------------------------------------------------
// 'QuickDeltaPtExtractor.C'
// Derek Anderson
// 04.25.2023
//
// Quickly apply cuts to and extract plots of the track
// DeltaPt/Pt from the track evaluator.
// ----------------------------------------------------------------------------

// standard c includes
#include <iostream>
// root includes
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"
#include "TError.h"
#include "TString.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const Ssiz_t  NTxt     = 3;
static const Ssiz_t  NPad     = 2;
static const Ssiz_t  NVtx     = 4;
static const Ssiz_t  NCuts    = 7;
static const Ssiz_t  NRange   = 2;
static const Ssiz_t  NTrkCuts = 6;
static const TString SInTrack = "ntp_track";
static const TString SInTruth = "ntp_gtrack";

// default parameters
static const TString SInDef  = "input/embed_only/final_merge/sPhenixG4_forSectorCheck_embedScanOn_embedOnly.pt020n5pim.d11m4y2023.root";
static const TString SOutDef = "varyDeltaPtCut.withInttCutAndPtDeltaVsTrack.pt020n5pim.d4m5y2023.root";



void QuickDeltaPtExtractor(const TString sInput = SInDef, const TString sOutput = SOutDef) {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning delta-pt extractor script..." << endl;

  // cut parameters
  const UInt_t   nInttTrkMin       = 1;
  const UInt_t   nMVtxTrkMin       = 2;
  const UInt_t   nTpcTrkMin        = 35;
  const Double_t qualTrkMax        = 10.;
  const Double_t vzTrkMax          = 10.;
  const Double_t ptTrkMin          = 0.1;
  const Double_t ptDeltaMax[NCuts] = {0.5, 0.25, 0.1, 0.05, 0.03, 0.02, 0.01};
  const Double_t normRange[NRange] = {0.2, 1.2};

  // histogram parameters
  const TString sPtTrueBase       = "PtTrue";
  const TString sPtRecoBase       = "PtReco";
  const TString sPtFracBase       = "PtFrac";
  const TString sPtDeltaBase      = "DeltaPt";
  const TString sPtTrkTruBase     = "PtTrkTruth";
  const TString sRejectBase       = "Reject";
  const TString sEffBase          = "Efficiency";
  const TString sDPtSuffix[NCuts] = {"_dPt50", "_dPt25", "_dPt10", "_dPt05", "_dPt03", "_dPt02", "_dPt01"};

  // histogram text parameters
  const TString sTitle        = "";
  const TString sCounts       = "counts";
  const TString sPtTrueAxis   = "p_{T}^{true} [GeV/c]";
  const TString sPtRecoAxis   = "p_{T}^{reco} [GeV/c]";
  const TString sPtFracAxis   = "p_{T}^{reco} / p_{T}^{true}";
  const TString sPtDeltaAxis  = "#Deltap_{T} / p_{T}^{reco}";
  const TString sDeltaCutAxis = "max #Deltap_{T} / p_{T}^{reco}";
  const TString sRejectAxis   = "rejection factor";
  const TString sEffAxis      = "#epsilon_{trk}";

  // histogram style parameters
  const UInt_t  fColTrue(923);
  const UInt_t  fColPure(923);
  const UInt_t  fColTrk(809);
  const UInt_t  fMarTrue(20);
  const UInt_t  fMarPure(20);
  const UInt_t  fMarTrk(46);
  const UInt_t  fColCut[NCuts]      = {899, 909, 879, 889, 859, 869, 839};
  const UInt_t  fMarCut[NCuts]      = {24,  26,  32,  25,  27,  28,  30};
  const Float_t rPtRange[NRange]    = {0., 30.};
  const Float_t rFracRange[NRange]  = {0., 4.};
  const Float_t rDeltaRange[NRange] = {0., 1.};

  // legend parameters
  const TString sLegTrue("truth");
  const TString sLegTrack("tracks (w/ cuts)");
  const TString sInfo[NTxt]        = {"#bf{#it{sPHENIX}} Simulation",
                                      "100 #pi^{-}/event, p_{T} #in (0, 20) GeV/c",
                                      "#bf{Only #pi^{-}}"};
  const TString sLegCut[NCuts]     = {"#Deltap_{T} / p_{T} < 0.5",
                                      "#Deltap_{T} / p_{T} < 0.25",
                                      "#Deltap_{T} / p_{T} < 0.1",
                                      "#Deltap_{T} / p_{T} < 0.05",
                                      "#Deltap_{T} / p_{T} < 0.03",
                                      "#Deltap_{T} / p_{T} < 0.02",
                                      "#Deltap_{T} / p_{T} < 0.01"};
  const TString sTrkCuts[NTrkCuts] = {"|v_{z}| < 10 cm",
                                      "N_{hit}^{intt} #geq 1",
                                      "N_{hit}^{mvtx} > 2",
                                      "N_{hit}^{tpc} > 35",
                                      "p_{T}^{reco} > 0.1 GeV/c",
                                      "quality < 10"};

  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput  = new TFile(sInput.Data(),  "read");
  if (!fInput || !fOutput) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fInput  = " << fInput  << "\n"
         << "       fOutput = " << fOutput << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input tuples
  TNtuple *ntTrack = (TNtuple*) fInput -> Get(SInTrack.Data());
  TNtuple *ntTruth = (TNtuple*) fInput -> Get(SInTruth.Data());
  if (!ntTrack || !ntTruth) {
    cerr << "PANIC: couldn't grab aninput tuple!\n"
         << "       ntTrack = " << ntTrack << "\n"
         << "       ntTruth = " << ntTruth << "\n"
         << endl;
    return;
  }
  cout << "    Grabbed input tuples." << endl;

  // declare track tuple addresses
  Float_t trk_event;
  Float_t trk_seed;
  Float_t trk_trackID;
  Float_t trk_crossing;
  Float_t trk_px;
  Float_t trk_py;
  Float_t trk_pz;
  Float_t trk_pt;
  Float_t trk_eta;
  Float_t trk_phi;
  Float_t trk_deltapt;
  Float_t trk_deltaeta;
  Float_t trk_deltaphi;
  Float_t trk_charge;
  Float_t trk_quality;
  Float_t trk_chisq;
  Float_t trk_ndf;
  Float_t trk_nhits;
  Float_t trk_nmaps;
  Float_t trk_nintt;
  Float_t trk_ntpc;
  Float_t trk_nmms;
  Float_t trk_ntpc1;
  Float_t trk_ntpc11;
  Float_t trk_ntpc2;
  Float_t trk_ntpc3;
  Float_t trk_nlmaps;
  Float_t trk_nlintt;
  Float_t trk_nltpc;
  Float_t trk_nlmms;
  Float_t trk_layers;
  Float_t trk_vertexID;
  Float_t trk_vx;
  Float_t trk_vy;
  Float_t trk_vz;
  Float_t trk_dca2d;
  Float_t trk_dca2dsigma;
  Float_t trk_dca3dxy;
  Float_t trk_dca3dxysigma;
  Float_t trk_dca3dz;
  Float_t trk_dca3dzsigma;
  Float_t trk_pcax;
  Float_t trk_pcay;
  Float_t trk_pcaz;
  Float_t trk_gtrackID;
  Float_t trk_gflavor;
  Float_t trk_gnhits;
  Float_t trk_gnmaps;
  Float_t trk_gnintt;
  Float_t trk_gntpc;
  Float_t trk_gnmms;
  Float_t trk_gnlmaps;
  Float_t trk_gnlintt;
  Float_t trk_gnltpc;
  Float_t trk_gnlmms;
  Float_t trk_gpx;
  Float_t trk_gpy;
  Float_t trk_gpz;
  Float_t trk_gpt;
  Float_t trk_geta;
  Float_t trk_gphi;
  Float_t trk_gvx;
  Float_t trk_gvy;
  Float_t trk_gvz;
  Float_t trk_gvt;
  Float_t trk_gfpx;
  Float_t trk_gfpy;
  Float_t trk_gfpz;
  Float_t trk_gfx;
  Float_t trk_gfy;
  Float_t trk_gfz;
  Float_t trk_gembed;
  Float_t trk_gprimary;
  Float_t trk_nfromtruth;
  Float_t trk_nwrong;
  Float_t trk_ntrumaps;
  Float_t trk_ntruintt;
  Float_t trk_ntrutpc;
  Float_t trk_ntrumms;
  Float_t trk_ntrutpc1;
  Float_t trk_ntrutpc11;
  Float_t trk_ntrutpc2;
  Float_t trk_ntrutpc3;
  Float_t trk_layersfromtruth;
  Float_t trk_nhittpcall;
  Float_t trk_nhittpcin;
  Float_t trk_nhittpcmid;
  Float_t trk_nhittpcout;
  Float_t trk_nclusall;
  Float_t trk_nclustpc;
  Float_t trk_nclusintt;
  Float_t trk_nclusmaps;
  Float_t trk_nclusmms;

  // declare truth tuple addresses
  Float_t tru_event;
  Float_t tru_seed;
  Float_t tru_gntracks;
  Float_t tru_gtrackID;
  Float_t tru_gflavor;
  Float_t tru_gnhits;
  Float_t tru_gnmaps;
  Float_t tru_gnintt;
  Float_t tru_gnmms;
  Float_t tru_gnintt1;
  Float_t tru_gnintt2;
  Float_t tru_gnintt3;
  Float_t tru_gnintt4;
  Float_t tru_gnintt5;
  Float_t tru_gnintt6;
  Float_t tru_gnintt7;
  Float_t tru_gnintt8;
  Float_t tru_gntpc;
  Float_t tru_gnlmaps;
  Float_t tru_gnlintt;
  Float_t tru_gnltpc;
  Float_t tru_gnlmms;
  Float_t tru_gpx;
  Float_t tru_gpy;
  Float_t tru_gpz;
  Float_t tru_gpt;
  Float_t tru_geta;
  Float_t tru_gphi;
  Float_t tru_gvx;
  Float_t tru_gvy;
  Float_t tru_gvz;
  Float_t tru_gvt;
  Float_t tru_gfpx;
  Float_t tru_gfpy;
  Float_t tru_gfpz;
  Float_t tru_gfx;
  Float_t tru_gfy;
  Float_t tru_gfz;
  Float_t tru_gembed;
  Float_t tru_gprimary;
  Float_t tru_trackID;
  Float_t tru_px;
  Float_t tru_py;
  Float_t tru_pz;
  Float_t tru_pt;
  Float_t tru_eta;
  Float_t tru_phi;
  Float_t tru_deltapt;
  Float_t tru_deltaeta;
  Float_t tru_deltaphi;
  Float_t tru_charge;
  Float_t tru_quality;
  Float_t tru_chisq;
  Float_t tru_ndf;
  Float_t tru_nhits;
  Float_t tru_layers;
  Float_t tru_nmaps;
  Float_t tru_nintt;
  Float_t tru_ntpc;
  Float_t tru_nmms;
  Float_t tru_ntpc1;
  Float_t tru_ntpc11;
  Float_t tru_ntpc2;
  Float_t tru_ntpc3;
  Float_t tru_nlmaps;
  Float_t tru_nlintt;
  Float_t tru_nltpc;
  Float_t tru_nlmms;
  Float_t tru_vertexID;
  Float_t tru_vx;
  Float_t tru_vy;
  Float_t tru_vz;
  Float_t tru_dca2d;
  Float_t tru_dca2dsigma;
  Float_t tru_dca3dxy;
  Float_t tru_dca3dxysigma;
  Float_t tru_dca3dz;
  Float_t tru_dca3dzsigma;
  Float_t tru_pcax;
  Float_t tru_pcay;
  Float_t tru_pcaz;
  Float_t tru_nfromtruth;
  Float_t tru_nwrong;
  Float_t tru_ntrumaps;
  Float_t tru_ntruintt;
  Float_t tru_ntrutpc;
  Float_t tru_ntrumms;
  Float_t tru_ntrutpc1;
  Float_t tru_ntrutpc11;
  Float_t tru_ntrutpc2;
  Float_t tru_ntrutpc3;
  Float_t tru_layersfromtruth;
  Float_t tru_nhittpcall;
  Float_t tru_nhittpcin;
  Float_t tru_nhittpcmid;
  Float_t tru_nhittpcout;
  Float_t tru_nclusall;
  Float_t tru_nclustpc;
  Float_t tru_nclusintt;
  Float_t tru_nclusmaps;
  Float_t tru_nclusmms;

  // set track branch addresses
  ntTrack -> SetBranchAddress("event",           &trk_event);
  ntTrack -> SetBranchAddress("seed",            &trk_seed);
  ntTrack -> SetBranchAddress("trackID",         &trk_trackID);
  ntTrack -> SetBranchAddress("crossing",        &trk_crossing);
  ntTrack -> SetBranchAddress("px",              &trk_px);
  ntTrack -> SetBranchAddress("py",              &trk_py);
  ntTrack -> SetBranchAddress("pz",              &trk_pz);
  ntTrack -> SetBranchAddress("pt",              &trk_pt);
  ntTrack -> SetBranchAddress("eta",             &trk_eta);
  ntTrack -> SetBranchAddress("phi",             &trk_phi);
  ntTrack -> SetBranchAddress("deltapt",         &trk_deltapt);
  ntTrack -> SetBranchAddress("deltaeta",        &trk_deltaeta);
  ntTrack -> SetBranchAddress("deltaphi",        &trk_deltaphi);
  ntTrack -> SetBranchAddress("charge",          &trk_charge);
  ntTrack -> SetBranchAddress("quality",         &trk_quality);
  ntTrack -> SetBranchAddress("chisq",           &trk_chisq);
  ntTrack -> SetBranchAddress("ndf",             &trk_ndf);
  ntTrack -> SetBranchAddress("nhits",           &trk_nhits);
  ntTrack -> SetBranchAddress("nmaps",           &trk_nmaps);
  ntTrack -> SetBranchAddress("nintt",           &trk_nintt);
  ntTrack -> SetBranchAddress("ntpc",            &trk_ntpc);
  ntTrack -> SetBranchAddress("nmms",            &trk_nmms);
  ntTrack -> SetBranchAddress("ntpc1",           &trk_ntpc1);
  ntTrack -> SetBranchAddress("ntpc11",          &trk_ntpc11);
  ntTrack -> SetBranchAddress("ntpc2",           &trk_ntpc2);
  ntTrack -> SetBranchAddress("ntpc3",           &trk_ntpc3);
  ntTrack -> SetBranchAddress("nlmaps",          &trk_nlmaps);
  ntTrack -> SetBranchAddress("nlintt",          &trk_nlintt);
  ntTrack -> SetBranchAddress("nltpc",           &trk_nltpc);
  ntTrack -> SetBranchAddress("nlmms",           &trk_nlmms);
  ntTrack -> SetBranchAddress("layers",          &trk_layers);
  ntTrack -> SetBranchAddress("vertexID",        &trk_vertexID);
  ntTrack -> SetBranchAddress("vx",              &trk_vx);
  ntTrack -> SetBranchAddress("vy",              &trk_vy);
  ntTrack -> SetBranchAddress("vz",              &trk_vz);
  ntTrack -> SetBranchAddress("dca2d",           &trk_dca2d);
  ntTrack -> SetBranchAddress("dca2dsigma",      &trk_dca2dsigma);
  ntTrack -> SetBranchAddress("dca3dxy",         &trk_dca3dxy);
  ntTrack -> SetBranchAddress("dca3dxysigma",    &trk_dca3dxysigma);
  ntTrack -> SetBranchAddress("dca3dz",          &trk_dca3dz);
  ntTrack -> SetBranchAddress("dca3dzsigma",     &trk_dca3dzsigma);
  ntTrack -> SetBranchAddress("pcax",            &trk_pcax);
  ntTrack -> SetBranchAddress("pcay",            &trk_pcay);
  ntTrack -> SetBranchAddress("pcaz",            &trk_pcaz);
  ntTrack -> SetBranchAddress("gtrackID",        &trk_gtrackID);
  ntTrack -> SetBranchAddress("gflavor",         &trk_gflavor);
  ntTrack -> SetBranchAddress("gnhits",          &trk_gnhits);
  ntTrack -> SetBranchAddress("gnmaps",          &trk_gnmaps);
  ntTrack -> SetBranchAddress("gnintt",          &trk_gnintt);
  ntTrack -> SetBranchAddress("gntpc",           &trk_gntpc);
  ntTrack -> SetBranchAddress("gnmms",           &trk_gnmms);
  ntTrack -> SetBranchAddress("gnlmaps",         &trk_gnlmaps);
  ntTrack -> SetBranchAddress("gnlintt",         &trk_gnlintt);
  ntTrack -> SetBranchAddress("gnltpc",          &trk_gnltpc);
  ntTrack -> SetBranchAddress("gnlmms",          &trk_gnlmms);
  ntTrack -> SetBranchAddress("gpx",             &trk_gpx);
  ntTrack -> SetBranchAddress("gpy",             &trk_gpy);
  ntTrack -> SetBranchAddress("gpz",             &trk_gpz);
  ntTrack -> SetBranchAddress("gpt",             &trk_gpt);
  ntTrack -> SetBranchAddress("geta",            &trk_geta);
  ntTrack -> SetBranchAddress("gphi",            &trk_gphi);
  ntTrack -> SetBranchAddress("gvx",             &trk_gvx);
  ntTrack -> SetBranchAddress("gvy",             &trk_gvy);
  ntTrack -> SetBranchAddress("gvz",             &trk_gvz);
  ntTrack -> SetBranchAddress("gvt",             &trk_gvt);
  ntTrack -> SetBranchAddress("gfpx",            &trk_gfpx);
  ntTrack -> SetBranchAddress("gfpy",            &trk_gfpy);
  ntTrack -> SetBranchAddress("gfpz",            &trk_gfpz);
  ntTrack -> SetBranchAddress("gfx",             &trk_gfx);
  ntTrack -> SetBranchAddress("gfy",             &trk_gfy);
  ntTrack -> SetBranchAddress("gfz",             &trk_gfz);
  ntTrack -> SetBranchAddress("gembed",          &trk_gembed);
  ntTrack -> SetBranchAddress("gprimary",        &trk_gprimary);
  ntTrack -> SetBranchAddress("nfromtruth",      &trk_nfromtruth);
  ntTrack -> SetBranchAddress("nwrong",          &trk_nwrong);
  ntTrack -> SetBranchAddress("ntrumaps",        &trk_ntrumaps);
  ntTrack -> SetBranchAddress("ntruintt",        &trk_ntruintt);
  ntTrack -> SetBranchAddress("ntrutpc",         &trk_ntrutpc);
  ntTrack -> SetBranchAddress("ntrumms",         &trk_ntrumms);
  ntTrack -> SetBranchAddress("ntrutpc1",        &trk_ntrutpc1);
  ntTrack -> SetBranchAddress("ntrutpc11",       &trk_ntrutpc11);
  ntTrack -> SetBranchAddress("ntrutpc2",        &trk_ntrutpc2);
  ntTrack -> SetBranchAddress("ntrutpc3",        &trk_ntrutpc3);
  ntTrack -> SetBranchAddress("layersfromtruth", &trk_layersfromtruth);
  ntTrack -> SetBranchAddress("nhittpcall",      &trk_nhittpcall);
  ntTrack -> SetBranchAddress("nhittpcin",       &trk_nhittpcin);
  ntTrack -> SetBranchAddress("nhittpcmid",      &trk_nhittpcmid);
  ntTrack -> SetBranchAddress("nhittpcout",      &trk_nhittpcout);
  ntTrack -> SetBranchAddress("nclusall",        &trk_nclusall);
  ntTrack -> SetBranchAddress("nclustpc",        &trk_nclustpc);
  ntTrack -> SetBranchAddress("nclusintt",       &trk_nclusintt);
  ntTrack -> SetBranchAddress("nclusmaps",       &trk_nclusmaps);
  ntTrack -> SetBranchAddress("nclusmms",        &trk_nclusmms);

  // Set branch addresses.
  ntTruth -> SetBranchAddress("event",           &tru_event);
  ntTruth -> SetBranchAddress("seed",            &tru_seed);
  ntTruth -> SetBranchAddress("gntracks",        &tru_gntracks);
  ntTruth -> SetBranchAddress("gtrackID",        &tru_gtrackID);
  ntTruth -> SetBranchAddress("gflavor",         &tru_gflavor);
  ntTruth -> SetBranchAddress("gnhits",          &tru_gnhits);
  ntTruth -> SetBranchAddress("gnmaps",          &tru_gnmaps);
  ntTruth -> SetBranchAddress("gnintt",          &tru_gnintt);
  ntTruth -> SetBranchAddress("gnmms",           &tru_gnmms);
  ntTruth -> SetBranchAddress("gnintt1",         &tru_gnintt1);
  ntTruth -> SetBranchAddress("gnintt2",         &tru_gnintt2);
  ntTruth -> SetBranchAddress("gnintt3",         &tru_gnintt3);
  ntTruth -> SetBranchAddress("gnintt4",         &tru_gnintt4);
  ntTruth -> SetBranchAddress("gnintt5",         &tru_gnintt5);
  ntTruth -> SetBranchAddress("gnintt6",         &tru_gnintt6);
  ntTruth -> SetBranchAddress("gnintt7",         &tru_gnintt7);
  ntTruth -> SetBranchAddress("gnintt8",         &tru_gnintt8);
  ntTruth -> SetBranchAddress("gntpc",           &tru_gntpc);
  ntTruth -> SetBranchAddress("gnlmaps",         &tru_gnlmaps);
  ntTruth -> SetBranchAddress("gnlintt",         &tru_gnlintt);
  ntTruth -> SetBranchAddress("gnltpc",          &tru_gnltpc);
  ntTruth -> SetBranchAddress("gnlmms",          &tru_gnlmms);
  ntTruth -> SetBranchAddress("gpx",             &tru_gpx);
  ntTruth -> SetBranchAddress("gpy",             &tru_gpy);
  ntTruth -> SetBranchAddress("gpz",             &tru_gpz);
  ntTruth -> SetBranchAddress("gpt",             &tru_gpt);
  ntTruth -> SetBranchAddress("geta",            &tru_geta);
  ntTruth -> SetBranchAddress("gphi",            &tru_gphi);
  ntTruth -> SetBranchAddress("gvx",             &tru_gvx);
  ntTruth -> SetBranchAddress("gvy",             &tru_gvy);
  ntTruth -> SetBranchAddress("gvz",             &tru_gvz);
  ntTruth -> SetBranchAddress("gvt",             &tru_gvt);
  ntTruth -> SetBranchAddress("gfpx",            &tru_gfpx);
  ntTruth -> SetBranchAddress("gfpy",            &tru_gfpy);
  ntTruth -> SetBranchAddress("gfpz",            &tru_gfpz);
  ntTruth -> SetBranchAddress("gfx",             &tru_gfx);
  ntTruth -> SetBranchAddress("gfy",             &tru_gfy);
  ntTruth -> SetBranchAddress("gfz",             &tru_gfz);
  ntTruth -> SetBranchAddress("gembed",          &tru_gembed);
  ntTruth -> SetBranchAddress("gprimary",        &tru_gprimary);
  ntTruth -> SetBranchAddress("trackID",         &tru_trackID);
  ntTruth -> SetBranchAddress("px",              &tru_px);
  ntTruth -> SetBranchAddress("py",              &tru_py);
  ntTruth -> SetBranchAddress("pz",              &tru_pz);
  ntTruth -> SetBranchAddress("pt",              &tru_pt);
  ntTruth -> SetBranchAddress("eta",             &tru_eta);
  ntTruth -> SetBranchAddress("phi",             &tru_phi);
  ntTruth -> SetBranchAddress("deltapt",         &tru_deltapt);
  ntTruth -> SetBranchAddress("deltaeta",        &tru_deltaeta);
  ntTruth -> SetBranchAddress("deltaphi",        &tru_deltaphi);
  ntTruth -> SetBranchAddress("charge",          &tru_charge);
  ntTruth -> SetBranchAddress("quality",         &tru_quality);
  ntTruth -> SetBranchAddress("chisq",           &tru_chisq);
  ntTruth -> SetBranchAddress("ndf",             &tru_ndf);
  ntTruth -> SetBranchAddress("nhits",           &tru_nhits);
  ntTruth -> SetBranchAddress("layers",          &tru_layers);
  ntTruth -> SetBranchAddress("nmaps",           &tru_nmaps);
  ntTruth -> SetBranchAddress("nintt",           &tru_nintt);
  ntTruth -> SetBranchAddress("ntpc",            &tru_ntpc);
  ntTruth -> SetBranchAddress("nmms",            &tru_nmms);
  ntTruth -> SetBranchAddress("ntpc1",           &tru_ntpc1);
  ntTruth -> SetBranchAddress("ntpc11",          &tru_ntpc11);
  ntTruth -> SetBranchAddress("ntpc2",           &tru_ntpc2);
  ntTruth -> SetBranchAddress("ntpc3",           &tru_ntpc3);
  ntTruth -> SetBranchAddress("nlmaps",          &tru_nlmaps);
  ntTruth -> SetBranchAddress("nlintt",          &tru_nlintt);
  ntTruth -> SetBranchAddress("nltpc",           &tru_nltpc);
  ntTruth -> SetBranchAddress("nlmms",           &tru_nlmms);
  ntTruth -> SetBranchAddress("vertexID",        &tru_vertexID);
  ntTruth -> SetBranchAddress("vx",              &tru_vx);
  ntTruth -> SetBranchAddress("vy",              &tru_vy);
  ntTruth -> SetBranchAddress("vz",              &tru_vz);
  ntTruth -> SetBranchAddress("dca2d",           &tru_dca2d);
  ntTruth -> SetBranchAddress("dca2dsigma",      &tru_dca2dsigma);
  ntTruth -> SetBranchAddress("dca3dxy",         &tru_dca3dxy);
  ntTruth -> SetBranchAddress("dca3dxysigma",    &tru_dca3dxysigma);
  ntTruth -> SetBranchAddress("dca3dz",          &tru_dca3dz);
  ntTruth -> SetBranchAddress("dca3dzsigma",     &tru_dca3dzsigma);
  ntTruth -> SetBranchAddress("pcax",            &tru_pcax);
  ntTruth -> SetBranchAddress("pcay",            &tru_pcay);
  ntTruth -> SetBranchAddress("pcaz",            &tru_pcaz);
  ntTruth -> SetBranchAddress("nfromtruth",      &tru_nfromtruth);
  ntTruth -> SetBranchAddress("nwrong",          &tru_nwrong);
  ntTruth -> SetBranchAddress("ntrumaps",        &tru_ntrumaps);
  ntTruth -> SetBranchAddress("ntruintt",        &tru_ntruintt);
  ntTruth -> SetBranchAddress("ntrutpc",         &tru_ntrutpc);
  ntTruth -> SetBranchAddress("ntrumms",         &tru_ntrumms);
  ntTruth -> SetBranchAddress("ntrutpc1",        &tru_ntrutpc1);
  ntTruth -> SetBranchAddress("ntrutpc11",       &tru_ntrutpc11);
  ntTruth -> SetBranchAddress("ntrutpc2",        &tru_ntrutpc2);
  ntTruth -> SetBranchAddress("ntrutpc3",        &tru_ntrutpc3);
  ntTruth -> SetBranchAddress("layersfromtruth", &tru_layersfromtruth);
  ntTruth -> SetBranchAddress("nhittpcall",      &tru_nhittpcall);
  ntTruth -> SetBranchAddress("nhittpcin",       &tru_nhittpcin);
  ntTruth -> SetBranchAddress("nhittpcmid",      &tru_nhittpcmid);
  ntTruth -> SetBranchAddress("nhittpcout",      &tru_nhittpcout);
  ntTruth -> SetBranchAddress("nclusall",        &tru_nclusall);
  ntTruth -> SetBranchAddress("nclustpc",        &tru_nclustpc);
  ntTruth -> SetBranchAddress("nclusintt",       &tru_nclusintt);
  ntTruth -> SetBranchAddress("nclusmaps",       &tru_nclusmaps);
  ntTruth -> SetBranchAddress("nclusmms",        &tru_nclusmms);
  cout << "    Set track tuple branches." << endl;

  // declare histograms
  TH1D *hEff;
  TH1D *hPtTruth;
  TH1D *hPtDelta;
  TH1D *hPtTrack;
  TH1D *hPtFrac;
  TH1D *hPtTrkTru;
  TH1D *hPtDeltaCut[NCuts];
  TH1D *hPtTrackCut[NCuts];
  TH1D *hPtFracCut[NCuts];
  TH1D *hPtTrkTruCut[NCuts];
  TH1D *hEffCut[NCuts];

  TH2D *hPtDeltaVsFrac;
  TH2D *hPtDeltaVsTrue;
  TH2D *hPtDeltaVsTrack;
  TH2D *hPtTrueVsTrack;
  TH2D *hPtDeltaVsFracCut[NCuts];
  TH2D *hPtDeltaVsTrueCut[NCuts];
  TH2D *hPtDeltaVsTrackCut[NCuts];
  TH2D *hPtTrueVsTrackCut[NCuts];

  // histogram binning
  const UInt_t  nPtBins(500);
  const UInt_t  nFracBins(1000);
  const UInt_t  nDeltaBins(5000);
  const Float_t rPtBins[NRange]    = {0., 50.};
  const Float_t rFracBins[NRange]  = {0., 10.};
  const Float_t rDeltaBins[NRange] = {0., 5.};

  // create names
  TString sPtTruth("h");
  TString sPtDelta("h");
  TString sPtTrack("h");
  TString sPtFrac("h");
  TString sPtTrkTru("h");
  sPtTruth.Append(sPtTrueBase.Data());
  sPtDelta.Append(sPtDeltaBase.Data());
  sPtTrack.Append(sPtRecoBase.Data());
  sPtFrac.Append(sPtFracBase.Data());
  sPtTrkTru.Append(sPtTrkTruBase.Data());

  TString sPtDeltaVsFrac("h");
  TString sPtDeltaVsTrue("h");
  TString sPtDeltaVsTrack("h");
  TString sPtTrueVsTrack("h");
  sPtDeltaVsFrac.Append(sPtDeltaBase.Data());
  sPtDeltaVsTrue.Append(sPtDeltaBase.Data());
  sPtDeltaVsTrack.Append(sPtDeltaBase.Data());
  sPtTrueVsTrack.Append(sPtTrueBase.Data());
  sPtDeltaVsFrac.Append("Vs");
  sPtDeltaVsTrue.Append("Vs");
  sPtDeltaVsTrack.Append("Vs");
  sPtTrueVsTrack.Append("Vs");
  sPtDeltaVsFrac.Append(sPtFracBase.Data());
  sPtDeltaVsTrue.Append(sPtTrueBase.Data());
  sPtDeltaVsTrack.Append(sPtRecoBase.Data());
  sPtTrueVsTrack.Append(sPtRecoBase.Data());

  TString sPtDeltaCut[NCuts];
  TString sPtTrackCut[NCuts];
  TString sPtFracCut[NCuts];
  TString sPtTrkTruCut[NCuts];
  TString sPtDeltaVsFracCut[NCuts];
  TString sPtDeltaVsTrueCut[NCuts];
  TString sPtDeltaVsTrackCut[NCuts];
  TString sPtTrueVsTrackCut[NCuts];
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    sPtDeltaCut[iCut]  = "h";
    sPtTrackCut[iCut]  = "h";
    sPtFracCut[iCut]   = "h";
    sPtTrkTruCut[iCut] = "h";
    sPtDeltaCut[iCut].Append(sPtDeltaBase.Data());
    sPtTrackCut[iCut].Append(sPtRecoBase.Data());
    sPtFracCut[iCut].Append(sPtFracBase.Data());
    sPtTrkTruCut[iCut].Append(sPtTrkTruBase.Data());
    sPtDeltaCut[iCut].Append(sDPtSuffix[iCut].Data());
    sPtTrackCut[iCut].Append(sDPtSuffix[iCut].Data());
    sPtFracCut[iCut].Append(sDPtSuffix[iCut].Data());
    sPtTrkTruCut[iCut].Append(sDPtSuffix[iCut].Data());

    sPtDeltaVsFracCut[iCut]  = "h";
    sPtDeltaVsTrueCut[iCut]  = "h";
    sPtDeltaVsTrackCut[iCut] = "h";
    sPtTrueVsTrackCut[iCut]  = "h";
    sPtDeltaVsFracCut[iCut].Append(sPtDeltaBase.Data());
    sPtDeltaVsFracCut[iCut].Append(sPtDeltaBase.Data());
    sPtDeltaVsTrueCut[iCut].Append(sPtDeltaBase.Data());
    sPtDeltaVsTrackCut[iCut].Append(sPtDeltaBase.Data());
    sPtTrueVsTrackCut[iCut].Append(sPtTrueBase.Data());
    sPtDeltaVsFracCut[iCut].Append("Vs");
    sPtDeltaVsTrueCut[iCut].Append("Vs");
    sPtDeltaVsTrackCut[iCut].Append("Vs");
    sPtTrueVsTrackCut[iCut].Append("Vs");
    sPtDeltaVsFracCut[iCut].Append(sPtFracBase.Data());
    sPtDeltaVsTrueCut[iCut].Append(sPtTrueBase.Data());
    sPtDeltaVsTrackCut[iCut].Append(sPtRecoBase.Data());
    sPtTrueVsTrackCut[iCut].Append(sPtRecoBase.Data());
    sPtDeltaVsFracCut[iCut].Append(sDPtSuffix[iCut].Data());
    sPtDeltaVsTrueCut[iCut].Append(sDPtSuffix[iCut].Data());
    sPtDeltaVsTrackCut[iCut].Append(sDPtSuffix[iCut].Data());
    sPtTrueVsTrackCut[iCut].Append(sDPtSuffix[iCut].Data());
  }

  // initialize histograms
  hPtTruth  = new TH1D(sPtTruth.Data(),  "", nPtBins,    rPtBins[0],    rPtBins[1]);
  hPtDelta  = new TH1D(sPtDelta.Data(),  "", nDeltaBins, rDeltaBins[0], rDeltaBins[1]);
  hPtTrack  = new TH1D(sPtTrack.Data(),  "", nPtBins,    rPtBins[0],    rPtBins[1]);
  hPtFrac   = new TH1D(sPtFrac.Data(),   "", nFracBins,  rFracBins[0],  rFracBins[1]);
  hPtTrkTru = new TH1D(sPtTrkTru.Data(), "", nPtBins,    rPtBins[0],    rPtBins[1]);
  hPtTruth  -> Sumw2();
  hPtDelta  -> Sumw2();
  hPtTrack  -> Sumw2();
  hPtFrac   -> Sumw2();
  hPtTrkTru -> Sumw2();

  hPtDeltaVsFrac  = new TH2D(sPtDeltaVsFrac.Data(),  "", nFracBins, rFracBins[0], rFracBins[1], nDeltaBins, rDeltaBins[0], rDeltaBins[1]);
  hPtDeltaVsTrue  = new TH2D(sPtDeltaVsTrue.Data(),  "", nPtBins,   rPtBins[0],   rPtBins[1],   nDeltaBins, rDeltaBins[0], rDeltaBins[1]);
  hPtDeltaVsTrack = new TH2D(sPtDeltaVsTrack.Data(), "", nPtBins,   rPtBins[0],   rPtBins[1],   nDeltaBins, rDeltaBins[0], rDeltaBins[1]);
  hPtTrueVsTrack  = new TH2D(sPtTrueVsTrack.Data(),  "", nPtBins,   rPtBins[0],   rPtBins[1],   nPtBins,    rPtBins[0],    rPtBins[1]);
  hPtDeltaVsFrac  -> Sumw2();
  hPtDeltaVsTrue  -> Sumw2();
  hPtDeltaVsTrack -> Sumw2();
  hPtTrueVsTrack  -> Sumw2();

  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    hPtDeltaCut[iCut]  = new TH1D(sPtDeltaCut[iCut].Data(),  "", nDeltaBins, rDeltaBins[0], rDeltaBins[1]);
    hPtTrackCut[iCut]  = new TH1D(sPtTrackCut[iCut].Data(),  "", nPtBins,    rPtBins[0],    rPtBins[1]);
    hPtFracCut[iCut]   = new TH1D(sPtFracCut[iCut].Data(),   "", nFracBins,  rFracBins[0],  rFracBins[1]);
    hPtTrkTruCut[iCut] = new TH1D(sPtTrkTruCut[iCut].Data(), "", nPtBins,    rPtBins[0],    rPtBins[1]);
    hPtDeltaCut[iCut]  -> Sumw2();
    hPtTrackCut[iCut]  -> Sumw2();
    hPtFracCut[iCut]   -> Sumw2();
    hPtTrkTruCut[iCut] -> Sumw2();

    hPtDeltaVsFracCut[iCut]  = new TH2D(sPtDeltaVsFracCut[iCut].Data(),  "", nFracBins, rFracBins[0], rFracBins[1], nDeltaBins, rDeltaBins[0], rDeltaBins[1]);
    hPtDeltaVsTrueCut[iCut]  = new TH2D(sPtDeltaVsTrueCut[iCut].Data(),  "", nPtBins,   rPtBins[0],   rPtBins[1],   nDeltaBins, rDeltaBins[0], rDeltaBins[1]);
    hPtDeltaVsTrackCut[iCut] = new TH2D(sPtDeltaVsTrackCut[iCut].Data(), "", nPtBins,   rPtBins[0],   rPtBins[1],   nDeltaBins, rDeltaBins[0], rDeltaBins[1]);
    hPtTrueVsTrackCut[iCut]  = new TH2D(sPtTrueVsTrackCut[iCut].Data(),  "", nPtBins,   rPtBins[0],   rPtBins[1],   nPtBins,    rPtBins[0],    rPtBins[1]);
    hPtDeltaVsFracCut[iCut]  -> Sumw2();
    hPtDeltaVsTrueCut[iCut]  -> Sumw2();
    hPtDeltaVsTrackCut[iCut] -> Sumw2();
    hPtTrueVsTrackCut[iCut]  -> Sumw2();
  }

  // grab no. of entries
  const Long64_t nTrks = ntTrack -> GetEntries();
  const Long64_t nTrus = ntTruth -> GetEntries();
  cout << "    Beginning tuple loops: " << nTrks << " reco. tracks and " << nTrus << " particles to process..." << endl;

  // for reject calculation
  UInt_t   nNorm[NCuts];
  UInt_t   nWeird[NCuts];
  Double_t reject[NCuts];
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    nNorm[iCut]  = 0;
    nWeird[iCut] = 0;
    reject[iCut] = 0.;
  }

  // track loop
  Long64_t nBytesTrk = 0;
  for (Long64_t iTrk = 0; iTrk < nTrks; iTrk++) {

    // grab entry
    const Long64_t bytesTrk = ntTrack -> GetEntry(iTrk);
    if (bytesTrk < 0.) {
      cerr << "WARNING: something wrong with track #" << iTrk << "! Aborting loop!" << endl;
      break;
    }
    nBytesTrk += bytesTrk;

    // announce progress
    const Long64_t iProgTrk = iTrk + 1;
    if (iProgTrk == nTrks) {
      cout << "      Processing track " << iProgTrk << "/" << nTrks << "..." << endl;
    } else {
      cout << "      Processing track " << iProgTrk << "/" << nTrks << "...\r" << flush;
    }

    // do calculations
    const Double_t ptFrac  = trk_pt / trk_gpt;
    const Double_t ptDelta = trk_deltapt / trk_pt;

    // apply trk cuts
    const Bool_t isInZVtxCut = (abs(trk_vz) <  vzTrkMax);
    const Bool_t isInInttCut = (trk_nintt   >= nInttTrkMin);
    const Bool_t isInMVtxCut = (trk_nlmaps  >  nMVtxTrkMin);
    const Bool_t isInTpcCut  = (trk_ntpc    >  nTpcTrkMin);
    const Bool_t isInPtCut   = (trk_pt      >  ptTrkMin);
    const Bool_t isInQualCut = (trk_quality <  qualTrkMax);
    const Bool_t isGoodTrk   = (isInZVtxCut && isInInttCut && isInMVtxCut && isInTpcCut && isInPtCut && isInQualCut);
    if (!isGoodTrk) continue;

    // fill histograms
    hPtDelta        -> Fill(ptDelta);
    hPtTrack        -> Fill(trk_pt);
    hPtFrac         -> Fill(ptFrac);
    hPtTrkTru       -> Fill(trk_gpt);
    hPtDeltaVsFrac  -> Fill(ptFrac,  ptDelta);
    hPtDeltaVsTrue  -> Fill(trk_gpt, ptDelta);
    hPtDeltaVsTrack -> Fill(trk_pt,  ptDelta);
    hPtTrueVsTrack  -> Fill(trk_pt,  trk_gpt);

    // apply delta-pt cuts
    const Bool_t isNormalTrk = ((ptFrac > normRange[0]) && (ptFrac < normRange[1]));
    for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
      const Bool_t isInDeltaPtCut = (ptDelta < ptDeltaMax[iCut]);
      if (isInDeltaPtCut) {

        // fill histograms
        hPtDeltaCut[iCut]        -> Fill(ptDelta);
        hPtTrackCut[iCut]        -> Fill(trk_pt);
        hPtFracCut[iCut]         -> Fill(ptFrac);
        hPtTrkTruCut[iCut]       -> Fill(trk_gpt);
        hPtDeltaVsFracCut[iCut]  -> Fill(ptFrac,  ptDelta);
        hPtDeltaVsTrueCut[iCut]  -> Fill(trk_gpt, ptDelta);
        hPtDeltaVsTrackCut[iCut] -> Fill(trk_pt,  ptDelta);
        hPtTrueVsTrackCut[iCut]  -> Fill(trk_pt,  trk_gpt);

        // increment counters
        if (isNormalTrk) {
          ++nNorm[iCut];
        } else {
          ++nWeird[iCut];
        }
      }
    }  // end delta-pt cut
  }  // end track loop

  // calculate purities
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    reject[iCut] = (Double_t) nNorm[iCut] / (Double_t) nWeird[iCut];
  }

  // truth loop
  Long64_t nBytesTru = 0;
  for (Long64_t iTru = 0; iTru < nTrus; iTru++) {

    // grab entry
    const Long64_t bytesTru = ntTruth -> GetEntry(iTru);
    if (bytesTru < 0.) {
      cerr << "WARNING: something wrong with particle #" << iTru << "! Aborting loop!" << endl;
      break;
    }
    nBytesTru += bytesTru;

    // announce progress
    const Long64_t iProgTru = iTru + 1;
    if (iProgTru == nTrus) {
      cout << "      Processing particle " << iProgTru << "/" << nTrus << "..." << endl;
    } else {
      cout << "      Processing particle" << iProgTru << "/" << nTrus << "...\r" << flush;
    }

    // fill truth histogram
    const Bool_t isPrimary = (tru_gprimary == 1);
    if (isPrimary) {
      hPtTruth -> Fill(tru_gpt);
    }
  }  // end track loop

  // announce purities
  cout << "    Finished tuple loops! Calculated rejection factors:" << endl;
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    cout << "      n(Norm, Weird) = (" << nNorm[iCut] << ", " << nWeird[iCut] << "), rejection = " << reject[iCut] << endl;
  }

  // make rejection graph
  TString sReject("gr");
  sReject.Append(sRejectBase.Data());

  TGraph *grReject = new TGraph(NCuts, ptDeltaMax, reject);
  grReject -> SetName(sReject.Data());
  cout << "    Made rejection factor graph." << endl; 

  // calculate efficiencies
  TString sEff("h");
  sEff.Append(sEffBase.Data());

  TString sEffCut[NCuts];
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    sEffCut[iCut] = "h";
    sEffCut[iCut].Append(sEffBase.Data());
    sEffCut[iCut].Append(sDPtSuffix[iCut].Data());
  }

  hEff = (TH1D*) hPtTruth -> Clone();
  hEff -> SetName(sEff.Data());
  hEff -> Reset("ICES");
  hEff -> Divide(hPtTrkTru, hPtTruth, 1., 1.);
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    hEffCut[iCut] = (TH1D*) hPtTruth -> Clone();
    hEffCut[iCut] -> SetName(sEffCut[iCut].Data());
    hEffCut[iCut] -> Reset("ICES");
    hEffCut[iCut] -> Divide(hPtTrkTruCut[iCut], hPtTruth, 1., 1.);
  }
  cout << "    Calculated efficiencies." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1,   1.};
  const Float_t fOffY[NPad] = {0.7,   1.3};
  const Float_t fOffZ[NPad] = {1.1,   1.1};
  grReject        -> SetMarkerColor(fColTrue);
  grReject        -> SetMarkerStyle(fMarTrue);
  grReject        -> SetFillColor(fColTrue);
  grReject        -> SetFillStyle(fFil);
  grReject        -> SetLineColor(fColTrue);
  grReject        -> SetLineStyle(fLin);
  grReject        -> SetLineWidth(fWid);
  grReject        -> SetTitle(sTitle.Data());
  grReject        -> GetXaxis() -> SetRangeUser(rDeltaRange[0], rDeltaRange[1]);
  grReject        -> GetXaxis() -> SetTitle(sDeltaCutAxis.Data());
  grReject        -> GetXaxis() -> SetTitleFont(fTxt);
  grReject        -> GetXaxis() -> SetTitleSize(fTit[1]);
  grReject        -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  grReject        -> GetXaxis() -> SetLabelFont(fTxt);
  grReject        -> GetXaxis() -> SetLabelSize(fLab[1]);
  grReject        -> GetXaxis() -> CenterTitle(fCnt);
  grReject        -> GetYaxis() -> SetTitle(sRejectAxis.Data());
  grReject        -> GetYaxis() -> SetTitleFont(fTxt);
  grReject        -> GetYaxis() -> SetTitleSize(fTit[1]);
  grReject        -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  grReject        -> GetYaxis() -> SetLabelFont(fTxt);
  grReject        -> GetYaxis() -> SetLabelSize(fLab[1]);
  grReject        -> GetYaxis() -> CenterTitle(fCnt);
  hEff            -> SetMarkerColor(fColTrk);
  hEff            -> SetMarkerStyle(fMarTrk);
  hEff            -> SetFillColor(fColTrk);
  hEff            -> SetFillStyle(fFil);
  hEff            -> SetLineColor(fColTrk);
  hEff            -> SetLineStyle(fLin);
  hEff            -> SetLineWidth(fWid);
  hEff            -> SetTitle(sTitle.Data());
  hEff            -> SetTitleFont(fTxt);
  hEff            -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
  hEff            -> GetXaxis() -> SetTitle(sPtTrueAxis.Data());
  hEff            -> GetXaxis() -> SetTitleFont(fTxt);
  hEff            -> GetXaxis() -> SetTitleSize(fTit[0]);
  hEff            -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hEff            -> GetXaxis() -> SetLabelFont(fTxt);
  hEff            -> GetXaxis() -> SetLabelSize(fLab[0]);
  hEff            -> GetXaxis() -> CenterTitle(fCnt);
  hEff            -> GetYaxis() -> SetTitle(sEffAxis.Data());
  hEff            -> GetYaxis() -> SetTitleFont(fTxt);
  hEff            -> GetYaxis() -> SetTitleSize(fTit[0]);
  hEff            -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hEff            -> GetYaxis() -> SetLabelFont(fTxt);
  hEff            -> GetYaxis() -> SetLabelSize(fLab[0]);
  hEff            -> GetYaxis() -> CenterTitle(fCnt);
  hPtTruth        -> SetMarkerColor(fColTrue);
  hPtTruth        -> SetMarkerStyle(fMarTrue);
  hPtTruth        -> SetFillColor(fColTrue);
  hPtTruth        -> SetFillStyle(fFil);
  hPtTruth        -> SetLineColor(fColTrue);
  hPtTruth        -> SetLineStyle(fLin);
  hPtTruth        -> SetLineWidth(fWid);
  hPtTruth        -> SetTitle(sTitle.Data());
  hPtTruth        -> SetTitleFont(fTxt);
  hPtTruth        -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
  hPtTruth        -> GetXaxis() -> SetTitle(sPtTrueAxis.Data());
  hPtTruth        -> GetXaxis() -> SetTitleFont(fTxt);
  hPtTruth        -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPtTruth        -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPtTruth        -> GetXaxis() -> SetLabelFont(fTxt);
  hPtTruth        -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPtTruth        -> GetXaxis() -> CenterTitle(fCnt);
  hPtTruth        -> GetYaxis() -> SetTitle(sCounts.Data());
  hPtTruth        -> GetYaxis() -> SetTitleFont(fTxt);
  hPtTruth        -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPtTruth        -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPtTruth        -> GetYaxis() -> SetLabelFont(fTxt);
  hPtTruth        -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPtTruth        -> GetYaxis() -> CenterTitle(fCnt);
  hPtDelta        -> SetMarkerColor(fColTrk);
  hPtDelta        -> SetMarkerStyle(fMarTrk);
  hPtDelta        -> SetFillColor(fColTrk);
  hPtDelta        -> SetFillStyle(fFil);
  hPtDelta        -> SetLineColor(fColTrk);
  hPtDelta        -> SetLineStyle(fLin);
  hPtDelta        -> SetLineWidth(fWid);
  hPtDelta        -> SetTitle(sTitle.Data());
  hPtDelta        -> SetTitleFont(fTxt);
  hPtDelta        -> GetXaxis() -> SetRangeUser(rDeltaRange[0], rDeltaRange[1]);
  hPtDelta        -> GetXaxis() -> SetTitle(sPtDeltaAxis.Data());
  hPtDelta        -> GetXaxis() -> SetTitleFont(fTxt);
  hPtDelta        -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPtDelta        -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPtDelta        -> GetXaxis() -> SetLabelFont(fTxt);
  hPtDelta        -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPtDelta        -> GetXaxis() -> CenterTitle(fCnt);
  hPtDelta        -> GetYaxis() -> SetTitle(sCounts.Data());
  hPtDelta        -> GetYaxis() -> SetTitleFont(fTxt);
  hPtDelta        -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPtDelta        -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPtDelta        -> GetYaxis() -> SetLabelFont(fTxt);
  hPtDelta        -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPtDelta        -> GetYaxis() -> CenterTitle(fCnt);
  hPtTrack        -> SetMarkerColor(fColTrk);
  hPtTrack        -> SetMarkerStyle(fMarTrk);
  hPtTrack        -> SetFillColor(fColTrk);
  hPtTrack        -> SetFillStyle(fFil);
  hPtTrack        -> SetLineColor(fColTrk);
  hPtTrack        -> SetLineStyle(fLin);
  hPtTrack        -> SetLineWidth(fWid);
  hPtTrack        -> SetTitle(sTitle.Data());
  hPtTrack        -> SetTitleFont(fTxt);
  hPtTrack        -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
  hPtTrack        -> GetXaxis() -> SetTitle(sPtRecoAxis.Data());
  hPtTrack        -> GetXaxis() -> SetTitleFont(fTxt);
  hPtTrack        -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPtTrack        -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPtTrack        -> GetXaxis() -> SetLabelFont(fTxt);
  hPtTrack        -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPtTrack        -> GetXaxis() -> CenterTitle(fCnt);
  hPtTrack        -> GetYaxis() -> SetTitle(sCounts.Data());
  hPtTrack        -> GetYaxis() -> SetTitleFont(fTxt);
  hPtTrack        -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPtTrack        -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPtTrack        -> GetYaxis() -> SetLabelFont(fTxt);
  hPtTrack        -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPtTrack        -> GetYaxis() -> CenterTitle(fCnt);
  hPtFrac         -> SetMarkerColor(fColTrk);
  hPtFrac         -> SetMarkerStyle(fMarTrk);
  hPtFrac         -> SetFillColor(fColTrk);
  hPtFrac         -> SetFillStyle(fFil);
  hPtFrac         -> SetLineColor(fColTrk);
  hPtFrac         -> SetLineStyle(fLin);
  hPtFrac         -> SetLineWidth(fWid);
  hPtFrac         -> SetTitle(sTitle.Data());
  hPtFrac         -> SetTitleFont(fTxt);
  hPtFrac         -> GetXaxis() -> SetRangeUser(rFracRange[0], rFracRange[1]);
  hPtFrac         -> GetXaxis() -> SetTitle(sPtFracAxis.Data());
  hPtFrac         -> GetXaxis() -> SetTitleFont(fTxt);
  hPtFrac         -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPtFrac         -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPtFrac         -> GetXaxis() -> SetLabelFont(fTxt);
  hPtFrac         -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPtFrac         -> GetXaxis() -> CenterTitle(fCnt);
  hPtFrac         -> GetYaxis() -> SetTitle(sCounts.Data());
  hPtFrac         -> GetYaxis() -> SetTitleFont(fTxt);
  hPtFrac         -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPtFrac         -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPtFrac         -> GetYaxis() -> SetLabelFont(fTxt);
  hPtFrac         -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPtFrac         -> GetYaxis() -> CenterTitle(fCnt);
  hPtTrkTru       -> SetMarkerColor(fColTrk);
  hPtTrkTru       -> SetMarkerStyle(fMarTrk);
  hPtTrkTru       -> SetFillColor(fColTrk);
  hPtTrkTru       -> SetFillStyle(fFil);
  hPtTrkTru       -> SetLineColor(fColTrk);
  hPtTrkTru       -> SetLineStyle(fLin);
  hPtTrkTru       -> SetLineWidth(fWid);
  hPtTrkTru       -> SetTitle(sTitle.Data());
  hPtTrkTru       -> SetTitleFont(fTxt);
  hPtTrkTru       -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
  hPtTrkTru       -> GetXaxis() -> SetTitle(sPtTrueAxis.Data());
  hPtTrkTru       -> GetXaxis() -> SetTitleFont(fTxt);
  hPtTrkTru       -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPtTrkTru       -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPtTrkTru       -> GetXaxis() -> SetLabelFont(fTxt);
  hPtTrkTru       -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPtTrkTru       -> GetXaxis() -> CenterTitle(fCnt);
  hPtTrkTru       -> GetYaxis() -> SetTitle(sCounts.Data());
  hPtTrkTru       -> GetYaxis() -> SetTitleFont(fTxt);
  hPtTrkTru       -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPtTrkTru       -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPtTrkTru       -> GetYaxis() -> SetLabelFont(fTxt);
  hPtTrkTru       -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPtTrkTru       -> GetYaxis() -> CenterTitle(fCnt);
  hPtDeltaVsFrac  -> SetMarkerColor(fColTrk);
  hPtDeltaVsFrac  -> SetMarkerStyle(fMarTrk);
  hPtDeltaVsFrac  -> SetFillColor(fColTrk);
  hPtDeltaVsFrac  -> SetFillStyle(fFil);
  hPtDeltaVsFrac  -> SetLineColor(fColTrk);
  hPtDeltaVsFrac  -> SetLineStyle(fLin);
  hPtDeltaVsFrac  -> SetLineWidth(fWid);
  hPtDeltaVsFrac  -> SetTitle(sTitle.Data());
  hPtDeltaVsFrac  -> SetTitleFont(fTxt);
  hPtDeltaVsFrac  -> GetXaxis() -> SetRangeUser(rFracRange[0], rFracRange[1]);
  hPtDeltaVsFrac  -> GetXaxis() -> SetTitle(sPtFracAxis.Data());
  hPtDeltaVsFrac  -> GetXaxis() -> SetTitleFont(fTxt);
  hPtDeltaVsFrac  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPtDeltaVsFrac  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPtDeltaVsFrac  -> GetXaxis() -> SetLabelFont(fTxt);
  hPtDeltaVsFrac  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPtDeltaVsFrac  -> GetXaxis() -> CenterTitle(fCnt);
  hPtDeltaVsFrac  -> GetYaxis() -> SetRangeUser(rDeltaRange[0], rDeltaRange[1]);
  hPtDeltaVsFrac  -> GetYaxis() -> SetTitle(sPtDeltaAxis.Data());
  hPtDeltaVsFrac  -> GetYaxis() -> SetTitleFont(fTxt);
  hPtDeltaVsFrac  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPtDeltaVsFrac  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPtDeltaVsFrac  -> GetYaxis() -> SetLabelFont(fTxt);
  hPtDeltaVsFrac  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPtDeltaVsFrac  -> GetYaxis() -> CenterTitle(fCnt);
  hPtDeltaVsFrac  -> GetZaxis() -> SetTitle(sCounts.Data());
  hPtDeltaVsFrac  -> GetZaxis() -> SetTitleFont(fTxt);
  hPtDeltaVsFrac  -> GetZaxis() -> SetTitleSize(fTit[1]);
  hPtDeltaVsFrac  -> GetZaxis() -> SetTitleOffset(fOffZ[1]);
  hPtDeltaVsFrac  -> GetZaxis() -> SetLabelFont(fTxt);
  hPtDeltaVsFrac  -> GetZaxis() -> SetLabelSize(fLab[1]);
  hPtDeltaVsFrac  -> GetZaxis() -> CenterTitle(fCnt);
  hPtDeltaVsTrue  -> SetMarkerColor(fColTrk);
  hPtDeltaVsTrue  -> SetMarkerStyle(fMarTrk);
  hPtDeltaVsTrue  -> SetFillColor(fColTrk);
  hPtDeltaVsTrue  -> SetFillStyle(fFil);
  hPtDeltaVsTrue  -> SetLineColor(fColTrk);
  hPtDeltaVsTrue  -> SetLineStyle(fLin);
  hPtDeltaVsTrue  -> SetLineWidth(fWid);
  hPtDeltaVsTrue  -> SetTitle(sTitle.Data());
  hPtDeltaVsTrue  -> SetTitleFont(fTxt);
  hPtDeltaVsTrue  -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
  hPtDeltaVsTrue  -> GetXaxis() -> SetTitle(sPtTrueAxis.Data());
  hPtDeltaVsTrue  -> GetXaxis() -> SetTitleFont(fTxt);
  hPtDeltaVsTrue  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPtDeltaVsTrue  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPtDeltaVsTrue  -> GetXaxis() -> SetLabelFont(fTxt);
  hPtDeltaVsTrue  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPtDeltaVsTrue  -> GetXaxis() -> CenterTitle(fCnt);
  hPtDeltaVsTrue  -> GetYaxis() -> SetRangeUser(rDeltaRange[0], rDeltaRange[1]);
  hPtDeltaVsTrue  -> GetYaxis() -> SetTitle(sPtDeltaAxis.Data());
  hPtDeltaVsTrue  -> GetYaxis() -> SetTitleFont(fTxt);
  hPtDeltaVsTrue  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPtDeltaVsTrue  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPtDeltaVsTrue  -> GetYaxis() -> SetLabelFont(fTxt);
  hPtDeltaVsTrue  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPtDeltaVsTrue  -> GetYaxis() -> CenterTitle(fCnt);
  hPtDeltaVsTrue  -> GetZaxis() -> SetTitle(sCounts.Data());
  hPtDeltaVsTrue  -> GetZaxis() -> SetTitleFont(fTxt);
  hPtDeltaVsTrue  -> GetZaxis() -> SetTitleSize(fTit[1]);
  hPtDeltaVsTrue  -> GetZaxis() -> SetTitleOffset(fOffZ[1]);
  hPtDeltaVsTrue  -> GetZaxis() -> SetLabelFont(fTxt);
  hPtDeltaVsTrue  -> GetZaxis() -> SetLabelSize(fLab[1]);
  hPtDeltaVsTrue  -> GetZaxis() -> CenterTitle(fCnt);
  hPtDeltaVsTrack -> SetMarkerColor(fColTrk);
  hPtDeltaVsTrack -> SetMarkerStyle(fMarTrk);
  hPtDeltaVsTrack -> SetFillColor(fColTrk);
  hPtDeltaVsTrack -> SetFillStyle(fFil);
  hPtDeltaVsTrack -> SetLineColor(fColTrk);
  hPtDeltaVsTrack -> SetLineStyle(fLin);
  hPtDeltaVsTrack -> SetLineWidth(fWid);
  hPtDeltaVsTrack -> SetTitle(sTitle.Data());
  hPtDeltaVsTrack -> SetTitleFont(fTxt);
  hPtDeltaVsTrack -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
  hPtDeltaVsTrack -> GetXaxis() -> SetTitle(sPtRecoAxis.Data());
  hPtDeltaVsTrack -> GetXaxis() -> SetTitleFont(fTxt);
  hPtDeltaVsTrack -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPtDeltaVsTrack -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPtDeltaVsTrack -> GetXaxis() -> SetLabelFont(fTxt);
  hPtDeltaVsTrack -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPtDeltaVsTrack -> GetXaxis() -> CenterTitle(fCnt);
  hPtDeltaVsTrack -> GetYaxis() -> SetRangeUser(rDeltaRange[0], rDeltaRange[1]);
  hPtDeltaVsTrack -> GetYaxis() -> SetTitle(sPtDeltaAxis.Data());
  hPtDeltaVsTrack -> GetYaxis() -> SetTitleFont(fTxt);
  hPtDeltaVsTrack -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPtDeltaVsTrack -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPtDeltaVsTrack -> GetYaxis() -> SetLabelFont(fTxt);
  hPtDeltaVsTrack -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPtDeltaVsTrack -> GetYaxis() -> CenterTitle(fCnt);
  hPtDeltaVsTrack -> GetZaxis() -> SetTitle(sCounts.Data());
  hPtDeltaVsTrack -> GetZaxis() -> SetTitleFont(fTxt);
  hPtDeltaVsTrack -> GetZaxis() -> SetTitleSize(fTit[1]);
  hPtDeltaVsTrack -> GetZaxis() -> SetTitleOffset(fOffZ[1]);
  hPtDeltaVsTrack -> GetZaxis() -> SetLabelFont(fTxt);
  hPtDeltaVsTrack -> GetZaxis() -> SetLabelSize(fLab[1]);
  hPtDeltaVsTrack -> GetZaxis() -> CenterTitle(fCnt);
  hPtTrueVsTrack  -> SetMarkerColor(fColTrk);
  hPtTrueVsTrack  -> SetMarkerStyle(fMarTrk);
  hPtTrueVsTrack  -> SetFillColor(fColTrk);
  hPtTrueVsTrack  -> SetFillStyle(fFil);
  hPtTrueVsTrack  -> SetLineColor(fColTrk);
  hPtTrueVsTrack  -> SetLineStyle(fLin);
  hPtTrueVsTrack  -> SetLineWidth(fWid);
  hPtTrueVsTrack  -> SetTitle(sTitle.Data());
  hPtTrueVsTrack  -> SetTitleFont(fTxt);
  hPtTrueVsTrack  -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
  hPtTrueVsTrack  -> GetXaxis() -> SetTitle(sPtRecoAxis.Data());
  hPtTrueVsTrack  -> GetXaxis() -> SetTitleFont(fTxt);
  hPtTrueVsTrack  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hPtTrueVsTrack  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hPtTrueVsTrack  -> GetXaxis() -> SetLabelFont(fTxt);
  hPtTrueVsTrack  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hPtTrueVsTrack  -> GetXaxis() -> CenterTitle(fCnt);
  hPtTrueVsTrack  -> GetYaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
  hPtTrueVsTrack  -> GetYaxis() -> SetTitle(sPtTrueAxis.Data());
  hPtTrueVsTrack  -> GetYaxis() -> SetTitleFont(fTxt);
  hPtTrueVsTrack  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hPtTrueVsTrack  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hPtTrueVsTrack  -> GetYaxis() -> SetLabelFont(fTxt);
  hPtTrueVsTrack  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hPtTrueVsTrack  -> GetYaxis() -> CenterTitle(fCnt);
  hPtTrueVsTrack  -> GetZaxis() -> SetTitle(sCounts.Data());
  hPtTrueVsTrack  -> GetZaxis() -> SetTitleFont(fTxt);
  hPtTrueVsTrack  -> GetZaxis() -> SetTitleSize(fTit[1]);
  hPtTrueVsTrack  -> GetZaxis() -> SetTitleOffset(fOffZ[1]);
  hPtTrueVsTrack  -> GetZaxis() -> SetLabelFont(fTxt);
  hPtTrueVsTrack  -> GetZaxis() -> SetLabelSize(fLab[1]);
  hPtTrueVsTrack  -> GetZaxis() -> CenterTitle(fCnt);
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    hEffCut[iCut]            -> SetMarkerColor(fColCut[iCut]);
    hEffCut[iCut]            -> SetMarkerStyle(fMarCut[iCut]);
    hEffCut[iCut]            -> SetFillColor(fColCut[iCut]);
    hEffCut[iCut]            -> SetFillStyle(fFil);
    hEffCut[iCut]            -> SetLineColor(fColCut[iCut]);
    hEffCut[iCut]            -> SetLineStyle(fLin);
    hEffCut[iCut]            -> SetLineWidth(fWid);
    hEffCut[iCut]            -> SetTitle(sTitle.Data());
    hEffCut[iCut]            -> SetTitleFont(fTxt);
    hEffCut[iCut]            -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
    hEffCut[iCut]            -> GetXaxis() -> SetTitle(sPtTrueAxis.Data());
    hEffCut[iCut]            -> GetXaxis() -> SetTitleFont(fTxt);
    hEffCut[iCut]            -> GetXaxis() -> SetTitleSize(fTit[0]);
    hEffCut[iCut]            -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hEffCut[iCut]            -> GetXaxis() -> SetLabelFont(fTxt);
    hEffCut[iCut]            -> GetXaxis() -> SetLabelSize(fLab[0]);
    hEffCut[iCut]            -> GetXaxis() -> CenterTitle(fCnt);
    hEffCut[iCut]            -> GetYaxis() -> SetTitle(sEffAxis.Data());
    hEffCut[iCut]            -> GetYaxis() -> SetTitleFont(fTxt);
    hEffCut[iCut]            -> GetYaxis() -> SetTitleSize(fTit[0]);
    hEffCut[iCut]            -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hEffCut[iCut]            -> GetYaxis() -> SetLabelFont(fTxt);
    hEffCut[iCut]            -> GetYaxis() -> SetLabelSize(fLab[0]);
    hEffCut[iCut]            -> GetYaxis() -> CenterTitle(fCnt);
    hPtDeltaCut[iCut]        -> SetMarkerColor(fColCut[iCut]);
    hPtDeltaCut[iCut]        -> SetMarkerStyle(fMarCut[iCut]);
    hPtDeltaCut[iCut]        -> SetFillColor(fColCut[iCut]);
    hPtDeltaCut[iCut]        -> SetFillStyle(fFil);
    hPtDeltaCut[iCut]        -> SetLineColor(fColCut[iCut]);
    hPtDeltaCut[iCut]        -> SetLineStyle(fLin);
    hPtDeltaCut[iCut]        -> SetLineWidth(fWid);
    hPtDeltaCut[iCut]        -> SetTitle(sTitle.Data());
    hPtDeltaCut[iCut]        -> SetTitleFont(fTxt);
    hPtDeltaCut[iCut]        -> GetXaxis() -> SetRangeUser(rDeltaRange[0], rDeltaRange[1]);
    hPtDeltaCut[iCut]        -> GetXaxis() -> SetTitle(sPtDeltaAxis.Data());
    hPtDeltaCut[iCut]        -> GetXaxis() -> SetTitleFont(fTxt);
    hPtDeltaCut[iCut]        -> GetXaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaCut[iCut]        -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hPtDeltaCut[iCut]        -> GetXaxis() -> SetLabelFont(fTxt);
    hPtDeltaCut[iCut]        -> GetXaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaCut[iCut]        -> GetXaxis() -> CenterTitle(fCnt);
    hPtDeltaCut[iCut]        -> GetYaxis() -> SetTitle(sCounts.Data());
    hPtDeltaCut[iCut]        -> GetYaxis() -> SetTitleFont(fTxt);
    hPtDeltaCut[iCut]        -> GetYaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaCut[iCut]        -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hPtDeltaCut[iCut]        -> GetYaxis() -> SetLabelFont(fTxt);
    hPtDeltaCut[iCut]        -> GetYaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaCut[iCut]        -> GetYaxis() -> CenterTitle(fCnt);
    hPtTrackCut[iCut]        -> SetMarkerColor(fColCut[iCut]);
    hPtTrackCut[iCut]        -> SetMarkerStyle(fMarCut[iCut]);
    hPtTrackCut[iCut]        -> SetFillColor(fColCut[iCut]);
    hPtTrackCut[iCut]        -> SetFillStyle(fFil);
    hPtTrackCut[iCut]        -> SetLineColor(fColCut[iCut]);
    hPtTrackCut[iCut]        -> SetLineStyle(fLin);
    hPtTrackCut[iCut]        -> SetLineWidth(fWid);
    hPtTrackCut[iCut]        -> SetTitle(sTitle.Data());
    hPtTrackCut[iCut]        -> SetTitleFont(fTxt);
    hPtTrackCut[iCut]        -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
    hPtTrackCut[iCut]        -> GetXaxis() -> SetTitle(sPtRecoAxis.Data());
    hPtTrackCut[iCut]        -> GetXaxis() -> SetTitleFont(fTxt);
    hPtTrackCut[iCut]        -> GetXaxis() -> SetTitleSize(fTit[1]);
    hPtTrackCut[iCut]        -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hPtTrackCut[iCut]        -> GetXaxis() -> SetLabelFont(fTxt);
    hPtTrackCut[iCut]        -> GetXaxis() -> SetLabelSize(fLab[1]);
    hPtTrackCut[iCut]        -> GetXaxis() -> CenterTitle(fCnt);
    hPtTrackCut[iCut]        -> GetYaxis() -> SetTitle(sCounts.Data());
    hPtTrackCut[iCut]        -> GetYaxis() -> SetTitleFont(fTxt);
    hPtTrackCut[iCut]        -> GetYaxis() -> SetTitleSize(fTit[1]);
    hPtTrackCut[iCut]        -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hPtTrackCut[iCut]        -> GetYaxis() -> SetLabelFont(fTxt);
    hPtTrackCut[iCut]        -> GetYaxis() -> SetLabelSize(fLab[1]);
    hPtTrackCut[iCut]        -> GetYaxis() -> CenterTitle(fCnt);
    hPtFracCut[iCut]         -> SetMarkerColor(fColCut[iCut]);
    hPtFracCut[iCut]         -> SetMarkerStyle(fMarCut[iCut]);
    hPtFracCut[iCut]         -> SetFillColor(fColCut[iCut]);
    hPtFracCut[iCut]         -> SetFillStyle(fFil);
    hPtFracCut[iCut]         -> SetLineColor(fColCut[iCut]);
    hPtFracCut[iCut]         -> SetLineStyle(fLin);
    hPtFracCut[iCut]         -> SetLineWidth(fWid);
    hPtFracCut[iCut]         -> SetTitle(sTitle.Data());
    hPtFracCut[iCut]         -> SetTitleFont(fTxt);
    hPtFracCut[iCut]         -> GetXaxis() -> SetRangeUser(rFracRange[0], rFracRange[1]);
    hPtFracCut[iCut]         -> GetXaxis() -> SetTitle(sPtFracAxis.Data());
    hPtFracCut[iCut]         -> GetXaxis() -> SetTitleFont(fTxt);
    hPtFracCut[iCut]         -> GetXaxis() -> SetTitleSize(fTit[1]);
    hPtFracCut[iCut]         -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hPtFracCut[iCut]         -> GetXaxis() -> SetLabelFont(fTxt);
    hPtFracCut[iCut]         -> GetXaxis() -> SetLabelSize(fLab[1]);
    hPtFracCut[iCut]         -> GetXaxis() -> CenterTitle(fCnt);
    hPtFracCut[iCut]         -> GetYaxis() -> SetTitle(sCounts.Data());
    hPtFracCut[iCut]         -> GetYaxis() -> SetTitleFont(fTxt);
    hPtFracCut[iCut]         -> GetYaxis() -> SetTitleSize(fTit[1]);
    hPtFracCut[iCut]         -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hPtFracCut[iCut]         -> GetYaxis() -> SetLabelFont(fTxt);
    hPtFracCut[iCut]         -> GetYaxis() -> SetLabelSize(fLab[1]);
    hPtFracCut[iCut]         -> GetYaxis() -> CenterTitle(fCnt);
    hPtTrkTruCut[iCut]       -> SetMarkerColor(fColCut[iCut]);
    hPtTrkTruCut[iCut]       -> SetMarkerStyle(fMarCut[iCut]);
    hPtTrkTruCut[iCut]       -> SetFillColor(fColCut[iCut]);
    hPtTrkTruCut[iCut]       -> SetFillStyle(fFil);
    hPtTrkTruCut[iCut]       -> SetLineColor(fColCut[iCut]);
    hPtTrkTruCut[iCut]       -> SetLineStyle(fLin);
    hPtTrkTruCut[iCut]       -> SetLineWidth(fWid);
    hPtTrkTruCut[iCut]       -> SetTitle(sTitle.Data());
    hPtTrkTruCut[iCut]       -> SetTitleFont(fTxt);
    hPtTrkTruCut[iCut]       -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
    hPtTrkTruCut[iCut]       -> GetXaxis() -> SetTitle(sPtTrueAxis.Data());
    hPtTrkTruCut[iCut]       -> GetXaxis() -> SetTitleFont(fTxt);
    hPtTrkTruCut[iCut]       -> GetXaxis() -> SetTitleSize(fTit[1]);
    hPtTrkTruCut[iCut]       -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hPtTrkTruCut[iCut]       -> GetXaxis() -> SetLabelFont(fTxt);
    hPtTrkTruCut[iCut]       -> GetXaxis() -> SetLabelSize(fLab[1]);
    hPtTrkTruCut[iCut]       -> GetXaxis() -> CenterTitle(fCnt);
    hPtTrkTruCut[iCut]       -> GetYaxis() -> SetTitle(sCounts.Data());
    hPtTrkTruCut[iCut]       -> GetYaxis() -> SetTitleFont(fTxt);
    hPtTrkTruCut[iCut]       -> GetYaxis() -> SetTitleSize(fTit[1]);
    hPtTrkTruCut[iCut]       -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hPtTrkTruCut[iCut]       -> GetYaxis() -> SetLabelFont(fTxt);
    hPtTrkTruCut[iCut]       -> GetYaxis() -> SetLabelSize(fLab[1]);
    hPtTrkTruCut[iCut]       -> GetYaxis() -> CenterTitle(fCnt);
    hPtDeltaVsFracCut[iCut]  -> SetMarkerColor(fColCut[iCut]);
    hPtDeltaVsFracCut[iCut]  -> SetMarkerStyle(fMarCut[iCut]);
    hPtDeltaVsFracCut[iCut]  -> SetFillColor(fColCut[iCut]);
    hPtDeltaVsFracCut[iCut]  -> SetFillStyle(fFil);
    hPtDeltaVsFracCut[iCut]  -> SetLineColor(fColCut[iCut]);
    hPtDeltaVsFracCut[iCut]  -> SetLineStyle(fLin);
    hPtDeltaVsFracCut[iCut]  -> SetLineWidth(fWid);
    hPtDeltaVsFracCut[iCut]  -> SetTitle(sTitle.Data());
    hPtDeltaVsFracCut[iCut]  -> SetTitleFont(fTxt);
    hPtDeltaVsFracCut[iCut]  -> GetXaxis() -> SetRangeUser(rFracRange[0], rFracRange[1]);
    hPtDeltaVsFracCut[iCut]  -> GetXaxis() -> SetTitle(sPtFracAxis.Data());
    hPtDeltaVsFracCut[iCut]  -> GetXaxis() -> SetTitleFont(fTxt);
    hPtDeltaVsFracCut[iCut]  -> GetXaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaVsFracCut[iCut]  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hPtDeltaVsFracCut[iCut]  -> GetXaxis() -> SetLabelFont(fTxt);
    hPtDeltaVsFracCut[iCut]  -> GetXaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaVsFracCut[iCut]  -> GetXaxis() -> CenterTitle(fCnt);
    hPtDeltaVsFracCut[iCut]  -> GetYaxis() -> SetRangeUser(rDeltaRange[0], rDeltaRange[1]);
    hPtDeltaVsFracCut[iCut]  -> GetYaxis() -> SetTitle(sPtDeltaAxis.Data());
    hPtDeltaVsFracCut[iCut]  -> GetYaxis() -> SetTitleFont(fTxt);
    hPtDeltaVsFracCut[iCut]  -> GetYaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaVsFracCut[iCut]  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hPtDeltaVsFracCut[iCut]  -> GetYaxis() -> SetLabelFont(fTxt);
    hPtDeltaVsFracCut[iCut]  -> GetYaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaVsFracCut[iCut]  -> GetYaxis() -> CenterTitle(fCnt);
    hPtDeltaVsFracCut[iCut]  -> GetZaxis() -> SetTitle(sCounts.Data());
    hPtDeltaVsFracCut[iCut]  -> GetZaxis() -> SetTitleFont(fTxt);
    hPtDeltaVsFracCut[iCut]  -> GetZaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaVsFracCut[iCut]  -> GetZaxis() -> SetTitleOffset(fOffZ[1]);
    hPtDeltaVsFracCut[iCut]  -> GetZaxis() -> SetLabelFont(fTxt);
    hPtDeltaVsFracCut[iCut]  -> GetZaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaVsFracCut[iCut]  -> GetZaxis() -> CenterTitle(fCnt);
    hPtDeltaVsTrueCut[iCut]  -> SetMarkerColor(fColTrk);
    hPtDeltaVsTrueCut[iCut]  -> SetMarkerStyle(fMarTrk);
    hPtDeltaVsTrueCut[iCut]  -> SetFillColor(fColTrk);
    hPtDeltaVsTrueCut[iCut]  -> SetFillStyle(fFil);
    hPtDeltaVsTrueCut[iCut]  -> SetLineColor(fColTrk);
    hPtDeltaVsTrueCut[iCut]  -> SetLineStyle(fLin);
    hPtDeltaVsTrueCut[iCut]  -> SetLineWidth(fWid);
    hPtDeltaVsTrueCut[iCut]  -> SetTitle(sTitle.Data());
    hPtDeltaVsTrueCut[iCut]  -> SetTitleFont(fTxt);
    hPtDeltaVsTrueCut[iCut]  -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetXaxis() -> SetTitle(sPtTrueAxis.Data());
    hPtDeltaVsTrueCut[iCut]  -> GetXaxis() -> SetTitleFont(fTxt);
    hPtDeltaVsTrueCut[iCut]  -> GetXaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetXaxis() -> SetLabelFont(fTxt);
    hPtDeltaVsTrueCut[iCut]  -> GetXaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetXaxis() -> CenterTitle(fCnt);
    hPtDeltaVsTrueCut[iCut]  -> GetYaxis() -> SetRangeUser(rDeltaRange[0], rDeltaRange[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetYaxis() -> SetTitle(sPtDeltaAxis.Data());
    hPtDeltaVsTrueCut[iCut]  -> GetYaxis() -> SetTitleFont(fTxt);
    hPtDeltaVsTrueCut[iCut]  -> GetYaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetYaxis() -> SetLabelFont(fTxt);
    hPtDeltaVsTrueCut[iCut]  -> GetYaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetYaxis() -> CenterTitle(fCnt);
    hPtDeltaVsTrueCut[iCut]  -> GetZaxis() -> SetTitle(sCounts.Data());
    hPtDeltaVsTrueCut[iCut]  -> GetZaxis() -> SetTitleFont(fTxt);
    hPtDeltaVsTrueCut[iCut]  -> GetZaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetZaxis() -> SetTitleOffset(fOffZ[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetZaxis() -> SetLabelFont(fTxt);
    hPtDeltaVsTrueCut[iCut]  -> GetZaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaVsTrueCut[iCut]  -> GetZaxis() -> CenterTitle(fCnt);
    hPtDeltaVsTrackCut[iCut] -> SetMarkerColor(fColTrk);
    hPtDeltaVsTrackCut[iCut] -> SetMarkerStyle(fMarTrk);
    hPtDeltaVsTrackCut[iCut] -> SetFillColor(fColTrk);
    hPtDeltaVsTrackCut[iCut] -> SetFillStyle(fFil);
    hPtDeltaVsTrackCut[iCut] -> SetLineColor(fColTrk);
    hPtDeltaVsTrackCut[iCut] -> SetLineStyle(fLin);
    hPtDeltaVsTrackCut[iCut] -> SetLineWidth(fWid);
    hPtDeltaVsTrackCut[iCut] -> SetTitle(sTitle.Data());
    hPtDeltaVsTrackCut[iCut] -> SetTitleFont(fTxt);
    hPtDeltaVsTrackCut[iCut] -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
    hPtDeltaVsTrackCut[iCut] -> GetXaxis() -> SetTitle(sPtRecoAxis.Data());
    hPtDeltaVsTrackCut[iCut] -> GetXaxis() -> SetTitleFont(fTxt);
    hPtDeltaVsTrackCut[iCut] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaVsTrackCut[iCut] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hPtDeltaVsTrackCut[iCut] -> GetXaxis() -> SetLabelFont(fTxt);
    hPtDeltaVsTrackCut[iCut] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaVsTrackCut[iCut] -> GetXaxis() -> CenterTitle(fCnt);
    hPtDeltaVsTrackCut[iCut] -> GetYaxis() -> SetRangeUser(rDeltaRange[0], rDeltaRange[1]);
    hPtDeltaVsTrackCut[iCut] -> GetYaxis() -> SetTitle(sPtDeltaAxis.Data());
    hPtDeltaVsTrackCut[iCut] -> GetYaxis() -> SetTitleFont(fTxt);
    hPtDeltaVsTrackCut[iCut] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaVsTrackCut[iCut] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hPtDeltaVsTrackCut[iCut] -> GetYaxis() -> SetLabelFont(fTxt);
    hPtDeltaVsTrackCut[iCut] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaVsTrackCut[iCut] -> GetYaxis() -> CenterTitle(fCnt);
    hPtDeltaVsTrackCut[iCut] -> GetZaxis() -> SetTitle(sCounts.Data());
    hPtDeltaVsTrackCut[iCut] -> GetZaxis() -> SetTitleFont(fTxt);
    hPtDeltaVsTrackCut[iCut] -> GetZaxis() -> SetTitleSize(fTit[1]);
    hPtDeltaVsTrackCut[iCut] -> GetZaxis() -> SetTitleOffset(fOffZ[1]);
    hPtDeltaVsTrackCut[iCut] -> GetZaxis() -> SetLabelFont(fTxt);
    hPtDeltaVsTrackCut[iCut] -> GetZaxis() -> SetLabelSize(fLab[1]);
    hPtDeltaVsTrackCut[iCut] -> GetZaxis() -> CenterTitle(fCnt);
    hPtTrueVsTrackCut[iCut]  -> SetMarkerColor(fColCut[iCut]);
    hPtTrueVsTrackCut[iCut]  -> SetMarkerStyle(fMarCut[iCut]);
    hPtTrueVsTrackCut[iCut]  -> SetFillColor(fColCut[iCut]);
    hPtTrueVsTrackCut[iCut]  -> SetFillStyle(fFil);
    hPtTrueVsTrackCut[iCut]  -> SetLineColor(fColCut[iCut]);
    hPtTrueVsTrackCut[iCut]  -> SetLineStyle(fLin);
    hPtTrueVsTrackCut[iCut]  -> SetLineWidth(fWid);
    hPtTrueVsTrackCut[iCut]  -> SetTitle(sTitle.Data());
    hPtTrueVsTrackCut[iCut]  -> SetTitleFont(fTxt);
    hPtTrueVsTrackCut[iCut]  -> GetXaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
    hPtTrueVsTrackCut[iCut]  -> GetXaxis() -> SetTitle(sPtRecoAxis.Data());
    hPtTrueVsTrackCut[iCut]  -> GetXaxis() -> SetTitleFont(fTxt);
    hPtTrueVsTrackCut[iCut]  -> GetXaxis() -> SetTitleSize(fTit[1]);
    hPtTrueVsTrackCut[iCut]  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hPtTrueVsTrackCut[iCut]  -> GetXaxis() -> SetLabelFont(fTxt);
    hPtTrueVsTrackCut[iCut]  -> GetXaxis() -> SetLabelSize(fLab[1]);
    hPtTrueVsTrackCut[iCut]  -> GetXaxis() -> CenterTitle(fCnt);
    hPtTrueVsTrackCut[iCut]  -> GetYaxis() -> SetRangeUser(rPtRange[0], rPtRange[1]);
    hPtTrueVsTrackCut[iCut]  -> GetYaxis() -> SetTitle(sPtTrueAxis.Data());
    hPtTrueVsTrackCut[iCut]  -> GetYaxis() -> SetTitleFont(fTxt);
    hPtTrueVsTrackCut[iCut]  -> GetYaxis() -> SetTitleSize(fTit[1]);
    hPtTrueVsTrackCut[iCut]  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hPtTrueVsTrackCut[iCut]  -> GetYaxis() -> SetLabelFont(fTxt);
    hPtTrueVsTrackCut[iCut]  -> GetYaxis() -> SetLabelSize(fLab[1]);
    hPtTrueVsTrackCut[iCut]  -> GetYaxis() -> CenterTitle(fCnt);
    hPtTrueVsTrackCut[iCut]  -> GetZaxis() -> SetTitle(sCounts.Data());
    hPtTrueVsTrackCut[iCut]  -> GetZaxis() -> SetTitleFont(fTxt);
    hPtTrueVsTrackCut[iCut]  -> GetZaxis() -> SetTitleSize(fTit[1]);
    hPtTrueVsTrackCut[iCut]  -> GetZaxis() -> SetTitleOffset(fOffZ[1]);
    hPtTrueVsTrackCut[iCut]  -> GetZaxis() -> SetLabelFont(fTxt);
    hPtTrueVsTrackCut[iCut]  -> GetZaxis() -> SetLabelSize(fLab[1]);
    hPtTrueVsTrackCut[iCut]  -> GetZaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLe       = 0;
  const UInt_t  fFilLe       = 0;
  const UInt_t  fLinLe       = 0;
  const Float_t yObjLe       = 0.1 + ((NCuts + 2) * 0.05);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, yObjLe};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hPtTruth,  sLegTrue.Data(),  "pf");
  leg -> AddEntry(hPtTrkTru, sLegTrack.Data(), "pf");
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    leg -> AddEntry(hPtTrkTruCut[iCut], sLegCut[iCut].Data(), "pf");
  }
  cout << "    Made legend." << endl;

  // make text boxes
  const UInt_t  fColInf       = 0;
  const UInt_t  fFilInf       = 0;
  const UInt_t  fLinInf       = 0;
  const Float_t yObjInf       = 0.1 + (NTxt * 0.05);
  const Float_t yObjCut       = 0.1 + (NTrkCuts * 0.05);
  const Float_t fInfXY[NVtx] = {0.3, 0.1, 0.5, yObjInf};
  const Float_t fCutXY[NVtx] = {0.5, 0.1, 0.7, yObjCut};

  TPaveText *info = new TPaveText(fInfXY[0], fInfXY[1], fInfXY[2], fInfXY[3], "NDC NB");
  info -> SetFillColor(fColInf);
  info -> SetFillStyle(fFilInf);
  info -> SetLineColor(fColInf);
  info -> SetLineStyle(fLinInf);
  info -> SetTextFont(fTxt);
  info -> SetTextAlign(fAln);
  for (Ssiz_t iTxt = 0; iTxt < NTxt; iTxt++) {
    info -> AddText(sInfo[iTxt].Data());
  }

  TPaveText *cuts = new TPaveText(fCutXY[0], fCutXY[1], fCutXY[2], fCutXY[3], "NDC NB");
  cuts -> SetFillColor(fColInf);
  cuts -> SetFillStyle(fFilInf);
  cuts -> SetLineColor(fColInf);
  cuts -> SetLineStyle(fLinInf);
  cuts -> SetTextFont(fTxt);
  cuts -> SetTextAlign(fAln);
  for (Ssiz_t iTrkCut = 0; iTrkCut < NTrkCuts; iTrkCut++) {
    cuts -> AddText(sTrkCuts[iTrkCut].Data());
  }
  cout << "    Made text." << endl;

  // make line
  const UInt_t  fColLi       = 1;
  const UInt_t  fLinLi       = 9;
  const UInt_t  fWidLi       = 1;
  const Float_t fLinXY[NVtx] = {rPtRange[0], 1., rPtRange[1], 1.};

  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);
  cout << "    Made line." << endl;

  // make plots
  const UInt_t  width(750);
  const UInt_t  width2D(1500);
  const UInt_t  height(950);
  const UInt_t  heightNR(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY1(1);
  const UInt_t  fLogY2(1);
  const UInt_t  fLogYNR(0);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT1(0.005);
  const Float_t fMarginT2(0.02);
  const Float_t fMarginTNR(0.02);
  const Float_t fMarginB1(0.25);
  const Float_t fMarginB2(0.005);
  const Float_t fMarginBNR(0.15);
  const Float_t fEffXY[NVtx]       = {0.,  0.,   1.,  0.35};
  const Float_t fTrksXY[NVtx]      = {0.,  0.35, 1.,  1.};
  const Float_t fBeforeDPtXY[NVtx] = {0.,  0.,   0.5, 1.};
  const Float_t fAfterDPtXY[NVtx]  = {0.5, 0.,   1.,  1.};

  TCanvas *cEff  = new TCanvas("cEfficiency", "", width, height);
  TPad    *pEff  = new TPad("pEff",  "", fEffXY[0],  fEffXY[1],  fEffXY[2],  fEffXY[3]);
  TPad    *pTrks = new TPad("pTrks", "", fTrksXY[0], fTrksXY[1], fTrksXY[2], fTrksXY[3]);
  cEff  -> SetGrid(fGrid, fGrid);
  cEff  -> SetTicks(fTick, fTick);
  cEff  -> SetBorderMode(fMode);
  cEff  -> SetBorderSize(fBord);
  pEff  -> SetGrid(fGrid, fGrid);
  pEff  -> SetTicks(fTick, fTick);
  pEff  -> SetLogx(fLogX);
  pEff  -> SetLogy(fLogY1);
  pEff  -> SetBorderMode(fMode);
  pEff  -> SetBorderSize(fBord);
  pEff  -> SetFrameBorderMode(fFrame);
  pEff  -> SetLeftMargin(fMarginL);
  pEff  -> SetRightMargin(fMarginR);
  pEff  -> SetTopMargin(fMarginT1);
  pEff  -> SetBottomMargin(fMarginB1);
  pTrks -> SetGrid(fGrid, fGrid);
  pTrks -> SetTicks(fTick, fTick);
  pTrks -> SetLogx(fLogX);
  pTrks -> SetLogy(fLogY2);
  pTrks -> SetBorderMode(fMode);
  pTrks -> SetBorderSize(fBord);
  pTrks -> SetFrameBorderMode(fFrame);
  pTrks -> SetLeftMargin(fMarginL);
  pTrks -> SetRightMargin(fMarginR);
  pTrks -> SetTopMargin(fMarginT2);
  pTrks -> SetBottomMargin(fMarginB2);
  cEff  -> cd();
  pEff  -> Draw();
  pTrks -> Draw();
  pEff  -> cd();
  hEff  -> Draw();
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    hEffCut[iCut] -> Draw("SAME");
  }
  line      -> Draw();
  pTrks     -> cd();
  hPtTruth  -> Draw();
  hPtTrkTru -> Draw("SAME");
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    hPtTrkTruCut[iCut] -> Draw("SAME");
  }
  leg     -> Draw();
  info    -> Draw();
  cuts    -> Draw();
  fOutput -> cd();
  cEff    -> Write();
  cEff    -> Close();

  TCanvas *cPtTruVsTrk = new TCanvas("cPtTruthVsTrack", "", width2D, heightNR);
  TPad    *pBefore     = new TPad("pBeforeDPt", "", fBeforeDPtXY[0], fBeforeDPtXY[1], fBeforeDPtXY[2], fBeforeDPtXY[3]);
  TPad    *pAfter      = new TPad("pAfterDPt", "",  fAfterDPtXY[0],  fAfterDPtXY[1],  fAfterDPtXY[2],  fAfterDPtXY[3]);
  cPtTruVsTrk                  -> SetGrid(fGrid, fGrid);
  cPtTruVsTrk                  -> SetTicks(fTick, fTick);
  cPtTruVsTrk                  -> SetBorderMode(fMode);
  cPtTruVsTrk                  -> SetBorderSize(fBord);
  pBefore                      -> SetGrid(fGrid, fGrid);
  pBefore                      -> SetTicks(fTick, fTick);
  pBefore                      -> SetLogx(fLogX);
  pBefore                      -> SetLogy(fLogYNR);
  pBefore                      -> SetBorderMode(fMode);
  pBefore                      -> SetBorderSize(fBord);
  pBefore                      -> SetFrameBorderMode(fFrame);
  pAfter                       -> SetGrid(fGrid, fGrid);
  pAfter                       -> SetTicks(fTick, fTick);
  pAfter                       -> SetLogx(fLogX);
  pAfter                       -> SetLogy(fLogYNR);
  pAfter                       -> SetBorderMode(fMode);
  pAfter                       -> SetBorderSize(fBord);
  pAfter                       -> SetFrameBorderMode(fFrame);
  cPtTruVsTrk                  -> cd();
  pBefore                      -> Draw();
  pAfter                       -> Draw();
  pBefore                      -> cd();
  hPtTrueVsTrack               -> SetTitle("Before #Deltap_{T}/p_{T} cut");
  hPtTrueVsTrack               -> Draw("colz");
  cuts                         -> Draw();
  pAfter                       -> cd();
  hPtTrueVsTrackCut[NCuts - 3] -> SetTitle("After #Deltap_{T}/p_{T} < 0.03 cut");
  hPtTrueVsTrackCut[NCuts - 3] -> Draw("colz");
  info                         -> Draw();
  fOutput                      -> cd();
  cPtTruVsTrk                  -> Write();
  cPtTruVsTrk                  -> Close();

  TCanvas *cReject = new TCanvas("cReject", "", width, heightNR);
  cReject  -> SetGrid(fGrid, fGrid);
  cReject  -> SetTicks(fTick, fTick);
  cReject  -> SetBorderMode(fMode);
  cReject  -> SetBorderSize(fBord);
  cReject  -> SetFrameBorderMode(fFrame);
  cReject  -> SetLeftMargin(fMarginL);
  cReject  -> SetRightMargin(fMarginR);
  cReject  -> SetTopMargin(fMarginTNR);
  cReject  -> SetBottomMargin(fMarginBNR);
  cReject  -> SetLogx(fLogX);
  cReject  -> SetLogy(fLogYNR);
  cReject  -> cd();
  grReject -> Draw("ALP");
  info     -> Draw();
  cuts     -> Draw();
  fOutput  -> cd();
  cReject  -> Write();
  cReject  -> Close();

  TCanvas *cDeltaPt = new TCanvas("cDeltaPt", "", width, heightNR);
  cDeltaPt  -> SetGrid(fGrid, fGrid);
  cDeltaPt  -> SetTicks(fTick, fTick);
  cDeltaPt  -> SetBorderMode(fMode);
  cDeltaPt  -> SetBorderSize(fBord);
  cDeltaPt  -> SetFrameBorderMode(fFrame);
  cDeltaPt  -> SetLeftMargin(fMarginL);
  cDeltaPt  -> SetRightMargin(fMarginR);
  cDeltaPt  -> SetTopMargin(fMarginTNR);
  cDeltaPt  -> SetBottomMargin(fMarginBNR);
  cDeltaPt  -> SetLogx(fLogX);
  cDeltaPt  -> SetLogy(fLogYNR);
  cDeltaPt  -> cd();
  hPtDelta  -> Draw();
  info      -> Draw();
  cuts      -> Draw();
  fOutput   -> cd();
  cDeltaPt  -> Write();
  cDeltaPt  -> Close();
  cout << "    Made plots." << endl;

  // save histograms
  fOutput         -> cd();
  grReject        -> Write();
  hEff            -> Write();
  hPtTruth        -> Write();
  hPtDelta        -> Write();
  hPtTrack        -> Write();
  hPtFrac         -> Write();
  hPtTrkTru       -> Write();
  hPtDeltaVsFrac  -> Write();
  hPtDeltaVsTrue  -> Write();
  hPtDeltaVsTrack -> Write();
  hPtTrueVsTrack  -> Write();
  for (Ssiz_t iCut = 0; iCut < NCuts; iCut++) {
    hEffCut[iCut]            -> Write();
    hPtDeltaCut[iCut]        -> Write();
    hPtTrackCut[iCut]        -> Write();
    hPtFracCut[iCut]         -> Write();
    hPtTrkTruCut[iCut]       -> Write();
    hPtDeltaVsFracCut[iCut]  -> Write();
    hPtDeltaVsTrueCut[iCut]  -> Write();
    hPtDeltaVsTrackCut[iCut] -> Write();
    hPtTrueVsTrackCut[iCut]  -> Write();
  }

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Finished delta-pt extractor script!\n" << endl;

}

// end ------------------------------------------------------------------------
