// ----------------------------------------------------------------------------
// 'STrackCutStudy.h'
// Derek Anderson
// 12.15.2022
//
// Reads in the 'ntp_track' Ntuple
// generated by the SVtxEvaluator
// class and studies the impact
// of cutting on various quantities.
// ----------------------------------------------------------------------------

#ifndef STRACKCUTSTUDY_H
#define STRACKCUTSTUDY_H

// standard c includes
#include <cmath>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <utility>
#include <iostream>
// root includes
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TFile.h>
#include <TMath.h>
#include <TError.h>
#include <TNtuple.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TDirectory.h>

using namespace std;

// global constants
static const Ssiz_t NVtx(4);
static const Ssiz_t NType(9);
static const Ssiz_t NTrkVar(12);
static const Ssiz_t NPhysVar(6);
static const Ssiz_t NRange(2);
static const Ssiz_t NPanel(2);
static const UInt_t FTxt(42);



class STrackCutStudy {

  public:

    // enums
    enum TRKVAR {
      VX       = 0,
      VY       = 1,
      VZ       = 2,
      NMMS     = 3,
      NMAP     = 4,
      NINT     = 5,
      NTPC     = 6,
      QUAL     = 7,
      DCAXY    = 8,
      DCAZ     = 9,
      DELDCAXY = 10,
      DELDCAZ  = 11
    };
    enum PHYSVAR {
      PHI    = 0,
      ETA    = 1,
      PT     = 2,
      DELPHI = 3,
      DELETA = 4,
      DELPT  = 5
    };
    enum TYPE {
      TRACK     = 0,
      TRUTH     = 1,
      WEIRD_ALL = 2,
      WEIRD_SI  = 3,
      WEIRD_TPC = 4,
      NORMAL    = 5,
      PILEUP    = 6,
      PRIMARY   = 7,
      NONPRIM   = 8
    };

    // ctor/dtor
    STrackCutStudy();
    ~STrackCutStudy();

    // public methods
    void SetInputOutputFiles(const TString sEmbedOnlyInput, const TString sPileupInput, const TString sOutput);
    void SetInputTuples(const TString sEmbedOnlyTuple, const TString sPileupTuple);
    void SetStudyParameters(const Bool_t intNorm, const Double_t weirdFracMin, const Double_t weirdFracMax);
    void SetCutFlags(const Bool_t doPrimary, const Bool_t doMVtx, const Bool_t doVz, const Bool_t doDcaXY, const Bool_t doDcaZ, const Bool_t doQuality);
    void SetTrackCuts(const pair<UInt_t, UInt_t> nMVtxRange, const pair<Double_t, Double_t> vzRange, const pair<Double_t, Double_t> dcaXyRange, const pair <Double_t, Double_t> dcaZRange, const pair<Double_t, Double_t> qualityRange);
    void SetPlotText(const Ssiz_t nTxtE, const Ssiz_t nTxtP, const TString sTxtE[], const TString sTxtP[]);
    void Init();
    void Analyze();
    void End();

  private:

    // track type/variable names/labels
    const Bool_t  isPileup[NType]     = {false,   false,   false,      false,     false,      false,    true,        true,          true};
    const TString sTrkNames[NType]    = {"Track", "Truth", "AllWeird", "SiWeird", "TpcWeird", "Normal", "AllPileup", "PrimePileup", "NonPrimePileup"};
    const TString sTrkLabels[NType]   = {"All tracks", "Truth tracks", "Weird tracks (all)", "Weird tracks (Si seed)", "Weird tracks (TPC seed)",
                                         "Normal tracks", "Including pileup tracks (all)", "Including pileup tracks (only primary)", "Including pileup gracks (non-primary)"};
    const TString sTrkVars[NTrkVar]   = {"Vx", "Vy", "Vz", "NMms", "NMap", "NInt", "NTpc", "Qual", "DcaXY", "DcaZ", "DeltaDcaXY", "DeltaDcaZ"};
    const TString sPhysVars[NPhysVar] = {"Phi", "Eta", "Pt", "DeltaPhi", "DeltaEta", "DeltaPt"};

    // io members
    TFile   *fOut;
    TFile   *fInEO;
    TFile   *fInPU;
    TString  sInFileEO;
    TString  sInFilePU;
    TString  sInTupleEO;
    TString  sInTuplePU;
    TString  sOutfile;
    TNtuple *ntTrkEO;
    TNtuple *ntTrkPU;

    // track-variable histograms
    TH1D *hTrkVar[NType][NTrkVar];
    TH1D *hTrkVarDiff[NType][NTrkVar];
    TH1D *hTrkVarFrac[NType][NTrkVar];
    TH2D *hTrkVarVsNTpc[NType][NTrkVar];
    TH2D *hTrkVarVsPtReco[NType][NTrkVar];
    TH2D *hTrkVarVsPtTrue[NType][NTrkVar];
    TH2D *hTrkVarVsPtFrac[NType][NTrkVar];

    // physics-variable histograms
    TH1D *hPhysVar[NType][NPhysVar];
    TH1D *hPhysVarDiff[NType][NPhysVar];
    TH1D *hPhysVarFrac[NType][NPhysVar];
    TH2D *hPhysVarVsNTpc[NType][NPhysVar];
    TH2D *hPhysVarVsPtReco[NType][NPhysVar];
    TH2D *hPhysVarVsPtTrue[NType][NPhysVar];
    TH2D *hPhysVarVsPtFrac[NType][NPhysVar];

    // text parameters
    Ssiz_t           nTxtEO;
    Ssiz_t           nTxtPU;
    TPaveText       *ptCut;
    vector<TString>  sTxtEO;
    vector<TString>  sTxtPU;

    // study parameters
    Bool_t   doIntNorm;
    Double_t normalPtFracMin;
    Double_t normalPtFracMax;

    // track cuts
    Bool_t doPrimaryCut;
    Bool_t doMVtxCut;
    Bool_t doVzCut;
    Bool_t doDcaXyCut;
    Bool_t doDcaZCut;
    Bool_t doQualityCut;
    pair<Double_t, Double_t> nMVtxCut;
    pair<Double_t, Double_t> vzCut;
    pair<Double_t, Double_t> dcaXyCut;
    pair<Double_t, Double_t> dcaZCut;
    pair<Double_t, Double_t> qualityCut;

    // embed-only leaves
    Float_t event;
    Float_t seed;
    Float_t trackID;
    Float_t crossing;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t pt;
    Float_t eta;
    Float_t phi;
    Float_t deltapt;
    Float_t deltaeta;
    Float_t deltaphi;
    Float_t charge;
    Float_t quality;
    Float_t chisq;
    Float_t ndf;
    Float_t nhits;
    Float_t nmaps;
    Float_t nintt;
    Float_t ntpc;
    Float_t nmms;
    Float_t ntpc1;
    Float_t ntpc11;
    Float_t ntpc2;
    Float_t ntpc3;
    Float_t nlmaps;
    Float_t nlintt;
    Float_t nltpc;
    Float_t nlmms;
    Float_t layers;
    Float_t vertexID;
    Float_t vx;
    Float_t vy;
    Float_t vz;
    Float_t dca2d;
    Float_t dca2dsigma;
    Float_t dca3dxy;
    Float_t dca3dxysigma;
    Float_t dca3dz;
    Float_t dca3dzsigma;
    Float_t pcax;
    Float_t pcay;
    Float_t pcaz;
    Float_t gtrackID;
    Float_t gflavor;
    Float_t gnhits;
    Float_t gnmaps;
    Float_t gnintt;
    Float_t gntpc;
    Float_t gnmms;
    Float_t gnlmaps;
    Float_t gnlintt;
    Float_t gnltpc;
    Float_t gnlmms;
    Float_t gpx;
    Float_t gpy;
    Float_t gpz;
    Float_t gpt;
    Float_t geta;
    Float_t gphi;
    Float_t gvx;
    Float_t gvy;
    Float_t gvz;
    Float_t gvt;
    Float_t gfpx;
    Float_t gfpy;
    Float_t gfpz;
    Float_t gfx;
    Float_t gfy;
    Float_t gfz;
    Float_t gembed;
    Float_t gprimary;
    Float_t nfromtruth;
    Float_t nwrong;
    Float_t ntrumaps;
    Float_t ntruintt;
    Float_t ntrutpc;
    Float_t ntrumms;
    Float_t ntrutpc1;
    Float_t ntrutpc11;
    Float_t ntrutpc2;
    Float_t ntrutpc3;
    Float_t layersfromtruth;
    Float_t nhittpcall;
    Float_t nhittpcin;
    Float_t nhittpcmid;
    Float_t nhittpcout;
    Float_t nclusall;
    Float_t nclustpc;
    Float_t nclusintt;
    Float_t nclusmaps;
    Float_t nclusmms;

    // with-pileup leaves
    Float_t pu_event;
    Float_t pu_seed;
    Float_t pu_gntracks;
    Float_t pu_gtrackID;
    Float_t pu_gflavor;
    Float_t pu_gnhits;
    Float_t pu_gnmaps;
    Float_t pu_gnintt;
    Float_t pu_gnmms;
    Float_t pu_gnintt1;
    Float_t pu_gnintt2;
    Float_t pu_gnintt3;
    Float_t pu_gnintt4;
    Float_t pu_gnintt5;
    Float_t pu_gnintt6;
    Float_t pu_gnintt7;
    Float_t pu_gnintt8;
    Float_t pu_gntpc;
    Float_t pu_gnlmaps;
    Float_t pu_gnlintt;
    Float_t pu_gnltpc;
    Float_t pu_gnlmms;
    Float_t pu_gpx;
    Float_t pu_gpy;
    Float_t pu_gpz;
    Float_t pu_gpt;
    Float_t pu_geta;
    Float_t pu_gphi;
    Float_t pu_gvx;
    Float_t pu_gvy;
    Float_t pu_gvz;
    Float_t pu_gvt;
    Float_t pu_gfpx;
    Float_t pu_gfpy;
    Float_t pu_gfpz;
    Float_t pu_gfx;
    Float_t pu_gfy;
    Float_t pu_gfz;
    Float_t pu_gembed;
    Float_t pu_gprimary;
    Float_t pu_trackID;
    Float_t pu_px;
    Float_t pu_py;
    Float_t pu_pz;
    Float_t pu_pt;
    Float_t pu_eta;
    Float_t pu_phi;
    Float_t pu_deltapt;
    Float_t pu_deltaeta;
    Float_t pu_deltaphi;
    Float_t pu_charge;
    Float_t pu_quality;
    Float_t pu_chisq;
    Float_t pu_ndf;
    Float_t pu_nhits;
    Float_t pu_layers;
    Float_t pu_nmaps;
    Float_t pu_nintt;
    Float_t pu_ntpc;
    Float_t pu_nmms;
    Float_t pu_ntpc1;
    Float_t pu_ntpc11;
    Float_t pu_ntpc2;
    Float_t pu_ntpc3;
    Float_t pu_nlmaps;
    Float_t pu_nlintt;
    Float_t pu_nltpc;
    Float_t pu_nlmms;
    Float_t pu_vertexID;
    Float_t pu_vx;
    Float_t pu_vy;
    Float_t pu_vz;
    Float_t pu_dca2d;
    Float_t pu_dca2dsigma;
    Float_t pu_dca3dxy;
    Float_t pu_dca3dxysigma;
    Float_t pu_dca3dz;
    Float_t pu_dca3dzsigma;
    Float_t pu_pcax;
    Float_t pu_pcay;
    Float_t pu_pcaz;
    Float_t pu_nfromtruth;
    Float_t pu_nwrong;
    Float_t pu_ntrumaps;
    Float_t pu_ntruintt;
    Float_t pu_ntrutpc;
    Float_t pu_ntrumms;
    Float_t pu_ntrutpc1;
    Float_t pu_ntrutpc11;
    Float_t pu_ntrutpc2;
    Float_t pu_ntrutpc3;
    Float_t pu_layersfromtruth;
    Float_t pu_nhittpcall;
    Float_t pu_nhittpcin;
    Float_t pu_nhittpcmid;
    Float_t pu_nhittpcout;
    Float_t pu_nclusall;
    Float_t pu_nclustpc;
    Float_t pu_nclusintt;
    Float_t pu_nclusmaps;
    Float_t pu_nclusmms;

    // private methods
    void   InitFiles();
    void   InitTuples();
    void   InitHists();
    void   MakeCutText();
    void   NormalizeHists();
    void   SetHistStyles();
    void   CreatePlots();
    void   SaveHists();
    void   FillTrackHistograms(const Int_t type, const Double_t recoTrkVars[], const Double_t trueTrkVars[], const Double_t recoPhysVars[], const Double_t truePhysVars[]);
    void   FillTruthHistograms(const Int_t type, const Double_t recoTrkVars[], const Double_t trueTrkVars[], const Double_t recoPhysVars[], const Double_t truePhysVars[]);
    void   ConstructPlots(const Ssiz_t nToDraw, const Int_t typesToDraw[], const TString sDirToSaveTo, const TString sPlotLabel);
    Bool_t ApplyCuts(const Bool_t isPrimary, const UInt_t trkNMVtx, const Double_t trkVz, const Double_t trkDcaXY, const Double_t trkDcaZ, const Double_t trkQuality);

};  // end STrackCutStudy definition

#endif

// end ------------------------------------------------------------------------
