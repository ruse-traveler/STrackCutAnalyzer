// 'STrackCutStudy.h'
// Derek Anderson
// 12.15.2022
//
// Reads in the 'ntp_track' Ntuple
// generated by the SVtxEvaluator
// class and studies the impact
// of cutting on various quantities.

#ifndef STRACKCUTSTUDY_H
#define STRACKCUTSTUDY_H

// standard c includes
#include <cmath>
#include <vector>
#include <cassert>
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
static const Ssiz_t NTxt(3);
static const Ssiz_t NTyp(4);
static const Ssiz_t NVar(18);
static const Ssiz_t NRange(2);
static const Ssiz_t NPanel(2);
static const UInt_t FTxt(42);



class STrackCutStudy {

  public:

    // enums
    // TODO: streamline histogram declarations, styling, and saving
    enum VARIABLE {
      NMMS  = 0,
      NMAP  = 1,
      NINT  = 2,
      NTPC  = 3,
      NTOT  = 4,
      CHISQ = 5,
      NDF   = 6,
      QUAL  = 7,
      DCAT  = 8,
      DCAZ  = 9,
      PHI   = 10,
      ETA   = 11,
      PT    = 12,
      DDCAT = 13,
      DDCAZ = 14,
      DPHI  = 15,
      DETA  = 16,
      DPT   = 17
    };
    enum TYPE {
      TRACK  = 0,
      TRUTH  = 1,
      WEIRD  = 2,
      PILEUP = 3
    };

    // ctor/dtor
    STrackCutStudy();
    ~STrackCutStudy();

    // public methods
    void SetInputOutputFiles(const TString sEmbedOnlyInput, const TString sPileupInput, const TString sOutput);
    void SetInputTuples(const TString sEmbedOnlyTuple, const TString sPileupTuple);
    void SetWeirdFractionCuts(const Double_t weirdFracMin, const Double_t weirdFracMax);
    void Init();
    void Analyze();
    void End();

  private:

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

    // embed-only track output histograms
    TH1D *hTrackNMms;
    TH1D *hTrackNMap;
    TH1D *hTrackNInt;
    TH1D *hTrackNTpc;
    TH1D *hTrackNTot;
    TH1D *hTrackPerMms;
    TH1D *hTrackPerMap;
    TH1D *hTrackPerInt;
    TH1D *hTrackPerTpc;
    TH1D *hTrackPerTot;
    TH1D *hTrackChi2;
    TH1D *hTrackNDF;
    TH1D *hTrackQuality;
    TH1D *hTrackDCAxy;
    TH1D *hTrackDCAz;
    TH1D *hTrackEta;
    TH1D *hTrackPhi;
    TH1D *hTrackPt;
    TH1D *hDeltaDCAxy;
    TH1D *hDeltaDCAz;
    TH1D *hDeltaEta;
    TH1D *hDeltaPhi;
    TH1D *hDeltaPt;
    TH2D *hTrackPtVsNMms;
    TH2D *hTrackPtVsNMap;
    TH2D *hTrackPtVsNInt;
    TH2D *hTrackPtVsNTpc;
    TH2D *hTrackPtVsNTot;
    TH2D *hTrackPtVsPerMms;
    TH2D *hTrackPtVsPerMap;
    TH2D *hTrackPtVsPerInt;
    TH2D *hTrackPtVsPerTpc;
    TH2D *hTrackPtVsPerTot;
    TH2D *hTrackPtVsChi2;
    TH2D *hTrackPtVsNDF;
    TH2D *hTrackPtVsQuality;
    TH2D *hTrackPtVsDCAxy;
    TH2D *hTrackPtVsDCAz;
    TH2D *hDeltaDCAxyVsTrkPt;
    TH2D *hDeltaDCAzVsTrkPt;
    TH2D *hDeltaEtaVsTrkPt;
    TH2D *hDeltaPhiVsTrkPt;
    TH2D *hDeltaPtVsTrkPt;

    // embed-only truth output histograms
    TH1D *hTruthNMms;
    TH1D *hTruthNMap;
    TH1D *hTruthNInt;
    TH1D *hTruthNTpc;
    TH1D *hTruthNTot;
    TH1D *hTruthEta;
    TH1D *hTruthPhi;
    TH1D *hTruthPt;
    TH1D *hTruthEtaFrac;
    TH1D *hTruthPhiFrac;
    TH1D *hTruthPtFrac;
    TH1D *hTruthEtaDiff;
    TH1D *hTruthPhiDiff;
    TH1D *hTruthPtDiff;
    TH2D *hTruthVsTrackEta;
    TH2D *hTruthVsTrackPhi;
    TH2D *hTruthVsTrackPt;
    TH2D *hFracVsTruthEta;
    TH2D *hFracVsTruthPhi;
    TH2D *hFracVsTruthPt;
    TH2D *hDiffVsTruthEta;
    TH2D *hDiffVsTruthPhi;
    TH2D *hDiffVsTruthPt;
    TH2D *hTruthPtVsNMms;
    TH2D *hTruthPtVsNMap;
    TH2D *hTruthPtVsNInt;
    TH2D *hTruthPtVsNTpc;
    TH2D *hTruthPtVsNTot;
    TH2D *hTruthPtVsChi2;
    TH2D *hTruthPtVsNDF;
    TH2D *hTruthPtVsQuality;
    TH2D *hTruthPtVsDCAxy;
    TH2D *hTruthPtVsDCAz;
    TH2D *hDeltaDCAxyVsTruPt;
    TH2D *hDeltaDCAzVsTruPt;
    TH2D *hDeltaEtaVsTruPt;
    TH2D *hDeltaPhiVsTruPt;
    TH2D *hDeltaPtVsTruPt;

    // embed-only weird output histograms
    TH1D *hWeirdNMms;
    TH1D *hWeirdNMap;
    TH1D *hWeirdNInt;
    TH1D *hWeirdNTpc;
    TH1D *hWeirdNTot;
    TH1D *hWeirdPerMms;
    TH1D *hWeirdPerMap;
    TH1D *hWeirdPerInt;
    TH1D *hWeirdPerTpc;
    TH1D *hWeirdPerTot;
    TH1D *hWeirdChi2;
    TH1D *hWeirdNDF;
    TH1D *hWeirdQuality;
    TH1D *hWeirdDCAxy;
    TH1D *hWeirdDCAz;
    TH1D *hWeirdEta;
    TH1D *hWeirdPhi;
    TH1D *hWeirdPt;
    TH1D *hWeirdDeltaDCAxy;
    TH1D *hWeirdDeltaDCAz;
    TH1D *hWeirdDeltaEta;
    TH1D *hWeirdDeltaPhi;
    TH1D *hWeirdDeltaPt;

    // declare with-pileup track histograms
    TH1D *hTrackNMms_PU;
    TH1D *hTrackNMap_PU;
    TH1D *hTrackNInt_PU;
    TH1D *hTrackNTpc_PU;
    TH1D *hTrackNTot_PU;
    TH1D *hTrackPerMms_PU;
    TH1D *hTrackPerMap_PU;
    TH1D *hTrackPerInt_PU;
    TH1D *hTrackPerTpc_PU;
    TH1D *hTrackPerTot_PU;
    TH1D *hTrackChi2_PU;
    TH1D *hTrackNDF_PU;
    TH1D *hTrackQuality_PU;
    TH1D *hTrackDCAxy_PU;
    TH1D *hTrackDCAz_PU;
    TH1D *hTrackEta_PU;
    TH1D *hTrackPhi_PU;
    TH1D *hTrackPt_PU;
    TH1D *hDeltaDCAxy_PU;
    TH1D *hDeltaDCAz_PU;
    TH1D *hDeltaEta_PU;
    TH1D *hDeltaPhi_PU;
    TH1D *hDeltaPt_PU;
    TH2D *hTrackPtVsNMms_PU;
    TH2D *hTrackPtVsNMap_PU;
    TH2D *hTrackPtVsNInt_PU;
    TH2D *hTrackPtVsNTpc_PU;
    TH2D *hTrackPtVsNTot_PU;
    TH2D *hTrackPtVsPerMms_PU;
    TH2D *hTrackPtVsPerMap_PU;
    TH2D *hTrackPtVsPerInt_PU;
    TH2D *hTrackPtVsPerTpc_PU;
    TH2D *hTrackPtVsPerTot_PU;
    TH2D *hTrackPtVsChi2_PU;
    TH2D *hTrackPtVsNDF_PU;
    TH2D *hTrackPtVsQuality_PU;
    TH2D *hTrackPtVsDCAxy_PU;
    TH2D *hTrackPtVsDCAz_PU;
    TH2D *hDeltaDCAxyVsTrkPt_PU;
    TH2D *hDeltaDCAzVsTrkPt_PU;
    TH2D *hDeltaEtaVsTrkPt_PU;
    TH2D *hDeltaPhiVsTrkPt_PU;
    TH2D *hDeltaPtVsTrkPt_PU;

    // calculation parameters
    Double_t weirdPtFracMin;
    Double_t weirdPtFracMax;

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
    void InitFiles();
    void InitTuples();
    void InitHists();
    void SetHistStyles();
    void CreatePlots();
    void SaveHists();

};  // end STrackCutStudy definition

#endif

// end ------------------------------------------------------------------------
