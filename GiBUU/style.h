#ifndef STYLE_H
#define STYLE_H

#include <math.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>

#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

//#include "TASImage.h"
#include "TAxis.h"
#include "TColor.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGeoManager.h"
#include "TGeoGlobalMagField.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraphPolar.h"
#include "TGrid.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLinearFitter.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TPolyMarker.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TUUID.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TVirtualPad.h"
#include "THStack.h"

#ifndef EPSILON
#define EPSILON 1e-12
#endif

using namespace std;

class style
{
 public:
  static TH1D * GetMCDataNsigma(TH1D * hmc, TH1D * hdata, const TString tag);
  static void SetStackColorCB(THStack * hh, const int *cbcol);
  static void IniColorCB();
  static TGraphAsymmErrors* ScaleGraph(const TGraphAsymmErrors * gin, const  TGraphAsymmErrors * gref);
  static void GraphMinMaxXY(const TGraphAsymmErrors * gr, Double_t & xmin, Double_t & xmax, Double_t &ymin, Double_t &ymax );
  static void PadSetup(TPad *currentPad, const Double_t currentLeft=0.12, const Double_t currentTop=0.09, const Double_t currentRight=0.13, const Double_t currentBottom=0.14);
  static void ToNaturalScale(TGraph* gr);
  static void ToNaturalScale(TGraphErrors* gr);
  static void ToNaturalScale(TGraphAsymmErrors* gr);
  static void ToNaturalScale(TAxis *ax);
  static void ToNaturalScale(TH1 *hh);
  static void AxisStyle(TAxis *ax, Bool_t kcen=0);
  static void ResetStyle(TLegend *obj, Double_t mar=-999, const Double_t ff=0.8);
  static void ResetStyle(TLatex *obj, const Double_t ff=1, const bool kndc = kTRUE);
  static void ResetStyle(TPaveText *obj, const Double_t ff=1);
  static void ResetStyle(TGraph *obj);
  static void ResetStyle(TF1 * obj, Bool_t kcen=0);
  static void ResetStyle(TH1 * obj, TVirtualPad* cpad=0x0, Bool_t kcen=0);
  static TH1D * GetStackedSum(THStack *stk, const char * name=0x0, const int col=kRed, const int lsty=kSolid, const int lwid=2, const int fsty=0);
  static void GetBins(TH1D * hh, double *xb, const int noff, const int ndim);
  static TMatrixD CleanErrorMatrix(TH1D * rawhist, const TMatrixD rawcov, const Int_t noff, const Int_t ndim, TH2D * & covHist, const double unit=1);
  static TMatrixD GetDiagonal(const TMatrixD err, const bool kinvert= false);
  static void Corr2Err(const TMatrixD corrMatrix, TMatrixD & err);
  static void Err2Corr(const TMatrixD err, TMatrixD & corrMatrix);
  static double GetChi2(TH1D * rawdata, const TMatrixD rawcov, const int noff, const int ndim, TH1D * rawtheory, const double unit=1);
  static void ResetStyle(THStack * obj,  Bool_t kcen=0);
  static void SetGlobalStyle(const Int_t lStat=0, Bool_t kcolor = 1);
  static void SetColor();

  static void BinLogX(TAxis *axis);
  static void UpdateLogX(TH1 *obj);
  static void UpdateLogX(TGraphAsymmErrors *obj);
  static void UpdateLogX(TGraphErrors *obj);

  static void DividegPad(int nx, int ny, double l, double r, double t, double b);
  static void DividePad(TPad * pin, int nx, int ny, double l, double r, double t, double b, TList * lout, Int_t &npad, const Int_t lastny);

  static TGraphAsymmErrors * GetInverse(const TGraphAsymmErrors * graw);
  static TH1D * Func2Hist(TF1 * ff, const Bool_t klogx);

  static Double_t fgkTextSize;
  static Double_t fgkTitleSize;
  static Double_t fgkMarkerSize;
  static Double_t fgkLineWidth;
  static Int_t fgkTextFont;
  static Double_t fgkLabelOffset;
  static Double_t fgkXTitleOffset;
  static Double_t fgkYTitleOffset;
  static Double_t fgkTickLength;

 private:
};

#endif
