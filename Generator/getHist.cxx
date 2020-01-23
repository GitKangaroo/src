#include "stdlib.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TSpline.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#define CC0piNpID 1
#define CC1piNpID 2
#define LOWRECOIL 8
#define LOWRECOIL0piNp 9
#define NUBAR1PI 16
#define MMECCQE 32

const bool OPENLR = false;
//true;

using namespace std;

//===============================================>
class GraphSpline
{
 public:
  GraphSpline(TGraph* g)
    {
      spline=new TSpline3( (std::string("graphspline") + g->GetName()).c_str(),g);
    }

  ~GraphSpline() { delete spline; };

  double operator()(double* x, double*)
  {
    return spline->Eval(x[0]);
  }

 private:
  TSpline3* spline;
};

TF1* GENIEgetCCSpline(const int nuPDG, const int targetA);

double CHspline(double *x, double *p)
{
  TF1 * spC = GENIEgetCCSpline((int)(p[0]), 12);
  TF1 * spH = GENIEgetCCSpline((int)(p[0]), 1);
  return spC->Eval(x[0])+spH->Eval(x[0]);
}

TF1* GENIEgetCCSpline(const int nuPDG, const int targetA)
{
  if(targetA==13){
    TF1 * fCH = new TF1("CCsp13",CHspline, 0, 120, 1);
    fCH->SetParameter(0,nuPDG);
    return fCH;
  }
  else{
    const TString objname=Form("CCsp%d",targetA);
    TF1 * sp = (TF1*) gDirectory->FindObject(objname);
    if(!sp){
      printf("initiating spline for nuPDG %d targetA %d\n", nuPDG, targetA);

      const TString splineFileName("/data/t2k/xlu/software/GENIE_v3_00_06/genie_xsec/v3_00_04a/NULL/G1802a00000-k250-e1000/data/xsec_graphs.root");
      TFile splineFile(splineFileName);
      if(splineFile.IsZombie()){
        cerr << "Can't find GENIE spline file at:" << endl;
        cerr << splineFileName << endl;
        cerr << "Exiting" << endl;
        exit(1);
      }
      
      TString splineName = "nu_";
      switch(abs(nuPDG)){
      case 14:
        splineName+="mu";
        break;
      case 12:
        splineName+="e";
        break;
      default:
        cerr << "XSecLooper::getCCSpline(): nu PDG is set to " << nuPDG
             << ", which we don't have splines for. Bailing" << endl;
        exit(1);
      }
      if(nuPDG<0) splineName+="_bar";
      
      if(targetA==12){
        splineName+="_C12/tot_cc";
      }
      else if(targetA==1){
        splineName+="_H1/tot_cc";
      }
      else{
        printf("wrong targetA %d\n", targetA); exit(1);
      }
      
      TGraph* gr=(TGraph*)splineFile.Get(splineName);
      if(!gr){
        cerr << "Can't find spline " << splineName << " in file " << splineFileName << ". Bailing." << endl;
        exit(1);
      }
      
      GraphSpline * gsp=new GraphSpline(gr);
      sp = new TF1(objname, gsp, 0, 120, 0, "GraphSpline" );
    }
    return sp;
  }
}

double GENIEgetFluxIntegral(const TH1F * h_rate,  const int nuPDG, const int targetA)
{
  TF1 * spline = GENIEgetCCSpline(nuPDG, targetA);
  
  double fluxInt = 0;
  for(int ii=1; ii<=h_rate->GetNbinsX(); ii++){

    //only do it within Enu cut!
    if(h_rate->GetBinContent(ii)<1){
      continue;
    }

    const double binLowEdge = h_rate->GetBinLowEdge(ii);
    const double binWidth   = h_rate->GetBinWidth(ii);

    //dimension: cross section
    const double xsec = spline->Integral(binLowEdge, binLowEdge+binWidth)/binWidth;

    if( xsec < 1E-6 ){
      std::cout << " Warning, xsec is 0 for bin from " << binLowEdge << " - " <<  binLowEdge+binWidth << endl;
    }
    else{
      //has bin width in it
      fluxInt += h_rate->GetBinContent(ii)/xsec;
      //printf("binLowEdge %f binWidth %f xsec %f rate %f flux %f\n", binLowEdge, binWidth, xsec, h_rate->GetBinContent(ii), fluxInt);
    }
  }

  return fluxInt;
}

/*//too slow
void drawInSlice(const int ncutbin, double cutbins[], const TString cutvar, const int ndrawbin, const double drawbins[], const TString drawvar, TTree *tree,  const TString cut0, const TString modename, TList *lout)
{
  int totnd = 0;
  const double *binlow=cutbins;
  const double *binup=&(cutbins[1]);

  for(int ibin=0; ibin<ncutbin; ibin++){
    const double binwidth = binup[ibin]-binlow[ibin];
    const TString tmpcut = Form("perweight * %f * (%s %s) ", 1/binwidth, cut0.Data(), Form("&& (%f<%s && %s<%f)", binlow[ibin], cutvar.Data(), cutvar.Data(), binup[ibin]));

    const TString hname=Form("%s%d%s", drawvar.Data(), ibin, modename.Data());
    TH1D * hhii = new TH1D(hname, tmpcut, ndrawbin, drawbins); lout->Add(hhii);
    const int ndraw = tree->Project(hname, drawvar, tmpcut);
    totnd+=ndraw;
    cout<<"Drawing "<<hname<<" "<<drawvar<<" tmpcut "<<tmpcut<<" ndraw "<<ndraw<<" totnd "<<totnd<<endl;
  }
}
*/
void getSliceYDrawX(TH2D * h2d, TList * lout)
{
  for(int ii=1; ii<=h2d->GetNbinsY(); ii++){
    const int histid=ii-1;
    TString tmpname(h2d->GetName());
    //remove xxxVS
    tmpname=tmpname(tmpname.First("S")+1,tmpname.Length());
    tmpname.Insert(tmpname.First("_"),Form("%d",histid));
    TH1D * hpj= h2d->ProjectionX(tmpname, ii, ii);
    hpj->SetTitle(Form("%s %f %f", h2d->GetName(), h2d->GetYaxis()->GetBinLowEdge(ii), h2d->GetYaxis()->GetBinUpEdge(ii)));
    lout->Add(hpj);
  }
}
//===============================================<

void getHist(const TString fn, const TString tag, const int anaid)
{
  cout<<"please check "<<fn<<" "<<tag<<" "<<anaid<<endl;

  //===========================================================>
  TString nuExp("NONE");
  if(fn.Contains("MINERvA"))
    nuExp = "MINERvA";
  else if(fn.Contains("T2K"))
    nuExp = "T2K";
  else if(fn.Contains("DUNE"))
    nuExp = "DUNE";

  TString tmpfilename(fn); tmpfilename.ToUpper();
  int nuPDG = 14;
  if(tmpfilename.Contains("NUBAR")){
    nuPDG = -14;
  }
  cout<<"nuPDG "<<nuPDG<<endl;

  int tarA=12;
  if(tmpfilename.Contains("HYDROGEN")){
    tarA = 1;
  }
  else if(tmpfilename.Contains("CH")){
    tarA = 13;
  }
  else if(tmpfilename.Contains("ARGON")){
    tarA = 40;
  }
  printf("\n\n****************************************** Using tarA %d **********************************************\n\n\n", tarA);

  /*
  if(anaid == CC1piNpID && fn.Contains("GENIE")){
    //tarA = 12;//no hydrogen for CC1piNpID
    tarA = 13;
    printf("\n\n****************************************** Using tarA %d **********************************************\n\n\n", tarA);
  }
  printf("tarA %d\n", tarA);
  */

  TFile * fin = new TFile(fn);

  int GiBUUnormFactor=-999;
  //update GiBUUnormFactor from header tree
  TTree * theader = (TTree*) fin->Get("header");
  if(theader){
    int tmpanr;
    theader->SetBranchAddress("nrun", &tmpanr);
    theader->GetEntry(0);
    if(tmpanr>0){
      GiBUUnormFactor=tmpanr;
    }
  }
  cout<<"GiBUUnormFactor "<<GiBUUnormFactor<<endl;

  double GENIEnormFactor = -999;
  TH1F * GENIEhCCrate=(TH1F*) fin->Get("hCCrate");
  if(GENIEhCCrate){
    //make it flux*Nnucleon
    GENIEnormFactor = GENIEgetFluxIntegral(GENIEhCCrate, nuPDG, tarA);
  }
  cout<<"GENIEnormFactor "<<GENIEnormFactor<<endl;

  double histnormFactor=1;
  TString gen;
  if(fn.Contains("GENIE")){
    if(tarA!=13){
      printf("GENIE only has tarA=13! %d\n", tarA); exit(1);
    }

    gen="GENIE";
    histnormFactor=GENIEnormFactor*tarA;//spline xsec is for each nucleus, now norm with tarA gives xsec per nucleon
  }
  else if(fn.Contains("GiBUU")){
    if(tarA!=1 && tarA!=12 && tarA!=40){
      printf("GiBUU can't have other tarA than 1 or 12!! %d\n", tarA);exit(1);
    }

    gen="GiBUU";
    histnormFactor=GiBUUnormFactor/tarA;//perweight is for per nucleon, norm with 1/tarA gives xsec per nucleus, then get xsec per nucleon in doGetStack.sh
  }
  else{
    printf("wrong generator name in file name! %s\n", fn.Data());
    exit(1);
  }
  cout<<"histnormFactor "<<histnormFactor<<endl;

  //=======================================================================<

  TTree * tree = (TTree*) fin->Get("tree");

  TList *lout=new TList;

  const TString modes[]={"_all","_qe","_res","_dis","_2p2h", "_other"};
  /*
https://gibuu.hepforge.org/Documentation/code/analysis/neutrinoAnalysis_f90.html#robo1063
 contains info on the very first neutrino-interaction with the nucleus:

    1: nucleon (QE)
    2-31: non-strange baryon resonance (as in IdTable)
    32: pi neutron-background (e.g. nu + n -> mu + pi+ + n)
    33: pi proton-background (e.g. nu + n -> mu + pi0 + p)
    34: DIS
    35: 2p2h QE
    36: 2p2h Delta
    37: two pion background

   */
  /*
enum mode{
    kALL=0,
    kQE=1,
    kRES=2,
    kDIS=3,
    k2P2H=4,
    kOTHER=5
  };
   */
  //const TString cuts[]={"1","(prod==1)", "(prod>=2 && prod<=33)", "(prod==34)", "(prod==35)", "(prod>=36)"};
  const TString cuts[]={"1","(evtMode==1)", "(evtMode==2)", "(evtMode==3)", "(evtMode==4)", "(evtMode>4)"};
  //const int pdsmin[]={0,   1, 2,  34, 35, 36};
  //const int pdsmax[]={1000,1, 33, 34, 35, 1000};
  const int pdsmin[]={0,   1, 2, 3, 4, 5};
  const int pdsmax[]={1000,1, 2, 3, 4, 1000};
  const Int_t nmode = sizeof(modes)/sizeof(TString);

  TString basecut="&& (perweight < 4000)";
  int minnp = 0, maxnp = 1000;
  if(anaid == CC0piNpID){
    basecut += "&& ( npar>=101 && npar<110 )";
    minnp = 101;
    maxnp = 109;
  }
  else if(anaid == CC1piNpID){
    basecut += "&& ( npar>=111 && npar<120 ) && (targetZ!=1) ";//no hydrogen
    printf("\n\n************************************** Require NO Hydrogen in the basecut! %d \"%s\"\n\n\n", anaid, basecut.Data());
    minnp = 111;
    maxnp = 119;
  }
  else if(anaid == LOWRECOIL){
    basecut += "";
  }
  else if(anaid == LOWRECOIL0piNp){
    basecut += "&& ( npar>=101 && npar<110 )";
  }
  else if(anaid == NUBAR1PI){
    basecut += "";
  }
  else if(anaid == MMECCQE){
    basecut += "";
  }
  else{
    printf("anaid not known! %d\n", anaid); exit(1);
  }

  for(int imode=0; imode<nmode; imode++){

    //================================================================ COMMON FILLING ========================================================================
    //******************************************************** Definition *****************************************
    int npar, evtMode, targetZ; //prod;      
    double perweight;
    tree->SetBranchAddress("npar", &npar);
    tree->SetBranchAddress("evtMode", &evtMode);
    tree->SetBranchAddress("targetZ", &targetZ);
    tree->SetBranchAddress("perweight", &perweight);
    //tree->SetBranchAddress("prod", &prod);

    //******************************************************** Definition *****************************************
    double Q2, xBj, xrest, Wtrue, Wrest;
    //, energyCCH, xCCH;
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("xBj", &xBj);
    tree->SetBranchAddress("xrest", &xrest);
    tree->SetBranchAddress("Wtrue", &Wtrue);
    tree->SetBranchAddress("Wrest", &Wrest);
    //tree->SetBranchAddress("energyCCH", &energyCCH);
    //tree->SetBranchAddress("xCCH", &xCCH);
    TH1D * hQ2  = new TH1D("Q2"+modes[imode],"", 30, 0, 2); lout->Add(hQ2);
    TH1D * hxBj = new TH1D("xBj"+modes[imode],"", 30, 0, 2); lout->Add(hxBj);
    TH1D * hxrest = new TH1D("xrest"+modes[imode],"", 30, 0, 2); lout->Add(hxrest);
    TH1D * hWtrue = new TH1D("Wtrue"+modes[imode],"", 60, 0, 3); lout->Add(hWtrue);
    TH1D * hWrest = new TH1D("Wrest"+modes[imode],"", 60, 0, 3); lout->Add(hWrest);
    //TH1D *henergyCCH  = new TH1D("energyCCH"+modes[imode],"", 60, 0, 20); lout->Add(henergyCCH);
    //TH1D * hxCCH = new TH1D("xCCH"+modes[imode],"", 30, 0, 2); lout->Add(hxCCH);

    //******************************************************** Definition *****************************************
    double muonmomentum, muontheta, enu;
    TH1D * hmuonmomentum = 0x0, * hmuontheta = 0x0, * henu  = 0x0;
    //at NUBAR1PI all these will have specific binning, therefore not needed here
    if(anaid!=NUBAR1PI){
      tree->SetBranchAddress("muonmomentum", &muonmomentum);
      tree->SetBranchAddress("muontheta", &muontheta);
      tree->SetBranchAddress("enu", &enu);

      const double Abin[]={0, 0.25, 0.5, 0.75, 1, 1.25, 1.500000, 2.250000, 2.500000, 2.750000, 3.000000, 3.250000, 3.500000, 3.750000, 4.000000, 4.250000, 4.500000, 4.750000, 5.000000, 5.250000, 5.500000, 5.750000, 6.000000, 6.250000, 6.500000, 6.750000, 7.000000, 7.250000, 7.500000, 7.750000, 8.000000, 8.250000, 8.500000, 8.750000, 9.000000, 9.250000, 9.500000, 9.750000, 10.000000}; 
      hmuonmomentum = new TH1D("muonmomentum"+modes[imode],"", sizeof(Abin)/sizeof(double)-1, Abin); lout->Add(hmuonmomentum); 
    
      const double Bbin[]={0.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000, 10.000000, 11.000000, 12.000000, 13.000000, 14.000000, 15.000000, 16.000000, 17.000000, 18.000000, 19.000000, 20.000000, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130};
      hmuontheta = new TH1D("muontheta"+modes[imode],"", sizeof(Bbin)/sizeof(double)-1, Bbin); lout->Add(hmuontheta);

      henu  = new TH1D("enu"+modes[imode],"", 60, 0, 20); lout->Add(henu);
    }

    //******************************************************** Definition *****************************************
    double RESmass, adlerPhi, lrsign, w2, KNsrc;
    //pseudoPhi, pseudosign, wpseudo2, 
    TH1D * hRESmass = 0x0, * hadlerPhi = 0x0, * hoobphi = 0x0, * hlrsign = 0x0, * hw2 = 0x0, * hKNsrc = 0x0;
    //* hpseudoPhi = 0x0, * hpseudosign = 0x0, * hwpseudo2 = 0x0, 
    if(OPENLR){      
      tree->SetBranchAddress("RESmass", &RESmass);
      tree->SetBranchAddress("adlerPhi", &adlerPhi);
      tree->SetBranchAddress("lrsign", &lrsign);
      tree->SetBranchAddress("w2", &w2);
      tree->SetBranchAddress("KNsrc", &KNsrc);

      hRESmass = new TH1D("RESmass"+modes[imode],"", 60, 0, 3); lout->Add(hRESmass);
      hadlerPhi = new TH1D("adlerPhi"+modes[imode],"", 60, 0, 360); lout->Add(hadlerPhi);
      hoobphi = new TH1D("oobphi"+modes[imode],"", 60, 0, 360); lout->Add(hoobphi);
      hlrsign = new TH1D("lrsign"+modes[imode],"", 3, -1.5, 1.5); lout->Add(hlrsign);
      hw2 = new TH1D("w2"+modes[imode],"", 3, -1.5, 1.5); lout->Add(hw2); 
      hKNsrc = new TH1D("KNsrc"+modes[imode],"", 3, -1.5, 1.5); lout->Add(hKNsrc);

      /*
      tree->SetBranchAddress("pseudoPhi", &pseudoPhi);
      tree->SetBranchAddress("pseudosign", &pseudosign);
      tree->SetBranchAddress("wpseudo2", &wpseudo2);
      hpseudoPhi = new TH1D("pseudoPhi"+modes[imode],"", 60, 0, 360); lout->Add(hpseudoPhi);
      hpseudosign = new TH1D("pseudosign"+modes[imode],"", 3, -1.5, 1.5); lout->Add(hpseudosign);
      hwpseudo2 = new TH1D("wpseudo2"+modes[imode],"", 3, -1.5, 1.5); lout->Add(hwpseudo2); 
      */
    }

    //********************************************************************************************************

    int ient=0;
    int ndraw = 0;
    while(tree->GetEntry(ient++)){
      if(perweight>4000){
        printf("\n\n\nAlert!! Common Filling Super weight!! %d %f skiping...\n\n\n", ient, perweight); 
        continue;
      }

      if(evtMode<pdsmin[imode] || evtMode>pdsmax[imode]){
        continue;
      }

      //-------- mode-specific cuts ->
      bool passtopology = true;
      if(anaid==CC0piNpID || anaid == CC1piNpID){
        if(npar<minnp || npar>maxnp){
          passtopology = false;
        }
      }

      bool passtarget = true;
      if(anaid == CC1piNpID){
        if(targetZ == 1){//no hydrogen for GFS1
          passtarget = false;
        }
      }

      bool passwr = true;
      if(anaid==NUBAR1PI){
        if(Wrest>=1.8){
          passwr = false;
        }
      }
      //<-------------------------------

      if(passtopology && passtarget){
        //no cut on W for W
        hWtrue->Fill(Wtrue,perweight);
        hWrest->Fill(Wrest,perweight);

        if(passwr){
          hQ2->Fill(Q2,perweight);
          hxBj->Fill(xBj,perweight);
          hxrest->Fill(xrest,perweight);
          //henergyCCH->Fill(energyCCH,perweight);
          //hxCCH->Fill(xCCH,perweight);

          if(anaid!=NUBAR1PI){
            hmuonmomentum->Fill(muonmomentum,perweight);
            hmuontheta->Fill(muontheta,perweight);
            henu->Fill(enu,perweight);
          }
          
          if(OPENLR){
            hRESmass->Fill(RESmass,perweight);
            hadlerPhi->Fill(adlerPhi,w2*perweight);
            hoobphi->Fill(adlerPhi, perweight);
            hlrsign->Fill(lrsign,perweight);
            hKNsrc->Fill(KNsrc,perweight);
            //hpseudoPhi->Fill(pseudoPhi, wpseudo2*perweight);
            //hpseudosign->Fill(pseudosign,perweight);
          }
          
          ndraw++;
        }
      }
    }
      cout << "ient: " << ient << "ndraw: " <<  ndraw << endl;
    cout<<"\nCommon filling: imode "<<imode<<", "<<pdsmin[imode]<<" evtMode "<<pdsmax[imode]<<", ndraw "<<ndraw<<endl<<endl;
    //====================================================================================================================================================

    if(anaid==CC0piNpID || anaid == CC1piNpID){
      //const TString tmpcut = Form("perweight * (%s %s) ", cuts[imode].Data(), basecut.Data());
      //cout<<"imode "<<imode<<" modes "<<modes[imode]<<" cuts "<<tmpcut<<endl;

      double protonmomentum, protontheta, pionmomentum, pionEk,  piontheta, baryonmomentum, baryontheta, baryonmass, dalphat, dphit, dpt, protonTT, pionTT, dpTT, dpTy, neutronmomentum;
      tree->SetBranchAddress("protonmomentum", &protonmomentum);
      tree->SetBranchAddress("protontheta", &protontheta);
      tree->SetBranchAddress("pionmomentum", &pionmomentum);
      tree->SetBranchAddress("pionEk", &pionEk);
      tree->SetBranchAddress("piontheta", &piontheta);
      tree->SetBranchAddress("baryonmomentum", &baryonmomentum);
      tree->SetBranchAddress("baryontheta", &baryontheta);
      tree->SetBranchAddress("baryonmass", &baryonmass);
      tree->SetBranchAddress("dalphat", &dalphat);
      tree->SetBranchAddress("dphit", &dphit);
      tree->SetBranchAddress("dpt", &dpt);
      tree->SetBranchAddress("protonTT", &protonTT);
      tree->SetBranchAddress("pionTT", &pionTT);
      tree->SetBranchAddress("dpTT", &dpTT);
      tree->SetBranchAddress("dpTy", &dpTy);
      tree->SetBranchAddress("neutronmomentum", &neutronmomentum);

      const double Cbin[]={0.3, 0.35, 0.4, 0.450000, 0.600000, 0.625000, 0.650000, 0.675000, 0.700000, 0.725000, 0.750000, 0.775000, 0.800000, 0.825000, 0.850000, 0.875000, 0.900000, 0.925000, 0.950000, 0.975000, 1.000000, 1.025000, 1.050000, 1.075000, 1.100000, 1.125000, 1.150000, 1.175000, 1.200000}; 
      TH1D * hprotonmomentum = new TH1D("protonmomentum"+modes[imode],"", sizeof(Cbin)/sizeof(double)-1, Cbin); lout->Add(hprotonmomentum);
      
      const double Dbin[]={0.000000, 7.500000, 10.000000, 12.500000, 15.000000, 17.500000, 20.000000, 22.500000, 25.000000, 27.500000, 30.000000, 32.500000, 35.000000, 37.500000, 40.000000, 42.500000, 45.000000, 47.500000, 50.000000, 52.500000, 55.000000, 57.500000, 60.000000, 62.500000, 65.000000, 67.500000, 70.000000}; 
      TH1D * hprotontheta = new TH1D("protontheta"+modes[imode],"", sizeof(Dbin)/sizeof(double)-1, Dbin); lout->Add(hprotontheta);

      TH1D * hpionmomentum = new TH1D("pionmomentum"+modes[imode],"", 40, 0, 1); lout->Add(hpionmomentum);
      TH1D * hpionEk = new TH1D("pionEk"+modes[imode],"", 40, 0, 0.5); lout->Add(hpionEk);
      TH1D * hpiontheta = new TH1D("piontheta"+modes[imode],"", 35, 0, 70); lout->Add(hpiontheta);
      TH1D * hbaryonmomentum = new TH1D("baryonmomentum"+modes[imode],"", 30, 0.2, 6); lout->Add(hbaryonmomentum);
      TH1D * hbaryontheta = new TH1D("baryontheta"+modes[imode],"", 35, 0, 70); lout->Add(hbaryontheta);
      TH1D * hbaryonmass = new TH1D("baryonmass"+modes[imode],"", 60, 0, 3); lout->Add(hbaryonmass);

      TH1D * hdalphat = 0x0, * hdphit = 0x0, * hdpt = 0x0, * hneutronmomentum = 0x0, * hprotonTT = 0x0, * hpionTT = 0x0, * hdpTT = 0x0, * hdpTy = 0x0, * hlrdpt = 0x0, * hlrdpTT = 0x0, * hKNdpt = 0x0, * hsumdpt = 0x0;
      //* hpseudodpt = 0x0, 

      if(nuExp=="MINERvA" || nuExp=="NONE" || nuExp =="DUNE"){
        const double Ebin[]= {0.000000, 20.000000, 40.000000, 60.000000, 80.000000, 100.000000, 120.000000, 130.000000, 140.000000, 150.000000, 160.000000, 170.000000, 180.000000};
        hdalphat = new TH1D("dalphat"+modes[imode],"", sizeof(Ebin)/sizeof(double)-1, Ebin); 

        const double Fbin[]={0.000000, 2.500000, 5.000000, 7.500000, 10.000000, 12.500000, 15.000000, 17.500000, 20.000000, 22.500000, 25.000000, 27.500000, 30.000000, 35.000000, 40.000000, 45.000000, 50.000000, 55.000000, 60.000000, 70.000000, 85.000000, 105.000000, 130.000000, 180.000000};
        hdphit = new TH1D("dphit"+modes[imode],"", sizeof(Fbin)/sizeof(double)-1, Fbin);

        const double Gbin[]={0.000000, 0.025000, 0.050000, 0.075000, 0.100000, 0.125000, 0.150000, 0.175000, 0.200000, 0.225000, 0.250000, 0.275000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.800000, 1.000000, 1.200000, 2.000000};
        hdpt = new TH1D("dpt"+modes[imode],"", sizeof(Gbin)/sizeof(double)-1, Gbin);

        const double Hbin[]={0.000000, 0.025000, 0.050000, 0.075000, 0.100000, 0.125000, 0.150000, 0.175000, 0.200000, 0.225000, 0.250000, 0.275000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.800000, 1.000000, 1.200000, 2.000000};
        hneutronmomentum = new TH1D("neutronmomentum"+modes[imode],"", sizeof(Hbin)/sizeof(double)-1, Hbin); lout->Add(hneutronmomentum);

        const double Jbin[]={-2, -0.800000, -0.700000, -0.650000, -0.600000, -0.550000, -0.500000, -0.450000, -0.400000, -0.350000, -0.300000, -0.250000, -0.200000, -0.150000, -0.100000, -0.050000, 0.000000, 0.050000, 0.100000, 0.150000, 0.200000, 0.250000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.800000, 2};
        hdpTT = new TH1D("dpTT"+modes[imode],"", sizeof(Jbin)/sizeof(double)-1, Jbin);

        hprotonTT = new TH1D("protonTT"+modes[imode],"", sizeof(Jbin)/sizeof(double)-1, Jbin);
        hpionTT = new TH1D("pionTT"+modes[imode],"", sizeof(Jbin)/sizeof(double)-1, Jbin);

        const double Kbin[]={-2.000000, -1.200000, -1.000000, -0.800000, -0.700000, -0.650000, -0.600000, -0.550000, -0.500000, -0.450000, -0.400000, -0.350000, -0.300000, -0.250000, -0.200000, -0.150000, -0.100000, -0.050000, 0.000000, 0.050000, 0.100000, 0.150000, 0.200000, 0.250000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 2,};
        hdpTy = new TH1D("dpTy"+modes[imode],"", sizeof(Kbin)/sizeof(double)-1, Kbin);
      }
      else if(nuExp=="T2K"){
        if(anaid==CC0piNpID){//defined by T2K analysis
          const double Ebin[] = {0.00000, 26.929016, 58.441695, 88.235500, 113.445643, 134.072124, 151.260858, 165.584803, 180.000};
          hdalphat = new TH1D("dalphat"+modes[imode],"", sizeof(Ebin)/sizeof(double)-1, Ebin);
        }
        else{
          const double Ebin[]= {0.000000, 20.000000, 40.000000, 60.000000, 80.000000, 100.000000, 120.000000, 130.000000, 140.000000, 150.000000, 160.000000, 170.000000, 180.000000};
          hdalphat = new TH1D("dalphat"+modes[imode],"", sizeof(Ebin)/sizeof(double)-1, Ebin);
        }

        const double Fbin[]={0.00000, 3.838817, 8.021409, 12.891550, 19.480565, 29.793805, 48.701413, 85.943669, 180.00};
        hdphit = new TH1D("dphit"+modes[imode],"", sizeof(Fbin)/sizeof(double)-1, Fbin);

        const double Gbin[]={0.000000, 0.080000, 0.120000, 0.155000, 0.200000, 0.260000, 0.360000, 0.510000, 1.10000};
        hdpt = new TH1D("dpt"+modes[imode],"", sizeof(Gbin)/sizeof(double)-1, Gbin);

        const double Hbin[]={0.000000, 0.030000, 0.060000, 0.090000, 0.120000, 0.150000, 0.180000, 0.210000, 0.240000, 0.270000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.800000, 1.000000, 1.200000, 1.5, 2.000000};
        hneutronmomentum = new TH1D("neutronmomentum"+modes[imode],"", sizeof(Hbin)/sizeof(double)-1, Hbin); lout->Add(hneutronmomentum);

        const double Jbin[]={-2, -0.800000, -0.700000, -0.650000, -0.600000, -0.550000, -0.500000, -0.450000, -0.400000, -0.350000, -0.300000, -0.250000, -0.200000, -0.150000, -0.100000, -0.050000, 0.000000, 0.050000, 0.100000, 0.150000, 0.200000, 0.250000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 0.800000, 2};
        hdpTT = new TH1D("dpTT"+modes[imode],"", sizeof(Jbin)/sizeof(double)-1, Jbin);
        hprotonTT = new TH1D("protonTT"+modes[imode],"", sizeof(Jbin)/sizeof(double)-1, Jbin);
        hpionTT = new TH1D("pionTT"+modes[imode],"", sizeof(Jbin)/sizeof(double)-1, Jbin);

        const double Kbin[]={-2.000000, -1.200000, -1.000000, -0.800000, -0.700000, -0.650000, -0.600000, -0.550000, -0.500000, -0.450000, -0.400000, -0.350000, -0.300000, -0.250000, -0.200000, -0.150000, -0.100000, -0.050000, 0.000000, 0.050000, 0.100000, 0.150000, 0.200000, 0.250000, 0.300000, 0.350000, 0.400000, 0.450000, 0.500000, 0.550000, 0.600000, 0.650000, 0.700000, 2,};
        hdpTy = new TH1D("dpTy"+modes[imode],"", sizeof(Kbin)/sizeof(double)-1, Kbin);
      }

      if(OPENLR){
        hlrdpt = new TH1D("lrdpt"+modes[imode],"", 20, 0, 2);
        hlrdpTT = new TH1D("lrdpTT"+modes[imode],"", 40, -2, 2);
        hKNdpt = new TH1D("KNdpt"+modes[imode],"", 20, 0, 2);
        hsumdpt = new TH1D("sumdpt"+modes[imode],"", 20, 0, 2);
        //hpseudodpt = new TH1D("pseudodpt"+modes[imode],"", 20, 0, 2);
      }

      TH1D * hstmp[]={hdalphat, hdphit, hdpt, hprotonTT, hpionTT, hdpTT, hdpTy, hlrdpt, hlrdpTT, hKNdpt, hsumdpt};
      //hpseudodpt, 
      for(int itmp=0; itmp<(int) (sizeof(hstmp)/sizeof(TH1D*)); itmp++){
        if(hstmp[itmp]) lout->Add(hstmp[itmp]);
      }
 
      int ient=0;
      int ndraw = 0;
      while(tree->GetEntry(ient++)){
        if(perweight>4000){
          printf("\n\n\nAlert!! GFS Super weight!! %d %f skiping...\n\n\n", ient, perweight);
          continue;
        }
          
        if(npar<minnp || npar>maxnp){
          continue;
        }

        if(evtMode<pdsmin[imode] || evtMode>pdsmax[imode]){
          continue;
        }

        if(anaid == CC1piNpID){
          if(targetZ == 1){
            continue;
          }
        }

        hprotonmomentum->Fill(protonmomentum,perweight);
        hprotontheta->Fill(protontheta,perweight);
        hpionmomentum->Fill(pionmomentum,perweight);
        hpionEk->Fill(pionEk,perweight);
        hpiontheta->Fill(piontheta,perweight);
        hbaryonmomentum->Fill(baryonmomentum,perweight);
        hbaryontheta->Fill(baryontheta,perweight);
        hbaryonmass->Fill(baryonmass,perweight);
        hdalphat->Fill(dalphat,perweight);
        hdphit->Fill(dphit,perweight);
        hdpt->Fill(dpt,perweight);
        hprotonTT->Fill(protonTT,perweight);
        hpionTT->Fill(pionTT,perweight);
        hdpTT->Fill(dpTT,perweight);
        hdpTy->Fill(dpTy,perweight);
        hneutronmomentum->Fill(neutronmomentum,perweight);

        if(OPENLR){
          if( fabs(lrsign)<1E-6 ){
            printf("\n\n\nstrange lrsign! %d %f\n\n\n", ient, lrsign);exit(1);
          }

          hlrdpt->Fill(dpt, lrsign * perweight);
          hlrdpTT->Fill(dpTT, w2 * perweight);
          hKNdpt->Fill(dpt, KNsrc * perweight);
          hsumdpt->Fill(dpt, perweight);//sampe as dpt but with different binning
          //hpseudodpt->Fill(dpt, pseudosign * perweight);
        }

        ndraw++;
      }
      cout<<"anaid "<<anaid<<"imode "<<imode<<", "<<minnp<<" npar "<<maxnp<<", "<<pdsmin[imode]<<" evtMode "<<pdsmax[imode]<<", ndraw "<<ndraw<<endl;

        /*
      const int ndraw = tree->Project("muonmomentum"+modes[imode],"muonmomentum",tmpcut);
      cout<<"ndraw "<<ndraw<<endl;
      
      tree->Project("muontheta"+modes[imode],"muontheta",tmpcut);
      tree->Project("protonmomentum"+modes[imode],"protonmomentum",tmpcut);
      tree->Project("protontheta"+modes[imode],"protontheta",tmpcut);
      tree->Project("pionmomentum"+modes[imode],"pionmomentum",tmpcut);
      tree->Project("piontheta"+modes[imode],"piontheta",tmpcut);
      tree->Project("baryonmomentum"+modes[imode],"baryonmomentum",tmpcut);
      tree->Project("baryontheta"+modes[imode],"baryontheta",tmpcut);
      tree->Project("baryonmass"+modes[imode],"baryonmass",tmpcut);
      tree->Project("dalphat"+modes[imode],"dalphat",tmpcut);
      tree->Project("dphit"+modes[imode],"dphit",tmpcut);
      tree->Project("dpt"+modes[imode],"dpt",tmpcut);
      tree->Project("neutronmomentum"+modes[imode],"neutronmomentum",tmpcut);
      */
    }
    else{
      vector<TString> hhs;
      if(anaid==LOWRECOIL || anaid==LOWRECOIL0piNp){
        double bq3[]={0, 0.2 ,  0.3 ,  0.4 ,  0.5 ,  0.6 ,  0.8 };
        TH1D * hq3=new TH1D("q3"+modes[imode],"", sizeof(bq3)/sizeof(double)-1, bq3); lout->Add(hq3); hhs.push_back("q3");
        
        const double bEav[]={0, 0.0166667, 0.0333333, 0.05, 0.0666667, 0.0833333, 0.1, 0.116667, 0.133333, 0.15, 0.166667, 0.183333, 0.2, 0.216667, 0.233333, 0.25, 0.266667, 0.283333, 0.3, 0.316667, 0.333333, 0.35, 0.366667, 0.383333, 0.4, 0.416667, 0.433333, 0.45, 0.466667, 0.483333, 0.5};
        TH2D * hq3VSEav = new TH2D("q3VSEav"+modes[imode],"", sizeof(bEav)/sizeof(double)-1, bEav, sizeof(bq3)/sizeof(double)-1, bq3); lout->Add(hq3VSEav); hhs.push_back("q3VSEav");
        
        if(anaid==LOWRECOIL0piNp){
          const double bneutronmomentum[]={0, 0.0333333, 0.0666667, 0.1, 0.133333, 0.166667, 0.2, 0.233333, 0.266667, 0.3, 0.333333, 0.366667, 0.4, 0.433333, 0.466667, 0.5, 0.533333, 0.566667, 0.6, 0.633333, 0.666667, 0.7, 0.733333, 0.766667, 0.8, 0.833333, 0.866667, 0.9, 0.933333, 0.966667, 1, 1.03333, 1.06667, 1.1, 1.13333, 1.16667, 1.2, 1.23333, 1.26667, 1.3, 1.33333, 1.36667, 1.4, 1.43333, 1.46667, 1.5, 1.53333, 1.56667, 1.6, 1.63333, 1.66667, 1.7, 1.73333, 1.76667, 1.8, 1.83333, 1.86667, 1.9, 1.93333, 1.96667, 2};
          TH2D * hq3VSneutronmomentum = new TH2D("q3VSneutronmomentum"+modes[imode],"", sizeof(bneutronmomentum)/sizeof(double)-1, bneutronmomentum, sizeof(bq3)/sizeof(double)-1, bq3);lout->Add(hq3VSneutronmomentum); hhs.push_back("q3VSneutronmomentum");
        }
      }
      else if(anaid==MMECCQE){
        const double q2qebins[]={0 , 0.00625 , 0.0125 , 0.025 , 0.0375 , 0.05 , 0.1 , 0.15 , 0.2 , 0.3 , 0.4 , 0.6 , 0.8 , 1.0 , 1.2 , 2.0 , 4.0 , 6.0 , 8.0 , 10.0};
        TH1D *hq2qe = new TH1D("q2qe"+modes[imode],"", sizeof(q2qebins)/sizeof(double)-1, q2qebins); lout->Add(hq2qe); hhs.push_back("q2qe");
        
        double muonptbins[]={0 , 0.07 , 0.15 , 0.25 , 0.33 , 0.4 , 0.47 , 0.55 ,0.7 , 0.85 , 1.00 , 1.25 , 1.5 , 2.5};
        TH1D *hmuonpt = new TH1D("muonpt"+modes[imode],"", sizeof(muonptbins)/sizeof(double)-1, muonptbins); lout->Add(hmuonpt); hhs.push_back("muonpt");
        
        double mupzbins[]={1.5 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0 , 8.0 , 9.0 , 10.0 , 15.0 , 20.0};
        TH1D *hmupz = new TH1D("mupz"+modes[imode],"", sizeof(mupzbins)/sizeof(double)-1, mupzbins); lout->Add(hmupz); hhs.push_back("mupz");
        
        TH2D * hmuonptVSmupz = new TH2D("muonptVSmupz"+modes[imode],"", sizeof(mupzbins)/sizeof(double)-1, mupzbins, sizeof(muonptbins)/sizeof(double)-1, muonptbins); lout->Add(hmuonptVSmupz); hhs.push_back("muonptVSmupz");
        TH2D * hmupzVSmuonpt = new TH2D("mupzVSmuonpt"+modes[imode],"", sizeof(muonptbins)/sizeof(double)-1, muonptbins, sizeof(mupzbins)/sizeof(double)-1, mupzbins); lout->Add(hmupzVSmuonpt); hhs.push_back("mupzVSmuonpt");

        /*
          pt00.150.250.40.712.5
          pz1.53.5820
          recoil040801201602002402803203604006008001000
         */
        const double muonpt3dbins[]={0, 0.15, 0.25, 0.4, 0.7, 1, 2.5};
        const double mupz3dbins[]={1.5, 3.5, 8, 20};
        const double Erecoil3dbins[]={0, 40, 80, 120, 160, 200, 240, 280, 320, 360, 400, 600, 800, 1000};
        TH3D * hErecoilTIMES1E3VSmupzVSmuonpt = new TH3D("ErecoilTIMES1E3VSmupzVSmuonpt"+modes[imode], "", sizeof(muonpt3dbins)/sizeof(double)-1, muonpt3dbins, sizeof(mupz3dbins)/sizeof(double)-1, mupz3dbins, sizeof(Erecoil3dbins)/sizeof(double)-1, Erecoil3dbins); lout->Add(hErecoilTIMES1E3VSmupzVSmuonpt);; hhs.push_back("ErecoilTIMES1E3VSmupzVSmuonpt");
      }
      else if(anaid==NUBAR1PI){
        const double enubins[]={0, 1.5, 2.0,3.0,3.50,4.0,5.0,6.0,8.0,10.0, 12.0};
        henu = new TH1D("enu"+modes[imode],"", 10, enubins); lout->Add(henu); hhs.push_back("enu");
        
        hmuonmomentum = new TH1D("muonmomentum"+modes[imode],"", 30, 0, 10); lout->Add(hmuonmomentum); hhs.push_back("muonmomentum");
        hmuontheta = new TH1D("muontheta"+modes[imode],"", 25, 0, 25); lout->Add(hmuontheta); hhs.push_back("muontheta");

        TH1D * hpionmomentum = new TH1D("pionmomentum"+modes[imode],"", 40, 0, 1); lout->Add(hpionmomentum); hhs.push_back("pionmomentum");
        TH1D * hpionEk = new TH1D("pionEk"+modes[imode],"", 40, 0, 0.5); lout->Add(hpionEk); hhs.push_back("pionEk");
        TH1D * hpiontheta = new TH1D("piontheta"+modes[imode],"", 35, 0, 180); lout->Add(hpiontheta); hhs.push_back("piontheta");
      }
      else{
        printf("\n\nunknown anaid!! %d\n\n", anaid); exit(1);
      }

      const int nhist = hhs.size();
      for(int ii=0; ii<nhist; ii++){
        printf("\nii %d hhs %s\n\n", ii, hhs[ii].Data());

        TString kincut;
        if(anaid==NUBAR1PI && (hhs[ii]!="Wrest")){
          kincut = "&& (Wrest<1.8)";
        }

        const TString tmpcut = Form("perweight * (%s %s %s) ", cuts[imode].Data(), basecut.Data(), kincut.Data());
        TString varname=hhs[ii];
        varname.ReplaceAll("VS",":");  
        varname.ReplaceAll("TIMES","*");
        cout<<"test "<<ii<<" cut "<<tmpcut<<" name "<<hhs[ii]+modes[imode]<<" varname "<<varname<<endl;
        const int ndraw = tree->Project(hhs[ii]+modes[imode],varname,tmpcut);
        cout<<"tmpcut "<<tmpcut<<" ndraw "<<ndraw<<endl<<endl;
      }
    }
  }

  const int nh = lout->GetEntries();

  //https://gibuu.hepforge.org/Documentation2017/code/analysis/neutrinoAnalysis_f90.html#robo1052
  /*
    For eventtype=5 and process_ID=CC and NC: 10^{-38} cm^2/GeV 
    All x-sec per particle (1/A) 
  */
  //GENIE spline same xsec unit as GiBUU
  /*
    root [8] tot_cc->GetYaxis()->GetTitle()
    (const char* 0x98a2358)"#sigma_{nuclear} (10^{-38} cm^{2})"
  */

  TList *l2d = new TList;
  for(int ii=0; ii<nh; ii++){
    //normalization confirmed by benchmark
    TH1D * htmp = (TH1D*) lout->At(ii);    
    const TString hname(htmp->GetName());
    if(hname.Contains("enu") || hname.Contains("Erecoil")){
      printf("For %s, no width normalization!\n", hname.Data());
      htmp->Scale(1E-38/histnormFactor);//scale up to full nucleus, so that can be added up with other nuclei
    }
    else{
      //htmp->Scale(1E-38*12/13/allnrun, "width");//per nucleon per bin width
      htmp->Scale(1E-38/histnormFactor, "width");//scale up to full nucleus, so that can be added up with other nuclei
    }

    //not needed any more and don't confuse with TH3D
    /*
    TH2D * h2d = (TH2D*) lout->At(ii);
    const TString n2d=h2d->GetName();
    if(n2d.Contains("VS")){
      getSliceYDrawX(h2d, l2d);
    }
    */
  }

  TFile * fout=new TFile(Form("outplot/outHist%s_A%d_ana%d_tag%s.root",nuExp.Data(), tarA, anaid, tag.Data()),"recreate");
  lout->Write();
  l2d->Write();
  fout->Save();
  fout->Close();

  cout<<"done!"<<endl;
}

int main(int argc, char * argv[])
{

  //void getHist(const TString fn, const TString tag, const int anaid)
  getHist(argv[1], argv[2], atoi(argv[3]));

  return 0;
}
