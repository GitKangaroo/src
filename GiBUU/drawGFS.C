#include "TString.h"
#include "TH2D.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TSystem.h"
#include "TPaveStats.h"
#include "TLegend.h"

#include "include/style.h"

//test LOWRECOIL
//CugnonFSI
//CutR

enum{
  knuGFS0pi,
  knubarGFS0pi,
  knuGFS1pi,
  knubarGFS1pi,
  kGFSall,

  knuLOWRECOIL,
  knubarLOWRECOIL,
  knuLOWRECOILnoFSI,
  knubarLOWRECOILnoFSI,
  knuLOWRECOIL0piNp,

  knubarLOWRECOIL0piNp,
  kNUBAR1PI,
  kNUBAR1PInoFSI,
  kNUBAR1PIsigmaEnu,
  kNUBAR1PIsigmaEnunoFSI,

  knuLOWRECOILCutR,
  knubarLOWRECOILCutR
};

const Int_t cbbase = 1002;
const Int_t cbcol[]={cbbase+0, cbbase+1, cbbase+3, cbbase+2, cbbase+4, cbbase+5, cbbase+6, cbbase+7, 1009, 1007, 1005, 1003};

#define NO_OTHER 0

THStack * getTHStack(const TString filename, const TString histname)
{
  TFile *fin=new TFile(filename);
  THStack *stk = (THStack*)fin->Get(Form("%s", histname.Data()));

  if(!stk){
    TString tmph(histname);
    tmph.ReplaceAll("/","_T0/");
    stk = (THStack*)fin->Get(Form("%s", tmph.Data()));
    if(!stk){
      cout<<"getTHStack null! "<<filename<<" "<<histname<<endl; exit(1);
    }
  }

  cout<<"getTHStack filename "<<filename<<endl;
  cout<<"getTHStack histname "<<histname<<endl;
  return stk;
}

TH1D * getTH1D(const TString filename, const TString histname, const int color, const int linestyle, double &hmax)
{
  TH1D *htmp = style::GetStackedSum(getTHStack(filename, histname));

  style::ResetStyle(htmp);
  htmp->SetLineColor(color);
  htmp->SetLineWidth(3);
  htmp->SetLineStyle(linestyle);
  htmp->SetFillStyle(0);
  if(htmp->GetMaximum()>hmax){
    hmax = htmp->GetBinContent(htmp->GetMaximumBin());
  }
  

  return htmp;
}

TH1D * getMINERvA(const TString filename, const TString var, const Bool_t klist=0)
{
  TFile * fin = new TFile("external_files/"+filename);
  if(!fin->IsOpen()){
    printf("no data file!! %s\n", filename.Data()); exit(1);
  }
  TH1D * hh = 0x0;
  if(klist){
    TList *ll = (TList*) fin->Get(var);
    hh=(TH1D*) ll->At(0);
  }
  else{
    hh=(TH1D*) fin->Get(var);
  }
  if(!hh){
    printf("no MINERvA! %s %s\n", filename.Data(), var.Data()); //exit(1);
  }
  else{
    printf("getMINERvA got %s %s\n", hh->GetName(), hh->GetTitle());
  }
  return hh;
}

TString getLOWRECOILCut(const TString varname)
{
  int ib = atoi(varname(3,1).Data());
  if(varname.Contains("pn")){
    ib = atoi(varname(2,1).Data());
  }
  const double binlow[]={0, 0.2, 0.3, 0.4, 0.5, 0.6};
  const double binup[]={0.2, 0.3, 0.4, 0.5, 0.6, 0.8};
  return Form("%.1f < #it{q}_{3}/(GeV/#it{c}) < %.1f", binlow[ib], binup[ib]);
}

void SetLegend(const Int_t mode, TLegend *lg, const TString appendhead, const TString prependhead="")
{
  TString hder;
    if(mode==knuGFS0pi){
      hder = ("#nu0#piNp");
    }
    else if(mode==knuGFS1pi){
      hder = ("#nu1#piNp");
    }
    else if(mode==knubarGFS0pi){
      hder = ("#bar{#nu}0#piNp");
    }
    else if(mode==knubarGFS1pi){
      hder = ("#bar{#nu}1#piNp");
    }
    else if(mode==knuLOWRECOIL || mode==knuLOWRECOILCutR){
      hder = ("#nu");
    }
    else if(mode==knubarLOWRECOIL || mode==knubarLOWRECOILCutR){
      hder = ("#bar{#nu}");
    }
    else if(mode==knuLOWRECOILnoFSI){
      hder = ("#nu no FSI");
    }
    else if(mode==knubarLOWRECOILnoFSI){
      hder = ("#bar{#nu} no FSI");
    }
    else if(mode==knuLOWRECOIL0piNp){
      hder = ("CH #nu 0#piNp");
    }
    else if(mode==knubarLOWRECOIL0piNp){
      hder = ("CH #bar{#nu} 0#piNp");
    }
    else if(mode==kNUBAR1PI){
      hder = ("CH #bar{#nu}1#pi");
    }
    else if(mode==kNUBAR1PInoFSI){
      hder = ("CH #bar{#nu}1#pi no FSI");
    }
    else if(mode==kNUBAR1PIsigmaEnu){
      hder = ("CH #bar{#nu}1#pi");
    }
    else if(mode==kNUBAR1PIsigmaEnunoFSI){
      hder = ("CH #bar{#nu}1#pi no FSI");
    }

    hder+=" ";
    hder+=appendhead;

    hder.Prepend(prependhead);

    hder.ReplaceAll("/(GeV/#it{c})","");

    lg->SetHeader(hder);
}

TH1D* SetStack(THStack * &stk, const int mode, const int sumcol, const int lsty, const int lwid, const int fsty, TLegend *lg=0x0, const bool kshape=kFALSE, const TString tag="")
{
  const TString bul[]={"","","","","","","","","","","","","","","",""};//editor asks not to do so "(a) ","(a) ","(b) ","(b) "};
  static int gIentry = 0;

  TH1D * hsum =  style::GetStackedSum(stk, Form("%ssum", stk->GetName()), /*mode>kGFSall?kBlack:*/sumcol, lsty, lwid, fsty);
  double scale = -999;
  if(kshape){
    scale = hsum->Integral(0, 100000, "width");
    hsum->Scale(1./scale);
  }
  lg->AddEntry(hsum, bul[gIentry++]+tag, "l");

  bool keepbit[10]={0,0,0,0,0,0,0,0,0,0};
  if(NO_OTHER){
    for(int ii=0; ii<10; ii++) keepbit[ii]=1;
  }
  else{
  if(mode==knuGFS0pi){
    keepbit[0]=1;
  }
  else if(mode==knubarGFS0pi){
    keepbit[3]=1;
  }
  else if(mode==knuGFS1pi||mode==knubarGFS1pi){
    keepbit[1]=1;
  }
  else if(mode==knuLOWRECOIL || mode==knubarLOWRECOIL){
    keepbit[0]=1;
    keepbit[1]=1;
    keepbit[2]=1;
    keepbit[3]=1;
  }
  else if(mode==kNUBAR1PI){
    keepbit[1]=1;
  }
  else{
    for(int ii=0; ii<10; ii++) keepbit[ii]=1;
    keepbit[2]=0;
  }
  }

  TList * hgstack = stk->GetHists();
  TH1D * hRES = 0x0;
  TH1D * hDIS = 0x0;

  for(int ii=0; ii<hgstack->GetEntries(); ii++){
    TH1D * tmph=(TH1D*) hgstack->At(ii);
    TString tit = tmph->GetTitle();
    if(tit=="RES"){
      hRES = tmph;
      //test tit="RES+DIS";
    }
    else if(tit=="DIS"){
      hDIS = tmph;
    }

    tmph->SetFillStyle(fsty);
    tmph->SetLineStyle(kDashed);
    tmph->SetLineWidth(3);
    if(NO_OTHER){
      tmph->SetLineColor(cbcol[ii]);
      tmph->SetFillColor(cbcol[ii]);
    }
    else{
      tmph->SetLineColor(sumcol);
    }
    if(kshape){
      tmph->Scale(1./scale);
    }
    if(lg){
      if(keepbit[ii]){
        lg->AddEntry(tmph, bul[gIentry++]+tit, "l");
        cout<<"testtest "<<bul[gIentry++]+tit<<endl;
      }
    }
  }

  //test hRES->Add(hDIS);

  THStack * hnew = new THStack;
  for(int ii=0; ii<hgstack->GetEntries(); ii++){
    TH1D * tmph=(TH1D*) hgstack->At(ii);
    
    if(keepbit[ii]){
      hnew->Add(tmph);
    }
  }
  stk = hnew;
  return hsum;
}

void UpdateMax(TH1D* hobj, const TH1D* hnew)
{
  const double max1 = hnew->GetBinContent(hnew->GetMaximumBin())*1.1;
  const double max2 = hobj->GetMaximum();
  
  if(max1>max2){
    printf("test %e %e\n", max1, max2);
    hobj->SetMaximum(max1);
  }
}

TString getMainDir(TString mainname, const TString tag)
{
  if(tag.Contains("T2K")){
    mainname+="T2K";
  }
  else if(tag.Contains("MME")){
    mainname+="MME";//MINERvA ME
  }
  else{
    mainname+="MINERvA";
  }

  if(tag.Contains("NuWro")){
    mainname+="NuWro";
  }
  else if(tag.Contains("GENIEOOB")){
    mainname+="GENIEOOB";
  }
  else if(tag.Contains("GENIE")){
    mainname+="GENIE";
  }
  else{
    mainname+="GiBUUC";
  }

  if(tag.Contains("noFSI")){
    mainname+="noFSI";
  }
    
  return mainname;
}

TString getMainFile(TString mainname,const TString tag)
{
  mainname=getMainDir(mainname, tag);
  
  return "outStack/"+mainname+"/"+mainname+".root";
}

void SetMaximum(TH1D * hh, const int mode, const double hmax, const TString var, const TH1D * href, const Int_t ivar)
{
  hh->SetMaximum(hmax*1.1);

  if(var.Contains("eav")){
    if(mode==knuLOWRECOIL || mode==knuLOWRECOILnoFSI || mode==knuLOWRECOILCutR){
      hh->SetMaximum(6.3E-38);
    }
    else if(mode==knubarLOWRECOIL || mode==knubarLOWRECOILnoFSI || mode==knubarLOWRECOILCutR ){
      hh->SetMaximum(3E-38);
    }
  }
  else  if(var.Contains("q3")){
    if(href) hh->SetMaximum(href->GetBinContent(href->GetMaximumBin())*1.2);
  }
  
  /*
    else if(mode==knuLOWRECOIL0piNp){
    hh->SetMaximum(6.3E-39);
    }
    else if(mode==knubarLOWRECOIL0piNp){
    hh->SetMaximum(3E-39);
    }
  */
  else if(mode==kNUBAR1PI || mode==kNUBAR1PInoFSI || mode==kNUBAR1PIsigmaEnu || mode==kNUBAR1PIsigmaEnunoFSI){
    //"enu", "Wrest", "Q2", "muonmomentum","muontheta", "pionEk", "piontheta"//19-25
    //"enu", "Wrest", "Q2", "muonmomentum","muontheta", "pionmomentum", "pionEk", "piontheta"//26-33
    printf("\n\n ivar dependent index!!\n\n");
    double hmaxtmp[]={80E-40, hmax*1.1 /*7E-39*/, 7E-39, 1.5E-39, 3E-40, hmax*1.1, 1.2E-38, 5E-41};//start ivar==19
    hh->SetMaximum(hmaxtmp[ivar-27]);
  }
  
  hh->SetMinimum(0);
}

void SetRangeUser(TH1D * hh, const int mode, const TString var)
{
  double xmin=999, xmax=-999;
  if(mode<knuLOWRECOIL){
    if(var=="muontheta"){
      xmin = 0;
      //xmax = 20;
      xmax = 130;
    }
    else if(var=="protontheta"){
      xmin = 0;
      xmax = 70;
    }
    else if(var=="muonmomentum"){
      //xmin = 1.5;
      xmin = 0;
      xmax = 6;
    }
    else if(var=="protonmomentum"){
      xmin = 0.3;
      xmax = 1.2;
    }
    else if(var=="dalphat"){
      xmin = 0;
      xmax = 180;
    }
    else if(var=="dphit"){
      xmin = 0;
      xmax = 60;//180;//
      if(mode==knubarGFS0pi){
        xmax = 180;
      }
    }
    else if(var=="dpt"){
      xmin = 0;
      xmax = 1;
      if(mode==knubarGFS0pi){
        xmax = 1.2;
      }
    }
    else if(var=="neutronmomentum"){
      xmin = 0;
      xmax = 1;
      //hh->GetXaxis()->SetTitle("#it{p}_{n} (GeV/#it{c})");
      if(mode==knubarGFS0pi){
        xmax = 1.2;
      }
    }
    else if(var=="pionmomentum"){
      xmin = 0;
      xmax = 1.5;
    }
    else if(var=="piontheta"){
      xmin = 0;
      xmax = 70;
    }
    else if(var=="baryonmomentum"){
      xmin = 0.2;
      xmax = 4;
    }
    else if(var=="baryontheta"){
      xmin = 0;
      xmax = 70;
    }
    else if(var=="baryonmass"){
      xmin = 1;
      xmax = 2;
    }
    else if(var=="Wtrue"){
      xmin = 1;
      xmax = 2;
    }
    else if(var=="protonTT" || var=="pionTT" || var=="dpTT"){
      xmin = -0.4;
      xmax = 0.4;
    }
  }
  
  if( mode==knuLOWRECOIL0piNp || mode ==  knubarLOWRECOIL0piNp ){
    if(var.Contains("pn")){
      xmin = 0;
      xmax = 1.2;
    }
  }
  
  if( mode==kNUBAR1PI || mode==kNUBAR1PInoFSI || mode==kNUBAR1PIsigmaEnu || mode==kNUBAR1PIsigmaEnunoFSI){
    if(var.Contains("pionEk")){
      xmin = 0;
      xmax = 0.3;
    }
  }
  
  if(xmin<xmax){
    hh->GetXaxis()->SetRangeUser(xmin, xmax);
  }
}

TLegend * getLegend(const int opt, const int mode)
{
  const bool kleft = (opt==1 || opt==3 || opt==4 || opt==14);

  const double x0 = kleft? 0.15 : (/*mode>kGFSall?0.75:*/0.65)-0.05;
  double y0 = 0.6;

  if(mode==knuLOWRECOIL || mode ==  knubarLOWRECOIL || mode==knuLOWRECOIL0piNp || mode ==  knubarLOWRECOIL0piNp || mode==knuLOWRECOILnoFSI || mode ==  knubarLOWRECOILnoFSI || mode==knuLOWRECOILCutR || mode ==  knubarLOWRECOILCutR ){
    y0-=0.07;
    if(mode==knuLOWRECOIL || mode ==  knubarLOWRECOIL || mode==knuLOWRECOILnoFSI || mode ==  knubarLOWRECOILnoFSI || mode==knuLOWRECOILCutR || mode ==  knubarLOWRECOILCutR){
      //y0-=0.07;
      y0+=0.07;
    }
  }

  const double x1= x0+0.34;

  TLegend *lg = new TLegend(x0, y0, x1, 0.91);
  style::ResetStyle(lg,0.18,0.68); lg->SetTextAlign(12); //lg->SetBorderSize(1);

  return lg;
}

void drawGFS( const int mode, const int opt, TString tag)
{
  Bool_t kshape = (mode==4)? 1 : 0;
  //test
  if(tag.Contains("NuWro") && tag.Contains("T2K")){
  //if(tag.Contains("T2K")){
    kshape = mode==0? 0 : 1;
  }

  if(0/*mode>kGFSall*/){
    tag.ReplaceAll("MINERvA","");
    tag.ReplaceAll(" ","");

    if( mode==knuLOWRECOILCutR || mode ==  knubarLOWRECOILCutR ){
      tag.Append(" CutR");
    }
  }

  printf("===========> check kshape %d mode %d tag %s\n", kshape, mode, tag.Data());
  //====================================

  style::SetGlobalStyle();
  style::IniColorCB();
  Double_t currentLeft=0.14;
  Double_t currentTop=0.08;
  Double_t currentRight=0.035;
  Double_t currentBottom=0.18;
  style::fgkTextSize = /*mode>kGFSall?0.067:*/0.07;
  style::fgkTitleSize = /*mode>kGFSall?0.067:*/0.07;
  style::fgkYTitleOffset=0.9;
  style::fgkXTitleOffset = 1.15;

  TH1D * hnuGFS0pi=0x0;
  TH1D * hnuGFS1pi=0x0;
  TH1D * hnubarGFS0pi=0x0;
  TH1D * hnubarGFS1pi=0x0;
  THStack * hany = 0x0;
  TH1D * hdata = 0x0;

  THStack * hother = 0x0;

  const TString varname[]={"muonmomentum","muontheta","protonmomentum","protontheta",//0-3
                           "dalphat", "dphit","dpt", "neutronmomentum", //4-7
                           "pionmomentum", "pionEk", "piontheta", "baryonmomentum", "baryontheta", "baryonmass",//8-13
                           "q3", "eav0","eav1","eav2","eav3","eav4","eav5", "pn0","pn1","pn2","pn3","pn4","pn5",//14-26
                           "enu", "Wrest", "Q2", "muonmomentum","muontheta", "pionmomentum", "pionEk", "piontheta",//27-34
                           "protonTT", "pionTT", "dpTT"//35-37
  };
  const TString xtit[]={"#it{p}_{#mu} (GeV/#it{c})", "#theta_{#mu} (degree)", "#it{p}_{p} (GeV/#it{c})", "#theta_{p} (degree)",
                        "#delta#alpha_{T} (degree)", "#delta#phi_{T} (degree)", "#delta#it{p}_{T} (GeV/#it{c})", "#it{p}_{N} (GeV/#it{c})",
                        "#it{p}_{#pi} (GeV/#it{c})", "T_{#pi} (GeV)", "#theta_{#pi} (degree)", "#it{p}_{#tilde{N}} (GeV/#it{c})", "#theta_{#tilde{N}} (degree)", "W(baryonmass) (Gev/#it{c}^{2})",
                        "#it{q}_{3} (GeV/#it{c})", "#it{E}_{av} (GeV)", "#it{E}_{av} (GeV)","#it{E}_{av} (GeV)","#it{E}_{av} (GeV)","#it{E}_{av} (GeV)","#it{E}_{av} (GeV)", "#it{p}_{N} (GeV/#it{c})","#it{p}_{N} (GeV/#it{c})","#it{p}_{N} (GeV/#it{c})","#it{p}_{N} (GeV/#it{c})","#it{p}_{N} (GeV/#it{c})","#it{p}_{N} (GeV/#it{c})",
                        "E_{#nu} (GeV)", "W (Gev/#it{c}^{2})", "Q^{2} (GeV^{2})", "#it{p}_{#mu} (GeV/#it{c})", "#theta_{#mu} (degree)", "#it{p}_{#pi} (GeV/#it{c})", "T_{#pi} (GeV)", "#theta_{#pi} (degree)",
                        "#it{p}_{TT}^{p} (GeV/#it{c})", "#it{p}_{TT}^{#pi} (GeV/#it{c})", "#delta#it{p}_{TT} (GeV/#it{c})"
  };
  
  const TString ytit[]={"d#sigma/d#it{p}_{#mu} [cm^{2}/(GeV/#it{c}) nucleon]", "d#sigma/d#theta_{#mu} (cm^{2}/degree nucleon)", "d#sigma/d#it{p}_{p} [cm^{2}/(GeV/#it{c}) nucleon]", "d#sigma/d#theta_{p} (cm^{2}/degree nucleon)",
                        "d#sigma/d#delta#alpha_{T} (cm^{2}/degree nucleon)", "d#sigma/d#delta#phi_{T} (cm^{2}/degree nucleon)", "d#sigma/d#delta#it{p}_{T} [cm^{2}/(GeV/#it{c}) nucleon]", "d#sigma/d#it{p}_{N} [cm^{2}/(GeV/#it{c}) nucleon]",
                        "d#sigma/d#it{p}_{#pi} [cm^{2}/(GeV/#it{c}) nucleon]", "d#sigma/dT_{#pi} (cm^{2}/GeV nucleon)", "d#sigma/d#theta_{#pi} (cm^{2}/degree nucleon)", "d#sigma/d#it{p}_{#tilde{N}} [cm^{2}/(GeV/#it{c}) nucleon]", "d#sigma/d#theta_{#tilde{N}} (cm^{2}/degree nucleon)", "d#sigma/dW [cm^{2}/(GeV/#it{c})^{2} nucleon]",
                        "d#sigma/d#it{q}_{3} [cm^{2}/(GeV/#it{c}) nucleon]", "d^{2}#sigma/d#it{E}_{av}d#it{q}_{3} (cm^{2}/GeV^{2} nucleon)", "d^{2}#sigma/d#it{E}_{av}d#it{q}_{3} (cm^{2}/GeV^{2} nucleon)", "d^{2}#sigma/d#it{E}_{av}d#it{q}_{3} (cm^{2}/GeV^{2} nucleon)", "d^{2}#sigma/d#it{E}_{av}d#it{q}_{3} (cm^{2}/GeV^{2} nucleon)", "d^{2}#sigma/d#it{E}_{av}d#it{q}_{3} (cm^{2}/GeV^{2} nucleon)", "d^{2}#sigma/d#it{E}_{av}d#it{q}_{3} (cm^{2}/GeV^{2} nucleon)", "d^{2}#sigma/d#it{p}_{N}d#it{q}_{3}[cm^{2}/(GeV/#it{c}) GeV nucleon]","d^{2}#sigma/d#it{p}_{N}d#it{q}_{3}[cm^{2}/(GeV/#it{c}) GeV nucleon]","d^{2}#sigma/d#it{p}_{N}d#it{q}_{3}[cm^{2}/(GeV/#it{c}) GeV nucleon]","d^{2}#sigma/d#it{p}_{N}d#it{q}_{3}[cm^{2}/(GeV/#it{c}) GeV nucleon]","d^{2}#sigma/d#it{p}_{N}d#it{q}_{3}[cm^{2}/(GeV/#it{c}) GeV nucleon]","d^{2}#sigma/d#it{p}_{N}d#it{q}_{3}[cm^{2}/(GeV/#it{c}) GeV nucleon]",
                        "#sigma(E_{#nu}) (cm^{2}/nucleon)", "d#sigma/dW [cm^{2}/(GeV/#it{c})^{2} nucleon]", "d#sigma/dQ^{2} (cm^{2}/GeV^{2} nucleon)", "d#sigma/d#it{p}_{#mu} [cm^{2}/(GeV/#it{c}) nucleon]", "d#sigma/d#theta_{#mu} (cm^{2}/degree nucleon)", "d#sigma/d#it{p}_{#pi} [cm^{2}/(GeV/#it{c}) nucleon]", "d#sigma/dT_{#pi} (cm^{2}/GeV nucleon)", "d#sigma/d#theta_{#pi} (cm^{2}/degree nucleon)",
                        "d#sigma/d#it{p}_{TT}^{p} [cm^{2}/(GeV/#it{c}) nucleon]", "d#sigma/d#it{p}_{TT}^{#pi} [cm^{2}/(GeV/#it{c}) nucleon]", "d#sigma/d#delta#it{p}_{TT} [cm^{2}/(GeV/#it{c}) nucleon]",
};

  const int ivar = opt;

  double hmax = -999;
  const TString var = varname[ivar];

  TH1D * htmp = 0x0;

  TString legother;

  if(mode==kGFSall){
    if(ivar<8){
      hnuGFS0pi = getTH1D(getMainFile("GFS0pi", tag), getMainDir("GFS0pi", tag)+"nu/"+var, kBlue, kDashed, hmax);
      hnubarGFS0pi = getTH1D(getMainFile("GFS0pi", tag), getMainDir("GFS0pi", tag)+"nubar/"+var, kRed, kDashed, hmax);
    }
    hnuGFS1pi = getTH1D(getMainFile("GFS1pi", tag), getMainDir("GFS1pi", tag)+"nu/"+var, kBlue, kSolid, hmax);
    hnubarGFS1pi = getTH1D(getMainFile("GFS1pi", tag), getMainDir("GFS1pi", tag)+"nubar/"+var, kRed, kSolid, hmax);

    if(kshape){
      hmax = -999;
      TH1D *hhstmp[]={hnuGFS0pi, hnuGFS1pi, hnubarGFS0pi, hnubarGFS1pi};
      for(int ii=0; ii<4; ii++){
        if(hhstmp[ii]){
          hhstmp[ii]->Scale(1./hhstmp[ii]->Integral(0,10000, "width"));
          if(hhstmp[ii]->GetMaximum()>hmax){
            hmax = hhstmp[ii]->GetBinContent(hhstmp[ii]->GetMaximumBin());
          }
        }
      }
    }
    htmp = hnubarGFS1pi;//not important which one
  }
  else{
    legother=tag;
    if(legother.Contains("noFSI")){
      legother.ReplaceAll("noFSI","");
    }

    //test
    legother.ReplaceAll("GENIEOOB","GiBUU");
      cout << "legother: " << legother << endl;
      
      
    hnuGFS0pi = getTH1D(getMainFile("GFS0pi", tag), getMainDir("GFS0pi", tag)+"nu/"+var, kBlue, kDashed, hmax);

    //legother.ReplaceAll("MME","MINERvA");

    //test legother.ReplaceAll("GiBUU","NuWro");
    /*
    if(legother.Contains("T2K")){
      legother.ReplaceAll("T2K","MINERvA");
    }
    else{
      legother.ReplaceAll("GiBUU","NuWro");
    }
    */
/*
    if(mode==knuGFS0pi)    { hany = getTHStack(getMainFile("GFS0pi", tag), getMainDir("GFS0pi", tag)+"nu/"+var);
        cout << "get hany!!!" << endl;
        hother = getTHStack(getMainFile("GFS0pi", legother), getMainDir("GFS0pi", legother)+"nu/"+var);
        cout << "get hother!!!" << endl;
    }
    
    if(!hany){
      printf("\n\n\nhany null!!!!!!!!! %d %s\n\n\n\n", mode, var.Data()); exit(1);
    }
    style::ResetStyle(hany);

    htmp=style::GetStackedSum(hany); style::ResetStyle(htmp);
    hmax = htmp->GetBinContent(htmp->GetMaximumBin());
    
    /*
    //get data
    if( mode==knuLOWRECOIL || mode==knuLOWRECOILnoFSI || mode==knuLOWRECOILCutR) hdata = getMINERvA("nuLowRecoil1D.root", var);
    else if( mode==knubarLOWRECOIL || mode==knubarLOWRECOILnoFSI || mode==knubarLOWRECOILCutR) hdata = getMINERvA("nubarLowRecoil1D.root", var);
    else if( mode==kNUBAR1PI||mode==kNUBAR1PInoFSI||mode==kNUBAR1PIsigmaEnu||mode==kNUBAR1PIsigmaEnunoFSI) hdata = getMINERvA("MINERvANUBAR1PI_1906.08300.root", var);
    else if( mode==knuGFS0pi ) hdata = getMINERvA("MINERvA_1805.05486_v3.root", var, 1);
    else if( mode==knuGFS1pi ) hdata = getMINERvA("GFSPIZERO_Coplowe_20191009.root", var);
    */
    
  }

  //_________________________________________________________________________________
  //_________________________________________________________________________________

  TCanvas *cc=0x0;
  cc=new TCanvas("cc","",600,400);

  style::PadSetup(cc);
  //cc->GetCanvas()->SetGrayscale();
  gPad->SetLeftMargin(currentLeft);//square
  gPad->SetRightMargin(currentRight);//0.03//square
  gPad->SetTopMargin(currentTop);
  gPad->SetBottomMargin(currentBottom);

  TH1D * hfirstdraw = (TH1D*) htmp->Clone("firstdraw");
  hfirstdraw->Scale(0);
  SetMaximum(hfirstdraw, mode, hmax, varname[ivar], hdata, ivar);
  SetRangeUser(hfirstdraw, mode, varname[ivar]);
  hfirstdraw->SetXTitle(xtit[ivar]);
  hfirstdraw->SetYTitle(kshape?"p.d.f.":ytit[ivar]);

  TLegend * lg = getLegend(opt, mode);

  //======================================================================================
  //======================================================================================
    
  if(hany){
    //hfirstdraw->Scale(0);

    TString appendhead;
    if(var.Contains("eav")||var.Contains("pn")){
      appendhead = getLOWRECOILCut(varname[ivar]);
    }

    TString prependhead;
    if(tag.Contains("GiBUU")){
      if(var.Contains("dalphat")){
        prependhead="(b) ";
      }
      else if(var.Contains("neutronmomentum")){
        prependhead="(a) ";
      }
    }
    else{
      if(var.Contains("dalphat")){
        prependhead="(d) ";
      }
      else if(var.Contains("neutronmomentum")){
        prependhead="(c) ";
      }
    }
    SetLegend(mode,lg, appendhead, prependhead);
    


    //drawing...

    style::ResetStyle(hfirstdraw);
      cout << "1draw" << endl;
    hfirstdraw->Draw("sameaxis");
    //if(varname[ivar].Contains("q3") || varname[ivar].Contains("eav")){
    //hany->Draw("same hist");
    //style::ResetStyle(hany);
      cout << "draw hany" << endl;
    //hany->Draw(/*mode>kGFSall?"same nostack hist":*/"same nostack hist");
      hnuGFS0pi->Draw("hist same");
    //hany->Draw(/*mode>kGFSall?"same nostack hist":*/"same nostack hist");
    //style::ResetStyle(hanySum);
    //hanySum->SetLineWidth(3);
    //  cout << "draw hanySum" << endl;
    //hanySum->Draw(/*mode>kGFSall?"same C hist":*/"same hist");
    //hanySum->Draw(/*mode>kGFSall?"same C hist":*/"same hist");


    /*
      }
      else{
      hany->Draw("same hist C");
      hh->Draw("same hist C");
      }
    */
    
    //printf("hanySum total xsec %s %e\n", hanySum->GetName(), hanySum->Integral(0,100000,"width"));
/*
    if(hother){
      style::ResetStyle(hother);
      //hother->SetLineWidth(3);
        cout << "draw hother" << endl;
      hother->Draw("same nostack hist C");
      //hother->Draw("same nostack hist");
      style::ResetStyle(hanySum);
      hotherSum->SetLineWidth(3);
      hotherSum->Draw("same hist C");
      //hotherSum->Draw("same hist");
      
      printf("hotherSum total xsec %s %e\n", hotherSum->GetName(), hotherSum->Integral(0,100000,"width"));
    }
    */
  }
  //======================================================================================
  //======================================================================================
  else if(mode==kGFSall){
    lg->SetHeader(tag);

    hfirstdraw->Draw("sameaxis");

    if(hnuGFS0pi){
      lg->AddEntry(hnuGFS0pi,"#nu0#piNp", "l");
      hnuGFS0pi->Draw("same hist C");
    }

    lg->AddEntry(hnuGFS1pi,"#nu1#piNp", "l");
    hnuGFS1pi->Draw("same hist C");

    if(hnubarGFS0pi){
      lg->AddEntry(hnubarGFS0pi,"#bar{#nu}0#piNp", "l");
      hnubarGFS0pi->Draw("same hist C");
    }

    lg->AddEntry(hnubarGFS1pi,"#bar{#nu}1#piNp", "l");
    hnubarGFS1pi->Draw("same hist C");
  }
  //======================================================================================
  //======================================================================================
    
  lg->Draw();
  hfirstdraw->SetLineWidth(3);
  hfirstdraw->Draw("sameaxis");

  tag.ReplaceAll(" ","");
  cc->Print(Form("outplot/final/GFSmode%dopt%dkshape%d%s.eps", mode, opt, kshape, tag.Data()));
  cc->Print(Form("outplot/final/GFSmode%dopt%dkshape%d%s.png", mode, opt, kshape, tag.Data()));
}

int main(int argc, char * argv[])
{
  //void drawGFS( const int mode, const int opt, TString tag)
  drawGFS(atoi(argv[1]), atoi(argv[2]), argv[3]);
  return 0;
}

