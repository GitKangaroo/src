#include "stdlib.h"

#include <iostream>
#include <fstream>

#include "TString.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TChain.h"
#include "AnaFunctions.h"
#include "AnaUtils.h"

using namespace std;

void GENIEReadChain(TChain * ch, TTree * tout, TH1F * &hCCrate, const int nEntryToStop = -999)
{
  const int _ARRAYSIZE_ = 100;

  // Declaration of leaf types                                                                                                                                                                                                                                       
  TObjString      *gEvtCode;

  Int_t           G2NeutEvtCode;
  Int_t           EvtNum;
  Double_t        EvtXSec;
  Double_t        EvtDXSec;
  Double_t        EvtWght;

  Double_t        EvtProb;
  Double_t        EvtVtx[4];
  Int_t           StdHepN;
  Int_t           StdHepPdg[_ARRAYSIZE_];
  Int_t           StdHepStatus[_ARRAYSIZE_];

  Int_t           StdHepRescat[_ARRAYSIZE_];
  Double_t        StdHepX4[_ARRAYSIZE_][4];
  Double_t        StdHepP4[_ARRAYSIZE_][4];
  Double_t        StdHepPolz[_ARRAYSIZE_][3];
  Int_t           StdHepFd[_ARRAYSIZE_];

  Int_t           StdHepLd[_ARRAYSIZE_];
  Int_t           StdHepFm[_ARRAYSIZE_];
  Int_t           StdHepLm[_ARRAYSIZE_];

  // List of branches                                                                                                                                                                                                                                                
  TBranch        *b_gEvtCode;

  TBranch        *b_G2NeutEvtCode;
  TBranch        *b_EvtNum;
  TBranch        *b_EvtXSec;
  TBranch        *b_EvtDXSec;
  TBranch        *b_EvtWght;

  TBranch        *b_EvtProb;
  TBranch        *b_EvtVtx;
  TBranch        *b_StdHepN;
  TBranch        *b_StdHepPdg;
  TBranch        *b_StdHepStatus;

  TBranch        *b_StdHepRescat;
  TBranch        *b_StdHepX4;
  TBranch        *b_StdHepP4;
  TBranch        *b_StdHepPolz;
  TBranch        *b_StdHepFd;

  TBranch        *b_StdHepLd;
  TBranch        *b_StdHepFm;
  TBranch        *b_StdHepLm;

  //only 0 works!!                                                                                                                                                                                                                                                   
  ch->SetMakeClass(0);
  gEvtCode = new TObjString;

  ch->SetBranchAddress("EvtCode", &gEvtCode, &b_gEvtCode);

  ch->SetBranchAddress("G2NeutEvtCode", &G2NeutEvtCode, &b_G2NeutEvtCode);
  ch->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
  ch->SetBranchAddress("EvtXSec", &EvtXSec, &b_EvtXSec);
  ch->SetBranchAddress("EvtDXSec", &EvtDXSec, &b_EvtDXSec);
  ch->SetBranchAddress("EvtWght", &EvtWght, &b_EvtWght);

  ch->SetBranchAddress("EvtProb", &EvtProb, &b_EvtProb);
  ch->SetBranchAddress("EvtVtx", EvtVtx, &b_EvtVtx);
  ch->SetBranchAddress("StdHepN", &StdHepN, &b_StdHepN);
  ch->SetBranchAddress("StdHepPdg", StdHepPdg, &b_StdHepPdg);
  ch->SetBranchAddress("StdHepStatus", StdHepStatus, &b_StdHepStatus);

  ch->SetBranchAddress("StdHepRescat", StdHepRescat, &b_StdHepRescat);
  ch->SetBranchAddress("StdHepX4", StdHepX4, &b_StdHepX4);
  ch->SetBranchAddress("StdHepP4", StdHepP4, &b_StdHepP4);
  ch->SetBranchAddress("StdHepPolz", StdHepPolz, &b_StdHepPolz);
  ch->SetBranchAddress("StdHepFd", StdHepFd, &b_StdHepFd);

  ch->SetBranchAddress("StdHepLd", StdHepLd, &b_StdHepLd);
  ch->SetBranchAddress("StdHepFm", StdHepFm, &b_StdHepFm);
  ch->SetBranchAddress("StdHepLm", StdHepLm, &b_StdHepLm);

  int ientry = 0;
  while(ch->GetEntry(ientry)){
    if(ientry%100000==0){
      printf("myEntries %d\n", ientry);
    }

    if(nEntryToStop>0){
      if(ientry>=nEntryToStop){
        printf("\n\n\n************************  GENIE Breaking after %d entries ***********************************************\n\n", nEntryToStop);
        break;
      }
    }

    //do it before the loop continues for any reason
    ientry++;

    //===========================================================================

    const TString ecode=gEvtCode->GetString();
    if(ecode.Contains("Weak[NC]")){//skip NC
      continue;
    }
    else if(!ecode.Contains("Weak[CC]")){
      printf("not cc!!! %s\n", ecode.Data()); exit(1);
    }

    if(!hCCrate){
      hCCrate = new TH1F("hCCrate","", 100000, 0, 100);//1MeV per bin
    }

    const double tmpenu = StdHepP4[0][3];
    hCCrate->Fill(tmpenu);

    const int tmpevent = EvtNum;
    const int tmpprod= abs(G2NeutEvtCode);
    const double tmppw = 1;//to-do EvtXSec;

    const int tmpnp = StdHepN;

    //printf("\ntest ientry %d tmpevent %d tmpprod %d tmpenu %f tmppw %f tmpnp %d\n", ientry, tmpevent, tmpprod, tmpenu, tmppw, tmpnp);
    
    AnaUtils::Ini();

    const bool isHydrogen = ecode.Contains("1000010010");//special case

    int idxIni = -999;
    int idxRESnucleon = -999;
    int idxRESpi = -999;
    int idxDelta = -999;
    for(int ii=1; ii<tmpnp; ii++){
      const int tmpid = StdHepPdg[ii];

      if( abs(tmpid)==13 && StdHepFm[ii]!=0 ){//skip non CC muon
        continue;
      }

      //https://github.com/luketpickering/NeutToRooTracker
      //If the StdHepStatus != 1 then you didn't see that particle. Always check the StdHepStatus.
      //Int_t StdHepStatus[$i<$StdHepN] //Status code of particle i --- 0 == incoming, 11 == incoming nucleus (if not using -E), 1 == good final state, anything else == intermediate or unseeable particle.
      AnaUtils::dtype IniOrFinaltype = AnaUtils::kNULL;
      if(isHydrogen){
        if(ii==1){//for hydrogen, initial nucleon at ii=1
          IniOrFinaltype = AnaUtils::kINI;
        }
      }
      else{
        if(StdHepStatus[ii]==11){//initial nucleon has status 11
          IniOrFinaltype = AnaUtils::kINI;
        }
      }

      if(IniOrFinaltype==AnaUtils::kINI){
        //tmpid can be 2000000201 for MEC
        if(!isHydrogen && StdHepFm[ii]!=1){
          printf("wrong initial nucleon event %d -- %d %d\n", tmpevent, StdHepFm[ii], tmpid); exit(1);
        }
      }

      //====================================================================>
      AnaUtils::dtype RESdtype = AnaUtils::kNULL;
      if(ecode.Contains("RES")){
        //particle list ordered by history; idx settings will be used for next loop
        if(IniOrFinaltype==AnaUtils::kINI){//initial nucleon
          idxIni = ii;
        }
        else if(//resonance
                StdHepFm[ii]==idxIni 
                && StdHepStatus[ii]!=13 //13 has no daughters
                ){
          
          if(StdHepStatus[ii]!=3){
            printf("\n\n\nGENIE RES resonance is not status 3! Strange!\n\n\n"); exit(1);
          }
          if( (StdHepLd[ii]-StdHepFd[ii])!=1 ){
            printf("\n\n\nGENIE RES resonance daughters not next to ech other! Strange!\n\n\n"); exit(1);
          }

          //cout<<"test event "<<tmpevent<<" ii "<<ii<<" idxIni "<<idxIni<<" ecode "<<ecode<<" status "<<StdHepStatus[ii]<<" pdg "<<StdHepPdg[ii]<<" scat "<<StdHepRescat[ii]<<" StdHepFd "<<StdHepFd[ii]<<" StdHepLd "<<StdHepLd[ii]<<" "<<endl;

          idxDelta = ii;

          const int tmppdg = StdHepPdg[StdHepFd[ii]];
          if(tmppdg>999){//nucleon
            idxRESnucleon = StdHepFd[ii];
            idxRESpi    = StdHepLd[ii];
          }
          else{
            idxRESpi    = StdHepFd[ii];
            idxRESnucleon = StdHepLd[ii];
          }
        }
        else if(ii == idxRESpi){
          RESdtype = AnaUtils::kPION;
        }
        else if(ii == idxRESnucleon){
          RESdtype = AnaUtils::kNUCLEON;
        }
      }
      //test-->
      /*      
      if(RESdtype>AnaUtils::kNULL){
        cout<<"test2 event "<<tmpevent<<" ii "<<ii<<" idxIni "<<idxIni<<" ecode "<<ecode<<" status "<<StdHepStatus[ii]<<" pdg "<<StdHepPdg[ii]<<" scat "<<StdHepRescat[ii]<<" StdHepFd "<<StdHepFd[ii]<<" StdHepLd "<<StdHepLd[ii]<<" RESdtype "<<RESdtype<<" "<<endl;
      }
      */
      //<--
      //<====================================================================

      //important to check status==1 for all cases
      if(StdHepStatus[ii]==1){
        IniOrFinaltype=AnaUtils::kFINAL;
      }

      //==== identify knock-out origin ===>
      double tmpKNsrc = 0;//proton 1; pion -1

      //printf("testjj %d %d %d %d\n", IniOrFinaltype, AnaUtils::kFINAL, RESdtype, AnaUtils::kNUCLEON);

      if(IniOrFinaltype==AnaUtils::kFINAL && StdHepPdg[ii]==2212 &&idxDelta>=0){//final-state proton
        int cursor = ii;
        int motherid=StdHepFm[cursor];
        
        while(motherid>idxDelta){
          cursor = motherid;
          motherid = StdHepFm[cursor];
        }
        if(motherid<idxDelta){
          printf("\n\n\ntmpKNsrc not found %d %d\n\n\n", motherid, idxDelta); exit(1);
        }
        else{
          if(StdHepPdg[cursor]>999){//nucleon
            tmpKNsrc = 1;
          }
          else{//pion
            tmpKNsrc = -1;
          }
        }
        //cout<<"test3 event "<<tmpevent<<" ii "<<ii<<" idxIni "<<idxIni<<" ecode "<<ecode<<" status "<<StdHepStatus[ii]<<" pdg "<<StdHepPdg[ii]<<" scat "<<StdHepRescat[ii]<<" StdHepFd "<<StdHepFd[ii]<<" StdHepLd "<<StdHepLd[ii]<<" RESdtype "<<RESdtype<<" idxDelta "<<idxDelta<<" tmpKNsrc "<<tmpKNsrc<<" "<<endl;     if(tmpevent==15)   exit(1);
      }
      //<====

      const bool isOK = (IniOrFinaltype!=AnaUtils::kNULL || RESdtype!=AnaUtils::kNULL);

      if(isOK){
        const double tmpmom1 = StdHepP4[ii][0];
        const double tmpmom2 = StdHepP4[ii][1];
        const double tmpmom3 = StdHepP4[ii][2];
        const double tmptote = StdHepP4[ii][3];

        //test
        //const TVector3 tmpvec(tmpmom1, tmpmom2, tmpmom3); printf("test particle %d %f %f %f %f %d mom %f theta %f\n", ii, tmpmom1, tmpmom2, tmpmom3, tmptote, tmpid, tmpvec.Mag(), tmpvec.Theta()*TMath::RadToDeg());

        //cout<<"test event "<<tmpevent<<" ii "<<ii<<" idxIni "<<idxIni<<" ecode "<<ecode<<" status "<<StdHepStatus[ii]<<" pdg "<<StdHepPdg[ii]<<" scat "<<StdHepRescat[ii]<<" StdHepFd "<<StdHepFd[ii]<<" StdHepLd "<<StdHepLd[ii]<<" RESdtype "<<RESdtype<<" "<<endl;
        
        AnaUtils::GENIEProceed(IniOrFinaltype, RESdtype, ecode, tmpevent, tmpprod, tmpenu, tmppw, tmpmom1, tmpmom2, tmpmom3, tmptote, tmpid, tmpKNsrc);
      }
    }
    AnaUtils::DoFill(tout);
  }
  cout<<"All entries "<<ientry<<endl;
}

int GiBUUReadFile(const TString filelist, TTree * tout, const int nFileToStop)
{
  TString tmpfl4target(filelist);
  tmpfl4target.ToUpper();
  int tmptargetZ = 6;
  if(tmpfl4target.Contains("HYDROGEN")){
    tmptargetZ = 1;
  }
  if(tmpfl4target.Contains("ARGON")){
    tmptargetZ = 18;
  }
    
  ifstream fl;
  fl.open(filelist);
  if(!fl.is_open()){
    printf("filelist %s not opened!\n", filelist.Data()); exit(1);
  }

  int filecount=0;
  int totlinecount=0;
  int totnrun = 0;

  string singlefile;
  while(getline(fl,singlefile)){
    if(filecount%10==0){
      cout<<"File count "<<filecount<<" total line count "<<totlinecount<<" total nrun "<<totnrun<<endl;
    }

    if(nFileToStop>0){
      if(filecount>=nFileToStop){
        printf("\n\n\n************************  GiBUU Breaking after %d files ***********************************************\n\n", nFileToStop);
        break;
      }
    }

    filecount++;
  //====================================== core loop ======================================
  ifstream elf;
  elf.open(singlefile.c_str());
  if(!elf.is_open()){
    printf("singlefile %s not opened!\n", singlefile.c_str()); exit(1);
  }
  string fline;
  AnaUtils::Ini();
  getline(elf, fline);//after reading even the last line it is still not eof
  while(!elf.eof()){
    TString tmpline(fline);
    totlinecount++;
    if(tmpline.Contains("Event")){//simply ignore it
      getline(elf, fline);
      continue;
    }
    
    double tmppw;
    double tmpmom1, tmpmom2, tmpmom3, tmptote;                                                                                                     
    double tmppos1, tmppos2, tmppos3;//position
    int tmprun, tmpevent, tmpid, tmpcharge, tmpprod;
    double tmpenu;
    
    sscanf(tmpline.Data(),"%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %*d %d %lf", &tmprun, &tmpevent, &tmpid, &tmpcharge, &tmppw, &tmppos1, &tmppos2, &tmppos3, &tmptote, &tmpmom1, &tmpmom2, &tmpmom3, /*&history,*/ &tmpprod, &tmpenu);
    
    AnaUtils::GiBUUProceed(tmprun, tmpevent, tmpcharge, tmpprod, tmpenu, tmppw, tmpmom1, tmpmom2, tmpmom3, tmptote, tmpid, tmptargetZ);

    getline(elf, fline);

    bool kChangeEvent = false;
    bool kChangeRun = false;
    if(elf.eof()){
      kChangeEvent = true;
      kChangeRun = true;
    }
    else{
    //only read in new event to check for event switching
      int nextRun, nextEvent;
      sscanf(fline.c_str(),"%d %d", &nextRun,  &nextEvent);
      if(nextEvent!=tmpevent){
        kChangeEvent = true;
      }
      if(nextRun!=tmprun){
        kChangeRun = true;
      }
    }

    if(kChangeEvent){
      AnaUtils::DoFill(tout);
      AnaUtils::Ini();
    }

    if(kChangeRun){
      totnrun++;
    }

  }
  elf.close();
  //============================================= end core ==============================================
  }

  cout<<"Final file count "<<filecount<<" total line count "<<totlinecount<<" total nrun "<<totnrun<<endl;
  return totnrun;
}

void anaGenerator(const TString tag, const TString filelist, const int tmpana, const int nToStop=-999)
{
  cout<<"please check "<<tag<<" "<<filelist<<" "<<tmpana<<endl;

  AnaUtils::experiment kExp = AnaUtils::kNONE;
  if(filelist.Contains("MINERvA")){
    kExp = AnaUtils::kMINERvA;
  }
  else if(filelist.Contains("T2K")){
    kExp = AnaUtils::kT2K;
  }
  else if(filelist.Contains("DUNE")){
    kExp = AnaUtils::kDUNE;
  }
  else{
    cout<<"\n\nnon-experimental filelist "<<filelist<<endl<<endl;
  }
    cout << "Now :do the TFile outplot\n";
    
  //TFile *fout=new TFile("outplot/outAna.root","recreate");
  TFile *fout=new TFile(Form("outplot/outAna%d_%s.root", tmpana, tag.Data()),"recreate");
  TTree * tout = AnaUtils::GetTree(tmpana, kExp);

  //_________________________________________________________________________________________________
  //_________________________________________________________________________________________________
  TChain * genieinput = AnaUtils::InputROOTFiles(filelist, "gRooTracker");
  TH1F * hCCrate = 0x0; 
  int nrun = -999;
  if(genieinput){
    GENIEReadChain(genieinput, tout, hCCrate, nToStop);
  }
  else{
    nrun = GiBUUReadFile(filelist, tout, nToStop);
  }

  //_________________________________________________________________________________________________
  //_________________________________________________________________________________________________

  fout->cd();

  if(hCCrate){
    hCCrate->Write();
  }

  TTree *theader = new TTree("header","header");
  theader->Branch("nrun",&nrun);
  theader->Fill();
  theader->Write();

  tout->Write();

  fout->Save();
  fout->Close();

  return;
}

int main(int argc, char* argv[])
{
  //void anaGenerator(const TString tag, const TString filelist, const int tmpana, const int nToStop)
  if(argc==4){
    anaGenerator(argv[1], argv[2], atoi(argv[3]));
  }
  else if(argc==5){
    anaGenerator(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
  }
  else{
    printf("wrong argc %d\n", argc); return 1;
  }
  return 0;
}
