#ifndef _ANAUTILS_
#define _ANAUTILS_

using namespace std;

using namespace AnaFunctions;

class AnaUtils
{
 public:
  enum experiment{
    kNONE, //kNONE = 0
    kMINERvA, //kMINERvA = 1
    kT2K, //KT2K = 2
    kDUNE //kDUNE = 3
  };

  enum dtype{
    kNULL,
    kPION,
    kNUCLEON,
    kINI,
    kFINAL
  };

  static TTree * GetTree(const int ana, const experiment exp);
  static void DoFill(TTree *tt);
  static void Ini();
  //GiBUU Proceed
  static void GiBUUProceed(const int tmprun, const int tmpevent, const int tmpcharge, const int tmpprod, const double tmpenu, const double tmppw, const double tmpmom1, const double tmpmom2, const double tmpmom3, const double tmptote, const int tmpid, const int tmptargetZ);
  //GENIE Proceed
  static void GENIEProceed(const dtype IniOrFinaltype, const dtype RESdtype, const TString code, const int tmpevent, const int tmpprod, const double tmpenu, const double tmppw, const double tmpmom1, const double tmpmom2, const double tmpmom3, const double tmptote, const int tmpid, const double tmpKNsrc);

  static TChain * InputROOTFiles(const TString file, const TString tr);

 private:
  enum analysis{
    GFS       =0,  
    LOWRECOIL =1,
    NUBAR1PI  =2,
    GFSEXP    =3,
    MMECCQE   =4,
    CLR       =5,
    RESPS     =6,
    GFSPIZERO =7
  };
  
  enum particlebit{
    PROTONBIT  =1,
    PIONBIT    =10,
    MUONBIT    =100,
    KAONBIT    =1000,
    ELECTRONBIT=10000,
    PIZEROBIT  =100000,
    GAMMABIT   =1000000,
    NEUTRONBIT =10000000,
    BKGBIT     =100000000
  };

  enum mode{
    kALL=0,
    kQE=1,
    kRES=2,
    kDIS=3,
    k2P2H=4,
    kOTHER=5
  };

  //Setup
  static void SetGENIEMode(const TString code);
  static void SetGENIETarget(const TString evtcode);
  static void GENIESetID(const int pdg, const double tmptote);

  static bool GiBUUIsIniN(const double tmppw, const double tmpmom1, const double tmpmom2, const double tmpmom3, const double tmptote);
  static void GiBUUSetID(const int id, const double tmptote);
  
  static void MainProceed();

  //Proceed
  static void AddEav(const bool kprint=false);
  static void ProceedT2KGFS();
  static void ProceedT2KGFSEXP();
  static void ProceedMINERvAGFS();
  static void ProceedMINERvAGFSPIZERO();
  static void ProceedMINERvAMECCQE();
  static void ProceedMINERvALOWRECOIL();
  static void ProceedMINERvANUBAR1PI();
  static void ProceedCLR();
  static void ProceedRESPS();

  //DoFill
  static bool IsMuon();
  static bool IsProton();
  static bool IsNeutron();
  static bool IsPion();
  static bool IsKaon();
  static bool IsElectron();
  static bool IsPiZero();
  static bool IsGamma();
  static bool IsBKG();
  static int GetNMuons();
  static int GetNProtons();
  static int GetNNeutrons();
  static int GetNPions();
  static bool IsGood();
  static void Calc();
  static void SetGiBUUMode();

  //---------------------------  
  //set by ReadLine->
  static int run;
  static int event;
  static int targetZ;
  static int prod;
  static int evtMode;
  static double perweight;
  static double KNsrc;
  static double enu;
  static int lineCharge;
  static int globalMuonCharge;
  static TLorentzVector *lineFullMom;
  static double lineMass;
  static double lineKNsource;
  static int linePID;
  static bool lineIsBkgParticle;
  static bool lineIsMMECCQEBkg;

  //<-

  static int anamode;
  static experiment expmode;

  //need to Ini-->
  static int npar;
  static int parbit;
  static TLorentzVector * iniNfullp;
  static TLorentzVector * RESpifullp;
  static TLorentzVector * RESnucleonfullp;
  static TLorentzVector * muonfullp;
  static TLorentzVector * protonfullp;
  static TLorentzVector * neutronFSfullp;
  static TLorentzVector * pionfullp;
  static TLorentzVector * baryonfullp;
  static double Eav;
  //<--

  //set via Calc()
  static double muoncostheta;
  static double muonmomentum;
  static double muontheta;
  static double protonmomentum;
  static double protontheta;
  static double neutronFSmomentum;
  static double neutronFStheta;
  static double pionmomentum;
  static double piontheta;
  static double pionEk;
  static double baryonmomentum;
  static double baryontheta;
  static double baryonmass;
  static double dpt;
  static double dphit;
  static double dalphat;
  static double dpTT;
  static double protonTT;
  static double pionTT;
  static double dpTy;
  static double neutronmomentum;
  static double calcEnu;
  static double q3;
  static double Wrest;
  static double Wtrue;
  static double xBj;
  static double xCCH;
  static double energyCCH;
  static double xrest;
  static double Q2;
  static double q2qe;
  static double Erecoil;
  static double muonpt;
  static double mupz;
  static double RESmass;
  static double adlerPhi;
  static double pseudoPhi;
  static double cosNuIniNAngle;
  static double cosQIniNAngle;
  static double lrsign;
  static double w2;
  static double wpseudo2;
  static double pseudosign;
};
//Initialise variables
int AnaUtils::prod=-999;
int AnaUtils::evtMode=-999;
int AnaUtils::run=-999;
int AnaUtils::event=-999;
int AnaUtils::targetZ=-999;
double AnaUtils::perweight=-999;
double AnaUtils::KNsrc = -999;
double AnaUtils::enu=-999;
int AnaUtils::lineCharge=-999;
int AnaUtils::globalMuonCharge=-999;

TLorentzVector * AnaUtils::lineFullMom = new TLorentzVector;
double AnaUtils::lineMass = -999;
double AnaUtils::lineKNsource = -999;
int AnaUtils::linePID = -999;
bool AnaUtils::lineIsBkgParticle = false;
bool AnaUtils::lineIsMMECCQEBkg = false;

int AnaUtils::anamode=0;
AnaUtils::experiment AnaUtils::expmode=AnaUtils::kNONE;
int AnaUtils::npar=0;
int AnaUtils::parbit=0;

TLorentzVector * AnaUtils::iniNfullp = new TLorentzVector;
TLorentzVector * AnaUtils::RESpifullp = new TLorentzVector;
TLorentzVector * AnaUtils::RESnucleonfullp = new TLorentzVector;
TLorentzVector * AnaUtils::muonfullp = new TLorentzVector;
TLorentzVector * AnaUtils::protonfullp = new TLorentzVector;
TLorentzVector * AnaUtils::neutronFSfullp = new TLorentzVector;
TLorentzVector * AnaUtils::pionfullp = new TLorentzVector;
TLorentzVector * AnaUtils::baryonfullp = new TLorentzVector;
double AnaUtils::Eav=-999;

double AnaUtils::muoncostheta=-999;
double AnaUtils::muonmomentum=-999;
double AnaUtils::muontheta=-999;
double AnaUtils::protonmomentum=-999;
double AnaUtils::protontheta=-999;
double AnaUtils::neutronFSmomentum=-999;
double AnaUtils::neutronFStheta=-999;
double AnaUtils::pionmomentum=-999;
double AnaUtils::piontheta=-999;
double AnaUtils::pionEk=-999;
double AnaUtils::baryonmomentum=-999;
double AnaUtils::baryontheta=-999;
double AnaUtils::baryonmass=-999;
double AnaUtils::dpt=-999;
double AnaUtils::dphit=-999;
double AnaUtils::dalphat=-999;
double AnaUtils::dpTT=-999;
double AnaUtils::protonTT=-999;
double AnaUtils::pionTT=-999;
double AnaUtils::dpTy=-999;
double AnaUtils::neutronmomentum=-999;
double AnaUtils::calcEnu=-999;
double AnaUtils::q3=-999;
double AnaUtils::Wrest=-999;
double AnaUtils::Wtrue=-999;
double AnaUtils::xBj=-999;
double AnaUtils::xCCH=-999;
double AnaUtils::energyCCH=-999;
double AnaUtils::xrest=-999;
double AnaUtils::Q2=-999;
double AnaUtils::q2qe=-999;
double AnaUtils::Erecoil=-999;
double AnaUtils::muonpt=-999;
double AnaUtils::mupz=-999;
double AnaUtils::RESmass=-999;
double AnaUtils::adlerPhi=-999;
double AnaUtils::pseudoPhi=-999;
double AnaUtils::cosNuIniNAngle=-999;
double AnaUtils::cosQIniNAngle=-999;
double AnaUtils::lrsign=-999;
double AnaUtils::w2=-999;
double AnaUtils::wpseudo2=-999;
double AnaUtils::pseudosign=-999;

//Declare Functions
TChain * AnaUtils::InputROOTFiles(const TString file, const TString tr)
{
  TChain *ch=0x0;

  if(file.Contains(".root")){
    ch = new TChain(tr);
    ch->Add(file);
  }
  else{
    ifstream fin(file);
    if(!fin){
      printf("AnaUtils::InputROOTFiles file not found \n%s\n\n",file.Data()); exit(1);
    }

    TString buff;
    while(fin.good()){
      fin>>buff;
      if(buff!=""){
        if(!buff.Contains(".root")){
          return 0x0;
        }
        else{
          if(!ch){
            ch = new TChain(tr);
          }
          ch->Add(buff);
        }
      }
    }
  }

  //const Int_t ent = ch->GetEntries(); //takes infinity time!!
  printf("\t%d trees!\n",ch->GetNtrees());

  return ch;
}

void AnaUtils::SetGENIEMode(const TString code)
{
  evtMode = kOTHER;

  if(code.Contains("QES")){
    evtMode = kQE;
  }
  else if(code.Contains("RES")){
    evtMode = kRES;
  }
  else if(code.Contains("DIS")){
    evtMode = kDIS;
  }
  else if(code.Contains("MEC")){
    evtMode = k2P2H;
  }

  return;

  //ignroe the following

  //genie mode from NEUT code
  //https://internal.dunescience.org/doxygen/neut__code__from__rootracker_8C.html

  //if      (is_qel && !is_charm && is_cc && is_nu           ) evtype =   1;
  //if (is_dis && W_gt_2 && is_cc && is_nu   ) evtype =  26;
  /*
 else if (is_nu    && is_cc && is_p && np==1 && nn==0 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  11;
  462        else if (is_nu    && is_cc && is_n && np==1 && nn==0 && npip==0 && npim==0 && npi0==1 && neKL==0) evtype =  12;
  463        else if (is_nu    && is_cc && is_n && np==0 && nn==1 && npip==1 && npim==0 && npi0==0 && neKL==0) evtype =  13;
   */

  //all work-in-progress, can't be tursted yet!!
  const int outmode[]={kQE, kRES, kDIS, k2P2H};
  const int pdsmin[] ={1, 11, 26, 35};
  const int pdsmax[] ={1, 13, 26, 35};
  const Int_t nmode = sizeof(outmode)/sizeof(int);

  //prod has to be unsigned!
  evtMode=kOTHER;
  for(int ii=0; ii<nmode; ii++){
    if(prod>=pdsmin[ii] && prod<=pdsmax[ii]){
      evtMode=outmode[ii];
      return;
    }
  }

}

void AnaUtils::SetGENIETarget(const TString evtcode)
{
    targetZ = -999;

    if(evtcode.Contains("1000060120")){
      targetZ = 6;
    }
    else if(evtcode.Contains("1000010010")){
      targetZ = 1;
    }
    else if(evtcode.Contains("1000180400")){
      //argon
      //Atomic number18
      //Standard atomic weight (±)39.948(1)[1]
      targetZ = 18;
    }
    else if(evtcode.Contains("1000080160")){
      //oxygen
      targetZ = 8;
    }
    else if(evtcode.Contains("1000260560")){
      //Fe
      //Atomic number26
      //Standard atomic weight (±)55.845(2)[1]
      targetZ = 26;
    }
    else if(evtcode.Contains("1000822070")){
      //Pb
      //Atomic number82
      //Standard atomic weight (±)207.2(1)[1]
      targetZ = 82;
    }
    else{
      printf("unknown tmptargetZ @%s@\n", evtcode.Data());exit(1);
    }
}

void AnaUtils::GENIESetID(const int pdg, const double tmptote)
{
  //GiBUU: tmpid, lineCharge -> lineIsBkgParticle, lineMass, linePID, globalMuonCharge, parbit      
  //GENIE: pdg -> lineCharge,   lineIsBkgParticle, lineMass, linePID, globalMuonCharge, parbit

  const int tmppdg = abs(pdg);
  ///minerva/app/users/xlu/cmtuser/Minerva_v10r8p12/Ana/AnalysisFramework/External/GENIEXSecExtract/apps/runXSecLooper_transverseCCQETwoTrack.cpp  
  /*
          if( apd == 13 ) genie_n_muons++;
           else if( apd == 22 || apd == 11 )  genie_n_photons++;
           else if( apd == 2212 ) { genie_n_protons++;  genie_n_nucleons++; }
           else if( apd == 2112 ) { genie_n_neutrons++; genie_n_nucleons++; }
           else if( apd == 211 )  { genie_n_pions++;    genie_n_mesons++; }
           else if( apd == 111 )  { genie_n_pi_zeros++; genie_n_photons++; }
           else if( apd == 411 || apd == 421 || apd == 431 ) { genie_n_charms++; genie_n_mesons++; }
           else if( apd == 321 || apd == 311 || apd == 310 || apd == 130 ) { genie_n_kaons++; genie_n_mesons++; }
           else if( apd > 3000 && apd < 5000 ) genie_n_heavy_baryons++;
           else genie_n_others++;

         }

         if( genie_n_muons         == 1 &&
             genie_n_mesons        == 0 &&
             genie_n_heavy_baryons == 0 &&
             genie_n_photons       == 0 &&
             genie_n_protons       != 0  ) {
           signal = true;
         }
   */
  lineIsBkgParticle = false;
  //two quarks (meson) 3 digits http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf  The general form is a 7–digit number: ±n nr nL nq1 nq2 nq3 nJ
  if( 
     (tmppdg>99 && tmppdg<1000) || 
     ( tmppdg == 22 || tmppdg == 11 ) ||
     ( tmppdg > 3000 )
      ){
    lineIsBkgParticle = true;
  }

  //allow electron and photon below 10MeV
  lineIsMMECCQEBkg = false;
  if(
     (tmppdg>99 && tmppdg<1000) ||
     ( tmppdg == 22 && tmptote > 0.01 ) ||
     ( tmppdg > 3000 )
     ){
    lineIsMMECCQEBkg = true;
  }

  lineCharge = (pdg>0? 1 : -1);

  linePID = 0;
  lineMass = 0;
  if(tmppdg==13){
    lineCharge *= -1;
    linePID = MUONBIT;
    lineMass = MuonMass();

    //no need to reset muon charge for event, should be the same for the whole sample                                                                                                                                                                       
    if(globalMuonCharge == -999){
      globalMuonCharge = lineCharge;
    }
    else if(globalMuonCharge != lineCharge){
      cout<<"\n\n ************** Muon charge not consistent! "<<globalMuonCharge<<" "<<lineCharge<<endl<<endl;
      exit(1);
    }
  }
  else if(tmppdg==11){
    lineCharge *= -1;
    linePID = ELECTRONBIT;
    lineMass = ElectronMass();
    parbit += ELECTRONBIT;
  }
  else if(tmppdg==2212){
    linePID = PROTONBIT;
    lineMass = ProtonMass();
    parbit += PROTONBIT;
  }
  else if(tmppdg==211){
    linePID = PIONBIT;
    lineMass = PionMass();
    parbit += PIONBIT;
  }
  else if(tmppdg==321){
    linePID = KAONBIT;
    lineMass = KaonMass();
    parbit += KAONBIT;
  }
  //neutrals now
  else if(tmppdg==111){
    lineCharge = 0;
    linePID = PIZEROBIT;
    lineMass = PiZeroMass();
    parbit += PIZEROBIT;
  }
  else if(tmppdg==130||tmppdg==310||tmppdg==311){
    lineCharge = 0;
    linePID = KAONBIT;
    parbit += KAONBIT;
  }
  else if(tmppdg==2112){
    lineCharge = 0;
    linePID = NEUTRONBIT;
    lineMass = NeutronMass();
    parbit += NEUTRONBIT;
  }
  else if(tmppdg==22){
    lineCharge = 0;
    linePID = GAMMABIT;
    parbit += GAMMABIT;
  }
}

void AnaUtils::GENIEProceed(const dtype IniOrFinaltype, const dtype RESdtype, const TString code, const int tmpevent, const int tmpprod, const double tmpenu, const double tmppw, const double tmpmom1, const double tmpmom2, const double tmpmom3, const double tmptote, const int tmpid, const double tmpKNsrc)
{
  if(IniOrFinaltype==kINI){
    iniNfullp->SetXYZT(tmpmom1, tmpmom2, tmpmom3, tmptote);
  }
  else{
    //both RESpi and RESnucleon can be final-state particles
    if(RESdtype==kPION){
      RESpifullp->SetXYZT(tmpmom1, tmpmom2, tmpmom3, tmptote);
    }
    else if(RESdtype==kNUCLEON){
      RESnucleonfullp->SetXYZT(tmpmom1, tmpmom2, tmpmom3, tmptote);
    }

    if(IniOrFinaltype==kFINAL){
      lineFullMom->SetXYZT(tmpmom1, tmpmom2, tmpmom3, tmptote);
      lineKNsource = tmpKNsrc;
            
      SetGENIEMode(code);
      SetGENIETarget(code);
      
      event = tmpevent;
      prod = tmpprod;
      enu = tmpenu;
      
      perweight = tmppw;

      GENIESetID(tmpid, tmptote);
      MainProceed();
    }
  }

  //lineCharge,   lineIsBkgParticle, lineMass, linePID, globalMuonCharge, parbit
  //printf("======> test evtMode %d targetZ %d event %d prod %d enu %f perweight %f --- lineCharge %d lineIsBkgParticle %d lineMass %f linePID %d globalMuonCharge %d parbit %d -- npar %d\n", evtMode, targetZ, event, prod, enu, perweight, lineCharge, lineIsBkgParticle, lineMass, linePID, globalMuonCharge, parbit, npar);
  /*
  cout<<"Muon"<<endl;
  muonfullp->Print();
  cout<<"Proton"<<endl;
  protonfullp->Print();
  cout<<"Pion"<<endl;
  pionfullp->Print();
  */
}

//============================================================================================================================================================
//============================================================================================================================================================

void AnaUtils::SetGiBUUMode()
{
  //const TString modes[]={"all","qe","res","dis","2p2h", "other"};
  //const TString cuts[]={"1","(prod==1)", "(prod>=2 && prod<=33)", "(prod==34)", "(prod==35)", "(prod>=36)"};
  const int outmode[]={kQE, kRES, kDIS, k2P2H};
  const int pdsmin[] ={1, 2,  34, 35};
  const int pdsmax[] ={1, 33, 34, 35};
  const Int_t nmode = sizeof(outmode)/sizeof(int);

  evtMode=kOTHER;
  for(int ii=0; ii<nmode; ii++){
    if(prod>=pdsmin[ii] && prod<=pdsmax[ii]){
      evtMode=outmode[ii];
      return;
    }
  }
}

void AnaUtils::GiBUUProceed(const int tmprun, const int tmpevent, const int tmpcharge, const int tmpprod, const double tmpenu, const double tmppw, const double tmpmom1, const double tmpmom2, const double tmpmom3, const double tmptote, const int tmpid, const int tmptargetZ)
{
  if(GiBUUIsIniN(tmppw, tmpmom1, tmpmom2, tmpmom3, tmptote)){
    return;
  }
  else{
    targetZ = tmptargetZ;

    run = tmprun;
    event = tmpevent;
    lineCharge = tmpcharge;
    prod = tmpprod;
    SetGiBUUMode();
    enu = tmpenu;

    GiBUUSetID(tmpid, tmptote);
    
    //have to stay after IsIniN, because should not do it if it is initial
    MainProceed();
  }
}

void AnaUtils::GiBUUSetID(const int tmpid, const double tmptote)
{
  //tmpid, lineCharge -> lineIsBkgParticle, lineMass, linePID, globalMuonCharge, parbit

  lineIsBkgParticle = false;
  ///minerva/app/users/xlu/cmtuser/Minerva_v10r8p12/Ana/AnalysisFramework/External/GENIEXSecExtract/apps/runXSecLooper_transverseCCQETwoTrack.cpp
  /*
    const int apd = fabs(pdg);
           if( apd == 13 ) genie_n_muons++;
           else if( apd == 22 || apd == 11 )  genie_n_photons++;
           else if( apd == 2212 ) { genie_n_protons++;  genie_n_nucleons++; }
           else if( apd == 2112 ) { genie_n_neutrons++; genie_n_nucleons++; }
           else if( apd == 211 )  { genie_n_pions++;    genie_n_mesons++; }
           else if( apd == 111 )  { genie_n_pi_zeros++; genie_n_photons++; }
           else if( apd == 411 || apd == 421 || apd == 431 ) { genie_n_charms++; genie_n_mesons++; }
           else if( apd == 321 || apd == 311 || apd == 310 || apd == 130 ) { genie_n_kaons++; genie_n_mesons++; }
           else if( apd > 3000 && apd < 5000 ) genie_n_heavy_baryons++;
           else genie_n_others++;

         }

         if( genie_n_muons         == 1 &&
             genie_n_mesons        == 0 &&
             genie_n_heavy_baryons == 0 &&
             genie_n_photons       == 0 &&
             genie_n_protons       != 0  ) {
           signal = true;
         }
   */
  //all mesons (pi, ka, charmed, all in GiBUU meson) + photon (including electron) + all heavy baryons (all s, c, b) , https://gibuu.hepforge.org/trac/wiki/ParticleIDs
  if( 
     ( tmpid>100 && tmpid<200 ) || 
     ( tmpid==999 || tmpid==901 ) || 
     ( tmpid>=32 && tmpid<=61 ) 
      ){
    lineIsBkgParticle = true;
  }
  
  //allow electron and photon below 10MeV
  lineIsMMECCQEBkg = false;
  if(
     ( tmpid>100 && tmpid<200 ) ||
     ( tmpid==999 && tmptote > 0.01 ) ||
     ( tmpid>=32 && tmpid<=61 )
     ){
    lineIsMMECCQEBkg = true;
  }

  lineMass = 0;
  linePID = 0;
  ////E_av = \sum KE_p + \sum KE_piplus + \sum KE_piminus + sum E_electrons + \sum E_Kplus + \sum E_kminus + sum E_pizero + sum E_gamma  
  if(lineCharge){
    if(tmpid == 902){
      lineMass = MuonMass(); 
      linePID = MUONBIT;

      //no need to reset muon charge for event, should be the same for the whole sample
      if(globalMuonCharge == -999){
        globalMuonCharge = lineCharge;
      }
      else if(globalMuonCharge != lineCharge){
        cout<<"\n\n ************** Muon charge not consistent! "<<globalMuonCharge<<" "<<lineCharge<<" "<<event<<endl<<endl;
        exit(1);
      }

    }
    else if(tmpid==1){
      //there are in fact quite some antiprotons in DIS (mode=34)
      /*
      if(lineCharge<0){
        cout<<"strange antiproton!!! "<<tmpline<<endl;
      }
      */
      if(lineCharge>0){
        lineMass = ProtonMass(); 
        linePID = PROTONBIT;
        parbit += PROTONBIT;
      }
    }
    else if(tmpid==101){
      lineMass = PionMass();
      linePID = PIONBIT;
      parbit += PIONBIT;
    }
    else if(tmpid==110||tmpid==111){
      lineMass = KaonMass();
      linePID = KAONBIT;
      parbit += KAONBIT;
    }
    else if(tmpid==901){
      lineMass = ElectronMass();
      linePID = ELECTRONBIT;
      parbit += ELECTRONBIT;
    }
  }//all charged
  else{
    if(tmpid == 1){
      lineMass = NeutronMass();
      linePID = NEUTRONBIT;
      parbit += NEUTRONBIT;
    }
    else if(tmpid==101){
      lineMass = PiZeroMass();
      linePID = PIZEROBIT;
      parbit += PIZEROBIT;
    }
    else if(tmpid==999){
      linePID = GAMMABIT;
      parbit += GAMMABIT;
    }
  }

}


bool AnaUtils::GiBUUIsIniN(const double tmppw, const double tmpmom1, const double tmpmom2, const double tmpmom3, const double tmptote)
{
  if(fabs(tmppw)<1E-10){
    //iniNfullp->SetXYZM(tmpmom1, tmpmom2, tmpmom3, lineCharge? ProtonMass() : NeutronMass());
    iniNfullp->SetXYZT(tmpmom1, tmpmom2, tmpmom3, tmptote);//bare mass might not work in the potential

    return true;
  }
  else{
    //only save non-0 perweight so that the whole event after this is not affected!
    perweight = tmppw;
  
    //lineFullMom->SetXYZM(tmpmom1, tmpmom2, tmpmom3, mass);
    lineFullMom->SetXYZT(tmpmom1, tmpmom2, tmpmom3, tmptote);  //tmptote needed for internal calculation of q, W

    return false;  
  }
}

//=======================================================================================================
//Generator independent procedures
//=======================================================================================================

void AnaUtils::MainProceed()
{
  if(expmode==kMINERvA || expmode==kNONE || expmode==kDUNE){
    if(anamode == GFS){
      ProceedMINERvAGFS();
    }
    else if(anamode == LOWRECOIL){
      ProceedMINERvALOWRECOIL();
    }
    else if(anamode == NUBAR1PI){
      ProceedMINERvANUBAR1PI();
    }
    else if(anamode == MMECCQE){
      ProceedMINERvAMECCQE();
    }
    else if(anamode == CLR){
      ProceedCLR();
    }
    else if(anamode == RESPS){
      ProceedRESPS();
    }
    else if(anamode == GFSPIZERO){
      ProceedMINERvAGFSPIZERO();
    }
  }
  else if(expmode==kT2K){
    if(anamode == GFS){
      ProceedT2KGFS();
    }
    else if(anamode == GFSEXP){
      ProceedT2KGFSEXP();
    }
  }
}


void AnaUtils::ProceedT2KGFS()
{
  const double tmpmom = lineFullMom->P();
  const double tmptheta = lineFullMom->Theta()*TMath::RadToDeg();
  //const double tmpEk = lineFullMom->E()-lineFullMom->M();
  const double tmpEk = Ekin(lineFullMom, lineMass);//only use experiment 3 momentum

  if(
     IsMuon()
     && tmpmom> 0.25 && tmptheta< 126.87 //cos > -0.6
     ){//mu- or +
    (*muonfullp)=(*lineFullMom);
    npar+= MUONBIT;
  }
  else if(
          IsProton()
          ){
    if(
       tmpmom>protonfullp->P()
       && tmpmom> 0.45 && tmpmom< 1.0 && tmptheta< 66.42 //cos > 0.4
       ){
      (*protonfullp)=(*lineFullMom);
      npar += PROTONBIT;
      KNsrc = lineKNsource;
    }
  }
  else if(
          IsPion() && (lineCharge+globalMuonCharge) == 0
          //the following is updated to PRC
          && tmpEk > 0.075 && tmpEk < 0.4
          && tmptheta<70
          ){//pion
        
    (*pionfullp)=(*lineFullMom);
    npar += PIONBIT;
  }
  else if(lineIsBkgParticle){//all mesons https://gibuu.hepforge.org/trac/wiki/ParticleIDs                                                                                                                             
    npar += BKGBIT;
  }
}

void AnaUtils::ProceedT2KGFSEXP()
{
  const double tmpmom = lineFullMom->P();
  const double tmptheta = lineFullMom->Theta()*TMath::RadToDeg();

  if(
     IsMuon()
     && tmpmom> 0.25 && tmpmom< 7 && tmptheta< 70
     ){//mu- or +
    (*muonfullp)=(*lineFullMom);
    npar+= MUONBIT;
  }
  else if(
          IsProton()
          ){
    if(
       tmpmom>protonfullp->P()
       && tmpmom> 0.45 && tmpmom< 1.2 && tmptheta< 70
       ){
      (*protonfullp)=(*lineFullMom);
      npar += PROTONBIT;
      KNsrc = lineKNsource;
    }
  }
  else if(
          IsPion() && (lineCharge+globalMuonCharge) == 0
          && tmpmom> 0.15 && tmpmom< 1.2 && tmptheta< 70
          ){//pion
    (*pionfullp)=(*lineFullMom);
    npar += PIONBIT;
  }
  else if(lineIsBkgParticle){//all mesons https://gibuu.hepforge.org/trac/wiki/ParticleIDs
    npar += BKGBIT;
  }
}

void AnaUtils::ProceedMINERvAMECCQE()
{
  const double tmptheta = lineFullMom->Theta()*TMath::RadToDeg();
  const double tmpmupz = lineFullMom->Pz();

  if(
     IsMuon()
     && tmpmupz>1.5 && tmptheta<20
     ){
    (*muonfullp)=(*lineFullMom);
    npar+= MUONBIT;
  }
  else if(
          IsProton()
          ){
    Erecoil += Ekin(lineFullMom, lineMass);
  }
  else if(lineIsMMECCQEBkg){
    npar += BKGBIT;
  }
}

void AnaUtils::ProceedCLR()
{
  const double tmpmom = lineFullMom->P();
  const double tmptheta = lineFullMom->Theta()*TMath::RadToDeg();

  if(
     IsMuon()
     && tmpmom> 1.5 && tmpmom< 10 && tmptheta<20
     ){
    (*muonfullp)=(*lineFullMom);
    npar+= MUONBIT;
  }
  else if(
          IsProton()
          ){
    if(
       tmpmom>protonfullp->P()
       && tmpmom> 0.45 && tmptheta<70
       ){
      (*protonfullp)=(*lineFullMom);
      npar += PROTONBIT;
      KNsrc = lineKNsource;
    }
  }
  else if(linePID == PIONBIT || linePID == PIZEROBIT){//jus tag, no use of momentum for calculation
    npar += PIONBIT;
  }
  else if(lineIsBkgParticle){
    if(globalMuonCharge == -999){
      printf("globalMuonCharge not set before others! %d %d npar %d event %d\n", linePID, globalMuonCharge, npar, event); exit(1);
    }
    npar += BKGBIT;
  }
}

void AnaUtils::ProceedRESPS()
{
  const double tmpmom = lineFullMom->P();
  const double tmptheta = lineFullMom->Theta()*TMath::RadToDeg();

  if(
     IsMuon()
     && tmpmom> 1.5 && tmpmom< 10 && tmptheta<20
     ){
    (*muonfullp)=(*lineFullMom);
    npar+= MUONBIT;
  }
  else if(
          IsProton()
          ){
    if(
       tmpmom>protonfullp->P()
       ){
      (*protonfullp)=(*lineFullMom);
      npar += PROTONBIT;
      KNsrc = lineKNsource;
    }
  }
  else if(linePID == PIONBIT || linePID == PIZEROBIT){//jus tag, no use of momentum for calculation
    (*pionfullp)=(*lineFullMom);
    npar += PIONBIT;
  }
  else if(lineIsBkgParticle){
    if(globalMuonCharge == -999){
      printf("globalMuonCharge not set before others! %d %d npar %d event %d\n", linePID, globalMuonCharge, npar, event); exit(1);
    }
    npar += BKGBIT;
  }
}

void AnaUtils::ProceedMINERvAGFS()
{
  const double tmpmom = lineFullMom->P();
  const double tmptheta = lineFullMom->Theta()*TMath::RadToDeg();
  //const double tmpEk = lineFullMom->E()-lineFullMom->M();
  const double tmpEk = Ekin(lineFullMom, lineMass);//only use experiment 3 momentum

  if(
     //all muons, in fact just one charge, depending on the input
     IsMuon()
     && tmpmom> 1.5 && tmpmom< 10 && tmptheta<20
     ){//mu+ or mu-
    (*muonfullp)=(*lineFullMom);
    npar+= MUONBIT;
  }
  else if(
          IsProton()
          ){
    if(
       tmpmom>protonfullp->P()
       && tmpmom> 0.45 && tmpmom< 1.2 && tmptheta<70
       ){
      (*protonfullp)=(*lineFullMom);
      npar += PROTONBIT;
      KNsrc = lineKNsource;
    }
  }
  else if(
          IsPion() && (lineCharge+globalMuonCharge) == 0
          //&& tmpmom> 0.2 && tmpmom< 4.0 
          && tmpEk > 0.075 && tmpEk < 0.4 
          && tmptheta<70
          ){//pion
    (*pionfullp)=(*lineFullMom);
    npar += PIONBIT;
  }
  else if(
          IsNeutron()
          ){
    (*neutronFSfullp)=(*lineFullMom);
    //npar += NEUTRONBIT;
  }
  else if(lineIsBkgParticle){//all mesons https://gibuu.hepforge.org/trac/wiki/ParticleIDs                                                                                                                                                                           
    if(globalMuonCharge == -999){
      printf("globalMuonCharge not set before others! %d %d npar %d event %d\n", linePID, globalMuonCharge, npar, event); exit(1);
    }
    npar += BKGBIT;
  }
}

void AnaUtils::ProceedMINERvAGFSPIZERO()
{
  const double tmpmom = lineFullMom->P();
  const double tmptheta = lineFullMom->Theta()*TMath::RadToDeg();

  if(
     IsMuon()
     && tmpmom> 1.5 && tmpmom< 20 && tmptheta<25
     ){
    (*muonfullp)=(*lineFullMom);
    npar+= MUONBIT;
  }
  else if(
          IsProton()
          ){
    if(
       tmpmom>protonfullp->P()
       && tmpmom> 0.45
       ){
      (*protonfullp)=(*lineFullMom);
      npar += PROTONBIT;
      KNsrc = lineKNsource;
    }
  }
  else if(
          IsPiZero()
          ){
    if(
       tmpmom>pionfullp->P()
       ){
      (*pionfullp)=(*lineFullMom);
      npar += PIONBIT;
    }
  }
  else if(lineIsBkgParticle){
    npar += BKGBIT;
  }
}

void AnaUtils::ProceedMINERvALOWRECOIL()
{
  const double tmpmom = lineFullMom->P();
  const double tmptheta = lineFullMom->Theta()*TMath::RadToDeg();
  const double tmpE = lineFullMom->E();

  //static void AddEav(const int id, const int charge, const double energy);
  //static int icounter=0;
  AddEav();//(icounter++)<50);

  //add Enu cut in case of GENIE; for GiBUU the cut is granted by flux cut in job card
  if(IsMuon() && tmptheta<20 &&
     enu > 2 && enu < 6
     ){
    double kk = -999;
    if(lineCharge<0){//nu
      kk = tmpE;
    }
    else{//nubar
      kk = tmpmom;
    }
    if(kk > 1.5){
      (*muonfullp)=(*lineFullMom);

      const TLorentzVector dummyNu(0,0,enu, enu);
      const TLorentzVector lvq= dummyNu - (*muonfullp);
      q3 = lvq.P();

      //only save below 0.8 to save space and time
      if(q3<0.8){      
        npar+= MUONBIT;
      }
    }
  }
  //subsample
  else if(
          IsProton()
          ){
    if(
       tmpmom>protonfullp->P()
       && tmpmom> 0.45 && tmpmom< 1.2 && tmptheta<70
       ){
      (*protonfullp)=(*lineFullMom);
      npar += PROTONBIT;
      KNsrc = lineKNsource;
    }
  }
  else if(
          IsPion() && (lineCharge+globalMuonCharge) == 0
          && tmpmom> 0.2 && tmpmom< 4.0 && tmptheta<70
          ){//pion
    (*pionfullp)=(*lineFullMom);
    npar += PIONBIT;
  }
  else if(lineIsBkgParticle){//all mesons https://gibuu.hepforge.org/trac/wiki/ParticleIDs
    npar += BKGBIT;
  }
}

void AnaUtils::ProceedMINERvANUBAR1PI()
{
  /*
    1) antineutrinos in the Tracker, so CH target
    2) charged-current
->    3) exact 1pi- (no other mesons, allow any number of nucleons)
    4) kinematic cuts:
->    a) 1.5 < Ev < 10 GeV
    b) W < 1.8 GeV (W plot without the W cut of course)
->    c) muon theta < 25 degrees
  */

  const double tmptheta = lineFullMom->Theta()*TMath::RadToDeg();
  const double tmpmom = lineFullMom->P();

  //add Enu cut in case of GENIE; for GiBUU the cut is granted by flux cut in job card 
  if(IsMuon() && tmptheta<25 &&
     enu > 1.5 && enu < 10
     ){
    (*muonfullp)=(*lineFullMom);
    npar+= MUONBIT;
  }
  //->save any leading proton also, not affecting selection
  else if(
          IsProton()
          ){
    if(tmpmom>protonfullp->P()){
      (*protonfullp)=(*lineFullMom);
      npar += PROTONBIT;
      KNsrc = lineKNsource;
    }
  }  
  //<-
  else if(
          IsPion() && (lineCharge+globalMuonCharge) == 0
          ){//pion
    (*pionfullp)=(*lineFullMom);
    npar += PIONBIT;
  }
  else if(lineIsBkgParticle){//all mesons https://gibuu.hepforge.org/trac/wiki/ParticleIDs
    if(globalMuonCharge == -999){
      printf("globalMuonCharge not set before others! %d %d\n", linePID, globalMuonCharge); exit(1);
    }
    npar += BKGBIT;
  }
}

void AnaUtils::AddEav(const bool kprint)
{
  //E_av = \sum KE_p + \sum KE_piplus + \sum KE_piminus + sum E_electrons + \sum E_Kplus + \sum E_kminus + sum E_pizero + sum E_gamma 
  double mass = -999;
  const double EPS = 1e-10;  //to ensure positivity

  //const double energy = lineFullMom->E();//experimetal energy deposit -> will have negative Eav
  const double energy = Energy(lineFullMom, lineMass);

  if(IsProton()){
    mass = ProtonMass();
  }
  else if(IsPion()){
    mass = PionMass();
  }
  else if(IsKaon() || IsElectron() || IsPiZero() || IsGamma()){
    mass = EPS;
  }
  else{//ignore all other particles
    return;
  }

  if(mass>0){
    if(energy<mass){
      printf("AnaUtils::AddEav negative kinematic energy! id %d charge %d energy %f mass %f %f Eav %f\n", linePID, lineCharge, energy, mass, energy-mass, Eav); exit(1);
    }

    Eav += (energy - mass) ;

    if(kprint){
      printf("id %d charge %d energy %f mass %f %f Eav %f\n", linePID, lineCharge, energy, mass, energy-mass, Eav);
    }
  }
}

void AnaUtils::DoFill(TTree *tt)
{
  if(IsGood()){
    Calc();
    tt->Fill();
    static int count = 0;
    count++;
    if(count%10000==0){
      cout<<"Fill Count "<<count<<endl;
    }
  }

  //need to do it also before first event
  //do it outside
  //Ini();
}

void AnaUtils::Ini()
{
  npar=0;
  parbit=0;
  iniNfullp->SetXYZT(0,0,0,0);
  RESpifullp->SetXYZT(0,0,0,0);
  RESnucleonfullp->SetXYZT(0,0,0,0);
  muonfullp->SetXYZT(0,0,0,0);
  protonfullp->SetXYZT(0,0,0,0);
  neutronFSfullp->SetXYZT(0,0,0,0);
  pionfullp->SetXYZT(0,0,0,0);
  baryonfullp->SetXYZT(0,0,0,0);

  Eav = 0;
  Erecoil = 0;
  KNsrc = 0;
}


bool AnaUtils::IsMuon()
{
  return (linePID==MUONBIT);
}

bool AnaUtils::IsProton()
{
  return (linePID==PROTONBIT);
}

bool AnaUtils::IsNeutron()
{
  return (linePID==NEUTRONBIT);
}

bool AnaUtils::IsPion()
{
  return (linePID==PIONBIT);
}

bool AnaUtils::IsKaon()
{
  return (linePID==KAONBIT);
}

bool AnaUtils::IsElectron()
{
  return (linePID==ELECTRONBIT);
}

bool AnaUtils::IsPiZero()
{
  return (linePID==PIZEROBIT);
}

bool AnaUtils::IsGamma()
{
  return (linePID==GAMMABIT);
}

bool AnaUtils::IsBKG()
{
  return (npar>=BKGBIT);
}

int AnaUtils::GetNProtons()
{
  //because of adding bit is in the if() of new leading proton, the nprotons is only such new leading protons in FV
  return (npar%(10*PROTONBIT)) /PROTONBIT;
}

int AnaUtils::GetNNeutrons()
{
  return (npar%(10*NEUTRONBIT)) /NEUTRONBIT;
}

int AnaUtils::GetNMuons()
{
  return (npar%(10*MUONBIT)) /MUONBIT;
}

int AnaUtils::GetNPions()
{
  return (npar%(10*PIONBIT)) /PIONBIT;
}

bool AnaUtils::IsGood()
{

  if(anamode==GFS || anamode==GFSEXP || anamode==CLR || anamode==GFSPIZERO ){
    //1mu
    if(GetNMuons()!=1){
      return false;
    }
    
    //Np at least one proton needed for baryon kinematic calculation
    if(GetNProtons()==0){
      return false;
    }
    

    //0 or 1pi
    if(GetNPions()>1){
      return false;
    }

    if(IsBKG()){
      return false;
    }

    return true;

    
    //1mu
    //0pi at least 1 p
    //1pi any p
    //return (npar >= 101 && npar < 120);
    
  }
  else if(anamode==RESPS){
    //1mu
    if(GetNMuons()!=1){
      return false;
    }

    //0 or 1pi
    if(GetNPions()>1){
      return false;
    }

    if(IsBKG()){
      return false;
    }

    //only care about RES
    if(evtMode!=2){
      return false;
    }

    return true;
  }
  else if(anamode==MMECCQE){
    if(GetNMuons()!=1){
      return false;
    }

    if(IsBKG()){
      return false;
    }

    return true;
  }
  else if(anamode==LOWRECOIL){
    return (GetNMuons()==1);
  }
  else if(anamode==NUBAR1PI){
    if(GetNMuons()!=1)
      return false;

    if(GetNPions()!=1)
      return false;

    if(IsBKG()){
      return false;
    }

    return true;
  }
  else{
    printf("unknown anamode %d\n", anamode); exit(1);
  }
}

TTree * AnaUtils::GetTree(const int ana, const experiment exp)
{
  anamode = ana;
  expmode = exp;
  cout<<"TTree * AnaUtils::GetTree anamode "<<anamode<<" expmode "<<expmode<<endl;

  //============
  TTree * tout = new TTree("tree","tree");
  //const Int_t spl = 1;

  tout->Branch("enu",&enu);
  tout->Branch("muonmomentum",&muonmomentum);
  tout->Branch("muontheta",&muontheta);
  tout->Branch("npar",&npar);

  tout->Branch("run",&run);
  tout->Branch("event",&event);
  tout->Branch("targetZ",&targetZ);
  tout->Branch("prod",&prod);
  tout->Branch("evtMode",&evtMode);

  tout->Branch("perweight",&perweight);

  //tout->Branch("muoncostheta",&muoncostheta);

  tout->Branch("Wtrue",&Wtrue);
  tout->Branch("Wrest",&Wrest);
  tout->Branch("xBj",&xBj);
  //tout->Branch("xCCH",&xCCH);
  //tout->Branch("energyCCH",&energyCCH);
  tout->Branch("xrest",&xrest);
  tout->Branch("Q2",&Q2);
  tout->Branch("muonpt",&muonpt);

  tout->Branch("RESmass",&RESmass);
  tout->Branch("adlerPhi",&adlerPhi);
  tout->Branch("lrsign",&lrsign);
  tout->Branch("KNsrc",&KNsrc);
  tout->Branch("w2",&w2);

  if(anamode==GFS||anamode==GFSEXP||anamode==CLR||anamode==RESPS||anamode==GFSPIZERO){

    if(anamode==CLR||anamode==RESPS){
      tout->Branch("pseudoPhi",&pseudoPhi);
      tout->Branch("pseudosign",&pseudosign);
      tout->Branch("wpseudo2",&wpseudo2);
      tout->Branch("cosNuIniNAngle",&cosNuIniNAngle);
      tout->Branch("cosQIniNAngle",&cosQIniNAngle);
    }

    tout->Branch("protonmomentum",&protonmomentum);
    tout->Branch("protontheta",&protontheta);
    tout->Branch("neutronFSmomentum",&neutronFSmomentum);
    tout->Branch("neutronFStheta",&neutronFStheta);
    tout->Branch("pionmomentum",&pionmomentum);
    tout->Branch("pionEk",&pionEk);
    tout->Branch("piontheta",&piontheta);
    tout->Branch("baryonmomentum",&baryonmomentum);
    tout->Branch("baryontheta",&baryontheta);
    tout->Branch("baryonmass",&baryonmass);
    tout->Branch("dpt",&dpt);
    tout->Branch("dphit",&dphit);
    tout->Branch("dalphat",&dalphat);
    tout->Branch("dpTT",&dpTT);
    tout->Branch("protonTT",&protonTT);
    tout->Branch("pionTT",&pionTT);
    tout->Branch("dpTy",&dpTy);
    tout->Branch("neutronmomentum",&neutronmomentum);
    //tout->Branch("calcEnu",&calcEnu);
    
    /*
      tout->Branch("muonfullp","TLorentzVector", &muonfullp, 128000, spl);
      tout->Branch("protonfullp","TLorentzVector", &protonfullp, 128000, spl);
      tout->Branch("pionfullp","TLorentzVector", &pionfullp, 128000, spl);
      tout->Branch("baryonfullp","TLorentzVector", &baryonfullp, 128000, spl);
    */
  }
  else if(anamode==MMECCQE){
    tout->Branch("q2qe",&q2qe);
    tout->Branch("Erecoil",&Erecoil);
    tout->Branch("mupz",&mupz);
  }
  else if(anamode==LOWRECOIL){
    tout->Branch("parbit",&parbit);
    tout->Branch("q3",&q3);
    tout->Branch("Eav",&Eav);
    tout->Branch("baryonmass",&baryonmass);
    tout->Branch("neutronmomentum",&neutronmomentum);
  }
  else if(anamode==NUBAR1PI){
    tout->Branch("pionmomentum",&pionmomentum);
    tout->Branch("piontheta",&piontheta);
    tout->Branch("pionEk",&pionEk);
  }

  return tout;
}


void AnaUtils::Calc()
{
  muoncostheta = TMath::Cos( muonfullp->Theta()  );
  muonmomentum = muonfullp->P();
  muontheta = muonfullp->Theta()*TMath::RadToDeg();
  muonpt = muonfullp->Pt();
  mupz = muonfullp->Pz();
  q2qe = GetTrueCCQEQ2(muonmomentum, muonfullp->Theta());
  protonmomentum = protonfullp->P();
  protontheta = protonfullp->Theta()*TMath::RadToDeg();
  neutronFSmomentum = neutronFSfullp->P();
  neutronFStheta = neutronFSfullp->Theta()*TMath::RadToDeg();
  pionmomentum = pionfullp->P();
  piontheta = pionfullp->Theta()*TMath::RadToDeg();
  //pionEk = pionfullp->E()-PionMass();
  pionEk = Ekin(pionfullp, PionMass()); //only use experimental momentum

  //(*baryonfullp) = (*protonfullp) + (*pionfullp);
  baryonfullp->SetXYZT(protonfullp->X()+pionfullp->X(), protonfullp->Y()+pionfullp->Y(), protonfullp->Z()+pionfullp->Z(), Energy(protonfullp, ProtonMass())+Energy(pionfullp, pionfullp->P()>1E-10? PionMass():0));//need to use experimental momentum only
  baryonmomentum = baryonfullp->P();
  baryontheta = baryonfullp->Theta()*TMath::RadToDeg();

  baryonmass = baryonfullp->M();

  const TLorentzVector dummyNu(0,0,enu, enu);

  const TVector3 ztt = (dummyNu.Vect().Cross(muonfullp->Vect())).Unit();

  protonTT = protonfullp->Vect().Dot(ztt);
  pionTT = pionfullp->Vect().Dot(ztt);

  const double ares = -0.3; //0.2

  lrsign = 0;
  if(protonfullp->E()>1E-10){
    //right is -1, left is +1
    lrsign = protonfullp->Vect().Dot(ztt) > 0 ? -1 : 1;

    //just use the unweighted lrsign (original one)
    pseudosign = lrsign;
  }

  RESmass = -999;
  adlerPhi = -999;
  pseudoPhi = -999;
  w2 = 1;
  wpseudo2 = 1;
  if(RESpifullp->E()>1E-10 && RESnucleonfullp->E()>1E-10){
    const TLorentzVector ressum = (*RESpifullp) + (*RESnucleonfullp);
    RESmass = ressum.M();
    
    //double GetAdlerPhi(TLorentzVector nufull, TLorentzVector muonfull, TLorentzVector pifull, TLorentzVector nucleonfull, TLorentzVector iniNfull)
    adlerPhi = GetOneBoostAdlerPhi(dummyNu, *muonfullp, *RESpifullp, *RESnucleonfullp, *iniNfullp);

    /*
    //->test
    const double tbadler = GetTwoBoostAdlerPhi(dummyNu, *muonfullp, *RESpifullp, *RESnucleonfullp, *iniNfullp);
    if(fabs(adlerPhi-tbadler)>1e-10){
      printf("one-boost and two-boost not consistent!\n");
      exit(1);
    }
    exit(1);
    //<--
    */

    //by weighting the event with lrsign, the adlerPhi weight is not needed
    w2 = GetW2(adlerPhi, ares);
    lrsign *= w2;

    //============= repeat incorrect calculation in CLR paper first version =====>
    pseudoPhi = GetPseudoPhi(dummyNu, *muonfullp, *RESpifullp, *RESnucleonfullp);
    wpseudo2 = GetW2(pseudoPhi, ares);
    pseudosign *= wpseudo2;
  }
  //<=====

  //==================
  const TLorentzVector lvq= dummyNu - (*muonfullp);
  q3 = lvq.P();
  Q2 = -lvq.M2();
  //=========
  TLorentzVector dummyP, dummyW;
  //assum proton target
  dummyP.SetXYZM(0,0,0,ProtonMass());
  dummyW = lvq+dummyP;
  Wrest = dummyW.M();
  xrest = Q2/2/(dummyP.Dot(lvq));

  dummyW = lvq+(*iniNfullp);
  Wtrue = dummyW.M();
  xBj = Q2/2/(iniNfullp->Dot(lvq));

  //===================
  cosNuIniNAngle = ((dummyNu.Vect()).Unit()).Dot((iniNfullp->Vect()).Unit());
  cosQIniNAngle = ((lvq.Vect()).Unit()).Dot((iniNfullp->Vect()).Unit());


  //not useful
  /*
  //standalone calculation
  energyCCH = EnuCCH(muonfullp);
  const TLorentzVector nuCCH(0, 0, energyCCH, energyCCH);
  const TLorentzVector qCCH=nuCCH-(*muonfullp);
  xCCH = -qCCH.M2()/2/ProtonMass()/qCCH.E();
  */

  /*
  //test->
  //with prod 1 and 2, Wtrue and baryonmass are nearly identical
  //check momentum conservation

  //2p2h 4-momentum can't be checked at vertex because not getting the other proton here
  //FSI also makes checking difficult -> so use noFSI for GiBUU; GENIE is fast with FSI

  //checked GiBUU (carbon noFSI) and GENIE carbon (FSI enabled) with GFS RES:
  //./doAna.sh list/noEuCut/FinalEventsList_005_T2K_nu_freeDelta_T0_noFSInoFSInoFSInoFSI_Flux9_noEuCut_1kruns_1000jobs.log testT2K 0
  //./doAna.sh list/BeamNu_EnuMINERvA_TargetCH.list testMINERvA 0  
  if( evtMode==2 && targetZ == 6 && GetNProtons()==1 && GetNPions()==1 ){

  //checked GiBUU and GENIE hydrogen with GFS RES
  //./doAna.sh list/noEuCut/FinalEventsList_005_T2K_nu_freeDelta_T0_noFSI_Flux9_noEuCut_hydrogen_1kruns_1000jobs.log testT2K 0
  //./doAna.sh list/BeamNu_EnuMINERvA_TargetCH.list testMINERvA 0
  //if( evtMode==2 && targetZ == 1 && GetNProtons()==1 && GetNPions()==1 ){
    cout<<"\ntargetZ "<<targetZ<<" evtMode "<<evtMode<<" event "<<event<<endl;
    const TLorentzVector ptot = dummyNu+(*iniNfullp)-(*muonfullp)-(*protonfullp)-(*pionfullp);
    cout<<"Nu "<<dummyNu.M()<<endl;    
    dummyNu.Print();
    cout<<endl;
    cout<<"iniN "<<iniNfullp->M()<<endl;
    iniNfullp->Print();
    cout<<endl;
    cout<<"mu "<<muonfullp->M()<<endl;
    muonfullp->Print(); 
    cout<<endl;
    cout<<"proton "<<protonfullp->M()<<endl;
    protonfullp->Print();
    cout<<endl;
    cout<<"pion "<<pionfullp->M()<<endl;
    pionfullp->Print();
    cout<<endl;
    cout<<"total "<<ptot.M()<<endl;
    ptot.Print();

    if(fabs(ptot.E())<1E-6 && fabs(ptot.P())<1E-6)
    {
      cout<<"4-momentum conserved!"<<endl;
      exit(1);
    }
  }
  //<-test
  */

  //==================

  const TVector3 unitneutrino(0,0,1);

  //from 
  //http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/AnalysisFramework/External/GENIEXSecExtract/src/XSec.cxx?annotate=1.20;cvsroot=mnvsoft;sortby=date
  
  const TVector3 plmuon   = unitneutrino *   muonfullp->Vect().Dot(unitneutrino);
  const TVector3 plbaryon = unitneutrino * baryonfullp->Vect().Dot(unitneutrino);
  
  const TVector3 pTmuon   = muonfullp->Vect()   - plmuon;
  const TVector3 pTbaryon = baryonfullp->Vect() - plbaryon;
  const TVector3 vdPt     = pTmuon + pTbaryon;

  dpt = vdPt.Mag();

  const TVector3 unitqt = -pTmuon.Unit();

  dphit   = TMath::ACos(pTbaryon.Dot(unitqt)/pTbaryon.Mag())*TMath::RadToDeg(); //in Deg 

  //dpt cutoff for hydrogen res dpt is 1E-5 
  if(dpt>1E-5){
    dalphat = TMath::ACos(    vdPt.Dot(unitqt)/vdPt.Mag()    )*TMath::RadToDeg(); //in Deg
  }
  else{//hydrogen
    dalphat=-999;
  }

  dpTT = dpt * sin(dalphat*TMath::DegToRad());
  dpTy = dpt * cos(dalphat*TMath::DegToRad());
  const Double_t dotcross = baryonfullp->Vect().Dot(unitneutrino.Cross( muonfullp->Vect() ));
  if(dotcross<0){
    dpTT *= -1;
  }

  //modified from original codes by J. Sobczyk 14 Nov 2017
  //https://journals.aps.org/prc/abstract/10.1103/PhysRevC.95.065501
  //Eq. 11
  //all in GeV

  //const double Eprim  = muonfullp->E();
  const double Eprim  = Energy(muonfullp, MuonMass());//use experimental momentum
  const double kprimL = plmuon.Mag();

  const double Epprim = baryonfullp->E();//already experimental momentum
  const double pprimL = plbaryon.Mag();

  const double factor = MA() - Eprim - Epprim + kprimL + pprimL;
  const double pT = vdPt.Mag();

  const double pL = -(MAstar()*MAstar() + pT*pT-factor*factor)/2.0/factor;

  neutronmomentum = sqrt(pL*pL + pT*pT);

  /*
  //test-->
  if(npar>=101 && npar<110){//0pi
  //if(npar>=111 && npar<120){//1pi
  static TList * lout= new TList;
  static int testcounter = 0;
  const TString hname=Form("hh%s", GetGiBUUMode().Data());
  const int nx = 150;
  const double xmin = 0; 
  const double xmax = 15;
  const int ny = 800;
  const double ymin = 0;
  const double ymax = 8;
  TH2D * hh = (TH2D*) gDirectory->Get(hname);
  if(!hh){
    hh = new TH2D(hname, hname, nx, xmin, xmax, ny, ymin, ymax); lout->Add(hh);
  }

  printf("test mastar %f neutronmementum %f factor %f pT %f\n", MAstar(), neutronmomentum, factor, pT);

  for(int ii=0; ii<hh->GetNbinsX(); ii++){
    const double imas = hh->GetXaxis()->GetBinCenter(ii);
    const double tmppl = -(imas*imas + pT*pT-factor*factor)/2.0/factor;
    const double tmppn = sqrt(tmppl*tmppl + pT*pT);
    hh->Fill(imas, tmppn, perweight);
  }
  testcounter++;

  if(testcounter==10000){
    TFile * fout=new TFile("testpn.root","recreate");
    lout->Write();
    fout->Save();
    fout->Write();
    fout->Close();
    exit(1);
  }
  }
  //<--test
  */

  calcEnu = kprimL + pprimL - pL;
}


#endif
