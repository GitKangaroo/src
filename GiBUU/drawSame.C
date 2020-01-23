void SetRangeUser(TH1D * hh,const TString var);

void drawSame(const TString fin_Ar, const TString fin_C, const TString tag, const TString var) {
    
    cout<<" fin_Ar "<<fin_Ar<<" fin_C "<<fin_C<<" tag "<<tag<<" var "<<var<<endl;
    
    //Define the variable list and their units and Y titiles
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
    TString x_title, y_title;
    for(int x = 0; x < 38; x++){
        if(var == varname[x]){
            x_title = xtit[x];
            y_title = ytit[x];
        }
    }
    cout << "x_title: " << x_title << " y_title: " << y_title << endl;
    
    // Open the root file that was produced by Carbon and Argon target
    TFile *g = new TFile(fin_C);
    TFile *f = new TFile(fin_Ar);
    
    //Number of histograms
    const int nHist=5;
    
    //Get Histograms for each event mode
    TH1D * hall_Ar=(TH1D*) f->Get("GFS0piDUNEGiBUUArnu/"+var+"all");
    TH1D * hqe_Ar=(TH1D*) f->Get("GFS0piDUNEGiBUUArnu/"+var+"QE");
    TH1D * hres_Ar=(TH1D*) f->Get("GFS0piDUNEGiBUUArnu/"+var+"RES");
    TH1D * hdis_Ar=(TH1D*) f->Get("GFS0piDUNEGiBUUArnu/"+var+"DIS");
    TH1D * h2p2h_Ar=(TH1D*) f->Get("GFS0piDUNEGiBUUArnu/"+var+"2p2h");
    
    TH1D * hall_C=(TH1D*) g->Get("GFS0piDUNEGiBUUCnu/"+var+"all");
    TH1D * hqe_C=(TH1D*) g->Get("GFS0piDUNEGiBUUCnu/"+var+"QE");
    TH1D * hres_C=(TH1D*) g->Get("GFS0piDUNEGiBUUCnu/"+var+"RES");
    TH1D * hdis_C=(TH1D*) g->Get("GFS0piDUNEGiBUUCnu/"+var+"DIS");
    TH1D * h2p2h_C=(TH1D*) g->Get("GFS0piDUNEGiBUUCnu/"+var+"2p2h");
    
    TH1D *hList_Ar[]={hall_Ar, hqe_Ar, hres_Ar, hdis_Ar, h2p2h_Ar};
    TH1D *hList_C[]={hall_C, hqe_C, hres_C, hdis_C, h2p2h_C};
    TString tit[]={"all", "QE","RES", "DIS", "2p2h"};
    
    //gStyle->SetPadGridX(kTRUE);
    //gStyle->SetPadGridY(kTRUE);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    
    for(int i=0; i<nHist; i++){
        //Create the TCanvas
        TCanvas *cc = new TCanvas;
        
        int binmax_Ar =  hList_Ar[i]->GetMaximumBin();
        int binmax_C = hList_C[i]->GetMaximumBin();
        
        double x_Ar = hList_Ar[i]->GetBinContent(binmax_Ar);
        double x_C = hList_C[i]->GetBinContent(binmax_C);
        
        //Define the legend
        if(var == "dalphat") TLegend *lg=new TLegend(0.1,0.77,0.23,0.9);
        else TLegend *lg=new TLegend(0.77, 0.77, 0.9, 0.9);
        lg->SetName("lg");
        lg->SetFillStyle(0);
        lg->SetBorderSize(0);
        //Set maximum of histogram
        if(x_Ar < x_C) hList_C[i]->SetMaximum(x_C*1.1);
        else hList_C[i]->SetMaximum(x_Ar*1.1);
        //Set user range
        SetRangeUser(hList_C[i],var);
        SetRangeUser(hList_Ar[i],var);
            
        hList_C[i]->SetLineColor(kRed);
        hList_C[i]->SetLineWidth(1);
        hList_C[i]->SetFillStyle(3005);
        hList_C[i]->SetFillColor( kRed);
        //Set x and y axis title
        hList_C[i]->GetXaxis()->SetTitle(x_title);
        hList_C[i]->GetYaxis()->SetTitle(y_title);
        hList_C[i]->Draw("hist C");
        hList_Ar[i]->SetLineColor(kBlack);
        hList_Ar[i]->SetLineWidth(3);
        hList_Ar[i]->Draw("hist C same");
        lg->AddEntry(hList_C[i], tit[i]+"_C", "f");
        lg->AddEntry(hList_Ar[i], tit[i]+"_Ar", "l");
        lg->Draw();
        /*}
        else{
            //Set user range
            SetRangeUser(hList_C[i],var);
            SetRangeUser(hList_Ar[i],var);
            
            hList_Ar[i]->SetLineColor(kBlack);
            hList_Ar[i]->SetLineWidth(3);
            hList_Ar[i]->GetXaxis()->SetTitle(x_title);
            hList_Ar[i]->GetYaxis()->SetTitle(y_title);
            hList_Ar[i]->Draw("hist C");
            hList_C[i]->SetLineColor(kRed -1);
            hList_C[i]->SetLineWidth(2);
            hList_C[i]->SetFillStyle( 3005);
            hList_C[i]->SetFillColor(kRed);
            hList_C[i]->Draw("hist C same");
            lg->AddEntry(hList_C[i], tit[i]+"_C", "f");
            lg->AddEntry(hList_Ar[i], tit[i]+"_Ar", "l");
            lg->Draw();
        }
         */
        cc->Print(Form("outStackSame/test/%s_%s_%s.png", var.Data(), tag.Data(), tit[i].Data()));
    }
}


void SetRangeUser(TH1D * hh,const TString var){
    double xmin=999, xmax=-999;
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
        xmax = 180;//180;//
      }
      else if(var=="dpt"){
        xmin = 0;
        xmax = 1.2;
      }
      else if(var=="neutronmomentum"){
        xmin = 0;
        xmax = 1.2;
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
        xmin = 0;
        xmax = 2;
      }
      else if(var=="Wtrue"){
        xmin = 0;
        xmax = 2;
      }
      else if(var=="protonTT" || var=="pionTT" || var=="dpTT"){
        xmin = -0.4;
        xmax = 0.4;
      }
    if(xmin<xmax){
      hh->GetXaxis()->SetRangeUser(xmin, xmax);
    }
}
