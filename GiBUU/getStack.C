void getStack(const TString fin, const TString anaid, const TString tag, const TString var, const TString nuwrovar)
{
    cout<<" fin "<<fin<<" anaid "<<anaid<<" tag "<<tag<<" var "<<var<<" nuwrovar "<<nuwrovar<<endl;
    TFile::Open(fin);
    cout << " File: " << fin << " is opened successfully" << endl;
    
    //Number of histograms
    const int nHist=5;
    
    TList *lout=new TList;

    //Get Histograms for each event mode
    TH1D * hall=(TH1D*) gDirectory->Get(var+"_all");
    TH1D * hqe=(TH1D*) gDirectory->Get(var+"_qe");
    TH1D * hres=(TH1D*) gDirectory->Get(var+"_res");
    TH1D * hdis=(TH1D*) gDirectory->Get(var+"_dis");
    TH1D * h2p2h=(TH1D*) gDirectory->Get(var+"_2p2h");
    TH1D * hother=(TH1D*) gDirectory->Get(var+"_other");
    
    //Check the other
    if(hother){
        if(hother->GetEntries()){
          printf("\n\n\n\n\n%s has entry!\n\n\n\n\n\n", hother->GetName()); exit(1);
        }
        else cout << hother->GetName()<<" is good and has entries "<< hother->GetEntries() << endl;
    }
    
    //List of histograms
    TH1D *hList[]={hall, hqe, hres, hdis, h2p2h};
    TString tit[]={"all", "QE","RES", "DIS", "2p2h"};
    //Define filling colors for each mode
    const Int_t cols[]={kBlack, kRed-3,  kBlue-3, kGreen-3, kOrange-3};
    //Define THStack to collect all modes
    THStack * stk = new THStack(var,tag);
    lout->Add(stk);
    //Define the legend
    TLegend *lg=new TLegend(0.7, 0.7, 0.9,0.9);
    lg->SetName("lg");
    lg->SetFillStyle(0);
    lout->Add(lg);
    
    //Create the TCanvas
    TCanvas *cc=new TCanvas;
    
    //Loop over each mode
    for(int i=0; i<nHist; i++){
        hList[i]->Scale(1./12.);
        //hList[i]->SetMaximum(12.0e-42);
        //cout << "Scaling to per nucleon in CH,  C and H should have been full nucleus and combined !!\n";
        //Add hitograms to lout
        lout->Add(hList[i]);
        //Setting histograms for each mode
        hList[i]->SetLineWidth(3);
        hList[i]->SetLineStyle(kSolid);
        hList[i]->SetLineColor(cols[i]);
        hList[i]->SetFillColor(cols[i]);
        hList[i]->SetFillStyle(0);
        hList[i]->SetTitle(Form("%s %s",tit[i].Data(), hList[i]->GetTitle()));
        hList[i]->SetName(Form("%s%s",var.Data(), tit[i].Data()));
        hList[i]->GetXaxis()->SetTitle(var);
        hList[i]->SetTitle(var);
        //Except the "all" mode, add to the THStack
        stk->Add(hList[i]);
        //Add the legend entries
        lg->AddEntry(hList[i], tit[i], "fl");
        hList[i]->SetMinimum(0);
        hList[i]->SetStats(0);
        hList[i]->Draw("same");
    }
    
    //stk->GetXaxis()->SetTitle(var);
    //stk->SetTitle(var);
    //gStyle->SetTitleX(0.5);
    //stk->Draw("hist C");
    //Set the "all" mode histogram
    /*hall->SetLineWidth(3);
    hall->SetLineColor(kBlack);
    hall->SetFillStyle(0);
    hall->SetMinimum(0);
    hall->Draw("same");*/
    //Cout xsec
    printf("total integrated xsec %s %s %e\n", anaid.Data(), var.Data(), hall->Integral(0,10000,"width"));
    //Draw the legend
    lg->Draw();
    //Save the plots
    cc->Print(Form("outStack/%s/%s_%s.eps", anaid.Data(), var.Data(), tag.Data()));
    cc->Print(Form("outStack/%s/%s_%s.png", anaid.Data(), var.Data(), tag.Data()));

    TFile *fout=new TFile(Form("outStack/%s/%s_%s.root", anaid.Data(), var.Data(), tag.Data()),"recreate");
    ldirectory=gDirectory->mkdir(tag);
    ldirectory->cd();
    lout->Write();
    fout->Save();
    fout->Close();
}
