void plotBt() {
    // Open the ROOT file
    TFile *f = TFile::Open("/eos/ams/user/z/zuhao/yanzx/ams_data/amsd69n_IHEPQcutDST/1542996677_10.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open file!" << std::endl;
        return;
    }
    
    // Get the tree
    TTree *tree = (TTree*)f->Get("amstreea");
    if (!tree) {
        std::cerr << "Error: Cannot find tree 'amstreea'!" << std::endl;
        return;
    }
    
    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->SetLogy(); // Optional: set log scale on y-axis if needed
    
    // Create histograms
    TH1F *h1 = new TH1F("h1", "", 300, 0, 3);
    TH1F *h2 = new TH1F("h2", "", 300, 0, 3);
    TH1F *h3 = new TH1F("h3", "", 300, 0, 3);
    
    // Set colors
    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h3->SetLineColor(kBlack);
    
    // Set line width for better visibility
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h3->SetLineWidth(2);
    
    // Draw with different cuts
    tree->Draw("tk_rigidity1[1][2][1]/mcutoffi[1][1]>>h1", "btstat_new==1", "goff");
    tree->Draw("tk_rigidity1[1][2][1]/mcutoffi[1][1]>>h2", "btstat_new==2", "goff");
    tree->Draw("tk_rigidity1[1][2][1]/mcutoffi[1][1]>>h3", "btstat_new==3", "goff");
    
    // Set axis titles
    h1->GetXaxis()->SetTitle("Rigidity/30#circ IGRF Max CutoffRig");
    h1->GetYaxis()->SetTitle("Events");
    h1->SetTitle("");
    
    // Find maximum for y-axis range
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    double max3 = h3->GetMaximum();
    double maxVal = TMath::Max(max1, TMath::Max(max2, max3));
    h1->SetMaximum(maxVal * 1.2);
    
    // Draw histograms
    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    h3->Draw("HIST SAME");
    
    // Add legend
    TLegend *leg = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg->AddEntry(h1, "Primary", "l");
    leg->AddEntry(h2, "Secondary", "l");
    leg->AddEntry(h3, "Trapped", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();
    
    // Save to file
    c1->SaveAs("/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleAnalysis/PureIso/1.png");
    
    // Clean up
    delete leg;
    delete h1;
    delete h2;
    delete h3;
    delete c1;
    f->Close();
    delete f;
    
}