void simple_plot() {
    // 打开文件
    TFile* f1_w1 = TFile::Open("/eos/user/z/zixuan/Isotope/Add/Be10_all_w1.root");
    TFile* f1_all = TFile::Open("/eos/user/z/zixuan/Isotope/Add/Be10_all.root");
    TFile* f2 = TFile::Open("/eos/user/z/zixuan/Isotope/FluxSmooth/tf1_integrals.root");
    
    cout << "=== Files opened ===" << endl;
    
    // w1版本
    cout << "\n=== Processing w1 version ===" << endl;
    TH1* h_data_w1 = (TH1*)f1_w1->Get("UnbiasedL1Inner_MC_BKG_H3a_TOF_Z4_Mass10");
    TH1* h_ratio_w1 = (TH1*)f2->Get("h_ratio_R_inv");
    TH1* h_flux_w1 = (TH1*)f1_w1->Get("MC_FLUX_H3");
    
    cout << "w1 - Data integral before: " << h_data_w1->Integral() << endl;
    cout << "w1 - Ratio integral before: " << h_ratio_w1->Integral() << endl;
    
    double tot_w1 = h_flux_w1->GetBinContent(1);
    cout << "w1 - Tot (bin 1): " << tot_w1 << endl;
    
    h_ratio_w1->Scale(tot_w1);
    cout << "w1 - Ratio integral after scale: " << h_ratio_w1->Integral() << endl;
    
    h_data_w1->Divide(h_ratio_w1);
    cout << "w1 - Data integral after divide: " << h_data_w1->Integral() << endl;
    
    h_data_w1->Scale(3.9 * 3.9 * TMath::Pi());
    cout << "w1 - Data integral after final scale: " << h_data_w1->Integral() << endl;
    
    // all版本
    cout << "\n=== Processing all version ===" << endl;
    TH1* h_data_all = (TH1*)f1_all->Get("UnbiasedL1Inner_MC_BKG_H3a_TOF_Z4_Mass10");
    TH1* h_ratio_all = (TH1*)f2->Get("h_ratio_Be10_R");
    TH1* h_tot_all = (TH1*)f1_all->Get("MC_FLUX_H3");
    
    cout << "all - Data integral before: " << h_data_all->Integral() << endl;
    cout << "all - Ratio integral before: " << h_ratio_all->Integral() << endl;
    
    double tot_all = tot_w1;
    cout << "all - Tot (integral): " << tot_all << endl;
    
    h_ratio_all->Scale(tot_all);
    cout << "all - Ratio integral after scale: " << h_ratio_all->Integral() << endl;
    
    h_data_all->Divide(h_ratio_all);
    cout << "all - Data integral after divide: " << h_data_all->Integral() << endl;
    
    h_data_all->Scale(3.9 * 3.9 * TMath::Pi());
    cout << "all - Data integral after final scale: " << h_data_all->Integral() << endl;
    
    // 检查一些bin内容
    cout << "\n=== Bin content check ===" << endl;
    for (int i = 1; i <= 5; i++) {
        cout << "Bin " << i << " - w1: " << h_data_w1->GetBinContent(i) 
             << ", all: " << h_data_all->GetBinContent(i) << endl;
    }
    
    // 画图
    cout << "\n=== Plotting ===" << endl;
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    h_data_w1->SetLineColor(kRed);
    h_data_w1->SetLineWidth(2);
    h_data_all->SetLineColor(kBlue);
    h_data_all->SetLineWidth(2);
    
    h_data_w1->Draw("hist");
    h_data_all->Draw("hist same");
    c->SetLogx();
    
    TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(h_data_w1, "w1 (R^{-1})", "l");
    leg->AddEntry(h_data_all, "all (Be10_R)", "l");
    
    cout << "Saving plot..." << endl;
    c->SaveAs("/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleAnalysis/BkgValid/10.png");
    
    cout << "Closing files..." << endl;
    f1_w1->Close();
    f1_all->Close();
    f2->Close();
    
    cout << "=== Done ===" << endl;
}