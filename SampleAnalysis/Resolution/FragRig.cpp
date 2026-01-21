#include "TFile.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include <vector>
#include <iostream>

void FragRig() {
    TString inName  = "/eos/user/z/zixuan/Isotope/Add/B11_rew_frag4_withBkg.rsl.root";
    TString outName = "/eos/user/z/zixuan/Isotope/Resolution/B11_Be10checkTrueFragBeta.pdf";

    std::vector<TString> histKeys = {
        "L1Inner_IDH7_TOF",
        "L1Inner_IDH7_NaF",
        "L1Inner_IDH7_AGL",
        "L1Inner_IDH7_Tracker"
    };

    TFile *file = TFile::Open(inName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 700);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird); 
    
    // 调整边距：左侧留给Y轴标题，底部留给X轴
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetBottomMargin(0.15);

    c1->Print(outName + "[");

    for (const auto& key : histKeys) {
        TH2F *h2 = (TH2F*)file->Get(key);
        if (!h2) continue;

        // --- 参数设置 ---
        // AGL从4.0开始，其余从1.92开始
        double minRig = (key.Contains("AGL")) ? 4.0 : 0.7;
        double maxRig = 1.2;
        double theoryVal = -0.136;
        double biasMin = -0.;
        double biasMax = 0.5;

        // ===========================
        // 1. Draw 2D Colz
        // ===========================
        c1->Clear();
        
        // 范围与刻度 (Ndivisions = 508)
        h2->GetXaxis()->SetRangeUser(biasMin, biasMax); 
        h2->GetYaxis()->SetRangeUser(minRig, maxRig);
        h2->GetXaxis()->SetNdivisions(508);
        
        // 标题 (使用 LaTeX 格式)
        h2->SetTitle("Beta Change of B11#rightarrowBe10 Fragments at L2 (MC Truth)");
        h2->GetXaxis()->SetTitle("Relative Beta Change: (beta_{L1} - beta_{L2}) / beta_{L1}");
        h2->GetYaxis()->SetTitle("L1 Beta");
        
        // 字体大小 (0.05) & 偏移
        h2->SetTitleSize(0.07);
        h2->GetXaxis()->SetLabelSize(0.05); h2->GetXaxis()->SetTitleSize(0.05);
        h2->GetYaxis()->SetLabelSize(0.05); h2->GetYaxis()->SetTitleSize(0.05);
        h2->GetYaxis()->SetTitleOffset(1.2); 

        // Log与画图
        c1->SetLogx(0); c1->SetLogy(0); c1->SetLogz(1);
        h2->Draw("COLZ");

        // 红色虚线 (x = -0.136)
        TLine line2D(theoryVal, minRig, theoryVal, maxRig);
        line2D.SetLineStyle(2); line2D.SetLineWidth(3); line2D.SetLineColor(kRed);
        line2D.Draw();

        c1->Print(outName);

        // ===========================
        // 2. Draw ProfileY
        // ===========================
        c1->Clear();
        TProfile *prof = h2->ProfileY(key + "_prof");
        
        // 范围与刻度
        prof->SetTitleSize(0.05);
        prof->GetXaxis()->SetRangeUser(minRig, maxRig);
        prof->SetMinimum(biasMin); 
        prof->GetYaxis()->SetNdivisions(508);

        prof->SetTitle("Beta Change of B11#rightarrowBe10 Fragments at L2 (MC Truth)");
        prof->GetXaxis()->SetTitle("Gen Beta [GV]");
        prof->GetYaxis()->SetTitle("Mean of Relative Beta Change");

        // 字体大小 (0.05) & 偏移
        prof->GetXaxis()->SetLabelSize(0.05); prof->GetXaxis()->SetTitleSize(0.05);
        prof->GetYaxis()->SetLabelSize(0.05); prof->GetYaxis()->SetTitleSize(0.05);
        prof->GetYaxis()->SetTitleOffset(1.3); // 避免遮挡

        // 样式 (黑点)
        prof->SetLineColor(kBlack);
        prof->SetMarkerColor(kBlack);
        prof->SetMarkerStyle(20);
        prof->SetMarkerSize(1.0);

        // Log与画图 (PZ: 只画纵向误差)
        c1->SetLogx(0); c1->SetLogy(0); c1->SetLogz(1);
        prof->Draw("PZ"); 

        // 红色虚线 (y = -0.136)
        TLine lineProf(minRig, theoryVal, maxRig, theoryVal);
        lineProf.SetLineStyle(2); lineProf.SetLineWidth(3); lineProf.SetLineColor(kRed);
        lineProf.Draw();

        c1->Print(outName);
    }

    c1->Print(outName + "]");
    file->Close();
    delete c1;
}