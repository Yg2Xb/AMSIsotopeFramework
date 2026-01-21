#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <iostream>
#include <vector>

void cpacc() {
    gStyle->SetOptStat(0);
    //gStyle->SetGridcolor(kGray);
    
    TString dir = "/eos/user/z/zixuan/Isotope/Add/";
    TString fileBe = dir + "Be7_rew_frag4_withBkg_bkgest_eff.root";
    TString fileC12 = dir + "C12_rew_frag4_withBkg_bkgest_eff.root";

    std::vector<TString> detectors = {"TOF", "NaF", "AGL"};
    
    // 创建结果直方图，横轴 0.3 - 22 GeV
    // 这里借用 C12 的一个直方图来获取 Binning 结构
    TFile *fTpl = TFile::Open(fileC12);
    TH1D *h_tpl = (TH1D*)fTpl->Get("L1Inner_BKG_H3a_TOF_Z4_Mass7");
    
    TH1D *res_C12 = (TH1D*)h_tpl->Clone("res_C12"); res_C12->Reset();
    TH1D *res_Be7 = (TH1D*)h_tpl->Clone("res_Be7"); res_Be7->Reset();
    res_C12->SetDirectory(0); res_Be7->SetDirectory(0);
    fTpl->Close();

    auto FillRatio = [&](TString fileName, TH1D* h_out) {
        TFile *f = TFile::Open(fileName);
        if (!f || f->IsZombie()) { std::cout << "Open failed: " << fileName << std::endl; return; }

        for (auto det : detectors) {
            TH1D *h_L1 = (TH1D*)f->Get("L1Inner_BKG_H3a_" + det + "_Z4_Mass7");
            TH1D *h_Unb = (TH1D*)f->Get("UnbiasedL1Inner_BKG_H3a_" + det + "_Z4_Mass7");
            
            if (!h_L1 || !h_Unb) continue;

            for (int i = 1; i <= h_out->GetNbinsX(); ++i) {
                double x = h_out->GetBinCenter(i);
                // 能量分段选择逻辑
                bool use = false;
                if (det == "TOF" && x >= 0.3 && x < 1.3) use = true;
                else if (det == "NaF" && x >= 1.3 && x < 3.2) use = true;
                else if (det == "AGL" && x >= 3.2 && x <= 22.0) use = true;

                if (use) {
                    double valL1 = h_L1->GetBinContent(h_L1->FindBin(x));
                    double valUnb = h_Unb->GetBinContent(h_Unb->FindBin(x));
                    if (valUnb > 0) {
                        h_out->SetBinContent(i, valL1 / valUnb);
                        // 简单误差传播: r = a/b -> dr = r * sqrt((da/a)^2 + (db/b)^2)
                        double errL1 = h_L1->GetBinError(h_L1->FindBin(x));
                        double errUnb = h_Unb->GetBinError(h_Unb->FindBin(x));
                        double relErr = std::sqrt(std::pow(errL1/valL1, 2) + std::pow(errUnb/valUnb, 2));
                        h_out->SetBinError(i, (valL1 / valUnb) * relErr);
                    }
                }
            }
        }
        f->Close();
    };

    // 执行计算
    FillRatio(fileC12, res_C12);
    FillRatio(fileBe, res_Be7);

    // 绘图
    TCanvas *c = new TCanvas("c", "Chain Efficiency Ratio Comparison", 800, 600);
    c->SetGrid();
    
    res_C12->SetTitle("Ratio of Acceptance(L1Inner / UnbiasedL1Inner);E_{k}/n [GeV];Ratio");
    res_C12->GetXaxis()->SetRangeUser(0.3, 22.0);
    res_C12->GetYaxis()->SetRangeUser(0.8, 1.2); // 通常这个比例在 0.9 左右
    
    res_C12->SetMarkerStyle(20);
    res_C12->SetMarkerColor(kRed);
    res_C12->SetLineColor(kRed);
    
    res_Be7->SetMarkerStyle(20);
    res_Be7->SetMarkerColor(kBlack);
    res_Be7->SetLineColor(kBlack);

    res_C12->Draw("P");
    res_Be7->Draw("P SAME");

    TLegend *leg = new TLegend(0.2, 0.12, 0.85, 0.25);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);
    leg->AddEntry(res_C12, "C12 -> Be7 Above L1 FragAcc(recEk)", "ep");
    leg->AddEntry(res_Be7, "Be7 Acc(recEk)", "ep");
    leg->Draw();

    c->SaveAs("Chain_Ratio_Comparison.png");
    std::cout << "Done! Result saved in Chain_Ratio_Comparison.png" << std::endl;
}