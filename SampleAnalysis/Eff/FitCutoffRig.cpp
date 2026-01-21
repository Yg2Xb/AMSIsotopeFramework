#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TPaveText.h"

void FitCutoffRig() {
    // ================= 配置路径 =================
    std::string inFile = "/eos/user/z/zixuan/Isotope/Add/Be_frag4_NoBkg_Tune_full.root";
    std::string outDir = "/eos/user/z/zixuan/Isotope/Eff/";
    std::string outName = outDir + "Rigidity_Estimator_Fit_Corrected.png";

    gSystem->mkdir(outDir.c_str(), true);

    // ================= 1. 读取原始文件 =================
    TFile *f = TFile::Open(inFile.c_str());
    if (!f || f->IsZombie()) {
        std::cout << "Error: Cannot open input file: " << inFile << std::endl;
        return;
    }

    TH2F *hRaw = (TH2F*)f->Get("ISS_FLUX_H5_1");
    if (!hRaw) {
        std::cout << "Error: Cannot find histogram ISS_FLUX_H5_0" << std::endl;
        return;
    }
    // 先对原始数据进行 Rebin，减少循环次数
    hRaw->RebinX(5);
    hRaw->RebinY(5);

    TGraph *gMPV = new TGraph();
    int pointIdx = 0;

    for (int i = 1; i <= hRaw->GetNbinsX(); ++i) {
        // 当前 Cutoff 值 (X轴)
        double currentCutoff = hRaw->GetXaxis()->GetBinCenter(i);

        // 投影到 Y 轴 (Inner Rigidity)
        TH1D *hProj = hRaw->ProjectionY(Form("py_%d", i), i, i);

        if (hProj->GetEntries() > 50) {
            // 直接找最大值的 Bin
            int maxBin = hProj->GetMaximumBin();
            double peakVal = hProj->GetBinCenter(maxBin);

            // 简单过滤
            if (peakVal > 0 && peakVal < 100) {
                gMPV->SetPoint(pointIdx, currentCutoff, peakVal);
                pointIdx++;
            }
        }
        delete hProj;
    }

    // ================= 4. 线性拟合 =================
    double fitMin = 5.;
    double fitMax = 21.0;
    TF1 *funcFit = new TF1("linearFit", "pol1", fitMin, fitMax);
    funcFit->SetLineColor(kMagenta);
    funcFit->SetLineWidth(3);

    gMPV->Fit(funcFit, "R");

    double p0 = funcFit->GetParameter(0);
    double p1 = funcFit->GetParameter(1);
    std::cout << ">>> Fit Result: InnerRig = " << p0 << " + " << p1 << " * CutoffRig" << std::endl;

    // ================= 5. 画图与保存 =================
    TCanvas *c1 = new TCanvas("c1", "Rigidity Estimator", 1000, 800);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    // 设置正确的标题
    hRaw->SetTitle("Rigidity Estimator Fit;Max Cutoff Rigidity [GV];Inner Tracker Rigidity [GV]");
    
    // 设置显示范围 (根据之前代码的数值)
    hRaw->GetXaxis()->SetRangeUser(2, 25); // Cutoff Range
    hRaw->GetYaxis()->SetRangeUser(2, 25); // Inner Range
    hRaw->GetZaxis()->SetRangeUser(1, 2000); // Inner Range
    c1->SetLogz();

    hRaw->Draw("COLZ");

    gMPV->SetMarkerStyle(20);
    gMPV->SetMarkerSize(0.6);
    gMPV->SetMarkerColor(kBlack); // 用黑色点更明显一点
    gMPV->Draw("P SAME");

    // 画拟合线
    funcFit->Draw("SAME");

    // 标注结果
    TPaveText *pt = new TPaveText(0.15, 0.75, 0.45, 0.88, "NDC");
    pt->SetFillColor(kWhite);
    pt->SetBorderSize(1);
    pt->AddText("Linear Fit (MPV):");
    pt->AddText(Form("Inner = %.3f + %.3f * Cutoff", p0, p1));
    pt->Draw();

    c1->SaveAs(outName.c_str());
    std::cout << ">>> Plot saved to: " << outName << std::endl;
}