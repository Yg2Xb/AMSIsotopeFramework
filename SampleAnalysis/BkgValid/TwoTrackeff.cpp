#include "../Tool.h" // 包含 SplineFit
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLatex.h"

using namespace std;

// ---------------- 1. 配置路径 ----------------
static const string pathEqu2 = "/eos/user/z/zixuan/Isotope/BkgValid/Boron_to_Beryllium_UnbiasedL1Inner_Validation.root";
static const string pathISS  = "/eos/user/z/zixuan/Isotope/Add/EffISS.root";
static const string pathMC   = "/eos/user/z/zixuan/Isotope/Add/EffMC.root";
static const string pathBkg  = "/eos/user/z/zixuan/Isotope/BkgEst/Epsilon_Results.root";

// 输出配置
static const string outRoot = "FullCorrection_Be10_Result.root";
static const string outDir  = "."; 

// Equ3 探测器区间定义
struct DetectorRange {
    string name;       
    double minEk;
    double maxEk;
};
static const vector<DetectorRange> detectors = {
    {"TOF", 0.3, 1.3},
    {"NaF", 1.3, 3.1},
    {"AGL", 3.1, 21.5}
};

// ---------------- 2. 辅助函数 ----------------

// [关键] 将 TH1 转为 TGraphErrors 仅用于 SplineFit 输入
// SplineFit 通常需要 Graph 结构来确定节点和误差
TGraphErrors* histToGraph(TH1* h, double xmin, double xmax) {
    TGraphErrors* g = new TGraphErrors();
    int pt = 0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double x = h->GetBinCenter(i);
        double y = h->GetBinContent(i);
        double ey = h->GetBinError(i);
        if (x < xmin || x > xmax) continue;
        if (y == 0 && ey == 0) continue; // 跳过空点

        g->SetPoint(pt, x, y);
        g->SetPointError(pt, 0, ey);
        pt++;
    }
    return g;
}

// 生成节点
static vector<double> generateXPoints(const TGraphErrors* g, double xmin, double xmax) {
    vector<double> pts;
    if (!g) return pts;
    const int n = g->GetN();
    const double* xs = g->GetX();
    for (int i = 0; i < n; ++i) {
        if (xs[i] < xmin || xs[i] > xmax) continue;
        pts.push_back(xs[i]); 
    }
    return pts;
}

TH1* getHist(TFile* f, string name) {
    TH1* h = (TH1*)f->Get(name.c_str());
    if (!h) cout << "[Warning] Histogram not found: " << name << " in " << f->GetName() << endl;
    return h;
}

// ---------------- 3. 主程序 ----------------
void CorrBkg_Hist() {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);

    TFile* fEqu2 = TFile::Open(pathEqu2.c_str());
    TFile* fISS  = TFile::Open(pathISS.c_str());
    TFile* fMC   = TFile::Open(pathMC.c_str());
    TFile* fBkg  = TFile::Open(pathBkg.c_str());

    if (!fEqu2 || !fISS || !fMC || !fBkg) {
        cerr << "Error: Cannot open input files." << endl;
        return;
    }
    
    TFile* fOut = new TFile(outRoot.c_str(), "RECREATE");

    // 获取一个模版直方图用于定义 Binning (从 Equ2 Combine 中取)
    string nameIso = "Be10";
    TH1* h_template = getHist(fEqu2, "Combine/Equ2/h_iss_ratio_" + nameIso);
    if (!h_template) h_template = getHist(fEqu2, "Combine/Equ2/h_iss_ratio_Combine_" + nameIso);
    if (!h_template) { cerr << "Fatal: Template hist not found." << endl; return; }

    // =========================================================
    // Step 1: Equ3 Efficiency Ratio (MC / ISS)
    // =========================================================
    cout << "--- Step 1: Processing Equ3 Efficiency (Hist Mode) ---" << endl;
    
    // 创建 Equ3 比值直方图
    TH1D* h_ratio_eq3 = (TH1D*)h_template->Clone("h_ratio_eq3_mc_over_iss");
    h_ratio_eq3->Reset(); // 清空内容，保留 Binning
    h_ratio_eq3->SetTitle("Equ3 Efficiency Ratio (MC/ISS)");

    for (const auto& det : detectors) {
        string hNameNum = "UnbiasedL1Inner_FLUXH1_BkgReduction_Num_" + det.name + "_Z6";
        string hNameDen = "UnbiasedL1Inner_FLUXH1_BkgReduction_Den_" + det.name + "_Z6";

        TH1* hNumISS = getHist(fISS, hNameNum);
        TH1* hDenISS = getHist(fISS, hNameDen);
        TH1* hNumMC  = getHist(fMC, hNameNum);
        TH1* hDenMC  = getHist(fMC, hNameDen);

        if (!hNumISS || !hDenISS || !hNumMC || !hDenMC) continue;

        // 遍历目标直方图的 Bin
        for (int i = 1; i <= h_ratio_eq3->GetNbinsX(); ++i) {
            double c = h_ratio_eq3->GetBinCenter(i);
            
            // 严格控制探测器区间拼接
            if (c < det.minEk || c >= det.maxEk) continue;

            // 在源直方图中找到对应的 Bin
            int binISS = hNumISS->FindBin(c);
            int binMC  = hNumMC->FindBin(c);

            double nI = hNumISS->GetBinContent(binISS); double dI = hDenISS->GetBinContent(binISS);
            double nM = hNumMC->GetBinContent(binMC);   double dM = hDenMC->GetBinContent(binMC);

            if (dI > 0 && dM > 0 && nI > 0 && nM > 0) {
                double effI = nI / dI;
                double effM = nM / dM;
                
                // Ratio = MC / ISS (Equ4 修正项)
                double r = effM / effI; 

                // 简单的误差传递 (假设二项分布误差)
                double relErrI = sqrt(1.0/nI - 1.0/dI); 
                double relErrM = sqrt(1.0/nM - 1.0/dM);
                if (std::isnan(relErrI)) relErrI = 0; // 防止 nI=dI 时的数值问题
                if (std::isnan(relErrM)) relErrM = 0;

                double err = r * sqrt(relErrI*relErrI + relErrM*relErrM);

                h_ratio_eq3->SetBinContent(i, r);
                h_ratio_eq3->SetBinError(i, err);
            }
        }
    }

    // 拟合 Step 1 (需要转 Graph 给 SplineFit)
    TGraphErrors* g_temp_eq3 = histToGraph(h_ratio_eq3, 0.3, 21.5);
    vector<double> nodesEq3 = generateXPoints(g_temp_eq3, 0.3, 21.5);
    TF1* fitEq3 = SplineFit((TGraphAsymmErrors*)g_temp_eq3, nodesEq3.data(), nodesEq3.size(), 0x38, "b2e2", "func_Equ3_Fit", 0.3, 21.5);
    if (fitEq3) fitEq3->SetLineColor(kRed);

    // Plot Step 1 

[Image of Efficiency Ratio Plot]

    TCanvas* c1 = new TCanvas("c1_Equ3", "Equ3 Ratio", 800, 600);
    c1->SetLogx();
    h_ratio_eq3->GetYaxis()->SetRangeUser(0.8, 1.2);
    h_ratio_eq3->SetMarkerStyle(20); h_ratio_eq3->SetLineColor(kBlack);
    h_ratio_eq3->Draw("P E");
    if(fitEq3) fitEq3->Draw("SAME");
    c1->SaveAs((outDir + "/Step1_Equ3_Ratio_Fit.png").c_str());

    // =========================================================
    // Step 2: Equ2 Ratio (ISS / MC) from Combine
    // =========================================================
    cout << "--- Step 2: Processing Equ2 (Combine) ---" << endl;
    
    TH1* h_iss_eq2 = getHist(fEqu2, "Combine/Equ2/h_iss_ratio_" + nameIso); 
    TH1* h_mc_eq2  = getHist(fEqu2, "Combine/Equ2/h_mc_ratio_" + nameIso);
    
    if (!h_iss_eq2) h_iss_eq2 = getHist(fEqu2, "Combine/Equ2/h_iss_ratio_Combine_" + nameIso);
    if (!h_mc_eq2)  h_mc_eq2  = getHist(fEqu2, "Combine/Equ2/h_mc_ratio_Combine_" + nameIso);

    if (!h_iss_eq2 || !h_mc_eq2) { cerr << "Fatal: Equ2 hists missing." << endl; return; }

    TH1D* h_ratio_eq2 = (TH1D*)h_iss_eq2->Clone("h_ratio_eq2_iss_over_mc");
    h_ratio_eq2->Divide(h_mc_eq2); // 计算 ISS / MC

    // Plot Step 2
    TCanvas* c2 = new TCanvas("c2_Equ2", "Equ2 Ratio", 800, 600);
    c2->SetLogx();
    h_ratio_eq2->SetTitle("Equ2 Combine Ratio (ISS/MC);E_{k}/n [GeV/n];Ratio");
    h_ratio_eq2->SetMarkerStyle(20); h_ratio_eq2->SetLineColor(kBlue);
    h_ratio_eq2->GetXaxis()->SetRangeUser(0.3, 21.5);
    h_ratio_eq2->GetYaxis()->SetRangeUser(0.5, 1.5);
    h_ratio_eq2->Draw("P E");
    c2->SaveAs((outDir + "/Step2_Equ2_Combine_Ratio.png").c_str());

    // =========================================================
    // Step 3: Combine & Final Fit (Hist Operation)
    // Correction = Equ2_Hist * Equ3_Fit_Eval
    // =========================================================
    cout << "--- Step 3: Calculating Final Correction ---" << endl;

    TH1D* h_final_corr = (TH1D*)h_ratio_eq2->Clone("h_final_corr");
    h_final_corr->Reset();
    h_final_corr->SetTitle("Total Correction Factor");

    for (int i = 1; i <= h_ratio_eq2->GetNbinsX(); ++i) {
        double x = h_ratio_eq2->GetBinCenter(i);
        if (x < 0.3 || x > 21.5) continue;

        double valEq2 = h_ratio_eq2->GetBinContent(i);
        double errEq2 = h_ratio_eq2->GetBinError(i);

        if (valEq2 <= 0) continue;

        double valEq3 = 1.0;
        if (fitEq3) valEq3 = fitEq3->Eval(x); // 取 Step 1 的拟合值
        
        // Final Value = Ratio_Eq2 * Ratio_Eq3
        double finalVal = valEq2 * valEq3;
        
        // Error = Err_Eq2 * Ratio_Eq3 (忽略 Fit 误差)
        double finalErr = errEq2 * valEq3; 

        h_final_corr->SetBinContent(i, finalVal);
        h_final_corr->SetBinError(i, finalErr);
    }

    // 分段拟合 (需要转 Graph)
    TGraphErrors* g_temp_final = histToGraph(h_final_corr, 0.3, 21.5);

    vector<double> nodesLow = generateXPoints(g_temp_final, 0.3, 6.1);
    TF1* fitFinalLow = SplineFit((TGraphAsymmErrors*)g_temp_final, nodesLow.data(), nodesLow.size(), 0x38, "b2e2", "fit_final_low", 0.3, 6.1);
    
    vector<double> nodesHigh = generateXPoints(g_temp_final, 6.1, 21.5);
    TF1* fitFinalHigh = SplineFit((TGraphAsymmErrors*)g_temp_final, nodesHigh.data(), nodesHigh.size(), 0x38, "b2e2", "fit_final_high", 6.1, 21.5);

    if(fitFinalLow) { fitFinalLow->SetLineColor(kRed); fitFinalLow->SetLineWidth(2); }
    if(fitFinalHigh) { fitFinalHigh->SetLineColor(kRed); fitFinalHigh->SetLineWidth(2); }

    // Plot Step 3 

[Image of Final Correction Plot]

    TCanvas* c3 = new TCanvas("c3_Final", "Final Correction", 800, 600);
    c3->SetLogx();
    TH1F* frame3 = new TH1F("f3", "Total Correction Factor;E_{k}/n [GeV/n];Correction", 100, 0.25, 25.0);
    frame3->SetMinimum(0.5); frame3->SetMaximum(1.5);
    frame3->Draw();
    
    h_final_corr->SetMarkerStyle(20); h_final_corr->SetMarkerColor(kBlack); h_final_corr->SetLineColor(kBlack);
    h_final_corr->Draw("P E SAME");
    
    if(fitFinalLow) fitFinalLow->Draw("SAME");
    if(fitFinalHigh) fitFinalHigh->Draw("SAME");
    
    TLegend* leg3 = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg3->AddEntry(h_final_corr, "Corr Points (Hist)", "lp");
    leg3->AddEntry(fitFinalLow, "Spline Fit", "l");
    leg3->Draw();
    c3->SaveAs((outDir + "/Step3_Final_Correction_Fit.png").c_str());

    // =========================================================
    // Step 4: Apply to Background Estimate
    // =========================================================
    cout << "--- Step 4: Applying to Background Est ---" << endl;

    TH1F* h_bkg_orig = (TH1F*)getHist(fBkg, "h_eps_B11_Be10");
    if (!h_bkg_orig) { cerr << "Missing Bkg Hist" << endl; fOut->Close(); return; }

    TH1F* h_bkg_corr = (TH1F*)h_bkg_orig->Clone("h_eps_B11_Be10_Corrected");
    h_bkg_corr->SetTitle("B11 #rightarrow Be10 Corrected");

    for (int i = 1; i <= h_bkg_corr->GetNbinsX(); ++i) {
        double x = h_bkg_corr->GetBinCenter(i);
        double content = h_bkg_corr->GetBinContent(i);
        double error   = h_bkg_corr->GetBinError(i);

        double factor = 1.0;
        if (x >= 0.3 && x < 6.1 && fitFinalLow) factor = fitFinalLow->Eval(x);
        else if (x >= 6.1 && x < 21.5 && fitFinalHigh) factor = fitFinalHigh->Eval(x);
        
        h_bkg_corr->SetBinContent(i, content * factor);
        h_bkg_corr->SetBinError(i, error * factor);
    }

    // Plot Step 4
    TCanvas* c4 = new TCanvas("c4_Bkg", "Background Comparison", 800, 600);
    c4->SetLogx();
    h_bkg_orig->SetLineColor(kBlue); h_bkg_orig->SetMarkerColor(kBlue); h_bkg_orig->SetMarkerStyle(24);
    h_bkg_corr->SetLineColor(kRed);  h_bkg_corr->SetMarkerColor(kRed);  h_bkg_corr->SetMarkerStyle(20);
    
    h_bkg_orig->GetXaxis()->SetRangeUser(0.25, 25.0);
    h_bkg_orig->GetYaxis()->SetTitle("Epsilon");
    
    h_bkg_orig->Draw("P E");
    h_bkg_corr->Draw("P E SAME");

    TLegend* leg4 = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg4->AddEntry(h_bkg_orig, "Original Est", "lp");
    leg4->AddEntry(h_bkg_corr, "Corrected Est", "lp");
    leg4->Draw();
    c4->SaveAs((outDir + "/Step4_Bkg_Comparison.png").c_str());

    // Save to ROOT
    fOut->cd();
    if(g_temp_eq3) g_temp_eq3->Write("g_temp_eq3");
    h_ratio_eq3->Write();
    if(fitEq3) fitEq3->Write();
    
    h_ratio_eq2->Write();
    
    h_final_corr->Write();
    if(fitFinalLow) fitFinalLow->Write();
    if(fitFinalHigh) fitFinalHigh->Write();
    
    h_bkg_corr->Write();

    fOut->Close();
    delete c1; delete c2; delete c3; delete c4;
    
    cout << "Done. All plots and root file saved." << endl;
}