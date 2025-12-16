#include "../Tool.h" // 包含 SplineFit
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"

using namespace std;

// ---------------- 配置部分 ----------------
static const string pathEqu2 = "/eos/user/z/zixuan/Isotope/BkgValid/Boron_to_Beryllium_UnbiasedL1Inner_Validation.root";
static const string pathISS  = "/eos/user/z/zixuan/Isotope/Add/EffISS.root";
static const string pathMC   = "/eos/user/z/zixuan/Isotope/Add/EffMC.root";
static const string pathBkg  = "/eos/user/z/zixuan/Isotope/BkgEst/Epsilon_Results.root";
static const string outDir   = "."; 

struct DetectorRange {
    string name;
    string histSuffix; 
    double minEk;
    double maxEk;
};

static const vector<DetectorRange> detectors = {
    {"TOF", "TOF", 0.3, 1.3},
    {"NaF", "NaF", 1.3, 3.1},
    {"AGL", "AGL", 3.1, 21.5}
};

// ---------------- 辅助函数 ----------------

// [关键] 智能设置 Y 轴范围 (Min*0.7, Max*1.3)
void SetSmartYRange(TH1* frame, TGraph* g = nullptr, TH1* h1 = nullptr, TH1* h2 = nullptr) {
    double ymin = 1.0e9;
    double ymax = -1.0e9;

    auto updateMinMax = [&](double val) {
        if (val == 0) return; // 忽略0值防止 range 错乱
        if (val < ymin) ymin = val;
        if (val > ymax) ymax = val;
    };

    // 扫描 Graph
    if (g) {
        for (int i = 0; i < g->GetN(); ++i) updateMinMax(g->GetY()[i]);
    }
    // 扫描 Hist1
    if (h1) {
        for (int i = 1; i <= h1->GetNbinsX(); ++i) {
            // 只考虑画图范围内的点
            if (h1->GetBinCenter(i) > 0.3 && h1->GetBinCenter(i) < 21.5) 
                updateMinMax(h1->GetBinContent(i));
        }
    }
    // 扫描 Hist2
    if (h2) {
        for (int i = 1; i <= h2->GetNbinsX(); ++i) {
            if (h2->GetBinCenter(i) > 0.3 && h2->GetBinCenter(i) < 21.5)
                updateMinMax(h2->GetBinContent(i));
        }
    }

    if (ymin > ymax) { ymin = 0.5; ymax = 1.5; } // Fallback

    double finalMin = ymin * 0.1;
    double finalMax = ymax * 1.2;
    
    frame->SetMinimum(finalMin);
    frame->SetMaximum(finalMax);
}

// [关键] 自适应样条拟合
// 从少节点开始，不满足 Chi2/NDF <= 2 就增加节点
TF1* RunAdaptiveSplineFit(TGraphErrors* g, string name, double xmin, double xmax) {
    TF1* bestFit = nullptr;
    int max_segments = 12; // 最大分段数，防止过拟合

    cout << "   >>> Adaptive Fit [" << name << "] Range: " << xmin << " - " << xmax << endl;

    for (int n_seg = 1; n_seg <= max_segments; ++n_seg) {
        vector<double> nodes;
        // 生成均匀分布的节点
        double step = (xmax - xmin) / (double)n_seg;
        for (int i = 0; i <= n_seg; ++i) {
            nodes.push_back(xmin + i * step);
        }

        string currentName = name + "_seg" + to_string(n_seg);
        TF1* f = SplineFit((TGraphAsymmErrors*)g, nodes.data(), nodes.size(), 0x38, "b2e2", currentName.c_str(), xmin, xmax);
        
        if (!f) continue;

        double chi2 = f->GetChisquare();
        double ndf = f->GetNDF();
        double ratio = (ndf > 0) ? chi2 / ndf : 999.0;

        bestFit = f; // 暂时保存当前结果

        cout << "      Segments: " << n_seg << " Nodes: " << nodes.size() << " Chi2/NDF: " << ratio << endl;

        if (ratio <= 2.0) {
            cout << "      [Converged] Satisfied Chi2/NDF <= 2.0" << endl;
            bestFit->SetLineColor(kRed); // 最终结果设为红色
            bestFit->SetLineWidth(2);
            return bestFit;
        }
    }
    
    cout << "      [Warning] Reached max segments without Chi2/NDF <= 2.0. Returning last fit." << endl;
    if (bestFit) {
        bestFit->SetLineColor(kRed);
        bestFit->SetLineWidth(2);
    }
    return bestFit;
}

TH1* getHist(TFile* f, string name) {
    if (!f) return nullptr;
    TH1* h = (TH1*)f->Get(name.c_str());
    if (!h) cout << " [ERROR] Histogram NOT FOUND: " << name << endl;
    return h;
}

// ---------------- 主程序 ----------------
void CorrBkg() {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    
    cout << "\n=== STARTING ADAPTIVE CORRECTION ===\n" << endl;

    TFile* fEqu2 = TFile::Open(pathEqu2.c_str());
    TFile* fISS  = TFile::Open(pathISS.c_str());
    TFile* fMC   = TFile::Open(pathMC.c_str());
    TFile* fBkg  = TFile::Open(pathBkg.c_str());

    if (!fEqu2 || !fISS || !fMC || !fBkg) {
        cerr << "[FATAL] Files open failed." << endl; return;
    }

    TFile* fOut = new TFile("FinalCorrection_Result_Adaptive.root", "RECREATE");

    // =========================================================================
    // Step 1: 效率修正
    // =========================================================================
    cout << "\n--- Step 1: Equ3 Efficiency Ratio ---" << endl;

    TGraphErrors* gRatio = new TGraphErrors(); gRatio->SetName("gRatio_MC_ISS");

    int pointIdx = 0;
    for (const auto& det : detectors) {
        string hNameNum = "UnbiasedL1Inner_FLUXH1_BkgReduction_Num_" + det.histSuffix + "_Z6";
        string hNameDen = "UnbiasedL1Inner_FLUXH1_BkgReduction_Den_" + det.histSuffix + "_Z6";

        TH1* hNumISS = getHist(fISS, hNameNum); TH1* hDenISS = getHist(fISS, hNameDen);
        TH1* hNumMC  = getHist(fMC, hNameNum);  TH1* hDenMC  = getHist(fMC, hNameDen);

        if (!hNumISS || !hDenISS || !hNumMC || !hDenMC) continue;

        for (int i = 1; i <= hNumISS->GetNbinsX(); ++i) {
            double center = hNumISS->GetBinCenter(i);
            if (center >= det.minEk && center < det.maxEk) {
                // ISS
                double N_ISS = hNumISS->GetBinContent(i); double D_ISS = hDenISS->GetBinContent(i);
                double eff_iss = (D_ISS > 0) ? N_ISS / D_ISS : 0;
                double err_iss = (D_ISS > 0) ? sqrt(eff_iss * (1.0 - eff_iss) / D_ISS) : 0;
                // MC
                double N_MC = hNumMC->GetBinContent(i); double D_MC = hDenMC->GetBinContent(i);
                double eff_mc = (D_MC > 0) ? N_MC / D_MC : 0;
                double err_mc = (D_MC > 0) ? sqrt(eff_mc * (1.0 - eff_mc) / D_MC) : 0;

                if (eff_iss > 0 && eff_mc > 0) {
                    double r = eff_mc / eff_iss;
                    double relErr = sqrt(pow(err_iss/eff_iss, 2) + pow(err_mc/eff_mc, 2));
                    gRatio->SetPoint(pointIdx, center, r);
                    gRatio->SetPointError(pointIdx, 0, r * relErr);
                    pointIdx++;
                }
            }
        }
    }

    // Adaptive Fit for Step 1
    TF1* fitFuncStep1 = RunAdaptiveSplineFit(gRatio, "fitRatioStep1", 0.3, 21.5);
    if(fitFuncStep1) fitFuncStep1->SetLineColor(kBlue); // 中间步骤用蓝色区分

    TCanvas* c1 = new TCanvas("c1", "Step1 Ratio", 800, 600); c1->SetLogx(1);
    TH1F* frame1 = new TH1F("frame1", ";E_{k}/n [GeV/n];Ratio (MC / ISS)", 100, 0.2, 30.0);
    frame1->GetYaxis()->SetRangeUser(0.99, 1.07);
    frame1->GetXaxis()->SetRangeUser(0.32, 30);
    frame1->Draw();
    gRatio->SetMarkerStyle(20); gRatio->Draw("P E SAME");
    if(fitFuncStep1) fitFuncStep1->Draw("SAME");
    
    // Legend
    TLegend* leg1 = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.05);
    leg1->AddEntry(gRatio, "MC/ISS Data", "lp");
    leg1->AddEntry(fitFuncStep1, "Spline Fit", "l");
    leg1->Draw();
    
    c1->SaveAs((outDir + "/Step1_Efficiency_Ratio.png").c_str());

    // =========================================================================
    // Step 2 & 3: 取 Equ2 Combine & 乘法
    // =========================================================================
    cout << "\n--- Step 2 & 3: Processing Equ2 Combine & Multiplying ---" << endl;

    string isoName = "Be10";
    TH1* h_iss_eq2 = getHist(fEqu2, "Combine/Equ2/h_iss_ratio_" + isoName);
    TH1* h_mc_eq2  = getHist(fEqu2, "Combine/Equ2/h_mc_ratio_" + isoName);
    if(!h_iss_eq2) h_iss_eq2 = getHist(fEqu2, "Combine/Equ2/h_iss_ratio_Combine_" + isoName);
    if(!h_mc_eq2)  h_mc_eq2  = getHist(fEqu2, "Combine/Equ2/h_mc_ratio_Combine_" + isoName);

    if(!h_iss_eq2 || !h_mc_eq2) { cerr << "[FATAL] Equ2 Histograms missing!" << endl; return; }

    TH1D* h_equ2_ratio = (TH1D*)h_iss_eq2->Clone("h_equ2_ratio");
    h_equ2_ratio->Divide(h_mc_eq2); 

    TGraphErrors* gFinalCorr = new TGraphErrors();
    gFinalCorr->SetName("gFinalCorrPoints");
    int ptFinal = 0;

    for(int i=1; i<=h_equ2_ratio->GetNbinsX(); ++i) {
        double x = h_equ2_ratio->GetBinCenter(i);
        if (x < 0.3 || x > 21.5) continue;
        double y_eq2 = h_equ2_ratio->GetBinContent(i);
        double e_eq2 = h_equ2_ratio->GetBinError(i);
        if (y_eq2 <= 0) continue;

        double correction_step1 = 1.0;
        if (fitFuncStep1) correction_step1 = fitFuncStep1->Eval(x);

        double final_val = y_eq2 * correction_step1;
        double final_err = e_eq2 * correction_step1;

        gFinalCorr->SetPoint(ptFinal, x, final_val);
        gFinalCorr->SetPointError(ptFinal, 0, final_err);
        ptFinal++;
    }

    // =========================================================================
    // Step 4: Final Fit (Adaptive)
    // =========================================================================
    cout << "\n--- Step 4: Final Fitting ---" << endl;

    TF1* fitFinalLow = RunAdaptiveSplineFit(gFinalCorr, "fit_final_low", 0.3, 6.5);
    TF1* fitFinalHigh = RunAdaptiveSplineFit(gFinalCorr, "fit_final_high", 6.5, 21.5);

    TCanvas* c2 = new TCanvas("c2", "Final Correction", 800, 600); c2->SetLogx(0);
    TH1F* frame2 = new TH1F("frame2", "Total Correction;E_{k}/n [GeV/n];Correction Factor", 100, 0.25, 21.5);
    SetSmartYRange(frame2, gFinalCorr); // 智能 Y 轴
    frame2->Draw();
    gFinalCorr->SetMarkerStyle(20); gFinalCorr->Draw("P E SAME");
    if(fitFinalLow) fitFinalLow->Draw("SAME");
    if(fitFinalHigh) fitFinalHigh->Draw("SAME");

    TLegend* leg2 = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.05);
    leg2->AddEntry(gFinalCorr, "Final Points", "lp");
    leg2->AddEntry(fitFinalLow, "Adaptive Fit", "l");
    //leg2->Draw();

    c2->SaveAs((outDir + "/Step4_Final_Correction.png").c_str());

    // =========================================================================
    // Step 5: Apply to Est
    // =========================================================================
    cout << "\n--- Step 5: Correcting Background ---" << endl;

    TH1F* h_bkg_orig = (TH1F*)getHist(fBkg, "h_eps_B11_Be10");
    if(!h_bkg_orig) return;

    TH1F* h_bkg_corr = (TH1F*)h_bkg_orig->Clone("h_eps_B11_Be10_Corrected");

    for (int i = 1; i <= h_bkg_corr->GetNbinsX(); ++i) {
        double x = h_bkg_corr->GetBinCenter(i);
        double content = h_bkg_corr->GetBinContent(i);
        
        double factor = 1.0;
        if (x >= 0.3 && x < 6.5 && fitFinalLow) factor = fitFinalLow->Eval(x);
        else if (x >= 6.5 && x < 21.5 && fitFinalHigh) factor = fitFinalHigh->Eval(x);
        
        h_bkg_corr->SetBinContent(i, content * factor);
    }

    TCanvas* c3 = new TCanvas("c3", "Bkg Comparison", 800, 600); c3->SetLogx(0);
    
    // 设置 Frame 和 Range
    TH1F* frame3 = new TH1F("frame3", ";E_{k}/n [GeV/n];Epsilon (Bkg Fraction)", 100, 0.25, 21.5);
    SetSmartYRange(frame3, nullptr, h_bkg_orig, h_bkg_corr); // 智能 Y 轴
    frame3->Draw();

    // 样式设置
    // Original: 黑色
    h_bkg_orig->SetLineColor(kBlack); 
    h_bkg_orig->SetMarkerColor(kBlack); 
    h_bkg_orig->SetMarkerStyle(20); 

    // Corrected: 红色
    h_bkg_corr->SetLineColor(kRed);  
    h_bkg_corr->SetMarkerColor(kRed);  
    h_bkg_corr->SetMarkerStyle(20);

    // Draw PZ
    h_bkg_orig->Draw("PZ SAME");
    h_bkg_corr->Draw("PZ SAME");

    // Legend: 无色透明大字
    TLegend* leg3 = new TLegend(0.55, 0.70, 0.88, 0.88);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.05); // 大字
    leg3->AddEntry(h_bkg_orig, "Original Est", "PZ");
    leg3->AddEntry(h_bkg_corr, "Corrected Est", "PZ");
    leg3->Draw();

    c3->SaveAs((outDir + "/Step5_Bkg_Comparison.png").c_str());

    // Save
    fOut->cd();
    gRatio->Write();
    if(fitFuncStep1) fitFuncStep1->Write();
    gFinalCorr->Write();
    if(fitFinalLow) fitFinalLow->Write();
    if(fitFinalHigh) fitFinalHigh->Write();
    h_bkg_corr->Write();
    fOut->Close();

    cout << "\n=== DONE ===" << endl;
}