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
#include "TSystem.h"

using namespace std;

// ---------------- 配置部分 ----------------

// 基础路径模板
static const string pathBase = "/eos/user/z/zixuan/Isotope/BkgValid/";
static const string suffixEqu2 = "_to_Beryllium_UnbiasedL1Inner_Validation.root";

// 全局效率文件
static const string pathISS  = "/eos/user/z/zixuan/Isotope/Add/EffISS.root";
static const string pathMC   = "/eos/user/z/zixuan/Isotope/Add/EffMC.root";
static const string pathBkg  = "/eos/user/z/zixuan/Isotope/BkgEst/Epsilon_Results.root";

// 输出目录
static const string outDir   = "/eos/user/z/zixuan/Isotope/BkgCorr"; 

// 定义用于修正计算的重核源 (B, C, N, O)
// 注意：Be9, Be10 不在此列，因为它们不通过 Equ2 修正
static const vector<string> sources = {"B10", "B11", "C12", "N14", "N15", "O16"};

// 定义用于最终绘图的所有源 (包含 Be)
// 顺序对应颜色：Be9, Be10, B10, B11, C12, N14, N15, O16
static const vector<string> all_sources = {"Be7", "Be9", "Be10", "B10", "B11", "C12", "N14", "N15", "O16"};
const vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta, kOrange-3, kCyan+1, kPink+7, kSpring-5};

static const vector<string> targets = {"Be7", "Be9", "Be10"};

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

// 根据同位素名称获取元素全名
string GetElementName(string isotope) {
    if (isotope.find("Be") == 0) return "Beryllium"; // 新增 Be 处理
    if (isotope.find("B") == 0) return "Boron";
    if (isotope.find("C") == 0) return "Carbon";
    if (isotope.find("N") == 0) return "Nitrogen";
    if (isotope.find("O") == 0) return "Oxygen";
    return "Unknown";
}

TH1* getHist(TFile* f, string name) {
    if (!f) return nullptr;
    TH1* h = (TH1*)f->Get(name.c_str());
    return h;
}

// [关键] 严格按照要求设置 Y 轴范围
void SetSmartYRange(TH1* frame, const vector<TH1*>& hists, TGraph* g = nullptr) {
    double globalMax = -1.0e9;
    double globalMinNonZero = 1.0e9;
    bool foundAny = false;

    auto checkVal = [&](double val) {
        if (val > globalMax) globalMax = val;
        if (val != 0 && val < globalMinNonZero) {
            globalMinNonZero = val;
            foundAny = true;
        }
    };

    // Check Hists
    for (auto h : hists) {
        if (!h) continue;
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double c = h->GetBinCenter(i);
            if (c > 0.3 && c < 21.5) {
                checkVal(h->GetBinContent(i));
            }
        }
    }

    // Check Graph
    if (g) {
        for (int i = 0; i < g->GetN(); ++i) {
            double y = g->GetY()[i];
            double x = g->GetX()[i];
            if (x > 0.3 && x < 21.5) checkVal(y);
        }
    }

    if (!foundAny) {
        frame->SetMinimum(0.1);
        frame->SetMaximum(1.0);
        return;
    }

    double diff = globalMax - globalMinNonZero;
    if (diff < 1e-5) diff = globalMax * 0.5; 
    if (diff == 0) diff = 1.0; 

    double finalMin = globalMinNonZero - 0.3 * diff; // 稍微扩大一点留白
    double finalMax = globalMax + 0.3 * diff;

    if (finalMin < 0 && globalMinNonZero >= 0) finalMin = 0;

    frame->SetMinimum(finalMin);
    frame->SetMaximum(finalMax);
}

// 自适应样条拟合
TF1* RunAdaptiveSplineFit(TGraphErrors* g, string name, double xmin, double xmax) {
    TF1* bestFit = nullptr;
    int max_segments = 5; 

    for (int n_seg = 1; n_seg <= max_segments; ++n_seg) {
        vector<double> nodes;
        double step = (xmax - xmin) / (double)n_seg;
        for (int i = 0; i <= n_seg; ++i) nodes.push_back(xmin + i * step);

        string currentName = name + "_seg" + to_string(n_seg);
        TF1* f = SplineFit((TGraphAsymmErrors*)g, nodes.data(), nodes.size(), 0x38, "b2e2", currentName.c_str(), xmin, xmax);
        
        if (!f) continue;
        double chi2 = f->GetChisquare();
        double ndf = f->GetNDF();
        double ratio = (ndf > 0) ? chi2 / ndf : 999.0;
        bestFit = f;

        if (ratio <= 2) { 
            bestFit->SetLineColor(kRed); 
            bestFit->SetLineWidth(2);
            return bestFit;
        }
    }
    if (bestFit) {
        bestFit->SetLineColor(kRed);
        bestFit->SetLineWidth(2);
    }
    return bestFit;
}

// ---------------- 计算模块 ----------------

// Step 1: 计算通用的效率修正 (MC/ISS Ratio)
TF1* CalculateEfficiencyCorrection(TFile* fISS, TFile* fMC, TFile* fOut) {
    cout << "\n--- Step 1: Calculating Efficiency Ratio (MC/ISS) ---" << endl;
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
                double eff_iss = (hDenISS->GetBinContent(i) > 0) ? hNumISS->GetBinContent(i)/hDenISS->GetBinContent(i) : 0;
                double err_iss = (hDenISS->GetBinContent(i) > 0) ? sqrt(eff_iss*(1-eff_iss)/hDenISS->GetBinContent(i)) : 0;
                
                double eff_mc = (hDenMC->GetBinContent(i) > 0) ? hNumMC->GetBinContent(i)/hDenMC->GetBinContent(i) : 0;
                double err_mc = (hDenMC->GetBinContent(i) > 0) ? sqrt(eff_mc*(1-eff_mc)/hDenMC->GetBinContent(i)) : 0;

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

    TF1* fitFunc = RunAdaptiveSplineFit(gRatio, "fitRatioStep1", 0.3, 21.5);
    if(fitFunc) fitFunc->SetLineColor(kBlue);

    // Save Plot
    TCanvas c1("c1_Eff", "Efficiency Ratio", 800, 600); c1.SetLogx(1);
    TH1F frame("frame1", ";E_{k}/n [GeV/n];Ratio (MC / ISS)", 100, 0.2, 30.0);
    SetSmartYRange(&frame, {}, gRatio);
    frame.Draw();
    gRatio->Draw("P E SAME");
    if(fitFunc) fitFunc->Draw("SAME");
    c1.SaveAs((outDir + "/Step1_Efficiency_Ratio.png").c_str());
    
    fOut->cd();
    gRatio->Write();
    if(fitFunc) fitFunc->Write();

    return fitFunc;
}

// ---------------- 主程序 ----------------
void CorrBkg_Multi() {
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gSystem->Exec(("mkdir -p " + outDir).c_str());

    // 打开全局文件
    TFile* fISS  = TFile::Open(pathISS.c_str());
    TFile* fMC   = TFile::Open(pathMC.c_str());
    TFile* fBkg  = TFile::Open(pathBkg.c_str());

    if (!fISS || !fMC || !fBkg) { cerr << "[FATAL] Global Files (ISS/MC/Bkg) missing." << endl; return; }

    TFile* fOut = new TFile((outDir + "/FinalCorrection_Results.root").c_str(), "RECREATE");

    // 1. 获取 Step 1 修正 (通用)
    TF1* step1Func = CalculateEfficiencyCorrection(fISS, fMC, fOut);

    // ----------------------------------------------------
    // 主循环：Target (Be7, Be9, Be10)
    // ----------------------------------------------------
    for (const string& target : targets) {
        cout << "\n================ Processing Target: " << target << " ================" << endl;

        TH1F* h_sum_orig = nullptr;
        TH1F* h_sum_corr = nullptr;

        // 用于存储该 Target 下所有分量的 Corrected Hist，用于最后的堆叠画图
        vector<TH1F*> hists_for_plot; 
        vector<string> legends_for_plot;

        // ----------------------------------------------------
        // 循环所有潜在源 (Be + B,C,N,O)
        // ----------------------------------------------------
        for (size_t srcIdx = 0; srcIdx < all_sources.size(); ++srcIdx) {
            string source = all_sources[srcIdx];
            int color = colors[srcIdx % colors.size()]; // 绑定颜色

            // [跳过逻辑]
            if (source == target) continue; // 源和目标相同
            if (source == "B10" && target == "Be10") {
                cout << "   [Skip] Skipping impossible channel: B10 -> Be10" << endl;
                continue;
            }
            if (source == "Be7") {
                continue;
            }

            // 获取原始背景估计
            string hNameBkg = "h_eps_" + source + "_" + target;
            TH1F* h_orig = (TH1F*)getHist(fBkg, hNameBkg);
            if (!h_orig) {
                // 如果某个物理上允许但文件中没有的道，跳过
                // cout << "   [Info] Histogram " << hNameBkg << " not found. Skipping." << endl;
                continue;
            }

            TH1F* h_corr = (TH1F*)h_orig->Clone(Form("h_eps_%s_%s_Corrected", source.c_str(), target.c_str()));
            
            // 判断源类型：是 Be (不修正) 还是 Heavy (Equ2 修正)
            bool isBeSource = (source.find("Be") == 0);

            TF1* fLow = nullptr; 
            TF1* fHigh = nullptr;

            if (isBeSource) {
                // Be 源：直接使用原始值，Factor = 1.0
                // 不需要做 Equ2 计算
                // h_corr 已经是 clone，内容一样，无需修改
                cout << "   [Info] " << source << " -> " << target << " is a Be-Be channel. No Equ2 correction applied." << endl;
            } 
            else {
                // Heavy 源 (B, C, N, O)：执行 Equ2 流程
                
                string elementName = GetElementName(source);
                string dynPathEqu2 = pathBase + elementName + suffixEqu2;
                TFile* fEqu2 = TFile::Open(dynPathEqu2.c_str());

                if (fEqu2) {
                    string hNameISS = "Combine/Equ2/h_iss_ratio_" + target;
                    string hNameMC  = "Combine/Equ2/h_mc_ratio_" + target;
                    TH1* h_iss = getHist(fEqu2, hNameISS);
                    TH1* h_mc  = getHist(fEqu2, hNameMC);

                    if (h_iss && h_mc) {
                        // 绘制 Equ2 分布对比图
                        TCanvas cDist(Form("cDist_%s_%s", source.c_str(), target.c_str()), "Equ2 Distributions", 800, 600);
                        TH1F frameDist("frameDist", (source + "->" + target + " Equ2 Check;R [GV] or E_{k};Events").c_str(), 100, 0.25, 21.5);
                        SetSmartYRange(&frameDist, {h_iss, h_mc}, nullptr);
                        frameDist.Draw();
                        h_iss->SetLineColor(kBlack); h_iss->SetMarkerColor(kBlack); h_iss->SetMarkerStyle(20);
                        h_mc->SetLineColor(kRed);    h_mc->SetMarkerColor(kRed);    h_mc->SetMarkerStyle(24);
                        h_iss->Draw("P E SAME"); h_mc->Draw("P E SAME");
                        cDist.SaveAs((outDir + "/Equ2_Dist_" + source + "_to_" + target + ".png").c_str());

                        // 计算 Ratio
                        TH1D* h_ratio = (TH1D*)h_iss->Clone(Form("h_ratio_%s_%s", source.c_str(), target.c_str()));
                        h_ratio->Divide(h_mc);

                        TGraphErrors* gCorr = new TGraphErrors();
                        gCorr->SetName(Form("gCorr_%s_%s", source.c_str(), target.c_str()));
                        int pt = 0;
                        for(int i=1; i<=h_ratio->GetNbinsX(); ++i) {
                            double x = h_ratio->GetBinCenter(i);
                            if (x < 0.3 || x > 21.5) continue;
                            double y = h_ratio->GetBinContent(i);
                            double e = h_ratio->GetBinError(i);
                            if (y <= 0) continue;
                            double s1 = (step1Func) ? step1Func->Eval(x) : 1.0;
                            gCorr->SetPoint(pt, x, y * s1);
                            gCorr->SetPointError(pt, 0, e * s1);
                            pt++;
                        }

                        // 拟合
                        fLow  = RunAdaptiveSplineFit(gCorr, "fit_" + source + "_" + target + "_low", 0.3, 6.5);
                        fHigh = RunAdaptiveSplineFit(gCorr, "fit_" + source + "_" + target + "_high", 6.5, 21.5);

                        // 绘制 Factor
                        TCanvas cCorr(Form("cCorr_%s_%s", source.c_str(), target.c_str()), "Correction Factor", 800, 600);
                        TH1F frameCorr("frameCorr", (source + "->" + target + " Correction Factor;E_{k}/n [GeV/n];Factor").c_str(), 100, 0.25, 21.5);
                        SetSmartYRange(&frameCorr, {}, gCorr);
                        frameCorr.Draw();
                        gCorr->SetMarkerStyle(20); gCorr->SetLineColor(kBlack); gCorr->Draw("P SAME");
                        if(fLow) fLow->Draw("SAME"); if(fHigh) fHigh->Draw("SAME");
                        cCorr.SaveAs((outDir + "/Factor_" + source + "_to_" + target + ".png").c_str());
                        
                        fOut->cd(); gCorr->Write();
                    }
                    fEqu2->Close();
                }

                // 应用修正
                for (int i = 1; i <= h_corr->GetNbinsX(); ++i) {
                    double x = h_corr->GetBinCenter(i);
                    double y = h_corr->GetBinContent(i);
                    double factor = 1.0;
                    if (x >= 0.3 && x < 6.5 && fLow) factor = fLow->Eval(x);
                    else if (x >= 6.5 && x < 21.5 && fHigh) factor = fHigh->Eval(x);
                    h_corr->SetBinContent(i, y * factor);
                }

                // 绘制 Est 对比 (Orig vs Corr)
                TCanvas cCh(Form("cCh_%s_%s", source.c_str(), target.c_str()), "Channel Comparison", 800, 600);
                TH1F frameCh("frameCh", (source + " -> " + target + " Est. Bkg;E_{k}/n [GeV/n];Fraction").c_str(), 100, 0.25, 21.5);
                SetSmartYRange(&frameCh, {h_orig, h_corr}, nullptr);
                frameCh.Draw();
                h_orig->SetLineColor(kBlack); h_orig->SetMarkerColor(kBlack); h_orig->SetMarkerStyle(20);
                h_corr->SetLineColor(kRed);   h_corr->SetMarkerColor(kRed);h_corr->SetMarkerStyle(20);
                h_orig->Draw("PZ SAME"); h_corr->Draw("PZ SAME");
                TLegend leg(0.4, 0.75, 0.88, 0.88);
                leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04);
                leg.AddEntry(h_orig, "Estimation before Correction", "lp");
                leg.AddEntry(h_corr, "Estimation after Correction", "lp");
                leg.Draw();
                cCh.SaveAs((outDir + "/Est_Comp_" + source + "_to_" + target + ".png").c_str());
            }

            // 累加到 Total Sum
            if (!h_sum_orig) {
                h_sum_orig = (TH1F*)h_orig->Clone(Form("h_Sum_Orig_%s", target.c_str()));
                h_sum_corr = (TH1F*)h_corr->Clone(Form("h_Sum_Corr_%s", target.c_str()));
                h_sum_orig->Reset(); h_sum_corr->Reset();
                h_sum_orig->SetTitle(("Total Est. Bkg for " + target).c_str());
            }
            h_sum_orig->Add(h_orig);
            h_sum_corr->Add(h_corr);

            // 准备最后的分量图数据
            // 设置颜色和样式
            h_corr->SetLineColor(color);
            h_corr->SetMarkerColor(color);
            h_corr->SetMarkerStyle(20);
            h_corr->SetMarkerSize(0.8);
            
            hists_for_plot.push_back((TH1F*)h_corr->Clone()); // 存入列表
            legends_for_plot.push_back(source + " -> " + target);

            fOut->cd();
            h_corr->Write();

        } // End Source Loop

        // ----------------------------------------------------
        // 绘制 Total Sum 对比图 (Orig vs Corr)
        // ----------------------------------------------------
        if (h_sum_orig && h_sum_corr) {
            TCanvas cSum(Form("cSum_%s", target.c_str()), ("Sum " + target).c_str(), 800, 600);
            TH1F frameSum("frSum", ("Total Background Estimation for Target: " + target + ";E_{k}/n [GeV/n];Background Fraction").c_str(), 100, 0.25, 21.5);
            SetSmartYRange(&frameSum, {h_sum_orig, h_sum_corr}, nullptr);
            frameSum.Draw();

            h_sum_orig->SetLineColor(kBlack); h_sum_orig->SetMarkerColor(kBlack); h_sum_orig->SetMarkerStyle(20);
            h_sum_corr->SetLineColor(kRed);   h_sum_corr->SetMarkerColor(kRed);   h_sum_corr->SetMarkerStyle(20);

            h_sum_orig->Draw("PZ SAME");
            h_sum_corr->Draw("PZ SAME");

            TLegend leg(0.4, 0.75, 0.88, 0.88);
            leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04);
            leg.AddEntry(h_sum_orig, "Estimation before Correction", "lp");
            leg.AddEntry(h_sum_corr, "Estimation after Correction", "lp");
            leg.Draw();

            cSum.SaveAs((outDir + "/Total_Sum_Comparison_" + target + ".png").c_str());
            
            fOut->cd();
            h_sum_orig->Write();
            h_sum_corr->Write();
        }

        // ----------------------------------------------------
        // [新增] 绘制 Breakdown 图 (各分量 + Total Corrected)
        // ----------------------------------------------------
        if (h_sum_corr && !hists_for_plot.empty()) {
            TCanvas cBreak(Form("cBreakdown_%s", target.c_str()), ("Breakdown " + target).c_str(), 1200, 600);
            cBreak.SetGrid();
            
            // 准备 Y 轴范围：包含 Total 和所有分量
            vector<TH1*> all_hists_ptr;
            all_hists_ptr.push_back(h_sum_corr);
            for(auto h : hists_for_plot) all_hists_ptr.push_back(h);

            TH1F frameBreak("frBreak", ("Fragmentation to " + target + ";E_{k}/n [GeV/n];Background Fraction").c_str(), 100, 0.2, 21.5);
            SetSmartYRange(&frameBreak, all_hists_ptr, nullptr);
            frameBreak.SetMaximum(h_sum_corr->GetMaximum()*1.25); 
            frameBreak.Draw();

            // 绘制 Total (黑色)
            h_sum_corr->SetLineColor(kBlack); 
            h_sum_corr->SetMarkerColor(kBlack); 
            h_sum_corr->SetMarkerStyle(20);
            h_sum_corr->Draw("PZ SAME");

            // 绘制分量
            for(auto h : hists_for_plot) {
                h->Draw("PZ SAME");
            }

            // Legend
            TLegend legBreak(0.65, 0.50, 0.88, 0.88);
            legBreak.SetBorderSize(1); 
            legBreak.SetFillStyle(0); 
            legBreak.SetTextSize(0.035);
            
            // 按列表顺序添加
            for(size_t i=0; i<hists_for_plot.size(); ++i) {
                legBreak.AddEntry(hists_for_plot[i], legends_for_plot[i].c_str(), "lp");
            }
            legBreak.AddEntry(h_sum_corr, "Total", "lp");
            legBreak.Draw();

            cBreak.SaveAs((outDir + "/Breakdown_" + target + ".png").c_str());
        }

    } // End Target Loop

    fOut->Close();
    cout << "\n=== All Done! Results saved to " << outDir << " ===" << endl;
}