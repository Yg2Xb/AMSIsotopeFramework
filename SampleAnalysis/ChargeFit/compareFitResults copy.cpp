#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLine.h>
#include <TError.h>
#include <TROOT.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "../Tool.h"

using namespace std;
using namespace AMS_Iso;

// 全局配置
struct Config {
    vector<string> elements = {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
    vector<string> detectors = {"TOF", "NaF", "AGL"};
    vector<string> chains = {"L1Inner", "UnbiasedL1Inner"};
    vector<string> templates = {"L1QTemplate", "L2QTemplate"};
    
    vector<pair<string, string>> paramList = {
        {"LG", "Width"}, {"LG", "MPV"}, {"LG", "Sigma"}, {"LG", "Chi2NDF"},
        {"EGE", "Peak"}, {"EGE", "SigmaL"}, {"EGE", "SigmaR"},
        {"EGE", "AlphaL"}, {"EGE", "AlphaR"}, {"EGE", "Chi2NDF"}
    };
    
    map<string, int> detColors = {{"TOF", kRed}, {"NaF", kBlue}, {"AGL", kGreen+2}};
    map<string, int> chainColors = {{"L1Inner", kBlack}, {"UnbiasedL1Inner", kRed}};
    map<string, int> templateColors = {{"L1QTemplate", kBlack}, {"L2QTemplate", kRed}};
    map<string, int> modelColors = {{"LG", kRed}, {"EGE", kBlue}};
    
    map<string, pair<double, double>> detRanges = {
        {"TOF", {0.35, 1.1}}, {"NaF", {0.9, 4.0}}, {"AGL", {3.0, 20.0}}
    };
};

Config gConfig;
vector<TF1*> gSplineFits;  // 存储spline fits用于保存

// 获取直方图
TH1* getHist(TFile* file, const string& name) {
    cout << "[DEBUG] Getting histogram: " << name << endl;
    auto* h = dynamic_cast<TH1*>(file->Get(name.c_str()));
    if (!h) {
        cout << "[WARN] Histogram not found: " << name << endl;
        return nullptr;
    }
    h->SetDirectory(nullptr);  // 让ROOT管理，但脱离文件
    h->SetStats(0);
    cout << "[DEBUG] Got histogram " << name << " with " << h->GetEntries() << " entries" << endl;
    return h;
}

// 检查直方图是否有效
bool isValidHist(TH1* h) {
    if (!h) return false;
    for (int b = 1; b <= h->GetNbinsX(); b++) {
        if (h->GetBinContent(b) > 0) return true;
    }
    return false;
}

// 获取Y轴范围
pair<double, double> getYRange(const vector<TH1*>& hists) {
    double minY = 1e30, maxY = -1e30;
    for (auto h : hists) {
        if (!isValidHist(h)) continue;
        for (int b = 1; b <= h->GetNbinsX(); b++) {
            double v = h->GetBinContent(b);
            if (v > 0 && isfinite(v)) {
                minY = min(minY, v);
                maxY = max(maxY, v);
            }
        }
    }
    if (minY >= maxY) return {0.1, 1.0};
    double margin = (maxY - minY) * 0.4;
    return {max(1e-9, minY - margin), maxY + margin};
}

// 执行智能Spline拟合
TF1* performSplineFit(TH1* hist, const string& detector, const string& baseName) {
    cout << "[DEBUG] Starting spline fit for: " << baseName << endl;
    if (!hist || !isValidHist(hist)) {
        cout << "[WARN] Invalid histogram for spline fit" << endl;
        return nullptr;
    }
    
    double fitStart = gConfig.detRanges[detector].first;
    double fitEnd = gConfig.detRanges[detector].second;
    
    // 尝试不同分段数
    for (int segments = 1; segments <= 4; segments++) {
        cout << "[DEBUG] Trying " << segments << " segments" << endl;
        vector<double> xpoints;
        for (int i = 0; i <= segments; i++) {
            xpoints.push_back(fitStart + i * (fitEnd - fitStart) / segments);
        }
        
        string fitName = baseName + "_seg" + to_string(segments);
        TF1* fit = nullptr;
        
        try {
            fit = SplineFit(hist, xpoints.data(), xpoints.size(), 0x38, "b2e2",
                          fitName.c_str(), 0.3, 20.0);
        } catch (...) {
            cout << "[ERROR] SplineFit threw exception for " << segments << " segments" << endl;
            continue;
        }
        
        if (!fit) {
            cout << "[WARN] SplineFit returned null for " << segments << " segments" << endl;
            continue;
        }
        
        int ndf = fit->GetNDF();
        if (ndf < 1) {
            cout << "[WARN] NDF < 1 for " << segments << " segments" << endl;
            delete fit;
            continue;
        }
        
        double chi2ndf = fit->GetChisquare() / ndf;
        cout << "[INFO] Fit with " << segments << " segments: chi2/ndf = " << chi2ndf << endl;
        
        if (chi2ndf < 1.8) {
            gSplineFits.push_back(fit);
            return fit;
        }
        
        if (segments == 4) {  // 最后一次尝试，返回最好的
            gSplineFits.push_back(fit);
            return fit;
        }
        delete fit;
    }
    
    cout << "[ERROR] All spline fit attempts failed for: " << baseName << endl;
    return nullptr;
}

// 画图基础函数
void drawPage(TCanvas* c, const vector<TH1*>& hists, const vector<int>& colors,
              const vector<string>& labels, const string& title, const string& yTitle,
              double xMin, double xMax, const string& pdfName, bool drawBoundaries = false,
              const vector<TF1*>& fits = {}) {
    
    cout << "[DEBUG] Drawing page: " << title << endl;
    c->Clear();
    c->SetLogx();
    
    // 检查有效性
    bool hasValid = false;
    for (auto h : hists) {
        if (isValidHist(h)) {
            hasValid = true;
            break;
        }
    }
    if (!hasValid) {
        cout << "[WARN] No valid histograms to draw" << endl;
        return;
    }
    
    // 创建框架
    auto yRange = getYRange(hists);
    cout << "[DEBUG] Y range: " << yRange.first << " to " << yRange.second << endl;
    
    TH1F* frame = c->DrawFrame(xMin, yRange.first, xMax, yRange.second);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
    frame->GetYaxis()->SetTitle(yTitle.c_str());
    
    // 画直方图
    for (size_t i = 0; i < hists.size(); i++) {
        if (!isValidHist(hists[i])) continue;
        hists[i]->SetLineColor(colors[i]);
        hists[i]->SetMarkerColor(colors[i]);
        hists[i]->SetLineWidth(2);
        hists[i]->SetMarkerStyle(20);
        hists[i]->Draw("E SAME");
        cout << "[DEBUG] Drew histogram " << i << endl;
    }
    
    // 画拟合线
    for (size_t i = 0; i < fits.size(); i++) {
        if (fits[i]) {
            fits[i]->SetLineColor(colors[i]);
            fits[i]->SetLineWidth(2);
            fits[i]->Draw("SAME");
            cout << "[DEBUG] Drew fit " << i << endl;
        }
    }
    
    // 画边界线
    if (drawBoundaries) {
        TLine l1(1.1, yRange.first, 1.1, yRange.second);
        l1.SetLineColor(kRed); l1.SetLineStyle(2); l1.SetLineWidth(2); l1.DrawClone();
        TLine l2(5.0, yRange.first, 5.0, yRange.second);
        l2.SetLineColor(kRed); l2.SetLineStyle(2); l2.SetLineWidth(2); l2.DrawClone();
    }
    
    // 标题和标签
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.10, 0.93, title.c_str());
    
    double y = 0.90;
    for (size_t i = 0; i < labels.size(); i++) {
        latex.SetTextColor(colors[i]);
        latex.DrawLatex(0.75, y, labels[i].c_str());
        y -= 0.06;
    }
    
    c->Update();
    c->Print(pdfName.c_str());
    cout << "[DEBUG] Page printed" << endl;
}

// 主要绘图函数
void plotDetectorComparison(TFile* fin, TCanvas* c, const string& pdfName) {
    cout << "\n[INFO] Starting Detector Comparison plots" << endl;
    
    for (const auto& element : gConfig.elements) {
        for (const auto& chain : gConfig.chains) {
            for (const auto& temp : gConfig.templates) {
                for (const auto& param : gConfig.paramList) {
                    vector<TH1*> hists;
                    for (const auto& det : gConfig.detectors) {
                        string name = chain + "_" + element + "_" + det + "_" + temp + "_" + param.first + "_" + param.second;
                        hists.push_back(getHist(fin, name));
                    }
                    
                    vector<int> colors = {kRed, kBlue, kGreen+2};
                    string title = element + " " + chain + " " + temp + " " + param.first + "_" + param.second;
                    drawPage(c, hists, colors, gConfig.detectors, title, param.second, 
                            0.4, 20.0, pdfName, true);
                }
            }
        }
    }
}

void plotTemplateComparison(TFile* fin, TCanvas* c, const string& pdfName) {
    cout << "\n[INFO] Starting Template Comparison plots" << endl;
    
    for (const auto& element : gConfig.elements) {
        for (const auto& chain : gConfig.chains) {
            for (const auto& det : gConfig.detectors) {
                for (const auto& param : gConfig.paramList) {
                    vector<TH1*> hists;
                    vector<TF1*> fits;
                    
                    for (const auto& temp : gConfig.templates) {
                        string name = chain + "_" + element + "_" + det + "_" + temp + "_" + param.first + "_" + param.second;
                        TH1* h = getHist(fin, name);
                        hists.push_back(h);
                        
                        // 如果不是Chi2参数，尝试拟合
                        if (param.second != "Chi2NDF" && isValidHist(h)) {
                            TF1* fit = performSplineFit(h, det, name);
                            fits.push_back(fit);
                        } else {
                            fits.push_back(nullptr);
                        }
                    }
                    
                    vector<int> colors = {kBlack, kRed};
                    vector<string> labels = {"L1Temp", "L2Temp"};
                    string title = element + " " + chain + " " + det + " " + param.first + "_" + param.second;
                    auto range = gConfig.detRanges[det];
                    drawPage(c, hists, colors, labels, title, param.second, 
                            range.first, range.second, pdfName, false, fits);
                }
            }
        }
    }
}

void plotChainComparison(TFile* fin, TCanvas* c, const string& pdfName) {
    cout << "\n[INFO] Starting Chain Comparison plots" << endl;
    
    for (const auto& element : gConfig.elements) {
        for (const auto& det : gConfig.detectors) {
            for (const auto& temp : gConfig.templates) {
                for (const auto& param : gConfig.paramList) {
                    vector<TH1*> hists;
                    for (const auto& chain : gConfig.chains) {
                        string name = chain + "_" + element + "_" + det + "_" + temp + "_" + param.first + "_" + param.second;
                        hists.push_back(getHist(fin, name));
                    }
                    
                    vector<int> colors = {kBlack, kRed};
                    string title = element + " " + det + " " + temp + " " + param.first + "_" + param.second;
                    auto range = gConfig.detRanges[det];
                    drawPage(c, hists, colors, gConfig.chains, title, param.second,
                            range.first, range.second, pdfName);
                }
            }
        }
    }
}

void plotModelChi2Comparison(TFile* fin, TCanvas* c, const string& pdfName) {
    cout << "\n[INFO] Starting Model Chi2 Comparison plots" << endl;
    
    vector<string> models = {"LG", "EGE"};
    vector<string> modelLabels = {"Landau-Gauss", "ExpGausExp"};
    
    for (const auto& element : gConfig.elements) {
        for (const auto& chain : gConfig.chains) {
            for (const auto& det : gConfig.detectors) {
                for (const auto& temp : gConfig.templates) {
                    vector<TH1*> hists;
                    for (const auto& model : models) {
                        string name = chain + "_" + element + "_" + det + "_" + temp + "_" + model + "_Chi2NDF";
                        hists.push_back(getHist(fin, name));
                    }
                    
                    vector<int> colors = {kRed, kBlue};
                    string title = element + " " + chain + " " + det + " " + temp + " Chi2/NDF";
                    auto range = gConfig.detRanges[det];
                    drawPage(c, hists, colors, modelLabels, title, "#chi^{2}/ndf",
                            range.first, range.second, pdfName);
                }
            }
        }
    }
}

void saveSplineFits(const string& outputFile) {
    cout << "\n[INFO] Saving " << gSplineFits.size() << " spline fits to " << outputFile << endl;
    TFile* outFile = TFile::Open(outputFile.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        cout << "[ERROR] Cannot create output file" << endl;
        return;
    }
    
    for (auto fit : gSplineFits) {
        if (fit) fit->Write();
    }
    outFile->Close();
    delete outFile;
}

void compareFitResults(
    const string& rootfile = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFitParams_BeToOxy_0.5.root",
    const string& outdir = "/eos/user/z/zixuan/Isotope/ChargeFit/comparison_plots/"
) {
    cout << "[INFO] Starting compareFitResults" << endl;
    cout << "[INFO] Input file: " << rootfile << endl;
    cout << "[INFO] Output directory: " << outdir << endl;
    
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    
    TFile* fin = TFile::Open(rootfile.c_str());
    if (!fin || fin->IsZombie()) {
        cout << "[ERROR] Cannot open input file: " << rootfile << endl;
        return;
    }
    cout << "[INFO] Input file opened successfully" << endl;
    
    TCanvas* c = new TCanvas("c", "c", 800, 400);
    c->SetGrid();
    
    // Template比较
    string pdfName = outdir + "FitParam_TemplateCompare_0.5.pdf";
    cout << "\n[INFO] Creating PDF: " << pdfName << endl;
    c->Print((pdfName + "[").c_str());
    plotTemplateComparison(fin, c, pdfName);
    c->Print((pdfName + "]").c_str());
    
    // Detector比较
    pdfName = outdir + "FitParam_DetectorCompare_0.5.pdf";
    cout << "\n[INFO] Creating PDF: " << pdfName << endl;
    c->Print((pdfName + "[").c_str());
    plotDetectorComparison(fin, c, pdfName);
    c->Print((pdfName + "]").c_str());
    
    // Chain比较
    pdfName = outdir + "FitParam_ChainCompare_0.5.pdf";
    cout << "\n[INFO] Creating PDF: " << pdfName << endl;
    c->Print((pdfName + "[").c_str());
    plotChainComparison(fin, c, pdfName);
    c->Print((pdfName + "]").c_str());
    
    // Model Chi2比较
    pdfName = outdir + "FitParam_Chi2ModelCompare_0.5.pdf";
    cout << "\n[INFO] Creating PDF: " << pdfName << endl;
    c->Print((pdfName + "[").c_str());
    plotModelChi2Comparison(fin, c, pdfName);
    c->Print((pdfName + "]").c_str());
    
    // 保存spline fits
    string splineOutFile = outdir + "allFitHistSplineSmooth_0.5.root";
    saveSplineFits(splineOutFile);
    
    delete c;
    fin->Close();
    delete fin;
    
    cout << "[INFO] All done!" << endl;
}