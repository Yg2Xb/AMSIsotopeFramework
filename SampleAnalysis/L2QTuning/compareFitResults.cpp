#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLine.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <memory>
#include "../Tool.h"

using namespace std;
using namespace AMS_Iso;

struct GlobalStorage {
    vector<TH1*> allHists;
    vector<TH1F*> allFrames;
    vector<pair<TF1*, string>> splineFits;
    
    void addHist(TH1* h) { if(h) allHists.push_back(h); }
    void addFrame(TH1F* f) { if(f) allFrames.push_back(f); }
    void addSplineFit(TF1* f, const string& name) { 
        if(f) splineFits.push_back(make_pair(f, name)); 
    }
};

GlobalStorage gStorage;

struct PlotConfig {
    vector<string> elements = {"Beryllium", "Boron", "Nitrogen"};
    vector<string> detectors = {"TOF", "NaF", "AGL"};
    vector<string> temps = {"L1Temp", "L2Temp"};
    vector<string> temps_unb = {"L1Temp", "L1Temp_unb", "L2Temp", "L2Temp_unb"};
    vector<string> L1_compare = {"L1Temp", "L1Temp_unb"};
    vector<string> L2_compare = {"L2Temp", "L2Temp_unb"};

    vector<pair<string, string>> paramList = {
        {"LG", "Width"}, {"LG", "MPV"}, {"LG", "Sigma"}, {"LG", "Chi2NDF"},
        {"EGE", "Peak"}, {"EGE", "SigmaL"}, {"EGE", "SigmaR"},
        {"EGE", "AlphaL"}, {"EGE", "AlphaR"}, {"EGE", "Chi2NDF"}
    };
    
    // NEW: Configuration for fixed spline segments
    map<string, int> splineSegments = {
        // Format: "Element_Detector_Template_Model_Parameter" -> segments
        // Beryllium, TOF, L1Temp
        {"Beryllium_TOF_L1Temp_LG_Width", 2}, {"Beryllium_TOF_L1Temp_LG_MPV", 2}, {"Beryllium_TOF_L1Temp_LG_Sigma", 3},
        {"Beryllium_TOF_L1Temp_EGE_Peak", 2}, {"Beryllium_TOF_L1Temp_EGE_SigmaL", 3}, {"Beryllium_TOF_L1Temp_EGE_SigmaR", 2},
        {"Beryllium_TOF_L1Temp_EGE_AlphaL", 2}, {"Beryllium_TOF_L1Temp_EGE_AlphaR", 2},
        // Beryllium, TOF, L2Temp
        {"Beryllium_TOF_L2Temp_LG_Width", 3}, {"Beryllium_TOF_L2Temp_LG_MPV", 2}, {"Beryllium_TOF_L2Temp_LG_Sigma", 3},
        {"Beryllium_TOF_L2Temp_EGE_Peak", 2}, {"Beryllium_TOF_L2Temp_EGE_SigmaL", 2}, {"Beryllium_TOF_L2Temp_EGE_SigmaR", 2},
        {"Beryllium_TOF_L2Temp_EGE_AlphaL", 2}, {"Beryllium_TOF_L2Temp_EGE_AlphaR", 2},
        // Beryllium, NaF, L1Temp
        {"Beryllium_NaF_L1Temp_LG_Width", 1}, {"Beryllium_NaF_L1Temp_LG_MPV", 1}, {"Beryllium_NaF_L1Temp_LG_Sigma", 1},
        {"Beryllium_NaF_L1Temp_EGE_Peak", 1}, {"Beryllium_NaF_L1Temp_EGE_SigmaL", 1}, {"Beryllium_NaF_L1Temp_EGE_SigmaR", 1},
        {"Beryllium_NaF_L1Temp_EGE_AlphaL", 1}, {"Beryllium_NaF_L1Temp_EGE_AlphaR", 1},
        // Beryllium, NaF, L2Temp
        {"Beryllium_NaF_L2Temp_LG_Width", 1}, {"Beryllium_NaF_L2Temp_LG_MPV", 1}, {"Beryllium_NaF_L2Temp_LG_Sigma", 1},
        {"Beryllium_NaF_L2Temp_EGE_Peak", 1}, {"Beryllium_NaF_L2Temp_EGE_SigmaL", 1}, {"Beryllium_NaF_L2Temp_EGE_SigmaR", 1},
        {"Beryllium_NaF_L2Temp_EGE_AlphaL", 1}, {"Beryllium_NaF_L2Temp_EGE_AlphaR", 1},
        // Beryllium, AGL, L1Temp
        {"Beryllium_AGL_L1Temp_LG_Width", 1}, {"Beryllium_AGL_L1Temp_LG_MPV", 1}, {"Beryllium_AGL_L1Temp_LG_Sigma", 1},
        {"Beryllium_AGL_L1Temp_EGE_Peak", 1}, {"Beryllium_AGL_L1Temp_EGE_SigmaL", 1}, {"Beryllium_AGL_L1Temp_EGE_SigmaR", 2},
        {"Beryllium_AGL_L1Temp_EGE_AlphaL", 2}, {"Beryllium_AGL_L1Temp_EGE_AlphaR", 2},
        // Beryllium, AGL, L2Temp
        {"Beryllium_AGL_L2Temp_LG_Width", 1}, {"Beryllium_AGL_L2Temp_LG_MPV", 2}, {"Beryllium_AGL_L2Temp_LG_Sigma", 1},
        {"Beryllium_AGL_L2Temp_EGE_Peak", 1}, {"Beryllium_AGL_L2Temp_EGE_SigmaL", 1}, {"Beryllium_AGL_L2Temp_EGE_SigmaR", 1},
        {"Beryllium_AGL_L2Temp_EGE_AlphaL", 1}, {"Beryllium_AGL_L2Temp_EGE_AlphaR", 1},
        // Boron, TOF, L1Temp
        {"Boron_TOF_L1Temp_LG_Width", 3}, {"Boron_TOF_L1Temp_LG_MPV", 3}, {"Boron_TOF_L1Temp_LG_Sigma", 2},
        {"Boron_TOF_L1Temp_EGE_Peak", 3}, {"Boron_TOF_L1Temp_EGE_SigmaL", 3}, {"Boron_TOF_L1Temp_EGE_SigmaR", 3},
        {"Boron_TOF_L1Temp_EGE_AlphaL", 2}, {"Boron_TOF_L1Temp_EGE_AlphaR", 2},
        // Boron, TOF, L2Temp
        {"Boron_TOF_L2Temp_LG_Width", 2}, {"Boron_TOF_L2Temp_LG_MPV", 2}, {"Boron_TOF_L2Temp_LG_Sigma", 3},
        {"Boron_TOF_L2Temp_EGE_Peak", 3}, {"Boron_TOF_L2Temp_EGE_SigmaL", 3}, {"Boron_TOF_L2Temp_EGE_SigmaR", 3},
        {"Boron_TOF_L2Temp_EGE_AlphaL", 3}, {"Boron_TOF_L2Temp_EGE_AlphaR", 3},
        // Boron, NaF, L1Temp
        {"Boron_NaF_L1Temp_LG_Width", 1}, {"Boron_NaF_L1Temp_LG_MPV", 1}, {"Boron_NaF_L1Temp_LG_Sigma", 1},
        {"Boron_NaF_L1Temp_EGE_Peak", 1}, {"Boron_NaF_L1Temp_EGE_SigmaL", 1}, {"Boron_NaF_L1Temp_EGE_SigmaR", 1},
        {"Boron_NaF_L1Temp_EGE_AlphaL", 1}, {"Boron_NaF_L1Temp_EGE_AlphaR", 1},
        // Boron, NaF, L2Temp
        {"Boron_NaF_L2Temp_LG_Width", 1}, {"Boron_NaF_L2Temp_LG_MPV", 2}, {"Boron_NaF_L2Temp_LG_Sigma", 1},
        {"Boron_NaF_L2Temp_EGE_Peak", 1}, {"Boron_NaF_L2Temp_EGE_SigmaL", 1}, {"Boron_NaF_L2Temp_EGE_SigmaR", 1},
        {"Boron_NaF_L2Temp_EGE_AlphaL", 1}, {"Boron_NaF_L2Temp_EGE_AlphaR", 1},
        // Boron, AGL, L1Temp
        {"Boron_AGL_L1Temp_LG_Width", 1}, {"Boron_AGL_L1Temp_LG_MPV", 1}, {"Boron_AGL_L1Temp_LG_Sigma", 1},
        {"Boron_AGL_L1Temp_EGE_Peak", 2}, {"Boron_AGL_L1Temp_EGE_SigmaL", 2}, {"Boron_AGL_L1Temp_EGE_SigmaR", 1},
        {"Boron_AGL_L1Temp_EGE_AlphaL", 1}, {"Boron_AGL_L1Temp_EGE_AlphaR", 1},
        // Boron, AGL, L2Temp
        {"Boron_AGL_L2Temp_LG_Width", 1}, {"Boron_AGL_L2Temp_LG_MPV", 1}, {"Boron_AGL_L2Temp_LG_Sigma", 2},
        {"Boron_AGL_L2Temp_EGE_Peak", 1}, {"Boron_AGL_L2Temp_EGE_SigmaL", 1}, {"Boron_AGL_L2Temp_EGE_SigmaR", 1},
        {"Boron_AGL_L2Temp_EGE_AlphaL", 1}, {"Boron_AGL_L2Temp_EGE_AlphaR", 1}
    };

    map<string, int> detColors = {{"TOF", kRed}, {"NaF", kBlue}, {"AGL", kGreen+2}};
    map<string, int> tempColors = {{"L1Temp", kBlack}, {"L2Temp", kRed}, 
                                  {"L1Temp_unb", kBlue}, {"L2Temp_unb", kGreen+2}};
    map<string, int> modelColors = {{"LG", kRed}, {"EGE", kBlue}};
    map<string, pair<double, double>> detRanges = {
        {"TOF", {0.42, 1.17}}, {"NaF", {1.17, 4.00}}, {"AGL", {3.23, 16.3}}
    };

    pair<double, double> getXRange(const string& det, bool isMulti = false) const {
        if (isMulti) return {0.42, 16.3};
        if (det == "TOF") return {0.42, 1.55};
        if (det == "NaF") return {0.86, 4.91};
        if (det == "AGL") return {2.88, 16.3};
        return {0.2, 16.3};
    }
};

class HistManager {
public:
    static vector<TH1*> getHists(TFile* file, const vector<string>& names) {
        vector<TH1*> result;
        for (const auto& name : names) {
            if (auto h = dynamic_cast<TH1*>(file->Get(name.c_str()))) {
                auto clone = static_cast<TH1*>(h->Clone((name + "_clone").c_str()));
                clone->SetDirectory(0);
                result.push_back(clone);
                gStorage.addHist(clone);
            }
        }
        return result;
    }

    static vector<TH1*> getHistsMultiFile(const vector<TFile*>& files, const string& name) {
        vector<TH1*> result;
        for (size_t i = 0; i < files.size(); i++) {
            if (files[i] && !files[i]->IsZombie()) {
                if (auto h = dynamic_cast<TH1*>(files[i]->Get(name.c_str()))) {
                    auto clone = static_cast<TH1*>(h->Clone((name + "_clone_" + to_string(i)).c_str()));
                    clone->SetDirectory(0);
                    result.push_back(clone);
                    gStorage.addHist(clone);
                }
            }
        }
        return result;
    }

    static pair<double, double> getYRange(const vector<TH1*>& hists) {
        double minY = 1e9, maxY = -1e9;
        for (auto h : hists) {
            if (!h) continue;
            for (int bin = 1; bin <= h->GetNbinsX(); bin++) {
                double val = h->GetBinContent(bin);
                if (val > 0) {
                    minY = min(minY, val);
                    maxY = max(maxY, val);
                }
            }
        }
        if (minY > maxY || minY > 1e8) return {0, 1};
        double margin = (maxY - minY) * 0.4;
        return {max(0.0, minY - margin), maxY + margin};
    }
};

class FittingUtils {
public:
    // Generate xpoints for spline fitting based on segments
    static vector<double> generateXPoints(const string& detector, int segments) {
        PlotConfig config;
        auto range = config.detRanges.at(detector);
        vector<double> xpoints;
        
        for (int i = 0; i <= segments; i++) {
            double x = range.first + i * (range.second - range.first) / segments;
            xpoints.push_back(x);
        }
        return xpoints;
    }

    // MODIFIED: Perform a single spline fit with a fixed number of segments
    static TF1* performSplineFit(TH1* hist, const string& detector, 
                               const string& baseName, const string& temp, int segments) {
        if (!hist || segments <= 0) return nullptr;
        
        auto xpoints = generateXPoints(detector, segments);
        string fitName = baseName;// + "_spline" + to_string(segments);
        
        try {
            auto fit = SplineFit(hist, xpoints.data(), xpoints.size(), 0x38, "b2e1",
                               fitName.c_str(), xpoints[0], xpoints.back());
            
            fit->SetLineColor(temp.find("L1") != string::npos ? kBlack : kRed); // L1 black, L2 red
            fit->SetLineWidth(2);
            fit->SetLineStyle(1);
            gStorage.addSplineFit(fit, fitName);
            return fit;
        } catch (...) {
            // Return nullptr if spline fit fails
            return nullptr;
        }
    }

    // Combined spline fit for detectorCompare
    static TF1* createCombinedSplineFit(const vector<TH1*>& hists, const vector<string>& detectors,
                                       const string& paramName, const string& element, const string& temp) {
        if (paramName == "Chi2NDF") return nullptr;
        
        TH1* tofHist = nullptr;
        for (size_t i = 0; i < hists.size(); ++i) {
            if (detectors[i] == "TOF" && hists[i]) {
                tofHist = hists[i];
                break;
            }
        }
        if (!tofHist) return nullptr;

        string combinedName = "combined_" + to_string(gStorage.allHists.size());
        TH1D* combinedHist = static_cast<TH1D*>(tofHist->Clone(combinedName.c_str()));
        combinedHist->SetDirectory(0);
        combinedHist->Reset();
        gStorage.addHist(combinedHist);

        // Fill combined histogram
        map<string, pair<double, double>> ranges = {
            {"TOF", {0.4, 1.17}}, {"NaF", {1.17, 3.23}}, {"AGL", {3.23, 16.3}}
        };

        for (size_t i = 0; i < hists.size(); ++i) {
            if (!hists[i] || ranges.find(detectors[i]) == ranges.end()) continue;
            auto [rangeMin, rangeMax] = ranges[detectors[i]];
            
            for (int bin = 1; bin <= hists[i]->GetNbinsX(); bin++) {
                double x = hists[i]->GetBinCenter(bin);
                if (x >= rangeMin && x < rangeMax) {
                    int targetBin = combinedHist->FindBin(x);
                    if (targetBin >= 1 && targetBin <= combinedHist->GetNbinsX()) {
                        combinedHist->SetBinContent(targetBin, hists[i]->GetBinContent(bin));
                        combinedHist->SetBinError(targetBin, max(hists[i]->GetBinError(bin), 
                                                               hists[i]->GetBinContent(bin) * 0.01));
                    }
                }
            }
        }

        // Get custom x-points
        vector<double> xpoints = getCustomXPoints(element, temp, paramName);
        try {
            auto fit = SplineFit(combinedHist, xpoints.data(), xpoints.size(), 0x38, "b2e1",
                               "combinedSpline", 0.4, 16.3);
            if (fit) {
                fit->SetLineColor(kBlack);
                fit->SetLineStyle(2);
                fit->SetLineWidth(2);
            }
            return fit;
        } catch (...) {
            return nullptr;
        }
    }

private:
    static vector<double> getCustomXPoints(const string& element, const string& temp, const string& paramName) {
        if (element == "Beryllium" && temp == "L1Temp" && (paramName == "Peak" || paramName == "AlphaL"))
            return {0.4, 0.75, 1.06, 1.4, 1.9, 6, 16.3};
        if (element == "Beryllium" && temp == "L2Temp" && (paramName == "SigmaR" || paramName == "AlphaL"))
            return {0.4, 0.7, 1., 2, 4, 16.3};
        if (element == "Beryllium" && temp == "L2Temp" && paramName == "AlphaR")
            return {0.4, 0.7, 1., 2, 4, 16.3};
        if (element == "Boron" && temp == "L1Temp" && (paramName == "MPV" || paramName == "Sigma"))
            return {0.4, 0.9, 1.5, 3, 7, 16.3};
        if (element == "Boron" && temp == "L1Temp" && 
            (paramName == "Peak" || paramName == "SigmaL" || paramName == "SigmaR" || paramName == "AlphaR"))
            return {0.4, 0.9, 1.4, 2.2, 6, 16.3};
        if (element == "Boron" && temp == "L1Temp" && paramName == "AlphaL")
            return {0.4, 2.2, 4.0, 8.0, 16.3};
        if (element == "Boron" && temp == "L2Temp" && 
            (paramName == "Peak" || paramName == "SigmaL" || paramName == "SigmaR" || 
             paramName == "AlphaL" || paramName == "AlphaR"))
            return {0.4, 0.75, 1.06, 1.4, 1.9, 6, 16.3};
        
        vector<double> defaultPoints;
        for (int i = 0; i < FluxAnalys::xpointsEkcombine.size(); ++i) {
            defaultPoints.push_back(FluxAnalys::xpointsEkcombine[i]);
        }
        return defaultPoints;
    }
};

class PlotMaker {
private:
    TCanvas* canvas;
    const PlotConfig& config;

    void setupFrame(double xMin, double xMax, const vector<TH1*>& hists, const string& yTitle) {
        canvas->Clear();
        string frameName = "frame_" + to_string(gStorage.allFrames.size());
        TH1F* frame = new TH1F(frameName.c_str(), "", 100, xMin, xMax);
        gStorage.addFrame(frame);
        
        auto yRange = HistManager::getYRange(hists);
        frame->SetStats(0);
        frame->SetTitle("");
        frame->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
        frame->GetYaxis()->SetTitle(yTitle.c_str());
        frame->GetYaxis()->SetRangeUser(yRange.first, yRange.second);
        frame->Draw();
    }

    void drawHistograms(const vector<TH1*>& hists, const vector<int>& colors, double xMin, double xMax) {
        for (size_t i = 0; i < hists.size(); i++) {
            if (!hists[i]) continue;
            hists[i]->GetXaxis()->SetRangeUser(xMin, xMax);
            hists[i]->SetLineColor(colors[i]);
            hists[i]->SetMarkerColor(colors[i]);
            hists[i]->SetLineWidth(2);
            hists[i]->SetMarkerStyle(20);
            hists[i]->Draw("E SAME");
        }
    }

    void addTitleAndLabels(const string& title, const vector<string>& labels, 
                          const vector<int>& colors, double labelPos = 0.75) {
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.045);
        latex.SetTextFont(42);
        latex.DrawLatex(0.0, 0.93, title.c_str());

        double yPos = 0.9;
        for (size_t i = 0; i < labels.size(); i++) {
            latex.SetTextColor(colors[i]);
            latex.DrawLatex(labelPos, yPos, labels[i].c_str());
            yPos -= 0.06;
        }
    }

public:
    PlotMaker(TCanvas* c, const PlotConfig& cfg) : canvas(c), config(cfg) {}

    void plotDetectorComparison(const vector<TH1*>& hists, const vector<int>& colors,
                               const vector<string>& labels, const string& title,
                               const string& yTitle, double xMin, double xMax, TF1* splineFit = nullptr) {
        if (hists.empty()) return;

        setupFrame(xMin, xMax, hists, yTitle);
        drawHistograms(hists, colors, xMin, xMax);

        if (splineFit) {
            splineFit->Draw("SAME");
            TLatex chi2Latex;
            chi2Latex.SetNDC();
            chi2Latex.SetTextSize(0.035);
            chi2Latex.SetTextFont(42);
            chi2Latex.SetTextColor(kBlack);
            double chi2 = splineFit->GetChisquare();
            int ndf = splineFit->GetNDF();
            double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;
            chi2Latex.DrawLatex(0.15, 0.85, Form("Chi2/NDF = %.2f/%d = %.2f", chi2, ndf, chi2ndf));
        }

        // Draw vertical lines
        auto yRange = HistManager::getYRange(hists);
        TLine* line1 = new TLine(1.17, yRange.first, 1.17, yRange.second);
        line1->SetLineColor(kRed);
        line1->SetLineStyle(2);
        line1->SetLineWidth(2);
        line1->Draw();

        TLine* line2 = new TLine(3.23, yRange.first, 3.23, yRange.second);
        line2->SetLineColor(kRed);
        line2->SetLineStyle(2);
        line2->SetLineWidth(2);
        line2->Draw();

        addTitleAndLabels(title, labels, colors, 0.75);
        canvas->SetLogx();
        canvas->Print(canvas->GetName());
    }

    // MODIFIED: plotTempComparison uses fixed segment spline fits
    void plotTempComparison(const vector<TH1*>& hists, const vector<int>& colors,
                           const vector<string>& labels, const string& title,
                           const string& yTitle, double xMin, double xMax,
                           const string& element, const string& detector, const pair<string, string>& param) {
        if (hists.empty()) return;

        // Skip Chi2 fitting
        if (title.find("Chi2NDF") != string::npos) {
            setupFrame(xMin, xMax, hists, yTitle);
            drawHistograms(hists, colors, xMin, xMax);
            addTitleAndLabels(title, labels, colors, 0.75);
            canvas->SetLogx();
            canvas->Print(canvas->GetName());
            return;
        }

        setupFrame(xMin, xMax, hists, yTitle);
        drawHistograms(hists, colors, xMin, xMax);

        TLatex chi2Latex;
        chi2Latex.SetNDC();
        chi2Latex.SetTextSize(0.035);
        chi2Latex.SetTextFont(42);

        // Perform spline fits with fixed segments and draw them
        for (size_t i = 0; i < hists.size(); i++) {
            if (!hists[i]) continue;
            
            // Construct key to get segment count
            string key = element + "_" + detector + "_" + labels[i] + "_" + param.first + "_" + param.second;
            int segments = config.splineSegments.count(key) ? config.splineSegments.at(key) : 1; // Default to 1 segment

            string baseName = key + "_splinefit";
            TF1* fit = FittingUtils::performSplineFit(hists[i], detector, baseName, labels[i], segments);
            
            if (fit) {
                fit->Draw("SAME");
                
                // Display Chi2/NDF for the fit
                double chi2 = fit->GetChisquare();
                int ndf = fit->GetNDF();
                double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;
                
                chi2Latex.SetTextColor(fit->GetLineColor());
                chi2Latex.DrawLatex(0.15, 0.85 - i * 0.05, 
                                   Form("%s (%d Seg): #chi^{2}/NDF = %.1f/%d = %.2f", 
                                        labels[i].c_str(), segments, chi2, ndf, chi2ndf));
            }
        }

        addTitleAndLabels(title, labels, colors, 0.75);
        canvas->SetLogx();
        canvas->Print(canvas->GetName());
    }

    void plotComparison(const vector<TH1*>& hists, const vector<int>& colors,
                       const vector<string>& labels, const string& title,
                       const string& yTitle, double xMin, double xMax,
                       TF1* extraFit = nullptr, double labelPos = 0.75) {
        if (hists.empty()) return;

        setupFrame(xMin, xMax, hists, yTitle);
        drawHistograms(hists, colors, xMin, xMax);

        if (extraFit) extraFit->Draw("SAME");

        addTitleAndLabels(title, labels, colors, labelPos);
        canvas->SetLogx();
        canvas->Print(canvas->GetName());
    }

    vector<string> generateHistNames(const string& element, const vector<string>& variables,
                                   const string& temp, const pair<string, string>& param) {
        vector<string> names;
        for (const auto& var : variables) {
            names.push_back(element + "_" + var + "_" + temp + "_" + param.first + "_" + param.second);
        }
        return names;
    }
};

void saveSplineFits(const string& outputFile) {
    
    unique_ptr<TFile> outFile(TFile::Open(outputFile.c_str(), "RECREATE"));
    if (!outFile || outFile->IsZombie()) return;
    
    for (const auto& fitPair : gStorage.splineFits) {
        if (fitPair.first) {
            TF1* clonedFit = static_cast<TF1*>(fitPair.first->Clone(fitPair.second.c_str()));
            clonedFit->Write();
        }
    }
    outFile->Close();
}

// Main plotting functions
void plotDetectorComparison(TFile* fin, TCanvas* c, const PlotConfig& config) {
    PlotMaker plotter(c, config);
    
    for (const auto& element : config.elements) {
        for (const auto& temp : config.temps) {
            for (const auto& param : config.paramList) {
                auto histNames = plotter.generateHistNames(element, config.detectors, temp, param);
                auto hists = HistManager::getHists(fin, histNames);
                if (hists.empty()) continue;

                auto colors = vector<int>{config.detColors.at("TOF"), config.detColors.at("NaF"), config.detColors.at("AGL")};
                auto xRange = config.getXRange("", true);
                
                auto splineFit = FittingUtils::createCombinedSplineFit(hists, config.detectors, param.second, element, temp);

                if (splineFit) {
                    string splineName = element + "_Combined_" + temp + "_" + param.first + "_" + param.second + "_splinefit";
                    gStorage.addSplineFit(splineFit, splineName);
                }

                plotter.plotDetectorComparison(hists, colors, config.detectors,
                    element + "  " + temp + "  " + param.first + "_" + param.second,
                    param.second, xRange.first, xRange.second, splineFit);
            }
        }
    }
}

// MODIFIED: This function's call to plotter.plotTempComparison now passes more arguments
void plotTempComparison(TFile* fin, TCanvas* c, const PlotConfig& config) {
    PlotMaker plotter(c, config);
    
    for (const auto& element : config.elements) {
        for (const auto& det : config.detectors) {
            for (const auto& param : config.paramList) {
                vector<string> tempHistNames;
                for (const auto& temp : config.temps) {
                    tempHistNames.push_back(element + "_" + det + "_" + temp + "_" + param.first + "_" + param.second);
                }
                
                auto hists = HistManager::getHists(fin, tempHistNames);
                if (hists.empty()) continue;

                vector<int> colors = {config.tempColors.at("L1Temp"), config.tempColors.at("L2Temp")};
                auto xRange = config.getXRange(det);

                plotter.plotTempComparison(hists, colors, config.temps,
                    element + "  " + det + "  " + param.first + "_" + param.second,
                    param.second, xRange.first, xRange.second, 
                    element, det, param); // Pass extra info for segment lookup
            }
        }
    }
}

void plotTemplateComparison(TFile* fin, TCanvas* c, const PlotConfig& config, 
                           const vector<string>& compareTemplates, const string& title) {
    PlotMaker plotter(c, config);
    
    for (const auto& element : config.elements) {
        for (const auto& det : config.detectors) {
            for (const auto& param : config.paramList) {
                vector<string> histNames;
                vector<int> colors;
                
                for (const auto& temp : compareTemplates) {
                    histNames.push_back(element + "_" + det + "_" + temp + "_" + param.first + "_" + param.second);
                    colors.push_back(config.tempColors.at(temp));
                }
                
                auto hists = HistManager::getHists(fin, histNames);
                if (hists.empty()) continue;

                auto xRange = config.getXRange(det);
                plotter.plotComparison(hists, colors, compareTemplates,
                    element + "  " + det + "  " + title, param.second, xRange.first, xRange.second);
            }
        }
    }
}

void plotL1TemplateUnbComparison(TFile* fin, TCanvas* c, const PlotConfig& config) {
    plotTemplateComparison(fin, c, config, config.L1_compare, "L1 unb vs L1");
}

void plotL2TemplateUnbComparison(TFile* fin, TCanvas* c, const PlotConfig& config) {
    plotTemplateComparison(fin, c, config, config.L2_compare, "L2 unb vs L2");
}

void plotCutComparison(const vector<TFile*>& fins, TCanvas* c, const PlotConfig& config) {
    PlotMaker plotter(c, config);
    vector<string> cutLabels = {"cut coe = 0.8", "cut coe = 0.5", "cut coe = 0.2"};
    vector<int> cutColors = {kBlack, kRed, kBlue};
    
    for (const auto& element : config.elements) {
        for (const auto& det : config.detectors) {
            for (const auto& temp : config.temps_unb) {
                for (const auto& param : config.paramList) {
                    string histName = element + "_" + det + "_" + temp + "_" + param.first + "_" + param.second;
                    auto hists = HistManager::getHistsMultiFile(fins, histName);
                    if (hists.empty()) continue;

                    auto xRange = config.getXRange(det);
                    plotter.plotComparison(hists, cutColors, cutLabels,
                        element + "  " + det + "  " + temp + "  " + param.first + "_" + param.second,
                        param.second, xRange.first, xRange.second);
                }
            }
        }
    }
}

void plotModelChi2Comparison(TFile* fin, TCanvas* c, const PlotConfig& config) {
    PlotMaker plotter(c, config);
    vector<string> fitTypes = {"LG", "EGE"};
    vector<string> modelLabels = {"Landau-Gauss", "ExpGausExp"};
    
    for (const auto& element : config.elements) {
        for (const auto& det : config.detectors) {
            for (const auto& temp : config.temps_unb) {
                vector<string> histNames;
                vector<int> colors;
                
                for (const auto& fitType : fitTypes) {
                    histNames.push_back(element + "_" + det + "_" + temp + "_" + fitType + "_Chi2NDF");
                    colors.push_back(config.modelColors.at(fitType));
                }
                
                auto hists = HistManager::getHists(fin, histNames);
                if (hists.empty()) continue;

                auto xRange = config.getXRange(det);
                plotter.plotComparison(hists, colors, modelLabels,
                    element + "  " + det + "  " + temp + "  Chi2/NDF",
                    "#chi^{2}/ndf", xRange.first, xRange.second, nullptr, 0.65);
            }
        }
    }
}

template<typename PlotFunc>
void createPDF(TFile* fin, TCanvas* c, const string& pdfName, const PlotConfig& config, PlotFunc plotFunc) {
    c->Print((pdfName + "[").c_str());
    plotFunc(fin, c, config);
    c->Print((pdfName + "]").c_str());
}

template<typename PlotFunc>
void createPDFMultiFile(const vector<TFile*>& fins, TCanvas* c, const string& pdfName,
                       const PlotConfig& config, PlotFunc plotFunc) {
    c->Print((pdfName + "[").c_str());
    plotFunc(fins, c, config);
    c->Print((pdfName + "]").c_str());
}

void compareFitResults(
    const string& rootfile02 = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/allFitHist_0p2.root",
    const string& rootfile05 = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/allFitHist_0p5.root",
    const string& rootfile08 = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/allFitHist_0p8.root",
    const string& outdir = "/eos/user/z/zixuan/Isotope/chargeTempFit/"
) {
    gStyle->SetOptStat(0);
    PlotConfig config;

    const string primaryFileToProcess = rootfile08;
    const string versionTag = "0p8";
    
    unique_ptr<TFile> fin02(TFile::Open(rootfile02.c_str()));
    unique_ptr<TFile> fin05(TFile::Open(rootfile05.c_str()));
    unique_ptr<TFile> fin08(TFile::Open(rootfile08.c_str()));
    
    if (!fin02 || fin02->IsZombie() || !fin05 || fin05->IsZombie() || !fin08 || fin08->IsZombie()) {
        cout << "ERROR: Cannot open one or more input files" << endl;
        return;
    }

    TFile* primaryFilePtr = nullptr;
    if (primaryFileToProcess == rootfile08) primaryFilePtr = fin08.get();
    else if (primaryFileToProcess == rootfile05) primaryFilePtr = fin05.get();
    else if (primaryFileToProcess == rootfile02) primaryFilePtr = fin02.get();

    auto c = make_unique<TCanvas>("c", "", 800, 400);
    c->SetGrid();
    
    vector<TFile*> cutFiles = {fin08.get(), fin05.get(), fin02.get()};
    
    try {
        string pdfName;

        pdfName = outdir + "FitParam_TempCompare_" + versionTag + ".pdf";
        c->SetName(pdfName.c_str());
        createPDF(primaryFilePtr, c.get(), pdfName, config, plotTempComparison);
       
        pdfName = outdir + "FitParam_DetectorCompare_" + versionTag + ".pdf";
        c->SetName(pdfName.c_str());
        createPDF(primaryFilePtr, c.get(), pdfName, config, plotDetectorComparison);
       
        /*
        pdfName = outdir + "FitParam_L1UnbCompare_" + versionTag + ".pdf";
        c->SetName(pdfName.c_str());
        createPDF(primaryFilePtr, c.get(), pdfName, config, plotL1TemplateUnbComparison);

        pdfName = outdir + "FitParam_L2UnbCompare_" + versionTag + ".pdf";
        c->SetName(pdfName.c_str());
        createPDF(primaryFilePtr, c.get(), pdfName, config, plotL2TemplateUnbComparison);

        pdfName = outdir + "FitParam_Chi2ModelCompare_" + versionTag + ".pdf";
        c->SetName(pdfName.c_str());
        createPDF(primaryFilePtr, c.get(), pdfName, config, plotModelChi2Comparison);

        pdfName = outdir + "FitParam_CutCompare.pdf";
        c->SetName(pdfName.c_str());
        createPDFMultiFile(cutFiles, c.get(), pdfName, config, plotCutComparison);
        
        */
        string splineOutFile = outdir + "allFitHistSplineSmooth_" + versionTag + ".root";
        saveSplineFits(splineOutFile);

    } catch (const exception& e) {
        cout << "ERROR: Exception during PDF generation: " << e.what() << endl;
    }
}