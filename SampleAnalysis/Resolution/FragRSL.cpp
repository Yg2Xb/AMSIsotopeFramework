#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iomanip>
#include <sstream>
#include <TStyle.h>
#include <TLine.h>
#include "../Tool.h"

using namespace std;
using namespace AMS_Iso;

struct HistInfo {
    string suffix;
    string title_part;
    string x_axis_label;
    int bin_min;
    int bin_max;
    double min_bin_width;
};

struct FitResult {
    double mean, mean_err;
    double sigma, sigma_err;
};

// Modified to include Element name for histogram lookup
struct FileConfig {
    string filename;
    string element; // "Boron", "Carbon", etc.
    double weight;
};

// Signal files (*2, *3): All are Beryllium isotopes
const vector<FileConfig> FILES_SIG = {
    {"Be7_rew_frag4.root",  "Beryllium", 1.8},
    {"Be9_rew_frag4.root",  "Beryllium", 0.9},
    {"Be10_rew_frag4.root", "Beryllium", 0.3}
};

// Background/Frag files (*4): Mixed elements
const vector<FileConfig> FILES_BKG = {
    {"B10_rew_frag4.root",  "Boron",    0.15},
    {"B11_rew_frag4.root",  "Boron",    0.35},
    {"C12_rew_frag4.root",  "Carbon",   2.50},
    {"N14_rew_frag4.root",  "Nitrogen", 0.25},
    {"N15_rew_frag4.root",  "Nitrogen", 0.25},
    {"O16_rew_frag4.root",  "Oxygen",   2.50}
};

void SetGraphStyle(TGraphErrors* g, int color, int marker) {
    g->SetMarkerStyle(marker); 
    g->SetMarkerSize(1.2);
    g->SetLineColor(color); 
    g->SetMarkerColor(color);
}

void AutoRebin(TH1* h, double min_width) {
    for(int r = 1; r <= 10; ++r) {
        if (h->GetMaximum() >= 50 && h->GetBinWidth(1) > min_width) break;
        h->Rebin(2);
    }
}

void findFitRange(TH1* hist, double coverage, double center_x, double& x_min, double& x_max) {
    if (!hist || hist->GetEntries() == 0) return;
    double total = hist->Integral();
    if (total <= 0) return;
    
    int center_bin = hist->GetXaxis()->FindFixBin(center_x);
    if (center_bin < 1 || center_bin > hist->GetNbinsX()) {
        center_bin = hist->GetMaximumBin(); 
        center_x = hist->GetXaxis()->GetBinCenter(center_bin);
    }
    
    double required = total * coverage;
    double current = hist->GetBinContent(center_bin);
    int low = center_bin - 1, high = center_bin + 1;
    x_min = hist->GetXaxis()->GetBinLowEdge(center_bin);
    x_max = hist->GetXaxis()->GetBinUpEdge(center_bin);

    int n_bins = hist->GetNbinsX();
    while (current < required) {
        double c_low = (low >= 1) ? hist->GetBinContent(low) : 0;
        double c_high = (high <= n_bins) ? hist->GetBinContent(high) : 0;
        
        bool extended = false;
        if (low >= 1 && high <= n_bins) {
            current += c_low + c_high;
            x_min = hist->GetXaxis()->GetBinLowEdge(low--);
            x_max = hist->GetXaxis()->GetBinUpEdge(high++);
            extended = true;
        } else if (low >= 1) {
            current += c_low;
            x_min = hist->GetXaxis()->GetBinLowEdge(low--);
            extended = true;
        } else if (high <= n_bins) {
            current += c_high;
            x_max = hist->GetXaxis()->GetBinUpEdge(high++);
            extended = true;
        }
        if (!extended) break;
    }
    
    double final_center = (x_min + x_max) / 2.0;
    double half_range = max(abs(x_max - final_center), abs(x_min - final_center));
    x_min = final_center - half_range;
    x_max = final_center + half_range;
}

void getGlobalYRange(const vector<TGraphErrors*>& graphs, double& y_min, double& y_max) {
    y_min = 1e10; y_max = -1e10;
    bool found = false;
    for (const auto& g : graphs) {
        if(!g) continue;
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            y_min = min(y_min, y); 
            y_max = max(y_max, y);
            found = true;
        }
    }
    if (found) {
        double r = (y_max - y_min == 0) ? (abs(y_min)*0.1 ? abs(y_min)*0.1 : 0.01) : (y_max - y_min);
        y_min -= r * 0.1;
        y_max += r * 0.1;
    } else {
        y_min = -0.01; y_max = 0.01;
    }
}

HistInfo getHistInfo(const string& suffix) {
    if (suffix == "ID_H5a2") return {"a2", "NaF-Tracker #Delta(1/#beta)", "gene Rigidity [GV]", 10, 25, 25e-5};
    if (suffix == "ID_H5b2") return {"b2", "AGL-Tracker #Delta(1/#beta)", "gene Rigidity [GV]", 20, 47, 5e-5};
    
    if (suffix == "ID_H5a3") return {"a3", "NaF-Gene #Delta(1/#beta)", "gene Beta", 10, 25, 25e-5};
    if (suffix == "ID_H5b3") return {"b3", "AGL-Gene #Delta(1/#beta)", "gene Beta", 20, 47, 5e-5};
    if (suffix == "ID_H5c3") return {"c3", "TOF-Gene #Delta(1/#beta)", "gene Beta", 7, 14, 5e-5};
    if (suffix == "ID_H5d3") return {"d3", "Track-Gene #Delta(1/rig)", "gene Rigidity [GV]", 6, 55, 25e-5};

    if (suffix == "ID_H5a4") return {"a4", "frag NaF-L2True #Delta(1/#beta)", "L2True Beta", 10, 25, 5e-5};
    if (suffix == "ID_H5b4") return {"b4", "frag AGL-L2True #Delta(1/#beta)", "L2True Beta", 20, 47, 5e-5};
    if (suffix == "ID_H5c4") return {"c4", "frag TOF-L2True #Delta(1/#beta)", "L2True Beta", 7, 14, 5e-5};
    if (suffix == "ID_H5d4") return {"d4", "frag Track-L2True #Delta(1/rig)", "L2True Rigidity [GV]", 6, 55, 5e-5};

    return {"", "", "", 0, 0, 0};
}

// Updated function to construct histogram name dynamically based on element
TH2F* LoadCombinedHist(const string& suffix, const vector<FileConfig>& files, const string& input_dir) {
    TH2F* h_total = nullptr;

    for (const auto& fc : files) {
        TFile* f = TFile::Open((input_dir + fc.filename).c_str(), "READ");
        if (!f || f->IsZombie()) {
            if (f) delete f;
            continue;
        }

        // Dynamically construct name: UnbiasedL1Inner_{Element}_{Suffix}
        string full_name = "UnbiasedL1Inner_" + fc.element + "_" + suffix;
        
        TH2F* h = dynamic_cast<TH2F*>(f->Get(full_name.c_str()));
        if (h) {
            if (!h_total) {
                h_total = (TH2F*)h->Clone(("comb_" + suffix).c_str());
                h_total->SetDirectory(nullptr);
                h_total->Scale(fc.weight);
            } else {
                h_total->Add(h, fc.weight);
            }
        } else {
            // Optional debug output
            // cout << "Warning: " << full_name << " not found in " << fc.filename << endl;
        }
        f->Close(); delete f;
    }
    return h_total;
}

FitResult twoStepGaussianFit(TH1* hist, const string& title, double center_x, double min_width) {
    FitResult result = {0, 0, 0, 0};
    if (!hist || hist->GetEntries() < 10) return result;
    hist->Sumw2();
    
    AutoRebin(hist, min_width);

    double x1 = 0, x2 = 0;
    findFitRange(hist, 0.85, center_x, x1, x2);
    if (x2 <= x1) return result;

    TF1* f1 = new TF1("f1", "gaus", x1, x2);
    f1->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS());
    if (hist->Fit(f1, "QRS") != 0) { delete f1; return result; }
    
    double mean0 = f1->GetParameter(1);
    double sigma0 = f1->GetParameter(2);
    delete f1;

    double x_min = mean0 - 2 * abs(sigma0);
    double x_max = mean0 + 2 * abs(sigma0);
    if (x_max <= x_min) return result;

    TCanvas* c = (TCanvas*)gROOT->FindObject("c_fit"); 
    hist->SetTitle(title.c_str());
    double buf = 0.5 * (x_max - x_min);
    hist->GetXaxis()->SetRangeUser(x_min - buf, x_max + buf);

    vector<vector<double>> data = DoGausPlusAsymGausFit(hist, x_min, x_max, c, false);

    if (data.size() >= 2 && data[1][4] > 0) { 
        result = {data[0][0], data[1][0], data[0][1], data[1][1]};
        
        c->cd();
        TLine *l1 = new TLine(x_min, 0, x_min, hist->GetMaximum()), *l2 = new TLine(x_max, 0, x_max, hist->GetMaximum());
        l1->SetLineStyle(2); l1->SetLineColor(kRed); l1->Draw("same");
        l2->SetLineStyle(2); l2->SetLineColor(kRed); l2->Draw("same");
    } else {
        c->cd(); hist->Draw("hist");
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
        lat.DrawLatex(0.6, 0.85, "Fit Failed");
    }
    c->Update();
    return result;
}

void DrawAndSave(TGraphErrors* g_mean, TGraphErrors* g_sigma, const HistInfo& info, const string& output_dir, TFile* root_file) {
    if (!g_mean || g_mean->GetN() == 0) return;

    string base_name = "Combined_" + info.suffix;
    string canvas_name = "c_" + base_name;
    
    TCanvas* c = new TCanvas(canvas_name.c_str(), base_name.c_str(), 1400, 500);
    c->Divide(2, 1);

    SetGraphStyle(g_mean, kRed, 20);
    SetGraphStyle(g_sigma, kBlue, 20);

    vector<TGraphErrors*> v_mean = {g_mean};
    vector<TGraphErrors*> v_sigma = {g_sigma};

    double ym_min, ym_max, ys_min, ys_max;
    getGlobalYRange(v_mean, ym_min, ym_max);
    getGlobalYRange(v_sigma, ys_min, ys_max);

    c->cd(1);
    gPad->SetGrid();
    gPad->SetLogx();
    g_mean->SetTitle((info.title_part + " Mean").c_str());
    g_mean->GetXaxis()->SetTitle(info.x_axis_label.c_str());
    g_mean->GetYaxis()->SetTitle("Mean");
    g_mean->GetYaxis()->SetRangeUser(ym_min, ym_max);
    g_mean->Draw("APZ");

    c->cd(2);
    gPad->SetGrid();
    gPad->SetLogx();
    g_sigma->SetTitle((info.title_part + " Sigma").c_str());
    g_sigma->GetXaxis()->SetTitle(info.x_axis_label.c_str());
    g_sigma->GetYaxis()->SetTitle("Sigma");
    g_sigma->GetYaxis()->SetRangeUser(ys_min, ys_max);
    g_sigma->Draw("APZ");

    c->SaveAs((output_dir + base_name + ".png").c_str());
    
    if (root_file) {
        root_file->cd();
        g_mean->SetName(("g_mean_" + base_name).c_str());
        g_mean->Write();
        g_sigma->SetName(("g_sigma_" + base_name).c_str());
        g_sigma->Write();
    }
    delete c;
}

void DrawCombined(const map<string, TGraphErrors*>& graphs, const string& detector, const string& output_dir) {
    TCanvas* c = new TCanvas(("c_comb_" + detector).c_str(), "Combined", 1400, 500);
    c->Divide(2, 1);

    string s3 = "ID_H5" + detector + "3";
    string s4 = "ID_H5" + detector + "4";
    
    if (graphs.count(s3 + "_Mean") == 0 || graphs.count(s4 + "_Mean") == 0) { delete c; return; }

    TGraphErrors* gm3 = graphs.at(s3 + "_Mean");
    TGraphErrors* gm4 = graphs.at(s4 + "_Mean");
    TGraphErrors* gs3 = graphs.at(s3 + "_Sigma");
    TGraphErrors* gs4 = graphs.at(s4 + "_Sigma");

    SetGraphStyle(gm3, kBlue, 20); SetGraphStyle(gs3, kBlue, 20);
    SetGraphStyle(gm4, kRed, 20);  SetGraphStyle(gs4, kRed, 20);

    vector<TGraphErrors*> vm = {gm3, gm4}, vs = {gs3, gs4};
    double ym_min, ym_max, ys_min, ys_max;
    getGlobalYRange(vm, ym_min, ym_max);
    getGlobalYRange(vs, ys_min, ys_max);

    HistInfo info = getHistInfo(s3);

    auto drawPad = [&](int pad, TGraphErrors* g3, TGraphErrors* g4, const string& type, double min, double max) {
        c->cd(pad); gPad->SetGrid(); gPad->SetLogx();
        g3->SetTitle((info.title_part + " " + type + " Comparison").c_str());
        g3->GetXaxis()->SetTitle(info.x_axis_label.c_str());
        g3->GetYaxis()->SetTitle(type.c_str());
        g3->GetYaxis()->SetRangeUser(min, max);
        g3->Draw("AP"); 
        g4->Draw("P same"); 
        
        TLegend* leg = new TLegend(0.75, 0.8, 0.9, 0.9);
        leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->AddEntry(g3, "Non-Frag (Be)", "p");
        leg->AddEntry(g4, "Frag (Mix)", "p");
        leg->Draw();
    };

    drawPad(1, gm3, gm4, "Mean", ym_min, ym_max);
    drawPad(2, gs3, gs4, "Sigma", ys_min, ys_max);

    c->SaveAs((output_dir + "Combined_Compare_" + detector + ".png").c_str());
    delete c;
}

void Analyze(const string& suffix, const string& input_dir, const string& output_dir, TFile* out_root, map<string, TGraphErrors*>& all_graphs) {
    HistInfo info = getHistInfo(suffix);
    if (info.suffix.empty()) return;

    cout << "Analyzing " << suffix << " ..." << endl;
    
    // Choose file list based on suffix ending (3/2 -> Signal, 4 -> Background)
    // Assuming *2 and *3 are signal (Be), and *4 is background (Frag)
    const vector<FileConfig>& files = (suffix.back() == '4') ? FILES_BKG : FILES_SIG;

    TH2F* h2 = LoadCombinedHist(suffix, files, input_dir);
    if (!h2) {
        cerr << "Warning: Failed to load combined histogram for " << suffix << endl;
        return;
    }

    string pdf_name = output_dir + "FitResults_Combined_" + suffix + ".pdf";
    TCanvas* c_fit = new TCanvas("c_fit", "Fits", 800, 600);
    c_fit->Print((pdf_name + "[").c_str(), "pdf");

    TGraphErrors* g_mean = new TGraphErrors();
    TGraphErrors* g_sigma = new TGraphErrors();

    int nBinsY = h2->GetNbinsY();
    
    for (int i = 1; i <= nBinsY; ++i) {
        if (i < info.bin_min || (info.bin_max != 999 && i > info.bin_max)) continue;

        double x_val = h2->GetYaxis()->GetBinCenter(i);
        
        TH1D* h1 = h2->ProjectionX(Form("px_%s_%d", suffix.c_str(), i), i, i);
        h1->SetDirectory(nullptr);
        
        if (h1->GetEntries() > 50) {
            stringstream ss; ss << fixed << setprecision(3) << x_val;
            string label = "Val: " + ss.str();
            
            FitResult res = twoStepGaussianFit(h1, info.title_part + " " + label, 0.0, info.min_bin_width);
            
            if (res.mean_err > 0 || res.sigma_err > 0) {
                int n = g_mean->GetN();
                g_mean->SetPoint(n, x_val, res.mean);
                g_mean->SetPointError(n, 0.0, res.mean_err);
                g_sigma->SetPoint(n, x_val, res.sigma);
                g_sigma->SetPointError(n, 0.0, res.sigma_err);
                c_fit->Print(pdf_name.c_str(), "pdf");
            }
        }
        delete h1;
    }

    DrawAndSave(g_mean, g_sigma, info, output_dir, out_root);
    
    all_graphs[suffix + "_Mean"] = g_mean;
    all_graphs[suffix + "_Sigma"] = g_sigma;

    c_fit->Print((pdf_name + "]").c_str(), "pdf");
    delete h2; 
    delete c_fit;
}

void FragRSL() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetErrorX(0); gStyle->SetOptFit(0);
    
    const string in_dir = "/eos/user/z/zixuan/Isotope/Add/";
    const string out_dir = "/eos/user/z/zixuan/Isotope/Resolution/FragRSL/";
    const string out_root = out_dir + "FragRSL_Combined_Final.root";

    if (gSystem->AccessPathName(out_dir.c_str())) gSystem->mkdir(out_dir.c_str(), kTRUE);

    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) { cerr << "Error output file" << endl; return; }

    map<string, TGraphErrors*> all_graphs;

    vector<string> suffixes = {
        "ID_H5a2", "ID_H5b2",
        "ID_H5a3", "ID_H5b3", "ID_H5c3", "ID_H5d3",
        "ID_H5a4", "ID_H5b4", "ID_H5c4", "ID_H5d4"
    };

    for (const auto& s : suffixes) {
        Analyze(s, in_dir, out_dir, fout, all_graphs);
    }

    vector<string> detectors = {"a", "b", "c", "d"};
    for(const auto& det : detectors) {
        DrawCombined(all_graphs, det, out_dir);
    }

    fout->Close(); delete fout;
    gROOT->SetBatch(kFALSE);
    cout << "\nDone. Results in: " << out_dir << endl;
}