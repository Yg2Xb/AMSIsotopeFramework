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

// --- Data Structures ---

struct HistInfo {
    string y_axis_label;
    string output_suffix;
    string title_description;
};

struct FitResult {
    double mean, mean_err;
    double sigma, sigma_err;
    double chi2, ndf;
    double LR, LR_err;
    double RR, RR_err;
};

// Unified configuration for MC files
struct MCConfig {
    string filename; // Input filename
    string nuclide;  // Label (e.g., "B10")
    string element;  // Element name for hist lookup (e.g., "Boron")
    double Z;        // Charge for plotting
};

// --- Constants & Config ---

const vector<double> H5_RIG_BINS_EDGES = {30.0, 50.0, 80.0, 120.0, 160.0, 240.0};
const vector<string> H5_RIG_LABELS = {"30-50 GV", "50-80 GV", "80-120 GV", "120-160 GV", "160-240 GV"};
const int N_RIG_BINS = H5_RIG_BINS_EDGES.size() - 1;

const vector<int> COLORS = {kYellow+2, kMagenta, kBlack, kRed, kBlue, kGreen + 2, kOrange + 1, kViolet, kCyan};

// Master list of MC definitions
const vector<MCConfig> MC_DEFINITIONS = {
    {"B10_rew_frag4.root",  "B10",  "Boron",     5.0},
    {"B11_rew_frag4.root",  "B11",  "Boron",     5.0},
    {"Be10_rew_frag4.root", "Be10", "Beryllium", 4.0},
    {"Be7_rew_frag4.root",  "Be7",  "Beryllium", 4.0},
    {"Be9_rew_frag4.root",  "Be9",  "Beryllium", 4.0},
    {"C12_rew_frag4.root",  "C12",  "Carbon",    6.0},
    {"N15_rew_frag4.root",  "N15",  "Nitrogen",  7.0},
    {"O16_rew_frag4.root",  "O16",  "Oxygen",    8.0}
};

// --- Helper Functions ---

void SetGraphStyle(TGraphErrors* g, int color, int style_offset = 0) {
    g->SetMarkerStyle(20 + style_offset); 
    g->SetMarkerSize(1.0);
    g->SetLineColor(color); 
    g->SetMarkerColor(color);
}

// Custom Rebin logic for MC (different from ISS)
void AutoRebin(TH1* h) {
    for(int r = 1; r <= 10; ++r) {
        if (h->GetMaximum() >= 60 || h->GetNbinsX() < 20) break;
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

void getGlobalYRange(const vector<TGraphErrors*>& graphs, double x_min, double x_max, double& y_min, double& y_max) {
    y_min = 1e10; y_max = -1e10;
    bool found = false;
    for (const auto& g : graphs) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            if (x >= x_min && x <= x_max) {
                y_min = min(y_min, y - g->GetErrorY(i));
                y_max = max(y_max, y + g->GetErrorY(i));
                found = true;
            }
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
    if (suffix == "ID_H5a") return {"Rigidity [GV]", "Rigidity", "NaF-Tracker #Delta(1/#beta)"};
    if (suffix == "ID_H5b") return {"Rigidity [GV]", "Rigidity", "AGL-Tracker #Delta(1/#beta)"};
    if (suffix == "ID_H4a") return {"Rigidity [GV]", "Rigidity", "NaF 1/#beta (Rig > 80GV)"};
    if (suffix == "ID_H4b") return {"Rigidity [GV]", "Rigidity", "AGL 1/#beta (Rig > 150GV)"};
    return {"Y Variable", "YVariable", "Unknown Delta Beta"};
}

// Core drawing function
void DrawAndSaveGraphs(vector<TGraphErrors*>& graphs, const vector<string>& legend_labels,
                       const string& x_title, const string& y_title, const string& title_prefix, 
                       const string& output_dir, const string& output_base, 
                       TFile* root_file, const string& suffix) {
    if (graphs.empty()) return;

    string canvas_name = "c_" + output_base + suffix;
    TCanvas* c = new TCanvas(canvas_name.c_str(), (title_prefix + suffix).c_str(), 700, 500);
    c->SetGrid();

    // MC specific X-axis range logic
    double x_min, x_max;
    if (suffix.find("Charge") != string::npos) {
        x_min = 1.5; x_max = 8.5; 
    } else {
        x_min = H5_RIG_BINS_EDGES.front();
        x_max = H5_RIG_BINS_EDGES.back() * 1.5;
    }

    double y_min, y_max;
    getGlobalYRange(graphs, x_min, x_max, y_min, y_max);

    TLegend* leg = new TLegend(0.8, 0.8, 0.99, 0.99);
    leg->SetFillStyle(0); leg->SetBorderSize(1);

    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->Draw(i == 0 ? "APZ" : "PZ same");
        if (i == 0) {
            graphs[i]->SetTitle((title_prefix + " vs " + x_title).c_str());
            graphs[i]->GetXaxis()->SetTitle(x_title.c_str());
            graphs[i]->GetYaxis()->SetTitle(y_title.c_str());
            graphs[i]->GetXaxis()->SetRangeUser(x_min, x_max);
            graphs[i]->GetYaxis()->SetRangeUser(y_min, y_max);
        }
        if (i < legend_labels.size()) leg->AddEntry(graphs[i], legend_labels[i].c_str(), "p");

        if (root_file) {
            root_file->cd();
            graphs[i]->SetName(("g_" + output_base + suffix + "_" + to_string(i)).c_str());
            graphs[i]->Write();
        }
    }
    leg->Draw();
    c->SaveAs((output_dir + output_base + suffix + ".png").c_str());
    delete leg; delete c;
}

// Combine H4 results for all nuclides
void DrawMC_H4_Combined(const map<string, FitResult>& all_results, const string& output_dir, const HistInfo& info, 
                      const string& output_base, const string& x_label, 
                      TFile* root_file) {
    
    TGraphErrors* g_mean = new TGraphErrors();
    TGraphErrors* g_sigma = new TGraphErrors();
    SetGraphStyle(g_mean, kRed);
    SetGraphStyle(g_sigma, kBlue);

    // Use MC_DEFINITIONS directly to determine Z for each nuclide
    for (const auto& conf : MC_DEFINITIONS) {
        string suffix = (info.output_suffix == "Rigidity") ? "ID_H4a" : "ID_H4b";
        string key = conf.nuclide + "_" + suffix;

        if (all_results.count(key)) {
            const auto& res = all_results.at(key);
            if (res.mean_err > 0 || res.sigma_err > 0) {
                int n = g_mean->GetN();
                g_mean->SetPoint(n, conf.Z, res.mean);
                g_mean->SetPointError(n, 0.0, res.mean_err);
                g_sigma->SetPoint(n, conf.Z, res.sigma);
                g_sigma->SetPointError(n, 0.0, res.sigma_err);
            }
        }
    }

    vector<string> legend = {"MC"};
    vector<TGraphErrors*> v_mean = {g_mean}, v_sigma = {g_sigma};

    DrawAndSaveGraphs(v_mean, legend, "Charge (Z)", ("#mu_{" + x_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_base, root_file, "_Mean_vs_Charge");
    DrawAndSaveGraphs(v_sigma, legend, "Charge (Z)", ("#sigma_{" + x_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_base, root_file, "_Sigma_vs_Charge");

    delete g_mean; delete g_sigma;
}

// Combine H5 results for all nuclides
void DrawMC_H5_Combined(const map<string, TH2F*>& h2_mean_map, const map<string, TH2F*>& h2_sigma_map,
                        const string& output_dir, const HistInfo& info, const string& output_base, 
                        const string& x_label, TFile* root_file) {
    
    vector<TGraphErrors*> gm_rig, gs_rig;
    vector<string> labels;
    int color_idx = 0;

    for (const auto& conf : MC_DEFINITIONS) {
        string nuclide = conf.nuclide;
        if (!h2_mean_map.count(nuclide) || !h2_sigma_map.count(nuclide)) continue;
        
        TH2F* h2_m = h2_mean_map.at(nuclide);
        TH2F* h2_s = h2_sigma_map.at(nuclide);

        TGraphErrors* gm = new TGraphErrors();
        TGraphErrors* gs = new TGraphErrors();
        int color = COLORS[color_idx % COLORS.size()];
        SetGraphStyle(gm, color, color_idx); // Use index as marker style offset
        SetGraphStyle(gs, color, color_idx);

        // H5 MC results stored in 1st column of TH2
        for (int j = 1; j <= N_RIG_BINS; ++j) {
            double m = h2_m->GetBinContent(j, 1);
            double me = h2_m->GetBinError(j, 1);
            double s = h2_s->GetBinContent(j, 1);
            double se = h2_s->GetBinError(j, 1);

            if (me > 0 || se > 0) {
                double rig_center = (H5_RIG_BINS_EDGES[j-1] + H5_RIG_BINS_EDGES[j]) / 2.0;
                int n = gm->GetN();
                gm->SetPoint(n, rig_center, m);
                gm->SetPointError(n, 0.0, me);
                gs->SetPoint(n, rig_center, s);
                gs->SetPointError(n, 0.0, se);
            }
        }
        gm_rig.push_back(gm);
        gs_rig.push_back(gs);
        labels.push_back(nuclide);
        color_idx++;
    }

    DrawAndSaveGraphs(gm_rig, labels, "Rigidity [GV]", ("#mu_{" + x_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_base, root_file, "_Mean_vs_Rigidity");
    DrawAndSaveGraphs(gs_rig, labels, "Rigidity [GV]", ("#sigma_{" + x_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_base, root_file, "_Sigma_vs_Rigidity");

    for (auto g : gm_rig) delete g;
    for (auto g : gs_rig) delete g;
}

FitResult twoStepGaussianFit(TH1* hist, const string& title_prefix, double default_center_x, const string& rig_label) {
    FitResult result = {0};
    if (!hist || hist->GetEntries() < 10) return result;
    hist->Sumw2();

    // Step 1: Rough fit
    double x1 = 0, x2 = 0;
    findFitRange(hist, 0.80, default_center_x, x1, x2);
    if (x2 <= x1) return result;

    TF1* f1 = new TF1("f1", "gaus", x1, x2);
    f1->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS());
    if (hist->Fit(f1, "QRS") != 0) { delete f1; return result; }
    
    double mean0 = f1->GetParameter(1);
    double sigma0 = f1->GetParameter(2);
    delete f1;

    // Step 2: Fine fit
    double x_min = mean0 - 4.0 * abs(sigma0);
    double x_max = mean0 + 4.0 * abs(sigma0);
    
    // Special handling for Beryllium in low rigidity
    if(rig_label == "30-50 GV" && title_prefix.find("Be") != string::npos){
        x_min = mean0 - 2.5 * abs(sigma0);
        x_max = mean0 + 2.5 * abs(sigma0);
    }
    if (x_max <= x_min) return result;

    TCanvas* c = (TCanvas*)gROOT->FindObject("c_fit"); 
    hist->SetTitle(title_prefix.c_str());
    hist->GetYaxis()->SetTitle("Events");
    double buf = 0.6 * (x_max - x_min);
    hist->GetXaxis()->SetRangeUser(x_min - buf, x_max + buf);

    vector<vector<double>> data = DoGausPlusAsymGausFit(hist, x_min, x_max, c, true);

    if (data.size() >= 2 && data[1][4] > 0) { 
        result = {data[0][0], data[1][0], data[0][1], data[1][1], 
                  data[0][4], data[1][4], data[0][2], data[1][2], data[0][3], data[1][3]};
        
        c->cd();
        TLine *l1 = new TLine(x_min, 0, x_min, hist->GetMaximum()), *l2 = new TLine(x_max, 0, x_max, hist->GetMaximum());
        l1->SetLineStyle(2); l1->SetLineColor(kRed); l1->Draw("same");
        l2->SetLineStyle(2); l2->SetLineColor(kRed); l2->Draw("same");
    } else {
        c->cd(); hist->Draw("hist");
        TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
        lat.DrawLatex(0.6, 0.85, "Fit Failed");
        if (!rig_label.empty()) lat.DrawLatex(0.15, 0.85, ("Rig: " + rig_label).c_str());
    }
    c->Update();
    return result;
}

// --- Main Analysis Logic ---
void Analyze(TFile* file, const string& suffix, const string& output_dir,
             const MCConfig& conf, TFile* out_root, bool is_h5,
             map<string, TH2F*>& all_h2_mean, map<string, TH2F*>& all_h2_sigma,
             map<string, FitResult>& all_h4_results) {

    cout << "\n--- Analyzing " << (is_h5 ? "H5" : "H4") << " Series: " << suffix << " for " << conf.nuclide << " ---" << endl;
    
    HistInfo info = getHistInfo(suffix);
    string base_name = info.output_suffix + suffix.substr(3) + "_" + conf.nuclide;
    string pdf_name = output_dir + "FitResults_" + base_name + ".pdf";
    string x_label = is_h5 ? "#Delta(1/#beta)" : "1/#beta";
    
    TCanvas* c_fit = new TCanvas("c_fit", "Fits", 800, 600);
    c_fit->Print((pdf_name + "[").c_str(), "pdf");

    TH2F *h2_m = nullptr, *h2_s = nullptr;
    if (is_h5) {
        out_root->cd();
        h2_m = new TH2F(("h2_mean_" + suffix + "_" + conf.nuclide).c_str(), (info.title_description + " Mean").c_str(), 
                        N_RIG_BINS, &H5_RIG_BINS_EDGES[0], 1, 1.5, 2.5);
        h2_s = new TH2F(("h2_sigma_" + suffix + "_" + conf.nuclide).c_str(), (info.title_description + " Sigma").c_str(), 
                        N_RIG_BINS, &H5_RIG_BINS_EDGES[0], 1, 1.5, 2.5);
    }

    string full_name = "UnbiasedL1Inner_" + conf.element + "_" + suffix;
    TObject* obj = file->Get(full_name.c_str());

    if (obj) {
        if (is_h5) {
            TH2F* h2 = dynamic_cast<TH2F*>(obj);
            if (h2) {
                for (int i = 0; i < N_RIG_BINS; ++i) {
                    int bin_y_max = (i == N_RIG_BINS - 1) ? h2->GetNbinsY() : h2->GetYaxis()->FindFixBin(H5_RIG_BINS_EDGES[i+1] - 1e-6);
                    TH1D* h1 = h2->ProjectionX(("px_" + full_name + conf.nuclide + to_string(i)).c_str(), h2->GetYaxis()->FindFixBin(H5_RIG_BINS_EDGES[i]), bin_y_max);
                    h1->SetDirectory(nullptr);
                    
                    if (h1->GetEntries() > 0) {
                        AutoRebin(h1);
                        string title = conf.nuclide + " - " + info.output_suffix + " [" + H5_RIG_LABELS[i] + "]";
                        FitResult res = twoStepGaussianFit(h1, title, 0.0, H5_RIG_LABELS[i]);
                        
                        if (res.mean_err > 0 || res.sigma_err > 0) {
                            h2_m->SetBinContent(i+1, 1, res.mean); h2_m->SetBinError(i+1, 1, res.mean_err);
                            h2_s->SetBinContent(i+1, 1, res.sigma); h2_s->SetBinError(i+1, 1, res.sigma_err);
                            c_fit->Print(pdf_name.c_str(), "pdf");
                        }
                    }
                    delete h1;
                }
                all_h2_mean[conf.nuclide] = h2_m;
                all_h2_sigma[conf.nuclide] = h2_s;
            }
        } else {
            TH1F* h1 = dynamic_cast<TH1F*>(obj);
            if (h1) {
                TH1F* h_clone = (TH1F*)h1->Clone(("cl_" + full_name + conf.nuclide).c_str());
                h_clone->SetDirectory(nullptr);
                h_clone->SetMarkerStyle(20); h_clone->SetMarkerSize(1.2);
                AutoRebin(h_clone);
                
                FitResult res = twoStepGaussianFit(h_clone, conf.nuclide + " - " + info.title_description, 1.0, "");
                if (res.mean_err > 0 || res.sigma_err > 0) {
                    all_h4_results[conf.nuclide + "_" + suffix] = res;
                    c_fit->Print(pdf_name.c_str(), "pdf");
                }
                delete h_clone;
            }
        }
    } else {
        cerr << "Warning: " << full_name << " not found in " << conf.filename << endl;
    }

    c_fit->Print((pdf_name + "]").c_str(), "pdf");
    delete c_fit;
}

// --- Entry Point ---
void BetaFit_MC() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetErrorX(0); gStyle->SetOptFit(0);
    
    const string in_dir = "/eos/user/z/zixuan/Isotope/Add/";
    const string out_dir = "/eos/user/z/zixuan/Isotope/Beta/MC/";
    const string out_root = out_dir + "MC_RICHBetaStudy.root";
    
    if (gSystem->AccessPathName(out_dir.c_str())) gSystem->mkdir(out_dir.c_str(), kTRUE);

    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) { cerr << "Error creating output: " << out_root << endl; return; }

    map<string, TH2F*> map_h5_mean, map_h5_sigma;
    map<string, FitResult> map_h4;

    // 1. Process all files defined in MC_DEFINITIONS
    for (const auto& conf : MC_DEFINITIONS) {
        TFile* fin = TFile::Open((in_dir + conf.filename).c_str(), "READ");
        if (!fin || fin->IsZombie()) { cerr << "Error input: " << conf.filename << endl; continue; }
        
        for (const auto& s : {"ID_H5a", "ID_H5b"}) Analyze(fin, s, out_dir, conf, fout, true, map_h5_mean, map_h5_sigma, map_h4);
        for (const auto& s : {"ID_H4a", "ID_H4b"}) Analyze(fin, s, out_dir, conf, fout, false, map_h5_mean, map_h5_sigma, map_h4);
        
        fin->Close(); delete fin;
    }

    // 2. Draw combined results
    fout->cd();
    // Save TH2s first
    for (const auto& pair : map_h5_mean) pair.second->Write();
    for (const auto& pair : map_h5_sigma) pair.second->Write();

    for (const string& s : {"ID_H5a", "ID_H5b"}) {
        HistInfo info = getHistInfo(s);
        DrawMC_H5_Combined(map_h5_mean, map_h5_sigma, out_dir, info, 
                           info.output_suffix + s.substr(3) + "_MC_Combined", "#Delta(1/#beta)", fout);
    }
    
    for (const string& s : {"ID_H4a", "ID_H4b"}) {
        HistInfo info = getHistInfo(s);
        DrawMC_H4_Combined(map_h4, out_dir, info, 
                           info.output_suffix + s.substr(3) + "_MC_Combined", "1/#beta", fout);
    }

    fout->Close(); delete fout;
    gROOT->SetBatch(kFALSE);
    cout << "\nDone. Results in: " << out_dir << endl;
}