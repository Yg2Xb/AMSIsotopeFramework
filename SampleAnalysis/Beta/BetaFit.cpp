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

using H4ResultsMap = map<string, FitResult>;

const vector<double> H5_RIG_BINS_EDGES = {30.0, 50.0, 80.0, 120.0, 160.0, 240.0};
const vector<string> H5_RIG_LABELS = {"30-50 GV", "50-80 GV", "80-120 GV", "120-160 GV", "160-240 GV"};
const int N_RIG_BINS = H5_RIG_BINS_EDGES.size() - 1;

const vector<string> H5_CHARGE_LABELS = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
const vector<double> H5_CHARGE_BINS_EDGES = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
const int N_CHARGE_BINS = H5_CHARGE_BINS_EDGES.size() - 1;

const vector<int> COLORS_Z = {kYellow+2, kMagenta, kBlack, kRed, kBlue, kGreen + 2, kOrange + 1};
const vector<int> COLORS_RIG = {kBlack, kRed, kBlue, kGreen+2, kMagenta};


void SetGraphStyle(TGraphErrors* g, int color) {
    g->SetMarkerStyle(20); 
    g->SetMarkerSize(1.0);
    g->SetLineColor(color); 
    g->SetMarkerColor(color);
}

void AutoRebin(TH1* h, double threshold = 80.0) {
    for(int r = 1; r <= 10; ++r) {
        if (h->GetMaximum() >= threshold && h->GetBinWidth(1) > 1e-6) break;
        h->Rebin(2);
    }
}

vector<double> GetBinCenters(const vector<double>& edges) {
    vector<double> centers;
    for (size_t i = 0; i < edges.size() - 1; ++i) 
        centers.push_back((edges[i] + edges[i+1]) / 2.0);
    return centers;
}

void findFitRange(TH1* hist, double coverage, double center_x, double& x_min, double& x_max) {
    if (!hist || hist->GetEntries() == 0) return;
    double total_integral = hist->Integral();
    if (total_integral <= 0) return;
    
    int center_bin = hist->GetXaxis()->FindFixBin(center_x);
    if (center_bin < 1 || center_bin > hist->GetNbinsX()) {
        center_bin = hist->GetMaximumBin(); 
        center_x = hist->GetXaxis()->GetBinCenter(center_bin);
    }
    
    double required = total_integral * coverage;
    double current = hist->GetBinContent(center_bin);
    int low = center_bin - 1;
    int high = center_bin + 1;
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

void DrawAndSaveGraphs(vector<TGraphErrors*>& graphs, const vector<string>& legend_labels,
                       const string& x_title, const string& y_title, const string& title_prefix, 
                       const string& output_dir, const string& output_base, 
                       TFile* root_file, const string& suffix) {
    if (graphs.empty()) return;

    string canvas_name = "c_" + output_base + suffix;
    TCanvas* c = new TCanvas(canvas_name.c_str(), (title_prefix + suffix).c_str(), 700, 500);
    c->SetGrid();

    double x_min = graphs[0]->GetXaxis()->GetXmin();
    double x_max = graphs[0]->GetXaxis()->GetXmax();
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

void DrawH4Results(const H4ResultsMap& results, const string& output_dir, const HistInfo& info,
                   const string& output_base, const string& x_label, 
                   const vector<string>& particles, TFile* root_file) {
    
    TGraphErrors* g_mean = new TGraphErrors();
    TGraphErrors* g_sigma = new TGraphErrors();
    SetGraphStyle(g_mean, kBlack);
    SetGraphStyle(g_sigma, kBlack);

    vector<double> z_centers = GetBinCenters(H5_CHARGE_BINS_EDGES);

    for (size_t p = 0; p < particles.size(); ++p) {
        if (p >= z_centers.size()) break;
        if (results.count(particles[p])) {
            const auto& res = results.at(particles[p]);
            if (res.mean_err > 0 || res.sigma_err > 0) {
                int n = g_mean->GetN();
                g_mean->SetPoint(n, z_centers[p], res.mean);
                g_mean->SetPointError(n, 0.0, res.mean_err);
                g_sigma->SetPoint(n, z_centers[p], res.sigma);
                g_sigma->SetPointError(n, 0.0, res.sigma_err);
            }
        }
    }

    double min_z = H5_CHARGE_BINS_EDGES.front();
    double max_z = H5_CHARGE_BINS_EDGES.back();
    g_mean->GetXaxis()->SetRangeUser(min_z, max_z);
    g_sigma->GetXaxis()->SetRangeUser(min_z, max_z);

    vector<string> legend = {"ISS"};
    vector<TGraphErrors*> v_mean = {g_mean}, v_sigma = {g_sigma};

    DrawAndSaveGraphs(v_mean, legend, "Charge (Z)", ("#mu_{" + x_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_base, root_file, "_Mean_vs_Charge");
    
    DrawAndSaveGraphs(v_sigma, legend, "Charge (Z)", ("#sigma_{" + x_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_base, root_file, "_Sigma_vs_Charge");

    delete g_mean; delete g_sigma;
}

void SliceAndExtractGraphs(TH2F* h2_mean, TH2F* h2_sigma, 
                          bool slice_along_charge, // true: X=Rigidity(fix charge), false: X=Charge(fix rig)
                          vector<TGraphErrors*>& g_means, vector<TGraphErrors*>& g_sigmas) {
    
    int n_outer = slice_along_charge ? N_CHARGE_BINS : N_RIG_BINS;
    int n_inner = slice_along_charge ? N_RIG_BINS : N_CHARGE_BINS;
    const vector<int>& colors = slice_along_charge ? COLORS_Z : COLORS_RIG;
    
    vector<double> x_centers = GetBinCenters(slice_along_charge ? H5_RIG_BINS_EDGES : H5_CHARGE_BINS_EDGES);
    double x_min = (slice_along_charge ? H5_RIG_BINS_EDGES : H5_CHARGE_BINS_EDGES).front();
    double x_max = (slice_along_charge ? H5_RIG_BINS_EDGES : H5_CHARGE_BINS_EDGES).back() * (slice_along_charge ? 1.5 : 1.0);

    for (int i = 0; i < n_outer; ++i) {
        TGraphErrors* gm = new TGraphErrors();
        TGraphErrors* gs = new TGraphErrors();
        int color = colors[i % colors.size()];
        SetGraphStyle(gm, color);
        SetGraphStyle(gs, color);

        for (int j = 0; j < n_inner; ++j) {
            int bin_x = slice_along_charge ? j + 1 : i + 1;
            int bin_y = slice_along_charge ? i + 1 : j + 1;

            double m = h2_mean->GetBinContent(bin_x, bin_y);
            double me = h2_mean->GetBinError(bin_x, bin_y);
            double s = h2_sigma->GetBinContent(bin_x, bin_y);
            double se = h2_sigma->GetBinError(bin_x, bin_y);

            if (me > 0 || se > 0) {
                int n = gm->GetN();
                gm->SetPoint(n, x_centers[j], m);
                gm->SetPointError(n, 0.0, me);
                gs->SetPoint(n, x_centers[j], s);
                gs->SetPointError(n, 0.0, se);
            }
        }
        gm->GetXaxis()->SetRangeUser(x_min, x_max);
        gs->GetXaxis()->SetRangeUser(x_min, x_max);
        g_means.push_back(gm);
        g_sigmas.push_back(gs);
    }
}

void DrawH5ResultsFromTH2(TH2F* h2_mean, TH2F* h2_sigma, const string& output_dir, const HistInfo& info,
                          const string& output_base, const string& x_label, TFile* root_file) {
    
    // 1. Plot vs Rigidity 
    vector<TGraphErrors*> gm_rig, gs_rig;
    SliceAndExtractGraphs(h2_mean, h2_sigma, true, gm_rig, gs_rig);
    
    DrawAndSaveGraphs(gm_rig, H5_CHARGE_LABELS, "Rigidity [GV]", ("#mu_{" + x_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_base, root_file, "_Mean_vs_Rigidity");
    DrawAndSaveGraphs(gs_rig, H5_CHARGE_LABELS, "Rigidity [GV]", ("#sigma_{" + x_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_base, root_file, "_Sigma_vs_Rigidity");

    // 2. Plot vs Charge 
    vector<TGraphErrors*> gm_z, gs_z;
    SliceAndExtractGraphs(h2_mean, h2_sigma, false, gm_z, gs_z);

    DrawAndSaveGraphs(gm_z, H5_RIG_LABELS, "Charge (Z)", ("#mu_{" + x_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_base, root_file, "_Mean_vs_Charge");
    DrawAndSaveGraphs(gs_z, H5_RIG_LABELS, "Charge (Z)", ("#sigma_{" + x_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_base, root_file, "_Sigma_vs_Charge");

    for (auto g : gm_rig) delete g; for (auto g : gs_rig) delete g;
    for (auto g : gm_z) delete g;   for (auto g : gs_z) delete g;
}

FitResult twoStepGaussianFit(TH1* hist, const string& title_prefix, double default_center_x, const string& rig_label) {
    FitResult result = {0};
    if (!hist || hist->GetEntries() < 10) return result;

    // Step 1
    double x1 = 0, x2 = 0;
    findFitRange(hist, 0.80, default_center_x, x1, x2);
    if (x2 <= x1) return result;

    TF1* f1 = new TF1("f1", "gaus", x1, x2);
    f1->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS());
    if (hist->Fit(f1, "QRS") != 0) { delete f1; return result; }
    
    double mean0 = f1->GetParameter(1);
    double sigma0 = f1->GetParameter(2);
    delete f1;

    // Step 2
    double x_min = mean0 - 4.0 * abs(sigma0);
    double x_max = mean0 + 4.0 * abs(sigma0);
    
    if(rig_label == "30-50 GV" && title_prefix.find("Beryllium") != string::npos){
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

    if (data[1][4] > 0) { // Success
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

void Analyze(TFile* file, const string& suffix, const string& output_dir,
             const vector<string>& particles, TFile* out_root, bool is_h5) {

    cout << "\n--- Analyzing " << (is_h5 ? "H5" : "H4") << " Series: " << suffix << " ---" << endl;
    
    HistInfo info = getHistInfo(suffix);
    string base_name = info.output_suffix + suffix.substr(3);
    string pdf_name = output_dir + "FitResults_" + base_name + ".pdf";
    string x_label = is_h5 ? "#Delta(1/#beta)" : "1/#beta";
    
    TCanvas* c_fit = new TCanvas("c_fit", "Fits", 800, 600);
    c_fit->Print((pdf_name + "[").c_str(), "pdf");

    H4ResultsMap h4_res;
    TH2F *h2_m = nullptr, *h2_s = nullptr;
    
    if (is_h5) {
        h2_m = new TH2F(("h2_mean_" + suffix).c_str(), (info.title_description + " Mean").c_str(), 
                        N_RIG_BINS, &H5_RIG_BINS_EDGES[0], N_CHARGE_BINS, &H5_CHARGE_BINS_EDGES[0]);
        h2_s = new TH2F(("h2_sigma_" + suffix).c_str(), (info.title_description + " Sigma").c_str(), 
                        N_RIG_BINS, &H5_RIG_BINS_EDGES[0], N_CHARGE_BINS, &H5_CHARGE_BINS_EDGES[0]);
    }

    for (size_t p = 0; p < particles.size(); ++p) {
        string name = "UnbiasedL1Inner_" + particles[p] + "_" + suffix;
        TObject* obj = file->Get(name.c_str());
        if (!obj) continue;

        if (is_h5) { 
            TH2F* h2 = dynamic_cast<TH2F*>(obj);
            if (!h2) continue;
            for (int i = 0; i < N_RIG_BINS; ++i) {
                int bin_y_max = (i == N_RIG_BINS - 1) ? h2->GetNbinsY() : h2->GetYaxis()->FindFixBin(H5_RIG_BINS_EDGES[i+1] - 1e-6);
                TH1D* h1 = h2->ProjectionX(("px_" + name + to_string(i)).c_str(), h2->GetYaxis()->FindFixBin(H5_RIG_BINS_EDGES[i]), bin_y_max);
                
                if (h1->GetEntries() > 0) {
                    AutoRebin(h1);
                    string title = particles[p] + " - " + info.output_suffix + " [" + H5_RIG_LABELS[i] + "]";
                    FitResult res = twoStepGaussianFit(h1, title, 0.0, H5_RIG_LABELS[i]);
                    
                    if (res.mean_err > 0 || res.sigma_err > 0) {
                        h2_m->SetBinContent(i+1, p+1, res.mean); h2_m->SetBinError(i+1, p+1, res.mean_err);
                        h2_s->SetBinContent(i+1, p+1, res.sigma); h2_s->SetBinError(i+1, p+1, res.sigma_err);
                        c_fit->Print(pdf_name.c_str(), "pdf");
                    }
                }
                delete h1;
            }
        } else { 
            TH1F* h1 = dynamic_cast<TH1F*>(obj);
            if (!h1) continue;
            TH1F* h_clone = (TH1F*)h1->Clone(("cl_" + name).c_str());
            h_clone->SetMarkerStyle(20); h_clone->SetMarkerSize(1.2);
            AutoRebin(h_clone);
            
            FitResult res = twoStepGaussianFit(h_clone, particles[p] + " - " + info.title_description, 1.0, "");
            if (res.mean_err > 0 || res.sigma_err > 0) {
                h4_res[particles[p]] = res;
                c_fit->Print(pdf_name.c_str(), "pdf");
            }
            delete h_clone;
        }
    }

    c_fit->Print((pdf_name + "]").c_str(), "pdf");
    delete c_fit;

    out_root->cd();
    if (is_h5) {
        h2_m->Write(); h2_s->Write();
        DrawH5ResultsFromTH2(h2_m, h2_s, output_dir, info, base_name, x_label, out_root);
        delete h2_m; delete h2_s;
    } else {
        DrawH4Results(h4_res, output_dir, info, base_name, x_label, particles, out_root);
    }
}

void BetaFit() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetErrorX(0); gStyle->SetOptFit(0);

    const string in_path = "/eos/user/z/zixuan/Isotope/Add/Be_frag4.root";
    const string out_dir = "/eos/user/z/zixuan/Isotope/Beta/ISS/";
    const string out_root = out_dir + "ISS_RICHBetaStudy.root";
    
    if (gSystem->AccessPathName(out_dir.c_str())) gSystem->mkdir(out_dir.c_str(), kTRUE);

    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie()) { cerr << "Error opening input: " << in_path << endl; return; }
    
    TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) { cerr << "Error creating output: " << out_root << endl; fin->Close(); return; }

    vector<string> parts = {"Helium","Lithium","Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
    
    for (const auto& s : {"ID_H5a", "ID_H5b"}) Analyze(fin, s, out_dir, parts, fout, true);
    for (const auto& s : {"ID_H4a", "ID_H4b"}) Analyze(fin, s, out_dir, parts, fout, false);

    fout->Close(); fin->Close();
    delete fout; delete fin;
    gROOT->SetBatch(kFALSE);
    cout << "\nDone. Results in: " << out_dir << endl;
}