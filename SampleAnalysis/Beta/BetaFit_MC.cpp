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
    double mean;
    double mean_err;
    double sigma;
    double sigma_err;
    double chi2;
    double ndf;
    double LR; 
    double LR_err; 
    double RR; 
    double RR_err; 
};

using H4ResultsMap = map<string, FitResult>;
using MCAllH4Results = map<string, FitResult>;

const vector<double> H5_RIG_BINS_EDGES = {30.0, 50.0, 80.0, 120.0, 160.0, 240.0};
const vector<string> H5_RIG_LABELS = {"30-50 GV", "50-80 GV", "80-120 GV", "120-160 GV", "160-240 GV"};
const int N_RIG_BINS = H5_RIG_BINS_EDGES.size() - 1;

const vector<string> H5_CHARGE_LABELS = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
const vector<double> H5_CHARGE_BINS_EDGES = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
const int N_CHARGE_BINS = H5_CHARGE_BINS_EDGES.size() - 1;

string getNuclideName(const string& filename) {
    size_t start = filename.find_last_of('/') == string::npos ? 0 : filename.find_last_of('/') + 1;
    size_t end = filename.find("_rew_frag4.root");
    if (end == string::npos) {
        end = filename.find_last_of('.');
    }
    if (end == string::npos || end <= start) return "Unknown";
    return filename.substr(start, end - start);
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
    
    double required_integral = total_integral * coverage;
    double current_integral = hist->GetBinContent(center_bin);
    int low_bin = center_bin - 1;
    int high_bin = center_bin + 1;
    x_min = hist->GetXaxis()->GetBinLowEdge(center_bin);
    x_max = hist->GetXaxis()->GetBinUpEdge(center_bin);

    while (current_integral < required_integral) {
        bool extended = false;
        double content_low = (low_bin >= 1) ? hist->GetBinContent(low_bin) : 0;
        double content_high = (high_bin <= hist->GetNbinsX()) ? hist->GetBinContent(high_bin) : 0;

        if (low_bin >= 1 && high_bin <= hist->GetNbinsX()) {
            current_integral += content_low + content_high;
            x_min = hist->GetXaxis()->GetBinLowEdge(low_bin);
            x_max = hist->GetXaxis()->GetBinUpEdge(high_bin);
            low_bin--;
            high_bin++;
            extended = true;
        } else if (low_bin >= 1) {
            current_integral += content_low;
            x_min = hist->GetXaxis()->GetBinLowEdge(low_bin);
            low_bin--;
            extended = true;
        } else if (high_bin <= hist->GetNbinsX()) {
            current_integral += content_high;
            x_max = hist->GetXaxis()->GetBinUpEdge(high_bin);
            high_bin++;
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
    y_min = 1e10;
    y_max = -1e10;
    bool found_point = false;
    for (const auto& g : graphs) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            if (x >= x_min && x <= x_max) {
                double y_low = y - g->GetErrorY(i);
                double y_high = y + g->GetErrorY(i);
                y_min = min(y_min, y_low);
                y_max = max(y_max, y_high);
                found_point = true;
            }
        }
    }
    if (found_point) {
        double range = y_max - y_min;
        if (range == 0.0) {
            range = abs(y_min * 0.1);
            if (range == 0) range = 0.01;
        }
        double buffer = range * 0.10; 
        y_min -= buffer;
        y_max += buffer;
    } else {
        y_min = -0.01;
        y_max = 0.01;
    }
}

HistInfo getHistInfo(const string& suffix) {
    HistInfo info;
    if (suffix == "ID_H5a") { info = {"Rigidity [GV]", "Rigidity", "NaF-Tracker #Delta(1/#beta)"}; }
    else if (suffix == "ID_H5b") { info = {"Rigidity [GV]", "Rigidity", "AGL-Tracker #Delta(1/#beta)"}; }
    else if (suffix == "ID_H4a") { info = {"Rigidity [GV]", "Rigidity", "NaF 1/#beta (Rig > 80GV)"}; }
    else if (suffix == "ID_H4b") { info = {"Rigidity [GV]", "Rigidity", "AGL 1/#beta (Rig > 150GV)"}; }
    else { info = {"Y Variable", "YVariable", "Unknown Delta Beta"}; }
    return info;
}

void DrawAndSaveGraphs(vector<TGraphErrors*>& graphs, const vector<string>& legend_labels, 
                         const string& x_axis_title, const string& y_axis_title, 
                         const string& canvas_title_prefix, const string& output_dir, 
                         const string& output_name_base, bool logx, TFile* save_to_root, 
                         const string& graph_name_suffix) {
    
    if (graphs.empty()) return;
    
    gROOT->cd();
    string canvas_name = "c_" + output_name_base + graph_name_suffix;
    TCanvas* c = new TCanvas(canvas_name.c_str(), (canvas_title_prefix + graph_name_suffix).c_str(), 700, 500);
    c->SetGrid();

    double x_min = H5_RIG_BINS_EDGES.front();
    double x_max = H5_RIG_BINS_EDGES.back() * 1.5;
    
    if (graph_name_suffix.find("Charge") != string::npos) {
        x_min = 3.5; 
        x_max = 8.5; 
    } else {
        x_min = graphs[0]->GetXaxis()->GetXmin();
        x_max = graphs[0]->GetXaxis()->GetXmax();
    }
    
    double y_min, y_max;
    getGlobalYRange(graphs, x_min, x_max, y_min, y_max);

    TLegend* leg = new TLegend(0.8, 0.8, 0.99, 0.99); 
    leg->SetFillStyle(0); 
    leg->SetBorderSize(1);

    for (size_t i = 0; i < graphs.size(); ++i) {
        string draw_opt = (i == 0) ? "APZ" : "PZ same";
        
        graphs[i]->Draw(draw_opt.c_str());
        
        if (i == 0) {
            graphs[i]->SetTitle((canvas_title_prefix + " vs " + x_axis_title).c_str());
            graphs[i]->GetXaxis()->SetTitle(x_axis_title.c_str());
            graphs[i]->GetYaxis()->SetTitle(y_axis_title.c_str());
            graphs[i]->GetXaxis()->SetRangeUser(x_min, x_max);
            graphs[i]->GetYaxis()->SetRangeUser(y_min, y_max);
        }

        if (i < legend_labels.size()) {
            leg->AddEntry(graphs[i], legend_labels[i].c_str(), "p");
        }

        if (save_to_root) {
            save_to_root->cd();
            graphs[i]->SetName(("g_" + output_name_base + graph_name_suffix + "_" + to_string(i)).c_str());
            graphs[i]->Write();
        }
    }

    leg->Draw();
    c->SaveAs((output_dir + output_name_base + graph_name_suffix + ".png").c_str());
}

void DrawH4Results_MC(const map<string, FitResult>& all_results, const string& output_dir, const HistInfo& info, 
                      const string& output_name_base, const string& x_axis_label, 
                      const vector<string>& all_nuclides, TFile* save_to_root) {
    
    gROOT->cd(); 
    TGraphErrors* g_mean_z = new TGraphErrors();
    TGraphErrors* g_sigma_z = new TGraphErrors();
    
    g_mean_z->SetMarkerStyle(20); g_sigma_z->SetMarkerStyle(20);
    g_mean_z->SetMarkerSize(1.0); g_sigma_z->SetMarkerSize(1.0);
    g_mean_z->SetMarkerColor(kRed); g_sigma_z->SetMarkerColor(kBlue);
    g_mean_z->SetLineColor(kRed); g_sigma_z->SetLineColor(kBlue);

    for (size_t p = 0; p < all_nuclides.size(); ++p) {
        const string& nuclide = all_nuclides[p];
        
        double z_center = 0.0;
        if (nuclide.size() > 0) {
            char first_char = nuclide[0];
            if (first_char == 'H') z_center = 2.0;
            else if (first_char == 'L') z_center = 3.0;
            else if (first_char == 'B' && nuclide.size() > 1 && (nuclide[1] == 'e' || nuclide[1] == 'E')) z_center = 4.0; 
            else if (first_char == 'B') z_center = 5.0; 
            else if (first_char == 'C') z_center = 6.0;
            else if (first_char == 'N') z_center = 7.0;
            else if (first_char == 'O') z_center = 8.0;
            else continue;
        }

        string suffix = (info.output_suffix == "Rigidity") ? "ID_H4a" : "ID_H4b";
        string key = nuclide + "_" + suffix;

        if (all_results.count(key)) {
            const auto& res = all_results.at(key);
            if (res.mean_err > 0 || res.sigma_err > 0) {
                int n_mean = g_mean_z->GetN();
                g_mean_z->SetPoint(n_mean, z_center, res.mean);
                g_mean_z->SetPointError(n_mean, 0.0, res.mean_err);
                int n_sigma = g_sigma_z->GetN();
                g_sigma_z->SetPoint(n_sigma, z_center, res.sigma);
                g_sigma_z->SetPointError(n_sigma, 0.0, res.sigma_err);
            }
        }
    }
    
    vector<string> single_legend = {"MC"}; 
    
    vector<TGraphErrors*> graphs_mean_z = {g_mean_z};
    DrawAndSaveGraphs(graphs_mean_z, single_legend, "Charge (Z)", ("#mu_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_name_base, false, 
                      save_to_root, "_Mean_vs_Charge");

    vector<TGraphErrors*> graphs_sigma_z = {g_sigma_z};
    DrawAndSaveGraphs(graphs_sigma_z, single_legend, "Charge (Z)", ("#sigma_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_name_base, false, 
                      save_to_root, "_Sigma_vs_Charge");
}

void DrawH5ResultsFromTH2_MC(const map<string, TH2F*>& h2_mean_map, const map<string, TH2F*>& h2_sigma_map,
                             const string& output_dir, const HistInfo& info, const string& output_name_base, 
                             const string& x_axis_label, TFile* save_to_root, const vector<string>& all_nuclides) {
    
    vector<int> colors = {kYellow+2, kMagenta, kBlack, kRed, kBlue, kGreen + 2, kOrange + 1, kViolet, kCyan};
    
    vector<TGraphErrors*> graphs_mean_rig, graphs_sigma_rig;
    vector<string> legend_labels_rig;

    double rig_min = H5_RIG_BINS_EDGES.front();
    double rig_max = H5_RIG_BINS_EDGES.back() * 1.5;

    int color_idx = 0;
    for (const string& nuclide : all_nuclides) {
        if (!h2_mean_map.count(nuclide) || !h2_sigma_map.count(nuclide)) continue;

        TH2F* h2_mean = h2_mean_map.at(nuclide);
        TH2F* h2_sigma = h2_sigma_map.at(nuclide);

        gROOT->cd(); 
        TGraphErrors* g_mean = new TGraphErrors();
        TGraphErrors* g_sigma = new TGraphErrors();
        int color = colors[color_idx % colors.size()];
        
        g_mean->SetMarkerStyle(20 + color_idx); g_sigma->SetMarkerStyle(20 + color_idx);
        g_mean->SetMarkerSize(1.0); g_sigma->SetMarkerSize(1.0);
        g_mean->SetMarkerColor(color); g_sigma->SetMarkerColor(color);
        g_mean->SetLineColor(color); g_sigma->SetLineColor(color);

        int p_bin = 1; 
        for (int j = 1; j <= N_RIG_BINS; ++j) {
            double mean = h2_mean->GetBinContent(j, p_bin);
            double mean_err = h2_mean->GetBinError(j, p_bin);
            double sigma = h2_sigma->GetBinContent(j, p_bin);
            double sigma_err = h2_sigma->GetBinError(j, p_bin);
            
            if (mean_err > 0 || sigma_err > 0) {
                double rig_center = (H5_RIG_BINS_EDGES[j-1] + H5_RIG_BINS_EDGES[j]) / 2.0;
                g_mean->SetPoint(g_mean->GetN(), rig_center, mean);
                g_mean->SetPointError(g_mean->GetN()-1, 0.0, mean_err); 
                g_sigma->SetPoint(g_sigma->GetN(), rig_center, sigma);
                g_sigma->SetPointError(g_sigma->GetN()-1, 0.0, sigma_err); 
            }
        }
        
        g_mean->GetXaxis()->SetRangeUser(rig_min, rig_max);
        g_sigma->GetXaxis()->SetRangeUser(rig_min, rig_max);
        graphs_mean_rig.push_back(g_mean);
        graphs_sigma_rig.push_back(g_sigma);
        legend_labels_rig.push_back(nuclide);
        color_idx++;
    }

    DrawAndSaveGraphs(graphs_mean_rig, legend_labels_rig, "Rigidity [GV]", ("#mu_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_name_base, false, 
                      save_to_root, "_Mean_vs_Rigidity");

    DrawAndSaveGraphs(graphs_sigma_rig, legend_labels_rig, "Rigidity [GV]", ("#sigma_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_name_base, false, 
                      save_to_root, "_Sigma_vs_Rigidity");
}


FitResult twoStepGaussianFit(TH1* hist, const string& x_axis_label, TCanvas* c_fit, const string& title_prefix, double default_center_x, const string& rig_label) {
    FitResult result = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (!hist || hist->GetEntries() < 10) return result;
    hist->Sumw2();

    double x_min_fit1 = 0.0, x_max_fit1 = 0.0;
    findFitRange(hist, 0.80, default_center_x, x_min_fit1, x_max_fit1);
    if (x_max_fit1 <= x_min_fit1) return result;

    gROOT->cd(); 
    TF1* f_gaus1 = new TF1("f_gaus1", "gaus", x_min_fit1, x_max_fit1);
    f_gaus1->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS());
    int fit_status1 = hist->Fit(f_gaus1, "QRS");
    
    if (fit_status1 != 0) {
        return result; 
    }

    double mean0 = f_gaus1->GetParameter(1);
    double sigma0 = f_gaus1->GetParameter(2);

    double x_min_fit2 = mean0 - 4.0 * abs(sigma0);
    double x_max_fit2 = mean0 + 4.0 * abs(sigma0);
    
    if (rig_label == "30-50 GV" && title_prefix.find("Be") != string::npos) {
        x_min_fit2 = mean0 - 2.5 * abs(sigma0);
        x_max_fit2 = mean0 + 2.5 * abs(sigma0);
    }
    
    if (x_max_fit2 <= x_min_fit2) return result;
    
    hist->SetTitle(title_prefix.c_str());
    hist->GetYaxis()->SetTitle("Events");
    double buffer = 0.6 * (x_max_fit2 - x_min_fit2);
    hist->GetXaxis()->SetRangeUser(x_min_fit2 - buffer, x_max_fit2 + buffer);
    
    vector<vector<double>> fit_data = DoGausPlusAsymGausFit(hist, x_min_fit2, x_max_fit2, c_fit, true);
    
    if (fit_data.size() == 2 && fit_data[1][4] > 0) { 
        result.mean = fit_data[0][0];
        result.mean_err = fit_data[1][0];
        result.sigma = fit_data[0][1];
        result.sigma_err = fit_data[1][1];
        result.LR = fit_data[0][2];
        result.LR_err = fit_data[1][2]; 
        result.RR = fit_data[0][3];
        result.RR_err = fit_data[1][3]; 
        result.chi2 = fit_data[0][4];
        result.ndf = fit_data[1][4];
        c_fit->cd();
        
        TLine* l_min = new TLine(x_min_fit2, 0, x_min_fit2, hist->GetMaximum() * 1.);
        TLine* l_max = new TLine(x_max_fit2, 0, x_max_fit2, hist->GetMaximum() * 1.);
        l_min->SetLineStyle(2); l_max->SetLineStyle(2);
        l_min->SetLineColor(kRed); l_max->SetLineColor(kRed);
        l_min->Draw("same"); l_max->Draw("same");

        c_fit->Update();
    } else {
        c_fit->cd();
        hist->Draw("hist");
        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
        latex.DrawLatex(0.6, 0.85, "Fit Failed");
        if (!rig_label.empty()) {
            stringstream ss; ss << "Rig: " << rig_label;
            latex.DrawLatex(0.15, 0.85, ss.str().c_str());
        }
        c_fit->Update();
    }
    
    return result;
}

void Analyze(TFile* file, const string& suffix, const string& output_dir,
             const string& nuclide, TFile* output_root_file, bool is_h5,
             map<string, TH2F*>& all_h2_mean, map<string, TH2F*>& all_h2_sigma,
             map<string, FitResult>& all_h4_results) {

    cout << "\n--- Analyzing " << (is_h5 ? "H5" : "H4") << " Series: " << suffix << " for " << nuclide << " ---" << endl;
    gStyle->SetOptFit(0);
    
    const string hist_particle_name = "Helium"; 
    const string prefix = "UnbiasedL1Inner_";
    const string x_axis_label = is_h5 ? "#Delta(1/#beta)" : "1/#beta";
    HistInfo info = getHistInfo(suffix);
    string output_name_base = info.output_suffix + suffix.substr(3) + "_" + nuclide;
    string fit_pdf_path = output_dir + "FitResults_" + output_name_base + ".pdf";
    
    gROOT->cd(); 
    TCanvas* c_fit = new TCanvas("c_fit", ("Gaussian Fit Results - " + nuclide).c_str(), 800, 600);
    c_fit->SetLogy(0);
    c_fit->Print((fit_pdf_path + "[").c_str(), "pdf");

    TH2F* h2_mean = nullptr;
    TH2F* h2_sigma = nullptr;
    if (is_h5) {
        output_root_file->cd();
        h2_mean = new TH2F(("h2_mean_" + suffix + "_" + nuclide).c_str(), (info.title_description + " Mean;Rigidity [GV];Charge (Z)").c_str(), 
                             N_RIG_BINS, &H5_RIG_BINS_EDGES[0], 1, 1.5, 2.5);
        h2_sigma = new TH2F(("h2_sigma_" + suffix + "_" + nuclide).c_str(), (info.title_description + " Sigma;Rigidity [GV];Charge (Z)").c_str(), 
                              N_RIG_BINS, &H5_RIG_BINS_EDGES[0], 1, 1.5, 2.5);
    }

    string full_name = prefix + hist_particle_name + "_" + suffix;
    TObject* obj = file->Get(full_name.c_str());

    if (!obj) {
        cerr << "Warning: Histogram " << full_name << " not found in " << file->GetName() << endl;
        return;
    }

    if (is_h5) {
        TH2F* h2 = dynamic_cast<TH2F*>(obj);
        if (!h2) { return; }

        for (size_t i = 0; i < N_RIG_BINS; ++i) {
            double rig_min = H5_RIG_BINS_EDGES[i];
            double rig_max = H5_RIG_BINS_EDGES[i+1];
            int rig_bin_index = i + 1;
            string rig_label = H5_RIG_LABELS[i];

            int bin_y_min = h2->GetYaxis()->FindFixBin(rig_min);
            int bin_y_max = (i == N_RIG_BINS - 1) ? h2->GetNbinsY() : h2->GetYaxis()->FindFixBin(rig_max - 1e-6);

            gROOT->cd(); 
            TH1D* h1_proj = h2->ProjectionX(("proj_" + full_name + "_" + nuclide + "_" + to_string(i)).c_str(), bin_y_min, bin_y_max);
            h1_proj->SetDirectory(nullptr); 

            if (h1_proj->GetEntries() == 0) { continue; }
            
            for(int r=1 ; r<=10; ++r) { if (h1_proj->GetMaximum() >= 60 || h1_proj->GetNbinsX() < 20) break; h1_proj->Rebin(2); }
            
            string title = nuclide + " - " + info.output_suffix + " [" + rig_label + "]";
            FitResult fit_res = twoStepGaussianFit(h1_proj, x_axis_label, c_fit, title, 0.0, rig_label);

            if (fit_res.mean_err > 0 || fit_res.sigma_err > 0) {
                h2_mean->SetBinContent(rig_bin_index, 1, fit_res.mean);
                h2_mean->SetBinError(rig_bin_index, 1, fit_res.mean_err);
                h2_sigma->SetBinContent(rig_bin_index, 1, fit_res.sigma);
                h2_sigma->SetBinError(rig_bin_index, 1, fit_res.sigma_err);
                c_fit->Print(fit_pdf_path.c_str(), "pdf");
            }
        }
        
        all_h2_mean[nuclide] = h2_mean;
        all_h2_sigma[nuclide] = h2_sigma;

    } else {
        TH1F* h1 = dynamic_cast<TH1F*>(obj);
        if (!h1) { return; }
        
        gROOT->cd(); 
        TH1F* h1_clone = (TH1F*)h1->Clone(("clone_" + full_name + "_" + nuclide).c_str());
        h1_clone->SetDirectory(nullptr);
        h1_clone->SetMarkerStyle(20);
        h1_clone->SetMarkerSize(1.2);

        for(int r=1 ; r<=10; ++r) { if (h1_clone->GetMaximum() >= 60 || h1_clone->GetNbinsX() < 20) break; h1_clone->Rebin(2); }
        
        string title = nuclide + " - " + info.title_description;
        FitResult fit_res = twoStepGaussianFit(h1_clone, x_axis_label, c_fit, title, 1.0, "");

        if (fit_res.mean_err > 0 || fit_res.sigma_err > 0) {
            all_h4_results[nuclide + "_" + suffix] = fit_res;
            c_fit->Print(fit_pdf_path.c_str(), "pdf");
        }
    }
    
    c_fit->Print((fit_pdf_path + "]").c_str(), "pdf");
}

void BetaFit_MC() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetErrorX(0); 

    vector<string> mc_files = {
        "B10_rew_frag4.root", "B11_rew_frag4.root", "Be10_rew_frag4.root", "Be7_rew_frag4.root",
        "Be9_rew_frag4.root", "C12_rew_frag4.root", "N15_rew_frag4.root", "O16_rew_frag4.root"
    };
    
    const string input_dir = "/eos/user/z/zixuan/Isotope/Add/";
    const string output_dir = "/eos/user/z/zixuan/Isotope/Beta/MC/";
    const string output_root_path = output_dir + "MC_RICHBetaStudy.root";
    
    vector<string> all_nuclides;
    for (const string& filename : mc_files) {
        all_nuclides.push_back(getNuclideName(filename));
    }

    if (gSystem->AccessPathName(output_dir.c_str())) {
        cout << "Creating output directory: " << output_dir << endl;
        gSystem->mkdir(output_dir.c_str(), kTRUE);
    }

    TFile* output_root_file = TFile::Open(output_root_path.c_str(), "RECREATE");
    if (!output_root_file || output_root_file->IsZombie()) {
        cerr << "ERROR: Cannot create output ROOT file: " << output_root_path << endl;
        return;
    }
    
    map<string, TH2F*> mc_h5_mean_map;
    map<string, TH2F*> mc_h5_sigma_map;
    map<string, FitResult> mc_h4_results;

    for (const string& filename : mc_files) {
        string full_path = input_dir + filename;
        string nuclide = getNuclideName(filename);
        
        if (nuclide == "Unknown") {
            cerr << "Warning: Skipping unknown file format: " << filename << endl;
            continue;
        }

        TFile* file = TFile::Open(full_path.c_str(), "READ");
        if (!file || file->IsZombie()) {
            cerr << "ERROR: Cannot open input file: " << full_path << endl;
            continue;
        }
        
        vector<string> h5_suffixes = {"ID_H5a", "ID_H5b"};
        for (const string& suffix : h5_suffixes) {
            Analyze(file, suffix, output_dir, nuclide, output_root_file, true, 
                    mc_h5_mean_map, mc_h5_sigma_map, mc_h4_results);
        }

        vector<string> h4_suffixes = {"ID_H4a", "ID_H4b"};
        for (const string& suffix : h4_suffixes) {
            Analyze(file, suffix, output_dir, nuclide, output_root_file, false, 
                    mc_h5_mean_map, mc_h5_sigma_map, mc_h4_results);
        }

        file->Close();
    }

    output_root_file->cd();

    vector<string> h5_suffixes = {"ID_H5a", "ID_H5b"};
    for (const string& suffix : h5_suffixes) {
        HistInfo info = getHistInfo(suffix);
        string output_name_base = info.output_suffix + suffix.substr(3) + "_MC_Combined";
        string x_axis_label = "#Delta(1/#beta)";
        
        map<string, TH2F*> current_mean_map, current_sigma_map;
        for (const string& nuclide : all_nuclides) {
            if (mc_h5_mean_map.count(nuclide)) {
                current_mean_map[nuclide] = mc_h5_mean_map[nuclide];
                current_sigma_map[nuclide] = mc_h5_sigma_map[nuclide];
            }
        }

        for(const auto& pair : current_mean_map) pair.second->Write();
        for(const auto& pair : current_sigma_map) pair.second->Write();

        DrawH5ResultsFromTH2_MC(current_mean_map, current_sigma_map, output_dir, info, 
                                 output_name_base, x_axis_label, output_root_file, all_nuclides);
    }
    
    vector<string> h4_suffixes = {"ID_H4a", "ID_H4b"};
    for (const string& suffix : h4_suffixes) {
        HistInfo info = getHistInfo(suffix);
        string output_name_base = info.output_suffix + suffix.substr(3) + "_MC_Combined";
        string x_axis_label = "1/#beta";

        map<string, FitResult> current_h4_results;
        for (const string& nuclide : all_nuclides) {
            string key = nuclide + "_" + suffix;
            if (mc_h4_results.count(key)) {
                current_h4_results[nuclide] = mc_h4_results[key];
            }
        }

        DrawH4Results_MC(current_h4_results, output_dir, info, 
                          output_name_base, x_axis_label, all_nuclides, output_root_file);
    }

    output_root_file->Close();

    gROOT->SetBatch(kFALSE);
    cout << "\nAnalysis finished. Results saved to: " << output_dir << endl;
}