#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TError.h>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <numeric>

// 假设 SplineFit 函数在此头文件中定义
#include "../Tool.h"

using namespace std;
using namespace AMS_Iso;

// =================================================================================
// ValueWithError struct for robust error propagation
// =================================================================================
struct ValueWithError {
    double val = 0.0;
    double err = 0.0;

    ValueWithError() = default;
    ValueWithError(double v, double e) : val(v), err(e) {}

    ValueWithError operator+(const ValueWithError& o) const {
        return {val + o.val, std::hypot(err, o.err)};
    }
    ValueWithError operator-(const ValueWithError& o) const {
        return {val - o.val, std::hypot(err, o.err)};
    }
    ValueWithError operator*(const ValueWithError& o) const {
        double new_val = val * o.val;
        if (new_val == 0.0) return {0.0, 0.0};
        double rel_err_sq = (val != 0 ? (err/val)*(err/val) : 0) + (o.val != 0 ? (o.err/o.val)*(o.err/o.val) : 0);
        return {new_val, std::abs(new_val) * std::sqrt(rel_err_sq)};
    }
    ValueWithError operator/(const ValueWithError& o) const {
        if (o.val == 0.0) return {0.0, 0.0}; // Or handle as error
        double new_val = val / o.val;
        if (new_val == 0.0) return {0.0, 0.0};
        double rel_err_sq = (val != 0 ? (err/val)*(err/val) : 0) + (o.val != 0 ? (o.err/o.val)*(o.err/o.val) : 0);
        return {new_val, std::abs(new_val) * std::sqrt(rel_err_sq)};
    }
};

// --- 全局配置 ---
const string FLUX_FILE_PATH = "/eos/user/z/zixuan/Isotope/FluxSmooth/FluxSmooth.root";
const string ACC_FILE_DIR = "/eos/user/z/zixuan/Isotope/Add/";
const string OUTPUT_DIR = "/eos/user/z/zixuan/Isotope/BkgEst/";
std::vector<double> xpoints = {0.33, 0.77, 1.41, 3.16, 5.62, 10.00, 30, 60, 100};
std::vector<double> xpoints2 = {0.33, 1.41, 5.62, 10.00, 30, 60, 100};

const bool DO_REBIN = true;

const vector<string> ALL_PARTICLES = {"Be7", "Be9", "Be10", "B10", "B11", "C12", "N14", "N15","O16"};
const vector<string> ALL_SECONDARY_ISOTOPES = {"Be7", "Be9", "Be10"};

const vector<string> DETECTORS = {"TOF", "NaF", "AGL"};

map<string, pair<double, double>> DETECTOR_RANGES = {
    {"TOF", {0.33, 1.1133}},
    {"NaF", {1.1133, 3.05913}},
    {"AGL", {3.05913, 22.0}}
};

const vector<int> PLOT_COLORS = {kRed, kBlue, kGreen+2, kMagenta, kOrange-3, kCyan+1, kPink+7, kSpring-5};

// --- 工具函数 ---

int getMassNumber(const string& particleName) {
    size_t first_digit = particleName.find_first_of("0123456789");
    if (first_digit == string::npos) {
        throw runtime_error("getMassNumber: Could not find mass number in particle name: " + particleName);
    }
    return stoi(particleName.substr(first_digit));
}

TH1* getHistClone(TFile* file, const string& histName) {
    if (!file || file->IsZombie()) {
        throw runtime_error("getHistClone: Invalid file provided.");
    }
    TH1* h_orig = dynamic_cast<TH1*>(file->Get(histName.c_str()));
    if (!h_orig) {
        throw runtime_error("Histogram not found: " + histName + " in file " + file->GetName());
    }
    TH1* h_clone = (TH1*)h_orig->Clone((histName + "_clone").c_str());
    h_clone->SetDirectory(nullptr);
    
    if (DO_REBIN) {
        h_clone->Rebin(2);
    }
    h_clone->Sumw2();
    return h_clone;
}

// =================================================================================
// Acceptance calculation updated to use Binomial errors.
// =================================================================================
map<int, ValueWithError> calculateAcceptance(const string& accFilePath, const string& eventHistName) {
    unique_ptr<TFile> file(TFile::Open(accFilePath.c_str()));
    if (!file || file->IsZombie()) {
        throw runtime_error("Cannot open acceptance file: " + accFilePath);
    }
    // 获取原始直方图，不使用Sumw2()，因为我们需要原始计数
    // FIX: Use dynamic_cast to safely convert TObject* to TH1*
    unique_ptr<TH1> h_events_raw(dynamic_cast<TH1*>(file->Get(eventHistName.c_str())));
    unique_ptr<TH1> h_mc_flux_raw(dynamic_cast<TH1*>(file->Get("MC_FLUX_H3")));

    if (!h_events_raw || !h_mc_flux_raw) {
        throw runtime_error("Could not get raw histograms from " + accFilePath);
    }
    
    // 克隆并进行rebin
    TH1* h_events = (TH1*)h_events_raw->Clone((eventHistName + "_clone_acc").c_str());
    h_events->SetDirectory(nullptr);
    TH1* h_mc_flux = (TH1*)h_mc_flux_raw->Clone("mc_flux_clone_acc");
    h_mc_flux->SetDirectory(nullptr);

    if (DO_REBIN) {
        h_events->Rebin(2);
        h_mc_flux->Rebin(2);
    }

    if (h_events->GetNbinsX() != h_mc_flux->GetNbinsX()) {
        delete h_events;
        delete h_mc_flux;
        throw runtime_error("Bin mismatch in " + accFilePath);
    }

    map<int, ValueWithError> acceptance_map;
    double scale = TMath::Power(3.9, 2) * TMath::Pi();

    for (int i = 1; i <= h_events->GetNbinsX(); ++i) {
        double k = h_events->GetBinContent(i);   // Number of passed events
        double N = h_mc_flux->GetBinContent(i);  // Total number of generated events

        if (N > 0) {
            double efficiency = k / N;
            // Binomial error for the efficiency
            double efficiency_err = sqrt(efficiency * (1.0 - efficiency) / N);
            
            // The final acceptance is efficiency * scale factor
            double acceptance_val = efficiency * scale;
            double acceptance_err = efficiency_err * scale; // Error also scales
            
            acceptance_map[i] = ValueWithError(acceptance_val, acceptance_err);
        } else {
            // If no events were generated, acceptance is zero with zero error.
            acceptance_map[i] = ValueWithError(0.0, 0.0);
        }
    }
    
    // Clean up the cloned histograms
    delete h_events;
    delete h_mc_flux;

    return acceptance_map;
}

TF1* smoothRatio(TH1* h_ratio, const string& fitName) {
    TF1* fit = nullptr;
    try {
        if(fitName == "fit_acc_ratio_combined_C12_to_Be7" || fitName == "fit_acc_ratio_combined_C12_to_B11" || fitName == "fit_acc_ratio_combined_C12_to_Be10" || fitName == "fit_acc_ratio_combined_N15_to_Be7" || fitName == "fit_acc_ratio_combined_O16_to_Be10" || fitName == "fit_acc_ratio_combined_O16_to_B11"){
            fit = SplineFit(h_ratio, xpoints2.data(), xpoints2.size(), 0x38, "b1e1", fitName.c_str(), 0.3, 20.5);
        }
        else{
            fit = SplineFit(h_ratio, xpoints.data(), xpoints.size(), 0x38, "b1e1", fitName.c_str(), 0.3, 20.5);
        }
    } catch (const std::exception& e) {
        cerr << "SplineFit failed for " << fitName << ": " << e.what() << endl; return nullptr;
    }
    if (!fit) {
        cerr << "SplineFit returned nullptr for " << fitName << endl; return nullptr;
    }
    return fit;
}

// Helper to set standard style for histograms
void setHistStyle(TH1* h, int color) {
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.2);
}


// --- 主执行函数 ---
void est_frag() {
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0); 
    gErrorIgnoreLevel = kWarning;

    unique_ptr<TFile> fluxFile(TFile::Open(FLUX_FILE_PATH.c_str()));
    if (!fluxFile || fluxFile->IsZombie()) {
        cerr << "FATAL: Cannot open flux file: " << FLUX_FILE_PATH << endl;
        return;
    }
    map<string, TF1*> fluxes;
    for (const auto& particle : ALL_PARTICLES) {
        fluxes[particle] = dynamic_cast<TF1*>(fluxFile->Get((particle + "_spline_Ek").c_str()));
        if (!fluxes[particle]) {
            cerr << "FATAL: Flux function for " << particle << " not found!" << endl;
            return;
        }
    }

    map<string, vector<TH1D*>> results_by_secondary;
    map<string, vector<string>> labels_by_secondary;
    
    map<string, bool> debug_printed_for_channel;

    // --- Calculation Loop ---
    for (const auto& secondary : ALL_SECONDARY_ISOTOPES) {
        vector<string> primaries_to_use;
        int secondary_mass = getMassNumber(secondary);
        for (const auto& potential_primary : ALL_PARTICLES) {
            if (getMassNumber(potential_primary) > secondary_mass) {
                primaries_to_use.push_back(potential_primary);
            }
        }

        if (primaries_to_use.empty()) continue;

        for (const auto& primary : primaries_to_use) {
            if(primary == "Be7") continue;
            
            cout << "\n>>> Processing channel: " << primary << " -> " << secondary << "..." << endl;
            string channel_label = primary + " #rightarrow " + secondary;
            
            string pdf_path = OUTPUT_DIR + "rew_results_" + primary + "_to_" + secondary + ".pdf";
            unique_ptr<TCanvas> intermediate_canvas(new TCanvas("intermediate_canvas", "Intermediate Results", 900, 700));
            intermediate_canvas->Print((pdf_path + "[").c_str());
            
            TF1* flux_X = fluxes.at(primary);
            TF1* flux_Y = fluxes.at(secondary);

            // --- PDF Page 1: Fluxes ---
            intermediate_canvas->Clear();
            intermediate_canvas->SetLogy(true);
            flux_X->SetLineColor(kRed);
            flux_Y->SetLineColor(kBlue);
            flux_X->SetTitle(Form("Input Fluxes for %s -> %s;E_{k}/n [GeV/n];Flux", primary.c_str(), secondary.c_str()));
            flux_X->GetXaxis()->SetRangeUser(0.3, 100.0);
            flux_X->Draw();
            flux_Y->Draw("SAME");
            {
                unique_ptr<TLegend> leg(new TLegend(0.78, 0.78, 0.88, 0.88));
                leg->SetBorderSize(0); leg->SetFillStyle(0);
                leg->AddEntry(flux_X, Form("Flux(%s)", primary.c_str()), "l");
                leg->AddEntry(flux_Y, Form("Flux(%s)", secondary.c_str()), "l");
                leg->Draw();
                intermediate_canvas->Print(pdf_path.c_str());
            }
            
            auto flux_ratio_func = [&](double ek) -> ValueWithError {
                double valX = flux_X->Eval(ek);
                double valY = flux_Y->Eval(ek);
                // Assume 1% relative error on flux values as a reasonable estimate
                ValueWithError flux_vw_X(valX, valX * 0.01);
                ValueWithError flux_vw_Y(valY, valY * 0.01);
                return flux_vw_X / flux_vw_Y;
            };
            
            unique_ptr<TH1> h_acc_ratio_combined(nullptr);
            TH1D* h_final_epsilon = nullptr;

            for (const auto& det : DETECTORS) {
                cout << "  -- Detector: " << det << " --" << endl;
                try {
                    string channel_suffix;
                    if (secondary.rfind("Be", 0) == 0) { channel_suffix = "_rew_frag4.root"; } 
                    else if (secondary.rfind("B", 0) == 0) { channel_suffix = "_rew_frag5.root"; } 
                    else { throw runtime_error("Unknown secondary isotope family: " + secondary); }
                    
                    string chargePart, massPart;
                    if (secondary.rfind("Be", 0) == 0) { chargePart = "_Z4_"; massPart = secondary.substr(2); } 
                    else { chargePart = "_Z5_"; massPart = secondary.substr(1); }
                    
                    string eventHistName = "L1Inner_MC_BKG_H3a_" + det + chargePart + "Mass" + massPart;
                    string acc_path_Y    = ACC_FILE_DIR + secondary + channel_suffix;
                    string acc_path_XfragY = ACC_FILE_DIR + primary + channel_suffix;
                    
                    map<int, ValueWithError> acc_map_Y = calculateAcceptance(acc_path_Y, eventHistName);
                    map<int, ValueWithError> acc_map_XfragY = calculateAcceptance(acc_path_XfragY, eventHistName);
                    
                    unique_ptr<TH1> h_temp_for_bins(getHistClone(TFile::Open(acc_path_Y.c_str()), eventHistName));
                    unique_ptr<TH1D> h_acc_ratio_det( (TH1D*)h_temp_for_bins->Clone(Form("h_acc_ratio_%s_%s_%s", primary.c_str(), secondary.c_str(), det.c_str())) );
                    h_acc_ratio_det->Reset();

                    for (int i = 1; i <= h_acc_ratio_det->GetNbinsX(); ++i) {
                        ValueWithError acc_ratio = acc_map_XfragY[i] / acc_map_Y[i];
                        h_acc_ratio_det->SetBinContent(i, acc_ratio.val);
                        h_acc_ratio_det->SetBinError(i, acc_ratio.err);
                    }

                    if (!h_final_epsilon) {
                        h_final_epsilon = (TH1D*)h_acc_ratio_det->Clone(Form("h_epsilon_%s_to_%s", primary.c_str(), secondary.c_str()));
                        h_final_epsilon->SetTitle(channel_label.c_str());
                        h_final_epsilon->Reset();
                    }

                    if (!h_acc_ratio_combined) {
                        h_acc_ratio_combined.reset((TH1*)h_acc_ratio_det->Clone(Form("h_acc_ratio_combined_%s_to_%s", primary.c_str(), secondary.c_str())));
                        h_acc_ratio_combined->Reset();
                        h_acc_ratio_combined->SetTitle("Combined Acceptance Ratio (All Detectors)");
                        setHistStyle(h_acc_ratio_combined.get(), kBlack);
                    }

                    for (int i = 1; i <= h_acc_ratio_det->GetNbinsX(); ++i) {
                        double ek_n = h_acc_ratio_det->GetXaxis()->GetBinCenter(i);
                        if (ek_n >= DETECTOR_RANGES.at(det).first && ek_n < DETECTOR_RANGES.at(det).second) {
                            h_acc_ratio_combined->SetBinContent(i, h_acc_ratio_det->GetBinContent(i));
                            h_acc_ratio_combined->SetBinError(i, h_acc_ratio_det->GetBinError(i));
                        }
                    }

                } catch (const std::runtime_error& e) {
                    cerr << "    ERROR processing " << det << " for " << channel_label << ": " << e.what() << endl;
                }
            } // Detector loop
            
            if (!h_acc_ratio_combined) {
                cerr << "ERROR: Combined acceptance ratio histogram could not be created for " << channel_label << ". Skipping." << endl;
                intermediate_canvas->Print((pdf_path + "]").c_str());
                continue;
            }
            
            intermediate_canvas->Clear();
            intermediate_canvas->SetLogy(false);
            h_acc_ratio_combined->SetTitle(Form("Acceptance Ratio & Smooth (%s);E_{k}/n [GeV/n];Acc(X#rightarrowY)/Acc(Y)", "Combined"));
            h_acc_ratio_combined->GetXaxis()->SetRangeUser(0.3, 20.5);
            h_acc_ratio_combined->Draw("PZ");
            
            unique_ptr<TF1> fit_acc_ratio_combined(smoothRatio(h_acc_ratio_combined.get(), Form("fit_acc_ratio_combined_%s_to_%s", primary.c_str(), secondary.c_str())));
            if (!fit_acc_ratio_combined) {
                cerr << "    WARNING: Global spline fit failed for " << channel_label << ". Skipping final calculation." << endl;
                intermediate_canvas->Print((pdf_path + "]").c_str());
                continue;
            }
            fit_acc_ratio_combined->SetLineColor(kRed);
            fit_acc_ratio_combined->Draw("SAME");
            intermediate_canvas->Print(pdf_path.c_str());
            
            for (int i = 1; i <= h_final_epsilon->GetNbinsX(); ++i) {
                double ek_n = h_final_epsilon->GetBinCenter(i);
                
                string current_det = "";
                for (const auto& det_pair : DETECTOR_RANGES) {
                    if (ek_n >= det_pair.second.first && ek_n < det_pair.second.second) {
                        current_det = det_pair.first;
                        break;
                    }
                }

                if (!current_det.empty()) {
                    ValueWithError flux_ratio_vw = flux_ratio_func(ek_n);
                    ValueWithError acc_ratio_vw(fit_acc_ratio_combined->Eval(ek_n), h_acc_ratio_combined->GetBinError(i));
                    
                    ValueWithError epsilon_vw = flux_ratio_vw * acc_ratio_vw;
                    h_final_epsilon->SetBinContent(i, epsilon_vw.val);
                    h_final_epsilon->SetBinError(i, epsilon_vw.err);
                }
            }

            if (h_final_epsilon) {
                results_by_secondary[secondary].push_back(h_final_epsilon);
                labels_by_secondary[secondary].push_back(channel_label);

                intermediate_canvas->Clear();
                intermediate_canvas->SetLogy(false);
                h_final_epsilon->SetTitle(Form("Final #epsilon for %s;E_{k}/n [GeV/n];#epsilon", channel_label.c_str()));
                setHistStyle(h_final_epsilon, kBlack);
                h_final_epsilon->SetMinimum(0);
                h_final_epsilon->GetXaxis()->SetRangeUser(0.3, 20.5);
                h_final_epsilon->Draw("P");
                intermediate_canvas->Print(pdf_path.c_str());
            }
            
            intermediate_canvas->Print((pdf_path + "]").c_str());
            cout << "    INFO: Intermediate results saved to " << pdf_path << endl;

        } // Primary loop
    } // Secondary loop

    cout << "\n--- All calculations finished. Starting to generate plots. ---\n" << endl;

    unique_ptr<TFile> outFile(TFile::Open((OUTPUT_DIR + "Epsilon_Results.root").c_str(), "RECREATE"));

    for (auto& pair : results_by_secondary) {
        const string& secondary_name = pair.first;
        vector<TH1D*>& hists_to_plot = pair.second;
        const vector<string>& labels = labels_by_secondary.at(secondary_name);

        if (hists_to_plot.empty()) continue;
        
        TH1D* h_sum = (TH1D*)hists_to_plot[0]->Clone(Form("h_epsilon_sum_%s", secondary_name.c_str()));
        h_sum->SetTitle(Form("Total Fragmentation to %s", secondary_name.c_str()));
        h_sum->Reset();

        for (TH1D* h : hists_to_plot) {
            h_sum->Add(h);
        }

        double y_min = 1.0e10, y_max = -1.0e10;
        for (const auto& h : hists_to_plot) {
            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                double content = h->GetBinContent(i);
                if (content > 0) {
                    y_min = std::min(y_min, content);
                    y_max = std::max(y_max, content);
                }
            }
        }
        y_max = std::max(y_max, h_sum->GetMaximum());

        y_min = (y_min < 1.0e9) ? y_min * 0.5 : 1e-6;
        y_max = (y_max > 0) ? y_max * 1.4 : 1.0;

        TCanvas* c1 = new TCanvas(Form("c_%s", secondary_name.c_str()), Form("Epsilon for %s", secondary_name.c_str()), 800, 400);
        c1->SetGrid();
        
        TH1* frame = c1->DrawFrame(0.0, 0, 20.5, y_max);
        frame->SetTitle(Form("Fragmentation Contribution to %s;E_{k}/n [GeV/n];#epsilon (Fragmentation Fraction)", secondary_name.c_str()));
        frame->GetYaxis()->SetTitleOffset(1.2);

        TLegend* legend = new TLegend(0.7, 0.56, 0.95, 0.86);
        legend->SetTextSize(0.04); 
        legend->SetBorderSize(0); legend->SetFillStyle(0);

        for (size_t i = 0; i < hists_to_plot.size(); ++i) {
            TH1D* h = hists_to_plot[i];
            int color = PLOT_COLORS[i % PLOT_COLORS.size()];
            setHistStyle(h, color);
            h->Draw("P SAME");
            legend->AddEntry(h, labels[i].c_str(), "p");
        }
        
        setHistStyle(h_sum, kBlack);
        h_sum->Draw("P SAME");
        legend->AddEntry(h_sum, "Total", "p");

        legend->Draw();
        c1->Update();

        string output_filename = OUTPUT_DIR + "rew_frag_epsilon_" + secondary_name + ".png";
        c1->SaveAs(output_filename.c_str());
        cout << ">>> Successfully created plot: " << output_filename << endl;
        
        outFile->cd();
        for (TH1D* h : hists_to_plot) {
            h->Write();
        }
        h_sum->Write();

        delete c1;
        delete legend;
        delete h_sum;
    }

    outFile->Close();
    cout << "\n--- All plots generated and results saved to " << outFile->GetName() << ". Script finished. ---\n" << endl;
}