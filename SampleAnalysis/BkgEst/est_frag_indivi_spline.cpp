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
#include <algorithm> // For std::min/max
#include <iomanip>   // For std::setprecision

// 假设 SplineFit 函数在此头文件中定义
#include "../Tool.h"

using namespace std;
using namespace AMS_Iso; // 假设 SplineFit 在此命名空间

// --- 全局配置 ---
const string FLUX_FILE_PATH = "/eos/user/z/zixuan/Isotope/FluxSmooth/FluxSmooth.root";
const string ACC_FILE_DIR = "/eos/user/z/zixuan/Isotope/Add/";
const string OUTPUT_DIR = "/eos/user/z/zixuan/Isotope/BkgEst/";
std::vector<double> xpoints = {0.26, 0.5, 0.7, 0.9, 1.41, 2.00, 3.16, 5.62, 10.00, 30, 60, 100};
// =================================================================================
// MODIFICATION 1: Add a global boolean to control rebinning.
// =================================================================================
const bool DO_REBIN = true;

const vector<string> ALL_PARTICLES = {"C12", "N14", "O16", "B10", "B11", "Be7", "Be9", "Be10"};
const vector<string> ALL_SECONDARY_ISOTOPES = {"B10", "B11", "Be7", "Be9", "Be10"};

const vector<string> DETECTORS = {"TOF", "NaF", "AGL"};

map<string, pair<double, double>> DETECTOR_RANGES = {
    {"TOF", {0.33, 1.1133}},
    {"NaF", {1.1133, 3.05913}},
    {"AGL", {3.05913, 20.0}}
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
    
    // Apply rebinning based on the global flag
    if (DO_REBIN) {
        h_clone->Rebin(2);
    }

    return h_clone;
}

TH1* calculateAcceptance(const string& accFilePath, const string& eventHistName) {
    unique_ptr<TFile> file(TFile::Open(accFilePath.c_str()));
    if (!file || file->IsZombie()) {
        throw runtime_error("Cannot open acceptance file: " + accFilePath);
    }
    unique_ptr<TH1> h_events(getHistClone(file.get(), eventHistName));
    unique_ptr<TH1> h_mc_flux(getHistClone(file.get(), "MC_FLUX_H3"));
    if (h_events->GetNbinsX() != h_mc_flux->GetNbinsX()) {
        throw runtime_error("Bin mismatch in " + accFilePath);
    }
    h_events->Divide(h_mc_flux.get());
    h_events->Scale(TMath::Power(3.9, 2) * TMath::Pi());
    return h_events.release();
}

// =================================================================================
// MODIFICATION 2: Change knot selection to be based on point count instead of range.
// =================================================================================
TF1* smoothRatio(TH1* h_ratio, const string& detector, const string& fitName) {
    vector<double> knots;
    double fitStart = 0.35;
    double fitEnd = 100;
    int segments = 8;

    // Collect x-values of all points within the fitting range
    vector<double> valid_points_x;
    for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
        double bin_center = h_ratio->GetXaxis()->GetBinCenter(i);
        if (bin_center >= fitStart && bin_center <= fitEnd) {
            valid_points_x.push_back(bin_center);
        }
    }

    if (valid_points_x.size() < 2) {
        cerr << "Not enough points (" << valid_points_x.size() << ") in range [" << fitStart << ", " << fitEnd << "] to perform spline fit for " << fitName << endl;
        return nullptr;
    }

    // Distribute knots based on the index of points, not their x-value range
    for (int i = 0; i <= segments; ++i) {
        // Calculate the index of the point to use for the knot
        int point_index = static_cast<int>(i * (valid_points_x.size() - 1.0) / segments);
        knots.push_back(valid_points_x[point_index]);
    }
    
    // Ensure knots are unique and sorted, as SplineFit might require this
    sort(knots.begin(), knots.end());
    knots.erase(unique(knots.begin(), knots.end()), knots.end());
    if (knots.size() < 2) {
         cerr << "Not enough unique knots generated for " << fitName << endl;
         return nullptr;
    }

    TF1* fit = nullptr;
    try {
        fit = SplineFit(h_ratio, xpoints.data(), xpoints.size(), 0x38, "b2e2", fitName.c_str(), 0.3, 100.0);
    } catch (const std::exception& e) {
        cerr << "SplineFit failed for " << fitName << ": " << e.what() << endl; return nullptr;
    }
    if (!fit) {
        cerr << "SplineFit returned nullptr for " << fitName << endl; return nullptr;
    }
    return fit;
}


// --- 主执行函数 ---
void est_frag() {
    gStyle->SetOptStat(0);
    // Set ErrorX to 0 to prevent ROOT from drawing horizontal error bars with "HIST" option
    gStyle->SetErrorX(0); 
    gErrorIgnoreLevel = kWarning;

    unique_ptr<TFile> fluxFile(TFile::Open(FLUX_FILE_PATH.c_str()));
    if (!fluxFile || fluxFile->IsZombie()) {
        cerr << "FATAL: Cannot open flux file: " << FLUX_FILE_PATH << endl;
        return;
    }
    map<string, TF1*> fluxes;
    for (const auto& particle : ALL_PARTICLES) {
        fluxes[particle] = dynamic_cast<TF1*>(fluxFile->Get((particle + "_Ek").c_str()));
        if (!fluxes[particle]) {
            cerr << "FATAL: Flux function for " << particle << " not found!" << endl;
            return;
        }
    }

    map<string, vector<TH1*>> results_by_secondary;
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

        if (primaries_to_use.empty()) {
            cout << "No heavier primary particles found for " << secondary << ", skipping." << endl;
            continue;
        }

        for (const auto& primary : primaries_to_use) {
            cout << "\n>>> Processing channel: " << primary << " -> " << secondary << "..." << endl;
            string channel_label = primary + " #rightarrow " + secondary;
            
            string pdf_path = OUTPUT_DIR + "intermediate_results_" + primary + "_to_" + secondary + ".pdf";
            unique_ptr<TCanvas> intermediate_canvas(new TCanvas("intermediate_canvas", "Intermediate Results", 900, 700));
            intermediate_canvas->SaveAs((pdf_path + "[").c_str());
            
            TF1* flux_X = fluxes.at(primary);
            TF1* flux_Y = fluxes.at(secondary);

            // --- PDF Page 1: Fluxes ---
            intermediate_canvas->Clear();
            intermediate_canvas->SetLogx(true);
            intermediate_canvas->SetLogy(true); // Only this plot has LogY
            flux_X->SetLineColor(kRed);
            flux_Y->SetLineColor(kBlue);
            flux_X->SetTitle(Form("Input Fluxes for %s -> %s;E_{k}/n [GeV/n];Flux", primary.c_str(), secondary.c_str()));
            flux_X->Draw();
            flux_Y->Draw("SAME");
            {
                unique_ptr<TLegend> leg(new TLegend(0.78, 0.78, 0.88, 0.88));
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(flux_X, Form("Flux(%s)", primary.c_str()), "l");
                leg->AddEntry(flux_Y, Form("Flux(%s)", secondary.c_str()), "l");
                //leg->Draw();
                intermediate_canvas->SaveAs(pdf_path.c_str());
            }

            // =================================================================================
            // MODIFICATION 3: Create a single TF1 for the flux ratio.
            // =================================================================================
            auto flux_ratio_func = [&](double *x, double *p) {
                double fy = flux_Y->Eval(x[0]);
                return (fy > 0) ? flux_X->Eval(x[0]) / fy : 0.0;
            };
            unique_ptr<TF1> tf1_flux_ratio(new TF1(Form("tf1_flux_ratio_%s_to_%s", primary.c_str(), secondary.c_str()), flux_ratio_func, 0.3, 100.0, 0));
            tf1_flux_ratio->SetTitle(Form("Flux Ratio for %s -> %s;E_{k}/n [GeV/n];Flux(%s) / Flux(%s)", primary.c_str(), secondary.c_str(), primary.c_str(), secondary.c_str()));
            tf1_flux_ratio->SetLineColor(kGreen+2);
            
            // --- PDF Page 2: Flux Ratio TF1 ---
            intermediate_canvas->Clear();
            intermediate_canvas->SetLogx(true);
            intermediate_canvas->SetLogy(false);
            tf1_flux_ratio->Draw();
            intermediate_canvas->SaveAs(pdf_path.c_str());

            TH1* h_final_epsilon = nullptr;

            for (const auto& det : DETECTORS) {
                cout << "  -- Detector: " << det << " --" << endl;
                try {
                    string channel_suffix;
                    if (secondary.rfind("Be", 0) == 0) { channel_suffix = "_all_w1_frag4.root"; } 
                    else if (secondary.rfind("B", 0) == 0) { channel_suffix = "_all_w1_frag5.root"; } 
                    else { throw runtime_error("Unknown secondary isotope family: " + secondary); }
                    
                    string chargePart, massPart;
                    if (secondary.rfind("Be", 0) == 0) { chargePart = "_Z4_"; massPart = secondary.substr(2); } 
                    else { chargePart = "_Z5_"; massPart = secondary.substr(1); }
                    
                    string eventHistName = "UnbiasedL1Inner_MC_BKG_H3a_" + det + chargePart + "Mass" + massPart;
                    string acc_path_Y    = ACC_FILE_DIR + secondary + channel_suffix;
                    string acc_path_XfragY = ACC_FILE_DIR + primary + channel_suffix;
                    
                    // Plotting raw histograms for acceptance calculation
                    {
                        unique_ptr<TFile> file_Y(TFile::Open(acc_path_Y.c_str()));
                        if (file_Y && !file_Y->IsZombie()) {
                            unique_ptr<TH1> h_events_Y_raw(getHistClone(file_Y.get(), eventHistName));
                            unique_ptr<TH1> h_mc_flux_Y_raw(getHistClone(file_Y.get(), "MC_FLUX_H3"));
                            
                            intermediate_canvas->Clear();
                            intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                            h_events_Y_raw->SetTitle(Form("Event Count for %s (%s);E_{k}/n [GeV/n];Counts", secondary.c_str(), det.c_str()));
                            h_events_Y_raw->SetMarkerColor(kBlue); h_events_Y_raw->SetLineColor(kBlue);
                            for(int i=0; i<=h_events_Y_raw->GetNbinsX()+1; ++i) h_events_Y_raw->SetBinError(i, 0);
                            h_events_Y_raw->Draw("HIST");
                            intermediate_canvas->SaveAs(pdf_path.c_str());

                            intermediate_canvas->Clear();
                            intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                            h_mc_flux_Y_raw->SetTitle(Form("MC Generation for %s (%s);E_{k}/n [GeV/n];MC Generated", secondary.c_str(), det.c_str()));
                            h_mc_flux_Y_raw->SetMarkerColor(kGreen+2); h_mc_flux_Y_raw->SetLineColor(kGreen+2);
                            for(int i=0; i<=h_mc_flux_Y_raw->GetNbinsX()+1; ++i) h_mc_flux_Y_raw->SetBinError(i, 0);
                            h_mc_flux_Y_raw->Draw("HIST");
                            intermediate_canvas->SaveAs(pdf_path.c_str());
                        }
                        unique_ptr<TFile> file_XfragY(TFile::Open(acc_path_XfragY.c_str()));
                         if (file_XfragY && !file_XfragY->IsZombie()) {
                            unique_ptr<TH1> h_events_XfragY_raw(getHistClone(file_XfragY.get(), eventHistName));
                            unique_ptr<TH1> h_mc_flux_XfragY_raw(getHistClone(file_XfragY.get(), "MC_FLUX_H3"));
                            
                            intermediate_canvas->Clear();
                            intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                            h_events_XfragY_raw->SetTitle(Form("Event Count for %s#rightarrow%s (%s);E_{k}/n [GeV/n];Counts", primary.c_str(), secondary.c_str(), det.c_str()));
                            h_events_XfragY_raw->SetMarkerColor(kRed); h_events_XfragY_raw->SetLineColor(kRed);
                            for(int i=0; i<=h_events_XfragY_raw->GetNbinsX()+1; ++i) h_events_XfragY_raw->SetBinError(i, 0);
                            h_events_XfragY_raw->Draw("HIST");
                            intermediate_canvas->SaveAs(pdf_path.c_str());

                            intermediate_canvas->Clear();
                            intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                            h_mc_flux_XfragY_raw->SetTitle(Form("MC Generation for %s#rightarrow%s (%s);E_{k}/n [GeV/n];MC Generated", primary.c_str(), secondary.c_str(), det.c_str()));
                            h_mc_flux_XfragY_raw->SetMarkerColor(kMagenta); h_mc_flux_XfragY_raw->SetLineColor(kMagenta);
                            for(int i=0; i<=h_mc_flux_XfragY_raw->GetNbinsX()+1; ++i) h_mc_flux_XfragY_raw->SetBinError(i, 0);
                            h_mc_flux_XfragY_raw->Draw("HIST");
                            intermediate_canvas->SaveAs(pdf_path.c_str());
                        }
                    }

                    unique_ptr<TH1> h_acc_Y(calculateAcceptance(acc_path_Y, eventHistName));
                    unique_ptr<TH1> h_acc_XfragY(calculateAcceptance(acc_path_XfragY, eventHistName));

                    intermediate_canvas->Clear();
                    intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                    h_acc_Y->SetTitle(Form("Acceptance for %s (%s);E_{k}/n [GeV/n];Acceptance", secondary.c_str(), det.c_str()));
                    h_acc_Y->SetMarkerColor(kBlue); h_acc_Y->SetLineColor(kBlue);
                    for(int i=0; i<=h_acc_Y->GetNbinsX()+1; ++i) h_acc_Y->SetBinError(i, 0);
                    h_acc_Y->Draw("HIST");
                    intermediate_canvas->SaveAs(pdf_path.c_str());
                    
                    intermediate_canvas->Clear();
                    intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                    h_acc_XfragY->SetTitle(Form("Acceptance for %s#rightarrow%s (%s);E_{k}/n [GeV/n];Acceptance", primary.c_str(), secondary.c_str(), det.c_str()));
                    h_acc_XfragY->SetMarkerColor(kRed); h_acc_XfragY->SetLineColor(kRed);
                    for(int i=0; i<=h_acc_XfragY->GetNbinsX()+1; ++i) h_acc_XfragY->SetBinError(i, 0);
                    h_acc_XfragY->Draw("HIST");
                    intermediate_canvas->SaveAs(pdf_path.c_str());

                    if (!h_final_epsilon) {
                        h_final_epsilon = (TH1*)h_acc_Y->Clone(Form("h_epsilon_%s_to_%s", primary.c_str(), secondary.c_str()));
                        h_final_epsilon->SetTitle(channel_label.c_str());
                        h_final_epsilon->Reset();
                    }

                    unique_ptr<TH1> h_acc_ratio((TH1*)h_acc_XfragY->Clone(Form("h_acc_ratio_%s_%s_%s", primary.c_str(), secondary.c_str(), det.c_str())));
                    h_acc_ratio->Divide(h_acc_Y.get());

                    unique_ptr<TF1> fit_acc_ratio(smoothRatio(h_acc_ratio.get(), det, Form("fit_acc_ratio_%s_%s_%s", primary.c_str(), secondary.c_str(), det.c_str())));
                    if (!fit_acc_ratio) { cout << "    WARNING: Skipping detector " << det << " due to failed spline fit." << endl; continue; }

                    intermediate_canvas->Clear();
                    intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                    h_acc_ratio->SetTitle(Form("Acceptance Ratio & Smooth (%s);E_{k}/n [GeV/n];Acc(X#rightarrowY)/Acc(Y)", det.c_str()));
                    h_acc_ratio->SetMarkerStyle(20); h_acc_ratio->SetMarkerColor(kBlack);
                    for(int i=0; i<=h_acc_ratio->GetNbinsX()+1; ++i) h_acc_ratio->SetBinError(i, 0);
                    h_acc_ratio->Draw("HIST");
                    fit_acc_ratio->SetLineColor(kRed);
                    fit_acc_ratio->Draw("SAME");
                    {
                        unique_ptr<TLegend> leg(new TLegend(0.2, 0.7, 0.5, 0.88));
                        leg->AddEntry(h_acc_ratio.get(), "Raw Ratio", "HIST");
                        leg->AddEntry(fit_acc_ratio.get(), "Smooth", "l");
                        //leg->Draw();
                        intermediate_canvas->SaveAs(pdf_path.c_str());
                    }

                    for (int i = 1; i <= h_final_epsilon->GetNbinsX(); ++i) {
                        double ek_n = h_final_epsilon->GetBinCenter(i);
                        if (ek_n >= DETECTOR_RANGES.at(det).first && ek_n < DETECTOR_RANGES.at(det).second) {
                            // Use the pre-calculated TF1 for flux ratio
                            double flux_ratio_val = tf1_flux_ratio->Eval(ek_n);
                            double acc_ratio_spline_val = fit_acc_ratio->Eval(ek_n);
                            
                            double epsilon_val = flux_ratio_val * acc_ratio_spline_val;
                            h_final_epsilon->SetBinContent(i, epsilon_val);
                            h_final_epsilon->SetBinError(i, 0.0);

                            // --- The original debugging block remains unchanged in its logic ---
                            string debug_key = channel_label + det;
                            if (true) { // Kept 'true' as in the user's last version
                                int DEBUG_BIN_INDEX = h_final_epsilon->FindBin(ek_n);
                                cout << "\n\n/********** DEBUGGING OUTPUT (BIN " << DEBUG_BIN_INDEX << ") **********/" << endl;
                                cout << "Channel: " << primary << " -> " << secondary << ", Detector: " << det << endl;
                                cout << "Energy (ek/n): " << ek_n << " GeV/n" << endl;
                                cout << "----------------------------------------" << endl;
                                cout << fixed << setprecision(10); 

                                cout << "--- Acceptance Calculation Details ---" << endl;
                                {
                                    unique_ptr<TFile> file_Y(TFile::Open(acc_path_Y.c_str()));
                                    unique_ptr<TFile> file_XfragY(TFile::Open(acc_path_XfragY.c_str()));
                                    if(file_Y && file_XfragY){
                                        unique_ptr<TH1> h_events_Y(getHistClone(file_Y.get(), eventHistName));
                                        unique_ptr<TH1> h_mc_flux_Y(getHistClone(file_Y.get(), "MC_FLUX_H3"));
                                        unique_ptr<TH1> h_events_XfragY(getHistClone(file_XfragY.get(), eventHistName));
                                        unique_ptr<TH1> h_mc_flux_XfragY(getHistClone(file_XfragY.get(), "MC_FLUX_H3"));

                                        if(h_events_Y && h_mc_flux_Y && h_events_XfragY && h_mc_flux_XfragY){
                                            double scale = TMath::Power(3.9, 2) * TMath::Pi();
                                            cout << "Scale Factor (3.9^2 * Pi): " << scale << endl;
                                            
                                            double raw_events_Y = h_events_Y->GetBinContent(DEBUG_BIN_INDEX);
                                            double raw_mc_flux_Y = h_mc_flux_Y->GetBinContent(DEBUG_BIN_INDEX);
                                            double calc_acc_Y = (raw_mc_flux_Y > 0) ? (raw_events_Y / raw_mc_flux_Y) * scale : 0.0;
                                            cout << "Acc(Y): " << endl;
                                            cout << "  Raw Events: " << raw_events_Y << " (from " << acc_path_Y << ")" << endl;
                                            cout << "  MC Flux:    " << raw_mc_flux_Y << endl;
                                            cout << "  => Calculated Acc(Y) = (Events/Flux)*Scale = " << calc_acc_Y << endl;
                                            cout << "  (Value from function: " << h_acc_Y->GetBinContent(DEBUG_BIN_INDEX) << ")" << endl;

                                            double raw_events_XfragY = h_events_XfragY->GetBinContent(DEBUG_BIN_INDEX);
                                            double raw_mc_flux_XfragY = h_mc_flux_XfragY->GetBinContent(DEBUG_BIN_INDEX);
                                            double calc_acc_XfragY = (raw_mc_flux_XfragY > 0) ? (raw_events_XfragY / raw_mc_flux_XfragY) * scale : 0.0;
                                            cout << "Acc(XfragY): " << endl;
                                            cout << "  Raw Events: " << raw_events_XfragY << " (from " << acc_path_XfragY << ")" << endl;
                                            cout << "  MC Flux:    " << raw_mc_flux_XfragY << endl;
                                            cout << "  => Calculated Acc(XfragY) = (Events/Flux)*Scale = " << calc_acc_XfragY << endl;
                                            cout << "  (Value from function: " << h_acc_XfragY->GetBinContent(DEBUG_BIN_INDEX) << ")" << endl;
                                        } else { cout << "Could not retrieve raw histograms for acceptance debug." << endl; }
                                    } else { cout << "Could not open files for acceptance debug." << endl; }
                                }
                                cout << "----------------------------------------" << endl;

                                cout << "--- Epsilon Calculation Details ---" << endl;
                                cout << "Flux(" << primary << ") at " << ek_n << " GeV/n:   " << flux_X->Eval(ek_n) << endl;
                                cout << "Flux(" << secondary << ") at " << ek_n << " GeV/n: " << flux_Y->Eval(ek_n) << endl;
                                cout << "=> Flux Ratio (X/Y) from TF1:  " << flux_ratio_val << endl;
                                cout << "----------------------------------------" << endl;
                                cout << "Acceptance Ratio (from raw hist): " << h_acc_ratio->GetBinContent(DEBUG_BIN_INDEX) << endl;
                                cout << "Acceptance Ratio (from spline):   " << acc_ratio_spline_val << endl;
                                cout << "----------------------------------------" << endl;
                                cout << "FINAL Epsilon = Flux Ratio * Splined Acc Ratio" << endl;
                                cout << "              = " << flux_ratio_val << " * " << acc_ratio_spline_val << endl;
                                cout << "              = " << epsilon_val << endl;
                                cout << "/**************************************/\n\n" << endl;
                                
                                debug_printed_for_channel[debug_key] = true;
                            }
                        }
                    }
                } catch (const std::runtime_error& e) {
                    cerr << "    ERROR processing " << det << " for " << channel_label << ": " << e.what() << endl;
                }
            } // Detector loop
            
            if (h_final_epsilon) {
                results_by_secondary[secondary].push_back(h_final_epsilon);
                labels_by_secondary[secondary].push_back(channel_label);

                // =================================================================================
                // MODIFICATION 4: Add final epsilon result to the PDF.
                // =================================================================================
                intermediate_canvas->Clear();
                intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                h_final_epsilon->SetTitle(Form("Final #epsilon for %s;E_{k}/n [GeV/n];#epsilon", channel_label.c_str()));
                h_final_epsilon->SetMarkerColor(kBlack); h_final_epsilon->SetLineColor(kBlack);
                h_final_epsilon->SetMarkerStyle(20);
                // Errors are already set to 0 during calculation
                h_final_epsilon->SetMinimum(0);
                h_final_epsilon->Draw("HIST");
                intermediate_canvas->SaveAs(pdf_path.c_str());
            }
            
            intermediate_canvas->SaveAs((pdf_path + "]").c_str());
            cout << "    INFO: Intermediate results saved to " << pdf_path << endl;

        } // Primary loop
    } // Secondary loop

    cout << "\n--- All calculations finished. Starting to generate plots. ---\n" << endl;

    // --- Plotting Loop (UNCHANGED) ---
    for (const auto& pair : results_by_secondary) {
        const string& secondary_name = pair.first;
        const vector<TH1*>& hists_to_plot = pair.second;
        const vector<string>& labels = labels_by_secondary.at(secondary_name);

        if (hists_to_plot.empty()) continue;

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

        y_min = (y_min < 1.0e9) ? y_min * 0.5 : 1e-6;
        y_max = (y_max > 0) ? y_max * 2.0 : 1.0;

        TCanvas* c1 = new TCanvas(Form("c_%s", secondary_name.c_str()), Form("Epsilon for %s", secondary_name.c_str()), 1000, 750);
        c1->SetGrid(); c1->SetLogx(); //c1->SetLogy();

        TH1* frame = c1->DrawFrame(0.2, y_min, 25.0, y_max);
        frame->SetTitle(Form("Fragmentation Contribution to %s;E_{k}/n [GeV/n];#epsilon (Fragmentation Fraction)", secondary_name.c_str()));
        frame->GetYaxis()->SetTitleOffset(1.2);
        frame->GetYaxis()->SetRangeUser(0, y_max);

        TLegend* legend = new TLegend(0.7, 0.65, 0.9, 0.88);
        legend->SetTextSize(0.03); 
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);

        for (size_t i = 0; i < hists_to_plot.size(); ++i) {
            TH1* h = hists_to_plot[i];
            int color = PLOT_COLORS[i % PLOT_COLORS.size()];
            h->SetLineColor(color); h->SetMarkerColor(color);
            h->SetMarkerStyle(20); h->SetMarkerSize(1.0);
            h->Draw("P SAME");
            legend->AddEntry(h, labels[i].c_str(), "HIST");
        }

        legend->Draw();
        c1->Update();

        string output_filename = OUTPUT_DIR + "frag_epsilon_" + secondary_name + ".png";
        c1->SaveAs(output_filename.c_str());
        cout << ">>> Successfully created plot: " << output_filename << endl;

        delete c1;
        delete legend;
    }

    cout << "\n--- All plots generated. Script finished. ---\n" << endl;
}