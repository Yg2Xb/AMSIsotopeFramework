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
std::vector<double> xpoints = {0.33, 0.77, 1.41, 3.16, 5.62, 10.00, 30, 60, 100};
std::vector<double> xpoints2 = {0.33, 1.41, 5.62, 10.00, 30, 60, 100};
// =================================================================================
// MODIFICATION 1: Add a global boolean to control rebinning.
// =================================================================================
const bool DO_REBIN = true;

const vector<string> ALL_PARTICLES = {"Be7", "Be9", "Be10", "B10", "B11", "C12", "N14", "N15","O16"};
//const vector<string> ALL_PARTICLES = {"Be7", "Be9", "Be10", "C12"};
const vector<string> ALL_SECONDARY_ISOTOPES = {"B10", "B11", "Be7", "Be9", "Be10"};
//const vector<string> ALL_SECONDARY_ISOTOPES = {"Be7", "Be9", "Be10"};

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
    
    // Apply rebinning based on the global flag
    if (DO_REBIN) {
        h_clone->Rebin(2);
    }
    h_clone->Sumw2(); // Ensure sum of squares of weights is stored for error calculation

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
        // Rebinning is done in getHistClone, so both should be rebinned consistently.
        // This check remains valid.
        throw runtime_error("Bin mismatch in " + accFilePath);
    }
    h_events->Divide(h_events.get(), h_mc_flux.get(), 1, 1, "B");
    h_events->Scale(TMath::Power(3.9, 2) * TMath::Pi());
    return h_events.release();
}

// =================================================================================
// MODIFICATION: The 'detector' argument is no longer needed as the fit is global.
// The function body already used the global 'xpoints', so this simplifies the signature.
// =================================================================================
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
        fluxes[particle] = dynamic_cast<TF1*>(fluxFile->Get((particle + "_spline_Ek").c_str()));
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
            if(primary == "Be7") {
                cout << "Skipping channel: " << primary << " -> " << secondary << " (both are secondaries)." << endl;
                continue;
            }
            cout << "\n>>> Processing channel: " << primary << " -> " << secondary << "..." << endl;
            string channel_label = primary + " #rightarrow " + secondary;
            
            string pdf_path = OUTPUT_DIR + "rew_results_" + primary + "_to_" + secondary + ".pdf";
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
            flux_X->GetXaxis()->SetRangeUser(0.3, 100.0);
            flux_X->Draw();
            flux_Y->Draw("SAME");
            {
                unique_ptr<TLegend> leg(new TLegend(0.78, 0.78, 0.88, 0.88));
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(flux_X, Form("Flux(%s)", primary.c_str()), "l");
                leg->AddEntry(flux_Y, Form("Flux(%s)", secondary.c_str()), "l");
                leg->Draw(); // Re-enabled legend drawing as it was likely commented out by mistake
                intermediate_canvas->SaveAs(pdf_path.c_str());
            }

            auto flux_ratio_func = [&](double *x, double *p) {
                double fy = flux_Y->Eval(x[0]);
                return (fy > 0) ? flux_X->Eval(x[0]) / fy : 0.0;
            };
            unique_ptr<TF1> tf1_flux_ratio(new TF1(Form("tf1_flux_ratio_%s_to_%s", primary.c_str(), secondary.c_str()), flux_ratio_func, 0.3, 20.5, 0));
            tf1_flux_ratio->SetTitle(Form("Flux Ratio for %s -> %s;E_{k}/n [GeV/n];Flux(%s) / Flux(%s)", primary.c_str(), secondary.c_str(), primary.c_str(), secondary.c_str()));
            tf1_flux_ratio->SetLineColor(kGreen+2);
            cout<<" flux ratio at 1 GeV/n: "<<tf1_flux_ratio->Eval(1.0)<<endl;
            
            // --- PDF Page 2: Flux Ratio TF1 ---
            intermediate_canvas->Clear();
            intermediate_canvas->SetLogx(true);
            intermediate_canvas->SetLogy(false);
            //tf1_flux_ratio->GetXaxis()->SetRangeUser(0.3, 20.5);
            tf1_flux_ratio->Draw();
            intermediate_canvas->SaveAs(pdf_path.c_str());

            // =================================================================================
            // NEW LOGIC: Prepare a combined histogram for acceptance ratios.
            // =================================================================================
            unique_ptr<TH1> h_acc_ratio_combined(nullptr);
            TH1* h_final_epsilon = nullptr;

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
                    
                    // The block for plotting raw histograms is preserved as requested.
                    {
                        unique_ptr<TFile> file_Y(TFile::Open(acc_path_Y.c_str()));
                        if (file_Y && !file_Y->IsZombie()) {
                            unique_ptr<TH1> h_events_Y_raw(getHistClone(file_Y.get(), eventHistName));
                            unique_ptr<TH1> h_mc_flux_Y_raw(getHistClone(file_Y.get(), "MC_FLUX_H3"));
                            
                            intermediate_canvas->Clear();
                            intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                            h_events_Y_raw->SetTitle(Form("Event Count for %s (%s);E_{k}/n [GeV/n];Counts", secondary.c_str(), det.c_str()));
                            h_events_Y_raw->SetMarkerColor(kBlue); h_events_Y_raw->SetLineColor(kBlue);
                            //for(int i=0; i<=h_events_Y_raw->GetNbinsX()+1; ++i) h_events_Y_raw->SetBinError(i, 0);
                            h_events_Y_raw->GetXaxis()->SetRangeUser(0.3, 20.5);
                            h_events_Y_raw->Draw("PZ");
                            intermediate_canvas->SaveAs(pdf_path.c_str());

                            intermediate_canvas->Clear();
                            intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                            h_mc_flux_Y_raw->SetTitle(Form("MC Generation for %s (%s);E_{k}/n [GeV/n];MC Generated", secondary.c_str(), det.c_str()));
                            h_mc_flux_Y_raw->SetMarkerColor(kGreen+2); h_mc_flux_Y_raw->SetLineColor(kGreen+2);
                            //for(int i=0; i<=h_mc_flux_Y_raw->GetNbinsX()+1; ++i) h_mc_flux_Y_raw->SetBinError(i, 0);
                            h_mc_flux_Y_raw->GetXaxis()->SetRangeUser(0.3, 20.5);
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
                            //for(int i=0; i<=h_events_XfragY_raw->GetNbinsX()+1; ++i) h_events_XfragY_raw->SetBinError(i, 0);
                            h_events_XfragY_raw->GetXaxis()->SetRangeUser(0.3, 20.5);
                            h_events_XfragY_raw->Draw("PZ");
                            intermediate_canvas->SaveAs(pdf_path.c_str());

                            intermediate_canvas->Clear();
                            intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                            h_mc_flux_XfragY_raw->SetTitle(Form("MC Generation for %s#rightarrow%s (%s);E_{k}/n [GeV/n];MC Generated", primary.c_str(), secondary.c_str(), det.c_str()));
                            h_mc_flux_XfragY_raw->SetMarkerColor(kMagenta); h_mc_flux_XfragY_raw->SetLineColor(kMagenta);
                            //for(int i=0; i<=h_mc_flux_XfragY_raw->GetNbinsX()+1; ++i) h_mc_flux_XfragY_raw->SetBinError(i, 0);
                            h_mc_flux_XfragY_raw->GetXaxis()->SetRangeUser(0.3, 20.5);
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
                    //for(int i=0; i<=h_acc_Y->GetNbinsX()+1; ++i) h_acc_Y->SetBinError(i, 0);
                    h_acc_Y->GetXaxis()->SetRangeUser(0.3, 20.5);
                    h_acc_Y->Draw("PZ");
                    intermediate_canvas->SaveAs(pdf_path.c_str());
                    
                    intermediate_canvas->Clear();
                    intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                    h_acc_XfragY->SetTitle(Form("Acceptance for %s#rightarrow%s (%s);E_{k}/n [GeV/n];Acceptance", primary.c_str(), secondary.c_str(), det.c_str()));
                    h_acc_XfragY->SetMarkerColor(kRed); h_acc_XfragY->SetLineColor(kRed);
                    //for(int i=0; i<=h_acc_XfragY->GetNbinsX()+1; ++i) h_acc_XfragY->SetBinError(i, 0);
                    h_acc_XfragY->GetXaxis()->SetRangeUser(0.3, 20.5);
                    h_acc_XfragY->Draw("PZ");
                    intermediate_canvas->SaveAs(pdf_path.c_str());

                    if (!h_final_epsilon) {
                        h_final_epsilon = (TH1*)h_acc_Y->Clone(Form("h_epsilon_%s_to_%s", primary.c_str(), secondary.c_str()));
                        h_final_epsilon->SetTitle(channel_label.c_str());
                        h_final_epsilon->Reset();
                    }

                    unique_ptr<TH1> h_acc_ratio_det((TH1*)h_acc_XfragY->Clone(Form("h_acc_ratio_%s_%s_%s", primary.c_str(), secondary.c_str(), det.c_str())));
                    h_acc_ratio_det->Divide(h_acc_Y.get());

                    // =================================================================================
                    // NEW LOGIC: Combine detector-specific ratios into one histogram.
                    // =================================================================================
                    if (!h_acc_ratio_combined) {
                        h_acc_ratio_combined.reset((TH1*)h_acc_ratio_det->Clone(Form("h_acc_ratio_combined_%s_to_%s", primary.c_str(), secondary.c_str())));
                        h_acc_ratio_combined->Reset();
                        h_acc_ratio_combined->SetTitle("Combined Acceptance Ratio (All Detectors)");
                        h_acc_ratio_combined->SetMarkerSize(1.2);
                        h_acc_ratio_combined->SetMarkerStyle(20);
                    }

                    for (int i = 1; i <= h_acc_ratio_det->GetNbinsX(); ++i) {
                        double ek_n = h_acc_ratio_det->GetXaxis()->GetBinCenter(i);
                        if (ek_n >= DETECTOR_RANGES.at(det).first && ek_n < DETECTOR_RANGES.at(det).second) {
                            h_acc_ratio_combined->SetBinContent(i, h_acc_ratio_det->GetBinContent(i));
                            h_acc_ratio_combined->SetBinError(i, h_acc_ratio_det->GetBinError(i));
                        }
                    }

                    // The individual fit and epsilon calculation are now removed from this loop.
                    // The debug block is also removed from here as its context is now global.

                } catch (const std::runtime_error& e) {
                    cerr << "    ERROR processing " << det << " for " << channel_label << ": " << e.what() << endl;
                }
            } // Detector loop
            
            // =================================================================================
            // NEW LOGIC: Perform a single fit on the combined histogram after the loop.
            // =================================================================================
            if (!h_acc_ratio_combined) {
                cerr << "ERROR: Combined acceptance ratio histogram could not be created for " << channel_label << ". Skipping." << endl;
                intermediate_canvas->SaveAs((pdf_path + "]").c_str());
                continue;
            }

            // --- PDF Page: Combined Ratio and Final Fit ---
            intermediate_canvas->Clear();
            intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
            h_acc_ratio_combined->SetTitle(Form("Acceptance Ratio & Smooth (%s);E_{k}/n [GeV/n];Acc(X#rightarrowY)/Acc(Y)", "Combined"));
            h_acc_ratio_combined->SetMarkerStyle(20); h_acc_ratio_combined->SetMarkerColor(kBlack);
            //for(int i=0; i<=h_acc_ratio_combined->GetNbinsX()+1; ++i) h_acc_ratio_combined->SetBinError(i, 0);
            //h_acc_ratio_combined->GetXaxis()->SetRangeUser(0.3, 20.5);
            h_acc_ratio_combined->Draw("PZ");
            
            unique_ptr<TF1> fit_acc_ratio_combined(smoothRatio(h_acc_ratio_combined.get(), Form("fit_acc_ratio_combined_%s_to_%s", primary.c_str(), secondary.c_str())));
            if (!fit_acc_ratio_combined) {
                cerr << "    WARNING: Global spline fit failed for " << channel_label << ". Skipping final calculation." << endl;
                intermediate_canvas->SaveAs((pdf_path + "]").c_str());
                continue;
            }
            fit_acc_ratio_combined->SetLineColor(kRed);
            fit_acc_ratio_combined->GetXaxis()->SetRangeUser(0.3, 20.5);
            fit_acc_ratio_combined->Draw("SAME");
            {
                unique_ptr<TLegend> leg(new TLegend(0.3, 0.7, 0.5, 0.88));
                leg->AddEntry(h_acc_ratio_combined.get(), "Combined Raw Ratio", "l"); // Use "l" for HIST
                leg->AddEntry(fit_acc_ratio_combined.get(), "Global Smooth", "l");
                //leg->Draw(); // Re-enabled legend drawing
                intermediate_canvas->SaveAs(pdf_path.c_str());
            }

            // =================================================================================
            // NEW LOGIC: Calculate final epsilon using the single global fit.
            // The original debug block is now placed here, adapted to the new logic.
            // =================================================================================
            for (int i = 1; i <= h_final_epsilon->GetNbinsX(); ++i) {
                double ek_n = h_final_epsilon->GetBinCenter(i);
                
                // Find which detector this energy bin belongs to
                string current_det = "";
                for (const auto& det_pair : DETECTOR_RANGES) {
                    if (ek_n >= det_pair.second.first && ek_n < det_pair.second.second) {
                        current_det = det_pair.first;
                        break;
                    }
                }

                if (!current_det.empty()) {
                    double flux_ratio_val = tf1_flux_ratio->Eval(ek_n);
                    double acc_ratio_spline_val = fit_acc_ratio_combined->Eval(ek_n);
                    
                    double epsilon_val = flux_ratio_val * acc_ratio_spline_val;
                    double epsilon_err = 0;// need study
                    h_final_epsilon->SetBinContent(i, epsilon_val);
                    h_final_epsilon->SetBinError(i, epsilon_err);

                    // --- The original debugging block, adapted for the new global context ---
                    string debug_key = channel_label + current_det;
                    if (true) { // Kept 'true' as in the user's last version
                        int DEBUG_BIN_INDEX = i;
                        cout << "\n\n/********** DEBUGGING OUTPUT (BIN " << DEBUG_BIN_INDEX << ") **********/" << endl;
                        cout << "Channel: " << primary << " -> " << secondary << ", Detector: " << current_det << endl;
                        cout << "Energy (ek/n): " << ek_n << " GeV/n" << endl;
                        cout << "----------------------------------------" << endl;
                        cout << fixed << setprecision(10); 

                        cout << "--- Acceptance Calculation Details ---" << endl;
                        {
                            string acc_path_Y    = ACC_FILE_DIR + secondary + (secondary.rfind("Be", 0) == 0 ? "_rew_frag4.root" : "_rew_frag5.root");
                            string acc_path_XfragY = ACC_FILE_DIR + primary + (secondary.rfind("Be", 0) == 0 ? "_rew_frag4.root" : "_rew_frag5.root");
                            string chargePart = (secondary.rfind("Be", 0) == 0) ? "_Z4_" : "_Z5_";
                            string massPart = (secondary.rfind("Be", 0) == 0) ? secondary.substr(2) : secondary.substr(1);
                            string eventHistName = "L1Inner_MC_BKG_H3a_" + current_det + chargePart + "Mass" + massPart;

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

                                    double raw_events_XfragY = h_events_XfragY->GetBinContent(DEBUG_BIN_INDEX);
                                    double raw_mc_flux_XfragY = h_mc_flux_XfragY->GetBinContent(DEBUG_BIN_INDEX);
                                    double calc_acc_XfragY = (raw_mc_flux_XfragY > 0) ? (raw_events_XfragY / raw_mc_flux_XfragY) * scale : 0.0;
                                    cout << "Acc(XfragY): " << endl;
                                    cout << "  Raw Events: " << raw_events_XfragY << " (from " << acc_path_XfragY << ")" << endl;
                                    cout << "  MC Flux:    " << raw_mc_flux_XfragY << endl;
                                    cout << "  => Calculated Acc(XfragY) = (Events/Flux)*Scale = " << calc_acc_XfragY << endl;

                                } else { cout << "Could not retrieve raw histograms for acceptance debug." << endl; }
                            } else { cout << "Could not open files for acceptance debug." << endl; }
                        }
                        cout << "----------------------------------------" << endl;

                        cout << "--- Epsilon Calculation Details ---" << endl;
                        cout << "Flux(" << primary << ") at " << ek_n << " GeV/n:   " << flux_X->Eval(ek_n) << endl;
                        cout << "Flux(" << secondary << ") at " << ek_n << " GeV/n: " << flux_Y->Eval(ek_n) << endl;
                        cout << "=> Flux Ratio (X/Y) from TF1:  " << flux_ratio_val << endl;
                        cout << "----------------------------------------" << endl;
                        cout << "Acceptance Ratio (from combined hist): " << h_acc_ratio_combined->GetBinContent(DEBUG_BIN_INDEX) << endl;
                        cout << "Acceptance Ratio (from global spline): " << acc_ratio_spline_val << endl;
                        cout << "----------------------------------------" << endl;
                        cout << "FINAL Epsilon = Flux Ratio * Splined Acc Ratio" << endl;
                        cout << "              = " << flux_ratio_val << " * " << acc_ratio_spline_val << endl;
                        cout << "              = " << epsilon_val << endl;
                        cout << "/**************************************/\n\n" << endl;
                        
                        debug_printed_for_channel[debug_key] = true;
                    }
                }
            }

            if (h_final_epsilon) {
                results_by_secondary[secondary].push_back(h_final_epsilon);
                labels_by_secondary[secondary].push_back(channel_label);

                intermediate_canvas->Clear();
                intermediate_canvas->SetLogx(true); intermediate_canvas->SetLogy(false);
                h_final_epsilon->SetTitle(Form("Final #epsilon for %s;E_{k}/n [GeV/n];#epsilon", channel_label.c_str()));
                h_final_epsilon->SetMarkerColor(kBlack); h_final_epsilon->SetLineColor(kBlack);
                h_final_epsilon->SetMarkerStyle(20);
                h_final_epsilon->SetMinimum(0);
                h_final_epsilon->GetXaxis()->SetRangeUser(0.3, 20.5);
                h_final_epsilon->Draw("P");
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
        y_max = (y_max > 0) ? y_max * 1.2 : 1.0;

        TCanvas* c1 = new TCanvas(Form("c_%s", secondary_name.c_str()), Form("Epsilon for %s", secondary_name.c_str()), 1000, 750);
        c1->SetGrid(); c1->SetLogx(true); //c1->SetLogy();

        TH1* frame = c1->DrawFrame(0.3, y_min, 20.5, y_max);
        frame->SetTitle(Form("Fragmentation Contribution to %s;E_{k}/n [GeV/n];#epsilon (Fragmentation Fraction)", secondary_name.c_str()));
        frame->GetYaxis()->SetTitleOffset(1.2);
        frame->GetYaxis()->SetRangeUser(0, y_max);

        TLegend* legend = new TLegend(0.7, 0.65, 0.9, 0.88);
        legend->SetTextSize(0.03); 
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);

        for (size_t i = 0; i < hists_to_plot.size(); ++i) {
            TH1* h = hists_to_plot[i];
            h->GetXaxis()->SetRangeUser(0.3, 20.5);
            int color = PLOT_COLORS[i % PLOT_COLORS.size()];
            h->SetLineColor(color); h->SetMarkerColor(color);
            h->SetMarkerStyle(20); h->SetMarkerSize(1.);
            h->Draw("P SAME");
            legend->AddEntry(h, labels[i].c_str(), "p"); // Corrected from "HIST" to "p" to match Draw option
        }

        legend->Draw();
        c1->Update();

        string output_filename = OUTPUT_DIR + "rew_frag_epsilon_" + secondary_name + ".png";
        c1->SaveAs(output_filename.c_str());
        cout << ">>> Successfully created plot: " << output_filename << endl;

        delete c1;
        delete legend;
    }

    cout << "\n--- All plots generated. Script finished. ---\n" << endl;
}