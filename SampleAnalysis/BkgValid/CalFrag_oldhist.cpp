#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <TPad.h>
#include <TLine.h>
#include <TF1.h>
#include <TBox.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <iomanip>

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <map>
#include <string>
#include <numeric>
#include <algorithm>
#include <sstream>

// ===================================================================
// ======================= DEBUG CONFIGURATION =======================
// ===================================================================
const bool DEBUG_MODE = false;
const std::string DEBUG_SOURCE = "Boron";
const std::string DEBUG_FRAGMENT = "Beryllium";
const std::string DEBUG_CHAIN = "L1Inner";
// ===================================================================

// --- Helper: ValueWithError (Error Propagation) ---
struct ValueWithError {
    double val = 0.0;
    double err = 0.0;
    ValueWithError() = default;
    ValueWithError(double v, double e) : val(v), err(e) {}
    ValueWithError operator+(const ValueWithError& o) const { return {val + o.val, std::hypot(err, o.err)}; }
    ValueWithError operator-(const ValueWithError& o) const { return {val - o.val, std::hypot(err, o.err)}; }
    ValueWithError operator*(double scalar) const { return {val * scalar, std::abs(err * scalar)}; }
    ValueWithError operator*(const ValueWithError& o) const {
        double p_val = val * o.val;
        double p_err = std::hypot(err * o.val, val * o.err);
        return {p_val, p_err};
    }
    ValueWithError operator/(const ValueWithError& o) const {
        if (o.val == 0.0) return {std::nan(""), std::nan("")};
        double q_val = val / o.val;
        double q_err = std::hypot(err / o.val, (val * o.err) / (o.val * o.val));
        return {q_val, q_err};
    }
    ValueWithError operator/(double scalar) const {
        if (scalar == 0.0) return {std::nan(""), std::nan("")};
        return {val / scalar, err / std::abs(scalar)};
    }
};
ValueWithError operator*(double scalar, const ValueWithError& ve) { return ve * scalar; }

// --- Helper: File and Histogram Operations ---
const std::map<std::string, std::string> particleNameToAbbr = {
    {"Beryllium", "Be"}, {"Boron", "B"}, {"Carbon", "C"}, {"Nitrogen", "N"}, {"Oxygen", "O"}
};

std::unique_ptr<TFile> openFile(const TString& filename, const char* option = "READ") {
    auto file = std::unique_ptr<TFile>(TFile::Open(filename.Data(), option));
    if (!file || file->IsZombie()) {
        throw std::runtime_error("Error opening file: " + std::string(filename.Data()));
    }
    if (DEBUG_MODE || std::string(option) == "RECREATE") {
        std::cout << "[INFO] Opened file: " << filename << std::endl;
    }
    return file;
}

// ADAPTIVE CHANGE: Template remains, but generally we will use TH1* or TH1D*
// This allows reading TH1F or TH1D transparently.
template<typename HistType = TH1>
HistType* getHistFromFile(TFile* file, const TString& histName) {
    if (!file) {
         throw std::runtime_error("Input TFile pointer is null when trying to get histogram: " + std::string(histName.Data()));
    }
    HistType* hist = nullptr;
    file->GetObject(histName, hist);
    if (!hist) {
        std::cout << "[WARNING] Could not find histogram: " << std::string(histName.Data()) << " in file " << file->GetName() << ". Assuming it's zero." << std::endl;
        return nullptr;
    }
    hist->SetDirectory(0);
    return hist;
}

void fillHistBin(TH1* hist, int bin_idx, const ValueWithError& vw) {
    if(hist) {
        hist->SetBinContent(bin_idx, vw.val);
        hist->SetBinError(bin_idx, vw.err);
    }
}

// --- Configuration for a specific analysis process ---
struct AnalysisConfig {
    std::string sourceParticle;
    std::string fragmentParticle;
    bool InIsoLevel;
    std::vector<std::string> primaryFitIsotopes;
    int useMass;
    std::string fragFileID;
    std::map<std::string, double> sourceIsotopeFractions;
    std::vector<std::string> fragmentIsotopes;
    std::map<std::string, int> fragZandMass;
    std::map<std::string, std::string> mcSourceFiles;
    std::string mcSourceGenHist = "MC_FLUX_H3";
};

// --- 68% Uncertainty Calculation (from plotR.cpp) ---
double calculate68PercentUncertainty_Linear(TH1* hist, TF1* fit_func, double fit_min, double fit_max) {
    if (!hist || !fit_func) return 0.0;
    
    int total_points = 0;
    for (int bin = hist->GetXaxis()->FindBin(fit_min); bin <= hist->GetXaxis()->FindBin(fit_max); ++bin) {
        double y = hist->GetBinContent(bin);
        double y_err = hist->GetBinError(bin);
        if (y != 0 && !std::isnan(y) && !std::isinf(y) && y_err > 0) {
            total_points++;
        }
    }
    
    if (total_points == 0) return 0.0;
    
    double delta = 0.0;
    double step = 0.001;
    int target_points = static_cast<int>(std::ceil(0.68 * total_points));
    
    while (delta < 10.0) {
        int points_within = 0;
        for (int bin = hist->GetXaxis()->FindBin(fit_min); bin <= hist->GetXaxis()->FindBin(fit_max); ++bin) {
            double x = hist->GetBinCenter(bin);
            double y = hist->GetBinContent(bin);
            double y_err = hist->GetBinError(bin);
            
            if (y != 0 && !std::isnan(y) && !std::isinf(y) && y_err > 0) {
                double fit_val = fit_func->Eval(x);
                if (y >= (fit_val - delta) && y <= (fit_val + delta)) {
                    points_within++;
                }
            }
        }
        
        if (points_within >= target_points) {
            return delta;
        }
        delta += step;
    }
    return delta;
}

// --- Plotting Helper (Modified to create two separate plots and add fits to comparison) ---
void createComparisonPlots(TH1* h_iss, TH1* h_mc, const std::string& title, 
                          const std::string& y_axis_title, const std::string& output_path) {
    if (!h_iss) {
        std::cerr << "[PLOT ERROR] ISS histogram is null for plot: " << title << std::endl;
        return;
    }

    // Renamed ranges to match the ratio plot's usage of 1.27-6.1 and 6.1-21.5
    const double x_min = 0.4;
    const double x_max = 21.0;
    // USER REQUESTED CHANGE 1: Low E fit range starts at 0.4 instead of 1.27
    const double fit_linear_min_low = 0.4; 
    const double fit_linear_max_low = 6.1;  
    const double fit_linear_min_high = 6.1; 
    const double fit_linear_max_high = 21.5; 
    
    // =================================================================================
    // ========== PLOT 1: ISS vs MC Comparison (New: Fits for both ISS and MC) ==========
    // =================================================================================
    TCanvas* c1 = new TCanvas("c1", title.c_str(), 800, 600);
    c1->SetBottomMargin(0.15);
    c1->SetTopMargin(0.12);
    c1->SetLeftMargin(0.18);
    c1->SetRightMargin(0.03);
    c1->SetGrid();
    
    // --- 1. Draw ISS Data ---
    h_iss->SetMarkerStyle(20);
    h_iss->SetMarkerColor(kBlack);
    h_iss->SetLineColor(kBlack);
    h_iss->SetTitle(title.c_str());
    h_iss->GetYaxis()->SetTitle(y_axis_title.c_str());
    h_iss->GetXaxis()->SetTitle("E_{k}/n [GeV]");
    h_iss->GetXaxis()->SetRangeUser(x_min, x_max);
    h_iss->SetStats(0);
    
    double max_y = 0;
    for(int i = 1; i <= h_iss->GetNbinsX(); ++i) {
        double x = h_iss->GetBinCenter(i);
        if (x > x_max || x < x_min) continue;
        max_y = std::max(max_y, h_iss->GetBinContent(i));
    }
    if (h_mc) {
        for(int i = 1; i <= h_mc->GetNbinsX(); ++i) {
            double x = h_mc->GetBinCenter(i);
            if (x > x_max || x < x_min) continue;
            max_y = std::max(max_y, h_mc->GetBinContent(i));
        }
    }
    h_iss->GetYaxis()->SetRangeUser(0, max_y * 1.2);
    h_iss->GetYaxis()->SetTitleOffset(1.7);
    h_iss->Draw("PZ");
    
    // --- Legend for Data/MC (Top Right) ---
    TLegend* legend1 = new TLegend(0.6, 0.70, 0.9, 0.88);
    legend1->SetBorderSize(1);
    legend1->SetFillStyle(0);
    legend1->AddEntry(h_iss, "ISS Data", "ep");
    
    if (h_mc) {
        h_mc->SetMarkerStyle(20);
        h_mc->SetMarkerColor(kRed);
        h_mc->SetLineColor(kRed);
        h_mc->Draw("PZ SAME");
        legend1->AddEntry(h_mc, "MC Simulation", "ep");
    }
    legend1->Draw();
    
    // --- Legend for Fits (Bottom Left) ---
    TLegend* legend_fits = new TLegend(0.3, 0.15, 0.8, 0.3); // Adjusted to bottom left
    legend_fits->SetBorderSize(1);
    legend_fits->SetFillStyle(0);
    legend_fits->SetTextSize(0.025);

    // Lambda to perform the fit and add to the legend
    auto perform_and_draw_fit = [&](TH1* hist, const char* name_prefix, double min_e, double max_e, 
                                    Color_t color, const char* label_prefix) {
        if (!hist) return (TF1*)nullptr;

        // Check for data in range
        bool hasData = false;
        for (int bin = hist->GetXaxis()->FindBin(min_e); 
             bin <= hist->GetXaxis()->FindBin(max_e); ++bin) {
            if (hist->GetBinContent(bin) != 0) {
                hasData = true;
                break;
            }
        }
        
        if (!hasData) return (TF1*)nullptr;
        
        TString fit_name = TString::Format("%s_%s_%.1f_%.1f", name_prefix, hist->GetName(), min_e, max_e);
        TF1* fit_func = new TF1(fit_name, "pol0", min_e, max_e);
        fit_func->SetLineColor(color);
        fit_func->SetLineWidth(2);
        
        // Use R option for range, Q for quiet, S for store fit results
        TFitResultPtr r = hist->Fit(fit_func, "SRQ+"); 
        
        if (r->IsValid() && r->Status() == 0) {
            double p0 = fit_func->GetParameter(0);
            
            // USER REQUESTED CHANGE 2: Display standard fit error (GetParError(0)) 
            // and use %.3f for precision.
            double fit_err = fit_func->GetParError(0); 
            // Calculate 68% uncertainty for the visual band width (as originally intended)
            double delta_for_band = calculate68PercentUncertainty_Linear(hist, fit_func, min_e, max_e); 
            delta_for_band = delta_for_band;
            
            // Draw uncertainty band
            const int n_points = 100;
            double x_arr[n_points * 2];
            double y_arr[n_points * 2];
            
            for (int i = 0; i < n_points; ++i) {
                double x = min_e + i * (max_e - min_e) / (n_points - 1);
                double y_center = fit_func->Eval(x);
                // Use the 68% uncertainty for the band width
                x_arr[i] = x;
                y_arr[i] = y_center + delta_for_band; 
                x_arr[2 * n_points - 1 - i] = x;
                y_arr[2 * n_points - 1 - i] = y_center - delta_for_band;
            }
            
            TGraph* uncertainty_band = new TGraph(2 * n_points, x_arr, y_arr);
            uncertainty_band->SetFillColorAlpha(color, 0.15);
            uncertainty_band->SetFillStyle(1001);
            uncertainty_band->SetLineColor(color);
            uncertainty_band->SetLineWidth(1);
            uncertainty_band->Draw("F SAME");
            fit_func->Draw("SAME");
            
            // Update label format to use fit_err and 3 decimal places
            TString fit_label = TString::Format("%s (%.1f-%.1f GeV): %.4f #pm %.4f", 
                                               label_prefix, min_e, max_e, p0, delta_for_band);
            cout<<"Nuc: "<<title<<", Fit Range: "<<min_e<<"-"<<max_e<<", Fit Value: "<<p0<<" +/- "<<fit_err<<", 68% Uncertainty Band: +/- "<<delta_for_band<<endl;
            legend_fits->AddEntry(uncertainty_band, fit_label, "lf");
            
            return fit_func;
        }
        delete fit_func;
        return (TF1*)nullptr;
    };


    // --- 2. ISS Fits ---
    perform_and_draw_fit(h_iss, "iss_low", fit_linear_min_low, fit_linear_max_low, kBlack, "ISS Data");
    perform_and_draw_fit(h_iss, "iss_high", fit_linear_min_high, fit_linear_max_high, kGray+2, "ISS Data");

    // --- 3. MC Fits ---
    if (h_mc) {
        perform_and_draw_fit(h_mc, "mc_low", fit_linear_min_low, fit_linear_max_low, kRed, "MC Sim.");
        perform_and_draw_fit(h_mc, "mc_high", fit_linear_min_high, fit_linear_max_high, kOrange+7, "MC Sim.");
    }
    
    // Redraw data points to be on top of the bands
    h_iss->Draw("PZ SAME");
    if (h_mc) h_mc->Draw("PZ SAME");
    
    legend_fits->Draw();
    
    c1->SaveAs((output_path + "_comparison.png").c_str());
    delete c1;
    
    // =================================================================================
    // ========== PLOT 2: ISS/MC Ratio with Fits (The original ratio plot) ==========
    // =================================================================================
    // 保持 Ratio Plot 的拟合范围不变 (1.27-6.1)
    if (!h_mc) {
        std::cout << "[INFO] No MC histogram, skipping ratio plot for: " << title << std::endl;
        return;
    }
    
    TCanvas* c2 = new TCanvas("c2", "Ratio", 800, 600);
    c2->SetTopMargin(0.05);
    c2->SetGrid();
    
    TH1* h_ratio = (TH1*)h_iss->Clone("h_ratio");
    h_ratio->Divide(h_mc);
    h_ratio->SetTitle("");
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetMarkerColor(kBlue);
    h_ratio->SetLineColor(kBlue);
    h_ratio->GetXaxis()->SetTitle("Measured E_{k}/n [GeV]");
    h_ratio->GetYaxis()->SetTitle("ISS/MC");
    h_ratio->GetYaxis()->SetTitleOffset(1.3);
    h_ratio->GetXaxis()->SetRangeUser(x_min, x_max);
    h_ratio->SetStats(0);
    
    TLegend* legend2 = new TLegend(0.15, 0.70, 0.76, 0.95);
    legend2->SetFillStyle(0);
    legend2->SetBorderSize(1);
    legend2->SetTextSize(0.028);
    legend2->AddEntry(h_ratio, (title+" ISS/MC").c_str(), "pe");
    
    // First linear fit (1.27 - 6.1 GeV)
    const double fit_linear_min_ratio_orig = 1.27; 
    const double fit_linear_max_ratio_orig = 6.1;
    
    bool hasDataLinear1 = false;
    for (int bin = h_ratio->GetXaxis()->FindBin(fit_linear_min_ratio_orig); 
         bin <= h_ratio->GetXaxis()->FindBin(fit_linear_max_ratio_orig); ++bin) {
        if (h_ratio->GetBinContent(bin) != 0) {
            hasDataLinear1 = true;
            break;
        }
    }
    
    TF1* fit_linear1 = nullptr;
    TGraph* uncertainty_band1 = nullptr;
    
    if (hasDataLinear1) {
        fit_linear1 = new TF1("fit_linear1", "pol0", fit_linear_min_ratio_orig, fit_linear_max_ratio_orig); 
        fit_linear1->SetLineColor(kMagenta);
        fit_linear1->SetLineWidth(2);
        
        TFitResultPtr r_linear1 = h_ratio->Fit(fit_linear1, "SRQ");
        
        if (r_linear1->IsValid() && r_linear1->Status() == 0) {
            double p0 = fit_linear1->GetParameter(0);
            double chi2 = r_linear1->Chi2();
            int ndf = r_linear1->Ndf();
            
            double delta_linear1 = calculate68PercentUncertainty_Linear(h_ratio, fit_linear1, 
                                                                        fit_linear_min_ratio_orig, fit_linear_max_ratio_orig); 
            
            const int n_points = 100;
            double x_arr[n_points * 2];
            double y_arr[n_points * 2];
            
            for (int i = 0; i < n_points; ++i) {
                double x = fit_linear_min_ratio_orig + i * (fit_linear_max_ratio_orig - fit_linear_min_ratio_orig) / (n_points - 1); 
                double y_center = fit_linear1->Eval(x);
                x_arr[i] = x;
                y_arr[i] = y_center + delta_linear1;
                x_arr[2 * n_points - 1 - i] = x;
                y_arr[2 * n_points - 1 - i] = y_center - delta_linear1;
            }
            
            uncertainty_band1 = new TGraph(2 * n_points, x_arr, y_arr);
            uncertainty_band1->SetFillColorAlpha(kMagenta, 0.25);
            uncertainty_band1->SetFillStyle(1001);
            uncertainty_band1->SetLineColor(kMagenta - 7);
            uncertainty_band1->SetLineWidth(1);
            
            TString fit_label = TString::Format("constant fit (1.27-6.1 GeV): %.2f #pm %.2f", 
                                               p0, delta_linear1);
            legend2->AddEntry(uncertainty_band1, fit_label, "lf");
            
            if (ndf > 0) {
                TString chi2_label = TString::Format("#chi^{2}/NDF = %.1f / %d = %.2f", 
                                                    chi2, ndf, chi2 / ndf);
                legend2->AddEntry((TObject*)0, chi2_label, "");
            }
        }
    }
    
    // Second linear fit (6.1 - 21.5 GeV)
    const double fit_second_min_ratio_orig = 6.1;
    const double fit_second_max_ratio_orig = 21.5;
    
    bool hasDataLinear2 = false;
    for (int bin = h_ratio->GetXaxis()->FindBin(fit_second_min_ratio_orig); 
         bin <= h_ratio->GetXaxis()->FindBin(fit_second_max_ratio_orig); ++bin) { 
        if (h_ratio->GetBinContent(bin) != 0) {
            hasDataLinear2 = true;
            break;
        }
    }
    
    TF1* fit_linear2 = nullptr;
    TGraph* uncertainty_band2 = nullptr;
    double min_y_second_fit = 1e9;
    
    if (hasDataLinear2) {
        fit_linear2 = new TF1("fit_linear2", "pol0", fit_second_min_ratio_orig, fit_second_max_ratio_orig); 
        fit_linear2->SetLineColor(kOrange + 1);
        fit_linear2->SetLineWidth(2);
        
        TFitResultPtr r_linear2 = h_ratio->Fit(fit_linear2, "SRQ+");
        
        if (r_linear2->IsValid() && r_linear2->Status() == 0) {
            double p0 = fit_linear2->GetParameter(0);
            double chi2 = r_linear2->Chi2();
            int ndf = r_linear2->Ndf();
            
            double delta_linear2 = calculate68PercentUncertainty_Linear(h_ratio, fit_linear2,
                                                                        fit_second_min_ratio_orig, fit_second_max_ratio_orig); 
            
            const int n_points = 100;
            for (int i = 0; i < n_points; ++i) {
                double x = fit_second_min_ratio_orig + i * (fit_second_max_ratio_orig - fit_second_min_ratio_orig) / (n_points - 1); 
                double y_center = fit_linear2->Eval(x);
                min_y_second_fit = std::min(min_y_second_fit, y_center - delta_linear2);
            }
            
            double x_arr[n_points * 2];
            double y_arr[n_points * 2];
            
            for (int i = 0; i < n_points; ++i) {
                double x = fit_second_min_ratio_orig + i * (fit_second_max_ratio_orig - fit_second_min_ratio_orig) / (n_points - 1); 
                double y_center = fit_linear2->Eval(x);
                x_arr[i] = x;
                y_arr[i] = y_center + delta_linear2;
                x_arr[2 * n_points - 1 - i] = x;
                y_arr[2 * n_points - 1 - i] = y_center - delta_linear2;
            }
            
            uncertainty_band2 = new TGraph(2 * n_points, x_arr, y_arr);
            uncertainty_band2->SetFillColorAlpha(kOrange + 1, 0.25);
            uncertainty_band2->SetFillStyle(1001);
            uncertainty_band2->SetLineColor(kOrange - 6);
            uncertainty_band2->SetLineWidth(1);
            
            TString fit_label = TString::Format("constant fit (6.1-21.5 GeV): %.2f #pm %.2f",
                                               p0, delta_linear2);
            legend2->AddEntry(uncertainty_band2, fit_label, "lf");
            
            if (ndf > 0) {
                TString chi2_label = TString::Format("#chi^{2}/NDF = %.1f / %d = %.2f", 
                                                    chi2, ndf, chi2 / ndf);
                legend2->AddEntry((TObject*)0, chi2_label, "");
            }
        }
    }
    
    double min_r = 1e9, max_r = -1e9;
    for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
        double x = h_ratio->GetBinCenter(i);
        if (x > x_max || x < x_min) continue;
        double y = h_ratio->GetBinContent(i);
        if (y != 0 && !std::isnan(y) && !std::isinf(y)) {
            min_r = std::min(min_r, y);
            max_r = std::max(max_r, y);
        }
    }
    
    if (min_y_second_fit < 1e9) {
        min_r = std::min(min_r, min_y_second_fit * 0.8);
    }
    
    if (min_r > max_r) { min_r = 0.2; max_r = 1.8; }
    h_ratio->GetYaxis()->SetRangeUser(min_r * 0.75, max_r * 1.25);
    
    h_ratio->Draw("PZ");
    
    if (uncertainty_band1) {
        uncertainty_band1->Draw("F SAME");
        fit_linear1->Draw("SAME");
    }
    
    if (uncertainty_band2) {
        uncertainty_band2->Draw("F SAME");
        fit_linear2->Draw("SAME");
    }
    
    legend2->Draw();
    c2->SaveAs((output_path + "_ratio.png").c_str());
    delete c2;
    delete h_ratio;
}

// --- MC Calculation Helpers ---
// ADAPTIVE CHANGE: Accept map of TH1* (generic), return TH1D* (Double precision)
TH1D* stitchHistograms(const std::map<std::string, TH1*>& hists, TH1* refHist) {
    // Create histogram using refHist bins, but ensure it is TH1D
    auto combined = new TH1D(TString::Format("%s_stitched", refHist->GetName()), 
                             refHist->GetTitle(), 
                             refHist->GetNbinsX(), 
                             refHist->GetXaxis()->GetXbins()->GetArray());
    combined->Reset();
    const double boundary1 = 1.28, boundary2 = 3.06;
    for (int i = 1; i <= combined->GetNbinsX(); ++i) {
        double binCenter = combined->GetBinCenter(i);
        std::string det = (binCenter < boundary1) ? "TOF" : (binCenter < boundary2) ? "NaF" : "AGL";
        if (hists.count(det) && hists.at(det) != nullptr) {
            TH1* h = hists.at(det);
            int source_bin = h->FindBin(binCenter);
            // GetBinContent works for both TH1F and TH1D, returns double
            combined->SetBinContent(i, h->GetBinContent(source_bin));
            combined->SetBinError(i, h->GetBinError(source_bin));
        }
    }
    return combined;
}

// ADAPTIVE CHANGE: Use TH1D for internal calculation logic
std::map<std::string, TH1D*> calculateMCRatios(const AnalysisConfig& config, const std::string& chain, 
                                                TH1* refHist, bool is_debug_target) {
    std::cout << "\n--- Calculating MC Ratios for " << config.sourceParticle << " -> " 
              << config.fragmentParticle << " ---\n";
    std::map<std::string, TH1D*> mc_ratios;
    const std::string mcBaseDir = "/eos/user/z/zixuan/Isotope/Add/";
    const std::array<std::string, 3> detectors = {"TOF", "NaF", "AGL"};
    const int nBins = refHist->GetNbinsX();

    // Use TH1D for accumulation
    TH1D* h_mc_total_source = new TH1D("h_mc_total_source_stitched", "", nBins, refHist->GetXaxis()->GetXbins()->GetArray());
    TH1D* h_mc_total_frag = new TH1D("h_mc_total_frag_stitched", "", nBins, refHist->GetXaxis()->GetXbins()->GetArray());
    
    std::map<std::string, TH1D*> h_mc_frag_isos;
    for (const auto& frag_iso : config.fragmentIsotopes) {
        h_mc_frag_isos[frag_iso] = new TH1D(TString::Format("h_mc_%s_stitched", frag_iso.c_str()), "", nBins, refHist->GetXaxis()->GetXbins()->GetArray());
    }

    std::map<std::string, std::unique_ptr<TH1>> h_gen_map;
    for (const auto& pair : config.mcSourceFiles) {
        const std::string& source_iso = pair.first;
        auto mc_file = openFile(mcBaseDir + pair.second);
        // Use generic TH1 GetObject
        h_gen_map[source_iso].reset(getHistFromFile<TH1>(mc_file.get(), config.mcSourceGenHist.c_str()));
        if(h_gen_map[source_iso]) h_gen_map[source_iso]->Rebin(2);
    }

    for (const auto& pair : config.mcSourceFiles) {
        const std::string& current_source_iso = pair.first;
        auto mc_file = openFile(mcBaseDir + pair.second);

        std::map<std::string, TH1*> source_hists_per_det;
        for(const auto& det : detectors) {
            source_hists_per_det[det] = getHistFromFile<TH1>(mc_file.get(), 
                TString::Format("%s_MC_BKG_H1_%s", chain.c_str(), det.c_str()));
            if(source_hists_per_det[det]) source_hists_per_det[det]->Rebin(2);
        }
        std::unique_ptr<TH1D> stitched_source(stitchHistograms(source_hists_per_det, refHist));

        std::map<std::string, std::unique_ptr<TH1D>> stitched_frags;
        for (const auto& frag_iso : config.fragmentIsotopes) {
            std::map<std::string, TH1*> frag_hists_per_det;
            int frag_Z = (config.fragmentParticle == "Beryllium") ? 4 : 5;
            for(const auto& det : detectors) {
                TString histname = TString::Format("%s_MC_BKG_H2_%s_Z%d_Mass%d", 
                    chain.c_str(), det.c_str(), frag_Z, config.fragZandMass.at(frag_iso));
                frag_hists_per_det[det] = getHistFromFile<TH1>(mc_file.get(), histname);
                if(frag_hists_per_det[det]) frag_hists_per_det[det]->Rebin(2);
            }
            stitched_frags[frag_iso].reset(stitchHistograms(frag_hists_per_det, refHist));
        }

        for (int i_bin = 1; i_bin <= nBins; ++i_bin) {
            if (refHist->GetBinCenter(i_bin) > 21) continue;

            double total_gen_in_bin = 0.0;
            for (const auto& inner_pair : config.mcSourceFiles) {
                const std::string& iso = inner_pair.first;
                if (h_gen_map.count(iso) && h_gen_map.at(iso)) {
                    total_gen_in_bin += h_gen_map.at(iso)->GetBinContent(i_bin);
                }
            }
            
            double n_gen_current = (h_gen_map.count(current_source_iso) && h_gen_map.at(current_source_iso)) ? 
                h_gen_map.at(current_source_iso)->GetBinContent(i_bin) : 0.0;
            
            double weight = 0.0;
            if (n_gen_current > 0) {
                weight = (total_gen_in_bin / n_gen_current) * config.sourceIsotopeFractions.at(current_source_iso);
            }

            ValueWithError raw_source_vw(stitched_source->GetBinContent(i_bin), stitched_source->GetBinError(i_bin));
            ValueWithError weighted_source_vw = raw_source_vw * weight;

            h_mc_total_source->SetBinContent(i_bin, h_mc_total_source->GetBinContent(i_bin) + weighted_source_vw.val);
            h_mc_total_source->SetBinError(i_bin, std::hypot(h_mc_total_source->GetBinError(i_bin), weighted_source_vw.err));
            
            for (const auto& frag_iso : config.fragmentIsotopes) {
                ValueWithError raw_frag_vw = {stitched_frags[frag_iso]->GetBinContent(i_bin), 
                                              stitched_frags[frag_iso]->GetBinError(i_bin)};
                ValueWithError weighted_frag_vw = raw_frag_vw * weight;
                
                h_mc_frag_isos[frag_iso]->SetBinContent(i_bin, 
                    h_mc_frag_isos[frag_iso]->GetBinContent(i_bin) + weighted_frag_vw.val);
                h_mc_frag_isos[frag_iso]->SetBinError(i_bin, 
                    std::hypot(h_mc_frag_isos[frag_iso]->GetBinError(i_bin), weighted_frag_vw.err));
                
                h_mc_total_frag->SetBinContent(i_bin, h_mc_total_frag->GetBinContent(i_bin) + weighted_frag_vw.val);
                h_mc_total_frag->SetBinError(i_bin, std::hypot(h_mc_total_frag->GetBinError(i_bin), weighted_frag_vw.err));
            }
        }
    }

    if (config.InIsoLevel) {
        for (const auto& frag_iso : config.fragmentIsotopes) {
            TH1D* ratio = (TH1D*)h_mc_frag_isos[frag_iso]->Clone(TString::Format("h_mc_ratio_%s", frag_iso.c_str()));
            ratio->Divide(h_mc_total_source);
            mc_ratios[frag_iso] = ratio;
        }
    }
    
    TH1D* total_ratio_hist = (TH1D*)h_mc_total_frag->Clone("h_mc_ratio_total");
    total_ratio_hist->Divide(h_mc_total_source);
    mc_ratios["Total"] = total_ratio_hist;

    delete h_mc_total_source;
    delete h_mc_total_frag;
    for(auto const& [key, val] : h_mc_frag_isos) delete val;
    
    return mc_ratios;
}

// CHANGED: Helper function to check if we should zero out specific MC points for specific isotopes
bool shouldZeroPoint(const AnalysisConfig& config, int bin_number_in_range, const std::string& isotope = "Total") {
    // For C to Be9 specifically: zero out 8th and 18th point
    if (config.sourceParticle == "Carbon" && isotope == "Be9") {
        if (bin_number_in_range == 8 || bin_number_in_range == 18) {
            return true;
        }
    }
    
    // For O to Be (Total only): zero out 7th point
    if (config.sourceParticle == "Oxygen" && config.fragmentParticle == "Beryllium" && isotope == "Total") {
        if (bin_number_in_range == 7) {
            return true;
        }
    }
    
    // For O to B (Total only): zero out 5th point
    if (config.sourceParticle == "Oxygen" && config.fragmentParticle == "Boron" && isotope == "Total") {
        if (bin_number_in_range == 5) {
            return true;
        }
    }
    
    return false;
}

// --- Main Calculation Function ---
void runAnalysis(const std::string& chain, const AnalysisConfig& config) {
    bool is_debug_target = (DEBUG_MODE && config.sourceParticle == DEBUG_SOURCE && 
                           config.fragmentParticle == DEBUG_FRAGMENT && chain == DEBUG_CHAIN);

    std::cout << "\n===================================================================\n";
    std::cout << "Starting Analysis: " << config.sourceParticle << " -> " << config.fragmentParticle
              << " for chain: " << chain << " (InIsoLevel Z: " << (config.InIsoLevel ? "Yes" : "No") << ")" << std::endl;
    std::cout << "===================================================================\n";

    const std::array<std::string, 3> detectors = {"TOF", "NaF", "AGL"};
    
    // ADAPTIVE CHANGE 2: Direct L1 Yield Path and Counts Path
    // Removed QFit path.
    TString l1YieldFilePath = "/eos/user/z/zixuan/Isotope/ChargeTemp/PureL1TempFit_AllL1QUnbiasedL1Inner.root";
    TString countsFilePath = TString::Format("/eos/user/z/zixuan/Isotope/Add/%s.root", config.fragFileID.c_str());
    TString outputFilePath = TString::Format("/eos/user/z/zixuan/Isotope/BkgValid/%s_to_%s_%s_Validation.root", 
        config.sourceParticle.c_str(), config.fragmentParticle.c_str(), chain.c_str());

    auto l1YieldFile = openFile(l1YieldFilePath);
    auto countsFile = openFile(countsFilePath);
    auto outputFile = openFile(outputFilePath, "RECREATE");

    // Use TH1D for calculation precision, regardless of file content type
    TH1F* refHist = getHistFromFile<TH1F>(countsFile.get(), 
        TString::Format("%s_ISS_BKG_H1_%s_TOF", chain.c_str(), config.sourceParticle.c_str()));
    if (!refHist) {
        throw std::runtime_error("Reference histogram for binning not found!");
    }
    refHist->Rebin(2);
    const TAxis* xAxis = refHist->GetXaxis();
    const int nBins = xAxis->GetNbins();
    const Double_t* binEdges = xAxis->GetXbins()->GetArray();

    // ADAPTIVE CHANGE: Use TH1* to hold whatever is read, allows for F or D
    std::map<std::string, TH1*> rawFragCountsHists;
    std::map<std::string, TH1*> l1SignalYieldHists; // Source from Source
    std::map<std::string, TH1*> l1ContamYieldHists; // Source from Fragment

    for (const auto& det : detectors) {
        // 1. Raw L2 Counts (Fragment candidates in L2)
        rawFragCountsHists[det] = getHistFromFile<TH1>(countsFile.get(), 
            TString::Format("%s_ISS_BKG_H3_%s_%s", chain.c_str(), config.sourceParticle.c_str(), det.c_str()));
        if(rawFragCountsHists[det]) rawFragCountsHists[det]->Rebin(2);
        
        // 2. Direct L1 Yields (Requirement 2)
        // Format: h_yield_in_{Source}_from_{Source}_{Det}
        TString sigName = TString::Format("h_yield_in_%s_from_%s_%s", 
                                          config.sourceParticle.c_str(), config.sourceParticle.c_str(), det.c_str());
        l1SignalYieldHists[det] = getHistFromFile<TH1>(l1YieldFile.get(), sigName);
        l1SignalYieldHists[det]->Rebin(2);

        // Format: h_yield_in_{Source}_from_{Fragment}_{Det}
        TString contamName = TString::Format("h_yield_in_%s_from_%s_%s", 
                                             config.fragmentParticle.c_str(), config.sourceParticle.c_str(), det.c_str());
        l1ContamYieldHists[det] = getHistFromFile<TH1>(l1YieldFile.get(), contamName);
        l1ContamYieldHists[det]->Rebin(2);
    }
    
    std::unique_ptr<TFile> fragmentFitFile, normalFitFile;
    if (config.InIsoLevel) {
        std::string fragAbbr = particleNameToAbbr.at(config.fragmentParticle);
        TString fragmentFitFilePath = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d_FragFrom%s.root", fragAbbr.c_str(), chain.c_str(), config.useMass, config.sourceParticle.c_str());
        TString normalFitFilePath = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d.root", fragAbbr.c_str(), chain.c_str(), config.useMass);
        fragmentFitFile = openFile(fragmentFitFilePath);
        normalFitFile   = openFile(normalFitFilePath);
    }

    // ADAPTIVE CHANGE: Always create TH1D for output calculations to maintain precision
    auto createHist = [&](const char* name, const char* title) {
        auto hist = new TH1D(name, title, nBins, binEdges);
        hist->SetStats(0);
        return hist;
    };
    
    std::map<std::string, TH1D*> h_iss_ratios;
    if (config.InIsoLevel) {
        for (const auto& iso : config.fragmentIsotopes) {
            h_iss_ratios[iso] = createHist(TString::Format("h_iss_ratio_%s", iso.c_str()), "");
        }
    }
    h_iss_ratios["Total"] = createHist("h_iss_ratio_total", "");

    // Get MC Ratios (returns TH1D* now)
    std::map<std::string, TH1D*> h_mc_ratios = calculateMCRatios(config, chain, refHist, is_debug_target);

    int bin_count_in_range = 0;
    
    for (int i_bin = 1; i_bin <= nBins; ++i_bin) {
        double binCenter = xAxis->GetBinCenter(i_bin);
        if (binCenter > 21.5 || binCenter < 0.2) continue;
        
        bin_count_in_range++;
        
        std::string det = (binCenter < 1.28) ? "TOF" : (binCenter < 3.06) ? "NaF" : "AGL";

        ValueWithError rawFragCounts_vw(0,0), final_source_yield_vw(0,0), final_contam_yield_vw(0,0);

        // Safe extraction from TH1* (supports F and D via GetBinContent)
        if (rawFragCountsHists.count(det) && rawFragCountsHists[det]) 
            rawFragCounts_vw = {rawFragCountsHists[det]->GetBinContent(i_bin), rawFragCountsHists[det]->GetBinError(i_bin)};
        
        if (l1SignalYieldHists.count(det) && l1SignalYieldHists[det])
            final_source_yield_vw = {l1SignalYieldHists[det]->GetBinContent(i_bin), l1SignalYieldHists[det]->GetBinError(i_bin)};
            
        if (l1ContamYieldHists.count(det) && l1ContamYieldHists[det])
            final_contam_yield_vw = {l1ContamYieldHists[det]->GetBinContent(i_bin), l1ContamYieldHists[det]->GetBinError(i_bin)};

        // Logic change: Denominator is simply the signal yield obtained from file
        if (final_source_yield_vw.val <= 0) continue;

        if (is_debug_target) {
            std::cout << std::string(125, '-') << "\n";
            std::cout << "Bin " << i_bin << " (Point " << bin_count_in_range << ", Center: " 
                     << std::fixed << std::setprecision(4) << binCenter << " GeV, Det: " << det << ")\n\n";
            std::cout << "--- ISS CALCULATION ---\n";
            std::cout << "Final Source (Denominator): " << final_source_yield_vw.val 
                     << " +/- " << final_source_yield_vw.err << "\n";
            std::cout << "Final Contamination (L1): " << final_contam_yield_vw.val 
                     << " +/- " << final_contam_yield_vw.err << "\n\n";
        }

        ValueWithError pure_total_frag_bin(0,0);
        if (config.InIsoLevel) {
            std::map<std::string, ValueWithError> frag_fracs_vw, norm_fracs_vw;
            ValueWithError frag_frac_sum(0,0), norm_frac_sum(0,0);

            auto get_fractions = [&](TFile* file, const std::vector<std::string>& isotopes, 
                                    std::map<std::string, ValueWithError>& vw_map, ValueWithError& sum_vw) {
                for (const auto& iso : isotopes) {
                    // Use TH1 for generic read
                    auto h = getHistFromFile<TH1>(file, TString::Format("h_best_%s_frac_%s", iso.c_str(), det.c_str()));
                    if(h) {
                        vw_map[iso] = {h->GetBinContent(i_bin), h->GetBinError(i_bin)};
                        sum_vw = sum_vw + vw_map[iso];
                        delete h;
                    }
                }
                const auto& last_iso = config.fragmentIsotopes.back();
                vw_map[last_iso] = ValueWithError(1.0, 0.0) - sum_vw;
            };
            
            get_fractions(fragmentFitFile.get(), config.primaryFitIsotopes, frag_fracs_vw, frag_frac_sum);
            get_fractions(normalFitFile.get(), config.primaryFitIsotopes, norm_fracs_vw, norm_frac_sum);
            
            for (const auto& iso : config.fragmentIsotopes) {
                // Calculation: (Raw_L2 * Frac) - (L1_Contam * Norm_Frac) -> All normalized by L1_Signal
                ValueWithError raw_iso_counts_vw = rawFragCounts_vw * frag_fracs_vw[iso];
                ValueWithError contamination_term_for_iso = final_contam_yield_vw * norm_fracs_vw[iso];
                ValueWithError pure_iso_counts_vw = raw_iso_counts_vw - contamination_term_for_iso;
                ValueWithError ratio_vw = pure_iso_counts_vw / final_source_yield_vw;
                fillHistBin(h_iss_ratios[iso], i_bin, ratio_vw);
                pure_total_frag_bin = pure_total_frag_bin + pure_iso_counts_vw;
            }
        } else {
            // If not doing isotopic analysis, Total Fragment - L1 Contamination
            pure_total_frag_bin = rawFragCounts_vw - final_contam_yield_vw;
        }

        ValueWithError total_ratio_vw = pure_total_frag_bin / final_source_yield_vw;
        fillHistBin(h_iss_ratios["Total"], i_bin, total_ratio_vw);
    }
    
    outputFile->cd();
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    std::string plot_dir = "/eos/user/z/zixuan/Isotope/BkgValid/Plots/";
    
    int bin_count = 0;
    for (int i_bin = 1; i_bin <= nBins; ++i_bin) {
        double binCenter = xAxis->GetBinCenter(i_bin);
        if (binCenter > 21.5 || binCenter < 0.4) continue;
        
        bin_count++;
        
        if (config.InIsoLevel) {
            for (const auto& iso : config.fragmentIsotopes) {
                if (shouldZeroPoint(config, bin_count, iso)) {
                    if (h_mc_ratios.count(iso) && h_mc_ratios[iso]) {
                        h_mc_ratios[iso]->SetBinContent(i_bin, 0);
                        h_mc_ratios[iso]->SetBinError(i_bin, 0);
                        std::cout << "[INFO] Zeroing MC " << iso << " point " << bin_count 
                                 << " (bin " << i_bin << ", E=" << binCenter << " GeV)" << std::endl;
                    }
                }
            }
        }
        
        if (shouldZeroPoint(config, bin_count, "Total")) {
            if (h_mc_ratios.count("Total") && h_mc_ratios["Total"]) {
                h_mc_ratios["Total"]->SetBinContent(i_bin, 0);
                h_mc_ratios["Total"]->SetBinError(i_bin, 0);
                std::cout << "[INFO] Zeroing MC Total point " << bin_count 
                         << " (bin " << i_bin << ", E=" << binCenter << " GeV)" << std::endl;
            }
        }
    }
    
    if (config.InIsoLevel) {
        for (const auto& iso : config.fragmentIsotopes) {
            h_iss_ratios[iso]->Write();
            if (h_mc_ratios.count(iso)) h_mc_ratios[iso]->Write();
            
            createComparisonPlots(h_iss_ratios[iso], h_mc_ratios.count(iso) ? h_mc_ratios[iso] : nullptr,
                TString::Format("%s#rightarrow^{%d}%s (%s)", config.sourceParticle.c_str(), 
                    config.fragZandMass.at(iso), config.fragmentParticle.c_str(), chain.c_str()).Data(),
                TString::Format("L2 frag ^{%d}%s / L1%s ", config.fragZandMass.at(iso), 
                    config.fragmentParticle.c_str(), config.sourceParticle.c_str()).Data(),
                plot_dir + TString::Format("%s_to_%s_%s", config.sourceParticle.c_str(), iso.c_str(), chain.c_str()).Data());
        }
    }
    
    h_iss_ratios["Total"]->Write();
    if (h_mc_ratios.count("Total")) h_mc_ratios["Total"]->Write();
    
    createComparisonPlots(h_iss_ratios["Total"], h_mc_ratios.count("Total") ? h_mc_ratios["Total"] : nullptr,
        TString::Format("%s#rightarrow%s (%s)", config.sourceParticle.c_str(), 
            config.fragmentParticle.c_str(), chain.c_str()).Data(),
        TString::Format("L2 frag %s / L1 %s", config.fragmentParticle.c_str(), 
            config.sourceParticle.c_str()).Data(),
        plot_dir + TString::Format("%s_to_%s_Total_%s", config.sourceParticle.c_str(), 
            config.fragmentParticle.c_str(), chain.c_str()).Data());

    std::cout << "Analysis complete. Results saved to " << outputFile->GetName() << std::endl;
    
    for(auto const& [key, val] : h_iss_ratios) delete val;
    for(auto const& [key, val] : h_mc_ratios) delete val;
}

// --- Main entry point ---
void CalFrag() {
    AnalysisConfig config_B_to_Be;
    config_B_to_Be.sourceParticle = "Boron";
    config_B_to_Be.fragmentParticle = "Beryllium";
    config_B_to_Be.InIsoLevel = true;
    config_B_to_Be.primaryFitIsotopes = {"Be7", "Be9"};
    config_B_to_Be.useMass = 7;
    config_B_to_Be.fragFileID = "Be_frag4";
    config_B_to_Be.sourceIsotopeFractions = {{"B10", 0.3}, {"B11", 0.7}};
    config_B_to_Be.fragmentIsotopes = {"Be7", "Be9", "Be10"};
    config_B_to_Be.fragZandMass = {{"Be7", 7}, {"Be9", 9}, {"Be10", 10}};
    config_B_to_Be.mcSourceFiles = {{"B10", "B10_rew_frag4.root"}, {"B11", "B11_rew_frag4.root"}};
    
    AnalysisConfig config_C_to_Be;
    config_C_to_Be.sourceParticle = "Carbon";
    config_C_to_Be.fragmentParticle = "Beryllium";
    config_C_to_Be.InIsoLevel = true;
    config_C_to_Be.primaryFitIsotopes = {"Be7", "Be9"};
    config_C_to_Be.useMass = 7;
    config_C_to_Be.fragFileID = "Be_frag4";
    config_C_to_Be.sourceIsotopeFractions = {{"C12", 1.0}};
    config_C_to_Be.fragmentIsotopes = {"Be7", "Be9", "Be10"};
    config_C_to_Be.fragZandMass = {{"Be7", 7}, {"Be9", 9}, {"Be10", 10}};
    config_C_to_Be.mcSourceFiles = {{"C12", "C12_rew_frag4.root"}};
    
    AnalysisConfig config_N_to_Be;
    config_N_to_Be.sourceParticle = "Nitrogen";
    config_N_to_Be.fragmentParticle = "Beryllium";
    config_N_to_Be.InIsoLevel = true;
    config_N_to_Be.primaryFitIsotopes = {"Be7", "Be9"};
    config_N_to_Be.useMass = 7;
    config_N_to_Be.fragFileID = "Be_frag4";
    config_N_to_Be.sourceIsotopeFractions = {{"N14", 0.5}, {"N15", 0.5}};
    config_N_to_Be.fragmentIsotopes = {"Be7", "Be9", "Be10"};
    config_N_to_Be.fragZandMass = {{"Be7", 7}, {"Be9", 9}, {"Be10", 10}};
    config_N_to_Be.mcSourceFiles = {{"N14", "N14_rew_frag4.root"}, {"N15", "N15_rew_frag4.root"}};
    
    AnalysisConfig config_O_to_Be;
    config_O_to_Be.sourceParticle = "Oxygen";
    config_O_to_Be.fragmentParticle = "Beryllium";
    config_O_to_Be.InIsoLevel = true;
    config_O_to_Be.primaryFitIsotopes = {"Be7", "Be9"};
    config_O_to_Be.useMass = 7;
    config_O_to_Be.fragFileID = "Be_frag4";
    config_O_to_Be.sourceIsotopeFractions = {{"O16", 1.0}};
    config_O_to_Be.fragmentIsotopes = {"Be7", "Be9", "Be10"};
    config_O_to_Be.fragZandMass = {{"Be7", 7}, {"Be9", 9}, {"Be10", 10}};
    config_O_to_Be.mcSourceFiles = {{"O16", "O16_rew_frag4.root"}};

    AnalysisConfig config_C_to_B;
    config_C_to_B.sourceParticle = "Carbon";
    config_C_to_B.fragmentParticle = "Boron";
    config_C_to_B.InIsoLevel = true;
    config_C_to_B.primaryFitIsotopes = {"B10"};
    config_C_to_B.useMass = 10;
    config_C_to_B.fragFileID = "B_frag5";
    config_C_to_B.sourceIsotopeFractions = {{"C12", 1.0}};
    config_C_to_B.fragmentIsotopes = {"B10", "B11"};
    config_C_to_B.fragZandMass = {{"B10", 10}, {"B11", 11}};
    config_C_to_B.mcSourceFiles = {{"C12", "C12_rew_frag5.root"}};
    
    AnalysisConfig config_N_to_B;
    config_N_to_B.sourceParticle = "Nitrogen";
    config_N_to_B.fragmentParticle = "Boron";
    config_N_to_B.InIsoLevel = true;
    config_N_to_B.primaryFitIsotopes = {"B10"};
    config_N_to_B.useMass = 10;
    config_N_to_B.fragFileID = "B_frag5";
    config_N_to_B.sourceIsotopeFractions = {{"N14", 0.5}, {"N15", 0.5}};
    config_N_to_B.fragmentIsotopes = {"B10", "B11"};
    config_N_to_B.fragZandMass = {{"B10", 10}, {"B11", 11}};
    config_N_to_B.mcSourceFiles = {{"N14", "N14_rew_frag5.root"}, {"N15", "N15_rew_frag5.root"}};
    
    AnalysisConfig config_O_to_B;
    config_O_to_B.sourceParticle = "Oxygen";
    config_O_to_B.fragmentParticle = "Boron";
    config_O_to_B.InIsoLevel = true;
    config_O_to_B.primaryFitIsotopes = {"B10"};
    config_O_to_B.useMass = 10;
    config_O_to_B.fragFileID = "B_frag5";
    config_O_to_B.sourceIsotopeFractions = {{"O16", 1.0}};
    config_O_to_B.fragmentIsotopes = {"B10", "B11"};
    config_O_to_B.fragZandMass = {{"B10", 10}, {"B11", 11}};
    config_O_to_B.mcSourceFiles = {{"O16", "O16_rew_frag5.root"}};
    
    try {
        const std::vector<std::string> chains = {"UnbiasedL1Inner"};
        const std::vector<AnalysisConfig> all_configs = {
            config_B_to_Be, config_C_to_Be,
            config_N_to_Be, config_O_to_Be
        };
        
        for (const auto& chain : chains) {
            for (const auto& config : all_configs) {
                runAnalysis(chain, config);
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "[FATAL ERROR] An exception occurred: " << e.what() << std::endl;
    }
}