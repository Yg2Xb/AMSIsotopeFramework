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

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <map>
#include <string>
#include <numeric>

// --- Helper: ValueWithError (Error Propagation) ---
// (这部分代码与之前版本完全相同，保持不变)
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
    ValueWithError operator*(double scalar) const {
        return {val * scalar, std::abs(err * scalar)};
    }
    ValueWithError operator*(const ValueWithError& o) const {
        double p_val = val * o.val;
        if (val == 0.0 || o.val == 0.0) return {p_val, 0.0};
        double rel_err_sq = (err/val)*(err/val) + (o.err/o.val)*(o.err/o.val);
        return {p_val, std::abs(p_val) * std::sqrt(rel_err_sq)};
    }
    ValueWithError operator/(const ValueWithError& o) const {
        if (o.val == 0.0) {
            return {std::nan(""), std::nan("")};
        }
        if (val == 0.0) return {0.0, 0.0};
        double q_val = val / o.val;
        double rel_err_sq = (err/val)*(err/val) + (o.err/o.val)*(o.err/o.val);
        return {q_val, std::abs(q_val) * std::sqrt(rel_err_sq)};
    }
};
ValueWithError operator*(double scalar, const ValueWithError& ve) {
    return ve * scalar;
}

// --- Helper: File and Histogram Operations ---
std::unique_ptr<TFile> openFile(const TString& filename, const char* option = "READ") {
    auto file = std::unique_ptr<TFile>(TFile::Open(filename.Data(), option));
    if (!file || file->IsZombie()) {
        throw std::runtime_error("Error opening file: " + std::string(filename.Data()));
    }
    std::cout << "[INFO] Opened file: " << filename << std::endl;
    return file;
}

template<typename HistType>
HistType* getHistFromFile(TFile* file, const TString& histName) {
    if (!file) {
         throw std::runtime_error("Input TFile pointer is null when trying to get histogram: " + std::string(histName.Data()));
    }
    HistType* hist = nullptr;
    file->GetObject(histName, hist);
    if (!hist) {
        throw std::runtime_error("Could not find histogram: " + std::string(histName.Data()) + " in file " + file->GetName());
    }
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
    std::vector<std::string> fragmentIsotopes;
    std::vector<std::string> primaryFitIsotopes;
    int useMass;
    std::string countsFileID;
    std::map<std::string, double> sourceIsotopeFractions;
};

// --- Main Calculation Function (Refactored for Dynamic Binning) ---
void runAnalysis(const std::string& chain, const AnalysisConfig& config) {
    std::cout << "\n===================================================================\n";
    std::cout << "Starting Analysis: " << config.sourceParticle << " -> " << config.fragmentParticle
              << " for chain: " << chain << std::endl;
    std::cout << "===================================================================\n";

    const std::array<std::string, 3> detectors = {"TOF", "NaF", "AGL"};

    // --- Dynamically build file paths ---
    TString fragmentFitFilePath = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d_FragFrom%s.root",
        config.fragmentParticle.c_str(), chain.c_str(), config.useMass, config.sourceParticle.c_str());
    TString normalFitFilePath = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d.root",
        config.fragmentParticle.c_str(), chain.c_str(), config.useMass);
    TString countsFilePath = TString::Format("/eos/user/z/zixuan/Isotope/Add/%s.root", config.countsFileID.c_str());
    TString yieldFilePath = TString::Format("/eos/user/z/zixuan/Isotope/ChargeTemp/QFit_%s_to_%s_%s.root",
        config.sourceParticle.c_str(), config.fragmentParticle.c_str(), chain.c_str());
    TString outputFilePath = TString::Format("/eos/user/z/zixuan/Isotope/BkgValid/%s_to_%s_%s_Validation.root",
        config.sourceParticle.c_str(), config.fragmentParticle.c_str(), chain.c_str());

    // --- Open all files ---
    auto fragmentFitFile = openFile(fragmentFitFilePath);
    auto normalFitFile   = openFile(normalFitFilePath);
    auto countsFile      = openFile(countsFilePath);
    auto yieldFile       = openFile(yieldFilePath);
    auto outputFile      = openFile(outputFilePath, "RECREATE");

    // --- Get Binning Information Dynamically ---
    // We use one of the yield histograms as the reference for our "wide bin" structure.
    TH1D* refHist = getHistFromFile<TH1D>(yieldFile.get(), TString::Format("h_yield_%s_TOF", config.sourceParticle.c_str()));
    const TAxis* xAxis = refHist->GetXaxis();
    const int nBins = xAxis->GetNbins();
    const Double_t* binEdges = xAxis->GetXbins()->GetArray(); // Get the array of bin edges
    std::cout << "[INFO] Dynamic binning retrieved: " << nBins << " bins from " << binEdges[0] << " to " << binEdges[nBins] << std::endl;

    // --- Prepare input histograms ---
    std::map<std::string, TH1F*> fragFracHists, normFracHists;
    std::map<std::string, TH1D*> rawFragCountsHists, sourceYieldHists, contamYieldHists;

    for (const auto& det : detectors) {
        rawFragCountsHists[det] = getHistFromFile<TH1D>(countsFile.get(), TString::Format("%s_ISS_BKG_H1_%s_%s", chain.c_str(), config.sourceParticle.c_str(), det.c_str()));
        rawFragCountsHists[det]->Rebin(2);
        sourceYieldHists[det] = getHistFromFile<TH1D>(yieldFile.get(), TString::Format("h_yield_%s_%s", config.sourceParticle.c_str(), det.c_str()));
        contamYieldHists[det] = getHistFromFile<TH1D>(yieldFile.get(), TString::Format("h_yield_%s_%s", config.fragmentParticle.c_str(), det.c_str()));

        for (const auto& iso : config.primaryFitIsotopes) {
            fragFracHists[iso + "_" + det] = getHistFromFile<TH1F>(fragmentFitFile.get(), TString::Format("h_best_%s_frac_%s", iso.c_str(), det.c_str()));
            normFracHists[iso + "_" + det] = getHistFromFile<TH1F>(normalFitFile.get(), TString::Format("h_best_%s_frac_%s", iso.c_str(), det.c_str()));
        }
        if (config.fragmentParticle == "Boron") {
             fragFracHists["B10_" + det] = getHistFromFile<TH1F>(fragmentFitFile.get(), TString::Format("h_best_B10_frac_%s", det.c_str()));
             normFracHists["B10_" + det] = getHistFromFile<TH1F>(normalFitFile.get(), TString::Format("h_best_B10_frac_%s", det.c_str()));
        }
    }

    // --- Prepare output histograms using the dynamic binning ---
    std::vector<TObject*> allObjectsToSave;
    auto createHist = [&](const char* name, const char* title) {
        // Use the dynamically retrieved nBins and binEdges
        auto hist = new TH1F(name, title, nBins, binEdges);
        hist->SetStats(0);
        allObjectsToSave.push_back(hist);
        return hist;
    };

    std::map<std::string, TH1F*> h_pure_counts, h_ratio_vs_total_source;
    std::map<std::string, std::map<std::string, TH1F*>> h_ratio_vs_source_iso;

    for (const auto& iso : config.fragmentIsotopes) {
        h_pure_counts[iso] = createHist(TString::Format("h_pure_counts_%s", iso.c_str()), TString::Format("Pure Fragment Counts of ^{%s};E_{k}/n [GeV/n];Counts", iso.c_str()));
        h_ratio_vs_total_source[iso] = createHist(TString::Format("h_ratio_%s_vs_Total_%s", iso.c_str(), config.sourceParticle.c_str()), TString::Format("Ratio: Pure ^{%s} / Total %s;E_{k}/n [GeV/n];Ratio", iso.c_str(), config.sourceParticle.c_str()));
        
        if (!config.sourceIsotopeFractions.empty()) {
            for (const auto& pair : config.sourceIsotopeFractions) {
                const std::string& source_iso = pair.first;
                h_ratio_vs_source_iso[iso][source_iso] = createHist(TString::Format("h_ratio_%s_vs_%s", iso.c_str(), source_iso.c_str()), TString::Format("Ratio: Pure ^{%s} / Source ^{%s};E_{k}/n [GeV/n];Ratio", iso.c_str(), source_iso.c_str()));
            }
        }
    }
    TH1F* h_total_raw_frag = createHist("h_total_raw_fragment_counts", TString::Format("Total Raw Counts of Selected %s;E_{k}/n [GeV/n];Counts", config.fragmentParticle.c_str()));
    TH1F* h_total_contam = createHist("h_total_contamination_counts", TString::Format("Total Contamination Counts from Primary %s;E_{k}/n [GeV/n];Counts", config.fragmentParticle.c_str()));
    TH1F* h_total_source = createHist("h_total_source_counts", TString::Format("Total Source Counts of %s;E_{k}/n [GeV/n];Counts", config.sourceParticle.c_str()));

    // --- Main loop over energy bins ---
    for (int i_bin = 1; i_bin <= nBins; ++i_bin) {
        // Get bin center directly from the reference axis
        double binCenter = xAxis->GetBinCenter(i_bin);
        if(binCenter > 22) break; // Limit to 22 GeV/n as before
        std::string det = (binCenter < 1.28) ? "TOF" : (binCenter < 3.06) ? "NaF" : "AGL";


        ValueWithError rawFragCounts_vw(rawFragCountsHists[det]->GetBinContent(i_bin), rawFragCountsHists[det]->GetBinError(i_bin));
        ValueWithError sourceYield_vw(sourceYieldHists[det]->GetBinContent(i_bin), sourceYieldHists[det]->GetBinError(i_bin));
        ValueWithError contamYield_vw(contamYieldHists[det]->GetBinContent(i_bin), contamYieldHists[det]->GetBinError(i_bin));

        std::map<std::string, ValueWithError> frag_fracs_vw, norm_fracs_vw;
        ValueWithError frag_frac_sum(0,0), norm_frac_sum(0,0);

        auto get_fractions = [&](std::map<std::string, TH1F*>& hist_map, std::map<std::string, ValueWithError>& vw_map, ValueWithError& sum_vw) {
            for (const auto& iso : config.primaryFitIsotopes) {
                double val = hist_map[iso + "_" + det]->GetBinContent(i_bin);
                double err = hist_map[iso + "_" + det]->GetBinError(i_bin);
                vw_map[iso] = {val, err};
                sum_vw = sum_vw + vw_map[iso];
            }
            const auto& last_iso = config.fragmentIsotopes.back();
            vw_map[last_iso] = ValueWithError(1.0, 0.0) - sum_vw;
        };
        
        get_fractions(fragFracHists, frag_fracs_vw, frag_frac_sum);
        get_fractions(normFracHists, norm_fracs_vw, norm_frac_sum);

        // --- Perform calculations (This core logic is unchanged) ---
        for (const auto& iso : config.fragmentIsotopes) {
            ValueWithError raw_iso_counts_vw = rawFragCounts_vw * frag_fracs_vw[iso];
            ValueWithError contamination_term_for_iso = contamYield_vw * norm_fracs_vw[iso];
            ValueWithError pure_iso_counts_vw = raw_iso_counts_vw - contamination_term_for_iso;
            fillHistBin(h_pure_counts[iso], i_bin, pure_iso_counts_vw);

            ValueWithError ratio_vs_total_source_vw = pure_iso_counts_vw / sourceYield_vw;
            fillHistBin(h_ratio_vs_total_source[iso], i_bin, ratio_vs_total_source_vw);

            if (!config.sourceIsotopeFractions.empty()) {
                for (const auto& pair : config.sourceIsotopeFractions) {
                    const std::string& source_iso = pair.first;
                    double source_iso_frac = pair.second;
                    ValueWithError source_iso_yield_vw = sourceYield_vw * source_iso_frac;
                    ValueWithError ratio_vs_source_iso_vw = pure_iso_counts_vw / source_iso_yield_vw;
                    fillHistBin(h_ratio_vs_source_iso[iso][source_iso], i_bin, ratio_vs_source_iso_vw);
                }
            }
        }
        
        fillHistBin(h_total_raw_frag, i_bin, rawFragCounts_vw);
        fillHistBin(h_total_contam, i_bin, contamYield_vw);
        fillHistBin(h_total_source, i_bin, sourceYield_vw);
    }

    // --- Save all histograms and canvases to file (Unchanged) ---
    outputFile->cd();
    for (TObject* obj : allObjectsToSave) {
        obj->Write(obj->GetName(), TObject::kOverwrite);
    }
    
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);

    auto createCanvas = [&](const char* name, const char* title, int divX, int divY) {
        auto canvas = new TCanvas(name, title, 400 * divX, 350 * divY);
        canvas->Divide(divX, divY);
        return canvas;
    };

    std::map<std::string, EColor> iso_colors = {
        {"Be7", kGreen+2}, {"Be9", kBlue}, {"Be10", kMagenta},
        {"B10", kOrange+1}, {"B11", kCyan+2}
    };

    TCanvas* c_pure = createCanvas("c_pure_counts", "Pure Fragment Counts", 1, config.fragmentIsotopes.size());
    for (size_t i = 0; i < config.fragmentIsotopes.size(); ++i) {
        const auto& iso = config.fragmentIsotopes[i];
        c_pure->cd(i + 1);
        gPad->SetLogy();
        h_pure_counts[iso]->SetLineColor(iso_colors[iso]);
        h_pure_counts[iso]->SetMarkerColor(iso_colors[iso]);
        h_pure_counts[iso]->SetMarkerStyle(20);
        h_pure_counts[iso]->Draw("E1");
    }
    c_pure->Write();

    if (!config.sourceIsotopeFractions.empty()) {
        for (const auto& pair : config.sourceIsotopeFractions) {
            const std::string& source_iso = pair.first;
            TCanvas* c_ratio_iso = createCanvas(TString::Format("c_ratio_vs_%s", source_iso.c_str()), TString::Format("Ratio vs Source %s", source_iso.c_str()), 1, config.fragmentIsotopes.size());
            for (size_t i = 0; i < config.fragmentIsotopes.size(); ++i) {
                const auto& frag_iso = config.fragmentIsotopes[i];
                c_ratio_iso->cd(i + 1);
                gPad->SetLogy();
                h_ratio_vs_source_iso[frag_iso][source_iso]->SetLineColor(iso_colors[frag_iso]);
                h_ratio_vs_source_iso[frag_iso][source_iso]->SetMarkerColor(iso_colors[frag_iso]);
                h_ratio_vs_source_iso[frag_iso][source_iso]->SetMarkerStyle(20);
                h_ratio_vs_source_iso[frag_iso][source_iso]->Draw("E1");
            }
            c_ratio_iso->Write();
        }
    } else {
        TCanvas* c_ratio_total = createCanvas("c_ratio_vs_total_source", TString::Format("Ratio vs Total %s", config.sourceParticle.c_str()), 1, config.fragmentIsotopes.size());
        for (size_t i = 0; i < config.fragmentIsotopes.size(); ++i) {
            const auto& iso = config.fragmentIsotopes[i];
            c_ratio_total->cd(i + 1);
            gPad->SetLogy();
            h_ratio_vs_total_source[iso]->SetLineColor(iso_colors[iso]);
            h_ratio_vs_total_source[iso]->SetMarkerColor(iso_colors[iso]);
            h_ratio_vs_total_source[iso]->SetMarkerStyle(20);
            h_ratio_vs_total_source[iso]->Draw("E1");
        }
        c_ratio_total->Write();
    }
    
    std::cout << "Analysis complete. Results saved to " << outputFile->GetName() << std::endl;
}

// --- Main entry point (Unchanged) ---
void CalFrag() {
    AnalysisConfig config_B_to_Be;
    config_B_to_Be.sourceParticle = "Boron";
    config_B_to_Be.fragmentParticle = "Beryllium";
    config_B_to_Be.fragmentIsotopes = {"Be7", "Be9", "Be10"};
    config_B_to_Be.primaryFitIsotopes = {"Be7", "Be9"};
    config_B_to_Be.useMass = 7;
    config_B_to_Be.countsFileID = "Be_frag4";
    config_B_to_Be.sourceIsotopeFractions = {{"B10", 0.3}, {"B11", 0.7}};

    AnalysisConfig config_C_to_B;
    config_C_to_B.sourceParticle = "Carbon";
    config_C_to_B.fragmentParticle = "Boron";
    config_C_to_B.fragmentIsotopes = {"B10", "B11"};
    config_C_to_B.primaryFitIsotopes = {"B10"};
    config_C_to_B.useMass = 10;
    config_C_to_B.countsFileID = "B_frag5";
    
    try {
        const std::vector<std::string> chains = {"L1Inner", "UnbiasedL1Inner"};
        for (const auto& chain : chains) {
            runAnalysis(chain, config_B_to_Be);
            runAnalysis(chain, config_C_to_B);
        }
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] An exception occurred: " << e.what() << std::endl;
    }
}