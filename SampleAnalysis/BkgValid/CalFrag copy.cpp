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
const bool DEBUG_MODE = true; // MODIFIED: Set back to true to enable all debug outputs
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

template<typename HistType>
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
    bool isAdjacentZ;
    std::vector<std::string> primaryFitIsotopes;
    int useMass;
    std::string fragFileID;
    std::map<std::string, double> sourceIsotopeFractions;
    std::vector<std::string> fragmentIsotopes;
    std::map<std::string, int> fragZandMass;
    std::map<std::string, std::string> mcSourceFiles;
    std::string mcSourceGenHist = "MC_FLUX_H3";
};


// --- Plotting Helper ---
void createComparisonPlot(TH1* h_iss, TH1* h_mc, const std::string& title, const std::string& y_axis_title, const std::string& output_path) {
    if (!h_iss) {
        std::cerr << "[PLOT ERROR] ISS histogram is null for plot: " << title << std::endl;
        return;
    }
    TCanvas* c = new TCanvas(TString::Format("c_%s", output_path.c_str()), title.c_str(), 1000, 800);
    TPad* top_pad = new TPad("top_pad", "top_pad", 0, 0.3, 1, 1.0);
    top_pad->SetBottomMargin(0.02); top_pad->SetTopMargin(0.08); top_pad->SetGridx(); top_pad->SetGridy();
    top_pad->Draw();
    TPad* bot_pad = new TPad("bot_pad", "bot_pad", 0, 0.05, 1, 0.3);
    bot_pad->SetTopMargin(0.05); bot_pad->SetBottomMargin(0.3); bot_pad->SetGridx(); bot_pad->SetGridy();
    bot_pad->Draw();
    top_pad->cd();
    h_iss->SetMarkerStyle(20); h_iss->SetMarkerColor(kBlack); h_iss->SetLineColor(kBlack);
    h_iss->SetTitle(title.c_str());
    h_iss->GetYaxis()->SetTitle(y_axis_title.c_str());
    h_iss->GetXaxis()->SetLabelSize(0); h_iss->GetXaxis()->SetTitleSize(0);
    double max_y = 0;
    for(int i = 1; i <= h_iss->GetNbinsX(); ++i) {
        if (h_iss->GetBinCenter(i) > 21 || h_iss->GetBinCenter(i) < 0.3) continue;
        max_y = std::max(max_y, h_iss->GetBinContent(i) + h_iss->GetBinError(i));
    }
    if (h_mc) {
        for(int i = 1; i <= h_mc->GetNbinsX(); ++i) {
            if (h_mc->GetBinCenter(i) > 21 || h_mc->GetBinCenter(i) < 0.3) continue;
            max_y = std::max(max_y, h_mc->GetBinContent(i) + h_mc->GetBinError(i));
        }
    }
    h_iss->GetYaxis()->SetRangeUser(0, max_y * 1.4);
    h_iss->GetXaxis()->SetRangeUser(0.3, 21);
    h_iss->GetYaxis()->SetTitleOffset(1.4);
    h_iss->Draw("PZ");
    TLegend* legend = new TLegend(0.55, 0.70, 0.88, 0.88);
    legend->SetBorderSize(0); legend->SetFillStyle(0);
    legend->AddEntry(h_iss, "ISS Data", "ep");
    if (h_mc) {
        h_mc->SetMarkerStyle(20); h_mc->SetMarkerColor(kRed); h_mc->SetLineColor(kRed);
        h_mc->Draw("PZ SAME");
        legend->AddEntry(h_mc, "MC Simulation", "ep");
    }
    legend->Draw();
    bot_pad->cd();
    if (h_mc) {
        TH1F* h_ratio = (TH1F*)h_iss->Clone("h_ratio");
        h_ratio->Divide(h_mc);
        h_ratio->SetTitle("");

        double min_r = 1e9, max_r = -1e9;
        for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
            if (h_ratio->GetBinCenter(i) > 21 || h_ratio->GetBinCenter(i) < 0.3) continue;
            if (h_ratio->GetBinContent(i) != 0) {
                min_r = std::min(min_r, h_ratio->GetBinContent(i) - h_ratio->GetBinError(i));
                max_r = std::max(max_r, h_ratio->GetBinContent(i) + h_ratio->GetBinError(i));
            }
        }
        if (min_r > max_r) { min_r = 0.2; max_r = 1.8; }

        h_ratio->GetYaxis()->SetTitle("ISS / MC");
        h_ratio->GetYaxis()->SetRangeUser(min_r * 0.7, max_r * 1.);
        h_ratio->GetYaxis()->SetNdivisions(504);
        h_ratio->GetYaxis()->CenterTitle();
        h_ratio->GetYaxis()->SetTitleSize(0.1);
        h_ratio->GetYaxis()->SetLabelSize(0.12);
        h_ratio->GetYaxis()->SetTitleOffset(0.5);

        h_ratio->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetLabelSize(0.12);
        h_ratio->GetXaxis()->SetTitleOffset(1.0);
        h_ratio->GetXaxis()->SetRangeUser(0.3, 21);
        h_ratio->Draw("PZ");
        
        TLine* unityLine = new TLine(0.3, 1.0,21.5, 1.0);
        unityLine->SetLineStyle(2);
        unityLine->SetLineColor(kGray + 2);
        unityLine->SetLineWidth(2);
        unityLine->Draw("SAME");
    }
    c->SaveAs((output_path + ".png").c_str());
    delete c;
}


// --- MC Calculation Helpers ---
TH1F* stitchHistograms(const std::map<std::string, TH1*>& hists, TH1* refHist) {
    auto combined = (TH1F*)refHist->Clone(TString::Format("%s_stitched", refHist->GetName()));
    combined->Reset();
    const double boundary1 = 1.28, boundary2 = 3.06;
    for (int i = 1; i <= combined->GetNbinsX(); ++i) {
        double binCenter = combined->GetBinCenter(i);
        std::string det = (binCenter < boundary1) ? "TOF" : (binCenter < boundary2) ? "NaF" : "AGL";
        if (hists.count(det) && hists.at(det) != nullptr) {
            TH1* h = hists.at(det);
            int source_bin = h->FindBin(binCenter);
            combined->SetBinContent(i, h->GetBinContent(source_bin));
            combined->SetBinError(i, h->GetBinError(source_bin));
        }
    }
    return combined;
}

std::map<std::string, TH1F*> calculateMCRatios(const AnalysisConfig& config, const std::string& chain, TH1* refHist, bool is_debug_target) {
    std::cout << "\n--- Calculating MC Ratios for " << config.sourceParticle << " -> " << config.fragmentParticle << " ---\n";
    std::map<std::string, TH1F*> mc_ratios;
    const std::string mcBaseDir = "/eos/user/z/zixuan/Isotope/Add/";
    const std::array<std::string, 3> detectors = {"TOF", "NaF", "AGL"};
    const int nBins = refHist->GetNbinsX();

    TH1F* h_mc_total_source = (TH1F*)refHist->Clone("h_mc_total_source_stitched"); h_mc_total_source->Reset();
    TH1F* h_mc_total_frag = (TH1F*)refHist->Clone("h_mc_total_frag_stitched"); h_mc_total_frag->Reset();
    std::map<std::string, TH1F*> h_mc_frag_isos;
    for (const auto& frag_iso : config.fragmentIsotopes) {
        h_mc_frag_isos[frag_iso] = (TH1F*)refHist->Clone(TString::Format("h_mc_%s_stitched", frag_iso.c_str()));
        h_mc_frag_isos[frag_iso]->Reset();
    }

    std::map<std::string, std::unique_ptr<TH1D>> h_gen_map;
    for (const auto& pair : config.mcSourceFiles) {
        const std::string& source_iso = pair.first;
        auto mc_file = openFile(mcBaseDir + pair.second);
        h_gen_map[source_iso].reset(getHistFromFile<TH1D>(mc_file.get(), config.mcSourceGenHist.c_str()));
        if(h_gen_map[source_iso]) h_gen_map[source_iso]->Rebin(2);
    }
    
    if (is_debug_target) {
        std::cout << "\n" << std::string(125, '*') << "\n";
        std::cout << "STARTING DETAILED MC BIN-BY-BIN CALCULATION\n";
        std::cout << std::string(125, '*') << "\n";
    }

    for (const auto& pair : config.mcSourceFiles) {
        const std::string& current_source_iso = pair.first;
        auto mc_file = openFile(mcBaseDir + pair.second);

        std::map<std::string, TH1*> source_hists_per_det;
        for(const auto& det : detectors) {
            source_hists_per_det[det] = getHistFromFile<TH1>(mc_file.get(), TString::Format("%s_MC_BKG_H1_%s", chain.c_str(), det.c_str()));
            if(source_hists_per_det[det]) source_hists_per_det[det]->Rebin(2);
        }
        std::unique_ptr<TH1F> stitched_source(stitchHistograms(source_hists_per_det, refHist));

        std::map<std::string, std::unique_ptr<TH1F>> stitched_frags;
        for (const auto& frag_iso : config.fragmentIsotopes) {
            std::map<std::string, TH1*> frag_hists_per_det;
            int frag_Z = (config.fragmentParticle == "Beryllium") ? 4 : 5;
            for(const auto& det : detectors) {
                TString histname = TString::Format("%s_MC_BKG_H2_%s_Z%d_Mass%d", chain.c_str(), det.c_str(), frag_Z, config.fragZandMass.at(frag_iso));
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
            double n_gen_current = (h_gen_map.count(current_source_iso) && h_gen_map.at(current_source_iso)) ? h_gen_map.at(current_source_iso)->GetBinContent(i_bin) : 0.0;
            double weight = 0.0;
            if (n_gen_current > 0) {
                weight = (total_gen_in_bin / n_gen_current) * config.sourceIsotopeFractions.at(current_source_iso);
            }

            ValueWithError raw_source_vw(stitched_source->GetBinContent(i_bin), stitched_source->GetBinError(i_bin));
            ValueWithError weighted_source_vw = raw_source_vw * weight;
            
            // MODIFIED: Re-introduced debug block
            if (is_debug_target) {
                if (i_bin == 1 || (refHist->GetBinCenter(i_bin-1) <= 21 && refHist->GetBinCenter(i_bin) < 21)) {
                    std::cout << std::string(125, '-') << "\n";
                    std::cout << "MC DEBUG FOR Bin " << i_bin << " (Center: " << std::fixed << std::setprecision(4) << refHist->GetBinCenter(i_bin) << " GeV/n)\n";
                    std::cout << "Contribution from: " << current_source_iso << "\n";
                    std::cout << "  - Weight Calc: N_gen(" << n_gen_current << "), N_total(" << total_gen_in_bin << "), Frac(" << config.sourceIsotopeFractions.at(current_source_iso) << ") ==> Weight: " << weight << "\n";
                    std::cout << "  - Source (H1): " << raw_source_vw.val << " -> " << weighted_source_vw.val << "\n";
                }
            }

            h_mc_total_source->SetBinContent(i_bin, h_mc_total_source->GetBinContent(i_bin) + weighted_source_vw.val);
            h_mc_total_source->SetBinError(i_bin, std::hypot(h_mc_total_source->GetBinError(i_bin), weighted_source_vw.err));
            
            for (const auto& frag_iso : config.fragmentIsotopes) {
                ValueWithError raw_frag_vw = {stitched_frags[frag_iso]->GetBinContent(i_bin), stitched_frags[frag_iso]->GetBinError(i_bin)};
                ValueWithError weighted_frag_vw = raw_frag_vw * weight;
                h_mc_frag_isos[frag_iso]->SetBinContent(i_bin, h_mc_frag_isos[frag_iso]->GetBinContent(i_bin) + weighted_frag_vw.val);
                h_mc_frag_isos[frag_iso]->SetBinError(i_bin, std::hypot(h_mc_frag_isos[frag_iso]->GetBinError(i_bin), weighted_frag_vw.err));
                
                h_mc_total_frag->SetBinContent(i_bin, h_mc_total_frag->GetBinContent(i_bin) + weighted_frag_vw.val);
                h_mc_total_frag->SetBinError(i_bin, std::hypot(h_mc_total_frag->GetBinError(i_bin), weighted_frag_vw.err));
            }
        }
    }

    if (config.isAdjacentZ) {
        for (const auto& frag_iso : config.fragmentIsotopes) {
            TH1F* ratio = (TH1F*)h_mc_frag_isos[frag_iso]->Clone(TString::Format("h_mc_ratio_%s", frag_iso.c_str()));
            ratio->Divide(h_mc_total_source);
            mc_ratios[frag_iso] = ratio;
        }
    }
    TH1F* total_ratio_hist = (TH1F*)h_mc_total_frag->Clone("h_mc_ratio_total");
    total_ratio_hist->Divide(h_mc_total_source);
    mc_ratios["Total"] = total_ratio_hist;

    delete h_mc_total_source; delete h_mc_total_frag;
    for(auto const& [key, val] : h_mc_frag_isos) delete val;
    return mc_ratios;
}


// --- Main Calculation Function ---
void runAnalysis(const std::string& chain, const AnalysisConfig& config) {
    // MODIFIED: This now correctly depends on the global DEBUG_MODE setting
    bool is_debug_target = true;//(DEBUG_MODE && config.sourceParticle == DEBUG_SOURCE && config.fragmentParticle == DEBUG_FRAGMENT && chain == DEBUG_CHAIN);

    std::cout << "\n===================================================================\n";
    std::cout << "Starting Analysis: " << config.sourceParticle << " -> " << config.fragmentParticle
              << " for chain: " << chain << " (Adjacent Z: " << (config.isAdjacentZ ? "Yes" : "No") << ")" << std::endl;
    std::cout << "===================================================================\n";

    const std::array<std::string, 3> detectors = {"TOF", "NaF", "AGL"};
    TString ratioFilePath = TString::Format("/eos/user/z/zixuan/Isotope/ChargeTemp/QFit_%s_to_%s_%s.root", config.sourceParticle.c_str(), config.fragmentParticle.c_str(), chain.c_str());
    TString countsFilePath = TString::Format("/eos/user/z/zixuan/Isotope/Add/%s.root", config.fragFileID.c_str());
    TString outputFilePath = TString::Format("/eos/user/z/zixuan/Isotope/BkgValid/%s_to_%s_%s_Validation.root", config.sourceParticle.c_str(), config.fragmentParticle.c_str(), chain.c_str());

    auto ratioFile = openFile(ratioFilePath);
    auto countsFile = openFile(countsFilePath);
    auto outputFile = openFile(outputFilePath, "RECREATE");

    TH1D* refHist = getHistFromFile<TH1D>(countsFile.get(), TString::Format("%s_ISS_BKG_H1_%s_TOF", chain.c_str(), config.sourceParticle.c_str()));
    if (!refHist) {
        throw std::runtime_error("Reference histogram for binning not found!");
    }
    refHist->Rebin(2);
    const TAxis* xAxis = refHist->GetXaxis();
    const int nBins = xAxis->GetNbins();
    const Double_t* binEdges = xAxis->GetXbins()->GetArray();

    std::map<std::string, TH1D*> l1SampleCountsHists, rawFragCountsHists, sourceRatioHists, contamRatioHists;
    for (const auto& det : detectors) {
        l1SampleCountsHists[det] = getHistFromFile<TH1D>(countsFile.get(), TString::Format("%s_ISS_BKG_H1_%s_%s", chain.c_str(), config.sourceParticle.c_str(), det.c_str()));
        if(l1SampleCountsHists[det]) l1SampleCountsHists[det]->Rebin(2);
        
        rawFragCountsHists[det] = getHistFromFile<TH1D>(countsFile.get(), TString::Format("%s_ISS_BKG_H3_%s_%s", chain.c_str(), config.sourceParticle.c_str(), det.c_str()));
        if(rawFragCountsHists[det]) rawFragCountsHists[det]->Rebin(2);
        
        sourceRatioHists[det] = getHistFromFile<TH1D>(ratioFile.get(), TString::Format("h_narrowfrac_%s_%s", config.sourceParticle.c_str(), det.c_str()));
        contamRatioHists[det] = getHistFromFile<TH1D>(ratioFile.get(), TString::Format("h_narrowfrac_%s_%s", config.fragmentParticle.c_str(), det.c_str()));
    }
    
    std::unique_ptr<TFile> fragmentFitFile, normalFitFile;
    if (config.isAdjacentZ) {
        std::string fragAbbr = particleNameToAbbr.at(config.fragmentParticle);
        TString fragmentFitFilePath = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d_FragFrom%s.root", fragAbbr.c_str(), chain.c_str(), config.useMass, config.sourceParticle.c_str());
        TString normalFitFilePath = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d.root", fragAbbr.c_str(), chain.c_str(), config.useMass);
        fragmentFitFile = openFile(fragmentFitFilePath);
        normalFitFile   = openFile(normalFitFilePath);
    }

    auto createHist = [&](const char* name, const char* title) {
        auto hist = new TH1F(name, title, nBins, binEdges); hist->SetStats(0); return hist;
    };
    std::map<std::string, TH1F*> h_iss_ratios;
    if (config.isAdjacentZ) {
        for (const auto& iso : config.fragmentIsotopes) { h_iss_ratios[iso] = createHist(TString::Format("h_iss_ratio_%s", iso.c_str()), ""); }
    }
    h_iss_ratios["Total"] = createHist("h_iss_ratio_total", "");

    std::map<std::string, TH1F*> h_mc_ratios = calculateMCRatios(config, chain, refHist, is_debug_target);

    for (int i_bin = 1; i_bin <= nBins; ++i_bin) {
        double binCenter = xAxis->GetBinCenter(i_bin);
        if (binCenter > 21) break;
        std::string det = (binCenter < 1.28) ? "TOF" : (binCenter < 3.06) ? "NaF" : "AGL";

        ValueWithError l1SampleCounts_vw(0,0), rawFragCounts_vw(0,0), sourceRatio_vw(0,0), contamRatio_vw(0,0);
        if (l1SampleCountsHists.count(det) && l1SampleCountsHists[det]) l1SampleCounts_vw = {l1SampleCountsHists[det]->GetBinContent(i_bin), l1SampleCountsHists[det]->GetBinError(i_bin)};
        if (rawFragCountsHists.count(det) && rawFragCountsHists[det]) rawFragCounts_vw = {rawFragCountsHists[det]->GetBinContent(i_bin), rawFragCountsHists[det]->GetBinError(i_bin)};
        if (sourceRatioHists.count(det) && sourceRatioHists[det]) sourceRatio_vw = {sourceRatioHists[det]->GetBinContent(i_bin), sourceRatioHists[det]->GetBinError(i_bin)};
        if (contamRatioHists.count(det) && contamRatioHists[det]) contamRatio_vw = {contamRatioHists[det]->GetBinContent(i_bin), contamRatioHists[det]->GetBinError(i_bin)};

        ValueWithError final_source_yield_vw = l1SampleCounts_vw * sourceRatio_vw;
        ValueWithError final_contam_yield_vw = l1SampleCounts_vw * contamRatio_vw;

        if (final_source_yield_vw.val <= 0) continue;

        // MODIFIED: Re-introduced debug block
        if (is_debug_target) {
            std::cout << std::string(125, '-') << "\n";
            std::cout << "Bin " << i_bin << " (Center: " << std::fixed << std::setprecision(4) << binCenter << " GeV/n, Det: " << det << ")\n\n";
            std::cout << "--- ISS CALCULATION ---\n";
            std::cout << std::left << std::setw(35) << "1. L1 Sample Counts:" << l1SampleCounts_vw.val << " +/- " << l1SampleCounts_vw.err << "\n";
            std::cout << std::left << std::setw(35) << "2. Source Ratio:" << sourceRatio_vw.val << " +/- " << sourceRatio_vw.err << "\n";
            std::cout << std::left << std::setw(35) << "   ==> Final Source (Denominator):" << final_source_yield_vw.val << " +/- " << final_source_yield_vw.err << "\n\n";
            std::cout << std::left << std::setw(35) << "3. Raw Fragment Counts:" << rawFragCounts_vw.val << " +/- " << rawFragCounts_vw.err << "\n";
            std::cout << std::left << std::setw(35) << "4. Contamination Ratio:" << contamRatio_vw.val << " +/- " << contamRatio_vw.err << "\n";
            std::cout << std::left << std::setw(35) << "   ==> Final Contamination:" << final_contam_yield_vw.val << " +/- " << final_contam_yield_vw.err << "\n\n";
        }

        ValueWithError pure_total_frag_bin(0,0);
        if (config.isAdjacentZ) {
            std::map<std::string, ValueWithError> frag_fracs_vw, norm_fracs_vw;
            ValueWithError frag_frac_sum(0,0), norm_frac_sum(0,0);

            auto get_fractions = [&](TFile* file, const std::vector<std::string>& isotopes, std::map<std::string, ValueWithError>& vw_map, ValueWithError& sum_vw) {
                for (const auto& iso : isotopes) {
                    auto h = getHistFromFile<TH1F>(file, TString::Format("h_best_%s_frac_%s", iso.c_str(), det.c_str()));
                    if(h) { vw_map[iso] = {h->GetBinContent(i_bin), h->GetBinError(i_bin)}; sum_vw = sum_vw + vw_map[iso]; delete h; }
                }
                const auto& last_iso = config.fragmentIsotopes.back();
                vw_map[last_iso] = ValueWithError(1.0, 0.0) - sum_vw;
            };
            get_fractions(fragmentFitFile.get(), config.primaryFitIsotopes, frag_fracs_vw, frag_frac_sum);
            get_fractions(normalFitFile.get(), config.primaryFitIsotopes, norm_fracs_vw, norm_frac_sum);
            
            if (is_debug_target) std::cout << "   --- Isotope Breakdown ---\n";
            for (const auto& iso : config.fragmentIsotopes) {
                ValueWithError raw_iso_counts_vw = rawFragCounts_vw * frag_fracs_vw[iso];
                ValueWithError contamination_term_for_iso = final_contam_yield_vw * norm_fracs_vw[iso];
                ValueWithError pure_iso_counts_vw = raw_iso_counts_vw - contamination_term_for_iso;
                ValueWithError ratio_vw = pure_iso_counts_vw / final_source_yield_vw;
                fillHistBin(h_iss_ratios[iso], i_bin, ratio_vw);
                pure_total_frag_bin = pure_total_frag_bin + pure_iso_counts_vw;

                if (is_debug_target) {
                    std::cout << "   For Isotope " << iso << ":\n";
                    std::cout << "     - Raw Counts Term: " << raw_iso_counts_vw.val << "\n";
                    std::cout << "     - Contam Term: " << contamination_term_for_iso.val << "\n";
                    std::cout << "     - ==> Pure Iso Counts: " << pure_iso_counts_vw.val << "\n";
                    std::cout << "     - ==> FINAL ISS RATIO (" << iso << "): " << ratio_vw.val << "\n";
                }
            }
        } else {
            pure_total_frag_bin = rawFragCounts_vw - final_contam_yield_vw;
        }

        ValueWithError total_ratio_vw = pure_total_frag_bin / final_source_yield_vw;
        fillHistBin(h_iss_ratios["Total"], i_bin, total_ratio_vw);
        
        if (is_debug_target) {
            std::cout << "\n   --- Totals ---\n";
            std::cout << std::left << std::setw(35) << "Total Pure Frag Counts:" << pure_total_frag_bin.val << " +/- " << pure_total_frag_bin.err << "\n";
            std::cout << std::left << std::setw(35) << "FINAL ISS RATIO (Total):" << total_ratio_vw.val << " +/- " << total_ratio_vw.err << "\n";
        }
    }
    
    outputFile->cd(); gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0);
    std::string plot_dir = "/eos/user/z/zixuan/Isotope/BkgValid/Plots/";
    if (config.isAdjacentZ) {
        for (const auto& iso : config.fragmentIsotopes) {
            h_iss_ratios[iso]->Write();
            if (h_mc_ratios.count(iso)) h_mc_ratios[iso]->Write();
            createComparisonPlot(h_iss_ratios[iso], h_mc_ratios.count(iso) ? h_mc_ratios[iso] : nullptr,
                TString::Format("%s#rightarrow^{%d}%s (%s)", config.sourceParticle.c_str(), config.fragZandMass.at(iso), config.fragmentParticle.c_str(), chain.c_str()).Data(),
                TString::Format("L2 frag ^{%d}%s / L1%s ", config.fragZandMass.at(iso), config.fragmentParticle.c_str(), config.sourceParticle.c_str()).Data(),
                plot_dir + TString::Format("%s_to_%s_%s", config.sourceParticle.c_str(), iso.c_str(), chain.c_str()).Data());
        }
    }
    h_iss_ratios["Total"]->Write();
    if (h_mc_ratios.count("Total")) h_mc_ratios["Total"]->Write();
    createComparisonPlot(h_iss_ratios["Total"], h_mc_ratios.count("Total") ? h_mc_ratios["Total"] : nullptr,
        TString::Format("%s#rightarrow%s (%s)", config.sourceParticle.c_str(), config.fragmentParticle.c_str(), chain.c_str()).Data(),
        TString::Format("L2 frag %s / L1 %s", config.fragmentParticle.c_str(), config.sourceParticle.c_str()).Data(),
        plot_dir + TString::Format("%s_to_%s_Total_%s", config.sourceParticle.c_str(), config.fragmentParticle.c_str(), chain.c_str()).Data());

    std::cout << "Analysis complete. Results saved to " << outputFile->GetName() << std::endl;
    for(auto const& [key, val] : h_iss_ratios) delete val;
    for(auto const& [key, val] : h_mc_ratios) delete val;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++ SPECIAL Be10 / B11 ANALYSIS MODULE +++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void runSpecialBe10B11Analysis(const std::string& chain) {
    std::cout << "\n===================================================================\n";
    std::cout << "Starting Special Analysis: B11 -> Be10 for chain: " << chain << std::endl;
    std::cout << "===================================================================\n";

    const AnalysisConfig config_B_to_Be = {
        "Boron", "Beryllium", true, {"Be7", "Be9"}, 7, "Be_frag4",
        {{"B10", 0.3}, {"B11", 0.7}}, {"Be7", "Be9", "Be10"},
        {{"Be7", 7}, {"Be9", 9}, {"Be10", 10}},
        {{"B11", "B11_rew_frag4.root"}}
    };
    const std::array<std::string, 3> detectors = {"TOF", "NaF", "AGL"};
    const std::string plot_dir = "/eos/user/z/zixuan/Isotope/BkgValid/Plots/";

    TString ratioFilePath = TString::Format("/eos/user/z/zixuan/Isotope/ChargeTemp/QFit_Boron_to_Beryllium_%s.root", chain.c_str());
    TString countsFilePath = TString::Format("/eos/user/z/zixuan/Isotope/Add/%s.root", config_B_to_Be.fragFileID.c_str());
    
    auto tempCountsFile = openFile(countsFilePath);
    TH1D* refHist = getHistFromFile<TH1D>(tempCountsFile.get(), TString::Format("%s_ISS_BKG_H1_Boron_TOF", chain.c_str()));
    if (!refHist) throw std::runtime_error("Reference histogram for binning not found in special analysis!");
    refHist->Rebin(2);
    const TAxis* xAxis = refHist->GetXaxis();
    const int nBins = xAxis->GetNbins();
    const Double_t* binEdges = xAxis->GetXbins()->GetArray();
    tempCountsFile.reset();

    auto createHist = [&](const char* name, const char* title) {
        return new TH1F(name, title, nBins, binEdges);
    };

    TH1F* h_iss_ratio_be10_b11 = createHist("h_iss_ratio_Be10_over_B11", "");
    {
        auto ratioFile = openFile(ratioFilePath);
        auto countsFile = openFile(countsFilePath);
        std::string fragAbbr = "Be";
        TString fragmentFitFilePath = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d_FragFrom%s.root", fragAbbr.c_str(), chain.c_str(), config_B_to_Be.useMass, config_B_to_Be.sourceParticle.c_str());
        TString normalFitFilePath = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d.root", fragAbbr.c_str(), chain.c_str(), config_B_to_Be.useMass);
        auto fragmentFitFile = openFile(fragmentFitFilePath);
        auto normalFitFile   = openFile(normalFitFilePath);

        for (int i_bin = 1; i_bin <= nBins; ++i_bin) {
            double binCenter = xAxis->GetBinCenter(i_bin);
            if (binCenter > 21) break;
            std::string det = (binCenter < 1.28) ? "TOF" : (binCenter < 3.06) ? "NaF" : "AGL";

            auto l1SampleCountsHist = getHistFromFile<TH1D>(countsFile.get(), TString::Format("%s_ISS_BKG_H1_Boron_%s", chain.c_str(), det.c_str()));
            auto rawFragCountsHist = getHistFromFile<TH1D>(countsFile.get(), TString::Format("%s_ISS_BKG_H3_Boron_%s", chain.c_str(), det.c_str()));
            auto sourceRatioHist = getHistFromFile<TH1D>(ratioFile.get(), TString::Format("h_narrowfrac_Boron_%s", det.c_str()));
            auto contamRatioHist = getHistFromFile<TH1D>(ratioFile.get(), TString::Format("h_narrowfrac_Beryllium_%s", det.c_str()));
            if(l1SampleCountsHist) l1SampleCountsHist->Rebin(2);
            if(rawFragCountsHist) rawFragCountsHist->Rebin(2);

            ValueWithError l1SampleCounts_vw(l1SampleCountsHist ? l1SampleCountsHist->GetBinContent(i_bin) : 0, l1SampleCountsHist ? l1SampleCountsHist->GetBinError(i_bin) : 0);
            ValueWithError rawFragCounts_vw(rawFragCountsHist ? rawFragCountsHist->GetBinContent(i_bin) : 0, rawFragCountsHist ? rawFragCountsHist->GetBinError(i_bin) : 0);
            ValueWithError sourceRatio_vw(sourceRatioHist ? sourceRatioHist->GetBinContent(i_bin) : 0, sourceRatioHist ? sourceRatioHist->GetBinError(i_bin) : 0);
            ValueWithError contamRatio_vw(contamRatioHist ? contamRatioHist->GetBinContent(i_bin) : 0, contamRatioHist ? contamRatioHist->GetBinError(i_bin) : 0);
            
            ValueWithError final_total_boron_yield_vw = l1SampleCounts_vw * sourceRatio_vw;
            ValueWithError final_b11_yield_vw = final_total_boron_yield_vw * 0.7;

            if (final_b11_yield_vw.val <= 0) continue;

            ValueWithError frag_frac_sum(0,0), norm_frac_sum(0,0);
            std::map<std::string, ValueWithError> frag_fracs_vw, norm_fracs_vw;
            auto get_fractions = [&](TFile* file, std::map<std::string, ValueWithError>& vw_map, ValueWithError& sum_vw) {
                for (const auto& iso : config_B_to_Be.primaryFitIsotopes) {
                    auto h = getHistFromFile<TH1F>(file, TString::Format("h_best_%s_frac_%s", iso.c_str(), det.c_str()));
                    if(h) { vw_map[iso] = {h->GetBinContent(i_bin), h->GetBinError(i_bin)}; sum_vw = sum_vw + vw_map[iso]; delete h; }
                }
                vw_map["Be10"] = ValueWithError(1.0, 0.0) - sum_vw;
            };
            get_fractions(fragmentFitFile.get(), frag_fracs_vw, frag_frac_sum);
            get_fractions(normalFitFile.get(), norm_fracs_vw, norm_frac_sum);

            ValueWithError final_contam_yield_vw = l1SampleCounts_vw * contamRatio_vw;
            ValueWithError raw_be10_counts_vw = rawFragCounts_vw * frag_fracs_vw["Be10"];
            ValueWithError contam_term_for_be10 = final_contam_yield_vw * norm_fracs_vw["Be10"];
            ValueWithError pure_be10_counts_vw = raw_be10_counts_vw - contam_term_for_be10;
            
            ValueWithError ratio_vw = pure_be10_counts_vw / final_b11_yield_vw;
            fillHistBin(h_iss_ratio_be10_b11, i_bin, ratio_vw);

            if(l1SampleCountsHist) delete l1SampleCountsHist;
            if(rawFragCountsHist) delete rawFragCountsHist;
            if(sourceRatioHist) delete sourceRatioHist;
            if(contamRatioHist) delete contamRatioHist;
        }
    }

    TH1F* h_mc_ratio_be10_b11 = createHist("h_mc_ratio_Be10_over_B11", "");
    {
        auto mc_file = openFile(TString::Format("/eos/user/z/zixuan/Isotope/Add/%s", config_B_to_Be.mcSourceFiles.at("B11").c_str()));

        std::map<std::string, TH1*> source_hists_per_det, frag_hists_per_det;
        for(const auto& det : detectors) {
            source_hists_per_det[det] = getHistFromFile<TH1>(mc_file.get(), TString::Format("%s_MC_BKG_H1_%s", chain.c_str(), det.c_str()));
            if(source_hists_per_det[det]) source_hists_per_det[det]->Rebin(2);
            
            frag_hists_per_det[det] = getHistFromFile<TH1>(mc_file.get(), TString::Format("%s_MC_BKG_H2_%s_Z4_Mass10", chain.c_str(), det.c_str()));
            if(frag_hists_per_det[det]) frag_hists_per_det[det]->Rebin(2);
        }
        
        std::unique_ptr<TH1F> stitched_source(stitchHistograms(source_hists_per_det, refHist));
        std::unique_ptr<TH1F> stitched_frag(stitchHistograms(frag_hists_per_det, refHist));
        
        h_mc_ratio_be10_b11->Add(stitched_frag.get());
        h_mc_ratio_be10_b11->Divide(stitched_source.get());
    }

    createComparisonPlot(h_iss_ratio_be10_b11, h_mc_ratio_be10_b11,
        TString::Format("^{11}Boron#rightarrow^{10}Beryllium (%s)", chain.c_str()).Data(),
        "L2 frag ^{10}Beryllium / L1^{11}Boron",
        plot_dir + TString::Format("Special_Be10_over_B11_%s", chain.c_str()).Data());

    std::cout << "Special B11->Be10 analysis for chain " << chain << " complete." << std::endl;
    delete h_iss_ratio_be10_b11;
    delete h_mc_ratio_be10_b11;
}


// --- Main entry point ---
void CalFrag() {
    AnalysisConfig config_B_to_Be;
    config_B_to_Be.sourceParticle = "Boron"; config_B_to_Be.fragmentParticle = "Beryllium";
    config_B_to_Be.isAdjacentZ = true;
    config_B_to_Be.primaryFitIsotopes = {"Be7", "Be9"};
    config_B_to_Be.useMass = 7;
    config_B_to_Be.fragFileID = "Be_frag4";
    config_B_to_Be.sourceIsotopeFractions = {{"B10", 0.3}, {"B11", 0.7}};
    config_B_to_Be.fragmentIsotopes = {"Be7", "Be9", "Be10"};
    config_B_to_Be.fragZandMass = {{"Be7", 7}, {"Be9", 9}, {"Be10", 10}};
    config_B_to_Be.mcSourceFiles = {{"B10", "B10_rew_frag4.root"}, {"B11", "B11_rew_frag4.root"}};
    
    AnalysisConfig config_C_to_B;
    config_C_to_B.sourceParticle = "Carbon"; config_C_to_B.fragmentParticle = "Boron";
    config_C_to_B.isAdjacentZ = true;
    config_C_to_B.primaryFitIsotopes = {"B10"};
    config_C_to_B.useMass = 10;
    config_C_to_B.fragFileID = "B_frag5";
    config_C_to_B.sourceIsotopeFractions = {{"C12", 1.0}};
    config_C_to_B.fragmentIsotopes = {"B10", "B11"};
    config_C_to_B.fragZandMass = {{"B10", 10}, {"B11", 11}};
    config_C_to_B.mcSourceFiles = {{"C12", "C12_rew_frag5.root"}};
    
    AnalysisConfig config_N_to_B;
    config_N_to_B.sourceParticle = "Nitrogen"; config_N_to_B.fragmentParticle = "Boron";
    config_N_to_B.isAdjacentZ = false;
    config_N_to_B.fragFileID = "B_frag5";
    config_N_to_B.sourceIsotopeFractions = {{"N14", 0.5}, {"N15", 0.5}};
    config_N_to_B.fragmentIsotopes = {"B10", "B11"};
    config_N_to_B.fragZandMass = {{"B10", 10}, {"B11", 11}};
    config_N_to_B.mcSourceFiles = {{"N14", "N14_rew_frag5.root"}, {"N15", "N15_rew_frag5.root"}};
    
    AnalysisConfig config_O_to_B;
    config_O_to_B.sourceParticle = "Oxygen"; config_O_to_B.fragmentParticle = "Boron";
    config_O_to_B.isAdjacentZ = false;
    config_O_to_B.fragFileID = "B_frag5";
    config_O_to_B.sourceIsotopeFractions = {{"O16", 1.0}};
    config_O_to_B.fragmentIsotopes = {"B10", "B11"};
    config_O_to_B.fragZandMass = {{"B10", 10}, {"B11", 11}};
    config_O_to_B.mcSourceFiles = {{"O16", "O16_rew_frag5.root"}};
    
    AnalysisConfig config_C_to_Be;
    config_C_to_Be.sourceParticle = "Carbon"; config_C_to_Be.fragmentParticle = "Beryllium";
    config_C_to_Be.isAdjacentZ = false;
    config_C_to_Be.fragFileID = "Be_frag4";
    config_C_to_Be.sourceIsotopeFractions = {{"C12", 1.0}};
    config_C_to_Be.fragmentIsotopes = {"Be7", "Be9", "Be10"};
    config_C_to_Be.fragZandMass = {{"Be7", 7}, {"Be9", 9}, {"Be10", 10}};
    config_C_to_Be.mcSourceFiles = {{"C12", "C12_rew_frag4.root"}};
    
    AnalysisConfig config_N_to_Be;
    config_N_to_Be.sourceParticle = "Nitrogen"; config_N_to_Be.fragmentParticle = "Beryllium";
    config_N_to_Be.isAdjacentZ = false;
    config_N_to_Be.fragFileID = "Be_frag4";
    config_N_to_Be.sourceIsotopeFractions = {{"N14", 0.5}, {"N15", 0.5}};
    config_N_to_Be.fragmentIsotopes = {"Be7", "Be9", "Be10"};
    config_N_to_Be.fragZandMass = {{"Be7", 7}, {"Be9", 9}, {"Be10", 10}};
    config_N_to_Be.mcSourceFiles = {{"N14", "N14_rew_frag4.root"}, {"N15", "N15_rew_frag4.root"}};
    
    AnalysisConfig config_O_to_Be;
    config_O_to_Be.sourceParticle = "Oxygen"; config_O_to_Be.fragmentParticle = "Beryllium";
    config_O_to_Be.isAdjacentZ = false;
    config_O_to_Be.fragFileID = "Be_frag4";
    config_O_to_Be.sourceIsotopeFractions = {{"O16", 1.0}};
    config_O_to_Be.fragmentIsotopes = {"Be7", "Be9", "Be10"};
    config_O_to_Be.fragZandMass = {{"Be7", 7}, {"Be9", 9}, {"Be10", 10}};
    config_O_to_Be.mcSourceFiles = {{"O16", "O16_rew_frag4.root"}};
    
    try {
        const std::vector<std::string> chains = {"L1Inner"};
        const std::vector<AnalysisConfig> all_configs = {
            config_B_to_Be, config_C_to_B,
            config_N_to_B, config_O_to_B,
            config_C_to_Be, config_N_to_Be, config_O_to_Be
        };
        for (const auto& chain : chains) {
            for (const auto& config : all_configs) {
                runAnalysis(chain, config);
            }
        }

        // Run the special analysis for each chain
        for (const auto& chain : chains) {
            runSpecialBe10B11Analysis(chain);
        }

    } catch (const std::exception& e) {
        std::cerr << "[FATAL ERROR] An exception occurred: " << e.what() << std::endl;
    }
}