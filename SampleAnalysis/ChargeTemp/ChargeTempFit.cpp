#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <cmath>

#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooArgList.h>
#include <RooMsgService.h>
#include <RooFormulaVar.h>
#include <RooAbsReal.h>

#include "../Tool.h"

using namespace AMS_Iso;
using namespace RooFit;
using namespace std;

// --- Global Configuration ---
const std::string inputFileName = "/eos/user/z/zixuan/Isotope/Add/Be_frag4.root";
const std::string outputDir = "/eos/user/z/zixuan/Isotope/ChargeTemp/"; 
const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};

// --- Detector-specific energy ranges ---
const std::map<std::string, std::pair<double, double>> detector_ek_ranges = {
    {"TOF", {0.4, 1.29}},
    {"NaF", {0.90, 5.10}},
    {"AGL", {2.90, 22.0}}
};

// --- Helper Utilities ---
struct ElementInfo { int Z; std::string type; };
const std::map<std::string, ElementInfo> element_db = {
    {"Beryllium", {4, "Secondary"}}, {"Boron",     {5, "Secondary"}},
    {"Carbon",    {6, "Primary"}},   {"Nitrogen",  {7, "Secondary"}},
    {"Oxygen",    {8, "Primary"}}
};
const std::map<int, std::string> Z_to_name = {
    {4, "Beryllium"}, {5, "Boron"}, {6, "Carbon"}, {7, "Nitrogen"}, {8, "Oxygen"}
};

// --- Result Data Structures ---
struct FitResult {
    double ekpernuc_center = 0.0;
    double ekpernuc_low = 0.0;
    double ekpernuc_up = 0.0;
    double ekpernuc_width = 0.0;
    int fit_status = -1;
    double chi2ndf = 0.0;
    int ndf = 0;
    double N_sig_narrow = 0.0;
    double fragment_yield = 0.0;
    double fragment_yield_err = 0.0;
    double source_yield = 0.0;
    double source_yield_err = 0.0;
    std::map<std::string, double> fit_fractions;
    std::map<std::string, double> fit_fractions_err;
    std::map<std::string, double> narrow_fractions;
    std::map<std::string, double> narrow_fractions_err;
};

// --- Core Fitting Class ---
class ChargeFitProcessor {
public:
    ChargeFitProcessor(
        const std::string& source, const std::string& fragment, const std::string& chain, const std::string& detector,
        int ekpernuc_bin, TH2F* signal_rebinned, const std::map<std::string, TH2F*>& templates_rebinned);

    bool initializeAndProject();
    bool runFit();
    FitResult calculateFinalYield();
    std::unique_ptr<RooPlot> generatePlotAndCalcChi2(FitResult& result);
    const std::vector<std::string>& getTemplateElements() const { return templateElements_; }

private:
    std::string sourceName_, fragmentName_, chainName_, detectorName_;
    int energyBin_;
    TH2F* h_signal_rebinned_;
    const std::map<std::string, TH2F*>& templates_rebinned_;
    std::vector<std::string> templateElements_;
    double fitMin_, fitMax_, narrowMin_, narrowMax_;
    std::unique_ptr<RooRealVar> charge_;
    std::unique_ptr<TH1D> h_signal_extended_;
    std::map<std::string, std::unique_ptr<TH1D>> h_templates_extended_;
    std::vector<std::unique_ptr<RooRealVar>> fractionParams_;
    std::unique_ptr<RooFormulaVar> lastFraction_;
    std::unique_ptr<RooAddPdf> total_pdf_;
    std::unique_ptr<RooFitResult> fitResult_;
    std::unique_ptr<RooDataHist> data_hist_;
    std::map<std::string, std::unique_ptr<RooHistPdf>> template_pdfs_;
    std::map<std::string, std::unique_ptr<RooDataHist>> template_data_hists_;
    void configure();
    std::unique_ptr<TH1D> projectSlice(TH2F* h2d, const char* name);
};

ChargeFitProcessor::ChargeFitProcessor(
    const std::string& source, const std::string& fragment, const std::string& chain, const std::string& detector,
    int ekpernuc_bin, TH2F* signal_rebinned, const std::map<std::string, TH2F*>& templates_rebinned)
    : sourceName_(source), fragmentName_(fragment), chainName_(chain), detectorName_(detector),
      energyBin_(ekpernuc_bin), h_signal_rebinned_(signal_rebinned), templates_rebinned_(templates_rebinned)
{
    configure();
}

void ChargeFitProcessor::configure() {
    int z_source = element_db.at(sourceName_).Z;

    // --- NEW TEMPLATE SELECTION STRATEGY ---
    templateElements_.clear();
    if (sourceName_ == "Oxygen") {
        // Special case for Oxygen: fit Carbon, Nitrogen, and Oxygen
        templateElements_.push_back("Carbon");
        templateElements_.push_back("Nitrogen");
        templateElements_.push_back("Oxygen");
        
        // Adjust fit range to include all three templates
        fitMin_ = 6.2; // Low enough to include Carbon (Z=6)
        fitMax_ = 8.6; // High enough to include Oxygen (Z=8) tail
    } else {
        // Default case: fit source, source-1, and source+1
        for (int z_offset = -1; z_offset <= 1; ++z_offset) {
            int current_z = z_source + z_offset;
            if (Z_to_name.count(current_z)) {
                templateElements_.push_back(Z_to_name.at(current_z));
            }
        }
        // Generic fit range for other sources
        fitMin_ = z_source - 1 - 0.2;
        fitMax_ = z_source + 1 + 0.5;
        if(sourceName_ == "Nitrogen") {
             fitMin_ = 5.8; 
             fitMax_ = 8.0;
        }
    }
    
    // Ensure the elements are sorted by Z for consistent parameter ordering
    std::sort(templateElements_.begin(), templateElements_.end(), [&](const std::string& a, const std::string& b) {
        return element_db.at(a).Z < element_db.at(b).Z;
    });

    // Configure narrow range for yield calculation (always based on source)
    const auto& sourceInfo = element_db.at(sourceName_);
    if (sourceInfo.type == "Primary") {
        narrowMin_ = sourceInfo.Z - 0.4;
        narrowMax_ = sourceInfo.Z + 0.4;
    } else {
        narrowMin_ = sourceInfo.Z - 0.2;
        narrowMax_ = sourceInfo.Z + 0.4;
    }
}

std::unique_ptr<TH1D> ChargeFitProcessor::projectSlice(TH2F* h2d, const char* name) {
    auto slice = std::unique_ptr<TH1D>(h2d->ProjectionX(name, energyBin_, energyBin_));
    slice->SetDirectory(nullptr);
    return slice;
}

bool ChargeFitProcessor::initializeAndProject() {
    auto h_signal_slice = projectSlice(h_signal_rebinned_, Form("h_signal_slice_bin%d", energyBin_));
    if (!h_signal_slice || h_signal_slice->GetEntries() < 20) return false;
    
    std::map<std::string, std::unique_ptr<TH1D>> h_templates_slices;
    for (const auto& el : templateElements_) {
        auto it = templates_rebinned_.find(el);
        if (it == templates_rebinned_.end()) {
            std::cerr << "      -> CRITICAL: Template for " << el << " not found." << std::endl;
            return false;
        }
        h_templates_slices[el] = projectSlice(it->second, Form("h_template_%s_slice_bin%d", el.c_str(), energyBin_));
        if (!h_templates_slices[el] || h_templates_slices[el]->GetEntries() < 5) return false;
    }

    h_signal_extended_ = extendHistogram(h_signal_slice.get(), fitMin_, fitMax_);
    for (const auto& el : templateElements_) {
        h_templates_extended_[el] = extendHistogram(h_templates_slices.at(el).get(), fitMin_, fitMax_);
    }

    charge_ = std::make_unique<RooRealVar>("charge", "Charge", fitMin_, fitMax_);
    data_hist_ = std::make_unique<RooDataHist>("data_hist", "Data", *charge_, h_signal_extended_.get());
    for (const auto& el : templateElements_) {
        template_data_hists_[el] = std::make_unique<RooDataHist>(Form("dhist_%s", el.c_str()), "", *charge_, h_templates_extended_.at(el).get());
        template_pdfs_[el] = std::make_unique<RooHistPdf>(Form("pdf_%s", el.c_str()), "", *charge_, *template_data_hists_[el]);
    }

    std::vector<std::string> free_params_elements;
    std::string constrained_element;
    if (templateElements_.size() > 1) {
        // Prioritize source and source-1 as free parameters
        free_params_elements.push_back(sourceName_);
        
        std::string z_minus_1_name = "";
        if (Z_to_name.count(element_db.at(sourceName_).Z - 1)) {
            z_minus_1_name = Z_to_name.at(element_db.at(sourceName_).Z - 1);
        }
        
        if (!z_minus_1_name.empty() && std::find(templateElements_.begin(), templateElements_.end(), z_minus_1_name) != templateElements_.end()) {
             if (free_params_elements.size() < templateElements_.size() - 1) {
                free_params_elements.push_back(z_minus_1_name);
             }
        }
        
        // Fill remaining free parameter slots
        for (const auto& el : templateElements_) {
            if (free_params_elements.size() >= templateElements_.size() - 1) break;
            if (std::find(free_params_elements.begin(), free_params_elements.end(), el) == free_params_elements.end()) {
                free_params_elements.push_back(el);
            }
        }
        // Find the one remaining element to be constrained
        for (const auto& el : templateElements_) {
            if (std::find(free_params_elements.begin(), free_params_elements.end(), el) == free_params_elements.end()) {
                constrained_element = el;
                break;
            }
        }
    } else if (!templateElements_.empty()) {
        free_params_elements.push_back(templateElements_[0]);
    }

    std::map<std::string, RooAbsReal*> frac_map;
    double initial_other_frac = (free_params_elements.size() > 1) ? 0.02 / (free_params_elements.size() - 1) : 0.0;
    for (const auto& el : free_params_elements) {
        double initial_val = (el == sourceName_) ? 0.98 : initial_other_frac;
        auto fracVar = std::make_unique<RooRealVar>(Form("frac_%s", el.c_str()), "", initial_val, 0.0, 1.0);
        frac_map[el] = fracVar.get();
        fractionParams_.push_back(std::move(fracVar));
    }
    if (!constrained_element.empty()) {
        std::string formula = "1.0";
        RooArgList formulaArgs;
        for (const auto& param : fractionParams_) {
            formula += " - @" + std::to_string(formulaArgs.getSize());
            formulaArgs.add(*param);
        }
        lastFraction_ = std::make_unique<RooFormulaVar>(Form("frac_%s", constrained_element.c_str()), "", formula.c_str(), formulaArgs);
        frac_map[constrained_element] = lastFraction_.get();
    }

    RooArgList pdfList, fracList;
    for (const auto& el : templateElements_) {
        pdfList.add(*template_pdfs_.at(el));
        fracList.add(*frac_map.at(el));
    }
    total_pdf_ = std::make_unique<RooAddPdf>("total_pdf", "Total PDF", pdfList, fracList, false);
    return true;
}

bool ChargeFitProcessor::runFit() {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

    constexpr double MAX_CHI2NDF = 90.0;
    const double minFitMax = element_db.at(sourceName_).Z + 0.5;
    double currentFitMax = fitMax_;

    while (currentFitMax >= minFitMax) {
        cout << "      -> Attempting fit with range [" << fitMin_ << ", " << currentFitMax << "]..." << endl;

        charge_->setRange("current_fit_range", fitMin_, currentFitMax);

        fitResult_ = std::unique_ptr<RooFitResult>(
            total_pdf_->fitTo(*data_hist_,
                Save(true),
                PrintLevel(-1),
                Range("current_fit_range"),
                Strategy(2),
                Minimizer("Minuit2", "migrad")
            )
        );

        if (fitResult_ && fitResult_->status() == 0) {
            auto tempFrame = std::unique_ptr<RooPlot>(charge_->frame(Range(fitMin_, currentFitMax)));
            data_hist_->plotOn(tempFrame.get(), Name("data_hist"));
            total_pdf_->plotOn(tempFrame.get(), Name("total_pdf")); 
            
            double chi2 = calculateChi2(tempFrame.get(), "data_hist", "total_pdf", fitMin_, currentFitMax);
            
            int nFreeParams = fitResult_->floatParsFinal().getSize();
            int nBinsInRange = 0;
            for(int i = 1; i <= h_signal_extended_->GetNbinsX(); ++i) {
                double binCenter = h_signal_extended_->GetBinCenter(i);
                if (binCenter >= fitMin_ && binCenter <= currentFitMax) {
                    nBinsInRange++;
                }
            }
            int ndf = nBinsInRange - nFreeParams;
            double chi2ndf = (ndf > 0) ? chi2 / ndf : 9999.0;

            if (chi2ndf < MAX_CHI2NDF) {
                cout << "      -> Fit Successful! Chi2/NDF = " << chi2ndf << " in range [" << fitMin_ << ", " << currentFitMax << "]" << endl;
                fitMax_ = currentFitMax;
                return true;
            } else {
                cout << "      -> Fit converged but Chi2/NDF is poor: " << chi2ndf << ". Shrinking range." << endl;
            }
        } else {
            cout << "      -> Fit failed to converge. Shrinking range." << endl;
        }

        currentFitMax -= 0.1;
    }

    cout << "      -> All fit attempts failed after shrinking range to its limit." << endl;
    return false;
}

FitResult ChargeFitProcessor::calculateFinalYield() {
    FitResult res;
    res.fit_status = fitResult_->status();

    const RooArgList& finalPars = fitResult_->floatParsFinal();
    for (const auto& el : templateElements_) {
        RooAbsReal* frac_param = (RooAbsReal*)finalPars.find(Form("frac_%s", el.c_str()));
        if (frac_param) {
            res.fit_fractions[el] = frac_param->getVal();
            RooRealVar* real_var = dynamic_cast<RooRealVar*>(frac_param);
            if (real_var) res.fit_fractions_err[el] = real_var->getError();
        } else {
            if (lastFraction_ && lastFraction_->GetName() == std::string("frac_" + el)) {
                res.fit_fractions[el] = lastFraction_->getVal();
                res.fit_fractions_err[el] = lastFraction_->getPropagatedError(*fitResult_);
            }
        }
    }

    int bin_low = h_signal_extended_->GetXaxis()->FindBin(narrowMin_);
    int bin_high = h_signal_extended_->GetXaxis()->FindBin(narrowMax_);
    res.N_sig_narrow = h_signal_extended_->Integral(bin_low, bin_high);
    double err_N_sig = sqrt(res.N_sig_narrow);

    charge_->setRange("narrow_range", narrowMin_, narrowMax_);
    std::map<std::string, double> shape_factors;
    for (const auto& el : templateElements_) {
        auto integral_obj = std::unique_ptr<RooAbsReal>(template_pdfs_.at(el)->createIntegral(*charge_, NormSet(*charge_), Range("narrow_range")));
        shape_factors[el] = integral_obj->getVal();
    }
    double total_model_narrow_integral = 0.0;
    for (const auto& el : templateElements_) {
        total_model_narrow_integral += res.fit_fractions.at(el) * shape_factors.at(el);
    }

    bool fragment_in_fit = std::find(templateElements_.begin(), templateElements_.end(), fragmentName_) != templateElements_.end();
    if (fragment_in_fit && total_model_narrow_integral > 1e-9) {
        double P_frag = (res.fit_fractions.at(fragmentName_) * shape_factors.at(fragmentName_)) / total_model_narrow_integral;
        double F_Y = res.fit_fractions.at(fragmentName_);
        double err_F_Y = res.fit_fractions_err.at(fragmentName_);
        double err_P_frag_approx = (F_Y > 1e-9) ? (err_F_Y / F_Y) * P_frag : 0.0;
        
        res.fragment_yield = res.N_sig_narrow * P_frag;
        if (res.fragment_yield > 0) {
            res.fragment_yield_err = sqrt(pow(err_N_sig * P_frag, 2) + pow(res.N_sig_narrow * err_P_frag_approx, 2));
        }
    } else {
        res.fragment_yield = 0.0;
        res.fragment_yield_err = 0.0;
    }

    if (total_model_narrow_integral > 1e-9) {
        double P_source = (res.fit_fractions.at(sourceName_) * shape_factors.at(sourceName_)) / total_model_narrow_integral;
        double F_S = res.fit_fractions.at(sourceName_);
        double err_F_S = res.fit_fractions_err.at(sourceName_);
        double err_P_source_approx = (F_S > 1e-9) ? (err_F_S / F_S) * P_source : 0.0;

        res.source_yield = res.N_sig_narrow * P_source;
        if (res.source_yield > 0) {
            res.source_yield_err = sqrt(pow(err_N_sig * P_source, 2) + pow(res.N_sig_narrow * err_P_source_approx, 2));
        }
    }

    if(total_model_narrow_integral > 1e-9){
        for(const auto& el : templateElements_){
            res.narrow_fractions[el] = (res.fit_fractions.at(el) * shape_factors.at(el)) / total_model_narrow_integral;
            res.narrow_fractions_err[el] = 0.0;
        }
    }
    
    charge_->setRange("full_range", charge_->getMin(), charge_->getMax());
    return res;
}

std::unique_ptr<RooPlot> ChargeFitProcessor::generatePlotAndCalcChi2(FitResult& result) {
    auto frame = std::unique_ptr<RooPlot>(charge_->frame(Range(fitMin_, fitMax_), Title(Form("Final Fit in [%.2f, %.2f]", fitMin_, fitMax_))));
    
    data_hist_->plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
    
    total_pdf_->plotOn(frame.get(), Name("total_pdf"), LineColor(kRed), LineWidth(2));
    
    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kOrange, kCyan, kYellow + 2};
    int color_idx = 0;
    for (const auto& el : templateElements_) {
        total_pdf_->plotOn(frame.get(), Components(*template_pdfs_.at(el)), Name(Form("comp_%s", el.c_str())), LineColor(colors[color_idx % colors.size()]), LineStyle(kDashed));
        color_idx++;
    }
    
    int nFreeParams = fitResult_->floatParsFinal().getSize();
    
    int nBinsInRange = 0;
    for(int i = 1; i <= h_signal_extended_->GetNbinsX(); ++i) {
        double binCenter = h_signal_extended_->GetBinCenter(i);
        if (binCenter >= fitMin_ && binCenter <= fitMax_) {
            nBinsInRange++;
        }
    }
    result.ndf = nBinsInRange - nFreeParams;

    double chi2 = calculateChi2(frame.get(), "data_hist", "total_pdf", fitMin_, fitMax_);
    result.chi2ndf = (result.ndf > 0) ? chi2 / result.ndf : 0.0;

    return frame;
}

void runFragmentationAnalysis(const std::string& source, const std::string& fragment, const std::string& chain) {
    if (element_db.find(source) == element_db.end() || element_db.find(fragment) == element_db.end()) return;
    std::cout << "Starting analysis for: " << source << " -> " << fragment << " (" << chain << ")" << std::endl;
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(inputFileName.c_str()));
    if (!inputFile || inputFile->IsZombie()) return;
    
    std::string pdf_filename = outputDir + "QFit_" + source + "_to_" + fragment + "_" + chain + ".pdf";
    std::string root_filename = outputDir + "QFit_" + source + "_to_" + fragment + "_" + chain + ".root";
    auto c_pdf = std::make_unique<TCanvas>("c_pdf", "PDF Canvas", 800, 600);
    c_pdf->Print((pdf_filename + "[").c_str());
    auto outputFile = std::make_unique<TFile>(root_filename.c_str(), "RECREATE");

    for (const auto& detector : detectors) {
        std::cout << "\n--- Processing Detector: " << detector << " ---" << std::endl;
        
        std::map<std::string, std::unique_ptr<TH2F>> templates_rebinned;
        
        // --- NEW TEMPLATE LOADING STRATEGY ---
        std::vector<std::string> required_elements;
        if (source == "Oxygen") {
            // Special case for Oxygen: load C, N, O templates
            required_elements.push_back("Carbon");
            required_elements.push_back("Nitrogen");
            required_elements.push_back("Oxygen");
        } else {
            // Default case: load source, source-1, and source+1
            int z_source = element_db.at(source).Z;
            for (int z_offset = -1; z_offset <= 1; ++z_offset) {
                int current_z = z_source + z_offset;
                if (Z_to_name.count(current_z)) {
                    required_elements.push_back(Z_to_name.at(current_z));
                }
            }
        }

        std::string signalHistName = chain + "_ISS_BKG_H2_" + source + "_L1QSignal_" + detector;
        TH2F* h_signal_raw = (TH2F*)inputFile->Get(signalHistName.c_str());
        if (!h_signal_raw) continue;
        auto h_signal_rebinned = std::unique_ptr<TH2F>((TH2F*)h_signal_raw->Clone(Form("%s_rebinned", signalHistName.c_str())));
        h_signal_rebinned->RebinY(2);

        bool all_templates_found = true;
        for (const auto& el : required_elements) {
            const auto& el_info = element_db.at(el);
            std::string template_type;
            if (source == "Nitrogen" || source == "Oxygen") {
                template_type = "L2QTemplate";
            } else {
                template_type = (el_info.type == "Primary") ? "L1QTemplate" : "L2QTemplate";
            }
            std::string templateHistName = chain + "_ISS_BKG_H2_" + el + "_" + template_type + "_" + detector;
            TH2F* h_template_raw = (TH2F*)inputFile->Get(templateHistName.c_str());
            if (!h_template_raw) { 
                std::cerr << "CRITICAL: Could not find required template: " << templateHistName << std::endl;
                all_templates_found = false; 
                break; 
            }
            templates_rebinned[el] = std::unique_ptr<TH2F>((TH2F*)h_template_raw->Clone(Form("%s_rebinned", templateHistName.c_str())));
            templates_rebinned[el]->RebinY(2);
        }
        if (!all_templates_found) continue;

        int n_bins_y = h_signal_rebinned->GetNbinsY();
        const TAxis* y_axis = h_signal_rebinned->GetYaxis();
        auto h_yield = std::make_unique<TH1D>(Form("h_yield_%s_%s", fragment.c_str(), detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
        auto h_source_yield = std::make_unique<TH1D>(Form("h_yield_%s_%s", source.c_str(), detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
        auto h_chi2ndf = std::make_unique<TH1D>(Form("h_chi2ndf_%s", detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
        std::map<std::string, std::unique_ptr<TH1D>> h_fitfracs, h_narrowfracs;
        
        // Create histograms for all potentially required elements to avoid runtime errors
        for (const auto& pair : element_db) {
            const std::string& el = pair.first;
            h_fitfracs[el] = std::make_unique<TH1D>(Form("h_fitfrac_%s_%s", el.c_str(), detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
            h_narrowfracs[el] = std::make_unique<TH1D>(Form("h_narrowfrac_%s_%s", el.c_str(), detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
        }


        for (int y_bin = 1; y_bin <= n_bins_y; ++y_bin) {
            double ek_center = y_axis->GetBinCenter(y_bin);
            cout<<"  [bin " << y_bin << "] Eklow = " << y_axis->GetBinLowEdge(y_bin) << " GeV/n" << endl;
            const auto& range = detector_ek_ranges.at(detector);
            if (ek_center < range.first || ek_center > range.second) continue;

            std::cout << "    Processing energy bin " << y_bin << " (Ek=" << ek_center << " GeV/n)" << std::endl;
            std::map<std::string, TH2F*> templates_raw_ptr;
            for(auto const& [key, val] : templates_rebinned) templates_raw_ptr[key] = val.get();
            
            ChargeFitProcessor processor(source, fragment, chain, detector, y_bin, h_signal_rebinned.get(), templates_raw_ptr);
            if (!processor.initializeAndProject()) continue;
            if (!processor.runFit()) continue;
            
            FitResult result = processor.calculateFinalYield();
            result.ekpernuc_center = ek_center;
            result.ekpernuc_low = y_axis->GetBinLowEdge(y_bin); 
            result.ekpernuc_up = y_axis->GetBinUpEdge(y_bin);
            result.ekpernuc_width = y_axis->GetBinWidth(y_bin);
            auto frame = processor.generatePlotAndCalcChi2(result);

            h_yield->SetBinContent(y_bin, result.fragment_yield);
            h_yield->SetBinError(y_bin, result.fragment_yield_err);
            h_source_yield->SetBinContent(y_bin, result.source_yield);
            h_source_yield->SetBinError(y_bin, result.source_yield_err);
            h_chi2ndf->SetBinContent(y_bin, result.chi2ndf);
            for (const auto& el : processor.getTemplateElements()) {
                if (h_fitfracs.count(el)) {
                    h_fitfracs.at(el)->SetBinContent(y_bin, result.fit_fractions[el]);
                    h_fitfracs.at(el)->SetBinError(y_bin, result.fit_fractions_err[el]);
                    h_narrowfracs.at(el)->SetBinContent(y_bin, result.narrow_fractions[el]);
                    h_narrowfracs.at(el)->SetBinError(y_bin, result.narrow_fractions_err[el]);
                }
            }

            c_pdf->Clear(); c_pdf->Divide(1, 2);
            TPad* pad1 = (TPad*)c_pdf->cd(1);
            pad1->SetPad(0, 0.3, 1, 1); pad1->SetLogy(); pad1->SetBottomMargin(0.02);
            
            frame->SetTitle(Form("%s #rightarrow %s (%s, E_{k}=%.2f-%.2f GeV/n)", source.c_str(), fragment.c_str(), detector.c_str(), result.ekpernuc_low, result.ekpernuc_up));
            frame->GetYaxis()->SetTitle("Events"); frame->GetXaxis()->SetLabelSize(0);
            frame->SetMinimum(0.5);
            frame->Draw();

            auto legend = std::make_unique<TLegend>(0.7, 0.55, 0.88, 0.88);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.03);
            legend->AddEntry("data_hist", "Data", "pe");
            legend->AddEntry("total_pdf", "Total Fit", "l");
            std::vector<int> colors = {kBlue + 3, kGreen + 3, kMagenta, kOrange + 3, kCyan, kYellow + 2};
            int color_idx = 0;
            for (const auto& el : processor.getTemplateElements()) {
                legend->AddEntry(Form("comp_%s", el.c_str()), el.c_str(), "l");
                color_idx++;
            }
            legend->Draw();

            auto info = std::make_unique<TPaveText>(0.15, 0.35, 0.6, 0.88, "NDC");
            info->SetFillStyle(0); info->SetBorderSize(0); info->SetTextAlign(12);
            info->AddText(Form("#chi^{2}/NDF = %.2f", result.chi2ndf));
            info->AddText("Fit Fractions:");
            for(const auto& el : processor.getTemplateElements()) info->AddText(Form("  F_{%s} = %.3f", el.substr(0,2).c_str(), result.fit_fractions.at(el)));
            info->AddText("Narrow Range Fractions:");
            for(const auto& el : processor.getTemplateElements()) info->AddText(Form("  P_{%s} = %.3f", el.substr(0,2).c_str(), result.narrow_fractions.at(el)));
            info->AddText(Form("N_{sig} (narrow) = %.0f", result.N_sig_narrow));
            info->AddText(Form("%s Yield = %.1f #pm %.1f", fragment.c_str(), result.fragment_yield, result.fragment_yield_err));
            info->AddText(Form("%s Yield = %.1f #pm %.1f", source.c_str(), result.source_yield, result.source_yield_err));
            info->Draw();
            
            TPad* pad2 = (TPad*)c_pdf->cd(2);
            pad2->SetPad(0, 0, 1, 0.3); pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.3); pad2->SetGridy();
            auto pullGraph = std::make_unique<TGraphErrors>();
            calculatePull(frame.get(), pullGraph.get(), frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
            setupPullPlot(pullGraph.get(), frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
            pullGraph->Draw("AP");
            TLine zeroLine(frame->GetXaxis()->GetXmin(), 0, frame->GetXaxis()->GetXmax(), 0);
            zeroLine.SetLineStyle(2); zeroLine.SetLineColor(kGray + 1);
            zeroLine.Draw("SAME");
            c_pdf->Print(pdf_filename.c_str());
        }
        outputFile->cd();
        h_yield->Write();
        h_source_yield->Write();
        h_chi2ndf->Write();
        // Write only the histograms for the elements that were actually used in this detector's fit
        for (const auto& el : required_elements) {
            if (h_fitfracs.count(el)) {
                h_fitfracs.at(el)->Write();
                h_narrowfracs.at(el)->Write();
            }
        }
    }
    c_pdf->Print((pdf_filename + "]").c_str());
    outputFile->Close();
}

void ChargeTempFit() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    //runFragmentationAnalysis("Boron", "Beryllium", "L1Inner");
    //runFragmentationAnalysis("Carbon", "Beryllium", "L1Inner");
    //runFragmentationAnalysis("Nitrogen", "Beryllium", "L1Inner");
    runFragmentationAnalysis("Oxygen", "Boron", "L1Inner");
    //runFragmentationAnalysis("Carbon", "Boron", "L1Inner");
    //runFragmentationAnalysis("Nitrogen", "Boron", "L1Inner");
    //runFragmentationAnalysis("Oxygen", "Boron", "L1Inner");
}