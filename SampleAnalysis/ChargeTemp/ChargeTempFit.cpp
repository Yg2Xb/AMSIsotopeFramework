#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <set>
#include <iomanip>

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
    bool should_fit_fragment_;
    
    std::string sourceName_, fragmentName_, chainName_, detectorName_;
    int energyBin_;
    TH2F* h_signal_rebinned_;
    const std::map<std::string, TH2F*>& templates_rebinned_;
    std::vector<std::string> templateElements_;
    double fitMin_, fitMax_, narrowMin_, narrowMax_; // UNIFIED NARROW RANGE
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
    int z_fragment = element_db.at(fragmentName_).Z;

    should_fit_fragment_ = (z_fragment >= z_source - 1);

    templateElements_.clear();
    std::set<std::string> elements_to_fit;

    elements_to_fit.insert(sourceName_);
    if (Z_to_name.count(z_source - 1)) {
        elements_to_fit.insert(Z_to_name.at(z_source - 1));
    }
    if (Z_to_name.count(z_source + 1)) {
        elements_to_fit.insert(Z_to_name.at(z_source + 1));
    }
    
    if (should_fit_fragment_) {
        elements_to_fit.insert(fragmentName_);
    }
    
    templateElements_.assign(elements_to_fit.begin(), elements_to_fit.end());
    std::sort(templateElements_.begin(), templateElements_.end(), [&](const std::string& a, const std::string& b) {
        return element_db.at(a).Z < element_db.at(b).Z;
    });

    if (!templateElements_.empty()) {
        int min_z = element_db.at(templateElements_.front()).Z;
        int max_z = element_db.at(templateElements_.back()).Z;
        fitMin_ = min_z - 0.4;
        fitMax_ = max_z + 0.5;
    } else {
        fitMin_ = z_source - 1.4;
        fitMax_ = z_source + 1.5;
    }

    const auto& sourceInfo = element_db.at(sourceName_);
    narrowMin_ = sourceInfo.Z - 0.2;
    narrowMax_ = sourceInfo.Z + 0.4;
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
        auto hist_clone = (TH1D*)h_templates_extended_.at(el)->Clone(Form("%s_clone_for_rdh", h_templates_extended_.at(el)->GetName()));
        template_data_hists_[el] = std::make_unique<RooDataHist>(Form("dhist_%s", el.c_str()), "", *charge_, hist_clone);
        template_pdfs_[el] = std::make_unique<RooHistPdf>(Form("pdf_%s", el.c_str()), "", *charge_, *template_data_hists_[el]);
    }

    std::vector<std::string> free_params_elements;
    std::string constrained_element;
    if (templateElements_.size() > 1) {
        free_params_elements.push_back(sourceName_);
        
        bool fragment_is_in_fit = std::find(templateElements_.begin(), templateElements_.end(), fragmentName_) != templateElements_.end();

        if (fragment_is_in_fit && fragmentName_ != sourceName_) {
            free_params_elements.push_back(fragmentName_);
        } else {
            std::string z_minus_1_name = Z_to_name.count(element_db.at(sourceName_).Z - 1) ? Z_to_name.at(element_db.at(sourceName_).Z - 1) : "";
            if (!z_minus_1_name.empty() && std::find(templateElements_.begin(), templateElements_.end(), z_minus_1_name) != templateElements_.end()) {
                 if (std::find(free_params_elements.begin(), free_params_elements.end(), z_minus_1_name) == free_params_elements.end()) {
                    free_params_elements.push_back(z_minus_1_name);
                 }
            }
        }
        
        for (const auto& el : templateElements_) {
            if (free_params_elements.size() >= templateElements_.size() - 1) break;
            if (std::find(free_params_elements.begin(), free_params_elements.end(), el) == free_params_elements.end()) {
                free_params_elements.push_back(el);
            }
        }

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
    for (const auto& el : free_params_elements) {
        double initial_val = 0.1;
        if (el == sourceName_) initial_val = 0.85;
        if (should_fit_fragment_ && el == fragmentName_) initial_val = 0.1;

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
            total_pdf_->fitTo(*data_hist_, Save(true), PrintLevel(-1), Range("current_fit_range"), Strategy(2), Minimizer("Minuit2", "migrad"))
        );

        if (fitResult_ && fitResult_->status() == 0) {
            auto tempFrame = std::unique_ptr<RooPlot>(charge_->frame(Range(fitMin_, currentFitMax)));
            data_hist_->plotOn(tempFrame.get(), Name("data_hist"));
            total_pdf_->plotOn(tempFrame.get(), Name("total_pdf")); 
            double chi2 = calculateChi2(tempFrame.get(), "data_hist", "total_pdf", fitMin_, currentFitMax);
            int nFreeParams = fitResult_->floatParsFinal().getSize();
            int nBinsInRange = 0;
            for(int i = 1; i <= h_signal_extended_->GetNbinsX(); ++i) {
                if (h_signal_extended_->GetBinCenter(i) >= fitMin_ && h_signal_extended_->GetBinCenter(i) <= currentFitMax) nBinsInRange++;
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
    const bool debug_output = true;
    FitResult res;

    if (!fitResult_ || fitResult_->status() != 0) {
        res.fit_status = fitResult_ ? fitResult_->status() : -1;
        return res;
    }
    res.fit_status = fitResult_->status();

    if (!should_fit_fragment_) {
        res.fragment_yield = 0.0;
        res.fragment_yield_err = 0.0;
    }

    if (debug_output) {
        cout << "\n\n<<<<<<<<<< DEBUGGING YIELD CALCULATION >>>>>>>>>>" << endl;
        cout << std::fixed << std::setprecision(5);
        cout << "Source: " << sourceName_ << ", Fragment: " << fragmentName_ << ", Energy Bin: " << energyBin_ << endl;
        if (!should_fit_fragment_) {
            cout << "NOTE: Fragment was not included in the fit model. Its yield is set to 0." << endl;
        }
        cout << "Unified Narrow Range: [" << narrowMin_ << ", " << narrowMax_ << "]" << endl;
    }

    const double N_total_fit = data_hist_->sumEntries();
    if (debug_output) cout << "N_total_fit (in fit range): " << N_total_fit << endl;

    int bin_low = h_signal_extended_->GetXaxis()->FindBin(narrowMin_);
    int bin_high = h_signal_extended_->GetXaxis()->FindBin(narrowMax_);
    res.N_sig_narrow = h_signal_extended_->Integral(bin_low, bin_high);
    double err_N_sig = (res.N_sig_narrow > 0) ? sqrt(res.N_sig_narrow) : 0.0;
    if (debug_output) cout << "N_sig_narrow (in unified range): " << res.N_sig_narrow << " +/- " << err_N_sig << endl;

    charge_->setRange("narrow_range_integral", narrowMin_, narrowMax_);
    auto integral_R_total_obj = std::unique_ptr<RooAbsReal>(total_pdf_->createIntegral(*charge_, NormSet(*charge_), Range("narrow_range_integral")));
    double R_total = integral_R_total_obj->getVal();
    double err_R_total = (N_total_fit > 0) ? sqrt(std::max(0.0, R_total * (1.0 - R_total) / N_total_fit)) : 0.0;
    if (debug_output) cout << "R_total (model integral in narrow range): " << R_total << " +/- " << err_R_total << endl;

    if (debug_output) cout << "\n--- Calculating values for each component ---" << endl;
    for (const auto& el : templateElements_) {
        if (debug_output) cout << "  Element: " << el << endl;
        
        RooAbsReal* frac_param = (RooAbsReal*)fitResult_->floatParsFinal().find(Form("frac_%s", el.c_str()));
        if (!frac_param) { if (lastFraction_ && lastFraction_->GetName() == std::string("frac_" + el)) frac_param = lastFraction_.get(); }
        if (!frac_param) continue;

        double F_A = frac_param->getVal();
        res.fit_fractions[el] = F_A;
        double err_F_A = (N_total_fit > 0) ? sqrt(std::max(0.0, F_A * (1.0 - F_A) / N_total_fit)) : 0.0;
        res.fit_fractions_err[el] = err_F_A;
        if (debug_output) cout << "    F_A (fit fraction): " << F_A << " +/- " << err_F_A << endl;

        auto integral_R_A_obj = std::unique_ptr<RooAbsReal>(template_pdfs_.at(el)->createIntegral(*charge_, NormSet(*charge_), Range("narrow_range_integral")));
        double R_A = integral_R_A_obj->getVal();

        const auto& h_template = h_templates_extended_.at(el);
        double N_template_total = h_template->Integral();
        double k_template_narrow = h_template->Integral(bin_low, bin_high);
        double R_A_for_err = (N_template_total > 0) ? k_template_narrow / N_template_total : 0.0;
        double err_R_A = (N_template_total > 0) ? sqrt(std::max(0.0, R_A_for_err * (1.0 - R_A_for_err) / N_template_total)) : 0.0;
        
        if (debug_output) {
            cout << "    R_A (shape factor, from PDF integral): " << R_A << " +/- " << err_R_A << endl;
        }

        if (R_total <= 1e-9) {
            res.narrow_fractions[el] = 0.0;
            res.narrow_fractions_err[el] = 0.0;
            if (debug_output) cout << "    R_total is zero, skipping P_A and Yield calculation." << endl;
            continue;
        }

        double P_A = F_A * R_A / R_total;
        res.narrow_fractions[el] = P_A;
        double rel_err_sq_P_A = 0.0;
        if (F_A > 0) rel_err_sq_P_A += pow(err_F_A / F_A, 2);
        if (R_A > 0) rel_err_sq_P_A += pow(err_R_A / R_A, 2);
        if (R_total > 0) rel_err_sq_P_A += pow(err_R_total / R_total, 2);
        double err_P_A = P_A * sqrt(rel_err_sq_P_A);
        res.narrow_fractions_err[el] = err_P_A;
        if (debug_output) cout << "    P_A (narrow fraction): " << P_A << " +/- " << err_P_A << endl;

        double Yield_A = res.N_sig_narrow * P_A;
        double rel_err_sq_Yield_A = 0.0;
        if (res.N_sig_narrow > 0) rel_err_sq_Yield_A += pow(err_N_sig / res.N_sig_narrow, 2);
        if (P_A > 0) rel_err_sq_Yield_A += pow(err_P_A / P_A, 2);
        double err_Yield_A = Yield_A * sqrt(rel_err_sq_Yield_A);
        if (debug_output) cout << "    Yield_A: " << Yield_A << " +/- " << err_Yield_A << endl;

        if (el == fragmentName_) {
            res.fragment_yield = Yield_A;
            res.fragment_yield_err = err_Yield_A;
        }
        if (el == sourceName_) {
            res.source_yield = Yield_A;
            res.source_yield_err = err_Yield_A;
        }
    }

    if (debug_output) cout << "\n<<<<<<<<<< END OF DEBUGGING >>>>>>>>>>>>\n" << endl;
    charge_->setRange("full_range", charge_->getMin(), charge_->getMax());
    return res;
}


std::unique_ptr<RooPlot> ChargeFitProcessor::generatePlotAndCalcChi2(FitResult& result) {
    if (!fitResult_) {
        // Create a dummy plot for visualization even if fit was skipped
        auto frame = std::unique_ptr<RooPlot>(charge_->frame(Title(Form("Fit Skipped for %s -> %s", sourceName_.c_str(), fragmentName_.c_str()))));
        data_hist_->plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
        return frame;
    }

    auto frame = std::unique_ptr<RooPlot>(charge_->frame(Range(fitMin_, fitMax_), Title(Form("Final Fit in [%.2f, %.2f]", fitMin_, fitMax_))));
    data_hist_->plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
    total_pdf_->plotOn(frame.get(), Name("total_pdf"), LineColor(kRed), LineWidth(2));
    
    std::vector<int> colors = {kBlue, kGreen + 2, kMagenta, kOrange, kCyan, kYellow + 2};
    for (size_t i = 0; i < templateElements_.size(); ++i) {
        total_pdf_->plotOn(frame.get(), Components(*template_pdfs_.at(templateElements_[i])), Name(Form("comp_%s", templateElements_[i].c_str())), LineColor(colors[i % colors.size()]), LineStyle(kDashed));
    }
    
    int nFreeParams = fitResult_->floatParsFinal().getSize();
    int nBinsInRange = 0;
    for(int i = 1; i <= h_signal_extended_->GetNbinsX(); ++i) {
        if (h_signal_extended_->GetBinCenter(i) >= fitMin_ && h_signal_extended_->GetBinCenter(i) <= fitMax_) nBinsInRange++;
    }
    result.ndf = nBinsInRange - nFreeParams;
    double chi2 = calculateChi2(frame.get(), "data_hist", "total_pdf", fitMin_, fitMax_);
    result.chi2ndf = (result.ndf > 0) ? chi2 / result.ndf : 0.0;

    return frame;
}

void runFragmentationAnalysis(const std::string& source, const std::string& fragment, const std::string& chain) {
    cout << "\n======================================================================\n";
    cout << "Starting analysis for: " << source << " -> " << fragment << " (" << chain << ")" << endl;
    cout << "======================================================================\n";
    
    std::string currentInputFileName;
    if (fragment == "Beryllium") currentInputFileName = "/eos/user/z/zixuan/Isotope/Add/Be_frag4.root";
    else if (fragment == "Boron") currentInputFileName = "/eos/user/z/zixuan/Isotope/Add/B_frag5.root";
    else { cerr << "CRITICAL: No input file defined for fragment: " << fragment << endl; return; }
    
    auto inputFile = std::unique_ptr<TFile>(TFile::Open(currentInputFileName.c_str())); 
    if (!inputFile || inputFile->IsZombie()) { cerr << "CRITICAL: Could not open input file: " << currentInputFileName << endl; return; }
    
    std::string pdf_filename = outputDir + "QFit_" + source + "_to_" + fragment + "_" + chain + ".pdf";
    std::string root_filename = outputDir + "QFit_" + source + "_to_" + fragment + "_" + chain + ".root";
    auto c_pdf = std::make_unique<TCanvas>("c_pdf", "PDF Canvas", 800, 600);
    c_pdf->Print((pdf_filename + "[").c_str());
    auto outputFile = std::make_unique<TFile>(root_filename.c_str(), "RECREATE");

    for (const auto& detector : detectors) {
        cout << "\n--- Processing Detector: " << detector << " ---" << endl;
        
        std::set<std::string> elements_to_load_set;
        elements_to_load_set.insert(source);
        elements_to_load_set.insert(fragment);
        int z_source = element_db.at(source).Z;
        if (Z_to_name.count(z_source - 1)) elements_to_load_set.insert(Z_to_name.at(z_source - 1));
        if (Z_to_name.count(z_source + 1)) elements_to_load_set.insert(Z_to_name.at(z_source + 1));
        std::vector<std::string> required_elements(elements_to_load_set.begin(), elements_to_load_set.end());

        std::string signalHistName = chain + "_ISS_BKG_H2_" + source + "_L1QSignal_" + detector;
        TH2F* h_signal_raw = (TH2F*)inputFile->Get(signalHistName.c_str());
        if (!h_signal_raw) { cout << "Signal histogram not found: " << signalHistName << ". Skipping detector." << endl; continue; }
        auto h_signal_rebinned = std::unique_ptr<TH2F>((TH2F*)h_signal_raw->Clone(Form("%s_rebinned", signalHistName.c_str())));
        h_signal_rebinned->RebinY(2);

        std::map<std::string, std::unique_ptr<TH2F>> templates_rebinned;
        bool all_templates_found = true;
        for (const auto& el : required_elements) {
            const auto& el_info = element_db.at(el);
            std::string template_type = (el_info.type == "Primary" || source == "Nitrogen" || source == "Oxygen") ? "L1QTemplate" : "L2QTemplate";
            std::string templateHistName = chain + "_ISS_BKG_H2_" + el + "_" + template_type + "_" + detector;
            TH2F* h_template_raw = (TH2F*)inputFile->Get(templateHistName.c_str());
            if (!h_template_raw) { 
                cerr << "CRITICAL: Could not find required template: " << templateHistName << endl;
                all_templates_found = false; 
                break; 
            }
            templates_rebinned[el] = std::unique_ptr<TH2F>((TH2F*)h_template_raw->Clone(Form("%s_rebinned", templateHistName.c_str())));
            templates_rebinned[el]->RebinY(2);
        }
        if (!all_templates_found) { cout << "Missing templates. Skipping detector." << endl; continue; }

        int n_bins_y = h_signal_rebinned->GetNbinsY();
        const TAxis* y_axis = h_signal_rebinned->GetYaxis();
        auto h_yield = std::make_unique<TH1D>(Form("h_yield_%s_%s", fragment.c_str(), detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
        auto h_source_yield = std::make_unique<TH1D>(Form("h_yield_%s_%s", source.c_str(), detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
        auto h_chi2ndf = std::make_unique<TH1D>(Form("h_chi2ndf_%s", detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
        std::map<std::string, std::unique_ptr<TH1D>> h_fitfracs, h_narrowfracs;
        for (const auto& pair : element_db) {
            h_fitfracs[pair.first] = std::make_unique<TH1D>(Form("h_fitfrac_%s_%s", pair.first.c_str(), detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
            h_narrowfracs[pair.first] = std::make_unique<TH1D>(Form("h_narrowfrac_%s_%s", pair.first.c_str(), detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
        }

        for (int y_bin = 1; y_bin <= n_bins_y; ++y_bin) {
            double ek_center = y_axis->GetBinCenter(y_bin);
            double ek_low = y_axis->GetBinLowEdge(y_bin);
            cout<<"  [bin " << y_bin << "] Eklow = " << ek_low << " GeV/n" << endl;
            const auto& range = detector_ek_ranges.at(detector);
            if (ek_center < range.first || ek_center > range.second) continue;

            cout << "    Processing energy bin " << y_bin << " (Ek=" << ek_center << " GeV/n)" << endl;
            std::map<std::string, TH2F*> templates_raw_ptr;
            for(auto const& [key, val] : templates_rebinned) templates_raw_ptr[key] = val.get();
            
            ChargeFitProcessor processor(source, fragment, chain, detector, y_bin, h_signal_rebinned.get(), templates_raw_ptr);
            if (!processor.initializeAndProject()) { cout << "    Initialization failed. Skipping bin." << endl; continue; }
            if (!processor.runFit()) { cout << "    Fit failed. Skipping bin." << endl; continue; }
            
            FitResult result = processor.calculateFinalYield();
            result.ekpernuc_center = ek_center;
            result.ekpernuc_low = ek_low; 
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
            
            c_pdf->Clear();
            if (result.fit_status != 0 && result.fit_status != -1) { // Handle failed fit visualization
                 auto text = std::make_unique<TPaveText>(0.1, 0.1, 0.9, 0.9);
                text->AddText(Form("%s #rightarrow %s (%s, E_{k}=%.2f-%.2f GeV/n)", source.c_str(), fragment.c_str(), detector.c_str(), result.ekpernuc_low, result.ekpernuc_up));
                text->AddText("FIT FAILED TO CONVERGE");
                text->Draw();
                c_pdf->Print(pdf_filename.c_str());
                continue;
            }

            c_pdf->Divide(1, 2);
            TPad* pad1 = (TPad*)c_pdf->cd(1);
            pad1->SetPad(0, 0.3, 1, 1); pad1->SetLogy(); pad1->SetBottomMargin(0.02);
            
            frame->SetTitle(Form("%s #rightarrow %s (%s, E_{k}=%.2f-%.2f GeV/n)", source.c_str(), fragment.c_str(), detector.c_str(), result.ekpernuc_low, result.ekpernuc_up));
            frame->GetYaxis()->SetTitle("Events"); frame->GetXaxis()->SetLabelSize(0);
            frame->SetMinimum(0.5);
            frame->Draw();

            auto legend = std::make_unique<TLegend>(0.7, 0.55, 0.88, 0.88);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.03);
            legend->AddEntry("data_hist", "Data", "pe");
            if(result.fit_status == 0) legend->AddEntry("total_pdf", "Total Fit", "l");
            std::vector<int> colors = {kBlue + 3, kGreen + 3, kMagenta, kOrange + 3, kCyan, kYellow + 2};
            if(result.fit_status == 0) {
                for (size_t i = 0; i < processor.getTemplateElements().size(); ++i) {
                    legend->AddEntry(Form("comp_%s", processor.getTemplateElements()[i].c_str()), processor.getTemplateElements()[i].c_str(), "l");
                }
            }
            legend->Draw();

            auto info = std::make_unique<TPaveText>(0.15, 0.35, 0.6, 0.88, "NDC");
            info->SetFillStyle(0); info->SetBorderSize(0); info->SetTextAlign(12);
            info->SetTextSize(0.025);
            if(result.fit_status == 0) {
                info->AddText(Form("#chi^{2}/NDF = %.2f", result.chi2ndf));
                info->AddText("Fit Fractions:");
                for(const auto& el : processor.getTemplateElements()) info->AddText(Form("  F_{%s} = %.3f #pm %.3f", el.substr(0,2).c_str(), result.fit_fractions.at(el), result.fit_fractions_err.at(el)));
                info->AddText("Narrow Range Fractions:");
                for(const auto& el : processor.getTemplateElements()) info->AddText(Form("  P_{%s} = %.3f #pm %.3f", el.substr(0,2).c_str(), result.narrow_fractions.at(el), result.narrow_fractions_err.at(el)));
                info->AddText(Form("N_{sig} (narrow) = %.1f", result.N_sig_narrow));
            } else {
                 info->AddText("Fit Skipped (Z_frag < Z_source - 1)");
            }
            info->AddText(Form("%s Yield = %.2f #pm %.2f", fragment.c_str(), result.fragment_yield, result.fragment_yield_err));
            info->AddText(Form("%s Yield = %.2f #pm %.2f", source.c_str(), result.source_yield, result.source_yield_err));
            info->Draw();
            
            TPad* pad2 = (TPad*)c_pdf->cd(2);
            pad2->SetPad(0, 0, 1, 0.3); pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.3); pad2->SetGridy();
            if(result.fit_status == 0) {
                auto pullGraph = std::make_unique<TGraphErrors>();
                calculatePull(frame.get(), pullGraph.get(), frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
                setupPullPlot(pullGraph.get(), frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
                pullGraph->Draw("AP");
                TLine zeroLine(frame->GetXaxis()->GetXmin(), 0, frame->GetXaxis()->GetXmax(), 0);
                zeroLine.SetLineStyle(2); zeroLine.SetLineColor(kGray + 1);
                zeroLine.Draw("SAME");
            }
            c_pdf->Print(pdf_filename.c_str());
        }
        outputFile->cd();
        h_yield->Write();
        h_source_yield->Write();
        h_chi2ndf->Write();
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
    
    runFragmentationAnalysis("Boron", "Beryllium", "L1Inner");
    runFragmentationAnalysis("Carbon", "Beryllium", "L1Inner");
}