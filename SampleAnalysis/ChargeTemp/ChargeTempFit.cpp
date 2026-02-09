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

const std::string outputDir = "/eos/user/z/zixuan/Isotope/ChargeTemp/";
const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
const std::vector<std::string> signal_types = {"L1Sig_Any", "L1Sig_PassLoose", "L1Sig_Pass"};

const std::map<std::string, std::pair<double, double>> detector_ek_ranges = {
  {"TOF", {0.25, 1.5}},
  {"NaF", {0.61, 6.10}},
  {"AGL", {2.70, 23.0}}
};

struct ElementInfo { int Z; std::string name; };
const std::map<int, ElementInfo> element_db_by_z = {
    {2, {2, "Helium"}}, {3, {3, "Lithium"}},
    {4, {4, "Beryllium"}}, {5, {5, "Boron"}}, {6, {6, "Carbon"}},
    {7, {7, "Nitrogen"}}, {8, {8, "Oxygen"}}
};

const std::map<std::string, int> element_db_by_name = {
    {"Helium", 2}, {"Lithium", 3},
    {"Beryllium", 4}, {"Boron", 5}, {"Carbon", 6},
    {"Nitrogen", 7}, {"Oxygen", 8}
};

struct UnifiedFitResult {
    double ekpernuc_center = 0.0;
    double ekpernuc_low = 0.0;
    double ekpernuc_up = 0.0;
    double ekpernuc_width = 0.0;
    
    int fit_status = -1;
    double chi2ndf = 0.0;
    int ndf = 0;

    std::map<std::string, double> fit_fractions;
    std::map<std::string, double> fit_fractions_err;

    std::map<std::string, double> yield_in_window;
    std::map<std::string, double> yield_in_window_err;
};

class UnifiedChargeFitter {
public:
    UnifiedChargeFitter(
        const std::string& chain, const std::string& detector, int ekpernuc_bin, 
        TH2F* signal_rebinned, const std::map<std::string, TH2F*>& templates_rebinned,
        double fitMin, double fitMax, const std::vector<std::string>& templates, const std::string& constrainedElement);

    bool initializeAndProject();
    bool runFit();
    
    std::pair<double, double> calculateYieldInWindow(const std::string& componentName, const std::string& windowCenterElement, bool isNoBkg);

    UnifiedFitResult getFitInfo(); 
    std::unique_ptr<RooPlot> generatePlotAndCalcChi2(UnifiedFitResult& result);
    const std::vector<std::string>& getTemplateElements() const { return templateElements_; }

private:
    std::string chainName_, detectorName_;
    int energyBin_;
    TH2F* h_signal_rebinned_;
    const std::map<std::string, TH2F*>& templates_rebinned_;
    
    std::vector<std::string> templateElements_;
    double fitMin_;
    double fitMax_;
    std::string constrainedElement_;

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

    std::unique_ptr<TH1D> projectSlice(TH2F* h2d, const char* name);
};

UnifiedChargeFitter::UnifiedChargeFitter(
    const std::string& chain, const std::string& detector, int ekpernuc_bin, 
    TH2F* signal_rebinned, const std::map<std::string, TH2F*>& templates_rebinned,
    double fitMin, double fitMax, const std::vector<std::string>& templates, const std::string& constrainedElement)
    : chainName_(chain), detectorName_(detector), energyBin_(ekpernuc_bin),
      h_signal_rebinned_(signal_rebinned), templates_rebinned_(templates_rebinned),
      fitMin_(fitMin), fitMax_(fitMax), templateElements_(templates), constrainedElement_(constrainedElement)
{
    TH1::AddDirectory(kFALSE);
}

std::unique_ptr<TH1D> UnifiedChargeFitter::projectSlice(TH2F* h2d, const char* name) {
    auto slice = std::unique_ptr<TH1D>(h2d->ProjectionX(name, energyBin_, energyBin_));
    return slice;
}

bool UnifiedChargeFitter::initializeAndProject() {
    auto h_signal_slice = projectSlice(h_signal_rebinned_, Form("h_signal_slice_bin%d", energyBin_));
    double N_signal_total = h_signal_slice ? h_signal_slice->GetEntries() : 0;
    
    if (!h_signal_slice || N_signal_total < 50) return false;
    
    std::map<std::string, std::unique_ptr<TH1D>> h_templates_slices;
    double N_templates_total_sum = 0.0;
    
    for (const auto& el : templateElements_) {
        auto it = templates_rebinned_.find(el);
        if (it == templates_rebinned_.end()) return false;
        h_templates_slices[el] = projectSlice(it->second, Form("h_template_%s_slice_bin%d", el.c_str(), energyBin_));
        double N_template = h_templates_slices[el] ? h_templates_slices[el]->GetEntries() : 0;
        if (!h_templates_slices[el] || N_template < 10) return false;
        N_templates_total_sum += N_template;
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

    std::map<std::string, RooAbsReal*> frac_map;
    for (const auto& el : templateElements_) {
        if (el == constrainedElement_) continue;
        double N_template_i = h_templates_slices.at(el)->GetEntries();
        double initial_guess = (N_templates_total_sum > 0) ? N_template_i / N_templates_total_sum : (1.0/templateElements_.size());
        initial_guess = std::max(0.001, std::min(0.999, initial_guess));
        auto fracVar = std::make_unique<RooRealVar>(Form("frac_%s", el.c_str()), "", initial_guess, 0.0, 1.0);
        frac_map[el] = fracVar.get();
        fractionParams_.push_back(std::move(fracVar));
    }

    std::string formula = "1.0";
    RooArgList formulaArgs;
    for (const auto& param : fractionParams_) { 
        formula += " - @" + std::to_string(formulaArgs.getSize());
        formulaArgs.add(*param);
    }
    lastFraction_ = std::make_unique<RooFormulaVar>(Form("frac_%s", constrainedElement_.c_str()), "", formula.c_str(), formulaArgs);
    frac_map[constrainedElement_] = lastFraction_.get();

    RooArgList pdfList, fracList;
    for (const auto& el : templateElements_) {
        pdfList.add(*template_pdfs_.at(el));
        fracList.add(*frac_map.at(el));
    }
    total_pdf_ = std::make_unique<RooAddPdf>("total_pdf", "Total PDF", pdfList, fracList, false); 
    return true;
}

bool UnifiedChargeFitter::runFit() {
    RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);     
    charge_->setRange("fit_range", fitMin_, fitMax_);
    fitResult_ = std::unique_ptr<RooFitResult>(
        total_pdf_->fitTo(*data_hist_, Save(true), PrintLevel(-1), Verbose(false), Range("fit_range"), Strategy(1), Minimizer("Minuit2"))
    );
    RooMsgService::instance().setSilentMode(true);
    return (fitResult_ != nullptr);
}

UnifiedFitResult UnifiedChargeFitter::getFitInfo() {
    UnifiedFitResult res;
    if (!fitResult_) return res;
    res.fit_status = fitResult_->status();

    for (const auto& el : templateElements_) { 
        RooAbsReal* frac_param = (RooAbsReal*)fitResult_->floatParsFinal().find(Form("frac_%s", el.c_str()));
        if (!frac_param) {      
            if (lastFraction_ && lastFraction_->GetName() == std::string("frac_" + el)) frac_param = lastFraction_.get();
        }
        if (!frac_param) continue;
        res.fit_fractions[el] = frac_param->getVal();
        double err = 0.0;
        if (frac_param->IsA()->InheritsFrom(RooRealVar::Class())) err = ((RooRealVar*)frac_param)->getError();
        else {
            double sum_err_sq = 0;
            for(const auto& p : fractionParams_) sum_err_sq += pow(p->getError(), 2);
            err = sqrt(sum_err_sq); 
        }
        res.fit_fractions_err[el] = err;
    }
    return res;
}

std::pair<double, double> UnifiedChargeFitter::calculateYieldInWindow(const std::string& componentName, const std::string& windowCenterElement, bool isNoBkg) {
    if (!fitResult_ || !element_db_by_name.count(windowCenterElement)) return {0.0, 0.0};
    
    if (std::find(templateElements_.begin(), templateElements_.end(), componentName) == templateElements_.end()) return {0.0, 0.0};

    int Z_window = element_db_by_name.at(windowCenterElement);
    
    double low = 0.5, up = 0.5;
    if(isNoBkg) {
        if(Z_window == 5) {low = 0.3; up = 0.4;}
        else if(Z_window == 7) {low = 0.4; up = 0.4;}
    } else {
        if(Z_window == 5 || Z_window == 7) {low = 0.3; up = 0.5;}
    }

    double narrowMin = Z_window - low;
    double narrowMax = Z_window + up;
    
    if (narrowMin < fitMin_ || narrowMax > fitMax_) return {0.0, 0.0};

    std::string range_name = Form("win_%s_for_%s", windowCenterElement.c_str(), componentName.c_str());
    charge_->setRange(range_name.c_str(), narrowMin, narrowMax);

    int bin_low = h_signal_extended_->GetXaxis()->FindBin(narrowMin);
    int bin_high = h_signal_extended_->GetXaxis()->FindBin(narrowMax);
    double N_sig_window = h_signal_extended_->Integral(bin_low, bin_high);
    double err_N_sig_window = (N_sig_window > 0) ? sqrt(N_sig_window) : 0.0;

    auto integral_R_total_obj = std::unique_ptr<RooAbsReal>(total_pdf_->createIntegral(*charge_, NormSet(*charge_), Range(range_name.c_str())));
    double R_total_window = integral_R_total_obj->getVal();

    if (R_total_window < 1e-9) return {0.0, 0.0};

    auto info = getFitInfo();
    double F_A = info.fit_fractions[componentName];
    double err_F_A = info.fit_fractions_err[componentName];

    auto integral_R_A_obj = std::unique_ptr<RooAbsReal>(template_pdfs_.at(componentName)->createIntegral(*charge_, NormSet(*charge_), Range(range_name.c_str())));
    double R_A = integral_R_A_obj->getVal();
    
    const auto& h_template = h_templates_extended_.at(componentName);
    double N_template_total = h_template->Integral();
    double k_template_narrow = h_template->Integral(h_template->GetXaxis()->FindBin(narrowMin), h_template->GetXaxis()->FindBin(narrowMax));
    double R_A_raw = (N_template_total > 0) ? k_template_narrow / N_template_total : 0.0;
    double err_R_A = (N_template_total > 0) ? sqrt(std::max(0.0, R_A_raw * (1.0 - R_A_raw) / N_template_total)) : 0.0;

    double P_A_window = F_A * R_A / R_total_window;
    
    double rel_err_sq = 0.0;
    if (F_A > 0) rel_err_sq += pow(err_F_A / F_A, 2);
    if (R_A > 0) rel_err_sq += pow(err_R_A / R_A, 2);
    double err_P_A_window = P_A_window * sqrt(rel_err_sq);

    double Yield = N_sig_window * P_A_window;
    double err_Yield = Yield * sqrt(pow(err_N_sig_window/N_sig_window, 2) + pow(err_P_A_window/P_A_window, 2));

    return {Yield, err_Yield};
}

std::unique_ptr<RooPlot> UnifiedChargeFitter::generatePlotAndCalcChi2(UnifiedFitResult& result) {
    auto frame = std::unique_ptr<RooPlot>(charge_->frame(Range(fitMin_, fitMax_), Title("Unified Charge Fit")));
    if (!fitResult_) {
        if(data_hist_) data_hist_->plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
        return frame;
    }
    data_hist_->plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
    total_pdf_->plotOn(frame.get(), Name("total_pdf"), LineColor(kRed), LineWidth(2));
    std::vector<int> colors = {28, kBlue, kOrange - 3, kGreen + 2, kAzure + 7, kMagenta, kCyan + 2}; 
    for (size_t i = 0; i < templateElements_.size(); ++i) { 
        total_pdf_->plotOn(frame.get(), Components(*template_pdfs_.at(templateElements_[i])), Name(Form("comp_%s", templateElements_[i].c_str())), LineColor(colors[i % colors.size()]), LineWidth(2));
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

void runUnifiedChargeAnalysis(const std::string& chain, const std::string& mode, const std::string& SourceNuc, bool isNoBkg) {
    std::string prefix = isNoBkg ? "NoBkg" : "withBkg";
    
    std::string signalInputFileName = Form("/eos/user/z/zixuan/Isotope/Add/Be_frag4_%s_Tune_full.root", prefix.c_str(), mode.c_str());
    auto signalFile = std::unique_ptr<TFile>(TFile::Open(signalInputFileName.c_str()));       
    if (!signalFile || signalFile->IsZombie()) { return; }

    std::string templateInputFileName = Form("/eos/user/z/zixuan/Isotope/PureChargeTemp/%s_PureChargeTemplates_%s.root", prefix.c_str(), chain.c_str());
    auto templateFile = std::unique_ptr<TFile>(TFile::Open(templateInputFileName.c_str()));       
    if (!templateFile || templateFile->IsZombie()) { return; }

    std::string pdf_filename = outputDir + prefix + "_" + mode + "QFit_Inner" + SourceNuc  + "_" + chain + ".pdf";
    std::string root_filename = outputDir + prefix + "_" + mode + "QFit_Inner" + SourceNuc + "_" + chain + ".root";
    
    TCanvas* c_pdf = new TCanvas("c_pdf", "PDF Canvas", 800, 600);
    c_pdf->Print((pdf_filename + "[").c_str());
    auto outputFile = std::make_unique<TFile>(root_filename.c_str(), "RECREATE");

    for (const auto& detector : detectors) {
        const std::vector<std::string> required_elements = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"}; 
        std::map<std::string, std::unique_ptr<TH2F>> templates_rebinned;
        bool all_templates_found = true;
        for (const auto& el : required_elements) {
            std::string useMode = (el == "Helium" || el == "Oxygen") ? "Pure" : mode;
            std::string templateHistName = Form("h2d_%sQTemp_%s_%s", useMode.c_str(), el.c_str(), detector.c_str());
            TH2F* h_template_raw = (TH2F*)templateFile->Get(templateHistName.c_str());
            if (!h_template_raw) { all_templates_found = false; break; }
            templates_rebinned[el] = std::unique_ptr<TH2F>((TH2F*)h_template_raw->Clone(Form("%s_rebinned", templateHistName.c_str())));
            templates_rebinned[el]->RebinX(1); 
            if(mode == "Tune") templates_rebinned[el]->Smooth(1,"G"); 
        }
        if (!all_templates_found) { continue; }

        for (const auto& sigType : signal_types) {
            if (SourceNuc != "Beryllium" && sigType == "L1Sig_Any") continue; 

            std::string effectiveSource = SourceNuc;
            if (SourceNuc == "Beryllium" && sigType == "L1Sig_Any") {
                effectiveSource = "Boron";
            }
            std::string signalHistName = chain + "_BKG_H4_" + effectiveSource + "_" + sigType + "_" + detector;
            
            TH2F* h_signal_raw = (TH2F*)signalFile->Get(signalHistName.c_str());
            if (!h_signal_raw) { continue; }
            auto h_signal_rebinned = std::unique_ptr<TH2F>((TH2F*)h_signal_raw->Clone(Form("%s_rebinned", signalHistName.c_str())));
            h_signal_rebinned->RebinX(2);

            int n_bins_y = h_signal_rebinned->GetNbinsY();
            const TAxis* y_axis = h_signal_rebinned->GetYaxis();
            auto h_chi2ndf = std::make_unique<TH1D>(Form("h_chi2ndf_%s_%s_%s", detector.c_str(), sigType.c_str(), chain.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
            
            std::map<std::string, std::unique_ptr<TH1D>> output_hists;
            
            if (SourceNuc == "Beryllium" && (sigType == "L1Sig_Pass" || sigType == "L1Sig_PassLoose")) {
                 for (int z = 4; z <= 8; ++z) {
                     std::string winName = element_db_by_z.at(z).name;
                     std::string rawHistName = Form("h_raw_count_in_%s_window_%s_%s_%s", winName.c_str(), detector.c_str(), sigType.c_str(), chain.c_str());
                     output_hists[rawHistName] = std::make_unique<TH1D>(rawHistName.c_str(), "", n_bins_y, y_axis->GetXbins()->GetArray());
                 }
            }

            if (SourceNuc == "Beryllium") {
                for (int compZ = 4; compZ <= 8; ++compZ) {
                    for (int winZ = 4; winZ <= 8; ++winZ) {
                        std::string compName = element_db_by_z.at(compZ).name;
                        std::string winName = element_db_by_z.at(winZ).name;
                        std::string histName = Form("h_yield_%s_in_%s_window_%s_%s_%s", compName.c_str(), winName.c_str(), detector.c_str(), sigType.c_str(), chain.c_str());
                        output_hists[histName] = std::make_unique<TH1D>(histName.c_str(), "", n_bins_y, y_axis->GetXbins()->GetArray());
                    }
                }
            } else {
                std::string histName = Form("h_yield_%s_in_%s_window_%s_%s_%s", SourceNuc.c_str(), SourceNuc.c_str(), detector.c_str(), sigType.c_str(), chain.c_str());
                output_hists[histName] = std::make_unique<TH1D>(histName.c_str(), "", n_bins_y, y_axis->GetXbins()->GetArray());
            }

            for (int y_bin = 1; y_bin <= n_bins_y; ++y_bin) {
                double ek_center = y_axis->GetBinCenter(y_bin);
                double ek_low = y_axis->GetBinLowEdge(y_bin);
                const auto& range = detector_ek_ranges.at(detector);
                if (ek_center < range.first || ek_center > range.second) continue;

                if (SourceNuc == "Beryllium" && (sigType == "L1Sig_Pass" || sigType == "L1Sig_PassLoose")) {
                    for (int z = 4; z <= 8; ++z) {
                        double winMin = z - 0.5;
                        double winMax = z + 0.5;
                        int binX_min = h_signal_raw->GetXaxis()->FindBin(winMin+0.001);
                        int binX_max = h_signal_raw->GetXaxis()->FindBin(winMax-0.001);
                        double rawCount = h_signal_raw->Integral(binX_min, binX_max, y_bin, y_bin);
                        double rawError = (rawCount > 0) ? sqrt(rawCount) : 0.0;

                        std::string winName = element_db_by_z.at(z).name;
                        std::string rawHistName = Form("h_raw_count_in_%s_window_%s_%s_%s", winName.c_str(), detector.c_str(), sigType.c_str(), chain.c_str());
                        if (output_hists.count(rawHistName)) {
                            output_hists[rawHistName]->SetBinContent(y_bin, rawCount);
                            output_hists[rawHistName]->SetBinError(y_bin, rawError);
                        }
                    }
                }

                std::map<std::string, TH2F*> templates_raw_ptr;
                for(auto const& [key, val] : templates_rebinned) templates_raw_ptr[key] = val.get();
                
                double fitMinVal = 2.6; double fitMaxVal = 8.5;
                std::vector<std::string> current_templates;
                int Z_source = element_db_by_name.at(SourceNuc);

                if (SourceNuc == "Beryllium") {
                    fitMinVal = 2.6; fitMaxVal = 8.5;
                    current_templates = required_elements; 
                } else {
                    fitMinVal = Z_source - 1.4; 
                    fitMaxVal = std::min(8.5, Z_source + 1.4);
                    for (int z = Z_source - 1; z <= Z_source + 1; ++z) {
                        if (element_db_by_z.count(z)) current_templates.push_back(element_db_by_z.at(z).name);
                    }
                }

                if (current_templates.empty()) continue;

                UnifiedChargeFitter fitter(chain, detector, y_bin, h_signal_rebinned.get(), templates_raw_ptr, fitMinVal, fitMaxVal, current_templates, SourceNuc);
                if (!fitter.initializeAndProject()) continue;
                if (!fitter.runFit()) continue;
                
                UnifiedFitResult plot_result = fitter.getFitInfo(); 
                
                if (SourceNuc == "Beryllium") {
                    for (int compZ = 4; compZ <= 8; ++compZ) {
                        std::string compName = element_db_by_z.at(compZ).name;
                        for (int winZ = 4; winZ <= 8; ++winZ) {
                            std::string winName = element_db_by_z.at(winZ).name;
                            
                            std::pair<double, double> res = fitter.calculateYieldInWindow(compName, winName, isNoBkg);
                            std::string histName = Form("h_yield_%s_in_%s_window_%s_%s_%s", compName.c_str(), winName.c_str(), detector.c_str(), sigType.c_str(), chain.c_str());
                            
                            if (output_hists.count(histName)) {
                                output_hists[histName]->SetBinContent(y_bin, res.first);
                                output_hists[histName]->SetBinError(y_bin, res.second);
                            }
                        }
                    }
                } else {
                    std::pair<double, double> res = fitter.calculateYieldInWindow(SourceNuc, SourceNuc, isNoBkg);
                    std::string histName = Form("h_yield_%s_in_%s_window_%s_%s_%s", SourceNuc.c_str(), SourceNuc.c_str(), detector.c_str(), sigType.c_str(), chain.c_str());
                     if (output_hists.count(histName)) {
                        output_hists[histName]->SetBinContent(y_bin, res.first);
                        output_hists[histName]->SetBinError(y_bin, res.second);
                    }
                }
                
                plot_result.ekpernuc_center = ek_center;
                plot_result.ekpernuc_low = ek_low;        
                plot_result.ekpernuc_up = y_axis->GetBinUpEdge(y_bin);
                auto frame = fitter.generatePlotAndCalcChi2(plot_result);
                h_chi2ndf->SetBinContent(y_bin, plot_result.chi2ndf);
                
                c_pdf->Clear();
                c_pdf->Divide(1, 2);
                TPad* pad1 = (TPad*)c_pdf->cd(1);
                pad1->SetPad(0, 0.3, 1, 1); pad1->SetLogy(); pad1->SetBottomMargin(0.02);
                
                frame->SetTitle(Form("%s %s %sFit (E_{k}=%.2f-%.2f GeV/n)", detector.c_str(), sigType.c_str(), mode.c_str(), plot_result.ekpernuc_low, plot_result.ekpernuc_up));
                frame->GetYaxis()->SetTitle("Events"); frame->GetXaxis()->SetLabelSize(0);
                if (sigType == "L1Sig_Any") {
                    frame->SetMinimum(20); 
                } else {
                    frame->SetMinimum(1); 
                }
                auto hproj = std::unique_ptr<TH1D>(h_signal_rebinned->ProjectionX(Form("hproj_%d", y_bin), y_bin, y_bin));
                double ymax = (sigType != "L1Sig_Any") ? hproj->GetMaximum() : hproj->GetBinContent(hproj->FindBin(2.6));
                frame->SetMaximum(5 * ymax); 
                frame->Draw();

                TLegend* legend = new TLegend(0.75, 0.55, 0.93, 0.88); 
                legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.03);
                legend->AddEntry(frame->findObject("data_hist"), "Data", "pe");
                legend->AddEntry(frame->findObject("total_pdf"), "Total Fit", "l");
                for (const auto& el : fitter.getTemplateElements()) {
                    legend->AddEntry(frame->findObject(Form("comp_%s", el.c_str())), el.c_str(), "l");
                }
                legend->Draw(); 

                TPaveText* info = new TPaveText(0.14, 0.48, 0.5, 0.88, "NDC"); 
                info->SetFillStyle(0); info->SetBorderSize(0); info->SetTextAlign(12); info->SetTextSize(0.03);
                info->AddText(Form("#chi^{2}/ndf = %.1f/%d = %.2f", plot_result.chi2ndf * plot_result.ndf, plot_result.ndf, plot_result.chi2ndf));
                info->AddText("Fractions:");
                for (const auto& el : current_templates) {
                    info->AddText(Form(" %s: %.4f#pm%.4f", el.c_str(), plot_result.fit_fractions[el], plot_result.fit_fractions_err[el]));
                }
                info->Draw();
                
                TPad* pad2 = (TPad*)c_pdf->cd(2);
                pad2->SetPad(0, 0, 1, 0.3); pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.3); pad2->SetGridy();
                TGraphErrors* pullGraph = new TGraphErrors();
                calculatePull(frame.get(), pullGraph, frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
                setupPullPlot(pullGraph, frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
                pullGraph->Draw("AP");
                
                c_pdf->Update();
                c_pdf->Print(pdf_filename.c_str());
            }
            
            outputFile->cd();
            h_chi2ndf->Write();
            for (auto& pair : output_hists) {
                if(pair.second->GetEntries() > 0) pair.second->Write();
            }
        }
    }
    c_pdf->Print((pdf_filename + "]").c_str());
    outputFile->Close();
}

void ChargeTempFit() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    
    const vector<string> sources = {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
    for (const auto& src : sources) runUnifiedChargeAnalysis("UnbiasedL1Inner", "Pure", src, true);
    for (const auto& src : sources) runUnifiedChargeAnalysis("UnbiasedL1Inner", "Pure", src, false);
}