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

// --- 全局常量和结构体 (已移除 N 和 O) ---

const std::string outputDir = "/eos/user/z/zixuan/Isotope/ChargeTemp/";
const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};

const std::map<std::string, std::pair<double, double>> detector_ek_ranges = {
    {"TOF", {0.33, 1.29}},
    {"NaF", {0.90, 5.10}},
    {"AGL", {2.90, 21.0}}
};

struct ElementInfo { int Z; std::string name; };
// 仅保留 Be, B, C
const std::map<int, ElementInfo> element_db_by_z = {
    {4, {4, "Beryllium"}}, {5, {5, "Boron"}}, {6, {6, "Carbon"}}
};
const std::map<std::string, int> element_db_by_name = {
    {"Beryllium", 4}, {"Boron", 5}, {"Carbon", 6}
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

    std::map<std::string, std::map<std::string, double>> yields_in_window;
    std::map<std::string, std::map<std::string, double>> yields_in_window_err;
    
    std::map<std::string, std::map<std::string, double>> fractions_in_window;
    std::map<std::string, std::map<std::string, double>> fractions_in_window_err;
};

// --- UnifiedChargeFitter 类定义 ---

class UnifiedChargeFitter {
public:
    UnifiedChargeFitter(
        const std::string& chain, const std::string& detector, int ekpernuc_bin, 
        TH2F* signal_rebinned, const std::map<std::string, TH2F*>& templates_rebinned);

    bool initializeAndProject();
    bool runFit();
    UnifiedFitResult calculateAllYields();
    std::unique_ptr<RooPlot> generatePlotAndCalcChi2(UnifiedFitResult& result);
    const std::vector<std::string>& getTemplateElements() const { return templateElements_; }
    TH1D* getSignalHist() const { return h_signal_extended_.get(); }
    const std::map<std::string, std::unique_ptr<TH1D>>& getTemplateHists() const { return h_templates_extended_; }
    double getFitMin() const { return fitMin_; }
    double getFitMax() const { return fitMax_; }

private:
    std::string chainName_, detectorName_;
    int energyBin_;
    TH2F* h_signal_rebinned_;
    const std::map<std::string, TH2F*>& templates_rebinned_;
    
    std::vector<std::string> templateElements_;
    double fitMin_ = 3.7, fitMax_ = 6.3;

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

    void configureFit();
    std::unique_ptr<TH1D> projectSlice(TH2F* h2d, const char* name);
};

// --- UnifiedChargeFitter 的实现 ---

UnifiedChargeFitter::UnifiedChargeFitter(
    const std::string& chain, const std::string& detector, int ekpernuc_bin, 
    TH2F* signal_rebinned, const std::map<std::string, TH2F*>& templates_rebinned)
    : chainName_(chain), detectorName_(detector), energyBin_(ekpernuc_bin),
      h_signal_rebinned_(signal_rebinned), templates_rebinned_(templates_rebinned)
{
    TH1::AddDirectory(kFALSE);
    configureFit();
}

/**
 * @brief 配置拟合中使用的元素模板列表。
 * 限制为 Beryllium, Boron, Carbon。
 */
void UnifiedChargeFitter::configureFit() {
    // 严格限制为 Be, B, C
    templateElements_ = {"Beryllium", "Boron", "Carbon"}; 
}

std::unique_ptr<TH1D> UnifiedChargeFitter::projectSlice(TH2F* h2d, const char* name) {
    auto slice = std::unique_ptr<TH1D>(h2d->ProjectionX(name, energyBin_, energyBin_));
    return slice;
}

bool UnifiedChargeFitter::initializeAndProject() {
    auto h_signal_slice = projectSlice(h_signal_rebinned_, Form("h_signal_slice_bin%d", energyBin_));
    double N_signal_total = h_signal_slice ? h_signal_slice->GetEntries() : 0;
    
    // --- 调试检查 1: 信号统计量 ---
    if (!h_signal_slice || N_signal_total < 50) { 
        std::cerr << "        -> CRITICAL: Signal count is too low (" << N_signal_total << "). Minimum 50 required." << std::endl;
        return false;
    }
    
    std::map<std::string, std::unique_ptr<TH1D>> h_templates_slices;
    double N_templates_total_sum = 0.0;
    
    for (const auto& el : templateElements_) {
        auto it = templates_rebinned_.find(el);
        if (it == templates_rebinned_.end()) {
            std::cerr << "        -> CRITICAL: Template for " << el << " not found." << std::endl;
            return false;
        }
        h_templates_slices[el] = projectSlice(it->second, Form("h_template_%s_slice_bin%d", el.c_str(), energyBin_));
        double N_template = h_templates_slices[el] ? h_templates_slices[el]->GetEntries() : 0;
        
        // --- 对每个模板的 1D 切片执行 Smooth(1) ---
        if (h_templates_slices[el]) {
             //h_templates_slices[el]->Smooth(1);
        }
        // ------------------------------------------

        // --- 调试检查 2: 模板统计量 ---
        if (!h_templates_slices[el] || N_template < 10) { 
            std::cerr << "        -> CRITICAL: Template count for " << el << " is too low (" << N_template << "). Minimum 10 required." << std::endl;
        }
        N_templates_total_sum += N_template;
    }

    // 扩展和准备直方图
    h_signal_extended_ = extendHistogram(h_signal_slice.get(), fitMin_, fitMax_);
    for (const auto& el : templateElements_) {
        h_templates_extended_[el] = extendHistogram(h_templates_slices.at(el).get(), fitMin_, fitMax_);
    }

    charge_ = std::make_unique<RooRealVar>("charge", "Charge", fitMin_, fitMax_);
    data_hist_ = std::make_unique<RooDataHist>("data_hist", "Data", *charge_, h_signal_extended_.get());
    for (const auto& el : templateElements_) {
        // Clone for RooDataHist ownership
        auto hist_clone = (TH1D*)h_templates_extended_.at(el)->Clone(Form("%s_clone_for_rdh", h_templates_extended_.at(el)->GetName()));
        template_data_hists_[el] = std::make_unique<RooDataHist>(Form("dhist_%s", el.c_str()), "", *charge_, hist_clone);
        template_pdfs_[el] = std::make_unique<RooHistPdf>(Form("pdf_%s", el.c_str()), "", *charge_, *template_data_hists_[el]);
    }

    // ************ 【修改点】设置拟合参数：Be 和 B 为自由参数，C 为约束参数。 ************
    std::vector<std::string> free_params_elements = {"Beryllium", "Boron"};
    std::string constrained_element = "Carbon";
    
    std::map<std::string, RooAbsReal*> frac_map;
    
    // --- 设置自由参数 Be 和 B ---
    for (const auto& el : free_params_elements) {
        double N_template_i = h_templates_slices.at(el)->GetEntries();
        double initial_guess = (N_templates_total_sum > 0) ? N_template_i / N_templates_total_sum : 0.3;
        initial_guess = std::max(0.001, std::min(0.999, initial_guess));
        
        // 保持开放的范围 [0, 1]
        auto fracVar = std::make_unique<RooRealVar>(Form("frac_%s", el.c_str()), "", initial_guess, 0.0, 1.0);
        frac_map[el] = fracVar.get();
        fractionParams_.push_back(std::move(fracVar));
    }

    // --- 约束 Carbon ---
    std::string formula = "1.0";
    RooArgList formulaArgs;
    for (const auto& param : fractionParams_) { // 循环 Be 和 B 的参数
        formula += " - @" + std::to_string(formulaArgs.getSize());
        formulaArgs.add(*param);
    }
    // C = 1 - frac_Be - frac_B
    lastFraction_ = std::make_unique<RooFormulaVar>(Form("frac_%s", constrained_element.c_str()), "", formula.c_str(), formulaArgs);
    frac_map[constrained_element] = lastFraction_.get();

    RooArgList pdfList, fracList;
    for (const auto& el : templateElements_) {
        pdfList.add(*template_pdfs_.at(el));
        fracList.add(*frac_map.at(el));
    }
    // 非扩展拟合 (Non-extended fit)
    total_pdf_ = std::make_unique<RooAddPdf>("total_pdf", "Total PDF", pdfList, fracList, false); 
    return true;
}

bool UnifiedChargeFitter::runFit() {
    RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);    
    cout << "        -> Attempting fit with fixed range [" << fitMin_ << ", " << fitMax_ << "]..." << endl;
    charge_->setRange("fit_range", fitMin_, fitMax_);
    fitResult_ = std::unique_ptr<RooFitResult>(
        total_pdf_->fitTo(*data_hist_, Save(true), PrintLevel(3), Range("fit_range"), Strategy(2), Minimizer("Minuit2", "migrad"))
    );
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);    

    if (fitResult_ && fitResult_->status() == 0) {
        cout << "        -> Fit Successful!" << endl;
        return true;
    } else {
        cout << "        -> Fit failed to converge (Status: " << (fitResult_ ? fitResult_->status() : -1) << ")." << endl;
        return false;
    }
}

UnifiedFitResult UnifiedChargeFitter::calculateAllYields() {
    UnifiedFitResult res;

    if (!fitResult_ || fitResult_->status() != 0) {
        res.fit_status = fitResult_ ? fitResult_->status() : -1;
        return res;
    }
    res.fit_status = fitResult_->status();

    // 计算拟合比例（Fraction）
    for (const auto& el : templateElements_) { // 循环 Be, B, C
        RooAbsReal* frac_param = (RooAbsReal*)fitResult_->floatParsFinal().find(Form("frac_%s", el.c_str()));
        if (!frac_param) {      
            if (lastFraction_ && lastFraction_->GetName() == std::string("frac_" + el)) {
                // 对于约束参数 'Carbon'
                frac_param = lastFraction_.get();
            }
        }
        if (!frac_param) continue;

        double F_A = frac_param->getVal();
        res.fit_fractions[el] = F_A;
        
        double err_F_A = 0.0;
        if (frac_param->IsA()->InheritsFrom(RooRealVar::Class())) {
            err_F_A = ((RooRealVar*)frac_param)->getError();
        } else if (lastFraction_ && lastFraction_.get() == frac_param) {
            err_F_A = 0.0; // 约束参数的误差处理
        }
        res.fit_fractions_err[el] = err_F_A;
    }

    // 计算窄窗口内的 Yields 和 Fractions
    for (const auto& window_el_pair : element_db_by_z) { // 循环 Be, B, C 窗口
        const std::string& window_element_name = window_el_pair.second.name;
        int Z_window = window_el_pair.first;
        // 窄窗口定义
        double narrowMin = Z_window - 0.2;
        double narrowMax = Z_window + 0.4;
        
        if (narrowMax < fitMin_ || narrowMin > fitMax_) continue;

        std::string range_name = "narrow_range_" + window_element_name;
        charge_->setRange(range_name.c_str(), narrowMin, narrowMax);

        int bin_low = h_signal_extended_->GetXaxis()->FindBin(narrowMin);
        int bin_high = h_signal_extended_->GetXaxis()->FindBin(narrowMax);
        double N_sig_window = h_signal_extended_->Integral(bin_low, bin_high);
        double err_N_sig_window = (N_sig_window > 0) ? sqrt(N_sig_window) : 0.0;

        auto integral_R_total_obj = std::unique_ptr<RooAbsReal>(total_pdf_->createIntegral(*charge_, NormSet(*charge_), Range(range_name.c_str())));
        double R_total_window = integral_R_total_obj->getVal();
        double err_R_total_window = 0.0; 

        if (R_total_window <= 1e-9) continue;

        for (const auto& component_el : templateElements_) { // 循环 Be, B, C 组分
            double F_A = res.fit_fractions.at(component_el);
            double err_F_A = res.fit_fractions_err.at(component_el);

            auto integral_R_A_obj = std::unique_ptr<RooAbsReal>(template_pdfs_.at(component_el)->createIntegral(*charge_, NormSet(*charge_), Range(range_name.c_str())));
            double R_A = integral_R_A_obj->getVal();
            
            // 简单形状误差估算
            const auto& h_template = h_templates_extended_.at(component_el);
            double N_template_total = h_template->Integral();
            double k_template_narrow = h_template->Integral(h_template->GetXaxis()->FindBin(narrowMin), h_template->GetXaxis()->FindBin(narrowMax));
            double R_A_for_err = (N_template_total > 0) ? k_template_narrow / N_template_total : 0.0;
            double err_R_A = (N_template_total > 0) ? sqrt(std::max(0.0, R_A_for_err * (1.0 - R_A_for_err) / N_template_total)) : 0.0;

            // 组分 A 在窗口内的比例 P_A_window
            double P_A_window = F_A * R_A / R_total_window;
            double rel_err_sq_P_A = 0.0;
            if (F_A > 0) rel_err_sq_P_A += pow(err_F_A / F_A, 2);
            if (R_A > 0) rel_err_sq_P_A += pow(err_R_A / R_A, 2);
            if (R_total_window > 0) rel_err_sq_P_A += pow(err_R_total_window / R_total_window, 2);
            double err_P_A_window = P_A_window * sqrt(rel_err_sq_P_A);

            res.fractions_in_window[window_element_name][component_el] = P_A_window;
            res.fractions_in_window_err[window_element_name][component_el] = err_P_A_window;

            // 产额 Yield_A 
            double Yield_A = N_sig_window * P_A_window;
            double rel_err_sq_Yield_A = 0.0;
            if (N_sig_window > 0) rel_err_sq_Yield_A += pow(err_N_sig_window / N_sig_window, 2);
            if (P_A_window > 0) rel_err_sq_Yield_A += pow(err_P_A_window / P_A_window, 2);
            double err_Yield_A = Yield_A * sqrt(rel_err_sq_Yield_A);

            res.yields_in_window[window_element_name][component_el] = Yield_A;
            res.yields_in_window_err[window_element_name][component_el] = err_Yield_A;
        }
    }

    charge_->setRange("full_range", charge_->getMin(), charge_->getMax());
    return res;
}

std::unique_ptr<RooPlot> UnifiedChargeFitter::generatePlotAndCalcChi2(UnifiedFitResult& result) {
    auto frame = std::unique_ptr<RooPlot>(charge_->frame(Range(fitMin_, fitMax_), Title("Unified Charge Fit")));
    
    if (!fitResult_) {
        if(data_hist_) data_hist_->plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
        return frame;
    }

    data_hist_->plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
    total_pdf_->plotOn(frame.get(), Name("total_pdf"), LineColor(kRed), LineWidth(2));
    
    std::vector<int> colors = {kAzure + 7, kOrange - 3, kGreen + 2}; // 只有3个颜色
    for (size_t i = 0; i < templateElements_.size(); ++i) { // 循环 Be, B, C
        total_pdf_->plotOn(frame.get(), Components(*template_pdfs_.at(templateElements_[i])), Name(Form("comp_%s", templateElements_[i].c_str())), LineColor(colors[i % colors.size()]), LineWidth(2));
    }
    
    int nFreeParams = fitResult_->floatParsFinal().getSize();
    int nBinsInRange = 0;
    for(int i = 1; i <= h_signal_extended_->GetNbinsX(); ++i) {
        if (h_signal_extended_->GetBinCenter(i) >= fitMin_ && h_signal_extended_->GetBinCenter(i) <= fitMax_) nBinsInRange++;
    }
    result.ndf = nBinsInRange - nFreeParams;
    // 假设 calculateChi2 函数是外部定义的
    double chi2 = calculateChi2(frame.get(), "data_hist", "total_pdf", fitMin_, fitMax_);
    result.chi2ndf = (result.ndf > 0) ? chi2 / result.ndf : 0.0;

    return frame;
}

// --- runUnifiedChargeAnalysis 函数 (已移除 N 和 O) ---

void runUnifiedChargeAnalysis(const std::string& chain) {
    cout << "\n======================================================================\n";
    cout << "Starting Unified Charge Analysis for chain: " << chain << endl;
    cout << "Using only Beryllium, Boron, and Carbon templates.\n";
    cout << "======================================================================\n";
    
    // --- 信号数据文件 (保持不变) ---
    std::string signalInputFileName = "/eos/user/z/zixuan/Isotope/Add_NoB/Be_frag4.root";
    auto signalFile = std::unique_ptr<TFile>(TFile::Open(signalInputFileName.c_str()));    
    if (!signalFile || signalFile->IsZombie()) {    
        cerr << "CRITICAL: Could not open signal input file: " << signalInputFileName << endl;    
        return;    
    }
    
    // --- 模板数据文件 (新路径) ---
    std::string templateInputFileName = "/eos/user/z/zixuan/Isotope/PureChargeTemp/PureChargeTemplates_UnbiasedL1Inner.root";
    auto templateFile = std::unique_ptr<TFile>(TFile::Open(templateInputFileName.c_str()));    
    if (!templateFile || templateFile->IsZombie()) {    
        cerr << "CRITICAL: Could not open template input file: " << templateInputFileName << endl;    
        return;    
    }

    std::string pdf_filename = outputDir + "UnifiedQFit_" + chain + ".pdf";
    std::string root_filename = outputDir + "UnifiedQFit_" + chain + ".root";
    auto c_pdf = std::make_unique<TCanvas>("c_pdf", "PDF Canvas", 800, 600);
    c_pdf->Print((pdf_filename + "[").c_str());
    auto outputFile = std::make_unique<TFile>(root_filename.c_str(), "RECREATE");

    for (const auto& detector : detectors) {
        cout << "\n--- Processing Detector: " << detector << " ---" << endl;
        
        // 严格限制为 Be, B, C
        const std::vector<std::string> required_elements = {"Beryllium", "Boron", "Carbon"}; 

        // --- 信号直方图 (从 signalFile 获取) ---
        std::string signalHistName = chain + "_ISS_BKG_H2_Boron_L1QSignal_" + detector;
        TH2F* h_signal_raw = (TH2F*)signalFile->Get(signalHistName.c_str());
        if (!h_signal_raw) {    
            cout << "Signal histogram not found: " << signalHistName << ". Skipping detector." << endl;    
            continue;    
        }
        auto h_signal_rebinned = std::unique_ptr<TH2F>((TH2F*)h_signal_raw->Clone(Form("%s_rebinned", signalHistName.c_str())));
        h_signal_rebinned->RebinX(2);

        // --- 模板直方图 (从 templateFile 获取，仅限 Be, B, C) ---
        std::map<std::string, std::unique_ptr<TH2F>> templates_rebinned;
        bool all_templates_found = true;
        for (const auto& el : required_elements) {
            std::string templateHistName = Form("h2d_PureQTemp_%s_%s", el.c_str(), detector.c_str());
            TH2F* h_template_raw = (TH2F*)templateFile->Get(templateHistName.c_str());
            if (!h_template_raw) {    
                cerr << "CRITICAL: Could not find required template: " << templateHistName << endl;
                all_templates_found = false;    
                break;    
            }
            templates_rebinned[el] = std::unique_ptr<TH2F>((TH2F*)h_template_raw->Clone(Form("%s_rebinned", templateHistName.c_str())));
            templates_rebinned[el]->RebinX(1); 
        }
        if (!all_templates_found) { cout << "Missing templates. Skipping detector." << endl; continue; }

        int n_bins_y = h_signal_rebinned->GetNbinsY();
        const TAxis* y_axis = h_signal_rebinned->GetYaxis();
        auto h_chi2ndf = std::make_unique<TH1D>(Form("h_chi2ndf_%s", detector.c_str()), "", n_bins_y, y_axis->GetXbins()->GetArray());
        
        // 仅为 Be, B, C 初始化直方图
        std::map<std::string, std::unique_ptr<TH1D>> h_fitfracs;
        for (const auto& el : required_elements) {
            h_fitfracs[el] = std::make_unique<TH1D>(Form("h_fitfrac_%s_%s", el.c_str(), detector.c_str()), Form("Fit Fraction of %s in %s;E_{k} [GeV/n];Fraction", el.c_str(), detector.c_str()), n_bins_y, y_axis->GetXbins()->GetArray());
        }

        std::map<std::string, std::map<std::string, std::unique_ptr<TH1D>>> h_yields, h_fracs_in_window;
        for (const auto& window_el : required_elements) {
            for (const auto& comp_el : required_elements) {
                h_yields[window_el][comp_el] = std::make_unique<TH1D>(Form("h_yield_in_%s_from_%s_%s", window_el.c_str(), comp_el.c_str(), detector.c_str()), Form("Yield from %s in %s window (%s);E_{k} [GeV/n];Yield", comp_el.c_str(), window_el.c_str(), detector.c_str()), n_bins_y, y_axis->GetXbins()->GetArray());
                h_fracs_in_window[window_el][comp_el] = std::make_unique<TH1D>(Form("h_frac_in_%s_from_%s_%s", window_el.c_str(), comp_el.c_str(), detector.c_str()), Form("Fraction from %s in %s window (%s);E_{k} [GeV/n];Fraction", comp_el.c_str(), window_el.c_str(), detector.c_str()), n_bins_y, y_axis->GetXbins()->GetArray());
            }
        }
        
        // 调试用：保存切片直方图的目录
        outputFile->mkdir(Form("Debug_Slices_%s", detector.c_str()));
        
        for (int y_bin = 1; y_bin <= n_bins_y; ++y_bin) {
            double ek_center = y_axis->GetBinCenter(y_bin);
            double ek_low = y_axis->GetBinLowEdge(y_bin);
            cout<<"  [bin " << y_bin << "] Ek_low = " << ek_low << " GeV/n" << endl;
            const auto& range = detector_ek_ranges.at(detector);
            if (ek_center < range.first || ek_center > range.second) continue;

            cout << "    Processing energy bin " << y_bin << " (Ek=" << ek_center << " GeV/n)" << endl;
            std::map<std::string, TH2F*> templates_raw_ptr;
            for(auto const& [key, val] : templates_rebinned) templates_raw_ptr[key] = val.get();
            
            UnifiedChargeFitter fitter(chain, detector, y_bin, h_signal_rebinned.get(), templates_raw_ptr);
            if (!fitter.initializeAndProject()) {    
                cout << "    Initialization failed. Skipping bin." << endl;    
                continue;    
            }
            
            // --- 调试改进 5: 写入切片直方图进行目视检查 ---
            outputFile->cd(Form("Debug_Slices_%s", detector.c_str()));
            fitter.getSignalHist()->Write(Form("h_signal_E%d", y_bin));
            for(auto const& [name, hist] : fitter.getTemplateHists()) {
                hist->Write(Form("h_template_%s_E%d", name.c_str(), y_bin));
            }
            outputFile->cd();
            
            if (!fitter.runFit()) {    
                cout << "    Fit failed. Skipping bin." << endl;    
                // 绘制失败的图
                UnifiedFitResult failed_result;
                failed_result.ekpernuc_low = ek_low; failed_result.ekpernuc_up = y_axis->GetBinUpEdge(y_bin);
                auto frame = fitter.generatePlotAndCalcChi2(failed_result);
                c_pdf->Clear();
                auto text = std::make_unique<TPaveText>(0.1, 0.1, 0.9, 0.9);
                text->AddText(Form("Unified Fit (%s, E_{k}=%.2f-%.2f GeV/n)", detector.c_str(), failed_result.ekpernuc_low, failed_result.ekpernuc_up));
                text->AddText("FIT FAILED TO CONVERGE - SEE ROOT OUTPUT FOR DETAILS");
                frame->Draw();
                text->Draw("SAME");
                c_pdf->SetLogy();
                c_pdf->Print(pdf_filename.c_str());
                continue;    
            }
            
            UnifiedFitResult result = fitter.calculateAllYields();
            result.ekpernuc_center = ek_center;
            result.ekpernuc_low = ek_low;    
            result.ekpernuc_up = y_axis->GetBinUpEdge(y_bin);
            result.ekpernuc_width = y_axis->GetBinWidth(y_bin);
            auto frame = fitter.generatePlotAndCalcChi2(result);

            h_chi2ndf->SetBinContent(y_bin, result.chi2ndf);
            for (const auto& el : required_elements) { // 仅 Be, B, C
                h_fitfracs.at(el)->SetBinContent(y_bin, result.fit_fractions[el]);
                h_fitfracs.at(el)->SetBinError(y_bin, result.fit_fractions_err[el]);
                for (const auto& window_el : required_elements) { // 仅 Be, B, C 窗口
                    h_yields.at(window_el).at(el)->SetBinContent(y_bin, result.yields_in_window[window_el][el]);
                    h_yields.at(window_el).at(el)->SetBinError(y_bin, result.yields_in_window_err[window_el][el]);
                    h_fracs_in_window.at(window_el).at(el)->SetBinContent(y_bin, result.fractions_in_window[window_el][el]);
                    h_fracs_in_window.at(window_el).at(el)->SetBinError(y_bin, result.fractions_in_window_err[window_el][el]);
                }
            }
            
            // --- 绘图逻辑 (已清理 N/O 依赖) ---
            c_pdf->Clear();
            c_pdf->Divide(1, 2);
            TPad* pad1 = (TPad*)c_pdf->cd(1);
            pad1->SetPad(0, 0.3, 1, 1); 
            pad1->SetLogy(); 
            pad1->SetBottomMargin(0.02);
            
            frame->SetTitle(Form("Unified Fit (%s, E_{k}=%.2f-%.2f GeV/n)", detector.c_str(), result.ekpernuc_low, result.ekpernuc_up));
            frame->GetYaxis()->SetTitle("Events"); frame->GetXaxis()->SetLabelSize(0);
            frame->SetMinimum(9);
            auto hmax = h_signal_rebinned.get()->ProjectionX(Form("hmax%d", y_bin), y_bin, y_bin);
            frame->SetMaximum(5 * hmax->GetMaximum());
            frame->Draw();

            auto legend = std::make_unique<TLegend>(0.7, 0.55, 0.88, 0.88); 
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.03);
            legend->AddEntry("data_hist", "Data", "pe");
            legend->AddEntry("total_pdf", "Total Fit", "l");
            for (const auto& el : fitter.getTemplateElements()) { // 仅 Be, B, C
                legend->AddEntry(Form("comp_%s", el.c_str()), el.c_str(), "l");
            }
            legend->Draw(); 

            auto info = std::make_unique<TPaveText>(0.15, 0.65, 0.65, 0.88, "NDC");
            info->SetFillStyle(0); info->SetBorderSize(0); info->SetTextAlign(12);
            info->SetTextSize(0.03);
            info->AddText(Form("#chi^{2}/NDF = %.2f", result.chi2ndf));
            string fracs_line = "";
            for(const auto& pair : result.fit_fractions) { // 仅 Be, B, C
                const std::string& el = pair.first;
                fracs_line += Form("F_{%s}=%.3f; ", el.c_str(), pair.second);
            }
            info->AddText(fracs_line.c_str());
            
            // 确保窗口计算仅使用 Be, B, C 的结果
            double yield_Be_in_Be = result.yields_in_window.at("Boron").at("Beryllium");
            double err_Be_in_Be = result.yields_in_window_err.at("Boron").at("Beryllium");
            double yield_B_in_Be = result.yields_in_window.at("Boron").at("Boron");
            double err_B_in_Be = result.yields_in_window_err.at("Boron").at("Boron");
            info->AddText(Form("In B window: N_{Be}=%.1f#pm%.1f, N_{B}=%.1f#pm%.1f", yield_Be_in_Be, err_Be_in_Be, yield_B_in_Be, err_B_in_Be));
            info->Draw();
            
            TPad* pad2 = (TPad*)c_pdf->cd(2);
            pad2->SetPad(0, 0, 1, 0.3); pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.3); pad2->SetGridy();
            auto pullGraph = std::make_unique<TGraphErrors>();
            // 假设 calculatePull 和 setupPullPlot 是外部定义的
            calculatePull(frame.get(), pullGraph.get(), frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
            setupPullPlot(pullGraph.get(), frame->GetXaxis()->GetXmin(), frame->GetXaxis()->GetXmax());
            pullGraph->Draw("AP");
            TLine zeroLine(frame->GetXaxis()->GetXmin(), 0, frame->GetXaxis()->GetXmax(), 0);
            zeroLine.SetLineStyle(2); zeroLine.SetLineColor(kGray + 1);
            zeroLine.Draw("SAME");
            c_pdf->Print(pdf_filename.c_str());
        }
        
        // 写入结果直方图 (仅 Be, B, C)
        outputFile->cd();
        outputFile->mkdir(detector.c_str())->cd();
        h_chi2ndf->Write();
        for (const auto& el : required_elements) {
            h_fitfracs.at(el)->Write();
            for (const auto& window_el : required_elements) {
                h_yields.at(window_el).at(el)->Write();
                h_fracs_in_window.at(window_el).at(el)->Write();
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
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    
    runUnifiedChargeAnalysis("UnbiasedL1Inner");
}