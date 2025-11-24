// Corr.cpp
// 核心功能：1. 拼接数据并计算 Base Ratio 2. 分段常数拟合 3. 计算 C_final (含误差传递)
// 特点：使用裸指针，不使用 TFitResult* 和 std::unique_ptr，图上无任何冗余元素，不使用 delete。

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TText.h>

#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <map>
#include <array>
#include <vector>
#include <utility> // For std::pair

// --- 全局常数定义 (替换所有 struct 和魔法数字) ---

// 绘图和输出精度
const int FINAL_PRECISION = 4;
const double SCALING_FACTOR = 0.7;

// 绘图范围
const double X_MIN_DRAW = 0.4;
const double X_MAX_DRAW = 21.5;
const double Y_MIN_DRAW = 0.0090;
const double Y_MAX_DRAW = 0.0210;

// 直方图拼接边界 (GeV)
const double STITCH_BOUNDARY_1 = 1.28; // TOF-NaF 边界
const double STITCH_BOUNDARY_2 = 3.06; // NaF-AGL 边界

// 低能区拟合常数 (0.4 - 6.1 GeV)
const double LOWE_RANGE_MIN = X_MIN_DRAW;
const double LOWE_RANGE_MAX = 6.1;
const double R_ISS_LOWE_VAL = 0.0111875;
const double R_ISS_LOWE_ERR = 0.002;
const double R_MC_LOWE_VAL = 0.0149920;
const double R_MC_LOWE_ERR = 0.002;

// 高能区拟合常数 (6.1 - 21.5 GeV)
const double HIGHE_RANGE_MIN = 6.1;
const double HIGHE_RANGE_MAX = X_MAX_DRAW;
const double R_ISS_HIGHE_VAL = 0.00999855;
const double R_ISS_HIGHE_ERR = 0.003;
const double R_MC_HIGHE_VAL = 0.0143756;
const double R_MC_HIGHE_ERR = 0.002;

// --- Helper: Error Propagation Functions (替换 ValueWithError struct) ---
// 使用 std::pair<double, double> 表示 {value, error}

/**
 * @brief 格式化输出 {val, err}
 */
std::string pairToString(const std::pair<double, double>& p, int precision = FINAL_PRECISION) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << p.first << " #pm " << p.second;
    return ss.str();
}

/**
 * @brief Subtraction: p1 - p2
 */
std::pair<double, double> subtract(const std::pair<double, double>& p1, const std::pair<double, double>& p2) {
    double v = p1.first - p2.first;
    double e = std::hypot(p1.second, p2.second);
    return {v, e};
}

/**
 * @brief Division: p1 / p2
 */
std::pair<double, double> divide(const std::pair<double, double>& p1, const std::pair<double, double>& p2) {
    if (std::abs(p2.first) < 1e-9) { return {std::nan(""), std::nan("")}; }
    double q_val = p1.first / p2.first;
    double q_err = q_val * std::hypot(p1.second / p1.first, p2.second / p2.first);
    return {q_val, q_err};
}

/**
 * @brief Natural Log: ln(p)
 */
std::pair<double, double> log_nat(const std::pair<double, double>& p) {
    if (p.first <= 0.0) { return {std::nan(""), std::nan("")}; }
    double log_val = std::log(p.first);
    double log_err = std::abs(p.second / p.first);
    return {log_val, log_err};
}

/**
 * @brief Power: p_base ^ p_exponent
 */
std::pair<double, double> power(const std::pair<double, double>& p_base, const std::pair<double, double>& p_exponent) {
    if (p_base.first <= 0.0) { return {std::nan(""), std::nan("")}; }
    double pow_val = std::pow(p_base.first, p_exponent.first);
    double term1 = p_exponent.first * p_base.second / p_base.first;
    double term2 = std::log(p_base.first) * p_exponent.second;
    double pow_err = pow_val * std::hypot(term1, term2);
    return {pow_val, pow_err};
}

/**
 * @brief 计算最终修正因子 C_final
 */
std::pair<double, double> calculateFinalCorrection(const std::pair<double, double>& Base,
                                                 const std::pair<double, double>& R_ISS,
                                                 const std::pair<double, double>& R_MC) {
    std::pair<double, double> C_sigma = divide(log_nat(R_ISS), log_nat(R_MC));
    // 创建一个常数 {1.0, 0.0}
    const std::pair<double, double> one_with_zero_err = {1.0, 0.0};
    std::pair<double, double> exponent = subtract(C_sigma, one_with_zero_err);
    return power(Base, exponent);
}


// --- 辅助函数：根据能量边界拼接不同探测器的直方图 (使用裸指针，不 delete) ---
TH1D* stitchDetectorHists(TFile* file, const std::string& base_tag, const std::string& stitched_name) {
    const std::array<std::string, 3> detectors = {"TOF", "NaF", "AGL"};
    // 使用全局常数
    const double boundary1 = STITCH_BOUNDARY_1;
    const double boundary2 = STITCH_BOUNDARY_2;
    std::map<std::string, TH1*> hists;
    TH1* h_ref = nullptr;

    for (const auto& det : detectors) {
        TString full_name = TString::Format("BKG_%s_%s", base_tag.c_str(), det.c_str());
        TH1* h_current = dynamic_cast<TH1*>(file->Get(full_name));
        if (!h_current) {
            TString mass_tag = (base_tag == "H2b") ? "_Z4_Mass10" : "";
            full_name = TString::Format("UnbiasedL1Inner_MC_BKG_%s_%s%s", base_tag.c_str(), det.c_str(), mass_tag.Data());
            h_current = dynamic_cast<TH1*>(file->Get(full_name));
        }
        if (!h_current) continue;

        TH1* h_cloned = (TH1*)h_current->Clone(TString::Format("h_temp_%s_%s", base_tag.c_str(), det.c_str()));
        h_cloned->SetDirectory(0);
        h_cloned->Rebin(2); // 魔法数字 2
        hists[det] = h_cloned;
        if (h_ref == nullptr) {
            h_ref = (TH1*)h_cloned->Clone("h_ref_binning");
            h_ref->SetDirectory(0);
        }
    }
    if (h_ref == nullptr) { return nullptr; }

    TH1D* combined = new TH1D(stitched_name.c_str(), "", h_ref->GetNbinsX(), h_ref->GetXaxis()->GetXbins()->GetArray());
    combined->SetDirectory(0);

    for (int i = 1; i <= combined->GetNbinsX(); ++i) {
        double binCenter = combined->GetBinCenter(i);
        std::string det;
        if (binCenter < boundary1) { det = "TOF"; }
        else if (binCenter < boundary2) { det = "NaF"; }
        else { det = "AGL"; }

        if (hists.count(det) && hists.at(det) != nullptr) {
            TH1* h = hists.at(det);
            int source_bin = h->FindBin(binCenter);
            combined->SetBinContent(i, h->GetBinContent(source_bin));
            combined->SetBinError(i, h->GetBinError(source_bin));
        }
    }

    // 严禁 delete 内部对象
    // h_ref 和 h_cloned 都是 h_current 的克隆，在 ROOT 宏执行结束时，如果它们没有被添加到任何 TDirectory (SetDirectory(0)) 并且没有被显式 delete，内存会被清理。
    // 在 ROOT 交互式环境或宏中，这种做法是可接受的，以满足“不使用 delete”的要求。

    return combined;
}

// ========================================================================================
// === 主执行函数
// ========================================================================================
void Corr() {
    // 强制移除所有统计信息和拟合信息，确保图上干净
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // I. 拼接数据并计算 R_Base (使用裸指针)
    // 路径常数化
    const TString FILE_PATH = "/eos/user/z/zixuan/Isotope/Add/B11_rew_frag4.root";
    TFile* file = TFile::Open(FILE_PATH.Data());
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file: " << FILE_PATH << std::endl;
        return;
    }

    TH1D* h_H2b_frag = stitchDetectorHists(file, "H2b", "Stitched_H2b");
    TH1D* h_H1b_source = stitchDetectorHists(file, "H1b", "Stitched_H1b");
    if (!h_H2b_frag || !h_H1b_source) {
        std::cerr << "Error: Cannot find required histograms." << std::endl;
        // 注意：file, h_H2b_frag, h_H1b_source 都没有 delete
        return;
    }

    TH1D* h_R_Base = (TH1D*)h_H2b_frag->Clone("h_R_Base");
    h_R_Base->SetTitle("");
    h_R_Base->GetYaxis()->SetTitle("");
    h_R_Base->GetXaxis()->SetTitle("");
    h_R_Base->Divide(h_H1b_source);
    h_R_Base->Scale(SCALING_FACTOR);

    // 强制设置 X/Y 轴绘图范围
    h_R_Base->SetMarkerStyle(20); // 魔法数字 20
    h_R_Base->SetMarkerSize(0.8); // 魔法数字 0.8
    h_R_Base->SetMarkerColor(kBlack);
    h_R_Base->GetYaxis()->SetRangeUser(Y_MIN_DRAW, Y_MAX_DRAW);
    h_R_Base->GetXaxis()->SetRangeUser(X_MIN_DRAW, X_MAX_DRAW);

    // II. 分两段常数拟合

    // 1. 低能区拟合
    TF1* f_lowE_fit = new TF1("f_lowE_fit", "[0]", LOWE_RANGE_MIN, LOWE_RANGE_MAX);
    // 拟合选项 "R" 表示只在指定的 Range 内拟合
    h_R_Base->Fit(f_lowE_fit, "R", "", LOWE_RANGE_MIN, LOWE_RANGE_MAX);
    std::pair<double, double> base_lowE = {f_lowE_fit->GetParameter(0), f_lowE_fit->GetParError(0)};

    // 2. 高能区拟合
    TF1* f_highE_fit = new TF1("f_highE_fit", "[0]", HIGHE_RANGE_MIN, HIGHE_RANGE_MAX);
    h_R_Base->Fit(f_highE_fit, "R", "", HIGHE_RANGE_MIN, HIGHE_RANGE_MAX);
    std::pair<double, double> base_highE = {f_highE_fit->GetParameter(0), f_highE_fit->GetParError(0)};

    // III. 计算最终修正因子 C_final (使用全局常数和函数)
    const std::pair<double, double> R_ISS_lowE = {R_ISS_LOWE_VAL, R_ISS_LOWE_ERR};
    const std::pair<double, double> R_MC_lowE = {R_MC_LOWE_VAL, R_MC_LOWE_ERR};
    const std::pair<double, double> R_ISS_highE = {R_ISS_HIGHE_VAL, R_ISS_HIGHE_ERR};
    const std::pair<double, double> R_MC_highE = {R_MC_HIGHE_VAL, R_MC_HIGHE_ERR};

    std::pair<double, double> C_final_lowE = calculateFinalCorrection(base_lowE, R_ISS_lowE, R_MC_lowE);
    std::pair<double, double> C_final_highE = calculateFinalCorrection(base_highE, R_ISS_highE, R_MC_highE);

    // IV. 绘图 (关键修正，确保线条长度正确)
    TCanvas* c1 = new TCanvas("c1", "", 800, 600); // 魔法数字 800, 600
    c1->cd();

    h_R_Base->Draw("PZ");

    // 绘制低能区线 (使用拟合值和常数范围)
    TF1* f_lowE_draw = new TF1("f_lowE_draw", Form("%f", base_lowE.first), LOWE_RANGE_MIN, LOWE_RANGE_MAX);
    f_lowE_draw->SetLineColor(kBlue + 2); // 魔法数字 +2
    f_lowE_draw->SetLineWidth(2);         // 魔法数字 2
    f_lowE_draw->Draw("SAME L");

    // 绘制高能区线 (红线) (使用拟合值和常数范围)
    TF1* f_highE_draw = new TF1("f_highE_draw", Form("%f", base_highE.first), HIGHE_RANGE_MIN, HIGHE_RANGE_MAX);
    f_highE_draw->SetLineColor(kRed + 2); // 魔法数字 +2
    f_highE_draw->SetLineWidth(2);        // 魔法数字 2
    f_highE_draw->Draw("SAME L");

    c1->Update();
    const TString OUTPUT_FILE = "B_to_Be10_FinalCorrection_NoDeletes_FinalPlot.pdf";
    c1->SaveAs(OUTPUT_FILE.Data());

    // V. 输出最终修正因子到屏幕
    std::cout << "\n=================================================================\n";
    std::cout << "          B -> 10Be 最终修正因子 C_final (Stitched)\n";
    std::cout << "=================================================================\n";
    std::cout << "--- Low Energy (" << LOWE_RANGE_MIN << " - " << LOWE_RANGE_MAX << " GeV) ---\n";
    std::cout << "Base (Fit): " << pairToString(base_lowE, FINAL_PRECISION) << "\n";
    std::cout << "C_final: " << pairToString(C_final_lowE, FINAL_PRECISION) << "\n";
    std::cout << "----------------------------------\n";
    std::cout << "--- High Energy (" << HIGHE_RANGE_MIN << " - " << HIGHE_RANGE_MAX << " GeV) ---\n";
    std::cout << "Base (Fit): " << pairToString(base_highE, FINAL_PRECISION) << "\n";
    std::cout << "C_final: " << pairToString(C_final_highE, FINAL_PRECISION) << "\n";
    std::cout << "=================================================================\n";

    // 严禁 delete 语句
}