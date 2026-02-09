#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdio>

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TVectorD.h>
#include <TString.h>
#include <TMath.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TSystem.h>

#include "../Tool.h" // 确保这里面包含了 funcExpGausExp 的定义

using namespace AMS_Iso;

// 函数声明
void buildCDFLookupTables(const std::string& nucleusName = "Beryllium");

// 全局常量
const double qmin_global = 1.5, qmax_global = 9.5;

// EGE 参数列表 (不再包含 LG)
const std::vector<std::pair<std::string, int>> EGE_PARAM_LIST = {
    {"Peak", 0}, {"SigmaL", 1}, {"AlphaL", 2},
    {"SigmaR", 3}, {"AlphaR", 4}, {"Norm", 5}, {"xmin", 6}, {"xmax", 7}
};

const std::map<std::string, int> ELEMENTS = {
    {"Helium", 2}, {"Lithium", 3}, {"Beryllium", 4}, {"Boron", 5}, 
    {"Carbon", 6}, {"Nitrogen", 7}, {"Oxygen", 8}
};
const std::vector<std::string> DETECTORS = {"TOF", "NaF", "AGL"};
const std::vector<std::string> CHAINS = {"L1Inner", "UnbiasedL1Inner"};
// 【确认修改 1】这里已经按照要求去掉了 Q，直接 L1Template, L2Template
const std::vector<std::string> TEMPLATES = {"L1Template", "L2Template"}; 

// 探测器能量范围
const std::map<std::string, std::pair<double, double>> DETECTOR_RANGES = {
    {"TOF", {0.2, 1.8}}, {"NaF", {0.61, 6.10}}, {"AGL", {2.70, 30.0}}
};

// 安全外推范围
const std::map<std::string, std::pair<double, double>> SAFE_DETECTOR_RANGES = {
    {"TOF", {0.25, 1.6}}, {"NaF", {0.7, 5.8}}, {"AGL", {2.8, 25.9}}
};

// 参数限制 (已移除 Width 和 Sigma，只保留 EGE 相关)
const std::map<std::string, std::pair<double, double>> PARAM_LIMITS = {
    {"SigmaL", {0.01,  0.5}},
    {"AlphaL", {0.1,   4.0}},
    {"SigmaR", {0.01,  0.5}},
    {"AlphaR", {0.1,   4.0}}
};

// 查找表结构体
struct LookupTable {
    std::vector<double> q_values;
    std::vector<double> cdf_l1;
    std::vector<double> cdf_l2;
    
    void build(TF1& f_l1, TF1& f_l2, double qmin, double qmax, const std::string& context, int npts = 4000) {
        if (npts <= 0) { std::cerr << "Error: npts must be positive." << std::endl; return; }
        q_values.resize(npts + 1); cdf_l1.resize(npts + 1); cdf_l2.resize(npts + 1);
        
        f_l1.SetRange(qmin, qmax); f_l2.SetRange(qmin, qmax);
        double total_l1 = f_l1.Integral(qmin, qmax, 1e-10);
        double total_l2 = f_l2.Integral(qmin, qmax, 1e-10);
        
        std::cout << "  [Debug] For context: " << context << std::endl;
        std::cout << "          Total Integral L1: " << total_l1 << ", Total Integral L2: " << total_l2 << std::endl;
        std::cout << "          Params for " << f_l1.GetName() << " (L1): [ ";
        for (int i = 0; i < f_l1.GetNpar(); ++i) std::cout << "p" << i << "=" << f_l1.GetParameter(i) << " ";
        std::cout << "]" << std::endl;
        std::cout << "          Params for " << f_l2.GetName() << " (L2): [ ";
        for (int i = 0; i < f_l2.GetNpar(); ++i) std::cout << "p" << i << "=" << f_l2.GetParameter(i) << " ";
        std::cout << "]" << std::endl;
        if (total_l1 <= 0 || total_l2 <= 0) {
            std::cerr << "  [FATAL WARNING] Integral of a PDF is zero or negative." << std::endl;
        }
        
        double dq = (qmax - qmin) / npts;
        for (int i = 0; i <= npts; ++i) {
            double q = qmin + i * dq;
            q_values[i] = q;
            double int1 = (total_l1 > 0) ? f_l1.Integral(qmin, q, 1e-9) : 0;
            double int2 = (total_l2 > 0) ? f_l2.Integral(qmin, q, 1e-9) : 0;
            cdf_l1[i] = (total_l1 > 0) ? int1 / total_l1 : 0;
            cdf_l2[i] = (total_l2 > 0) ? int2 / total_l2 : 0;
        }
    }
};

// 预加载 Splines (只加载 EGE)
void preloadSplines(std::map<std::string, std::unique_ptr<TF1>>& cache, TFile* fin, const std::string& nucleusName) {
    std::cout << "Pre-loading splines for " << nucleusName << " (EGE only)..." << std::endl;
    for (const auto& chain : CHAINS) {
        for (const auto& det : DETECTORS) {
            for (const auto& temp : TEMPLATES) {
                // 不再循环 model，直接指定 "EGE"
                std::string model = "EGE"; 
                for (const auto& p : EGE_PARAM_LIST) {
                    // 【确认修改 1】直接使用 vector 中的字符串，不加 Q
                    std::string tempNameInFile = temp; 
                    
                    std::string splineName = chain + "_" + nucleusName + "_" + det + "_" + tempNameInFile + "_" + model + "_" + p.first + "_spline";
                    TF1* spline = dynamic_cast<TF1*>(fin->Get(splineName.c_str()));
                    if (spline) {
                        cache[splineName] = std::unique_ptr<TF1>(static_cast<TF1*>(spline->Clone()));
                    }
                }
            }
        }
    }
    std::cout << "Finished pre-loading " << cache.size() << " splines into memory." << std::endl;
}

// 从 Spline 获取参数 (只针对 EGE)
void getParamsFromSpline(const std::map<std::string, std::unique_ptr<TF1>>& splineCache, 
                         const std::string& chain, const std::string& elem, const std::string& det, 
                         const std::string& tempNameInFile, double ek, double* params, int z = 0) {
    
    // 固定使用 EGE 模型
    std::string model = "EGE";
    const auto& det_range = DETECTOR_RANGES.at(det);
    bool is_extrapolating = (ek < det_range.first || ek > det_range.second);

    for (const auto& p : EGE_PARAM_LIST) {
        const std::string& param_name = p.first;
        int param_index = p.second;
        std::string splineName = chain + "_" + elem + "_" + det + "_" + tempNameInFile + "_" + model + "_" + param_name + "_spline";
        
        auto it = splineCache.find(splineName);
        if (it != splineCache.end()) {
            TF1* spline = it->second.get();
            double val = spline->Eval(ek);

            if (is_extrapolating || val <= 0) {
                double min_lim = 0, max_lim = 0;
                bool has_limits = false;

                if (param_name == "Peak") {
                    min_lim = z - 0.25;
                    max_lim = z + 0.25;
                    has_limits = true;
                } else {
                    auto limit_it = PARAM_LIMITS.find(param_name);
                    if (limit_it != PARAM_LIMITS.end()) {
                        min_lim = limit_it->second.first;
                        max_lim = limit_it->second.second;
                        has_limits = true;
                    }
                }

                if (has_limits && (val < min_lim || val > max_lim)) {
                    double old_val = val;
                    const auto& safe_range = SAFE_DETECTOR_RANGES.at(det);
                    double eval_ek = (ek < det_range.first) ? safe_range.first : safe_range.second;
                    val = spline->Eval(eval_ek);
                    printf("  [Info] Corrected param '%s' at Ek=%.3f. Extrapolated val %.4f was out of bounds. Using val %.4f from SAFE Ek boundary %.3f\n",
                           param_name.c_str(), ek, old_val, val, eval_ek);
                }
            }

            params[param_index] = val;

            if (val <= 0 && param_name != "xmin" && param_name != "xmax") { 
                printf("  [Value Warning] Param '%s' is non-positive (%.4f) at Ek=%.3f for %s/%s/%s/%s.\n",
                       param_name.c_str(), val, ek, chain.c_str(), elem.c_str(), det.c_str(), tempNameInFile.c_str());
            }

        } else {
            params[param_index] = 0;
        }
    }

    // EGE 特定固定参数
    params[5] = 1.0;          // Norm
    params[6] = qmin_global;  // xmin
    params[7] = qmax_global;  // xmax
}

void buildCDFLookupTables(const std::string& nucleusName) {
    if (ELEMENTS.find(nucleusName) == ELEMENTS.end()) {
        std::cerr << "Error: Invalid nucleus name '" << nucleusName << "'." << std::endl;
        return;
    }
    const int z = ELEMENTS.at(nucleusName);

    const std::clock_t start_clock = std::clock();
    const time_t start_time = time(nullptr);
    std::cout << "Starting performance monitoring..." << std::endl;

    const std::string paramFileName = "/eos/user/z/zixuan/Isotope/ChargeFit/smooth/withBkg_ChargeFitParamsSmooth_HeToOxy_NoTune_iter2.root";
    const std::string outFileName = "/eos/user/z/zixuan/Isotope/L2QTuning/NoBkg_CDFLookupTable_fromSpline_" + nucleusName + ".root";
    
    // 读取 Binning 信息 (保持原样)
    const std::string binningFileName = "/eos/user/z/zixuan/Isotope/ChargeFit/NoBkg_ChargeFitParams_HeToOxy_NoTune_iter2.root";
    const std::string binningHistName = "L1Inner_Beryllium_AGL_L1Template_EGE_Peak";
    std::vector<double> energyBins;
    auto finBinning = std::unique_ptr<TFile>(TFile::Open(binningFileName.c_str()));
    if (!finBinning || finBinning->IsZombie()) { std::cerr << "Error: Failed to open binning file: " << binningFileName << std::endl; return; }
    TH1D* hBinning = dynamic_cast<TH1D*>(finBinning->Get(binningHistName.c_str()));
    if (!hBinning) { std::cerr << "Error: Failed to get binning histogram: " << binningHistName << std::endl; return; }
    for (int i = 1; i <= hBinning->GetXaxis()->GetNbins() + 1; ++i) energyBins.push_back(hBinning->GetXaxis()->GetBinLowEdge(i));
    finBinning->Close();
    std::cout << "Successfully loaded " << energyBins.size() - 1 << " energy bins." << std::endl;

    auto finParam = std::unique_ptr<TFile>(TFile::Open(paramFileName.c_str()));
    if (!finParam || finParam->IsZombie()) { std::cerr << "Error: Failed to open parameter file: " << paramFileName << std::endl; return; }
    
    std::map<std::string, std::unique_ptr<TF1>> splineCache;
    preloadSplines(splineCache, finParam.get(), nucleusName);
    finParam->Close();

    auto fout = std::unique_ptr<TFile>(TFile::Open(outFileName.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) { std::cerr << "Error: Failed to create output file: " << outFileName << std::endl; return; }
    std::cout << "Starting to build CDF lookup tables (EGE ONLY) for: " << nucleusName << " -> " << outFileName << std::endl;
    
    for (const auto& chain : CHAINS) {
        for (const auto& det : DETECTORS) {
            std::cout << "\nProcessing: " << chain << " / " << nucleusName << " / " << det << " ..." << std::endl;
            const size_t total_bins_in_loop = energyBins.size() - 1;

            for (size_t iy = 0; iy < total_bins_in_loop; ++iy) {
                double ekLow = energyBins[iy];
                double ekCenter = (energyBins[iy] + energyBins[iy+1]) / 2.0;

                // 简单的能量切割
                if ((det == "TOF" && (ekCenter < 0.2 || ekCenter > 1.8)) ||
                    (det == "NaF" && (ekCenter < 0.61 || ekCenter > 6.1)) ||
                    (det == "AGL" && (ekCenter < 2.7 || ekCenter > 30.0))) continue;
                
                printf("--> Processing Bin %zu (Ek_low = %.4f, Ek_center = %.4f)\n", iy, ekLow, ekCenter);

                double par_ege_l1[8] = {0}, par_ege_l2[8] = {0};
                
                // 【确认修改 1】直接传递 "L1Template" / "L2Template"
                getParamsFromSpline(splineCache, chain, nucleusName, det, "L1Template", ekCenter, par_ege_l1, z);
                getParamsFromSpline(splineCache, chain, nucleusName, det, "L2Template", ekCenter, par_ege_l2, z);

                // 定义函数
                TF1 ege_l1("ege_l1", funcExpGausExp, qmin_global, qmax_global, 8); ege_l1.SetParameters(par_ege_l1); 
                TF1 ege_l2("ege_l2", funcExpGausExp, qmin_global, qmax_global, 8); ege_l2.SetParameters(par_ege_l2);
                
                LookupTable ege_lookup;
                std::string ege_context = Form("EGE model for %s, %s, %s, bin %zu", chain.c_str(), nucleusName.c_str(), det.c_str(), iy);
                ege_lookup.build(ege_l1, ege_l2, qmin_global, qmax_global, ege_context);
                
                fout->cd();
                std::string baseName = Form("%s_%s_%s_bin%zu", chain.c_str(), nucleusName.c_str(), det.c_str(), iy);
                
                // 保存 EGE CDF 数据
                TVectorD(ege_lookup.q_values.size(), ege_lookup.q_values.data()).Write((baseName + "_EGE_q").c_str());
                TVectorD(ege_lookup.cdf_l1.size(), ege_lookup.cdf_l1.data()).Write((baseName + "_EGE_cdf_l1").c_str());
                TVectorD(ege_lookup.cdf_l2.size(), ege_lookup.cdf_l2.data()).Write((baseName + "_EGE_cdf_l2").c_str());

                // 将 double 数组转换为 TVectorD 并保存
                TVectorD vec_par_l1(8, par_ege_l1);
                TVectorD vec_par_l2(8, par_ege_l2);
                vec_par_l1.Write((baseName + "_EGE_pars_l1").c_str());
                vec_par_l2.Write((baseName + "_EGE_pars_l2").c_str());

                std::clock_t now_clock = std::clock();
                time_t now_time = time(nullptr);
                double cpu_sec = double(now_clock - start_clock) / CLOCKS_PER_SEC;
                double wall_sec = difftime(now_time, start_time);
                double mem_mb = getCurrentRSS_MB(); 

                printf("  [Progress] Bin %3zu/%zu (%s/%s) done. | Total CPU: %8.2f s | Total Wall: %8.2f s | Mem: %8.2f MB\n",
                       iy + 1, total_bins_in_loop, chain.c_str(), det.c_str(), cpu_sec, wall_sec, mem_mb);
            }
        }
    }
    
    fout->Close();
    std::cout << "\nSuccessfully built and saved all CDF lookup tables (EGE Only) to: " << outFileName << std::endl;
}