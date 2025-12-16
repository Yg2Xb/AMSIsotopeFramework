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

#include "../Tool.h"

using namespace AMS_Iso;

void buildCDFLookupTables(const std::string& nucleusName = "Beryllium");

const double qmin_global = 1.5, qmax_global = 8.5;

const std::vector<std::pair<std::string, int>> LG_PARAM_LIST = {
    {"Width", 0}, {"MPV", 1}, {"Area", 2}, {"Sigma", 3}
};
const std::vector<std::pair<std::string, int>> EGE_PARAM_LIST = {
    {"Peak", 0}, {"SigmaL", 1}, {"AlphaL", 2},
    {"SigmaR", 3}, {"AlphaR", 4}, {"Norm", 5}, {"xmin", 6}, {"xmax", 7}
};
const std::map<std::string, int> ELEMENTS = {{"Helium", 2},{"Lithium", 3},{"Beryllium", 4}, {"Boron", 5}, {"Carbon", 6}, {"Nitrogen", 7}, {"Oxygen", 8}};
const std::vector<std::string> DETECTORS = {"TOF", "NaF", "AGL"};
const std::vector<std::string> CHAINS = {"L1Inner", "UnbiasedL1Inner"};
const std::vector<std::string> TEMPLATES = {"L1Template", "L2Template"};
const std::map<std::string, std::pair<double, double>> DETECTOR_RANGES = {
    {"TOF", {0.3, 1.28}}, {"NaF", {0.71, 5.1}}, {"AGL", {2.8, 20.0}}
};

const std::map<std::string, std::pair<double, double>> PARAM_LIMITS = {
    {"Width",  {0.005, 0.11}},
    {"Sigma",  {0.01,  0.5}},
    {"SigmaL", {0.01,  0.5}},
    {"AlphaL", {0.1,   4.0}},
    {"SigmaR", {0.01,  0.5}},
    {"AlphaR", {0.1,   4.0}}
};

struct LookupTable {
    std::vector<double> q_values;
    std::vector<double> cdf_l1;
    std::vector<double> cdf_l2;
    
    void build(TF1& f_l1, TF1& f_l2, double qmin, double qmax, const std::string& context, int npts = 1400) {
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
            std::cerr << "  [FATAL WARNING] Integral of a PDF is zero or negative. See details above." << std::endl;
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

void preloadSplines(std::map<std::string, std::unique_ptr<TF1>>& cache, TFile* fin, const std::string& nucleusName) {
    std::cout << "Pre-loading splines for " << nucleusName << "..." << std::endl;
    for (const auto& chain : CHAINS) {
        for (const auto& det : DETECTORS) {
            for (const auto& temp : TEMPLATES) {
                const std::vector<std::string> models = {"LG", "EGE"};
                for(const auto& model : models) {
                    const auto& param_list = (model == "LG") ? LG_PARAM_LIST : EGE_PARAM_LIST;
                    for (const auto& p : param_list) {
                        std::string splineName = chain + "_" + nucleusName + "_" + det + "_" + temp + "_" + model + "_" + p.first + "_spline";
                        TF1* spline = dynamic_cast<TF1*>(fin->Get(splineName.c_str()));
                        if (spline) {
                            cache[splineName] = std::unique_ptr<TF1>(static_cast<TF1*>(spline->Clone()));
                        }
                    }
                }
            }
        }
    }
    std::cout << "Finished pre-loading " << cache.size() << " splines into memory." << std::endl;
}

// 【新】用于修正参数的“安全”能量范围
const std::map<std::string, std::pair<double, double>> SAFE_DETECTOR_RANGES = {
    {"TOF", {0.33, 1.26}}, {"NaF", {0.75, 5.4}}, {"AGL", {2.95, 19.9}}
};

void getParamsFromSpline(const std::map<std::string, std::unique_ptr<TF1>>& splineCache, const std::string& chain, const std::string& elem, const std::string& det, 
                         const std::string& temp, const std::string& model, double ek, double* params, int z = 0) {
    const auto& param_list = (model == "LG") ? LG_PARAM_LIST : EGE_PARAM_LIST;
    
    const auto& det_range = DETECTOR_RANGES.at(det);
    bool is_extrapolating = (ek < det_range.first || ek > det_range.second);

    for (const auto& p : param_list) {
        const std::string& param_name = p.first;
        int param_index = p.second;
        std::string splineName = chain + "_" + elem + "_" + det + "_" + temp + "_" + model + "_" + param_name + "_spline";
        
        auto it = splineCache.find(splineName);
        if (it != splineCache.end()) {
            TF1* spline = it->second.get();
            double val = spline->Eval(ek);

            if (is_extrapolating || val <= 0) {
                double min_lim = 0, max_lim = 0;
                bool has_limits = false;

                if (param_name == "MPV" || param_name == "Peak") {
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

            if (val <= 0) {
                printf("  [!!!!!!!!!!!!!!!!!!!!!VALUE WARNING] Param '%s' is non-positive (%.4f) at Ek=%.3f for %s/%s/%s/%s.\n",
                       param_name.c_str(), val, ek, chain.c_str(), elem.c_str(), det.c_str(), temp.c_str());
            }

        } else {
            params[param_index] = 0;
        }
    }

    if (model == "LG") {
        params[2] = 1.0;
        params[4] = 1.0;

        // 检查 Width / Sigma 的比值
        double width = params[0]; // "Width" index 0
        double sigma = params[3]; // "Sigma" index 3
        
        if (sigma != 0 && (width / sigma > 20.0)) {
            printf("  [!!!!!!!!!!!RATIO WARNING] LG model Width/Sigma ratio > 20 (%.2f / %.2f = %.2f) at Ek=%.3f for %s/%s/%s/%s.\n",
                   width, sigma, width / sigma, ek, chain.c_str(), elem.c_str(), det.c_str(), temp.c_str());
        }

    } else if (model == "EGE") {
        params[5] = 1.0;
        params[6] = qmin_global;
        params[7] = qmax_global;
    }
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

    const std::string paramFileName = "/eos/user/z/zixuan/Isotope/ChargeFit/comparison_plots/allFitHistSplineSmooth_0.8_iter2.root";
    const std::string outFileName = "/eos/user/z/zixuan/Isotope/L2QTuning/CDFLookupTable_fromSpline_" + nucleusName + ".root";
    
    const std::string binningFileName = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFitParams_BeToOxy_0.8_iter2.root";
    const std::string binningHistName = "UnbiasedL1Inner_Beryllium_AGL_L1QTemplate_EGE_Peak";
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
    std::cout << "Starting to build CDF lookup tables for: " << nucleusName << " -> " << outFileName << std::endl;
    
    for (const auto& chain : CHAINS) {
        for (const auto& det : DETECTORS) {
            std::cout << "\nProcessing: " << chain << " / " << nucleusName << " / " << det << " ..." << std::endl;
            const size_t total_bins_in_loop = energyBins.size() - 1;

            for (size_t iy = 0; iy < total_bins_in_loop; ++iy) {
                double ekLow = energyBins[iy];
                double ekCenter = (energyBins[iy] + energyBins[iy+1]) / 2.0;

                if ((det == "TOF" && (ekCenter < 0.3 || ekCenter > 1.5)) ||
                    (det == "NaF" && (ekCenter < 0.6 || ekCenter > 6.0)) ||
                    (det == "AGL" && (ekCenter < 2.5 || ekCenter > 22.0))) continue;
                
                printf("--> Processing Bin %zu (Ek_low = %.4f, Ek_center = %.4f)\n", iy, ekLow, ekCenter);

                double par_lg_l1[5] = {0}, par_lg_l2[5] = {0};
                double par_ege_l1[8] = {0}, par_ege_l2[8] = {0};
                
                getParamsFromSpline(splineCache, chain, nucleusName, det, "L1QTemplate", "LG", ekCenter, par_lg_l1, z);
                getParamsFromSpline(splineCache, chain, nucleusName, det, "L2QTemplate", "LG", ekCenter, par_lg_l2, z);
                getParamsFromSpline(splineCache, chain, nucleusName, det, "L1QTemplate", "EGE", ekCenter, par_ege_l1, z);
                getParamsFromSpline(splineCache, chain, nucleusName, det, "L2QTemplate", "EGE", ekCenter, par_ege_l2, z);

                TF1 lg_l1("lg_l1", langaufun, qmin_global, qmax_global, 5); lg_l1.SetParameters(par_lg_l1); 
                TF1 lg_l2("lg_l2", langaufun, qmin_global, qmax_global, 5); lg_l2.SetParameters(par_lg_l2);
                TF1 ege_l1("ege_l1", funcExpGausExp, qmin_global, qmax_global, 8); ege_l1.SetParameters(par_ege_l1); 
                TF1 ege_l2("ege_l2", funcExpGausExp, qmin_global, qmax_global, 8); ege_l2.SetParameters(par_ege_l2);
                
                LookupTable lg_lookup, ege_lookup;
                std::string lg_context = Form("LG model for %s, %s, %s, bin %zu", chain.c_str(), nucleusName.c_str(), det.c_str(), iy);
                std::string ege_context = Form("EGE model for %s, %s, %s, bin %zu", chain.c_str(), nucleusName.c_str(), det.c_str(), iy);
                lg_lookup.build(lg_l1, lg_l2, qmin_global, qmax_global, lg_context);
                ege_lookup.build(ege_l1, ege_l2, qmin_global, qmax_global, ege_context);
                
                fout->cd();
                std::string baseName = Form("%s_%s_%s_bin%zu", chain.c_str(), nucleusName.c_str(), det.c_str(), iy);
                
                TVectorD(lg_lookup.q_values.size(), lg_lookup.q_values.data()).Write((baseName + "_LG_q").c_str());
                TVectorD(lg_lookup.cdf_l1.size(), lg_lookup.cdf_l1.data()).Write((baseName + "_LG_cdf_l1").c_str());
                TVectorD(lg_lookup.cdf_l2.size(), lg_lookup.cdf_l2.data()).Write((baseName + "_LG_cdf_l2").c_str());
                
                TVectorD(ege_lookup.q_values.size(), ege_lookup.q_values.data()).Write((baseName + "_EGE_q").c_str());
                TVectorD(ege_lookup.cdf_l1.size(), ege_lookup.cdf_l1.data()).Write((baseName + "_EGE_cdf_l1").c_str());
                TVectorD(ege_lookup.cdf_l2.size(), ege_lookup.cdf_l2.data()).Write((baseName + "_EGE_cdf_l2").c_str());

                std::clock_t now_clock = std::clock();
                time_t now_time = time(nullptr);
                double cpu_sec = double(now_clock - start_clock) / CLOCKS_PER_SEC;
                double wall_sec = difftime(now_time, start_time);
                double mem_mb = getCurrentRSS_MB();

                printf("  [Progress] Bin %3zu/%zu (%s/%s) done. | Total CPU: %8.2f s | Total Wall: %8.2f s | Mem: %8.2f MB\n",
                       iy + 1,
                       total_bins_in_loop,
                       chain.c_str(),
                       det.c_str(),
                       cpu_sec,
                       wall_sec,
                       mem_mb);
            }
        }
    }
    
    fout->Close();
    std::cout << "\nSuccessfully built and saved all CDF lookup tables to: " << outFileName << std::endl;
}