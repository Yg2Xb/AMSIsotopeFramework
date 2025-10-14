// calculate_all_integrals.C

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>

// ROOT headers
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"

// 包含工具
#include "../Tool.h"

using namespace AMS_Iso;

// 核素信息：名称 -> (Z, A)
std::map<std::string, std::pair<int, int>> getNucleiInfo() {
    return {
        {"Be7",  {4, 7}},
        {"Be9",  {4, 9}},
        {"Be10", {4, 10}},
        {"B10",  {5, 10}},
        {"B11",  {5, 11}},
        {"C12",  {6, 12}},
        {"N14",  {7, 14}},
        {"N15",  {7, 15}},
        {"O16",  {8, 16}}
    };
}

void cal_integral() {
    // --- 输入输出文件 ---
    const char* be10_fileName = "/eos/user/z/zixuan/Isotope/Add/Be10_all.root";
    const char* tf1_inFileName = "/eos/user/z/zixuan/Isotope/FluxSmooth/FluxSmooth.root";
    const char* outFileName = "/eos/user/z/zixuan/Isotope/FluxSmooth/tf1_integrals.root";

    // --- 获取MC_FLUX_H3的动能bin结构 ---
    TFile* be10_file = TFile::Open(be10_fileName);
    if (!be10_file || be10_file->IsZombie()) {
        std::cerr << "[FATAL] Cannot open Be10 file: " << be10_fileName << std::endl;
        return;
    }
    
    TH1* h_flux_h3 = (TH1*)be10_file->Get("MC_FLUX_H3");
    if (!h_flux_h3) {
        std::cerr << "[FATAL] Cannot find MC_FLUX_H3" << std::endl;
        be10_file->Close();
        return;
    }
    
    // 提取动能bin边界
    int nEkBins = h_flux_h3->GetNbinsX();
    std::vector<double> ekBins(nEkBins + 1);
    for (int i = 0; i <= nEkBins; i++) {
        ekBins[i] = h_flux_h3->GetBinLowEdge(i + 1);
    }
    ekBins[nEkBins] = h_flux_h3->GetBinLowEdge(nEkBins) + h_flux_h3->GetBinWidth(nEkBins);
    
    std::cout << "Extracted " << nEkBins << " Ek bins from MC_FLUX_H3" << std::endl;
    std::cout << "Ek range: " << ekBins[0] << " - " << ekBins[nEkBins] << " GeV/n" << std::endl;
    
    be10_file->Close();
    
    // --- 打开输出文件 ---
    TFile* outFile = new TFile(outFileName, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "[FATAL] Cannot create output file: " << outFileName << std::endl;
        return;
    }
    
    // =======================================================================
    //  第一部分: 计算解析函数的积分 (R^-1) - 用Be10 (Z=4, A=10)
    // =======================================================================
    std::cout << "\n--- Part 1: Calculating R^-1 with Be10 (Z=4, A=10) ---" << std::endl;
    
    // 将Be10的Ek bin转为刚度bin
    std::vector<double> rigBins_Be10(nEkBins + 1);
    for (int i = 0; i <= nEkBins; i++) {
        rigBins_Be10[i] = kineticEnergyToRigidity(ekBins[i], 4, 10); // Be10: Z=4, A=10
    }
    
    std::cout << "Be10 rigidity range: " << rigBins_Be10[0] << " - " << rigBins_Be10[nEkBins] << " GV" << std::endl;
    
    TH1D* h_integral_R_inv = new TH1D("h_integral_R_inv", "Integral of R^{-1}", nEkBins, ekBins.data());
    TH1D* h_ratio_R_inv = new TH1D("h_ratio_R_inv", "Ratio of R^{-1}", nEkBins, ekBins.data());
    
    const double r_min = rigBins_Be10[0];
    const double r_max = rigBins_Be10[nEkBins];
    const double total_inv = std::log(r_max / r_min);
    
    for (int i = 0; i < nEkBins; i++) {
        double integral = std::log(rigBins_Be10[i+1] / rigBins_Be10[i]);
        h_integral_R_inv->SetBinContent(i + 1, integral);
        h_ratio_R_inv->SetBinContent(i + 1, integral / total_inv);
    }
    
    outFile->cd();
    h_integral_R_inv->Write();
    h_ratio_R_inv->Write();
    
    delete h_integral_R_inv;
    delete h_ratio_R_inv;
    
    std::cout << "  -> R^-1 calculation completed with total integral = " << total_inv << std::endl;
    
    // =======================================================================
    //  第二部分: 处理各个核素的TF1
    // =======================================================================
    std::cout << "\n--- Part 2: Calculating for each nucleus TF1 ---" << std::endl;
    
    auto nucleiInfo = getNucleiInfo();
    
    TFile* inFile = TFile::Open(tf1_inFileName);
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "[FATAL] Cannot open input file: " << tf1_inFileName << std::endl;
        outFile->Close();
        return;
    }
    
    for (const auto& [name, info] : nucleiInfo) {
        std::string tf1_name = name + "_R";
        TF1* func = inFile->Get<TF1>(tf1_name.c_str());
        if (!func) {
            std::cerr << "[WARNING] TF1 '" << tf1_name << "' not found. Skipping." << std::endl;
            continue;
        }
        
        int Z = info.first;
        int A = info.second;
        std::cout << "Processing " << name << " (Z=" << Z << ", A=" << A << ")..." << std::endl;
        
        // 为该核素转换动能bin为刚度bin
        std::vector<double> rigBins_nucleus(nEkBins + 1);
        for (int i = 0; i <= nEkBins; i++) {
            rigBins_nucleus[i] = kineticEnergyToRigidity(ekBins[i], Z, A);
        }
        
        // 创建直方图（使用相同的Ek bin结构）
        TH1D* h_integral = new TH1D(Form("h_integral_%s", tf1_name.c_str()), 
                                   Form("Integral of %s", tf1_name.c_str()), 
                                   nEkBins, ekBins.data());
        TH1D* h_ratio = new TH1D(Form("h_ratio_%s", tf1_name.c_str()), 
                                 Form("Ratio of %s", tf1_name.c_str()), 
                                 nEkBins, ekBins.data());
        
        // 计算总积分（在该核素的刚度空间）
        double totalIntegral = func->Integral(rigBins_nucleus[0], rigBins_nucleus[nEkBins]);
        if (totalIntegral == 0) {
            std::cerr << "[WARNING] Total integral for '" << tf1_name << "' is zero." << std::endl;
        }
        
        // 逐bin计算积分
        for (int i = 0; i < nEkBins; i++) {
            // 在该核素的刚度空间积分
            double binIntegral = func->Integral(rigBins_nucleus[i], rigBins_nucleus[i+1]);
            
            // 存储到Ek bin中
            h_integral->SetBinContent(i + 1, binIntegral);
            h_ratio->SetBinContent(i + 1, (totalIntegral != 0) ? (binIntegral / totalIntegral) : 0.0);
        }
        
        // 写入文件
        outFile->cd();
        h_integral->Write();
        h_ratio->Write();
        
        delete h_integral;
        delete h_ratio;
        
        std::cout << "  -> " << name << " completed. Total integral = " << totalIntegral << std::endl;
        std::cout << "     Rigidity range: " << rigBins_nucleus[0] << " - " << rigBins_nucleus[nEkBins] << " GV" << std::endl;
    }
    
    // --- 清理 ---
    inFile->Close();
    outFile->Close();
    
    std::cout << "\n>>> All tasks finished. Results saved to " << outFileName << std::endl;
}