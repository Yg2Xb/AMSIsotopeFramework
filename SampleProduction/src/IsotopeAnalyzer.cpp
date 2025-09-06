#include "IsotopeAnalyzer.h"
#include "ProductRegistry.h"
#include "BinningManager.h"
#include "HistManager.h"      // 确保包含了 HistManager 的完整定义
#include "basic_var.h"        // 确保 getIsotopeVar 可用
#include <iostream>
#include <stdexcept>
#include <fstream>

namespace AMS_Iso {

IsotopeAnalyzer& IsotopeAnalyzer::Instance() {
    static IsotopeAnalyzer instance;
    return instance;
}

void IsotopeAnalyzer::setConfig(const TString& outDir, const TString& outName,
                               const TString& inData, const TString& inOptions,
                               int UseMass) {
    outDir_ = outDir;
    outName_ = outName;
    inData_ = inData;
    inOptions_ = inOptions;

    int charge = std::stoi(inOptions.Data());
    if (charge > 0 && charge <= Constants::ELEMENT_COUNT) {
        isotope = &getIsotopeVar(charge);
    } else {
        throw std::runtime_error("Invalid charge specified in options: " + std::to_string(charge));
    }
    
    isISS_ = !inData_.Contains("MC");
    UseMass_ = (UseMass != -1) ? UseMass : (isotope ? isotope->getMass(0) : -1);

    std::cout << "Analyzer configured for: " << isotope->getName()
              << (isISS_ ? " (ISS Data)" : " (MC Data)")
              << ", Z=" << charge << ", A=" << UseMass_ << std::endl;
}

void IsotopeAnalyzer::readDataFrom(TChain* chain, const TString& filename) {
    if (filename.EndsWith(".txt")) {
        std::ifstream filelist(filename.Data());
        if (!filelist.is_open()) {
            throw std::runtime_error("Could not open file list: " + std::string(filename.Data()));
        }
        std::string line;
        while (std::getline(filelist, line)) {
            if (!line.empty() && line[0] != '#') { // 忽略空行和注释
                chain->Add(line.c_str());
            }
        }
    } else {
        chain->Add(filename);
    }
}

void IsotopeAnalyzer::initialize() {
    // 1. 初始化 Binning 管理器，加载所有预定义的分箱
    // MODIFIED: 调用无参数版本的 Initialize
    BinningManager::GetInstance().Initialize();

    // 2. 创建 HistManager
    TString output_filename = outDir_ + "/" + outName_ + ".root";
    m_histManager = std::make_unique<HistManager>(output_filename.Data());
    
    // 3. 将自身信息传递给注册表
    ProductRegistry::GetInstance().SetAnalyzerInfo(this);

    // 4. 定义分析链并注册所有产品蓝图
    active_chains_ = {"L1Inner", "Unbiased"};
    ProductRegistry::GetInstance().RegisterProducts(!isISS_, active_chains_);

    // 5. 让 HistManager 根据蓝图创建所有直方图
    m_histManager->CreateProducts(ProductRegistry::GetInstance().GetBlueprints());

    // 6. 设置数据链
    dataChain = std::make_unique<TChain>("amstreea");
    readDataFrom(dataChain.get(), inData_);

    std::cout << "IsotopeAnalyzer initialized with " << dataChain->GetEntries() << " entries." << std::endl;
}

void IsotopeAnalyzer::write() {
    if (m_histManager) {
        std::cout << "Saving all histograms to file..." << std::endl;
        m_histManager->Save();
    }
}

void IsotopeAnalyzer::cleanup() {
    // 目前没有需要清理的动态资源，保留为空
    std::cout << "Analysis finished." << std::endl;
}

TH1* IsotopeAnalyzer::getHist(const std::string& name) const {
    if (!m_histManager) {
        throw std::runtime_error("FATAL: HistManager is not initialized in IsotopeAnalyzer.");
    }
    TH1* hist = m_histManager->GetHistogram(name);
    if (!hist) {
        std::string error_msg = "FATAL: Histogram with name '" + name + "' was requested but not found.";
        throw std::runtime_error(error_msg);
    }
    return hist;
}

TH1F* IsotopeAnalyzer::getHist1F(const std::string& name) const {
    return dynamic_cast<TH1F*>(getHist(name));
}

TH2F* IsotopeAnalyzer::getHist2F(const std::string& name) const {
    return dynamic_cast<TH2F*>(getHist(name));
}

} // namespace AMS_Iso

