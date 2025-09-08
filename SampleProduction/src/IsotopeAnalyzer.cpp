#include "IsotopeAnalyzer.h"
#include "BinningManager.h"
#include "HistManager.h"      // 使用我们新实现的 HistManager
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
    isISS_ = !inOptions.Contains("|MC");

    int charge = std::stoi(inOptions.Data());
    if (charge > 0 && charge <= Constants::ELEMENT_COUNT) {
        isotope = &getIsotopeVar(charge);
    } else {
        throw std::runtime_error("Invalid charge specified in options: " + std::to_string(charge));
    }

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
    BinningManager::GetInstance().Initialize();

    // 2. 定义分析链
    active_chains_ = {"UnbiasedL1Inner", "L1Inner"};

    // 3. 创建 HistManager （替代原来的 ProductRegistry）
    TString output_filename = outDir_ + "/" + outName_;
    m_histManager = std::make_unique<HistManager>(output_filename.Data(),
                                                  isISS_,
                                                  active_chains_,
                                                  isotope->getCharge(),
                                                  isotope,
                                                  UseMass_);

    // 4. 设置数据链
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
    std::cout << "Analysis finished." << std::endl;
}

TH1* IsotopeAnalyzer::getHist(const std::string& name) const {
    throw std::runtime_error("Direct getHist(name) is not supported anymore. "
                             "Please access histograms through HistManager arrays.");
}

TH1F* IsotopeAnalyzer::getHist1F(const std::string& name) const {
    return nullptr; // 不再支持
}

TH2F* IsotopeAnalyzer::getHist2F(const std::string& name) const {
    return nullptr; // 不再支持
}

} // namespace AMS_Iso