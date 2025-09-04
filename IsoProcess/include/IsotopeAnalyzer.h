/***********************************************************
 *  File: IsotopeAnalyzer.h
 *
 *  Header file for AMS Isotopes analysis framework.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/

#pragma once

// 1. 首先包含标准库 string_view
#include <string_view>

#include <vector>
#include <memory>
#include <map>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include "basic_var.h"
#include "HistogramManager.h"

namespace AMS_Iso {

class IsotopeAnalyzer {
public:
    static IsotopeAnalyzer& Instance() {
        static IsotopeAnalyzer instance;
        return instance;
    }

    // 配置和初始化
    void setConfig(const TString& outDir, const TString& outName,
                  const TString& inData, const TString& inOptions,
                  int UseMass = -1);
    void initialize();
    void write();
    void cleanup();

    // 访问器
    TFile* getOutputFile() { return fOutRoot.get(); }
    TChain* getDataChain() { return dataChain.get(); }
    TTree* getSaveTree() { return saveTree.get(); }
    const IsotopeVar* getIsotope() const { return isotope; }
    bool isISS() const { return isISS_; }
    // 添加这两个方法
    int getCharge() const { return isotope ? isotope->getCharge() : -1; }
    int getUseMass() const { return UseMass_; }
    const TString& getInOptions() const { return inOptions_; }  // 如果还需要访问原始选项
    // 直方图访问代理方法
    TH1* getHistogram(HistType type, const std::string& name) const {
        return histManager_->getHistogram(type, name);
    }
    TH1D* getHist1D(HistType type, const std::string& name) const {
        return histManager_->getHist1D(type, name);
    }
    TH2D* getHist2D(HistType type, const std::string& name) const {
        return histManager_->getHist2D(type, name);
    }
    TH1D* getBetaExposureHist(int isotopeIdx, int betaTypeIdx) const {
        return histManager_->getBetaExposureHist(isotopeIdx, betaTypeIdx);
    }
    TH1D* getEventHist(DetType det, int cutIdx) const {
        return histManager_->getEventHist(det, cutIdx);
    }

private:
    IsotopeAnalyzer() = default;
    void setupBranches();

    // 配置
    const IsotopeVar* isotope = nullptr;
    TString outDir_;
    TString outName_;
    TString inData_;
    TString inOptions_;
    int UseMass_ = -1;
    bool isISS_ = true;

    // ROOT对象
    std::unique_ptr<TFile> fOutRoot;
    std::unique_ptr<TChain> dataChain;
    std::unique_ptr<TTree> saveTree;

    void readDataFrom(TChain* chain, const TString& filename);
    
    // 直方图管理器
    std::unique_ptr<HistogramManager> histManager_;

    // 禁用拷贝和赋值
    IsotopeAnalyzer(const IsotopeAnalyzer&) = delete;
    IsotopeAnalyzer& operator=(const IsotopeAnalyzer&) = delete;

};

} // namespace AMS_Iso