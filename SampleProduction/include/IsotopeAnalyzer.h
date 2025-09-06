#ifndef ISOTOPE_ANALYZER_H
#define ISOTOPE_ANALYZER_H

#include <string>
#include <vector>
#include <memory>
#include "basic_var.h"

// ROOT类的正向声明
class TChain;
class TH1;
class TH1F;
class TH2F;

namespace AMS_Iso {

// 其他类的正向声明
class HistManager; 

class IsotopeAnalyzer {
public:
    static IsotopeAnalyzer& Instance();

    // 配置和初始化
    void setConfig(const TString& outDir, const TString& outName,
                  const TString& inData, const TString& inOptions,
                  int UseMass = -1);
    void initialize();
    void write();
    void cleanup();

    // 访问器
    TChain* getDataChain() { return dataChain.get(); }
    const IsotopeVar* getIsotope() const { return isotope; }
    bool isISS() const { return isISS_; }
    int getCharge() const { return isotope ? isotope->getCharge() : -1; }
    int getUseMass() const { return UseMass_; }
    const std::vector<std::string>& getActiveChains() const { return active_chains_; }

    // --- 关键接口 ---
    // 为 selectdata 提供按名称获取任何直方图的通用方法
    TH1* getHist(const std::string& name) const;
    TH1F* getHist1F(const std::string& name) const;
    TH2F* getHist2F(const std::string& name) const;
    HistManager* getHistManager() const { return m_histManager.get(); }


private:
    IsotopeAnalyzer() = default;

    // 配置
    const IsotopeVar* isotope = nullptr;
    TString outDir_;
    TString outName_;
    TString inData_;
    TString inOptions_;
    int UseMass_ = -1;
    bool isISS_ = true;
    std::vector<std::string> active_chains_; // <--- 存储分析链

    // ROOT对象
    std::unique_ptr<TChain> dataChain;
    
    // 新的 HistManager 作为成员
    std::unique_ptr<HistManager> m_histManager;

    void readDataFrom(TChain* chain, const TString& filename);
    
    // 禁用拷贝和赋值
    IsotopeAnalyzer(const IsotopeAnalyzer&) = delete;
    IsotopeAnalyzer& operator=(const IsotopeAnalyzer&) = delete;
};

} // namespace AMS_Iso

#endif // ISOTOPE_ANALYZER_H

