#pragma once

#include "IsoToolbox/AnalysisContext.h"
#include <vector>
#include <string>
#include <map>

namespace IsoToolbox {

// 存储单个同位素的分档方案
struct IsotopeBinning {
    std::vector<double> ek_bins;
    std::vector<double> velocity_bins;
};

class BinningManager {
public:
    explicit BinningManager(AnalysisContext* context);

    void Initialize();

    const IsotopeBinning& GetIsotopeBinning(const std::string& isotope_name) const;
    const std::vector<double>& GetRigidityBinning() const;
    
    // 智能获取分bin
    std::vector<double> GetBinning(const std::string& binning_type,
                                  const std::string& isotope_name = "",
                                  const std::string& detector = "") const;
                                  
    // 获取电荷分bin（用于背景分析）
    std::vector<double> GetChargeBinning(const std::string& source_type) const;

private:
    void LoadRigidityBins();
    void GenerateEkPerNucleonBins();
    
    // 生成默认的窄带Ek/n分bin（基于4/7质荷比）
    std::vector<double> GenerateDefaultNarrowEkBins() const;
    
    // 生成电荷分bin（Z±2，400个bin）
    std::vector<double> GenerateChargeBins(int central_z) const;
    
    // 生成1/M分bin（0.0-0.5，200个bin）
    std::vector<double> GenerateInvMassBins() const;

    AnalysisContext* m_context;
    std::vector<double> m_rigidity_bins;
    // 使用map来存储每种同位素的专属分bin方案
    std::map<std::string, IsotopeBinning> m_isotope_binnings;
    
    // 缓存默认的窄带Ek/n分bin
    mutable std::vector<double> m_default_narrow_ek_bins;
    
    // 缓存不同源粒子的电荷分bin
    mutable std::map<std::string, std::vector<double>> m_charge_binnings;
};

} // namespace IsoToolbox