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
    // 新增：智能获取分bin
    std::vector<double> GetBinning(const std::string& binning_type,
                                  const std::string& isotope_name = "",
                                  const std::string& detector = "") const;
                                  
    // 新增：获取电荷分bin（用于背景分析）
    std::vector<double> GetChargeBinning(const std::string& source_type) const;

private:
    void LoadRigidityBins();
    void GenerateEkPerNucleonBins();

    AnalysisContext* m_context;
    std::vector<double> m_rigidity_bins;
    // 使用map来存储每种同位素的专属分bin方案
    std::map<std::string, IsotopeBinning> m_isotope_binnings;
    // 新增：不同源粒子的电荷分bin缓存
    std::map<std::string, std::vector<double>> m_charge_binnings;
    void LoadChargeBinnings();
};

} // namespace IsoToolbox