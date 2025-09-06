#ifndef BINNING_MANAGER_H
#define BINNING_MANAGER_H

#include <vector>
#include <string>
#include <map>
#include "basic_var.h" // 确保包含了 basic_var.h

namespace AMS_Iso {

class BinningManager {
public:
    static BinningManager& GetInstance();

    void Initialize();
    
    // 按名称获取通用分箱
    const std::vector<double>& Get(const std::string& name) const;
    // 按质量数获取Beta分箱 (源于KineticEnergyBins，主要用于曝光时间)
    const std::vector<double>& GetIsotopeBetaBins(int mass) const;

    // ADDED: 新增查询函数，方便通过(Z,A)获取专属的Ek/n分箱
    const std::vector<double>& GetEkPerNucleonBins(int charge, int mass) const;
    // ADDED: 新增查询函数，方便通过(Z,A)获取专属的Beta分箱 (源于标准刚度)
    const std::vector<double>& GetBetaBins(int charge, int mass) const;

    // 为 HistManager 提供别名，保持兼容性
    const std::vector<double>& GetBinning(const std::string& name) const {
        return Get(name);
    }

    // --- 静态工具函数 ---
    static std::vector<double> ConvertRigidityToEk(const std::vector<double>& rig_bins, int charge, int mass);
    static std::vector<double> ConvertRigidityToBeta(const std::vector<double>& rig_bins, int charge, int mass);
    static std::vector<double> ConvertEkToBeta(const std::vector<double>& ek_bins);

private:
    BinningManager() = default;
    
    BinningManager(const BinningManager&) = delete;
    BinningManager& operator=(const BinningManager&) = delete;
    
    std::map<std::string, std::vector<double>> m_bin_map; 
    std::map<int, std::vector<double>> m_isotope_beta_bins;
};

} // namespace AMS_Iso

#endif // BINNING_MANAGER_H

