#include "IsoToolbox/BinningManager.h"
#include "IsoToolbox/Logger.h"
#include "IsoToolbox/PhysicsConstants.h"

namespace IsoToolbox {

BinningManager::BinningManager(AnalysisContext* context) : m_context(context) {}

void BinningManager::Initialize() {
    Logger::Info("Initializing BinningManager...");
    LoadRigidityBins();
    GenerateEkPerNucleonBins();
}

void BinningManager::LoadRigidityBins() {
    // 优先从YAML加载自定义刚度分档，如果不存在则使用C++中的默认值
    const auto& runSettings = m_context->GetRunSettings();
    if (runSettings["custom_rigidity_binning"]) {
        m_rigidity_bins = runSettings["custom_rigidity_binning"].as<std::vector<double>>();
        Logger::Info("Loaded custom rigidity binning from config file.");
    } else {
        m_rigidity_bins = PhysicsConstants::GetRigidityBins();
        Logger::Info("Using default rigidity binning from PhysicsConstants.");
    }
}

void BinningManager::GenerateEkPerNucleonBins() {
    const auto& particleInfo = m_context->GetParticleInfo();
    
    // 修正: 通过公共接口 GetCharge() 和 GetName() 访问数据
    int chargeZ = particleInfo.GetCharge();
    Logger::Debug("Generating Ek/n bins for nucleus: {}", particleInfo.GetName());

    // 修正: 遍历同位素质量信息需要通过 GetMasses() 和 GetIsotopeCount()
    const auto& masses = particleInfo.GetMasses();
    for (int i = 0; i < particleInfo.GetIsotopeCount(); ++i) {
        if (masses[i] == 0) continue; // 跳过无效的质量数
        
        int massA = masses[i];
        std::string isotope_name = particleInfo.GetName().substr(0, 2) + std::to_string(massA);
        
        Logger::Debug("... for isotope {}", isotope_name);

        m_isotope_binnings[isotope_name].ek_bins = BinningUtils::ConvertRigidityToEk(m_rigidity_bins, massA, chargeZ);
    }
}

const IsotopeBinning& BinningManager::GetIsotopeBinning(const std::string& isotope_name) const {
    auto it = m_isotope_binnings.find(isotope_name);
    if (it == m_isotope_binnings.end()) {
        throw std::runtime_error("Binning for isotope " + isotope_name + " not found.");
    }
    return it->second;
}

const std::vector<double>& BinningManager::GetRigidityBinning() const {
    return m_rigidity_bins;
}

// libs/IsoToolbox/src/BinningManager.cpp 添加实现
std::vector<double> BinningManager::GetBinning(const std::string& binning_type,
                                              const std::string& isotope_name,
                                              const std::string& detector) const {
    if (binning_type == "ek_per_nucleon") {
        if (!isotope_name.empty()) {
            const auto& isotope_binning = GetIsotopeBinning(isotope_name);
            return isotope_binning.ek_bins;
        } else {
            // 返回通用的narrow bins
            return {0.42, 0.51, 0.61, 0.73, 0.86, 1.00, 1.17, 1.35, 1.55, 1.77, 2.01, 2.28, 2.57, 2.88, 3.23, 3.60, 4.00, 4.44, 4.91, 5.42, 5.99, 6.56, 7.18, 7.86, 8.60, 9.40, 10.25, 11.16, 12.13, 14.35, 16.38};
        }
    } else if (binning_type == "rigidity") {
        return GetRigidityBinning();
    }
    
    return {};
}

std::vector<double> BinningManager::GetChargeBinning(const std::string& source_type) const {
    // 根据源粒子类型返回合适的电荷分bin
    if (source_type == "Carbon") {
        return {3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0};
    } else if (source_type == "Nitrogen") {
        return {4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0};
    } else if (source_type == "Oxygen") {
        return {5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0};
    }
    
    // 默认通用电荷分bin
    return {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0};
}

} // namespace IsoToolbox