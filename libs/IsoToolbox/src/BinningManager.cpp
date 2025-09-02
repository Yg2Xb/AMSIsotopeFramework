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
    
    int chargeZ = particleInfo.GetCharge();
    Logger::Debug("Generating Ek/n bins for nucleus: {}", particleInfo.GetName());

    // 为每个同位素生成专属的Ek/n分bin
    const auto& masses = particleInfo.GetMasses();
    for (int i = 0; i < particleInfo.GetIsotopeCount(); ++i) {
        if (masses[i] == 0) continue; // 跳过无效的质量数
        
        int massA = masses[i];
        std::string isotope_name = particleInfo.GetName().substr(0, 2) + std::to_string(massA);
        
        Logger::Debug("... for isotope {}", isotope_name);

        m_isotope_binnings[isotope_name].ek_bins = BinningUtils::ConvertRigidityToEk(m_rigidity_bins, massA, chargeZ);
    }
}

std::vector<double> BinningManager::GenerateDefaultNarrowEkBins() const {
    // 科研惯例：使用4/7质荷比生成默认的窄带Ek/n分bin
    const int massA = 7;
    const int chargeZ = 4;
    
    return BinningUtils::ConvertRigidityToEk(m_rigidity_bins, massA, chargeZ);
}

std::vector<double> BinningManager::GenerateChargeBins(int central_z) const {
    // 生成Z±2范围，400个bin的电荷分bin
    const double z_min = central_z - 2.0;
    const double z_max = central_z + 2.0;
    const int n_bins = 400;
    
    std::vector<double> charge_bins;
    for (int i = 0; i <= n_bins; ++i) {
        double z = z_min + (z_max - z_min) * i / n_bins;
        charge_bins.push_back(z);
    }
    
    return charge_bins;
}

std::vector<double> BinningManager::GenerateInvMassBins() const {
    // 1/M分bin：0.0到0.5，200个bin
    const double inv_mass_min = 0.0;
    const double inv_mass_max = 0.5;
    const int n_bins = 200;
    
    std::vector<double> inv_mass_bins;
    for (int i = 0; i <= n_bins; ++i) {
        double inv_mass = inv_mass_min + (inv_mass_max - inv_mass_min) * i / n_bins;
        inv_mass_bins.push_back(inv_mass);
    }
    
    return inv_mass_bins;
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

std::vector<double> BinningManager::GetBinning(const std::string& binning_type,
                                              const std::string& isotope_name,
                                              const std::string& detector) const {
    if (binning_type == "ek_per_nucleon") {
        if (!isotope_name.empty()) {
            // 使用特定同位素的分bin
            const auto& isotope_binning = GetIsotopeBinning(isotope_name);
            return isotope_binning.ek_bins;
        } else {
            // 使用默认的窄带分bin（4/7质荷比）
            if (m_default_narrow_ek_bins.empty()) {
                m_default_narrow_ek_bins = GenerateDefaultNarrowEkBins();
            }
            return m_default_narrow_ek_bins;
        }
    } else if (binning_type == "rigidity") {
        return GetRigidityBinning();
    } else if (binning_type == "invmass") {
        // 1/M分bin：0.0到0.5，200个bin
        static std::vector<double> inv_mass_bins = GenerateInvMassBins();
        return inv_mass_bins;
    }
    
    return {};
}

std::vector<double> BinningManager::GetChargeBinning(const std::string& source_type) const {
    // 检查缓存
    auto it = m_charge_binnings.find(source_type);
    if (it != m_charge_binnings.end()) {
        return it->second;
    }
    
    // 智能查找：利用PhysicsConstants获取电荷数
    int central_z = 6; // 默认值
    try {
        const auto& isotope = PhysicsConstants::GetIsotope(source_type);
        central_z = isotope.GetCharge();
    } catch (const std::exception& e) {
        Logger::Warn("Source type '{}' not found in PhysicsConstants, using default Z=6", source_type);
    }
    
    // 生成并缓存电荷分bin
    m_charge_binnings[source_type] = GenerateChargeBins(central_z);
    return m_charge_binnings[source_type];
}

} // namespace IsoToolbox