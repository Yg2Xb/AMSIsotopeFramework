#include "BinningManager.h"
#include "Tool.h"
#include "basic_var.h" 
#include <stdexcept>
#include <iostream>

namespace AMS_Iso {

BinningManager& BinningManager::GetInstance() {
    static BinningManager instance;
    return instance;
}

void BinningManager::Initialize() {
    if (!m_bin_map.empty()) return;

    m_bin_map.clear();
    m_isotope_beta_bins.clear();

    std::cout << "BinningManager: Pre-calculating bins for all isotopes (Z=1 to 8)..." << std::endl;

    const auto& rig_bins_arr = Binning::RigidityBins;
    std::vector<double> rigidity_bins_vec(rig_bins_arr.begin(), rig_bins_arr.end());
    m_bin_map["Rigidity"] = rigidity_bins_vec;

    const auto& ek_wide_arr = Binning::EkWideBin;
    std::vector<double> ek_wide_vec(ek_wide_arr.begin(), ek_wide_arr.end());
    m_bin_map["BetaFromEkWide"] = ConvertEkToBeta(ek_wide_vec);
    
    m_bin_map["EkPerNucleon"] = ConvertRigidityToEk(rigidity_bins_vec, 2, 4);
    m_bin_map["Beta"] = ConvertRigidityToBeta(rigidity_bins_vec, 2, 4);

    m_bin_map["InverseMass"] = {}; 
    m_bin_map["DeltaBeta"] = {};   
    m_bin_map["BetaGamma"] = {};   
    m_bin_map["Charge"] = {};      
    m_bin_map["EkGen"] = {};       

    for (int charge = 1; charge <= Constants::ELEMENT_COUNT; ++charge) {
        const auto& isotope_info = getIsotopeVar(charge);
        for (int i = 0; i < isotope_info.getIsotopeCount(); ++i) {
            int current_mass = isotope_info.getMass(i);
            std::string isotope_name = isotope_info.getName() + std::to_string(current_mass);
            
            // --- 计算并缓存该同位素专属的 EkPerNucleon bin ---
            std::string ek_bin_key = "EkPerNucleon_" + isotope_name;
            if (m_bin_map.find(ek_bin_key) == m_bin_map.end()) {
                 m_bin_map[ek_bin_key] = ConvertRigidityToEk(rigidity_bins_vec, charge, current_mass);
            }
            
            // --- 计算并缓存该同位素专属的 Beta bin (源于标准刚度) ---
            std::string beta_bin_key = "Beta_" + isotope_name;
            if (m_bin_map.find(beta_bin_key) == m_bin_map.end()) {
                 m_bin_map[beta_bin_key] = ConvertRigidityToBeta(rigidity_bins_vec, charge, current_mass);
            }

            // --- 计算并缓存该同位素专属的 Beta bin (源于KineticEnergyBins, 用于曝光时间) ---
            if (m_isotope_beta_bins.find(current_mass) == m_isotope_beta_bins.end()) {
                const auto& ek_bins_for_isotope_arr = Binning::KineticEnergyBins[charge - 1][i];
                std::vector<double> ek_bins_for_isotope_vec;
                ek_bins_for_isotope_vec.reserve(Constants::RIGIDITY_BINS + 1);
                for(int j = 0; j <= Constants::RIGIDITY_BINS; ++j) {
                    ek_bins_for_isotope_vec.push_back(ek_bins_for_isotope_arr[j]);
                }
                m_isotope_beta_bins[current_mass] = ConvertEkToBeta(ek_bins_for_isotope_vec);
            }
        }
    }
    
    std::cout << "BinningManager initialized with " << m_bin_map.size() << " general/isotope-specific bins and "
              << m_isotope_beta_bins.size() << " total isotope-specific beta binnings." << std::endl;
}

const std::vector<double>& BinningManager::Get(const std::string& name) const {
    auto it = m_bin_map.find(name);
    if (it != m_bin_map.end()) {
        return it->second;
    }
    throw std::runtime_error("Binning with name '" + name + "' not found in BinningManager.");
}

const std::vector<double>& BinningManager::GetIsotopeBetaBins(int mass) const {
    auto it = m_isotope_beta_bins.find(mass);
    if (it != m_isotope_beta_bins.end()) {
        return it->second;
    }
    throw std::runtime_error("Beta binning for mass " + std::to_string(mass) + " not found in BinningManager cache.");
}

const std::vector<double>& BinningManager::GetEkPerNucleonBins(int charge, int mass) const {
    const auto& isotope_info = getIsotopeVar(charge);
    std::string isotope_name = isotope_info.getName() + std::to_string(mass);
    std::string ek_bin_key = "EkPerNucleon_" + isotope_name;
    return Get(ek_bin_key);
}

const std::vector<double>& BinningManager::GetBetaBins(int charge, int mass) const {
    const auto& isotope_info = getIsotopeVar(charge);
    std::string isotope_name = isotope_info.getName() + std::to_string(mass);
    std::string beta_bin_key = "Beta_" + isotope_name;
    return Get(beta_bin_key);
}


// --- 静态工具函数实现 ---
std::vector<double> BinningManager::ConvertRigidityToEk(const std::vector<double>& rig_bins, int charge, int mass) {
    if (charge == 0 || mass == 0) throw std::runtime_error("Charge and mass cannot be zero for Ek conversion.");
    std::vector<double> ek_bins;
    ek_bins.reserve(rig_bins.size());
    for (double R : rig_bins) {
        ek_bins.push_back(Tools::rigidityToKineticEnergy(R, charge, mass));
    }
    return ek_bins;
}

std::vector<double> BinningManager::ConvertRigidityToBeta(const std::vector<double>& rig_bins, int charge, int mass) {
    if (charge == 0 || mass == 0) throw std::runtime_error("Charge and mass cannot be zero for Beta conversion.");
    std::vector<double> beta_bins;
    beta_bins.reserve(rig_bins.size());
    for (double R : rig_bins) {
        beta_bins.push_back(Tools::rigidityToBeta(R, charge, mass, false));
    }
    return beta_bins;
}

std::vector<double> BinningManager::ConvertEkToBeta(const std::vector<double>& ek_bins) {
    std::vector<double> beta_bins;
    beta_bins.reserve(ek_bins.size());
    for (double ek : ek_bins) {
        beta_bins.push_back(Tools::kineticEnergyToBeta(ek));
    }
    return beta_bins;
}

} // namespace AMS_Iso

