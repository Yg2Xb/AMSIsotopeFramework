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
    
    m_bin_map["EkPerNucleon"] = ConvertRigidityToEk(rigidity_bins_vec, 2, 4);
    m_bin_map["Beta"] = ConvertRigidityToBeta(rigidity_bins_vec, 2, 4);
    m_bin_map["BetaGamma"] = ConvertRigidityToBetaGamma(rigidity_bins_vec, 2, 4);

    m_bin_map["InverseMass"] = {}; 
    m_bin_map["1/NaFBeta"] = {};   
    m_bin_map["1/AGLBeta"] = {};   
    m_bin_map["DeltaNaFBeta"] = {};   
    m_bin_map["DeltaAGLBeta"] = {};   
    m_bin_map["DeltaTOFBeta"] = {};   
    m_bin_map["Charge"] = {};      

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
            std::string betagamma_bin_key = "BetaGamma_" + isotope_name;
            if (m_bin_map.find(betagamma_bin_key) == m_bin_map.end()) {
                 m_bin_map[betagamma_bin_key] = ConvertRigidityToBetaGamma(rigidity_bins_vec, charge, current_mass);
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

const std::vector<double>& BinningManager::GetBetaGammaBins(int charge, int mass) const {
    const auto& isotope_info = getIsotopeVar(charge);
    std::string isotope_name = isotope_info.getName() + std::to_string(mass);
    std::string betagamma_bin_key = "BetaGamma_" + isotope_name;
    return Get(betagamma_bin_key);
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

std::vector<double> BinningManager::ConvertRigidityToBetaGamma(const std::vector<double>& rig_bins, int charge, int mass) {
    if (charge == 0 || mass == 0) throw std::runtime_error("Charge and mass cannot be zero for Beta conversion.");
    std::vector<double> betagamma_bins;
    betagamma_bins.reserve(rig_bins.size());
    for (double R : rig_bins) {
        double converted_beta = Tools::rigidityToBeta(R, charge, mass, false);
        double beta_gamma = converted_beta / std::sqrt(1 - converted_beta * converted_beta);
        betagamma_bins.push_back(beta_gamma); 
    }
    return betagamma_bins;
}


} // namespace AMS_Iso

