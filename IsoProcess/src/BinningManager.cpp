// File: IsoProcess/src/BinningManager.cpp
// Purpose: Corrected version that calls the namespace functions from Tool.h.

#include "BinningManager.h"
#include "Tool.h"           // <-- Includes your namespace AMS_Iso::Tools
#include "basic_var.h"      // For AMS_Iso::Binning::RigidityBins etc.
#include <stdexcept>
#include <iostream>

// Constructor is now simpler
BinningManager::BinningManager() {}

void BinningManager::Initialize(int charge, int mass) {
    m_chargeZ = charge;
    m_massA = mass;

    const auto& rig_bins_arr = AMS_Iso::Binning::RigidityBins;
    std::vector<double> rigidity_bins_vec(rig_bins_arr.begin(), rig_bins_arr.end());

    m_bin_map["EkPerNucleon"] = ConvertRigidityToEk(rigidity_bins_vec);
    m_bin_map["Beta"] = ConvertRigidityToBeta(rigidity_bins_vec);
    m_bin_map["Beta_vs_Rigidity"] = m_bin_map["Beta"]; // For compatibility
}

const std::vector<double>& BinningManager::Get(const std::string& name) const {
    auto it = m_bin_map.find(name);
    if (it != m_bin_map.end()) {
        return it->second;
    }

    if (name == "Rigidity" || name == "RigidityBins") {
        const auto& bins = AMS_Iso::Binning::RigidityBins;
        m_bin_map[name] = std::vector<double>(bins.begin(), bins.end());
        return m_bin_map.at(name);
    }
    if (name == "EkWideBin") {
        const auto& bins = AMS_Iso::Binning::EkWideBin;
        m_bin_map[name] = std::vector<double>(bins.begin(), bins.end());
        return m_bin_map.at(name);
    }

    throw std::runtime_error("Binning with name '" + name + "' not found.");
}

std::vector<double> BinningManager::ConvertRigidityToEk(const std::vector<double>& rig_bins) {
    if (m_chargeZ == 0 || m_massA == 0) {
        throw std::runtime_error("BinningManager must be initialized (Ek).");
    }
    std::vector<double> ek_bins;
    ek_bins.reserve(rig_bins.size());
    for (double R : rig_bins) {
        // *** CORRECTED: Calling the namespace function directly ***
        ek_bins.push_back(AMS_Iso::Tools::rigidityToKineticEnergy(R, m_chargeZ, m_massA));
    }
    return ek_bins;
}

std::vector<double> BinningManager::ConvertRigidityToBeta(const std::vector<double>& rig_bins) {
    if (m_chargeZ == 0 || m_massA == 0) {
        throw std::runtime_error("BinningManager must be initialized (Beta).");
    }
    std::vector<double> beta_bins;
    beta_bins.reserve(rig_bins.size());
    for (double R : rig_bins) {
        // *** CORRECTED: Calling the namespace function directly ***
        beta_bins.push_back(AMS_Iso::Tools::rigidityToBeta(R, m_chargeZ, m_massA));
    }
    return beta_bins;
}