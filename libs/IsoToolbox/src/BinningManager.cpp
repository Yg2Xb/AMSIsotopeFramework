#include "IsoToolbox/BinningManager.h"
#include "IsoToolbox/AnalysisContext.h"
#include "IsoToolbox/PhysicsConstants.h"
#include "IsoToolbox/Tools.h"
#include "IsoToolbox/Logger.h"
#include <stdexcept>

BinningManager::BinningManager(const AnalysisContext* context) : m_context(context) {
    if (!m_context) {
        throw std::runtime_error("BinningManager requires a valid AnalysisContext upon construction.");
    }
}

void BinningManager::Initialize() {
    LOG_INFO("Initializing BinningManager...");
    GenerateEkPerNucleonBins();
    LOG_INFO("BinningManager initialized successfully.");
}

void BinningManager::GenerateEkPerNucleonBins() {
    const auto& rigidityBins = Constants::RigidityBinsVec;
    if (rigidityBins.empty()) {
        LOG_WARNING("Default Rigidity Bins vector is empty. Cannot generate Ek/n bins.");
        return;
    }

    const auto& isotopes = m_context->GetIsotopes();
    int chargeZ = m_context->GetChargeZ();

    LOG_DEBUG("Generating Ek/n bins for nucleus: %s", m_context->GetTargetNucleus().c_str());

    for (const auto& isotope : isotopes) {
        std::vector<double> ekBins;
        ekBins.reserve(rigidityBins.size());

        LOG_DEBUG("  - Isotope: %s (Z=%d, A=%d)", isotope.name.c_str(), chargeZ, isotope.mass);

        for (double rigidity : rigidityBins) {
            ekBins.push_back(Tools::RigidityToEkPerNucleon(rigidity, chargeZ, isotope.mass));
        }
        m_ekPerNucleonBins[isotope.name] = ekBins;
    }
}

const std::vector<double>& BinningManager::GetEkPerNucleonBins(const std::string& isotopeName) const {
    auto it = m_ekPerNucleonBins.find(isotopeName);
    if (it != m_ekPerNucleonBins.end()) {
        return it->second;
    }
    LOG_ERROR("Could not find Ek/n bins for isotope: %s. Returning empty vector.", isotopeName.c_str());
    return m_emptyBins;
}