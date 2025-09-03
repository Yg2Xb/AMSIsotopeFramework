#include "IsoToolbox/AnalysisContext.h"
#include "IsoToolbox/Logger.h"
#include <stdexcept>

namespace IsoToolbox {

AnalysisContext::AnalysisContext(const std::string& config_path) : m_particleInfo(nullptr) {
    try {
        m_configNode = YAML::LoadFile(config_path);
    } catch (const YAML::Exception& e) {
        throw std::runtime_error("Failed to load or parse config file '" + config_path + "': " + e.what());
    }
    parseConfig();
}

void AnalysisContext::parseConfig() {
    std::string particle_name = m_configNode["run_settings"]["target_nucleus"].as<std::string>();
    loadParticleData(particle_name);

    const auto& chain_node = m_configNode["analysis_chain"];
    m_analysisChain.chain_id           = chain_node["chain_id"].as<std::string>();
    m_analysisChain.rigidity_version = chain_node["rigidity_version"].as<std::string>();
    m_analysisChain.beta_version     = chain_node["beta_version"].as<std::string>(); // Changed from velocity_version
    m_analysisChain.cut_version      = chain_node["cut_version"].as<std::string>();

    const auto& samples_nodes = m_configNode["samples_to_process"];
    for (const auto& node : samples_nodes) {
        m_samples.push_back({
            node["name"].as<std::string>(),
            node["type"].as<std::string>()
        });
    }
}

void AnalysisContext::loadParticleData(const std::string& particleName) {
    m_particleInfo = std::make_unique<Isotope>(PhysicsConstants::GetIsotope(particleName));
}

const Isotope& AnalysisContext::GetParticleInfo() const {
    if (!m_particleInfo) {
        throw std::runtime_error("ParticleInfo has not been initialized.");
    }
    return *m_particleInfo;
}

const AnalysisChain& AnalysisContext::GetAnalysisChain() const {
    return m_analysisChain;
}

const std::vector<SampleInfo>& AnalysisContext::GetSamplesToProcess() const {
    return m_samples;
}

const YAML::Node AnalysisContext::GetRunSettings() const {
    return m_configNode;
}

} // namespace IsoToolbox
