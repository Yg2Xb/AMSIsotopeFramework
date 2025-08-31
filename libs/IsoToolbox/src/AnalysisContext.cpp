#include <IsoToolbox/AnalysisContext.h>
#include <stdexcept>

namespace IsoToolbox {

// BUG FIX: Implementation for the new convenience constructor.
AnalysisContext::AnalysisContext(const std::string& config_path) {
    // This smart pointer uses a no-op deleter because the ConfigManager is a singleton
    // and its lifetime is managed by its static instance, not by this shared_ptr.
    m_configManager = std::shared_ptr<ConfigManager>(&ConfigManager::GetInstance(), [](auto*){});
    if (!m_configManager->Load(config_path)) {
        throw std::runtime_error("AnalysisContext: Failed to load config file: " + config_path);
    }
    parseConfig();
}

AnalysisContext::AnalysisContext(std::shared_ptr<ConfigManager> config) : m_configManager(config) {
    parseConfig();
}

void AnalysisContext::parseConfig() {
    // Parse particle info
    std::string particle_name = m_configManager->Get<std::string>("current_particle");
    m_particleInfo.name = particle_name;
    m_particleInfo.charge = m_configManager->Get<int>("particle_definitions." + particle_name + ".charge");
    auto isotopes_nodes = m_configManager->GetNode("particle_definitions." + particle_name + ".isotopes");
    for (const auto& node : isotopes_nodes) {
        // BUG FIX: Use the renamed IsotopeDef struct
        m_particleInfo.isotopes.push_back({node["name"].as<std::string>(), node["mass"].as<int>()});
    }

    // Parse samples
    auto samples_nodes = m_configManager->GetNode("samples_to_process");
    for (const auto& node : samples_nodes) {
        m_samples.push_back({
            node["name"].as<std::string>(),
            node["type"].as<std::string>(),
            node["files"].as<std::vector<std::string>>(),
            node["histograms"].as<std::vector<std::string>>(),
        });
    }
}

const ParticleInfo& AnalysisContext::GetParticleInfo() const {
    return m_particleInfo;
}

const std::vector<Sample>& AnalysisContext::GetSamplesToProcess() const {
    return m_samples;
}

const Sample& AnalysisContext::GetSample(const std::string& name) const {
    for (const auto& sample : m_samples) {
        if (sample.name == name) {
            return sample;
        }
    }
    throw std::runtime_error("Sample with name '" + name + "' not found in configuration.");
}

} // namespace IsoToolbox