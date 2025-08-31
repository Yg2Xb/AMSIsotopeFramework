#include <IsoToolbox/AnalysisContext.h>
#include <stdexcept>

namespace IsoToolbox {

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

// BUG FIX: The implementation for GetParticleInfo() was missing.
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