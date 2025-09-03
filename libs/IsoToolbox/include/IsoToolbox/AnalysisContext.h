#pragma once

#include <string>
#include <vector>
#include <memory>
#include "IsoToolbox/PhysicsConstants.h"
#include "yaml-cpp/yaml.h"

namespace IsoToolbox {

// Defines the analysis chain's specific physics model versions
struct AnalysisChain {
    std::string chain_id;
    std::string rigidity_version;
    std::string beta_version; // Changed from velocity_version
    std::string cut_version;
};

// Defines information for a data sample to be processed
struct SampleInfo {
    std::string name;
    std::string type; // "ISS" or "MC"
};

class AnalysisContext {
public:
    explicit AnalysisContext(const std::string& config_path);
    
    const Isotope& GetParticleInfo() const;
    const AnalysisChain& GetAnalysisChain() const;
    const std::vector<SampleInfo>& GetSamplesToProcess() const;
    const YAML::Node GetRunSettings() const;

private:
    void parseConfig();
    void loadParticleData(const std::string& particleName);

    YAML::Node m_configNode;
    std::unique_ptr<Isotope> m_particleInfo; 
    AnalysisChain m_analysisChain;
    std::vector<SampleInfo> m_samples;
};

} // namespace IsoToolbox
