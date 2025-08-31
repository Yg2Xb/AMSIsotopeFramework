#pragma once

#include <string>
#include <vector>
#include <memory>
#include <IsoToolbox/ConfigManager.h>

namespace IsoToolbox {

// --- Data Structures that define the analysis context ---

// BUG FIX: Renamed from 'Isotope' to 'IsotopeDef' to avoid redefinition conflict
// with the class in PhysicsConstants.h
struct IsotopeDef {
    std::string name;
    int mass;
};

struct ParticleInfo {
    std::string name;
    int charge;
    std::vector<IsotopeDef> isotopes;
};

struct Sample {
    std::string name;
    std::string type;
    std::vector<std::string> files;
    std::vector<std::string> histograms;
};

class AnalysisContext {
public:
    // BUG FIX: Added a new constructor that takes the config file path directly.
    // This simplifies instantiation and resolves the constructor mismatch error in main.
    AnalysisContext(const std::string& config_path);
    
    AnalysisContext(std::shared_ptr<ConfigManager> config);

    const ParticleInfo& GetParticleInfo() const;
    const std::vector<Sample>& GetSamplesToProcess() const;
    const Sample& GetSample(const std::string& name) const;

private:
    void parseConfig();

    std::shared_ptr<ConfigManager> m_configManager;
    ParticleInfo m_particleInfo;
    std::vector<Sample> m_samples;
};

} // namespace IsoToolbox