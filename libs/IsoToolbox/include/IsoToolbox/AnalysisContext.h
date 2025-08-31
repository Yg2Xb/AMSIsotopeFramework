#pragma once

#include <string>
#include <vector>
#include <memory>
#include <IsoToolbox/ConfigManager.h>

namespace IsoToolbox {

// --- Data Structures that define the analysis context ---

struct Isotope {
    std::string name;
    int mass;
};

struct ParticleInfo {
    std.string name;
    int charge;
    std::vector<Isotope> isotopes;
};

// BUG FIX: The definition for 'Sample' was missing.
struct Sample {
    std::string name;
    std::string type;
    std::vector<std::string> files;
    std::vector<std::string> histograms;
};

class AnalysisContext {
public:
    AnalysisContext(std::shared_ptr<ConfigManager> config);

    // BUG FIX: The GetParticleInfo() method declaration was missing.
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