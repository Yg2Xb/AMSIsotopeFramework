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
    std::string name;
    int charge;
    std::vector<Isotope> isotopes;
};

// 修正: 添加了之前缺失的 Sample 结构体定义
struct Sample {
    std::string name;
    std::string type;
    std::vector<std::string> files;
    std::vector<std::string> histograms;
};

class AnalysisContext {
public:
    AnalysisContext(std::shared_ptr<ConfigManager> config);

    // 修正: 添加了之前缺失的 GetParticleInfo 方法
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
