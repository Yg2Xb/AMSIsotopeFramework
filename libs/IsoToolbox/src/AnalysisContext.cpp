#include "IsoToolbox/AnalysisContext.h"
#include <stdexcept>

namespace IsoToolbox {

AnalysisContext::AnalysisContext(const std::string& config_path) : m_particleInfo(nullptr) { // 初始化 unique_ptr
    try {
        m_configNode = YAML::LoadFile(config_path);
    } catch (const YAML::Exception& e) {
        throw std::runtime_error("Failed to load or parse config file '" + config_path + "': " + e.what());
    }
    parseConfig();
}

void AnalysisContext::parseConfig() {
    // 1. 解析要分析的目标粒子名称
    std::string particle_name = m_configNode["run_settings"]["target_nucleus"].as<std::string>();
    
    // 2. 从 C++ 代码中的 PhysicsConstants 加载该粒子的固定物理信息
    loadParticleData(particle_name);

    // 3. 解析用户选择的、灵活可变的“分析链”
    const auto& chain_node = m_configNode["analysis_chain"];
    m_analysisChain.chain_id           = chain_node["chain_id"].as<std::string>();
    m_analysisChain.rigidity_version = chain_node["rigidity_version"].as<std::string>();
    m_analysisChain.velocity_version = chain_node["velocity_version"].as<std::string>();
    m_analysisChain.cut_version      = chain_node["cut_version"].as<std::string>();

    // 4. 解析待处理的数据样本列表
    const auto& samples_nodes = m_configNode["samples_to_process"];
    for (const auto& node : samples_nodes) {
        m_samples.push_back({
            node["name"].as<std::string>(),
            node["type"].as<std::string>()
        });
    }
}

// 修正后的逻辑：直接获取并存储由 PhysicsConstants 提供的 Isotope 对象
void AnalysisContext::loadParticleData(const std::string& particleName) {
    // 创建 Isotope 类的一个新实例并将其所有权转移给 m_particleInfo
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
    return m_configNode["run_settings"];
}

} // namespace IsoToolbox