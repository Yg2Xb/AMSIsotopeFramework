#include "IsoToolbox/AnalysisContext.h"
#include "IsoToolbox/Logger.h"  // 添加Logger头文件
#include <stdexcept>

namespace IsoToolbox {

AnalysisContext::AnalysisContext(const std::string& config_path) : m_particleInfo(nullptr) {
    try {
        m_configNode = YAML::LoadFile(config_path);
        
        // 🔍 添加完整配置调试输出
        Logger::Debug("=== Full YAML Config Debug ===");
        YAML::Emitter emitter;
        emitter << m_configNode;
        Logger::Debug("Raw YAML content:\n{}", emitter.c_str());
        
        // 🔍 检查关键节点是否存在
        Logger::Debug("=== Key Nodes Check ===");
        Logger::Debug("Has 'run_settings': {}", m_configNode["run_settings"] ? "YES" : "NO");
        Logger::Debug("Has 'template_config': {}", m_configNode["template_config"] ? "YES" : "NO");
        Logger::Debug("Has 'active_templates': {}", m_configNode["active_templates"] ? "YES" : "NO");
        
        if (m_configNode["template_config"] && m_configNode["template_config"]["background_sources"]) {
            auto sources = m_configNode["template_config"]["background_sources"].as<std::vector<std::string>>();
            Logger::Debug("Found background_sources in AnalysisContext: [{}]", fmt::join(sources, ", "));
        } else {
            Logger::Warn("background_sources NOT found in AnalysisContext!");
        }
        Logger::Debug("=== End Config Debug ===");
        
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

    // 3. 解析用户选择的、灵活可变的"分析链"
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

// 🔧 关键修改：返回完整配置，而不只是run_settings子节点
const YAML::Node AnalysisContext::GetRunSettings() const {
    Logger::Debug("GetRunSettings() called, returning full config with {} top-level entries", m_configNode.size());
    return m_configNode;  // 返回完整配置！
}

} // namespace IsoToolbox