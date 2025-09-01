#pragma once

#include <string>
#include <vector>
#include <memory>
#include "IsoToolbox/PhysicsConstants.h" // 包含真实的类型定义
#include "yaml-cpp/yaml.h"

namespace IsoToolbox {

// 结构体：定义了分析链中使用的具体物理模型版本
struct AnalysisChain {
    std::string chain_id;
    std::string rigidity_version;
    std::string velocity_version;
    std::string cut_version;
};

// 结构体：定义了要处理的数据样本信息
struct SampleInfo {
    std::string name;
    std::string type; // "ISS" or "MC"
};

class AnalysisContext {
public:
    // 构造函数：接收配置文件路径，并加载所有分析上下文
    explicit AnalysisContext(const std::string& config_path);
    
    // 修正：返回类型为 Isotope, 这是 PhysicsConstants 中定义的真实类型
    const Isotope& GetParticleInfo() const;

    // 获取当前选择的分析链配置
    const AnalysisChain& GetAnalysisChain() const;
    
    // 获取待处理样本的信息
    const std::vector<SampleInfo>& GetSamplesToProcess() const;

    // 获取通用的运行设置
    const YAML::Node GetRunSettings() const;

private:
    void parseConfig();
    void loadParticleData(const std::string& particleName);

    YAML::Node m_configNode; // 持有整个配置文件的内容
    
    // 修正：成员变量的类型必须是 Isotope
    // 使用 std::unique_ptr 来存储，因为它没有默认构造函数
    std::unique_ptr<Isotope> m_particleInfo; 
    
    AnalysisChain m_analysisChain;
    std::vector<SampleInfo> m_samples;
};

} // namespace IsoToolbox