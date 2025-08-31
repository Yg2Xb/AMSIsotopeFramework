#include "IsoToolbox/ConfigManager.h"

namespace IsoToolbox {

ConfigManager& ConfigManager::GetInstance() {
    static ConfigManager instance;
    return instance;
}

bool ConfigManager::Load(const std::string& filepath) {
    try {
        m_config = YAML::LoadFile(filepath);
        return true;
    } catch (const YAML::Exception& e) {
        Logger::Error("Failed to load config file '{}': {}", filepath, e.what());
        return false;
    }
}

YAML::Node ConfigManager::GetNode(const std::string& key) const {
    return m_config[key];
}
const YAML::Node& ConfigManager::GetRootNode() const {
    return m_config;
}

} // namespace IsoToolbox