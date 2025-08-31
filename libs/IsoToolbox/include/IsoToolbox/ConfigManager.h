#pragma once

#include <string>
#include "yaml-cpp/yaml.h"
#include "IsoToolbox/Logger.h"

namespace IsoToolbox {

class ConfigManager {
public:
    bool Load(const std::string& filepath);

    // A simpler, more robust Get method
    template<typename T>
    T Get(const std::string& key) const;

    // Method to get the sub-node for nested parameters
    YAML::Node GetNode(const std::string& key) const;

    const YAML::Node& GetRootNode() const;

    static ConfigManager& GetInstance();

private:
    ConfigManager() = default;
    ConfigManager(const ConfigManager&) = delete;
    ConfigManager& operator=(const ConfigManager&) = delete;

    YAML::Node m_config;
};

// --- Template Implementation for simple keys ---
template<typename T>
T ConfigManager::Get(const std::string& key) const {
    try {
        return m_config[key].as<T>();
    } catch (const YAML::Exception& e) {
        Logger::Error("Error accessing config key '{}': {}", key, e.what());
        throw;
    }
}

} // namespace IsoToolbox