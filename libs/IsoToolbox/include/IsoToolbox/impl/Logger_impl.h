#pragma once

// This file provides the implementation for the templated static methods
// in Logger.h. It should only be included at the very end of Logger.h.

#include "IsoToolbox/Logger.h"
#include "spdlog/spdlog.h" // Include the actual spdlog header

namespace IsoToolbox {

template<typename... Args>
void Logger::Debug(const std::string& format, Args&&... args) {
    GetInstance().m_logger->debug(format, std::forward<Args>(args)...);
}

template<typename... Args>
void Logger::Info(const std::string& format, Args&&... args) {
    GetInstance().m_logger->info(format, std::forward<Args>(args)...);
}

template<typename... Args>
void Logger::Warn(const std::string& format, Args&&... args) {
    GetInstance().m_logger->warn(format, std::forward<Args>(args)...);
}

template<typename... Args>
void Logger::Error(const std::string& format, Args&&... args) {
    GetInstance().m_logger->error(format, std::forward<Args>(args)...);
}

} // namespace IsoToolbox