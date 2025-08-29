#pragma once

// This is an implementation detail file and should not be included directly by users.
// It uses the popular 'spdlog' library for high-performance logging.
#include "spdlog/spdlog.h"
#include "IsoToolbox/Logger.h"

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