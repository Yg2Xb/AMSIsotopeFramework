#include "IsoToolbox/Logger.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace IsoToolbox {

    Logger& Logger::GetInstance() {
        static Logger instance;
        return instance;
    }

    Logger::Logger() {
        // Create a thread-safe, color-coded console logger
        m_logger = spdlog::stdout_color_mt("console");
        // Default level
        SetLevel(Level::INFO); 
    }

    Logger::~Logger() {
        spdlog::drop_all();
    }

    void Logger::SetLevel(Level level) {
        switch (level) {
            case Level::DEBUG:
                m_logger->set_level(spdlog::level::debug);
                break;
            case Level::INFO:
                m_logger->set_level(spdlog::level::info);
                break;
            case Level::WARN:
                m_logger->set_level(spdlog::level::warn);
                break;
            case Level::ERROR:
                m_logger->set_level(spdlog::level::err);
                break;
        }
    }

} // namespace IsoToolbox