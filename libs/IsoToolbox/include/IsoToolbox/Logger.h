#pragma once

#include <memory>
#include <string>

// Forward declaration for the spdlog logger class
namespace spdlog {
    class logger;
}

namespace IsoToolbox {

    class Logger {
    public:
        enum class Level {
            DEBUG,
            INFO,
            WARN,
            ERROR
        };

        // Get the singleton instance of the logger
        static Logger& GetInstance();

        // Set the logging level
        void SetLevel(Level level);

        // Static logging methods for easy access
        template<typename... Args>
        static void Debug(const std::string& format, Args&&... args);

        template<typename... Args>
        static void Info(const std::string& format, Args&&... args);
        
        template<typename... Args>
        static void Warn(const std::string& format, Args&&... args);

        template<typename... Args>
        static void Error(const std::string& format, Args&&... args);


    private:
        Logger(); // Private constructor for singleton
        ~Logger();

        // Delete copy and move constructors/assignments
        Logger(const Logger&) = delete;
        Logger& operator=(const Logger&) = delete;
        Logger(Logger&&) = delete;
        Logger& operator=(Logger&&) = delete;

        std::shared_ptr<spdlog::logger> m_logger;
    };

} // namespace IsoToolbox

// Include the implementation of the templated methods
#include "impl/Logger_impl.h"