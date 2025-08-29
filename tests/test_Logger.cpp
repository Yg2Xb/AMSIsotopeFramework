#include "IsoToolbox/Logger.h"
#include <iostream>

using Logger = IsoToolbox::Logger;
// A minimal assertion macro for our tests

#define ASSERT(condition, message) \
    if (!(condition)) { \
        Logger::Error("--> ASSERT FAILED: {}", message); \
        return 1; \
    }

int main() {
    // We can use the logger even within our tests
    Logger::GetInstance().SetLevel(IsoToolbox::Logger::Level::INFO);
    Logger::Info("--- Running Logger Unit Test ---");

    // This is a "smoke test". It doesn't test complex logic, but
    // it verifies that we can instantiate and use the logger without crashing.
    // The real test is that this program compiles, links, and runs successfully.
    
    Logger::Info("Test: Calling logger with different levels.");
    Logger::Debug("This debug message should NOT be visible."); // Level is INFO
    Logger::Info("This info message should be visible.");

    std::cout << "[Test PASSED]" << std::endl;
    Logger::Info("--- Logger Unit Test PASSED ---");

    // Return 0 to indicate success to the CTest framework
    return 0;
}
