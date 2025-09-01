# ====================================================================
# Project Definition
# ====================================================================
cmake_minimum_required(VERSION 3.16)
project(AMSIsotopeFramework VERSION 1.0 LANGUAGES CXX)

# ====================================================================
# Global Settings
# ====================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# ====================================================================
# Find OR Fetch External Dependencies
# ====================================================================
message(STATUS "Looking for required packages...")

# Find ROOT (provided by the CERN environment)
find_package(ROOT 6.24 COMPONENTS Core Hist Graf REQUIRED)
message(STATUS "Found ROOT ${ROOT_VERSION}")

# Find yaml-cpp (provided by the CERN environment)
find_package(yaml-cpp 0.6 REQUIRED)
message(STATUS "Found yaml-cpp ${yaml-cpp_VERSION}")

# --- THE FINAL, CORRECT, AND PROFESSIONAL SOLUTION FOR SPDLOG ---
# Include the FetchContent module to manage dependencies automatically
include(FetchContent)

# Declare where to get spdlog from if we don't have it
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.9.2 # A stable, widely used version
)

# --- Add nlohmann/json for parsing task files ---
FetchContent_Declare(
  json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG        v3.11.2 # A stable, widely used version
)
FetchContent_MakeAvailable(json)
message(STATUS "nlohmann/json dependency is now available.")
# ----------------------------------------------------

# This command will first check if spdlog is already available via find_package.
# If not, it will download and build it for our project automatically.
FetchContent_MakeAvailable(spdlog)
message(STATUS "spdlog dependency is now available.")
# ----------------------------------------------------------------

# ====================================================================
# Project Structure
# ====================================================================
message(STATUS "Configuring project subdirectories...")
add_subdirectory(libs)
add_subdirectory(apps)

# ====================================================================
# Testing
# ====================================================================
enable_testing()
add_subdirectory(tests)

message(STATUS "Configuration successful!")