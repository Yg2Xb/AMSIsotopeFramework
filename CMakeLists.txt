# Root CMakeLists.txt for the AMSIsotopeFramework project

# 1. Set minimum required version of CMake and project name
cmake_minimum_required(VERSION 3.10)
project(AMSIsotopeFramework CXX)

# Ensure all targets are built with Position Independent Code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 2. Find ROOT package and its required components
find_package(ROOT REQUIRED
    COMPONENTS
        Core
        RIO
        Net
        Tree
        Graf
        Graf3d
        Gpad
        Hist
        RooFit
        RooFitCore
        MathCore
        Matrix
)

# Include ROOT's CMake utilities
include(${ROOT_USE_FILE})

# 3. Add crucial compile option for compatibility with CERN's ROOT builds
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)


# 4. Configure FetchContent module to manage external dependencies
include(FetchContent)

# --- spdlog ---
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.x
)
FetchContent_MakeAvailable(spdlog)

# --- yaml-cpp ---
FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG        0.8.0
)
FetchContent_MakeAvailable(yaml-cpp)

# --- nlohmann_json ---  <-- 新增下面这个代码块
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.2 # A recent, stable version
)
FetchContent_MakeAvailable(nlohmann_json)


# 5. Set the C++ standard for the entire project
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 6. Add the subdirectories that contain libraries and applications
add_subdirectory(libs)
add_subdirectory(apps)

# 7. Enable testing with CTest
enable_testing()
add_subdirectory(tests)
