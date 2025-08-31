# CMakeLists.txt for the IsoToolbox library

# 1. Add all source files for the library
set(SOURCES
    src/Logger.cpp
    src/ConfigManager.cpp
    src/PhysicsConstants.cpp # <-- 新增
    src/AnalysisContext.cpp # <-- 新增
    src/SplineFitter.cpp
    src/Tools.cpp
)

# 2. Add all public header files for the library
set(HEADERS
    include/IsoToolbox/Logger.h
    include/IsoToolbox/ConfigManager.h
    include/IsoToolbox/PhysicsConstants.h # <-- 新增
    include/IsoToolbox/AnalysisContext.h # <-- 新增
    include/IsoToolbox/SplineFitter.h
    include/IsoToolbox/Tools.h
)

# 3. Create the shared library target
add_library(IsoToolbox SHARED ${SOURCES} ${HEADERS})

# 4. Specify the public include directory.
target_include_directories(IsoToolbox
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# 5. Link against external dependencies
target_link_libraries(IsoToolbox
    PUBLIC
        spdlog::spdlog
        yaml-cpp
        ROOT::Core
        ROOT::MathCore # <-- ADD THIS LINE
)