# IsoToolbox - A library for common analysis tools

add_library(IsoToolbox)

# Public headers (available to other modules)
target_include_directories(IsoToolbox
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# Private headers (implementation details)
target_include_directories(IsoToolbox
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/IsoToolbox/impl>
)

# Source files
target_sources(IsoToolbox
    PRIVATE
        src/AnalysisContext.cpp
        src/ConfigManager.cpp
        src/Logger.cpp
        src/PhysicsConstants.cpp
        src/SplineFitter.cpp
        src/Tools.cpp
        src/BinningManager.cpp # <-- ADDED
)

# Dependencies
target_link_libraries(IsoToolbox
    PUBLIC
        ROOT::Core
        ROOT::Hist
        ROOT::Graf
        YAML::yaml-cpp
)