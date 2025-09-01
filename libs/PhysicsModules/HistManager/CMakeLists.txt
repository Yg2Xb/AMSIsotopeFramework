add_library(HistManager)

target_include_directories(HistManager
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_sources(HistManager
    PRIVATE
        src/HistManager.cpp
)

target_link_libraries(HistManager
    PUBLIC
        # Our headers include TH1.h, so ROOT::Hist must be public.
        # We also need IsoToolbox for ConfigManager, etc.
        IsoToolbox
        ROOT::Hist
)