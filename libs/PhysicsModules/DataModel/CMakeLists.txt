add_library(DataModel)

target_include_directories(DataModel
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_sources(DataModel
    PRIVATE
        src/AMSDstTreeA.cpp
)

target_link_libraries(DataModel
    PUBLIC
        # We depend on ROOT for TTree, and our header includes ROOT headers.
        ROOT::Core
        # --- FIX ---
        # Added ROOT::Tree, which is required for TTree functionality used in AMSDstTreeA
        ROOT::Tree
)