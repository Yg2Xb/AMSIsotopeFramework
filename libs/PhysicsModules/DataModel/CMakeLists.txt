add_library(DataModel)

target_include_directories(DataModel
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# --- MODIFIED ---
# We now compile both AMSDstTreeA.cpp (which contains the Init function)
# and the new AMSEventBuilder.cpp.
target_sources(DataModel
    PRIVATE
        src/AMSDstTreeA.cpp 
        src/AMSEventBuilder.cpp
)

target_link_libraries(DataModel
    PUBLIC
        ROOT::Core
        ROOT::Tree
        IsoToolbox
)

