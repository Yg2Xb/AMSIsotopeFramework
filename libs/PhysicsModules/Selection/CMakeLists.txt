add_library(Selection)

target_include_directories(Selection
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_sources(Selection
    PRIVATE
        src/CutBase.cpp
        src/RTICut.cpp
)

target_link_libraries(Selection
    PUBLIC
        # Selection modules will depend heavily on the data model and tools
        DataModel
        IsoToolbox
)