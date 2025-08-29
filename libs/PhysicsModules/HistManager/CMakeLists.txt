# CMakeLists.txt for the HistManager library
add_library(HistManager SHARED
    src/HistManager.cpp
)
target_include_directories(HistManager
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(HistManager
    PUBLIC
        DataModel
        IsoToolbox # <-- 新增
        ROOT::Hist
        ROOT::Core
)