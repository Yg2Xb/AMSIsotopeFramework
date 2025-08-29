# CMakeLists.txt for the Selection library
add_library(Selection SHARED
    src/CutBase.cpp # <-- ADD this new source file
    src/RTICut.cpp
)
target_include_directories(Selection
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(Selection
    PUBLIC
        DataModel
        IsoToolbox # <-- 新增
        ROOT::Core
)