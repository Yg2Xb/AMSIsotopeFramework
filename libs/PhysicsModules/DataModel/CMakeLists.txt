# CMakeLists.txt for the DataModel library

# 1. Define the library and its source files
add_library(DataModel SHARED
    src/AMSDstTreeA.cpp   
)

# 2. Specify the public include directory
target_include_directories(DataModel
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# 3. Link the library against ROOT components
target_link_libraries(DataModel
    PUBLIC
        ROOT::Tree
        ROOT::Core
)