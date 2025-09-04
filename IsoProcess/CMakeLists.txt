# Set minimum required CMake version
cmake_minimum_required(VERSION 3.20.2)

# Define project root directory
set(PROJECT_ROOT_DIR ${CMAKE_SOURCE_DIR}/../..)

# 设置编译器
set(CMAKE_CXX_COMPILER "/cvmfs/projects.cern.ch/intelsw/oneAPI/linux/x86_64/2023/compiler/2023.2.1/linux/bin/icpx")
set(CMAKE_C_COMPILER "/cvmfs/projects.cern.ch/intelsw/oneAPI/linux/x86_64/2023/compiler/2023.2.1/linux/bin/icx")

# Declaring the project
project(IsoProcess LANGUAGES CXX) # Changed project name for clarity

###########################################
# Basic Settings
###########################################
# 保持使用 C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ROOT flags
execute_process(
    COMMAND root-config --cflags
    OUTPUT_VARIABLE ROOT_CFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND root-config --libs
    OUTPUT_VARIABLE ROOT_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND root-config --libdir
    OUTPUT_VARIABLE ROOT_LIBDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND root-config --incdir
    OUTPUT_VARIABLE ROOT_INCDIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# 设置编译标志，匹配合作者的Makefile，但保持C++17
set(CMAKE_CXX_FLAGS "-O3 ${ROOT_CFLAGS} -std=c++17")

# 输出目录设置
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

###########################################
# Dependencies
###########################################
# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS
    Core RIO Net Hist Graf Graf3d Gpad Tree Matrix Physics MathCore Thread
)

include(${ROOT_USE_FILE})

###########################################
# GAMModel Library (UNCHANGED)
###########################################
# Generate dictionary
ROOT_GENERATE_DICTIONARY(G__GAMModel 
    ${CMAKE_SOURCE_DIR}/include/GAMModel.h
    LINKDEF ${CMAKE_SOURCE_DIR}/LinkDef.h
    MODULE GAMModel
    NOINSTALL
)

# GAMModel library
add_library(GAMModel SHARED
    src/GAMModel.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/G__GAMModel.cxx
)

target_include_directories(GAMModel 
    PUBLIC 
        ${CMAKE_SOURCE_DIR}/include
        ${ROOT_INCDIR}
)

target_link_libraries(GAMModel
    PUBLIC
        ${ROOT_LIBS}
)

###########################################
# Main Analysis Library (MODIFIED SECTION)
###########################################
# Collect source files
set(ANALYSIS_SOURCES
    src/basic_var.cpp
    # src/HistogramManager.cpp # <-- REMOVED
    src/IsotopeAnalyzer.cpp
    src/ModelManager.cpp
    src/RICHCut.cpp
    src/RTICut.cpp
    src/selectdata.cpp
    src/TOFCut.cpp
    src/Tool.cpp
    src/TrackerCut.cpp
    # --- ADDED NEW COMPONENTS ---
    src/ProductRegistry.cpp
    src/BinningManager.cpp
    src/HistManager.cpp
)

# Main analysis library
add_library(isotope_analysis SHARED ${ANALYSIS_SOURCES})

target_include_directories(isotope_analysis
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
        ${ROOT_INCDIR}
)

target_link_libraries(isotope_analysis
    PUBLIC
        GAMModel
        ${ROOT_LIBS}
)

###########################################
# Executable
###########################################
add_executable(IsoProcess main.cpp) # Renamed executable to IsoProcess

target_include_directories(IsoProcess
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${ROOT_INCDIR}
)

target_link_libraries(IsoProcess
    PRIVATE
        isotope_analysis
        GAMModel
        ${ROOT_LIBS}
)

# RPATH settings
set_target_properties(IsoProcess PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${ROOT_LIBDIR}"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

# 添加PIC标志
set_target_properties(GAMModel PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(isotope_analysis PROPERTIES POSITION_INDEPENDENT_CODE ON)

###########################################
# Debug Information
###########################################
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Project Root Dir: ${PROJECT_ROOT_DIR}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "ROOT flags: ${ROOT_CFLAGS}")
message(STATUS "ROOT libs: ${ROOT_LIBS}")
message(STATUS "ROOT libdir: ${ROOT_LIBDIR}")
message(STATUS "ROOT include dir: ${ROOT_INCDIR}")
message(STATUS "ROOT version: ${ROOT_VERSION}")