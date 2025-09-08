# Set minimum required CMake version
cmake_minimum_required(VERSION 3.20.2)

# Define project root directory
set(PROJECT_ROOT_DIR ${CMAKE_SOURCE_DIR}/../..)

# 设置编译器
#set(CMAKE_CXX_COMPILER "/cvmfs/projects.cern.ch/intelsw/oneAPI/linux/x86_64/2023/compiler/2023.2.1/linux/bin/icpx")
#set(CMAKE_C_COMPILER  "/cvmfs/projects.cern.ch/intelsw/oneAPI/linux/x86_64/2023/compiler/2023.2.1/linux/bin/icx")

# ⚠️ 建议在开启 ASan 时用 g++/gcc，而不是 Intel 编译器
# cmake -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc ..

# Declaring the project
project(SampleProduction LANGUAGES CXX)

###########################################
# Basic Settings
###########################################
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

# 基本编译标志
set(CMAKE_CXX_FLAGS "-O3 ${ROOT_CFLAGS} -std=c++17 -Wall -Wextra")

# 输出目录
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

###########################################
# Dependencies
###########################################
find_package(ROOT REQUIRED COMPONENTS
    Core RIO Net Hist Graf Graf3d Gpad Tree Matrix Physics MathCore Thread
)

if(NOT ROOT_FOUND)
    message(FATAL_ERROR "ROOT not found")
endif()
include(${ROOT_USE_FILE})

###########################################
# Enable AddressSanitizer in Debug
###########################################
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building with AddressSanitizer enabled")
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer -g)
    add_link_options(-fsanitize=address)
endif()

###########################################
# GAMModel Library
###########################################
ROOT_GENERATE_DICTIONARY(G__GAMModel 
    ${CMAKE_SOURCE_DIR}/include/GAMModel.h
    LINKDEF ${CMAKE_SOURCE_DIR}/LinkDef.h
    MODULE GAMModel
    NOINSTALL
)

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
# Main Analysis Library
###########################################
set(ANALYSIS_SOURCES
    src/basic_var.cpp
    src/IsotopeAnalyzer.cpp
    src/ModelManager.cpp
    src/RICHCut.cpp
    src/RTICut.cpp
    src/selectdata.cpp      
    src/TOFCut.cpp
    src/Tool.cpp
    src/TrackerCut.cpp
    src/BinningManager.cpp
    src/HistManager.cpp
)

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
add_executable(SampleProduction main.cpp) 

target_include_directories(SampleProduction
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${ROOT_INCDIR}
)

target_link_libraries(SampleProduction
    PRIVATE
        isotope_analysis
        GAMModel
        ${ROOT_LIBS}
)

# RPATH settings
set_target_properties(SampleProduction PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${ROOT_LIBDIR}"
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

# Add PIC
set_target_properties(GAMModel PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(isotope_analysis PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Testing support
enable_testing()

###########################################
# Debug Info
###########################################
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Compiler flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "ROOT libs: ${ROOT_LIBS}")