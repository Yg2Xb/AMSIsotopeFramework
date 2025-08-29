# CMakeLists.txt for the 'process' executable

# 1. Define the executable and its source files
add_executable(process
    src/main.cpp
)

# 2. Link the executable against our libraries.
#    We now link to the main PhysicsModules interface, which provides everything.
target_link_libraries(process
    PRIVATE
        IsoToolbox
        PhysicsModules  
)