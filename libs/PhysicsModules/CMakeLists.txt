# CMakeLists.txt for the PhysicsModules collection

# 1. Create a virtual INTERFACE library to group all physics modules.
#    It doesn't compile any code itself, it only collects properties.
add_library(PhysicsModules INTERFACE)

# 2. Add subdirectories for each module.
add_subdirectory(DataModel)
# add_subdirectory(Selection) # etc.

# 3. Link the sub-module libraries to the main interface.
#    This crucial step propagates the PUBLIC properties (like include paths)
#    from DataModel up to the PhysicsModules interface.
target_link_libraries(PhysicsModules
    PUBLIC
        DataModel
        # Selection # etc.
)