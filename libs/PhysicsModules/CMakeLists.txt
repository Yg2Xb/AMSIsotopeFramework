# CMakeLists.txt for the PhysicsModules collection

# 1. 创建一个“虚拟的”INTERFACE库来组织所有的物理模块。
add_library(PhysicsModules INTERFACE)

# 2. 添加子模块目录
add_subdirectory(DataModel)
# 未来我们会在这里添加更多模块

# 3. 将子模块的公开属性“注入”到总接口中
target_link_libraries(PhysicsModules
    INTERFACE       # <-- 这里是修正的地方，从PUBLIC改为INTERFACE
        DataModel
        # Selection
)