# CMakeLists.txt for the PhysicsModules collection

# 1. 创建“虚拟”接口库
add_library(PhysicsModules INTERFACE)

# 2. 添加子模块目录
add_subdirectory(DataModel)
add_subdirectory(Selection)      # <-- 新增
add_subdirectory(HistManager)    # <-- 新增

# 3. 将子模块链接到总接口上
target_link_libraries(PhysicsModules
    INTERFACE
        DataModel
        Selection            # <-- 新增
        HistManager          # <-- 新增
)