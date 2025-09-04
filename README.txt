核心类和功能:

数据处理与分析框架:
    IsotopeAnalyzer - 主分析框架类,管理整个分析流程
    HistogramManager - 管理所有直方图的创建、填充和保存
    selectdata - ROOT数据读取与处理类
探测器select类:
    TrackerCut - 径迹探测器select
    TOFCut - 飞行时间探测器select
    RICHCut - 切伦科夫探测器select
    RTICut - 运行时间信息select
RICH探测器beta修正:
    ModelManager - 管理RICH beta修正模型
    GAMModel - RICH beta修正的具体模型实现
基础工具和常量:
    Tool - 通用工具函数集合
    basic_var - 基本物理常量和数据结构
    rich_var - RICH探测器相关常量
    tracker_var - 径迹探测器相关常量

工作流程:
    main.cpp 创建 IsotopeAnalyzer 实例
    IsotopeAnalyzer 初始化直方图管理器和数据处理
    selectdata 读取数据并通过各个探测器select类进行事例选择
    结果通过 HistogramManager 保存到输出文件