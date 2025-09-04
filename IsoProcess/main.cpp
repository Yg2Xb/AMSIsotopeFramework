/***********************************************************
 *  File: main.cpp
 *  Main program for AMS Isotopes analysis
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/
#include <iostream>
#include <stdexcept>
#include "IsotopeAnalyzer.h"
#include "selectdata.h"
#include "basic_var.h"

using namespace AMS_Iso;

int main(int argc, char *argv[]) {
    try {
        // 参数检查
        if (argc < 4) {
            std::cerr << "Usage: " << argv[0] 
                     << " <outDir> <outName> <inData> <charge> [UseMass]" 
                     << std::endl;
            return 1;
        }

        // 获取分析器实例
        auto& analyzer = IsotopeAnalyzer::Instance();
        
        // 设置配置
        analyzer.setConfig(
            argv[1],                    // outDir
            argv[2],                    // outName
            argv[3],                    // inData
            argv[4],                    // charge
            argc >= 6 ? atoi(argv[5]) : -1  // UseMass 
        );
        
        // 初始化分析器
        analyzer.initialize();
        
        // 创建数据处理对象
        selectdata processor(analyzer.getDataChain());
        processor.SetAnalyzer(&analyzer);
        
        // 处理数据
        processor.Loop();
        
        // 保存结果
        analyzer.write();
        
        // 清理
        analyzer.cleanup();
        
        std::cout << "Analysis completed successfully" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error in analysis: " << e.what() << std::endl;
        return 1;
    }
}