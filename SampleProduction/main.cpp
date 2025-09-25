/***********************************************************
 * File: main.cpp
 * Main program for AMS Isotopes analysis
 *
 * Modified to support multiple root files (TChain)
 * while keeping compatibility with single-file scripts
 *
 * History:
 * 20241029 - created by ZX.Yan
 * 20250909 - modified by ZX.Yan for multi-file chain support
 ***********************************************************/

#include <iostream>
#include <string>
#include "IsotopeAnalyzer.h"
#include "selectdata.h"
#include "basic_var.h"

using namespace AMS_Iso;

int main(int argc, char* argv[]) {
    try {
        // 参数检查
        if (argc < 5) {
            std::cerr << "Usage for ISS: " << argv[0]
                      << " <outDir> <outName> <inData> <charge|MC> [UseMass]"
                      << std::endl;
            return 1;
        }

        // 获取分析器实例
        auto& analyzer = IsotopeAnalyzer::Instance();

        std::string outDir  = argv[1];
        std::string outName = argv[2];
        std::string inData  = argv[3]; // 可以是单个文件，或.txt列表文件
        std::string charge  = argv[4];
        int useMass = (argc > 5) ? atoi(argv[5]) : 0;

        // 配置与初始化（内部会根据 inData 自动 readDataFrom）
        analyzer.setConfig(outDir, outName, inData, charge, useMass);
        analyzer.initialize();

        // 获取数据链并处理
        TChain* chain = analyzer.getDataChain();
        selectdata processor(chain);
        processor.SetAnalyzer(&analyzer);
        processor.Loop();

        // 保存结果与清理
        analyzer.write();
        analyzer.cleanup();

        std::cout << "Analysis completed successfully" << std::endl;
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error in analysis: " << e.what() << std::endl;
        return 1;
    }
}