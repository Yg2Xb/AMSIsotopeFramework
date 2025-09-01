// apps/process/src/main.cpp
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <exception>
#include <fstream>

#include "TChain.h"
#include "nlohmann/json.hpp"

#include "IsoToolbox/AnalysisContext.h"
#include "IsoToolbox/BinningManager.h"
#include "IsoToolbox/Logger.h"
#include "IsoToolbox/ProductRegistry.h"  // 新增
#include "DataModel/AMSDstTreeA.h"
#include "HistManager/HistManager.h"

using json = nlohmann::json;

int main(int argc, char* argv[]) {
    IsoToolbox::Logger::GetInstance().SetLevel(IsoToolbox::Logger::Level::DEBUG);

    if (argc != 4 || std::string(argv[2]) != "--task") {
        IsoToolbox::Logger::Error("Usage: {} <config_file> --task <task_file.json>", argv[0]);
        return 1;
    }

    std::string configFile = argv[1];
    std::string taskFile = argv[3];

    try {
        IsoToolbox::Logger::Info("======= AMS Isotope Analysis Framework - IsoProcess =======");
        IsoToolbox::Logger::Info("Master Config: {}", configFile);
        IsoToolbox::Logger::Info("Task File: {}", taskFile);

        // 1. 解析任务文件 (JSON)
        std::ifstream f(taskFile);
        if (!f.is_open()) {
            throw std::runtime_error("Could not open task file: " + taskFile);
        }
        json task_data = json::parse(f);
        
        const std::vector<std::string> inputFiles = task_data["input_files"];
        const std::string outputFile = task_data["output_file"];
        const std::string sampleName = task_data["sample_name"];

        // 2. 创建核心分析组件
        auto context = std::make_shared<IsoToolbox::AnalysisContext>(configFile);
        auto binningManager = std::make_shared<IsoToolbox::BinningManager>(context.get());
        auto histManager = std::make_shared<IsoToolbox::HistManager>();
        
        // 3. 初始化组件
        binningManager->Initialize();
        
        // 确保ProductRegistry已注册模板
        auto& registry = IsoToolbox::ProductRegistry::GetInstance();
        
        // 使用新的简化接口预订直方图
        histManager->BookHistograms(context.get(), binningManager.get());
        
        const auto& particleInfo = context->GetParticleInfo();
        const auto& analysisChain = context->GetAnalysisChain();
        const auto& runSettings = context->GetRunSettings();

        IsoToolbox::Logger::Info("Processing Sample: {}", sampleName);
        IsoToolbox::Logger::Info("Analysis Target: {} (Z={})", particleInfo.GetName(), particleInfo.GetCharge());
        IsoToolbox::Logger::Info("Analysis Chain -> Rigidity: {}, Velocity: {}, Cuts: {}",
                                 analysisChain.rigidity_version, 
                                 analysisChain.velocity_version, 
                                 analysisChain.cut_version);

        // 4. 事件循环
        auto chain = std::make_unique<TChain>("amstreea");
        for (const auto& inputFile : inputFiles) {
            chain->Add(inputFile.c_str());
        }
        
        AMSDstTreeA dst(chain.get());
        long long nEntries = chain->GetEntries();
        if (nEntries == 0) {
            IsoToolbox::Logger::Warn("No events found in the input files for this task. Generating empty output.");
        } else {
            IsoToolbox::Logger::Info("Starting event loop for {} events...", nEntries);
            long long maxEvents = runSettings["max_events"].as<long long>(-1);
            if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;

            for (long long i = 0; i < nEntries; ++i) {
                chain->GetEntry(i);
                if (i > 0 && i % 100000 == 0) {
                    IsoToolbox::Logger::Info("  ... processing event {} / {}", i, nEntries);
                }
                
                // 示例：填充背景分析直方图
                // 这里应该根据你的事件选择逻辑来填充相应的直方图
                
                // 模拟一些事件数据填充
                for (const auto& source : {"Carbon", "Nitrogen", "Oxygen"}) {
                    for (const auto& detector : {"TOF", "NaF", "AGL"}) {
                        double fake_ek = (i % 1000) / 100.0 + 0.5;  // 0.5-10.5 GeV/n
                        
                        // 填充ISS.BKG.H2 (源粒子计数)
                        std::string hist_name = "ISS_BKG_H2_" + std::string(source) + "_" + std::string(detector);
                        histManager->Fill1D(hist_name, fake_ek);
                        
                        // 填充ISS.BKG.H3 (碎片产物计数)
                        hist_name = "ISS_BKG_H3_" + std::string(source) + "_" + std::string(detector);
                        histManager->Fill1D(hist_name, fake_ek);
                    }
                }
            }
        }
        
        IsoToolbox::Logger::Info("Event loop finished. Processed {} events.", nEntries);
        histManager->SaveHistograms(outputFile);

    } catch (const std::exception& e) {
        IsoToolbox::Logger::Error("An unrecoverable error occurred: {}", e.what());
        return 1;
    }

    IsoToolbox::Logger::Info("Task finished successfully.");
    return 0;
}