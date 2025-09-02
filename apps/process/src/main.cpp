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
#include "IsoToolbox/ProductRegistry.h"
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

        // 4. 从配置中读取背景源列表
        std::vector<std::string> background_sources;
        if (runSettings["template_config"] && runSettings["template_config"]["background_sources"]) {
            background_sources = runSettings["template_config"]["background_sources"].as<std::vector<std::string>>();
            IsoToolbox::Logger::Info("Using background sources from config: [{}]", fmt::join(background_sources, ", "));
        } else {
            // 默认值作为后备
            background_sources = {"Carbon", "Nitrogen", "Oxygen"};
            IsoToolbox::Logger::Warn("template_config.background_sources not found, using defaults: [{}]", 
                                   fmt::join(background_sources, ", "));
        }
        
        // 5. 探测器列表
        std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};

        IsoToolbox::Logger::Info("Processing Sample: {}", sampleName);
        IsoToolbox::Logger::Info("Analysis Target: {} (Z={})", particleInfo.GetName(), particleInfo.GetCharge());
        IsoToolbox::Logger::Info("Analysis Chain -> Rigidity: {}, Velocity: {}, Cuts: {}",
                                 analysisChain.rigidity_version, 
                                 analysisChain.velocity_version, 
                                 analysisChain.cut_version);

        // 6. 事件循环
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
            long long maxEvents = runSettings["run_settings"]["max_events"].as<long long>(-1);
            if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;

            for (long long i = 0; i < nEntries; ++i) {
                chain->GetEntry(i);
                if (i > 0 && i % 100000 == 0) {
                    IsoToolbox::Logger::Info("  ... processing event {} / {}", i, nEntries);
                }
                
                // 模拟事件数据
                double fake_ek = (i % 1000) / 100.0 + 0.5;  // 0.5-10.5 GeV/n
                double fake_charge = 4.0 + (i % 100) / 100.0; // 4.0-5.0
                double fake_invmass = 0.1 + (i % 50) / 500.0;  // 0.1-0.2
                
                for (const auto& source : background_sources) {
                    for (const auto& detector : detectors) {
                        
                        // 填充 ISS.BKG.H1 (2D: 电荷 vs Ek/n) - 需要3种电荷类型
                        std::vector<std::string> charge_types = {"Signal", "Template", "L2Template"};
                        for (const auto& charge_type : charge_types) {
                            std::string h1_name = "ISS_BKG_H1_" + source + "_" + charge_type + "_" + detector;
                            histManager->Fill2D(h1_name, fake_charge, fake_ek);  // x=charge, y=ek
                        }
                        
                        // 填充 ISS.BKG.H2 (1D: 源粒子计数)
                        std::string h2_name = "ISS_BKG_H2_" + source + "_" + detector;
                        histManager->Fill1D(h2_name, fake_ek);
                        
                        // 填充 ISS.BKG.H3 (1D: L2碎片计数)  
                        std::string h3_name = "ISS_BKG_H3_" + source + "_" + detector;
                        histManager->Fill1D(h3_name, fake_ek);
                        
                        // 填充 ISS.BKG.H4 (2D: 1/M vs Ek/n)
                        std::string h4_name = "ISS_BKG_H4_" + source + "_" + detector;
                        histManager->Fill2D(h4_name, fake_invmass, fake_ek);  // x=1/M, y=ek
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