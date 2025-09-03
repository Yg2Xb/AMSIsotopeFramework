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
#include "DataModel/StandardizedEvent.h"
#include "DataModel/AMSEventBuilder.h"
#include "HistManager/HistManager.h"
#include "DataModel/AMSDstTreeA.h"

using json = nlohmann::json;

int main(int argc, char* argv[]) {
    IsoToolbox::Logger::GetInstance().SetLevel(IsoToolbox::Logger::Level::INFO);

    if (argc != 4 || std::string(argv[2]) != "--task") {
        IsoToolbox::Logger::Error("Usage: {} <config_file> --task <task_file.json>", argv[0]);
        return 1;
    }

    std::string configFile = argv[1];
    std::string taskFile = argv[3];

    try {
        IsoToolbox::Logger::Info("======= AMS Isotope Analysis Framework - IsoProcess =======");
        auto context = std::make_shared<IsoToolbox::AnalysisContext>(configFile);
        auto binningManager = std::make_shared<IsoToolbox::BinningManager>(context.get());
        auto histManager = std::make_shared<IsoToolbox::HistManager>();
        
        binningManager->Initialize();
        histManager->BookHistograms(context.get(), binningManager.get());
        
        json task_data;
        std::ifstream f(taskFile);
        f >> task_data;
        
        const std::vector<std::string> inputFiles = task_data["input_files"];
        const std::string outputFile = task_data["output_file"];
        const auto& runSettings = context->GetRunSettings();
        const auto& analysisChain = context->GetAnalysisChain();

        IsoToolbox::Logger::Info("Analysis Chain -> Rigidity: {}, Beta: {}, Cuts: {}",
                                 analysisChain.rigidity_version, 
                                 analysisChain.beta_version, 
                                 analysisChain.cut_version);
        
        auto chain = std::make_unique<TChain>("amstreea");
        for (const auto& inputFile : inputFiles) {
            chain->Add(inputFile.c_str());
        }
        
        long long nEntries = chain->GetEntries();
        if (nEntries > 0) {
            PhysicsModules::AMSEventBuilder builder(chain.get(), context.get());
            StandardizedEvent event;

            IsoToolbox::Logger::Info("Starting event loop for {} events...", nEntries);
            long long maxEvents = runSettings["run_settings"]["max_events"].as<long long>(-1);
            if (maxEvents > 0 && maxEvents < nEntries) nEntries = maxEvents;

            for (long long i = 0; i < nEntries; ++i) {
                chain->GetEntry(i);
                if (i > 0 && i % 100000 == 0) {
                    IsoToolbox::Logger::Info("  ... processing event {} / {}", i, nEntries);
                }
                
                builder.Fill(event);
                const AMSDstTreeA* raw_event = builder.GetRawEvent();

                // --- Real Analysis Logic Placeholder ---
                if (event.rigidity > 0) {
                    histManager->Fill1D("ISS_BKG_H2_Carbon_TOF", event.ek_tof);
                }
            }
        }
        
        histManager->SaveHistograms(outputFile);

    } catch (const std::exception& e) {
        IsoToolbox::Logger::Error("An unrecoverable error occurred: {}", e.what());
        return 1;
    }

    return 0;
}

