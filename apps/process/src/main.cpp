#include "IsoToolbox/Logger.h"
#include "IsoToolbox/ConfigManager.h"
#include "IsoToolbox/AnalysisContext.h"
#include "DataModel/AMSDstTreeA.h"
#include "Selection/RTICut.h"
#include "HistManager/HistManager.h"

#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <vector>

using Logger = IsoToolbox::Logger;

int main(int argc, char** argv) {
    
    Logger::GetInstance().SetLevel(Logger::Level::INFO);
    Logger::Info("Starting IsoProcess application...");

    // 1. Load configuration
    IsoToolbox::ConfigManager::GetInstance().Load("config/config_v1.0.yaml");
    auto config = IsoToolbox::ConfigManager::GetInstance().GetNode("");

    // 2. Setup Analysis Context
    std::string target = config["analysis_target"].as<std::string>();
    IsoToolbox::AnalysisContext context(target);
    Logger::Info("Analysis target set to: {}", context.GetTargetIsotope().GetName());

    // 3. Prepare data input
    TChain chain("amstreea");
    auto files = config["input_files"].as<std::vector<std::string>>();
    for (const auto& file : files) {
        chain.Add(file.c_str());
    }
    Logger::Info("Chain created with {} entries.", chain.GetEntries());
    
    AMSDstTreeA eventReader;
    eventReader.Init(&chain);

    // 4. Prepare Cuts and Histograms
    Selection::RTICut rtiCut;
    // ... we will add more cuts here ...
    
    HistManager::HistManager histManager(context);
    
    // 5. Event Loop
    Logger::Info("Starting event loop...");
    for (Long64_t i = 0; i < chain.GetEntries(); ++i) {
        chain.GetEntry(i);

        if (i % 10000 == 0) {
            Logger::Debug("Processing entry {}", i);
        }

        // Apply Cuts
        if (!rtiCut.pass(eventReader)) continue;
        // ... more cuts ...

        // Fill Histograms (simplified logic)
        // We assume isotope 0 (Be7) and detector 0 (TOF) for now
        histManager.fill(eventReader, 0, 0); 
    }
    Logger::Info("Event loop finished.");

    // 6. Save results
    std::string outputFile = config["output_file"].as<std::string>();
    histManager.save(outputFile);
    Logger::Info("Histograms saved to {}", outputFile);

    return 0;
}