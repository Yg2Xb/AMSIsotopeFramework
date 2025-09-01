#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <exception>

#include "TFile.h"
#include "TTree.h"

#include "IsoToolbox/Logger.h"
#include "IsoToolbox/ConfigManager.h"
#include "IsoToolbox/AnalysisContext.h"
#include "IsoToolbox/BinningManager.h"
#include "PhysicsModules/DataModel/AMSDstTreeA.h"
#include "PhysicsModules/HistManager/HistManager.h"

int main(int argc, char* argv[]) {
    // Set logger level, e.g., to DEBUG for more verbose output
    Logger::GetInstance().SetLogLevel(LogLevel::INFO);

    if (argc != 3) {
        LOG_CRITICAL("Usage: %s <config_file> <input_file>", argv[0]);
        return 1;
    }

    std::string configFile = argv[1];
    std::string inputFile = argv[2];

    try {
        // 1. ====================== INITIALIZATION PHASE ======================
        LOG_INFO("======= AMS Isotope Analysis Framework =======");

        // -- Load Configuration
        ConfigManager config(configFile);

        // -- Setup Analysis Context (The brain of the operation)
        auto context = std::make_shared<AnalysisContext>();
        context->SetTargetNucleus(config.Get<std::string>("run_settings.target_nucleus"));
        context->Initialize();
        LOG_INFO("Analysis Target: %s (Z=%d)", context->GetTargetNucleus().c_str(), context->GetChargeZ());

        // -- Setup Binning Manager (The ruler maker)
        auto binningManager = std::make_shared<BinningManager>(context.get());
        binningManager->Initialize();

        // -- Book Histograms (The data containers)
        auto histManager = std::make_shared<HistManager>();
        histManager->BookHistograms(config, context.get(), binningManager.get());

        // -- Prepare Input Data
        auto inFile = std::make_unique<TFile>(inputFile.c_str(), "READ");
        if (!inFile || inFile->IsZombie()) {
            throw std::runtime_error("Error opening input file: " + inputFile);
        }
        TTree* tree = static_cast<TTree*>(inFile->Get("AMSDst"));
        if (!tree) {
            throw std::runtime_error("Could not find TTree 'AMSDst' in file: " + inputFile);
        }
        AMSDstTreeA dst(tree);

        // 2. ====================== EVENT LOOP PHASE ======================
        LOG_INFO("Starting event loop for %s...", inputFile.c_str());
        long long nEntries = tree->GetEntries();
        long long maxEvents = config.Get<long long>("run_settings.max_events", -1);
        if (maxEvents > 0 && maxEvents < nEntries) {
            nEntries = maxEvents;
        }

        for (long long i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);
            if (i > 0 && i % 100000 == 0) {
                LOG_INFO("Processing event %lld / %lld (%.1f%%)", i, nEntries, (double)i / nEntries * 100.0);
            }

            histManager->Fill1D("Event_Count_Total", 1.0);

            // ===============================================================
            // NEXT STEP: This is where the SELECTION logic will be added.
            // A 'SelectionManager' will be created and called here.
            // if (!selectionManager.ProcessEvent(dst)) continue;
            // ===============================================================

            // For now, fill histograms with fake data to test the system.
            for (const auto& isotope : context->GetIsotopes()) {
                // FAKE physics values for testing
                double fake_ek_per_nucleon = (i % 1000) / 100.0; // GeV/n
                double fake_inv_mass = (i % 200) / 400.0;

                std::string h1_name = "ISS_ID_H1_" + isotope.name;
                histManager->Fill1D(h1_name, fake_ek_per_nucleon);

                std::string h2_name = "ISS_ID_H2_" + isotope.name;
                histManager->Fill2D(h2_name, fake_inv_mass, fake_ek_per_nucleon);
            }
        }
        LOG_INFO("Event loop finished. Processed %lld events.", nEntries);

        // 3. ====================== FINALIZATION PHASE ======================
        std::string outputFile = config.Get<std::string>("run_settings.output_file");
        histManager->SaveHistograms(outputFile);

    } catch (const std::exception& e) {
        LOG_CRITICAL("An unrecoverable error occurred: %s", e.what());
        return 1;
    }

    LOG_INFO("Analysis finished successfully.");
    return 0;
}