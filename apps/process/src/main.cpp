#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <nlohmann/json.hpp> // This include is now found correctly

// ... The rest of the main.cpp file from my previous response was already correct and robust.
// It correctly implements the detector loop and task-driven workflow. No changes needed from that version.
// For clarity, I am including it again.

// Framework headers
#include <IsoToolbox/Logger.h>
#include <IsoToolbox/AnalysisContext.h>
#include <IsoToolbox/Tools.h>
#include <PhysicsModules/HistManager.h>
#include <DataModel/AMSDstTreeA.h>
#include <Selection/RTICut.h>

// ROOT headers
#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TMath.h>

// A simple structure to hold the task information
struct Task {
    std::string sample_name;
    std::vector<std::string> input_files;
    std::string output_file;
};

void show_usage(const std::string& name) {
    std::cerr << "Usage: " << name << " <physics_config.yaml> --task <task_description.json>" << std::endl;
}

int main(int argc, char** argv) {
    if (argc != 4 || std::string(argv[2]) != "--task") {
        show_usage(argv[0]);
        return 1;
    }

    std::string physics_config_path = argv[1];
    std::string task_path = argv[3];

    // --- 1. Parse Task File ---
    Task current_task;
    try {
        std::ifstream f(task_path);
        nlohmann::json data = nlohmann::json::parse(f);
        current_task.sample_name = data.at("sample_name");
        current_task.input_files = data.at("input_files").get<std::vector<std::string>>();
        current_task.output_file = data.at("output_file");
    } catch (const std::exception& e) {
        std::cerr << "Error parsing task file " << task_path << ": " << e.what() << std::endl;
        return 1;
    }

    // --- 2. Initialization with Physics Config ---
    IsoToolbox::Logger::Initialize();
    LOG_INFO("Starting job for sample '{}' from task file '{}'", current_task.sample_name, task_path);

    auto context = std::make_shared<IsoToolbox::AnalysisContext>(physics_config_path);
    const auto& sample_config = context->GetSample(current_task.sample_name);
    const auto& particle_info = context->GetParticleInfo();

    // --- 3. Setup Physics Modules ---
    auto hist_manager = std::make_unique<PhysicsModules::HistManager>(*context);
    auto rti_cut = std::make_unique<PhysicsModules::RTICut>(*context);
    
    hist_manager->initializeForSample(sample_config);

    // --- 4. Event Loop ---
    auto chain = std::make_unique<TChain>("AMSDstTreeA");
    for (const auto& file_path : current_task.input_files) {
        chain->Add(file_path.c_str());
    }

    auto data = std::make_unique<DataModel::AMSDstTreeA>(chain.get());
    long n_entries = data->fChain->GetEntries();
    LOG_INFO("Processing {} entries.", n_entries);

    for (long i = 0; i < n_entries; ++i) {
        data->fChain->GetEntry(i);

        if (!rti_cut->IsPass(data)) continue;
        if (data->rigidity[1] <= 0) continue;
        
        const std::vector<std::pair<std::string, float>> detectors = {
            {"TOF", data->beta},
            {"NaF", data->beta_rich[0]},
            {"AGL", data->beta_rich[1]}
        };

        for (const auto& det : detectors) {
            const std::string& detector_name = det.first;
            const float beta_val = det.second;

            if (!IsoToolbox::Tools::isValidBeta(beta_val)) continue;

            auto mass_result = IsoToolbox::Tools::calculateMass(beta_val, 1.0, data->rigidity[1], particle_info.charge);
            if (mass_result.invMass <= 0) continue;
            
            double ek_per_nucleon = mass_result.ek;
            int mass_number = TMath::Nint(1.0 / mass_result.invMass);

            std::string selected_isotope_name = "";
            for (const auto& iso : particle_info.isotopes) {
                if (iso.mass == mass_number) {
                    selected_isotope_name = iso.name;
                    break;
                }
            }

            if (!selected_isotope_name.empty()) {
                hist_manager->fill(Form("ISS.ID.H1.CountsVsEk.%s.%s", selected_isotope_name.c_str(), detector_name.c_str()), ek_per_nucleon);
                hist_manager->fill(Form("ISS.ID.H2.InvMassVsEk.%s.%s", selected_isotope_name.c_str(), detector_name.c_str()), ek_per_nucleon, mass_result.invMass);
            }
        }
    }
    
    // --- 5. Write Output ---
    hist_manager->write(current_task.output_file);

    LOG_INFO("Job finished successfully.");
    return 0;
}