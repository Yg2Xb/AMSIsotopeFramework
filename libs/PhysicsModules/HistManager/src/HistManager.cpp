#include <HistManager/HistManager.h>
#include <TFile.h>
#include <TString.h>
#include <IsoToolbox/Logger.h> // BUG FIX: Missing include for logging macros

namespace PhysicsModules {

HistManager::HistManager(const IsoToolbox::AnalysisContext& context) : m_context(context) {
    registerAllBlueprints();
}

HistManager::~HistManager() {
    for (auto const& [key, val] : m_hists) {
        delete val;
    }
}

void HistManager::registerAllBlueprints() {
    auto particle_info = m_context.GetParticleInfo(); // This call now works
    std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};

    for (const auto& isotope : particle_info.isotopes) {
        for (const auto& detector : detectors) {
            std::string h1_key = Form("ISS.ID.H1.CountsVsEk.%s.%s", isotope.name.c_str(), detector.c_str());
            // BUG FIX: Corrected lambda capture. The variables used inside must be captured.
            m_registry[h1_key] = [this, h1_key, isotope, detector]() {
                std::string title = Form("Event counts for %s (%s) vs. E_k/n; E_k/n (GeV/n); Counts", isotope.name.c_str(), detector.c_str());
                m_hists[h1_key] = new TH1F(h1_key.c_str(), title.c_str(), 100, 0.5, 10.0);
            };

            std::string h2_key = Form("ISS.ID.H2.InvMassVsEk.%s.%s", isotope.name.c_str(), detector.c_str());
            // BUG FIX: Corrected lambda capture.
            m_registry[h2_key] = [this, h2_key, isotope, detector]() {
                std::string title = Form("1/M for %s (%s) vs. E_k/n; E_k/n (GeV/n); 1/M (amu^-1)", isotope.name.c_str(), detector.c_str());
                m_hists[h2_key] = new TH2F(h2_key.c_str(), title.c_str(), 100, 0.5, 10.0, 100, 0.08, 0.15);
            };
        }
    }
}

void HistManager::bookFromBlueprint(const std::string& base_key) {
    auto particle_info = m_context.GetParticleInfo(); // This call now works
    std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};

    if (base_key == "ISS.ID.H1" || base_key == "ISS.ID.H2") {
        std::string desc = (base_key == "ISS.ID.H1") ? "CountsVsEk" : "InvMassVsEk";
        for (const auto& isotope : particle_info.isotopes) {
            for (const auto& detector : detectors) {
                std::string full_key = Form("%s.%s.%s.%s", base_key.c_str(), desc.c_str(), isotope.name.c_str(), detector.c_str());
                if (m_registry.count(full_key)) {
                    m_registry.at(full_key)();
                    if(m_hists.count(full_key)) m_hists.at(full_key)->SetDirectory(nullptr);
                } else {
                    SPDLOG_WARN("Blueprint for key '{}' not found in registry.", full_key);
                }
            }
        }
    } else {
        if (m_registry.count(base_key)) {
            m_registry.at(base_key)();
            if(m_hists.count(base_key)) m_hists.at(base_key)->SetDirectory(nullptr);
        } else {
            SPDLOG_WARN("Blueprint for key '{}' not found in registry.", base_key);
        }
    }
}

void HistManager::initializeForSample(const IsoToolbox::Sample& sample) {
    for (auto const& [key, val] : m_hists) { delete val; }
    m_hists.clear();

    for (const auto& hist_base_key : sample.histograms) {
        bookFromBlueprint(hist_base_key);
    }
    SPDLOG_INFO("Booked {} histograms for sample '{}'.", m_hists.size(), sample.name);
}
// ... the rest of the file (fill and write methods) is correct and remains unchanged.
void HistManager::Fill1D(const std::string& key, double value, double weight) {
    try {
        m_hists.at(key)->Fill(value, weight);
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Attempted to fill non-existent histogram with key: " + key);
    }
}

void HistManager::Fill2D(const std::string& key, double x_value, double y_value, double weight) {
    try {
        dynamic_cast<TH2*>(m_hists.at(key))->Fill(x_value, y_value, weight);
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Attempted to fill non-existent 2D histogram with key: " + key);
    }
}

void HistManager::write(const std::string& output_path) {
    TFile* output_file = TFile::Open(output_path.c_str(), "RECREATE");
    if (!output_file || output_file->IsZombie()) {
         throw std::runtime_error("Failed to open output file: " + output_path);
    }
    output_file->cd();
    for (auto const& [key, hist] : m_hists) {
        hist->Write();
    }
    output_file->Close();
    delete output_file;
    SPDLOG_INFO("Histograms successfully written to {}", output_path);
}

} // namespace PhysicsModules