#include "HistManager/HistManager.h"
#include "IsoToolbox/PhysicsConstants.h"
#include <TH1F.h>
#include <TFile.h>

namespace HistManager {

HistManager::HistManager(const IsoToolbox::AnalysisContext& context) {
    const auto& target = context.GetTargetIsotope();
    const auto& masses = target.GetMasses();
    const auto& rigidity_bins = IsoToolbox::PhysicsConstants::GetRigidityBins();
    
    h_event_counts_.resize(3); // 3 detectors
    for (int i_det = 0; i_det < 3; ++i_det) {
        h_event_counts_[i_det].resize(target.GetIsotopeCount());
        for (int i_iso = 0; i_iso < target.GetIsotopeCount(); ++i_iso) {
            
            // Auto-calculate Ek/n bins
            auto ek_bins = IsoToolbox::BinningUtils::ConvertRigidityToEk(rigidity_bins, masses[i_iso], target.GetCharge());
            
            std::string name = "h_event_counts_" + target.GetName() + std::to_string(masses[i_iso]) + "_det" + std::to_string(i_det);
            std::string title = "Event Counts for " + target.GetName() + std::to_string(masses[i_iso]);
            
            h_event_counts_[i_det][i_iso] = std::make_unique<TH1F>(name.c_str(), title.c_str(), ek_bins.size() - 1, ek_bins.data());
        }
    }
}

HistManager::~HistManager() = default;

void HistManager::fill(const AMSDstTreeA& event, int isotope_idx, int detector_idx) {
    // This is a simplified fill logic for now
    // We need to properly calculate Ek/n for the event
    h_event_counts_[detector_idx][isotope_idx]->Fill(1.0); // Placeholder
}

void HistManager::save(const std::string& output_path) {
    auto file = std::make_unique<TFile>(output_path.c_str(), "RECREATE");
    file->cd();
    for (const auto& vec : h_event_counts_) {
        for (const auto& hist : vec) {
            hist->Write();
        }
    }
    file->Close();
}

}