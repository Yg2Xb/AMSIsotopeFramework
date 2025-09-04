// File: IsoProcess/src/ProductRegistry.cpp
// Purpose: Implements the registration for multiple analysis chains (e.g., L1Inner, unbiasedL1Inner).

#include "ProductRegistry.h"
#include <iostream>
#include <algorithm>

ProductRegistry& ProductRegistry::GetInstance() {
    static ProductRegistry instance;
    return instance;
}

ProductRegistry::ProductRegistry() {}

const std::vector<HistogramBlueprint>& ProductRegistry::GetBlueprints() const {
    return m_blueprints;
}

// This helper now returns a vector of base histogram templates for ISS.
std::vector<HistogramBlueprint> ProductRegistry::GetISSBaseTemplates() {
    std::vector<HistogramBlueprint> templates;
    // === 1. Isotope Identification ===
    templates.push_back({"ISS.ID.H1", "Event counts vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Counts", "ID"});
    templates.push_back({"ISS.ID.H2", "1/M vs. E_k/n", "TH2F", "InverseMass", "EkPerNucleon", "1/Mass (1/GeV)", "E_k/n (GeV/n)", "ID"});
    templates.push_back({"ISS.ID.H3", "Heaviest Isotope 1/M vs. E_k/n (Cutoff)", "TH2F", "InverseMass", "EkPerNucleon", "1/Mass (1/GeV)", "E_k/n (GeV/n)", "ID"});
    templates.push_back({"ISS.ID.H4a", "1/beta RICH-NaF (beta~1)", "TH1F", "", "Beta", "1/#beta", "Counts", "ID"});
    templates.push_back({"ISS.ID.H4b", "1/beta RICH-AGL (beta~1)", "TH1F", "", "Beta", "1/#beta", "Counts", "ID"});
    templates.push_back({"ISS.ID.H5a", "d(1/beta) NaF-Trk vs. Rigidity", "TH2F", "Rigidity", "DeltaBeta", "Rigidity (GV)", "1/#beta_{NaF} - 1/#beta_{Trk}", "ID"});
    templates.push_back({"ISS.ID.H5b", "d(1/beta) AGL-Trk vs. Rigidity", "TH2F", "Rigidity", "DeltaBeta", "Rigidity (GV)", "1/#beta_{AGL} - 1/#beta_{Trk}", "ID"});
    templates.push_back({"ISS.ID.H6a", "d(1/beta) TOF-NaF vs. E_k/n", "TH2F", "EkPerNucleon", "DeltaBeta", "E_k/n (GeV/n)", "1/#beta_{TOF} - 1/#beta_{NaF}", "ID"});
    templates.push_back({"ISS.ID.H6b", "d(1/beta) TOF-AGL vs. E_k/n", "TH2F", "EkPerNucleon", "DeltaBeta", "E_k/n (GeV/n)", "1/#beta_{TOF} - 1/#beta_{AGL}", "ID"});
    templates.push_back({"ISS.ID.H7a", "d(1/beta) TOF-NaF vs. betagamma", "TH2F", "BetaGamma", "DeltaBeta", "#beta#gamma", "1/#beta_{TOF} - 1/#beta_{NaF}", "ID"});
    templates.push_back({"ISS.ID.H7b", "d(1/beta) TOF-AGL vs. betagamma", "TH2F", "BetaGamma", "DeltaBeta", "#beta#gamma", "1/#beta_{TOF} - 1/#beta_{AGL}", "ID"});

    // === 2. Background Subtraction ===
    templates.push_back({"ISS.BKG.H1", "L1/L2 Charge Template vs. E_k/n", "TH2F", "Charge", "EkPerNucleon", "Charge", "E_k/n (GeV/n)", "BKG"});
    templates.push_back({"ISS.BKG.H2", "Source Counts (L1Q) vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Counts", "BKG"});
    templates.push_back({"ISS.BKG.H3", "L2 Frag Counts (InnerQ) vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Counts", "BKG"});
    templates.push_back({"ISS.BKG.H4", "L2 Frag 1/M vs. E_k/n", "TH2F", "InverseMass", "EkPerNucleon", "1/Mass (1/GeV)", "E_k/n (GeV/n)", "BKG"});

    // === 3. Flux Calculation ===
    templates.push_back({"ISS.FLUX.H1", "Exposure Time vs. Rigidity", "TH1F", "", "Rigidity", "Rigidity (GV)", "Exposure Time (s)", "FLUX"});
    templates.push_back({"ISS.FLUX.H2", "Exposure Time vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Exposure Time (s)", "FLUX"});
    templates.push_back({"ISS.FLUX.H3", "Cut Efficiency Num/Den vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Counts", "FLUX"});
    return templates;
}

// This helper now returns a vector of base histogram templates for MC.
std::vector<HistogramBlueprint> ProductRegistry::GetMCBaseTemplates() {
    std::vector<HistogramBlueprint> templates;
    // === 1. Isotope Identification ===
    templates.push_back({"MC.ID.H1a", "Mass Template (Sel) 1/M vs. E_k/n", "TH2F", "InverseMass", "EkPerNucleon", "1/Mass (1/GeV)", "E_k/n (GeV/n)", "ID"});
    templates.push_back({"MC.ID.H1b", "Mass Template (Sel+NoTofFrag) 1/M vs. E_k/n", "TH2F", "InverseMass", "EkPerNucleon", "1/Mass (1/GeV)", "E_k/n (GeV/n)", "ID"});
    templates.push_back({"MC.ID.H2", "Heaviest Isotope 1/M vs. E_k/n (Cutoff)", "TH2F", "InverseMass", "EkPerNucleon", "1/Mass (1/GeV)", "E_k/n (GeV/n)", "ID"});
    templates.push_back({"MC.ID.H3a", "1/beta RICH-NaF (beta~1)", "TH1F", "", "Beta", "1/#beta", "Counts", "ID"});
    templates.push_back({"MC.ID.H3b", "1/beta RICH-AGL (beta~1)", "TH1F", "", "Beta", "1/#beta", "Counts", "ID"});
    templates.push_back({"MC.ID.H4a", "d(1/beta) NaF-Trk vs. Rigidity", "TH2F", "Rigidity", "DeltaBeta", "Rigidity (GV)", "1/#beta_{NaF} - 1/#beta_{Trk}", "ID"});
    templates.push_back({"MC.ID.H4b", "d(1/beta) AGL-Trk vs. Rigidity", "TH2F", "Rigidity", "DeltaBeta", "Rigidity (GV)", "1/#beta_{AGL} - 1/#beta_{Trk}", "ID"});
    templates.push_back({"MC.ID.H5a", "d(1/beta) TOF-NaF vs. E_k/n", "TH2F", "EkPerNucleon", "DeltaBeta", "E_k/n (GeV/n)", "1/#beta_{TOF} - 1/#beta_{NaF}", "ID"});
    templates.push_back({"MC.ID.H5b", "d(1/beta) TOF-AGL vs. E_k/n", "TH2F", "EkPerNucleon", "DeltaBeta", "E_k/n (GeV/n)", "1/#beta_{TOF} - 1/#beta_{AGL}", "ID"});
    templates.push_back({"MC.ID.H6a", "d(1/beta) TOF-NaF vs. betagamma", "TH2F", "BetaGamma", "DeltaBeta", "#beta#gamma", "1/#beta_{TOF} - 1/#beta_{NaF}", "ID"});
    templates.push_back({"MC.ID.H6b", "d(1/beta) TOF-AGL vs. betagamma", "TH2F", "BetaGamma", "DeltaBeta", "#beta#gamma", "1/#beta_{TOF} - 1/#beta_{AGL}", "ID"});

    // === 2. Background Subtraction ===
    templates.push_back({"MC.BKG.H1", "L1 Source Counts vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Counts", "BKG"});
    templates.push_back({"MC.BKG.H2", "L2 Frag Counts vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Counts", "BKG"});
    templates.push_back({"MC.BKG.H3a", "upTOF Frag Counts vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Counts", "BKG"});
    templates.push_back({"MC.BKG.H3b", "upTOF Frag (RICH Survived) vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Counts", "BKG"});

    // === 3. Flux Calculation ===
    templates.push_back({"MC.FLUX.H1", "Accumulated Passed Events vs. E_k", "TH1F", "", "Ek", "E_k (GeV)", "Counts", "FLUX"});
    templates.push_back({"MC.FLUX.H2", "Cut Efficiency Num/Den vs. E_k/n", "TH1F", "", "EkPerNucleon", "E_k/n (GeV/n)", "Counts", "FLUX"});
    templates.push_back({"MC.FLUX.H3", "Total Generated Events vs. E_k,gen", "TH1F", "", "EkGen", "E_{k,gen} (GeV/n)", "Counts", "FLUX"});
    return templates;
}

// The main registration function, now fully upgraded for multiple analysis chains.
void ProductRegistry::RegisterProducts(bool isMC,
                                     const std::vector<std::string>& active_chains,
                                     const std::vector<std::string>& active_groups) {
    if (m_isRegistered) return;
    if (active_chains.empty()) {
        std::cerr << "Warning: No active chains provided to ProductRegistry. No histograms will be registered." << std::endl;
        return;
    }

    // 1. Get the base templates (ISS or MC)
    std::vector<HistogramBlueprint> base_templates = isMC ? GetMCBaseTemplates() : GetISSBaseTemplates();

    // 2. Loop over each requested analysis chain (e.g., "unbiasedL1Inner", "L1Inner")
    for (const auto& chain : active_chains) {
        // 3. Loop over each base template
        for (auto base_bp : base_templates) { // Make a copy of the blueprint
            // 4. Apply filtering by category/group first (if any groups are specified)
            if (!active_groups.empty()) {
                if (std::find(active_groups.begin(), active_groups.end(), base_bp.category) == active_groups.end()) {
                    continue; // Skip this blueprint if its category is not in the active list
                }
            }

            // 5. Create the chain-specific blueprint by adding prefixes
            base_bp.uniqueID = chain + "_" + base_bp.uniqueID;
            base_bp.title = "[" + chain + "] " + base_bp.title;

            // 6. Add the final, fully-named blueprint to our list
            m_blueprints.push_back(base_bp);
        }
    }

    m_isRegistered = true;
    std::cout << "Registered " << m_blueprints.size() << " total histogram blueprints for " << active_chains.size() << " chain(s)." << std::endl;
}