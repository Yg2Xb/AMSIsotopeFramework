#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include <memory>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "basic_var.h"

// Namespace matches your project
namespace AMS_Iso {

using H1Ptr = std::unique_ptr<TH1F>;
using H2Ptr = std::unique_ptr<TH2F>;

class HistManager {
public:
    // Constructor
    HistManager(const std::string& output_filename,
                bool isISS,
                const std::vector<std::string>& chains,
                int charge,
                const IsotopeVar* iso,
                int UseMass,
                int FragmentZ); 

    // Prepare output tree (called before event loop)
    void PrepareFilteredTree(TChain* dataChain);
    
    // Get tree for filling
    TTree* GetFilteredTree() const { return m_filteredTree; }
    
    // Save method
    void Save(bool saveTree = true);

    // ==========================================
    // ======== ID AREA =========
    // ==========================================
    
    // ISS
    std::vector<std::vector<std::vector<H1Ptr>>> ISS_IDH1; // [chain][det][iso]
    
    // MC: Extended to Niso dim
    std::vector<std::vector<std::vector<H2Ptr>>> MC_IDH1; // [chain][det][iso]

    std::vector<std::vector<std::vector<H2Ptr>>> IDH2; // [chain][det][iso]
    std::vector<std::vector<H2Ptr>> IDH3; // [chain][det]

    // Beta study - with Z dimension
    std::vector<std::vector<H1Ptr>> IDH4a; // [z][chain]
    std::vector<std::vector<H1Ptr>> IDH4b; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5a; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5a2; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5a3; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b2; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b3; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH6a; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH6b; // [z][chain]
    //std::vector<std::vector<H2Ptr>> IDH7a; // [z][chain]
    //std::vector<std::vector<H2Ptr>> IDH7b; // [z][chain]

    // ==========================================
    // ======= BKG AREA =========
    // ==========================================

    // ---- H1 Series: frag corr den samples ----
    // Dimensions: [chain][source][det]
    // MC source dim is 1, ISS source dim is sources.size()
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H1a; // Equ1 Denom
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H1b; // Equ2 Denom
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H1c; // Equ3 Denom

    // ---- H2 Series: frag corr num samples ----
    
    // H2a: Fragment to Element Y counts
    // Dimensions: [chain][source][det]
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H2a; 

    // H2b: Mass distribution (1/Mass vs Ek) for fitting
    // Dimensions: [chain][source][det]
    std::vector<std::vector<std::vector<H2Ptr>>> BKG_H2b;

    // H2a2: MC Unique - Fragment Isotope Counts
    // Dimensions: [chain][det][iso] (Counts of specific isotopes from MC truth)
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H2a2;

    // ---- H3 Series: Background Estimation (MC Only) ----
    // Dimensions: [chain][det][iso]
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H3a;
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H3b;

    // ---- H4 Series: Charge Study ----
    // Dimensions: [chain][source][det][charge_type]
    std::vector<std::vector<std::vector<std::vector<H2Ptr>>>> BKG_H4;

    // ==========================================
    // ======= FLUX AREA ===========
    // ==========================================

    // FLUXH1: Efficiency samples (Equ4)
    // Dimensions: [chain][cut_group][num_den][det][iso]
    std::vector<std::vector<std::vector<std::vector<std::vector<H1Ptr>>>>> FLUXH1; 

    // ISS Flux Aux
    std::vector<H1Ptr> ISS_FLUXH2; // rig expoT
    std::vector<std::vector<H1Ptr>> ISS_FLUXH3; // ek expoT [det][iso]

    // MC Flux Aux
    std::vector<std::vector<std::vector<std::vector<H1Ptr>>>> MC_FLUXH2; // [chain][cut][det][gen/rec]
    std::vector<H1Ptr> MC_FLUXH3; // Generated counts

private:
    std::unique_ptr<TFile> m_outputFile;
    TTree* m_filteredTree = nullptr; 

    // Define branches to save
    static std::vector<std::string> GetActiveBranches();

    // Helper to create hists
    template <typename HistType, typename... Args>
    std::unique_ptr<HistType> createHist(Args&&... args) {
        auto hist = new HistType(std::forward<Args>(args)...);
        hist->SetDirectory(nullptr);
        return std::unique_ptr<HistType>(hist);
    }
};

} // namespace AMS_Iso

#endif