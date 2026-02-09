#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include <memory>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TChain.h"
#include "basic_var.h"

// Namespace matches your project
namespace AMS_Iso {

using H1Ptr = std::unique_ptr<TH1F>;
using H2Ptr = std::unique_ptr<TH2F>;
using H3Ptr = std::unique_ptr<TH3F>;

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
    std::vector<std::vector<H2Ptr>> IDH5a4; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b2; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b3; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b4; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5c3; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5c4; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5d3; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5d4; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH6a; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH6b; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH7; // [chain][det]
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
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H1b2; // Equ3 Denom
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H1c; // Equ4 Denom

    // ---- H2 Series: frag corr num samples ----
    
    // H2b: Mass distribution (1/Mass vs Ek) for fitting
    // Dimensions: [chain][source][det]
    std::vector<std::vector<std::vector<H2Ptr>>> BKG_H2b;
    // H2b: Mass distribution (1/Mass vs Ek) for fitting with L1Q window cut and standard selection except L1Q
    std::vector<std::vector<std::vector<H2Ptr>>> BKG_H2b2;

    // H2a2: MC Unique - Fragment Isotope Counts
    // Dimensions: [chain][det][iso] (Counts of specific isotopes from MC truth)
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H2a; 
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H2a2;

    // ---- H3 Series: Background Estimation (MC Only) ----
    // Dimensions: [chain][det][iso]
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H3a;
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_H3b;

    // ---- H4 Series: Charge Study ----
    // Dimensions: [chain][source][det][charge_type]
    std::vector<std::vector<std::vector<std::vector<H2Ptr>>>> BKG_H4;

    // ---- H5 Series: L1-Inner 2D Q Study ----
    std::vector<std::vector<std::vector<H3Ptr>>> BKG_H5; // [chain][det][type]

    // ==========================================
    // ======= FLUX AREA ===========
    // ==========================================

    // FLUXH1: Efficiency samples 
    // Dimensions: [cut_group][num_den][det][iz]
    std::vector<std::vector<std::vector<std::vector<H1Ptr>>>> FLUXH1; 

    // ISS Flux Aux
    std::vector<H1Ptr> ISS_FLUXH2; // rig expoT
    std::vector<std::vector<H1Ptr>> ISS_FLUXH3; // ek expoT [det][iso]
    //old bt check
    std::vector<H2Ptr> ISS_FLUXH4; // Generated counts
    std::vector<H2Ptr> ISS_FLUXH5; // cutoff rig vs measure rig

    // MC Flux Aux
    //std::vector<std::vector<std::vector<std::vector<H1Ptr>>>> MC_FLUXH2; // [chain][cut][det][gen/rec]
    std::vector<H1Ptr> MC_FLUXH3; // Generated counts

    // ==========================================
    // ======= FRAG STUDY AREA (Added) ==========
    // ==========================================
    
    // --- 1. 变量分布研究 (TH2F: Var vs Rigidity) ---
    // [type: 0=AboveL1, 1=BelowL1][geo: 0=TOF, 1=NaF, 2=AGL]
    std::vector<std::vector<H2Ptr>> BKG_FRAG_UTOFQ;

    // [type: 0=AboveL1, 1=BelowL1][geo_rich: 0=NaF, 1=AGL][var: 0..7]
    // Var Index: 0:LTOFQ, 1:richQ, 2:rich_pmt, 3:rich_pb, 4:rich_npe_ratio, 5:rich_used_ratio, 6:rich_good, 7:rich_clean
    std::vector<std::vector<std::vector<H2Ptr>>> BKG_FRAG_RICH;

    // [type: 0=AboveL1, 1=BelowL1][var: 0..1]
    // Var Index: 0:tof_chisc, 1:tof_chist
    std::vector<std::vector<H2Ptr>> BKG_FRAG_TOF;

    // --- 2. cutBackground 排除能力研究 (TH1F: Counts vs Rigidity) ---
    // [type: 0=AboveL1, 1=BelowL1][geo: 0=TOF, 1=NaF, 2=AGL][num_den: 0=Total, 1=PassBkg]
    std::vector<std::vector<std::vector<H1Ptr>>> BKG_FRAG_REJ;


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