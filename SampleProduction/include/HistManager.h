#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include <memory>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "basic_var.h"

namespace AMS_Iso {

class HistManager {
public:
    HistManager(const std::string& output_filename,
                bool isISS,
                const std::vector<std::string>& chains,
                int charge,
                const IsotopeVar* iso,
                int UseMass);

    ~HistManager();

    void Save();

    // ======== ID 区域 ========
    std::vector<std::vector<std::vector<TH1F*>>> ISS_IDH1; // [chain][det][iso]
    std::vector<std::vector<TH2F*>> MC_IDH1;               // [chain][det]

    std::vector<std::vector<std::vector<TH2F*>>> IDH2;     // [chain][det][iso/1]
    std::vector<std::vector<TH2F*>> IDH3;                  // [chain][det]
    std::vector<TH1F*> IDH4a; // [chain]
    std::vector<TH1F*> IDH4b; // [chain]
    std::vector<TH2F*> IDH5a; // [chain]
    std::vector<TH2F*> IDH5b; // [chain]
    std::vector<TH2F*> IDH6a; // [chain]
    std::vector<TH2F*> IDH6b; // [chain]
    std::vector<TH2F*> IDH7a; // [chain]
    std::vector<TH2F*> IDH7b; // [chain]

    // ======== BKG 区域 ========
    // ISS
    std::vector<std::vector<std::vector<TH1F*>>> ISS_BKGH1;                // [chain][source][det]
    std::vector<std::vector<std::vector<std::vector<TH2F*>>>> ISS_BKGH2;  // [chain][source][charge_type][det]
    std::vector<std::vector<std::vector<TH1F*>>> ISS_BKGH3;                // [chain][source][det]
    std::vector<std::vector<std::vector<TH2F*>>> ISS_BKGH4;                // [chain][source][det]
    // MC
    std::vector<std::vector<TH1F*>> MC_BKGH1;               // [chain][det]
    std::vector<std::vector<std::vector<TH1F*>>> MC_BKGH2;  // [chain][det][iso=1]
    std::vector<std::vector<std::vector<TH1F*>>> MC_BKGH3a; // [chain][det][iso=1]
    std::vector<std::vector<std::vector<TH1F*>>> MC_BKGH3b; // [chain][det][iso=1]

    // ======== FLUX 区域 ========
    std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*>>>>> FLUXH1; // [chain][cut_group][num_den][det][iso/1]
    std::vector<TH1F*> FLUXH2;  // [chain]
    std::vector<std::vector<std::vector<TH1F*>>> FLUXH3; // [chain][det][iso]

    // MC 独有
    std::vector<TH1F*> MC_FLUXH2; // [chain]
    std::vector<TH1F*> MC_FLUXH3; // [chain]

private:
    std::unique_ptr<TFile> m_outputFile;
};

} // namespace AMS_Iso

#endif