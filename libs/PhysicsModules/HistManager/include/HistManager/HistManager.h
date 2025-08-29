#pragma once
#include "IsoToolbox/AnalysisContext.h"
#include "DataModel/AMSDstTreeA.h"
#include <vector>
#include <memory>

class TH1F;
class TFile;

namespace HistManager {

class HistManager {
public:
    HistManager(const IsoToolbox::AnalysisContext& context);
    ~HistManager();
    
    void fill(const AMSDstTreeA& event, int isotope_idx, int detector_idx);
    void save(const std::string& output_path);

private:
    // 3 detectors x N isotopes
    std::vector<std::vector<std::unique_ptr<TH1F>>> h_event_counts_; 
};

}