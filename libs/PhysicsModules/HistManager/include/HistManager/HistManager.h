// libs/PhysicsModules/HistManager/include/HistManager/HistManager.h
#pragma once

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "IsoToolbox/AnalysisContext.h"
#include "IsoToolbox/BinningManager.h"
#include "IsoToolbox/ProductRegistry.h"
#include <map>
#include <string>
#include <memory>

namespace IsoToolbox {

class HistManager {
public:
    HistManager();
    ~HistManager();

    // 新的简化接口 - 移除ConfigManager依赖
    void BookHistograms(const AnalysisContext* context, const BinningManager* binningManager);
    
    void Fill1D(const std::string& name, double value, double weight = 1.0);
    void Fill2D(const std::string& name, double x_val, double y_val, double weight = 1.0);
    void SaveHistograms(const std::string& output_path);

private:
    std::map<std::string, std::unique_ptr<TH1>> m_hists;
    
    void CreateHistogram(const HistogramBlueprint& blueprint, 
                        const std::vector<double>& x_bins,
                        const std::vector<double>& y_bins = {});
};

} // namespace IsoToolbox