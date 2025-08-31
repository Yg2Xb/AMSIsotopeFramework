#pragma once

#include <string>
#include <vector>
#include <map>
#include <functional>
#include <stdexcept>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <IsoToolbox/AnalysisContext.h>

namespace PhysicsModules {

class HistManager {
public:
    HistManager(const IsoToolbox::AnalysisContext& context);
    ~HistManager();

    void initializeForSample(const IsoToolbox::Sample& sample);

    // FIX: Renamed methods to resolve ambiguity.
    void Fill1D(const std::string& key, double value, double weight = 1.0);
    void Fill2D(const std::string& key, double x_value, double y_value, double weight = 1.0);

    void write(const std::string& output_path);

private:
    void registerAllBlueprints();
    void bookFromBlueprint(const std::string& base_key);

    const IsoToolbox::AnalysisContext& m_context;
    std::map<std::string, TH1*> m_hists;
    
    std::map<std::string, std::function<void()>> m_registry;
};

} // namespace PhysicsModules