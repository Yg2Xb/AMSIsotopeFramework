// File: IsoProcess/include/ProductRegistry.h
// Purpose: Updated to support multiple analysis chains (e.g., L1Inner, unbiasedL1Inner).

#ifndef PRODUCT_REGISTRY_H
#define PRODUCT_REGISTRY_H

#include <string>
#include <vector>

// This struct holds all the information needed to define one histogram.
// No changes needed here.
struct HistogramBlueprint {
    std::string uniqueID;
    std::string title;
    std::string type; // "TH1F" or "TH2F"
    std::string binningX;
    std::string binningY; // Empty for 1D histograms
    std::string axisTitleX;
    std::string axisTitleY;
    std::string category; // "ID", "BKG", "FLUX"
};

class ProductRegistry {
public:
    static ProductRegistry& GetInstance();

    // The main registration function is updated to accept a list of active chains.
    void RegisterProducts(bool isMC,
                          const std::vector<std::string>& active_chains,
                          const std::vector<std::string>& active_groups = {});

    const std::vector<HistogramBlueprint>& GetBlueprints() const;

private:
    ProductRegistry();
    
    // These helpers now return a vector of templates instead of modifying the member variable.
    std::vector<HistogramBlueprint> GetISSBaseTemplates();
    std::vector<HistogramBlueprint> GetMCBaseTemplates();

    std::vector<HistogramBlueprint> m_blueprints;
    bool m_isRegistered = false;
};

#endif // PRODUCT_REGISTRY_H