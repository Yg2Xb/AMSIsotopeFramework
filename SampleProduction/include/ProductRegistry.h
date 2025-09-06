#ifndef PRODUCT_REGISTRY_H
#define PRODUCT_REGISTRY_H

#include <string>
#include <vector>
#include <map>
#include "IsotopeAnalyzer.h"

namespace AMS_Iso {

struct HistogramBlueprint {
    std::string name;
    std::string title;
    std::string type;
    std::string binning_x_key;
    std::string binning_y_key;
    std::string x_axis_title;
    std::string y_axis_title;
};

struct ProductTemplate {
    std::string template_id;
    std::string title_pattern;
    std::string name_pattern;
    std::string type;
    std::string binning_x_key;
    std::string binning_y_key;
    std::string x_title_pattern;
    std::string y_title_pattern;
    std::vector<std::string> multipliers;
    std::map<std::string, std::vector<std::string>> default_values;
};


class ProductRegistry {
public:
    static ProductRegistry& GetInstance();

    void SetAnalyzerInfo(const IsotopeAnalyzer* analyzer);
    void RegisterProducts(bool isMC, const std::vector<std::string>& active_chains);
    const std::vector<HistogramBlueprint>& GetBlueprints() const;

private:
    ProductRegistry(); // MODIFIED: Declare constructor
    void RegisterCoreTemplates_ISS();
    void RegisterCoreTemplates_MC();
    void ExpandTemplates(const std::vector<std::string>& active_chains);
    static std::string ReplacePatterns(const std::string& pattern, const std::map<std::string, std::string>& replacements);

    const IsotopeAnalyzer* m_analyzer = nullptr;
    std::map<std::string, ProductTemplate> m_templates;
    std::vector<HistogramBlueprint> m_blueprints;
    bool m_isRegistered = false;
};

} // namespace AMS_Iso

#endif // PRODUCT_REGISTRY_H
