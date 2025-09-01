// libs/IsoToolbox/src/ProductRegistry.cpp
#include "IsoToolbox/ProductRegistry.h"
#include "IsoToolbox/Logger.h"
#include <stdexcept>
#include <algorithm>

namespace IsoToolbox {

ProductRegistry& ProductRegistry::GetInstance() {
    static ProductRegistry instance;
    if (instance.m_templates.empty()) {
        instance.RegisterCoreTemplates();
    }
    return instance;
}

void ProductRegistry::RegisterCoreTemplates() {
    Logger::Info("Registering core product templates...");
    
    // ISS.BKG.H1: L1电荷分布 (3×3×S个直方图)
    m_templates["ISS.BKG.H1"] = {
        .template_id = "ISS.BKG.H1",
        .name_pattern = "ISS_BKG_H1_{source}_{charge_type}_{detector}",
        .title_pattern = "L1 charge for {source} {charge_type} vs E_k/n in {detector}",
        .type = "TH2F",
        .binning_type = "charge_ek_2d",
        .multipliers = {"detector", "charge_type", "source"},
        .default_values = {
            {"detector", {"TOF", "NaF", "AGL"}},
            {"charge_type", {"Signal", "Template", "L2Template"}},
            {"source", {"Carbon", "Nitrogen", "Oxygen"}}
        }
    };
    
    // ISS.BKG.H2: 源粒子事件计数 (3×S个直方图)
    m_templates["ISS.BKG.H2"] = {
        .template_id = "ISS.BKG.H2",
        .name_pattern = "ISS_BKG_H2_{source}_{detector}",
        .title_pattern = "Source {source} event counts vs E_k/n in {detector}",
        .type = "TH1F",
        .binning_type = "ek_per_nucleon",
        .multipliers = {"detector", "source"},
        .default_values = {
            {"detector", {"TOF", "NaF", "AGL"}},
            {"source", {"Carbon", "Nitrogen", "Oxygen"}}
        }
    };
    
    // ISS.BKG.H3: L2碎片产物计数 (3×S个直方图)
    m_templates["ISS.BKG.H3"] = {
        .template_id = "ISS.BKG.H3",
        .name_pattern = "ISS_BKG_H3_{source}_{detector}",
        .title_pattern = "L2 frag products from {source} vs E_k/n in {detector}",
        .type = "TH1F",
        .binning_type = "ek_per_nucleon",
        .multipliers = {"detector", "source"},
        .default_values = {
            {"detector", {"TOF", "NaF", "AGL"}},
            {"source", {"Carbon", "Nitrogen", "Oxygen"}}
        }
    };
    
    // ISS.BKG.H4: 1/M of L2 fragmentation products (3×S个直方图)
    m_templates["ISS.BKG.H4"] = {
        .template_id = "ISS.BKG.H4",
        .name_pattern = "ISS_BKG_H4_{source}_{detector}",
        .title_pattern = "1/M of L2 frag products from {source} vs E_k/n in {detector}",
        .type = "TH2F",
        .binning_type = "invmass_ek_2d",
        .multipliers = {"detector", "source"},
        .default_values = {
            {"detector", {"TOF", "NaF", "AGL"}},
            {"source", {"Carbon", "Nitrogen", "Oxygen"}}
        }
    };
    
    Logger::Info("Registered {} product templates", m_templates.size());
}

std::vector<HistogramBlueprint> ProductRegistry::ExpandTemplate(
    const std::string& template_id, 
    const AnalysisContext* context) const {
    
    auto it = m_templates.find(template_id);
    if (it == m_templates.end()) {
        throw std::runtime_error("Template not found: " + template_id);
    }
    
    const auto& tmpl = it->second;
    std::vector<HistogramBlueprint> blueprints;
    
    // 获取源粒子列表（优先从YAML，否则使用默认值）
    std::vector<std::string> sources = GetSourceList(context);
    std::vector<std::string> detectors = tmpl.default_values.at("detector");
    
    if (tmpl.multipliers.size() == 2) {
        // 2D展开：detector × source (如ISS.BKG.H2, H3, H4)
        for (const auto& detector : detectors) {
            for (const auto& source : sources) {
                HistogramBlueprint blueprint;
                blueprint.product_id = tmpl.template_id;
                blueprint.name = ReplacePatterns(tmpl.name_pattern, {
                    {"{detector}", detector},
                    {"{source}", source}
                });
                blueprint.title = ReplacePatterns(tmpl.title_pattern, {
                    {"{detector}", detector},
                    {"{source}", source}
                });
                blueprint.type = tmpl.type;
                blueprint.binning_type = tmpl.binning_type;
                blueprint.detector_name = detector;
                blueprint.source_name = source;
                
                blueprints.push_back(blueprint);
            }
        }
    } else if (tmpl.multipliers.size() == 3) {
        // 3D展开：detector × charge_type × source (如ISS.BKG.H1)
        std::vector<std::string> charge_types = tmpl.default_values.at("charge_type");
        
        for (const auto& detector : detectors) {
            for (const auto& charge_type : charge_types) {
                for (const auto& source : sources) {
                    HistogramBlueprint blueprint;
                    blueprint.product_id = tmpl.template_id;
                    blueprint.name = ReplacePatterns(tmpl.name_pattern, {
                        {"{detector}", detector},
                        {"{charge_type}", charge_type},
                        {"{source}", source}
                    });
                    blueprint.title = ReplacePatterns(tmpl.title_pattern, {
                        {"{detector}", detector},
                        {"{charge_type}", charge_type},
                        {"{source}", source}
                    });
                    blueprint.type = tmpl.type;
                    blueprint.binning_type = tmpl.binning_type;
                    blueprint.detector_name = detector;
                    blueprint.source_name = source;
                    blueprint.charge_type = charge_type;
                    
                    blueprints.push_back(blueprint);
                }
            }
        }
    }
    
    Logger::Debug("Expanded template {} into {} products", template_id, blueprints.size());
    return blueprints;
}

std::vector<HistogramBlueprint> ProductRegistry::GetActiveProducts(
    const std::vector<std::string>& active_templates,
    const AnalysisContext* context) const {
    
    std::vector<HistogramBlueprint> all_products;
    
    for (const auto& template_id : active_templates) {
        try {
            auto products = ExpandTemplate(template_id, context);
            all_products.insert(all_products.end(), products.begin(), products.end());
        } catch (const std::exception& e) {
            Logger::Error("Failed to expand template {}: {}", template_id, e.what());
        }
    }
    
    Logger::Info("Generated {} total products from {} templates", 
                 all_products.size(), active_templates.size());
    return all_products;
}

std::string ProductRegistry::ReplacePatterns(const std::string& pattern,
                                            const std::map<std::string, std::string>& replacements) const {
    std::string result = pattern;
    for (const auto& [placeholder, value] : replacements) {
        size_t pos = 0;
        while ((pos = result.find(placeholder, pos)) != std::string::npos) {
            result.replace(pos, placeholder.length(), value);
            pos += value.length();
        }
    }
    return result;
}

std::vector<std::string> ProductRegistry::GetSourceList(const AnalysisContext* context) const {
    // 尝试从YAML配置获取背景源列表
    try {
        const auto& runSettings = context->GetRunSettings();
        if (runSettings["template_config"] && runSettings["template_config"]["background_sources"]) {
            return runSettings["template_config"]["background_sources"].as<std::vector<std::string>>();
        }
    } catch (const std::exception& e) {
        Logger::Debug("Could not load background_sources from config, using defaults: {}", e.what());
    }
    
    // 返回默认值
    return {"Carbon", "Nitrogen", "Oxygen"};
}

} // namespace IsoToolbox