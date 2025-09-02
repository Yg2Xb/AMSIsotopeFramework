// libs/IsoToolbox/src/ProductRegistry.cpp
#include "IsoToolbox/ProductRegistry.h"
#include "IsoToolbox/Logger.h"
#include <stdexcept>
#include <algorithm>
#include <functional>

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
    
    // === ISS Background Templates ===
    // ISS.BKG.H1: L1电荷分布 (3×3×S个直方图)
    m_templates["ISS.BKG.H1"] = {
        .template_id = "ISS.BKG.H1",
        .name_pattern = "ISS_BKG_H1_{source}_{charge_type}_{detector}",
        .title_pattern = "L1 charge for {source} {charge_type} vs E_{k}/n in {detector}",
        .x_title = "{detector} E_{k}/n [GeV/nucleon]",
        .y_title = "{detector} Charge [a.u.]",
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
        .title_pattern = "Source {source} event counts vs E_{k}/n in {detector}",
        .x_title = "{detector} E_{k}/n [GeV/nucleon]",
        .y_title = "Counts",
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
        .title_pattern = "L2 frag products from {source} vs E_{k}/n in {detector}",
        .x_title = "{detector} E_{k}/n [GeV/nucleon]",
        .y_title = "Counts",
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
        .title_pattern = "1/M of L2 frag products from {source} vs E_{k}/n in {detector}",
        .x_title = "{detector} E_{k}/n [GeV/nucleon]",
        .y_title = "{detector} 1/M [nucleon/GeV]",
        .type = "TH2F",
        .binning_type = "invmass_ek_2d",
        .multipliers = {"detector", "source"},
        .default_values = {
            {"detector", {"TOF", "NaF", "AGL"}},
            {"source", {"Carbon", "Nitrogen", "Oxygen"}}
        }
    };
    
    // === MC Background Templates ===
    
    // MC.BKG.H1: L1 source nuclei counts vs. E_k/n (3个探测器)
    m_templates["MC.BKG.H1"] = {
        .template_id = "MC.BKG.H1",
        .name_pattern = "MC_BKG_H1_{detector}",
        .title_pattern = "L1 source nuclei counts vs E_{k}/n in {detector}",
        .x_title = "{detector} E_{k}/n [GeV/nucleon]",
        .y_title = "Counts",
        .type = "TH1F",
        .binning_type = "ek_per_nucleon",
        .multipliers = {"detector"},
        .default_values = {
            {"detector", {"TOF", "NaF", "AGL"}}
        }
    };
    
    // MC.BKG.H2: L2 fragmentation into target isotopes (3×N个)
    m_templates["MC.BKG.H2"] = {
        .template_id = "MC.BKG.H2",
        .name_pattern = "MC_BKG_H2_{detector}_{isotope}",
        .title_pattern = "L2 fragmentation into {isotope} vs E_{k}/n in {detector}",
        .x_title = "{detector} E_{k}/n [GeV/nucleon]",
        .y_title = "Counts",
        .type = "TH1F",
        .binning_type = "ek_per_nucleon",
        .multipliers = {"detector", "isotope"},
        .default_values = {
            {"detector", {"TOF", "NaF", "AGL"}},
            {"isotope", {"Be7", "Be9", "Be10"}}  // 从配置读取
        }
    };
    
    // MC.BKG.H3a: upTOF fragmented isotope counts (3×N个)
    m_templates["MC.BKG.H3a"] = {
        .template_id = "MC.BKG.H3a",
        .name_pattern = "MC_BKG_H3a_{detector}_{isotope}",
        .title_pattern = "upTOF fragmented {isotope} counts vs E_{k}/n in {detector}",
        .x_title = "{detector} E_{k}/n [GeV/nucleon]",
        .y_title = "Counts",
        .type = "TH1F",
        .binning_type = "ek_per_nucleon",
        .multipliers = {"detector", "isotope"},
        .default_values = {
            {"detector", {"TOF", "NaF", "AGL"}},
            {"isotope", {"Be7", "Be9", "Be10"}}
        }
    };
    
    // MC.BKG.H3b: upTOF fragmented and RICH survived (3×N个)
    m_templates["MC.BKG.H3b"] = {
        .template_id = "MC.BKG.H3b",
        .name_pattern = "MC_BKG_H3b_{detector}_{isotope}",
        .title_pattern = "upTOF fragmented & RICH survived {isotope} vs E_{k}/n in {detector}",
        .x_title = "{detector} E_{k}/n [GeV/nucleon]",
        .y_title = "Counts",
        .type = "TH1F",
        .binning_type = "ek_per_nucleon",
        .multipliers = {"detector", "isotope"},
        .default_values = {
            {"detector", {"TOF", "NaF", "AGL"}},
            {"isotope", {"Be7", "Be9", "Be10"}}
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
    
    // 构建所有维度的值列表
    std::vector<std::vector<std::string>> dimension_values;
    std::vector<std::string> dimension_names;
    
    for (const auto& multiplier : tmpl.multipliers) {
        dimension_names.push_back(multiplier);
        
        if (multiplier == "source") {
            dimension_values.push_back(GetSourceList(context));
        } else if (multiplier == "isotope") {
            dimension_values.push_back(GetIsotopeList(context));  // 🎯 新增
        } else {
            dimension_values.push_back(tmpl.default_values.at(multiplier));
        }
    }
    
    // 递归生成所有组合
    std::function<void(int, std::map<std::string, std::string>&)> generateCombinations = 
        [&](int dim_index, std::map<std::string, std::string>& current_replacements) {
            
            if (dim_index == dimension_names.size()) {
                // 到达叶子节点，创建blueprint
                HistogramBlueprint blueprint;
                blueprint.product_id = tmpl.template_id;
                blueprint.name = ReplacePatterns(tmpl.name_pattern, current_replacements);
                blueprint.title = ReplacePatterns(tmpl.title_pattern, current_replacements);
                blueprint.x_title = ReplacePatterns(tmpl.x_title, current_replacements);
                blueprint.y_title = ReplacePatterns(tmpl.y_title, current_replacements);
                blueprint.type = tmpl.type;
                blueprint.binning_type = tmpl.binning_type;
                
                // 设置特定字段
                if (current_replacements.find("{detector}") != current_replacements.end()) {
                    blueprint.detector_name = current_replacements["{detector}"];
                }
                if (current_replacements.find("{source}") != current_replacements.end()) {
                    blueprint.source_name = current_replacements["{source}"];
                }
                if (current_replacements.find("{charge_type}") != current_replacements.end()) {
                    blueprint.charge_type = current_replacements["{charge_type}"];
                }
                
                blueprints.push_back(blueprint);
                return;
            }
            
            // 递归处理当前维度
            const std::string& dim_name = dimension_names[dim_index];
            const auto& values = dimension_values[dim_index];
            
            for (const auto& value : values) {
                current_replacements["{" + dim_name + "}"] = value;
                generateCombinations(dim_index + 1, current_replacements);
            }
        };
    
    std::map<std::string, std::string> replacements;
    generateCombinations(0, replacements);
    
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
    Logger::Debug("ProductRegistry::GetSourceList() called");
    
    try {
        const auto& runSettings = context->GetRunSettings();
        Logger::Debug("runSettings size: {}", runSettings.size());
        
        if (runSettings["template_config"]) {
            Logger::Debug("template_config found");
            if (runSettings["template_config"]["background_sources"]) {
                auto sources = runSettings["template_config"]["background_sources"].as<std::vector<std::string>>();
                Logger::Debug("✅ Successfully read background_sources: [{}]", fmt::join(sources, ", "));
                return sources;
            } else {
                Logger::Warn("template_config.background_sources NOT found");
            }
        } else {
            Logger::Warn("template_config section NOT found");
        }
    } catch (const std::exception& e) {
        Logger::Error("Exception reading background_sources: {}", e.what());
    }
    
    Logger::Warn("Using hardcoded defaults");
    return {"Carbon", "Nitrogen", "Oxygen"};
}
    
    std::vector<std::string> ProductRegistry::GetIsotopeList(const AnalysisContext* context) const {
        try {
            const auto& runSettings = context->GetRunSettings();
            if (runSettings["template_config"]["target_isotopes"]) {
                auto isotopes = runSettings["template_config"]["target_isotopes"].as<std::vector<std::string>>();
                Logger::Debug("✅ Successfully read target_isotopes: [{}]", fmt::join(isotopes, ", "));
                return isotopes;
            }
        } catch (const std::exception& e) {
            Logger::Error("Exception reading target_isotopes: {}", e.what());
        }
        
        Logger::Warn("Using hardcoded isotope defaults");
        return {"Be7", "Be9", "Be10"};
    }

} // namespace IsoToolbox