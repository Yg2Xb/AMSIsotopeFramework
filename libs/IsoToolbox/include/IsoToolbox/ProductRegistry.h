// libs/IsoToolbox/include/IsoToolbox/ProductRegistry.h
#pragma once

#include "IsoToolbox/AnalysisContext.h"
#include <vector>
#include <string>
#include <map>
#include <memory>

namespace IsoToolbox {

struct ProductTemplate {
    std::string template_id;           // "ISS.BKG.H1"
    std::string name_pattern;          // "ISS_BKG_H1_{source}_{charge_type}_{detector}"
    std::string title_pattern;         // "L1 charge for {source} {charge_type} vs E_{k}/n in {detector}"  // 🎯 更新注释
    std::string x_title;               // "{detector} E_{k}/n [GeV/nucleon]"  // 🎯 更新注释示例
    std::string y_title;               // "{detector} Charge [a.u.]"          // 🎯 更新注释示例
    std::string type;                  // "TH1F", "TH2F"
    std::string binning_type;          // "ek_per_nucleon", "charge_ek_2d"
    
    // 扩展维度
    std::vector<std::string> multipliers;  // ["detector", "charge_type", "source"]
    std::map<std::string, std::vector<std::string>> default_values;
};

struct HistogramBlueprint {
    std::string product_id;      // "ISS.BKG.H1"
    std::string name;           // "ISS_BKG_H1_Carbon_Signal_TOF"
    std::string title;          // "L1 charge for Carbon Signal vs E_{k}/n in TOF"  // 🎯 更新注释
    std::string x_title;        // "TOF E_{k}/n [GeV/nucleon]"     // 🎯 更新注释示例
    std::string y_title;        // "TOF Charge [a.u.]"             // 🎯 更新注释示例
    std::string type;           // "TH1F"
    std::string binning_type;   // "ek_per_nucleon"
    
    // 用于BinningManager的参数
    std::string isotope_name;   // "Be10" (如果需要)
    std::string detector_name;  // "TOF"
    std::string source_name;    // "Carbon"
    std::string charge_type;    // "Signal"
    
    // 自定义bins（如果有YAML覆盖）
    std::vector<double> custom_x_bins;
    std::vector<double> custom_y_bins;
};

class ProductRegistry {
public:
    static ProductRegistry& GetInstance();
    
    void RegisterCoreTemplates();
    
    std::vector<HistogramBlueprint> ExpandTemplate(
        const std::string& template_id, 
        const AnalysisContext* context) const;
    
    // 获取激活的产品列表
    std::vector<HistogramBlueprint> GetActiveProducts(
        const std::vector<std::string>& active_templates,
        const AnalysisContext* context) const;

private:
    ProductRegistry() = default;
    std::map<std::string, ProductTemplate> m_templates;
    
    std::string ReplacePatterns(const std::string& pattern,
                               const std::map<std::string, std::string>& replacements) const;
                               
    std::vector<std::string> GetSourceList(const AnalysisContext* context) const;
    
    std::vector<std::string> GetIsotopeList(const AnalysisContext* context) const;

};

} // namespace IsoToolbox