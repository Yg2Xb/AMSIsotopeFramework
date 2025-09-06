#include "ProductRegistry.h"
#include "basic_var.h"
#include <iostream>
#include <functional>
#include <stdexcept>
#include <algorithm> // for std::replace

namespace AMS_Iso {

ProductRegistry::ProductRegistry()
    : m_analyzer(nullptr), m_isRegistered(false) {}

ProductRegistry& ProductRegistry::GetInstance() {
    static ProductRegistry instance;
    return instance;
}

void ProductRegistry::SetAnalyzerInfo(const IsotopeAnalyzer* analyzer) {
    m_analyzer = analyzer;
}

const std::vector<HistogramBlueprint>& ProductRegistry::GetBlueprints() const {
    return m_blueprints;
}

void ProductRegistry::RegisterProducts(bool isMC, const std::vector<std::string>& active_chains) {
    if (m_isRegistered) return;
    if (!m_analyzer) {
        throw std::runtime_error("ProductRegistry Error: IsotopeAnalyzer info was not set before registering products.");
    }

    if (isMC) {
        RegisterCoreTemplates_MC();
    } else {
        RegisterCoreTemplates_ISS();
    }
    
    ExpandTemplates(active_chains);
    m_isRegistered = true;
}

// MODIFIED: 为需要同位素专属分bin的模板，在binning_x_key中加入{isotope}占位符
void ProductRegistry::RegisterCoreTemplates_ISS() {
    m_templates.clear();
    const std::vector<std::string> sources = {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
    const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
    const std::vector<std::string> charge_types = {"L1QSignal", "L1QTemplate", "L2QTemplate"};
    const std::vector<std::string> cut_groups = {"CutGroup1", "CutGroup2", "CutGroup3"};
    const std::vector<std::string> num_den = {"Num", "Den"};

    m_templates["ISS.ID.H1"] = {"ISS.ID.H1", "{chain}: Event counts of {isotope} vs E_k/n in {detector}", "{chain}_{id}_{detector}_{isotope}", "TH1F", "EkPerNucleon_{isotope}", "", "{detector} E_{k}/n [GeV/nucleon]", "Counts", {"chain", "detector", "isotope"}, {{"detector", detectors}}};
    m_templates["ISS.ID.H2"] = {"ISS.ID.H2", "{chain}: 1/M of {isotope} vs E_k/n in {detector}", "{chain}_{id}_{detector}_{isotope}", "TH2F", "EkPerNucleon_{isotope}", "InverseMass", "{detector} E_{k}/n [GeV/nucleon]", "1/M [nucleon/GeV]", {"chain", "detector", "isotope"}, {{"detector", detectors}}};
    m_templates["ISS.ID.H3"] = {"ISS.ID.H3", "{chain}: Heaviest Isotope 1/M vs E_k/n (Cutoff) in {detector}", "{chain}_{id}_{detector}", "TH2F", "EkPerNucleon", "InverseMass", "{detector} E_{k}/n [GeV/nucleon]", "1/M [nucleon/GeV]", {"chain", "detector"}, {{"detector", detectors}}};
    m_templates["ISS.ID.H4a"] = {"ISS.ID.H4a", "{chain}: 1/beta RICH-NaF (beta~1)", "{chain}_{id}", "TH1F", "Beta", "", "1/#beta", "Counts", {"chain"}, {}};
    m_templates["ISS.ID.H4b"] = {"ISS.ID.H4b", "{chain}: 1/beta RICH-AGL (beta~1)", "{chain}_{id}", "TH1F", "Beta", "", "1/#beta", "Counts", {"chain"}, {}};
    m_templates["ISS.ID.H5a"] = {"ISS.ID.H5a", "{chain}: d(1/beta) NaF-Trk vs Rigidity", "{chain}_{id}", "TH2F", "Rigidity", "DeltaBeta", "Rigidity (GV)", "1/#beta_{NaF} - 1/#beta_{Trk}", {"chain"}, {}};
    m_templates["ISS.ID.H5b"] = {"ISS.ID.H5b", "{chain}: d(1/beta) AGL-Trk vs Rigidity", "{chain}_{id}", "TH2F", "Rigidity", "DeltaBeta", "Rigidity (GV)", "1/#beta_{AGL} - 1/#beta_{Trk}", {"chain"}, {}};
    m_templates["ISS.ID.H6a"] = {"ISS.ID.H6a", "{chain}: d(1/beta) TOF-NaF vs E_k/n", "{chain}_{id}", "TH2F", "EkPerNucleon", "DeltaBeta", "E_k/n (GeV/n)", "1/#beta_{TOF} - 1/#beta_{NaF}", {"chain"}, {}};
    m_templates["ISS.ID.H6b"] = {"ISS.ID.H6b", "{chain}: d(1/beta) TOF-AGL vs E_k/n", "{chain}_{id}", "TH2F", "EkPerNucleon", "DeltaBeta", "E_k/n (GeV/n)", "1/#beta_{TOF} - 1/#beta_{AGL}", {"chain"}, {}};
    m_templates["ISS.ID.H7a"] = {"ISS.ID.H7a", "{chain}: d(1/beta) TOF-NaF vs betagamma", "{chain}_{id}", "TH2F", "BetaGamma", "DeltaBeta", "#beta#gamma", "1/#beta_{TOF} - 1/#beta_{NaF}", {"chain"}, {}};
    m_templates["ISS.ID.H7b"] = {"ISS.ID.H7b", "{chain}: d(1/beta) TOF-AGL vs betagamma", "{chain}_{id}", "TH2F", "BetaGamma", "DeltaBeta", "#beta#gamma", "1/#beta_{TOF} - 1/#beta_{AGL}", {"chain"}, {}};
    m_templates["ISS.BKG.H1"] = {"ISS.BKG.H1", "{chain}: L1 charge for {source} {charge_type} vs E_k/n in {detector}", "{chain}_{id}_{source}_{charge_type}_{detector}", "TH2F", "EkPerNucleon", "Charge", "{detector} E_{k}/n [GeV/nucleon]", "{detector} Charge", {"chain", "detector", "charge_type", "source"}, {{"detector", detectors}, {"charge_type", charge_types}, {"source", sources}}};
    m_templates["ISS.BKG.H2"] = {"ISS.BKG.H2", "{chain}: Source {source} event counts (L1Q) vs E_k/n in {detector}", "{chain}_{id}_{source}_{detector}", "TH1F", "EkPerNucleon", "", "{detector} E_{k}/n [GeV/nucleon]", "Counts", {"chain", "detector", "source"}, {{"detector", detectors}, {"source", sources}}};
    m_templates["ISS.BKG.H3"] = {"ISS.BKG.H3", "{chain}: L2 frag from {source} (InnerQ) vs E_k/n in {detector}", "{chain}_{id}_{source}_{detector}", "TH1F", "EkPerNucleon", "", "{detector} E_{k}/n [GeV/nucleon]", "Counts", {"chain", "detector", "source"}, {{"detector", detectors}, {"source", sources}}};
    m_templates["ISS.BKG.H4"] = {"ISS.BKG.H4", "{chain}: 1/M of L2 frag from {source} vs E_k/n in {detector}", "{chain}_{id}_{source}_{detector}", "TH2F", "EkPerNucleon", "InverseMass", "{detector} E_{k}/n [GeV/nucleon]", "1/M [nucleon/GeV]", {"chain", "detector", "source"}, {{"detector", detectors}, {"source", sources}}};
    m_templates["ISS.FLUX.H1"] = {"ISS.FLUX.H1", "{chain}: Exposure time vs Rigidity", "{chain}_{id}", "TH1F", "Rigidity", "", "Rigidity [GV]", "Exposure Time [s]", {"chain"}, {}};
    m_templates["ISS.FLUX.H2"] = {"ISS.FLUX.H2", "{chain}: Exposure time for {isotope} vs E_k/n in {detector}", "{chain}_{id}_{detector}_{isotope}", "TH1F", "EkPerNucleon_{isotope}", "", "{detector} E_{k}/n [GeV/nucleon]", "Exposure Time [s]", {"chain", "detector", "isotope"}, {{"detector", detectors}}};
    m_templates["ISS.FLUX.H3"] = {"ISS.FLUX.H3", "{chain}: CutEff {cut_group} {num_den} for {isotope} vs E_k/n in {detector}", "{chain}_{id}_{cut_group}_{num_den}_{detector}_{isotope}", "TH1F", "EkPerNucleon_{isotope}", "", "{detector} E_{k}/n [GeV/nucleon]", "Counts", {"chain", "detector", "isotope", "cut_group", "num_den"}, {{"detector", detectors}, {"cut_group", cut_groups}, {"num_den", num_den}}};
}

void ProductRegistry::RegisterCoreTemplates_MC() {
    m_templates.clear();
    const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
    const std::vector<std::string> cut_groups = {"CutGroup1", "CutGroup2", "CutGroup3"};
    const std::vector<std::string> num_den = {"Num", "Den"};
    const std::vector<std::string> gen_rec = {"gen", "rec"};

    m_templates["MC.ID.H1a"] = {"MC.ID.H1a", "{chain}: Mass Template (Sel) 1/M of {isotope} vs E_k/n in {detector}", "{chain}_{id}_{detector}_{isotope}", "TH2F", "EkPerNucleon_{isotope}", "InverseMass", "{detector} E_{k}/n [GeV/nucleon]", "1/M [nucleon/GeV]", {"chain", "detector", "isotope"}, {{"detector", detectors}}};
    m_templates["MC.ID.H1b"] = {"MC.ID.H1b", "{chain}: Mass Template (Sel+NoTofFrag) 1/M of {isotope} vs E_k/n in {detector}", "{chain}_{id}_{detector}_{isotope}", "TH2F", "EkPerNucleon_{isotope}", "InverseMass", "{detector} E_{k}/n [GeV/nucleon]", "1/M [nucleon/GeV]", {"chain", "detector", "isotope"}, {{"detector", detectors}}};
    m_templates["MC.ID.H2"] = {"MC.ID.H2", "{chain}: Heaviest Isotope 1/M vs E_k/n (Cutoff) in {detector}", "{chain}_{id}_{detector}", "TH2F", "EkPerNucleon", "InverseMass", "{detector} E_{k}/n [GeV/nucleon]", "1/M [nucleon/GeV]", {"chain", "detector"}, {{"detector", detectors}}};
    m_templates["MC.ID.H3a"] = {"MC.ID.H3a", "{chain}: 1/beta RICH-NaF (beta~1)", "{chain}_{id}", "TH1F", "Beta", "", "1/#beta", "Counts", {"chain"}, {}};
    m_templates["MC.ID.H3b"] = {"MC.ID.H3b", "{chain}: 1/beta RICH-AGL (beta~1)", "{chain}_{id}", "TH1F", "Beta", "", "1/#beta", "Counts", {"chain"}, {}};
    m_templates["MC.ID.H4a"] = {"MC.ID.H4a", "{chain}: d(1/beta) NaF-Trk vs Rigidity", "{chain}_{id}", "TH2F", "Rigidity", "DeltaBeta", "Rigidity (GV)", "1/#beta_{NaF} - 1/#beta_{Trk}", {"chain"}, {}};
    m_templates["MC.ID.H4b"] = {"MC.ID.H4b", "{chain}: d(1/beta) AGL-Trk vs Rigidity", "{chain}_{id}", "TH2F", "Rigidity", "DeltaBeta", "Rigidity (GV)", "1/#beta_{AGL} - 1/#beta_{Trk}", {"chain"}, {}};
    m_templates["MC.ID.H5a"] = {"MC.ID.H5a", "{chain}: d(1/beta) TOF-NaF vs E_k/n", "{chain}_{id}", "TH2F", "EkPerNucleon", "DeltaBeta", "E_k/n (GeV/n)", "1/#beta_{TOF} - 1/#beta_{NaF}", {"chain"}, {}};
    m_templates["MC.ID.H5b"] = {"MC.ID.H5b", "{chain}: d(1/beta) TOF-AGL vs E_k/n", "{chain}_{id}", "TH2F", "EkPerNucleon", "DeltaBeta", "E_k/n (GeV/n)", "1/#beta_{TOF} - 1/#beta_{AGL}", {"chain"}, {}};
    m_templates["MC.ID.H6a"] = {"MC.ID.H6a", "{chain}: d(1/beta) TOF-NaF vs betagamma", "{chain}_{id}", "TH2F", "BetaGamma", "DeltaBeta", "#beta#gamma", "1/#beta_{TOF} - 1/#beta_{NaF}", {"chain"}, {}};
    m_templates["MC.ID.H6b"] = {"MC.ID.H6b", "{chain}: d(1/beta) TOF-AGL vs betagamma", "{chain}_{id}", "TH2F", "BetaGamma", "DeltaBeta", "#beta#gamma", "1/#beta_{TOF} - 1/#beta_{AGL}", {"chain"}, {}};
    m_templates["MC.BKG.H1"] = {"MC.BKG.H1", "{chain}: L1 Source Counts vs E_k/n in {detector}", "{chain}_{id}_{detector}", "TH1F", "EkPerNucleon", "", "{detector} E_{k}/n [GeV/nucleon]", "Counts", {"chain", "detector"}, {{"detector", detectors}}};
    m_templates["MC.BKG.H2"] = {"MC.BKG.H2", "{chain}: L2 frag to {isotope} vs E_k/n in {detector}", "{chain}_{id}_{detector}_{isotope}", "TH1F", "EkPerNucleon_{isotope}", "", "{detector} E_{k}/n [GeV/nucleon]", "Counts", {"chain", "detector", "isotope"}, {{"detector", detectors}}};
    m_templates["MC.BKG.H3a"] = {"MC.BKG.H3a", "{chain}: upTOF frag {isotope} vs E_k/n in {detector}", "{chain}_{id}_{detector}_{isotope}", "TH1F", "EkPerNucleon_{isotope}", "", "{detector} E_{k}/n [GeV/nucleon]", "Counts", {"chain", "detector", "isotope"}, {{"detector", detectors}}};
    m_templates["MC.BKG.H3b"] = {"MC.BKG.H3b", "{chain}: upTOF frag {isotope} (RICH survived) vs E_k/n in {detector}", "{chain}_{id}_{detector}_{isotope}", "TH1F", "EkPerNucleon_{isotope}", "", "{detector} E_{k}/n [GeV/nucleon]", "Counts", {"chain", "detector", "isotope"}, {{"detector", detectors}}};
    m_templates["MC.FLUX.H1"] = {"MC.FLUX.H1", "{chain}: Accumulated Passed Events for {cut_group} vs E_k_{gen_rec}", "{chain}_{id}_{cut_group}_{gen_rec}", "TH1F", "EkGen", "", "E_k [GeV/nucleon]", "Counts", {"chain", "cut_group", "gen_rec"}, {{"cut_group", cut_groups}, {"gen_rec", gen_rec}}};
    m_templates["MC.FLUX.H2"] = {"MC.FLUX.H2", "{chain}: CutEff {cut_group} {num_den} for {isotope} vs E_k/n in {detector}", "{chain}_{id}_{cut_group}_{num_den}_{detector}_{isotope}", "TH1F", "EkPerNucleon_{isotope}", "", "{detector} E_{k}/n [GeV/nucleon]", "Counts", {"chain", "detector", "isotope", "cut_group", "num_den"}, {{"detector", detectors}, {"cut_group", cut_groups}, {"num_den", num_den}}};
    m_templates["MC.FLUX.H3"] = {"MC.FLUX.H3", "{chain}: Total Generated Events vs E_k,gen", "{chain}_{id}", "TH1F", "EkGen", "", "E_{k,gen} [GeV/nucleon]", "Counts", {"chain"}, {}};
}


void ProductRegistry::ExpandTemplates(const std::vector<std::string>& active_chains) {
    m_blueprints.clear();

    const auto* iso_info = m_analyzer ? m_analyzer->getIsotope() : nullptr;
    if (!iso_info) {
        throw std::runtime_error("Isotope information is not available in ProductRegistry.");
    }
    
    std::vector<std::string> isotope_list;
    if (m_analyzer->isISS()) {
        for (int i = 0; i < iso_info->getIsotopeCount(); ++i) {
            isotope_list.push_back(iso_info->getName() + std::to_string(iso_info->getMass(i)));
        }
    } else { 
        isotope_list.push_back(iso_info->getName() + std::to_string(m_analyzer->getUseMass()));
    }
    
    for (auto const& [key, tmpl] : m_templates) {
        std::vector<std::vector<std::string>> dimension_values;
        std::vector<std::string> dimension_names;

        for (const auto& multiplier : tmpl.multipliers) {
            dimension_names.push_back(multiplier);
            if (multiplier == "chain") {
                dimension_values.push_back(active_chains);
            } else if (multiplier == "isotope") {
                dimension_values.push_back(isotope_list);
            } else if (tmpl.default_values.count(multiplier)) {
                dimension_values.push_back(tmpl.default_values.at(multiplier));
            } else {
                 throw std::runtime_error("Multiplier '" + multiplier + "' in template '" + tmpl.template_id + "' has no default values");
            }
        }

        std::function<void(int, std::map<std::string, std::string>)> generate_combinations =
            [&](int dim_index, std::map<std::string, std::string> current_replacements) {
            if (dim_index == (int)dimension_names.size()) {
                HistogramBlueprint bp;
                
                std::string sanitized_id = tmpl.template_id;
                std::replace(sanitized_id.begin(), sanitized_id.end(), '.', '_');
                current_replacements["id"] = sanitized_id;
                
                bp.name = ReplacePatterns(tmpl.name_pattern, current_replacements);
                bp.title = ReplacePatterns(tmpl.title_pattern, current_replacements);
                bp.x_axis_title = ReplacePatterns(tmpl.x_title_pattern, current_replacements);
                bp.y_axis_title = ReplacePatterns(tmpl.y_title_pattern, current_replacements);
                bp.type = tmpl.type;

                // --- MODIFIED: 这里是关键修正 ---
                // 动态替换binning key中的占位符
                bp.binning_x_key = ReplacePatterns(tmpl.binning_x_key, current_replacements);
                bp.binning_y_key = ReplacePatterns(tmpl.binning_y_key, current_replacements);
                
                m_blueprints.push_back(bp);
                return;
            }

            const std::string& dim_name = dimension_names[dim_index];
            const auto& values = dimension_values[dim_index];
            for (const auto& value : values) {
                current_replacements[dim_name] = value;
                generate_combinations(dim_index + 1, current_replacements);
            }
        };

        std::map<std::string, std::string> initial_replacements;
        generate_combinations(0, initial_replacements);
    }
}

std::string ProductRegistry::ReplacePatterns(const std::string& pattern, const std::map<std::string, std::string>& replacements) {
    std::string result = pattern;
    for (const auto& [placeholder, value] : replacements) {
        std::string token = "{" + placeholder + "}";
        size_t pos = 0;
        while ((pos = result.find(token, pos)) != std::string::npos) {
            result.replace(pos, token.length(), value);
            pos += value.length();
        }
    }
    return result;
}

} // namespace AMS_Iso