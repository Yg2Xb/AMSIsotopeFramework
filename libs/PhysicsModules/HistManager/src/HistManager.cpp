#include "HistManager/HistManager.h"
#include "IsoToolbox/Logger.h"
#include <stdexcept>

namespace IsoToolbox {

HistManager::HistManager() {}
HistManager::~HistManager() {}

void HistManager::BookHistograms(const AnalysisContext* context, const BinningManager* binningManager) {
    Logger::Info("Booking histograms using ProductRegistry...");
    
    auto& registry = ProductRegistry::GetInstance();
    
    // 从AnalysisContext获取激活的模板
    std::vector<std::string> active_templates;
    try {
        const auto& runSettings = context->GetRunSettings();
        if (runSettings["active_templates"]) {
            active_templates = runSettings["active_templates"].as<std::vector<std::string>>();
        } else {
            Logger::Warn("No 'active_templates' found in config, using defaults");
            active_templates = {"ISS.BKG.H2", "ISS.BKG.H3"}; // 默认激活简单的1D直方图
        }
    } catch (const std::exception& e) {
        Logger::Error("Failed to load active templates: {}", e.what());
        return;
    }
    
    // 获取所有激活的产品
    auto blueprints = registry.GetActiveProducts(active_templates, context);
    
    // 为每个产品创建直方图
    for (const auto& blueprint : blueprints) {
        std::vector<double> x_bins, y_bins;
        
        // 根据binning_type获取合适的分bin
        if (blueprint.binning_type == "ek_per_nucleon") {
            // 1D直方图：使用默认的窄带Ek/n分bin
            y_bins = binningManager->GetBinning("ek_per_nucleon");
            
        } else if (blueprint.binning_type == "charge_ek_2d") {
            // 2D直方图：x轴是电荷，y轴是ek/n
            x_bins = binningManager->GetChargeBinning(blueprint.source_name);
            y_bins = binningManager->GetBinning("ek_per_nucleon");
            
        } else if (blueprint.binning_type == "invmass_ek_2d") {
            // 2D直方图：x轴是1/M，y轴是ek/n
            x_bins = binningManager->GetBinning("invmass");
            y_bins = binningManager->GetBinning("ek_per_nucleon");
        }
        
        CreateHistogram(blueprint, x_bins, y_bins);
    }
    
    Logger::Info("Booked {} histograms", m_hists.size());
}

void HistManager::CreateHistogram(const HistogramBlueprint& blueprint, 
                                 const std::vector<double>& x_bins,
                                 const std::vector<double>& y_bins) {
    
    // 🎯 直方图创建函数映射表（只支持TH1F和TH2F）
    static const std::map<std::string, std::function<TH1*(const HistogramBlueprint&, 
                                                          const std::vector<double>&, 
                                                          const std::vector<double>&)>> creators = {
        {"TH1F", [](const auto& bp, const auto& x_bins, const auto& y_bins) -> TH1* {
            if (y_bins.empty()) {
                Logger::Error("No bins provided for 1D histogram: {}", bp.name);
                return nullptr;
            }
            return new TH1F(bp.name.c_str(), bp.title.c_str(), 
                           y_bins.size() - 1, y_bins.data());
        }},
        
        {"TH2F", [](const auto& bp, const auto& x_bins, const auto& y_bins) -> TH1* {
            if (x_bins.empty() || y_bins.empty()) {
                Logger::Error("Insufficient bins for 2D histogram: {}", bp.name);
                return nullptr;
            }
            return new TH2F(bp.name.c_str(), bp.title.c_str(), 
                           x_bins.size() - 1, x_bins.data(),
                           y_bins.size() - 1, y_bins.data());
        }}
    };
    
    // 🎯 查找并创建
    auto creator_it = creators.find(blueprint.type);
    if (creator_it == creators.end()) {
        Logger::Error("Unsupported histogram type: {} (only TH1F and TH2F supported)", blueprint.type);
        return;
    }
    
    TH1* hist = creator_it->second(blueprint, x_bins, y_bins);
    if (!hist) return;
    
    // 🎯 统一设置轴标题
    hist->GetXaxis()->SetTitle(blueprint.x_title.c_str());
    hist->GetYaxis()->SetTitle(blueprint.y_title.c_str());
    
    m_hists[blueprint.name] = std::unique_ptr<TH1>(hist);
    Logger::Debug("Created {} histogram: {} with axes ['{}', '{}']", 
                  blueprint.type, blueprint.name, blueprint.x_title, blueprint.y_title);
}

void HistManager::Fill1D(const std::string& name, double value, double weight) {
    auto it = m_hists.find(name);
    if (it != m_hists.end()) {
        it->second->Fill(value, weight);
    } else {
        Logger::Warn("Histogram '{}' not found for Fill1D.", name);
    }
}

void HistManager::Fill2D(const std::string& name, double x_val, double y_val, double weight) {
    auto it = m_hists.find(name);
    if (it != m_hists.end()) {
        auto* hist2d = dynamic_cast<TH2*>(it->second.get());
        if (hist2d) {
            hist2d->Fill(x_val, y_val, weight);
        } else {
            Logger::Error("Histogram '{}' is not a 2D histogram.", name);
        }
    } else {
        Logger::Warn("Histogram '{}' not found for Fill2D.", name);
    }
}

void HistManager::SaveHistograms(const std::string& output_path) {
    TFile* file = TFile::Open(output_path.c_str(), "RECREATE");
    if (!file || file->IsZombie()) {
        throw std::runtime_error("Failed to open output file: " + output_path);
    }
    
    for (const auto& pair : m_hists) {
        pair.second->Write();
    }
    
    file->Close();
    delete file;
    Logger::Info("Histograms saved to {}", output_path);
}

} // namespace IsoToolbox