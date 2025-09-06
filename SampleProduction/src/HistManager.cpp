#include "HistManager.h"
#include "ProductRegistry.h" // 包含完整的蓝图定义
#include "BinningManager.h"
#include <iostream>
#include <stdexcept>

namespace AMS_Iso {

HistManager::HistManager(const std::string& output_filename) {
    m_outputFile = std::make_unique<TFile>(output_filename.c_str(), "RECREATE");
    if (!m_outputFile || m_outputFile->IsZombie()) {
        throw std::runtime_error("Failed to create output ROOT file: " + output_filename);
    }
    std::cout << "HistManager: Output file '" << output_filename << "' opened." << std::endl;
}

HistManager::~HistManager() {
    if (m_outputFile && m_outputFile->IsOpen()) {
        m_outputFile->Close();
    }
}

void HistManager::CreateProducts(const std::vector<HistogramBlueprint>& blueprints) {
    m_outputFile->cd();
    const auto& binningManager = BinningManager::GetInstance();

    for (const auto& bp : blueprints) {
        if (m_histograms.count(bp.name)) {
            std::cerr << "WARNING: Duplicate histogram name detected, skipping creation: " << bp.name << std::endl;
            continue;
        }

        TH1* hist = nullptr;
        try {
            if (bp.type == "TH1F") {
                const auto& bins = binningManager.GetBinning(bp.binning_x_key);
                // MODIFIED: 如果binning数组为空，则使用默认分bin，增强鲁棒性
                if (bins.empty()) {
                    hist = new TH1F(bp.name.c_str(), bp.title.c_str(), 100, 0, 1);
                     std::cerr << "WARNING: Binning key '" << bp.binning_x_key << "' for TH1F '" << bp.name << "' is empty. Using default binning." << std::endl;
                } else {
                    hist = new TH1F(bp.name.c_str(), bp.title.c_str(), bins.size() - 1, bins.data());
                }
            } else if (bp.type == "TH2F") {
                const auto& xbins = binningManager.GetBinning(bp.binning_x_key);
                const auto& ybins = binningManager.GetBinning(bp.binning_y_key);
                // MODIFIED: 如果任一binning数组为空，则使用默认分bin
                if (xbins.empty() || ybins.empty()) {
                     hist = new TH2F(bp.name.c_str(), bp.title.c_str(), 100, 0, 1, 100, 0, 1);
                     std::cerr << "WARNING: Binning key '" << bp.binning_x_key << "' or '" << bp.binning_y_key << "' for TH2F '" << bp.name << "' is empty. Using default binning." << std::endl;
                } else {
                    hist = new TH2F(bp.name.c_str(), bp.title.c_str(),
                                    xbins.size() - 1, xbins.data(),
                                    ybins.size() - 1, ybins.data());
                }
            } else {
                std::cerr << "WARNING: Unknown histogram type '" << bp.type << "' for '" << bp.name << "'" << std::endl;
                continue;
            }
        } catch (const std::runtime_error& e) {
            std::cerr << "ERROR: Failed to create histogram '" << bp.name << "'. Reason: " << e.what() << std::endl;
            continue;
        }
        
        hist->SetXTitle(bp.x_axis_title.c_str());
        hist->SetYTitle(bp.y_axis_title.c_str());
        hist->Sumw2();
        
        m_histograms[bp.name] = hist;
    }
    std::cout << "HistManager created " << m_histograms.size() << " histograms." << std::endl;
}

void HistManager::Save() {
    if (m_outputFile && m_outputFile->IsOpen()) {
        m_outputFile->cd();
        for (auto const& [name, hist] : m_histograms) {
            hist->Write();
        }
        std::cout << "HistManager saved " << m_histograms.size() << " histograms to " << m_outputFile->GetName() << std::endl;
    }
}

TH1* HistManager::GetHistogram(const std::string& name) {
    auto it = m_histograms.find(name);
    if (it != m_histograms.end()) {
        return it->second;
    }
    return nullptr;
}

TH1F* HistManager::GetHist1F(const std::string& name) {
    return dynamic_cast<TH1F*>(GetHistogram(name));
}

TH2F* HistManager::GetHist2F(const std::string& name) {
    return dynamic_cast<TH2F*>(GetHistogram(name));
}

} // namespace AMS_Iso

