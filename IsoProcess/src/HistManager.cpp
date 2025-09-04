// File: IsoProcess/src/HistManager.cpp
// Purpose: The final version of HistManager, supporting multiple analysis chains.

#include "HistManager.h"
#include "BinningManager.h" // Needs the full definition of BinningManager
#include "ProductRegistry.h"
#include <iostream>
#include <stdexcept>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

// We pass the fully constructed BinningManager from main.cpp
HistManager::HistManager(BinningManager* binner) : m_binner(binner) {
    if (!m_binner) {
        throw std::runtime_error("HistManager received a null BinningManager pointer!");
    }
    // Tell ROOT to let us manage memory for histograms
    TH1::AddDirectory(kFALSE);
}

HistManager::~HistManager() {}

void HistManager::BookHistograms(bool isMC,
                               const std::vector<std::string>& active_chains,
                               const std::vector<std::string>& active_groups) {
    std::cout << "Booking histograms for " << active_chains.size() << " chain(s)..." << std::endl;
    
    ProductRegistry::GetInstance().RegisterProducts(isMC, active_chains, active_groups);

    const auto& blueprints = ProductRegistry::GetInstance().GetBlueprints();

    for (const auto& bp : blueprints) {
        CreateHistogram(bp);
    }
    std::cout << "Successfully booked " << m_hists.size() << " total histograms." << std::endl;
}

void HistManager::CreateHistogram(const HistogramBlueprint& bp) {
    TH1* hist = nullptr;
    try {
        if (bp.type == "TH1F") {
            const auto& ybins = m_binner->Get(bp.binningY);
            hist = new TH1F(bp.uniqueID.c_str(), bp.title.c_str(), ybins.size() - 1, ybins.data());
        } else if (bp.type == "TH2F") {
            const auto& xbins = m_binner->Get(bp.binningX);
            const auto& ybins = m_binner->Get(bp.binningY);
            hist = new TH2F(bp.uniqueID.c_str(), bp.title.c_str(), xbins.size() - 1, xbins.data(), ybins.size() - 1, ybins.data());
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Error creating histogram '" << bp.uniqueID << "': " << e.what() << std::endl;
        return;
    }

    if (hist) {
        hist->GetXaxis()->SetTitle(bp.axisTitleX.c_str());
        hist->GetYaxis()->SetTitle(bp.axisTitleY.c_str());
        m_hists[bp.uniqueID] = std::unique_ptr<TH1>(hist);
    } else {
        std::cerr << "Warning: Could not create histogram of unsupported type: " << bp.type << std::endl;
    }
}

// Fill function for histograms with dynamic suffixes (isotope, detector, etc.)
void HistManager::Fill1D(const std::string& chain, const std::string& base_name, const std::string& suffix, double value, double weight) {
    const std::string full_name = chain + "_" + base_name + suffix;
    auto it = m_hists.find(full_name);
    if (it != m_hists.end()) {
        it->second->Fill(value, weight);
    }
}

void HistManager::Fill2D(const std::string& chain, const std::string& base_name, const std::string& suffix, double x_val, double y_val, double weight) {
    const std::string full_name = chain + "_" + base_name + suffix;
    auto it = m_hists.find(full_name);
    if (it != m_hists.end()) {
        auto* hist2d = dynamic_cast<TH2*>(it->second.get());
        if (hist2d) {
            hist2d->Fill(x_val, y_val, weight);
        }
    }
}

// Simpler versions for histograms without a suffix.
void HistManager::Fill1D(const std::string& chain, const std::string& base_name, double value, double weight) {
    Fill1D(chain, base_name, "", value, weight);
}

void HistManager::Fill2D(const std::string& chain, const std::string& base_name, double x_val, double y_val, double weight) {
    Fill2D(chain, base_name, "", x_val, y_val, weight);
}

void HistManager::SaveAll(const std::string& output_path) {
    TFile* file = TFile::Open(output_path.c_str(), "RECREATE");
    if (!file || file->IsZombie()) {
        throw std::runtime_error("Failed to open output file: " + output_path);
    }
    for (const auto& pair : m_hists) {
        pair.second->Write();
    }
    file->Close();
    delete file;
    std::cout << "All histograms saved to " << output_path << std::endl;
}