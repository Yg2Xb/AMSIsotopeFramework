#include "HistManager/HistManager.h"
#include "IsoToolbox/ConfigManager.h"
#include "IsoToolbox/AnalysisContext.h"
#include "IsoToolbox/BinningManager.h"
#include "IsoToolbox/Logger.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <stdexcept>

void HistManager::BookHistograms(const ConfigManager& config, const AnalysisContext* context, const BinningManager* binningManager) {
    LOG_INFO("Booking histograms...");

    if (!context || !binningManager) {
        throw std::runtime_error("HistManager requires valid AnalysisContext and BinningManager.");
    }

    auto histConfigs = config.GetNode("histograms");
    if (!histConfigs) {
        LOG_WARNING("No 'histograms' section found in config file. No histograms will be booked.");
        return;
    }

    for (const auto& histNode : histConfigs) {
        std::string name = histNode["name"].as<std::string>();
        std::string type = histNode["type"].as<std::string>();
        std::string title = histNode["title"].as<std::string>();
        std::string category = histNode["category"] ? histNode["category"].as<std::string>() : "default";

        if (category == "per_isotope_ek_bin") {
            const auto& isotopes = context->GetIsotopes();
            for (const auto& isotope : isotopes) {
                std::string hist_name = name + "_" + isotope.name;
                std::string hist_title = title + " (" + isotope.name + ")";

                const auto& ek_bins = binningManager->GetEkPerNucleonBins(isotope.name);
                if (ek_bins.size() < 2) {
                    LOG_ERROR("Skipping histogram %s for isotope %s: requires at least 2 bin edges, but got %zu.", name.c_str(), isotope.name.c_str(), ek_bins.size());
                    continue;
                }

                if (type == "TH2D") {
                    int nbinsx = histNode["nbinsx"].as<int>();
                    double xlow = histNode["xlow"].as<double>();
                    double xup = histNode["xup"].as<double>();
                    m_histograms[hist_name] = std::make_shared<TH2D>(hist_name.c_str(), hist_title.c_str(), nbinsx, xlow, xup, ek_bins.size() - 1, ek_bins.data());
                    LOG_DEBUG("Booked per-isotope TH2D: %s", hist_name.c_str());
                } else if (type == "TH1D") {
                     m_histograms[hist_name] = std::make_shared<TH1D>(hist_name.c_str(), hist_title.c_str(), ek_bins.size() - 1, ek_bins.data());
                     LOG_DEBUG("Booked per-isotope TH1D: %s", hist_name.c_str());
                }
            }
        } else { // Default category for non-isotope-specific histograms
            if (type == "TH1D") {
                int nbins = histNode["nbins"].as<int>();
                double low = histNode["low"].as<double>();
                double up = histNode["up"].as<double>();
                m_histograms[name] = std::make_shared<TH1D>(name.c_str(), title.c_str(), nbins, low, up);
                LOG_DEBUG("Booked default TH1D: %s", name.c_str());
            }
            // Add other types like TH2D if needed for the 'default' category
        }
    }
    LOG_INFO("All histograms booked.");
}

template<typename T>
std::shared_ptr<T> HistManager::GetHist(const std::string& name) {
    auto it = m_histograms.find(name);
    if (it == m_histograms.end()) {
        LOG_WARNING("Histogram '%s' not found.", name.c_str());
        return nullptr;
    }
    auto hist = std::dynamic_pointer_cast<T>(it->second);
    if (!hist) {
        LOG_WARNING("Histogram '%s' found, but has incorrect type.", name.c_str());
        return nullptr;
    }
    return hist;
}

void HistManager::Fill1D(const std::string& name, double value, double weight) {
    if (auto hist = GetHist<TH1D>(name)) {
        hist->Fill(value, weight);
    }
}

void HistManager::Fill2D(const std::string& name, double x, double y, double weight) {
    if (auto hist = GetHist<TH2D>(name)) {
        hist->Fill(x, y, weight);
    }
}

void HistManager::SaveHistograms(const std::string& outputFileName) {
    auto outFile = std::make_unique<TFile>(outputFileName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        LOG_ERROR("Failed to create output file: %s", outputFileName.c_str());
        return;
    }
    outFile->cd();
    for (const auto& pair : m_histograms) {
        pair.second->Write();
    }
    outFile->Close();
    LOG_INFO("Histograms successfully saved to %s", outputFileName.c_str());
}