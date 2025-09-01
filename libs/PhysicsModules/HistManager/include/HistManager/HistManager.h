#ifndef HIST_MANAGER_H
#define HIST_MANAGER_H

#include <map>
#include <string>
#include <memory>
#include "TH1.h"

// Forward declarations
class ConfigManager;
class AnalysisContext;
class BinningManager;

class HistManager {
public:
    HistManager() = default;
    ~HistManager() = default;

    // Disable copy and move semantics
    HistManager(const HistManager&) = delete;
    HistManager& operator=(const HistManager&) = delete;
    HistManager(HistManager&&) = delete;
    HistManager& operator=(HistManager&&) = delete;

    /**
     * @brief Books all histograms defined in the configuration file.
     * It intelligently creates per-isotope histograms based on their category.
     */
    void BookHistograms(const ConfigManager& config, const AnalysisContext* context, const BinningManager* binningManager);

    void Fill1D(const std::string& name, double value, double weight = 1.0);
    void Fill2D(const std::string& name, double x, double y, double weight = 1.0);

    void SaveHistograms(const std::string& outputFileName);

private:
    template<typename T>
    std::shared_ptr<T> GetHist(const std::string& name);

    std::map<std::string, std::shared_ptr<TH1>> m_histograms;
};

#endif // HIST_MANAGER_H