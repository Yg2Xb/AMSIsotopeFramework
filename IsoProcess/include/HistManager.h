// File: IsoProcess/include/HistManager.h
// Purpose: The final version of HistManager, supporting multiple analysis chains.

#ifndef HIST_MANAGER_H
#define HIST_MANAGER_H

#include <string>
#include <vector>
#include <map>
#include <memory>

// Forward declarations to keep this header light
class TH1;
class BinningManager;

class HistManager {
public:
    // We now pass the BinningManager during construction
    HistManager(BinningManager* binner);
    ~HistManager();

    // The main function to create all histograms for all specified chains
    void BookHistograms(bool isMC,
                        const std::vector<std::string>& active_chains,
                        const std::vector<std::string>& active_groups = {});

    // The Fill function is now simpler. It takes the chain name as the first argument.
    void Fill1D(const std::string& chain, const std::string& base_name, const std::string& suffix, double value, double weight = 1.0);
    void Fill2D(const std::string& chain, const std::string& base_name, const std::string& suffix, double x_val, double y_val, double weight = 1.0);
    
    // Simpler versions for histograms without a suffix (like isotope, detector, etc.)
    void Fill1D(const std::string& chain, const std::string& base_name, double value, double weight = 1.0);
    void Fill2D(const std::string& chain, const std::string& base_name, double x_val, double y_val, double weight = 1.0);

    // Saves all created histograms to a file
    void SaveAll(const std::string& output_path);

private:
    // Helper function to create a single histogram instance
    void CreateHistogram(const struct HistogramBlueprint& blueprint);

    BinningManager* m_binner = nullptr; // Pointer to our binning provider
    
    // This map stores all our histogram objects, accessed by their full unique ID
    // (e.g., "L1Inner_ISS.ID.H1_TOF_Be7")
    std::map<std::string, std::unique_ptr<TH1>> m_hists;
};

#endif // HIST_MANAGER_H