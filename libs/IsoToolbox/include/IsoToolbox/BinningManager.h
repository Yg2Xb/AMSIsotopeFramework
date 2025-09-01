#ifndef BINNINIG_MANAGER_H
#define BINNINIG_MANAGER_H

#include <string>
#include <vector>
#include <map>

// Forward declarations to reduce header dependencies
class AnalysisContext;

/**
 * @class BinningManager
 * @brief Manages the creation of isotope-specific binning schemes.
 *
 * This class is responsible for generating and providing binning arrays
 * (e.g., Ek/n bins) based on a master set of rigidity bins. It fetches
 * information about the target nucleus and its isotopes from the AnalysisContext
 * and performs the necessary physics calculations.
 */
class BinningManager {
public:
    /**
     * @brief Constructor.
     * @param context A non-null pointer to the global AnalysisContext.
     */
    explicit BinningManager(const AnalysisContext* context);
    ~BinningManager() = default;

    // Disable copy and move semantics to prevent slicing and ownership issues
    BinningManager(const BinningManager&) = delete;
    BinningManager& operator=(const BinningManager&) = delete;
    BinningManager(BinningManager&&) = delete;
    BinningManager& operator=(BinningManager&&) = delete;

    /**
     * @brief Initializes the manager by generating all required bin schemes.
     * This must be called after the AnalysisContext is fully initialized.
     */
    void Initialize();

    /**
     * @brief Retrieves the Ek/n binning for a specific isotope.
     * @param isotopeName The name of the isotope (e.g., "Be7").
     * @return A const reference to the vector of bin edges. Returns an empty vector if not found.
     */
    const std::vector<double>& GetEkPerNucleonBins(const std::string& isotopeName) const;

private:
    void GenerateEkPerNucleonBins();

    const AnalysisContext* m_context;
    std::map<std::string, std::vector<double>> m_ekPerNucleonBins;
    std::vector<double> m_emptyBins; // Return this if isotope not found to avoid dangling refs
};

#endif // BINNINIG_MANAGER_H