// File: IsoProcess/include/BinningManager.h
// Purpose: Corrected version that does not take a Tool pointer.

#ifndef BINNING_MANAGER_H
#define BINNING_MANAGER_H

#include <vector>
#include <string>
#include <map>

class BinningManager {
public:
    // Constructor no longer needs a Tool pointer
    BinningManager();

    void Initialize(int charge, int mass);

    const std::vector<double>& Get(const std::string& name) const;

private:
    std::vector<double> ConvertRigidityToEk(const std::vector<double>& rig_bins);
    std::vector<double> ConvertRigidityToBeta(const std::vector<double>& rig_bins);

    int m_chargeZ = 0;
    int m_massA = 0;

    mutable std::map<std::string, std::vector<double>> m_bin_map;
};

#endif // BINNING_MANAGER_H