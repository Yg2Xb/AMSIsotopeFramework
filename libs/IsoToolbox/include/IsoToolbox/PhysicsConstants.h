#pragma once
#include <string>
#include <vector>
#include <array>
#include <TMath.h> 

// Forward declaration
namespace TMath {
    Double_t AMU();
}

namespace IsoToolbox {

namespace Constants {
    constexpr double AMU_GEV = 0.93149410242;
    constexpr int MAX_ISOTOPES = 3;
    constexpr int ELEMENT_COUNT = 8;
    constexpr int RIGIDITY_BINS = 80;
}

class Isotope {
public:
    Isotope(int charge, int n_isotopes, std::string name,
            const std::array<int, Constants::MAX_ISOTOPES>& masses,
            const std::array<int, Constants::MAX_ISOTOPES>& particle_ids);

    int GetCharge() const { return charge_; }
    int GetIsotopeCount() const { return n_isotopes_; }
    const std::string& GetName() const { return name_; }
    const std::array<int, Constants::MAX_ISOTOPES>& GetMasses() const { return masses_; }

private:
    int charge_;
    int n_isotopes_;
    std::string name_;
    std::array<int, Constants::MAX_ISOTOPES> masses_;
    std::array<int, Constants::MAX_ISOTOPES> particle_ids_;
};

class PhysicsConstants {
public:
    static const Isotope& GetIsotope(int charge);
    static const Isotope& GetIsotope(const std::string& name);
    static const std::vector<double>& GetRigidityBins();
};

class BinningUtils {
public:
    static std::vector<double> ConvertRigidityToEk(const std::vector<double>& rigidity_bins, int mass, int charge);
};

} // namespace IsoToolbox