#include "IsoToolbox/PhysicsConstants.h"
#include <stdexcept>

namespace IsoToolbox {

namespace {
    const std::array<Isotope, Constants::ELEMENT_COUNT> IsotopeData {{
        Isotope(1, 2, "Proton",     {{1, 2, 0}}, {{0, 0, 0}}),
        Isotope(2, 2, "Helium",     {{3, 4, 0}}, {{46, 47, 0}}),
        Isotope(3, 2, "Lithium",    {{6, 7, 0}}, {{61, 62, 0}}),
        Isotope(4, 3, "Beryllium",  {{7, 9, 10}}, {{63, 64, 114}}),
        Isotope(5, 2, "Boron",      {{10, 11, 0}}, {{65, 66, 0}}),
        Isotope(6, 2, "Carbon",     {{12, 13, 0}}, {{67, 117, 0}}),
        Isotope(7, 2, "Nitrogen",   {{14, 15, 0}}, {{68, 118, 0}}),
        Isotope(8, 3, "Oxygen",     {{16, 17, 18}}, {{69, 0, 0}})
    }};

    const std::vector<double> RigidityBinsVec = {
        0.8, 1.00, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.40, 2.67, 2.97, 3.29, 
        3.64, 4.02, 4.43, 4.88, 5.37, 5.90, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 
        11.0, 12.0, 13.0, 14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 
        28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 
        69.7, 74.9, 80.5, 86.5, 93.0, 100., 108., 116., 125., 135., 147., 160., 
        175., 192., 211., 233., 259., 291., 330., 379., 441., 525., 660., 880., 
        1300., 3300.
    };
}

Isotope::Isotope(int charge, int n_isotopes, std::string name,
                 const std::array<int, Constants::MAX_ISOTOPES>& masses,
                 const std::array<int, Constants::MAX_ISOTOPES>& particle_ids)
    : charge_(charge), n_isotopes_(n_isotopes), name_(std::move(name)), masses_(masses), particle_ids_(particle_ids) {}

const Isotope& PhysicsConstants::GetIsotope(int charge) {
    if (charge < 1 || charge > Constants::ELEMENT_COUNT) {
        throw std::out_of_range("Charge must be between 1 and " + std::to_string(Constants::ELEMENT_COUNT));
    }
    return IsotopeData[charge - 1];
}

const Isotope& PhysicsConstants::GetIsotope(const std::string& name) {
    for (const auto& iso : IsotopeData) {
        if (iso.GetName() == name) {
            return iso;
        }
    }
    throw std::out_of_range("Isotope with name " + name + " not found.");
}

const std::vector<double>& PhysicsConstants::GetRigidityBins() {
    return RigidityBinsVec;
}

std::vector<double> BinningUtils::ConvertRigidityToEk(const std::vector<double>& rigidity_bins, int mass, int charge) {
    std::vector<double> ek_bins;
    ek_bins.reserve(rigidity_bins.size());
    
    for (double R : rigidity_bins) {
        
        double factor = (mass * Constants::AMU_GEV) / charge;
        double term = std::pow(R / factor, 2);
        double E_kinetic_perNuclear = Constants::AMU_GEV * (std::sqrt(1 + term) - 1);
        
        ek_bins.push_back(E_kinetic_perNuclear);
    }
    return ek_bins;
}

} // namespace IsoToolbox