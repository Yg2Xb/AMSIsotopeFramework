/***********************************************************
 *  File: basic_var.cpp
 *
 *  Modern C++ implementation file for basic variables in AMS Isotopes analysis.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/
#include "basic_var.h"

namespace AMS_Iso {

// IsotopeVar implementation
IsotopeVar::IsotopeVar(int charge, int num, std::string name,
                       std::array<int, Constants::MAX_ISOTOPES> mass,
                       std::array<int, Constants::MAX_ISOTOPES> particle)
    : charge_(charge)
    , isotope_count_(num)
    , name_(std::move(name))
    , mass_(mass)
    , particle_(particle)
{
    if (charge < 1 || charge > Constants::ELEMENT_COUNT) {
        throw IsotopeError(
            "Invalid charge: " + std::to_string(charge) + 
            ". Must be between 1 and " + std::to_string(Constants::ELEMENT_COUNT));
    }
}

// BetaExpoT implementation
BetaExpoT::BetaExpoT(std::string name, double safety_factor,
                     double beta_limit, std::array<int, 2> bin_range)
    : name_(std::move(name))
    , safety_factor_(safety_factor)
    , beta_limit_(beta_limit)
    , bin_range_(bin_range)
{}

namespace {
// Define isotope data
const std::array<IsotopeVar, Constants::ELEMENT_COUNT> IsotopeData {{
    IsotopeVar(1, 2, "Proton",  std::array<int,3>{{1, 2, 0}}, std::array<int,3>{{0, 0, 0}}),
    IsotopeVar(2, 2, "Helium",  std::array<int,3>{{3, 4, 0}}, std::array<int,3>{{46, 47, 0}}),
    IsotopeVar(3, 2, "Lithium", std::array<int,3>{{6, 7, 0}}, std::array<int,3>{{61, 62, 0}}),
    IsotopeVar(4, 3, "Beryllium", std::array<int,3>{{7, 9, 10}}, std::array<int,3>{{63, 64, 114}}),
    IsotopeVar(5, 2, "Boron",   std::array<int,3>{{10, 11, 0}}, std::array<int,3>{{65, 66, 0}}),
    IsotopeVar(6, 2, "Carbon",  std::array<int,3>{{12, 13, 0}}, std::array<int,3>{{67, 117, 0}}),
    IsotopeVar(7, 2, "Nitrogen",std::array<int,3>{{14, 15, 0}}, std::array<int,3>{{68, 118, 0}}),
    IsotopeVar(8, 3, "Oxygen",  std::array<int,3>{{16, 17, 18}}, std::array<int,3>{{69, 0, 0}})
}};
} // anonymous namespace

// Global function implementation
const IsotopeVar& getIsotopeVar(int charge) {
    if (charge < 1 || charge > Constants::ELEMENT_COUNT) {
        throw IsotopeError(
            "Charge must be between 1 and " + std::to_string(Constants::ELEMENT_COUNT));
    }
    return IsotopeData[charge - 1];
}

int findIsotopeIndex(int mass, int charge) {
    if (charge < 1 || charge > Constants::ELEMENT_COUNT) {
        return -1;
    }
    
    const IsotopeVar& isotope = IsotopeData[charge - 1];
    const auto& masses = isotope.getMasses();
    
    // 在该元素的质量数组中查找对应的质量
    for (int i = 0; i < isotope.getIsotopeCount(); ++i) {
        if (masses[i] == mass) {
            return i;
        }
    }
    
    return -1;
}

namespace Binning {
    // Rigidity bins
    const std::array<double, Constants::RIGIDITY_BINS + 1> RigidityBins {{
        0.8, 1.00, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.40, 2.67, 2.97, 3.29, 
        3.64, 4.02, 4.43, 4.88, 5.37, 5.90, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 
        11.0, 12.0, 13.0, 14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 
        28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 
        69.7, 74.9, 80.5, 86.5, 93.0, 100., 108., 116., 125., 135., 147., 160., 
        175., 192., 211., 233., 259., 291., 330., 379., 441., 525., 660., 880., 
        1300., 3300.
    }};

    const std::array<double, 23> EkWideBin = {0.25, 0.42, 0.61, 0.86, 1.17, 1.55, 2.01, 2.57, 
                                            3.23, 4.00, 4.91, 5.99, 7.18, 8.60, 10.25, 12.13, 
                                            13.54, 15.85, 18.49, 21.53, 25.01, 29.04, 33.68};
    
    const std::array<double, 15> BkgEkWideBin {
        0.61, 0.86, 1.17, 1.55, 2.01, 2.57, 3.23, 4.00, 4.91, 5.99, 7.18, 8.60, 10.25, 12.13, 16.38
    };
    
} // namespace Binning

namespace Detector {
    // Beta ExpoTime types
    const std::array<BetaExpoT, Constants::BETA_TYPES> BetaTypes {{
        BetaExpoT("TOF",  1.06,   0.4,   {{1, 7}}),
        BetaExpoT("NaF",  1.005,  0.75,  {{8, 20}}),
        BetaExpoT("Aero", 1.0005, 0.953, {{21, 27}})
    }};


    // RICH detector parameters
    const std::array<std::array<int, Constants::ELEMENT_COUNT>, 2> RichBins {{
        {{400, 400, 240, 140, 200, 240, 240, 260}},  // NaF
        {{800, 800, 280, 140, 320, 400, 320, 336}}   // AGL
    }};
    
    const std::array<double, 2> RichAxis {{0.02, 0.012}};
    
    const std::array<std::array<int, Constants::ELEMENT_COUNT>, 2> TOFRichBins {{
        {{400, 400, 320, 320, 320, 400, 320, 400}},  // NaF
        {{400, 400, 320, 320, 320, 400, 320, 400}}   // AGL
    }};
    
    const std::array<double, 2> TOFAxis {{0.2, 0.2}};
    
    const std::array<std::string, 2> RichDetectorNames {{"NaF", "AGL"}};
} // namespace Detector

} // namespace AMS_Iso