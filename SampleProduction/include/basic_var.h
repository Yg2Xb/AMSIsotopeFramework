/***********************************************************
 *  File: basic_var.h
 *
 *  Modern C++ header file for basic variables in AMS Isotopes analysis.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/

#pragma once
#include <array>
#include <string>
#include <stdexcept>
#include <vector>
#include "TString.h"
#include "TChain.h"

namespace AMS_Iso {
	
// Version information
namespace Version {
    inline constexpr int MAJOR = 1;
    inline constexpr int MINOR = 0;
    inline constexpr int PATCH = 0;
    inline constexpr const char* DATE = "20241029";
    
    // 可选：提供一个获取版本字符串的函数
    inline std::string getVersionString() {
        return std::to_string(MAJOR) + "." + 
               std::to_string(MINOR) + "." + 
               std::to_string(PATCH);
    }
}
// Forward declarations
class IsotopeVar;
class BetaExpoT;

//!!!!!!!!!!!!!!!!!!!!
inline constexpr bool yanzx_dst = false; 



namespace Constants {
    inline constexpr int MAX_ISOTOPES = 3;
    inline constexpr int ELEMENT_COUNT = 8;
    inline constexpr int RIGIDITY_BINS = 73;
    inline constexpr int EFF_RIGIDITY_BINS = 40;
    inline constexpr int BETA_TYPES = 3;
    inline constexpr int TRACKER_CUTS = 5;
    inline constexpr double SAFE_FACTOR_RIG = 1.2;

    inline constexpr  int N_nuc = 13;
    inline constexpr  int nuclei_Z[N_nuc] = {2,2,3,3,4,4,4,5,5,6,7,7,8};
    inline constexpr  int nuclei_A[N_nuc] = {3,4,6,7,7,9,10,10,11,12,14,15,16};
}

	const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
	const std::vector<std::string> cut_groups = {
		"BasicAndFiducial", "Trigger", "InnerTracker", "L1BigZ", "L1PickUp", "L1UpperQ", "UpperTOFQ", "BkgReduction",
		"TOFGeo", "TOFBetaQuality", "NaFGeo", "NaFReconstruction", "AGLGeo", "AGLReconstruction"
	};
	const std::vector<std::string> num_den = {"Num", "Den"};
	const std::vector<std::string> charge_types = {"L1QSignal", "L1QTemplate", "L2QTemplate"};
	const std::vector<std::string> sources = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
	const std::vector<std::string> gene_rec = {"Gene", "Rec"};

const std::map<int, std::vector<std::pair<int, double>>> isotopeWeights = {
		{2, {{3, 0.1}, {4, 0.9}}},  // He
		{3, {{6, 0.5}, {7, 0.5}}},  // Li
		{4, {{7, 0.6}, {9, 0.3}, {10, 0.1}}},  // Be
		{5, {{10, 0.3}, {11, 0.7}}},           // B
		{6, {{12, 1.0}}},                      // C
		{7, {{14, 0.5}, {15, 0.5}}},           // N
		{8, {{16, 1.0}}}                       // O
};

// Exception classes
class IsotopeError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

class IsotopeVar {
public:
    // 构造函数
    IsotopeVar(int charge, int num, std::string name,
               std::array<int, Constants::MAX_ISOTOPES> mass,
               std::array<int, Constants::MAX_ISOTOPES> particle);

    // Getters
    int getCharge() const { return charge_; }
    int getIsotopeCount() const { return isotope_count_; }
    const std::string& getName() const { return name_; }
    const std::array<int, Constants::MAX_ISOTOPES>& getMasses() const { return mass_; }
    const std::array<int, Constants::MAX_ISOTOPES>& getParticles() const { return particle_; }
    int getParticleIDByMass(int mass) const;
    
    int getMass(int index) const {
        if (index < 0 || index >= isotope_count_) {
            throw IsotopeError("Invalid isotope index in getMass");
        }
        return mass_[index];
    }

private:
    int charge_;            // 电荷数
    int isotope_count_;     // 同位素数量
    std::string name_;      // 元素名称
    std::array<int, Constants::MAX_ISOTOPES> mass_;      // 质量数组
    std::array<int, Constants::MAX_ISOTOPES> particle_;  // 粒子ID数组
};

class BetaExpoT {
public:
    BetaExpoT(std::string name, double safety_factor, 
                   double beta_limit, std::array<int, 2> bin_range);

    const std::string& getName() const { return name_; }
    double getSafetyFactor() const { return safety_factor_; }
    double getBetaLimit() const { return beta_limit_; }
    const std::array<int, 2>& getBinRange() const { return bin_range_; }

private:
    std::string name_;
    double safety_factor_;
    double beta_limit_;
    std::array<int, 2> bin_range_;
};

namespace Binning {
    // Rigidity bins and related constants
    extern const std::array<double, Constants::RIGIDITY_BINS + 1> RigidityBins;
    extern const double RigiditySafetyFactor;
    
    // Efficiency estimation bins
    extern const std::array<double, Constants::EFF_RIGIDITY_BINS + 1> EffRigidityBins;

    // wei jh ek bin
    extern const std::array<double, 23> EkWideBin;
    extern const std::array<double, 15> BkgEkWideBin;
    
    extern const std::array<double, 38> common_RigidityBins;
    extern const std::array<double, 38> common_EkBins;
    extern const std::array<double, 74> old_EkBins;

    // Run parameters
    inline constexpr double RunBegin = 1304179200;
    inline constexpr double RunEnd = 1672502400;
    inline constexpr int RunBins = 3532;
}

enum class DetType {
    TOF = 0,
    NaF = 1,
    AGL = 2
};

namespace Detector {
    // Beta measurement types
    extern const std::array<BetaExpoT, Constants::BETA_TYPES> BetaTypes;
    
    // Tracker cuts
    extern const std::array<std::string, Constants::TRACKER_CUTS> TrackerCutNames;
    extern const std::array<std::string, 11> TOFCuts;
    extern const std::array<std::string, 11> RICHCuts;
    
    // RICH detector parameters
    extern const std::array<std::array<int, Constants::ELEMENT_COUNT>, 2> RichBins;
    extern const std::array<double, 2> RichAxis;
    extern const std::array<std::array<int, Constants::ELEMENT_COUNT>, 2> TOFRichBins;
    extern const std::array<double, 2> TOFAxis;
    extern const std::array<std::string, 2> RichDetectorNames;
}

// Global functions
const IsotopeVar& getIsotopeVar(int charge);
int findIsotopeIndex(int mass, int charge = 4);

} // namespace AMS_Iso