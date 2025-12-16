/***********************************************************
 * File: Tool.h
 *
 * Modern C++ header file for AMS analysis tools.
 *
 * History:
 * 20241029 - created by ZX.Yan
 ***********************************************************/

#pragma once

#include <array>
#include <algorithm>  // for std::all_of
#include <numeric>  // 为 std::accumulate
#include <optional>
#include <bitset>
#include <random>
#include <iostream>  // 为了使用 std::cerr
#include <TF1.h>  // 为了使用 TF1
#include <map>

namespace AMS_Iso {

template<size_t N>
struct CutResult {
    bool total;
    std::array<bool, N> details;

    CutResult() : total(false), details{} {}
    
    explicit CutResult(const std::array<bool, N>& cuts, bool calculateTotal = true) 
        : details(cuts) {
        total = calculateTotal ? 
                std::all_of(details.begin(), details.end(), [](bool b){ return b; }) 
                : cuts[0];
    }
};

namespace Tools {

// 物理常量
inline constexpr double MASS_UNIT = 0.9315;

// for reweight
extern const double geneRig_low;
extern const double geneRig_up;
extern TF1 f_MC;
extern TF1 f_Reweight;
extern const double MC_norm;
extern const double Reweight_norm;

// 初始化 AMS Flux TF1（线程安全，只执行一次）
void initFluxFunctions(const std::string& filename = "/eos/ams/group/ihep/zixuan/ForSampleProduction/FluxSmooth.root");

// 可选清理函数
void cleanupFluxFunctions();

// 获取 flux TF1 映射
std::map<std::string, std::shared_ptr<TF1>>& getFluxMap();
std::map<std::string, double>& getFluxNorm();

std::string selectFluxName(int charge, double mass);

// 计算事件权重
double calculateWeight(double mmom, int charge, double mass, bool isISS);



// 坐标计算结果类型
struct Point2D {
    double x;
    double y;
};

// 数组平均值计算
double calculateAverage(const double* values, 
                       int count,
                       double ignoreValue = -1000000.001);

// 位置修正
void modifyPositionByZ(double targetZ,
                       std::array<double, 3>& position,
                       double theta,
                       double phi);

// 动能和Beta转换
double betaToKineticEnergy(double beta);
double kineticEnergyToBeta(double kineticEnergy);
// 刚度和Beta转换
double rigidityToBeta(double rigidity, int charge, double mass, bool isElectron = false);
double betaToRigidity(double beta, int charge, double mass, bool isElectron = false);
//
double rigidityToKineticEnergy(double rig_gv, int z, double a);
double kineticEnergyToRigidity(double ek_per_nucleon, int z, double a);
double dR_dEk(double ek_per_nucleon, int z, double a);
void setCutStatus(std::bitset<32>& cutStatus, bool expectedValue, int bitPosition);
//
bool isValidBeta(double beta);
int findBin(std::vector<double> Rbins_beta, double beta);
bool isBeyondCutoff(double beta_low, double cutoffRig, double safetyFactor, 
                    int charge, double UseMass, bool isMC);

//carlos corr
double CorrectCalibrationBias(double beta, bool naf_rad);
//carlos corr 2025 for iss
double CorrectCalibrationBiasInData(double beta, bool naf_rad);

// 位置计算函数
std::optional<Point2D> calculateXYAtZ(const Float_t positions[9][3],
                                      const Float_t directions[9][3],
                                      double zpl,
                                      int trackIndex);

std::optional<Point2D> calculateXYAtNaFZ(const std::array<double, 3>& pos1,
                                         const std::array<double, 3>& pos2,
                                         double nafZ);

struct MassResult {
    double beta;
    double gamma;
    double invMass;
    double ek;
    
    bool isValid() const {
        return isValidBeta(beta) && gamma > 1.0 && 
               invMass > 0.0 && ek >= 0.0;
    }
};
MassResult calculateMass(double beta, double alpha, double innerRig, int charge);

// Function to get the RICH width based on particle charge and radiator type.
// iz: particle charge Z
// isNaF: true for NaF, false for AGL
double GetRichWidth(int iz, bool isNaF);
// Function to get the smeared RICH beta value.
// iz: particle charge Z
// beta: CIEMAT beta after applying corrections
// seed: a random number for each event (e.g., Run + Event number)
// isNaF: true for NaF, false for AGL
double GetSmearRichBeta(int iz, double beta, bool isNaF);
double GetSmearRigidity(double Rigidity, bool isISS, int idet);

// Initializes the charge tuning lookup tables. Thread-safe.
void initChargeTuning(const std::string& filename = "/eos/ams/group/ihep/zixuan/ForSampleProduction/CDFLookupTable_fromSpline.root");

// Performs L2 charge tuning using pre-loaded tables.
double tuneL2Charge(
    const std::string& chain, 
    const std::string& nucleusName, 
    const std::string& detectorName, 
    int ekBin, 
    double q_l2
);

} // namespace Tools
} // namespace AMS_Iso
