/***********************************************************
 *  File: Tool.h
 *
 *  Modern C++ header file for AMS analysis tools.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/

#pragma once

#include <array>
#include <algorithm>  // for std::all_of
#include <numeric>  // 为 std::accumulate
#include <optional>
#include <bitset>
#include <iostream>  // 为了使用 std::cerr
#include <TF1.h>  // 为了使用 TF1

namespace AMS_Iso {

template<size_t N>
struct CutResult {
    bool total;                  
    std::array<bool, N> details;

    // 默认构造函数
    CutResult() : total(false), details{} {}
    
    // 主构造函数
    // isPrecomputedTotal = true (默认): 计算所有details的AND结果
    // isPrecomputedTotal = false: 使用cuts[0]作为total
    explicit CutResult(const std::array<bool, N>& cuts, 
                      bool calculateTotal = true) 
        : details(cuts) {
        total = calculateTotal ? 
                std::all_of(details.begin(), details.end(), 
                           [](bool b) { return b; }) 
                : cuts[0];
    }
};

namespace Tools {

// 物理常量
inline constexpr double MASS_UNIT = 0.931;  // 核子平均质量
// for reweight
extern const double geneRig_low;     
extern const double geneRig_up;  
extern TF1 f_MC;
extern TF1 f_Reweight;
extern const double MC_norm;
extern const double Reweight_norm;

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
double rigidityToBeta(double rigidity, int charge, int mass, bool isElectron = false);
double betaToRigidity(double beta, int charge, int mass, bool isElectron = false);
//
double rigidityToKineticEnergy(double rig_gv, int z, int a);
double kineticEnergyToRigidity(double ek_per_nucleon, int z, int a);
double dR_dEk(double ek_per_nucleon, int z, int a);
void setCutStatus(std::bitset<32>& cutStatus, bool expectedValue, int bitPosition);
//
bool isValidBeta(double beta);
int findBin(std::vector<double> Rbins_beta, double beta);
bool isBeyondCutoff(double beta_low, double cutoffRig, double safetyFactor, 
                    int charge, int UseMass, bool isMC);

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

double calculateWeight(double mmom, int charge, bool isISS);
} // namespace Tools
} // namespace AMS_Iso