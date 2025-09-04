/***********************************************************
 *  File: Tool.cpp
 *
 *  Modern C++ implementation file for AMS analysis tools.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <bitset>
#include <optional>
#include <array>
#include <limits>
#include "Tool.h"

namespace AMS_Iso {
namespace Tools {

// 定义全局变量
const double geneRig_low = 1.0;     
const double geneRig_up = 2000.0;  
TF1 f_MC("f_MC", "1/x", geneRig_low, geneRig_up);
TF1 f_Reweight("f_Reweight", "pow(x, -2.7)", geneRig_low, geneRig_up);
const double MC_norm = f_MC.Integral(geneRig_low, geneRig_up);
const double Reweight_norm = f_Reweight.Integral(geneRig_low, geneRig_up);

// 函数实现...
double calculateWeight(double mmom, int charge, bool isISS) {
    if (isISS) return 1.0;
    if (charge == 0) {
        std::cerr << "Charge is zero, setting weight to 0." << std::endl;
        return 0.0;
    }
    double geneRig = mmom / charge;
    return (f_Reweight.Eval(geneRig) / Reweight_norm) / (f_MC.Eval(geneRig) / MC_norm);
}

double calculateAverage(const double* values, int count, double ignoreValue) {
    if (!values || count <= 0) return -1000000.0;

    double sum = 0.0;
    int validCount = 0;

    for (int i = 0; i < count; ++i) {
        if (values[i] > 0 && std::abs(values[i] - ignoreValue) > std::numeric_limits<double>::epsilon()) {
            sum += values[i];
            ++validCount;
        }
    }

    return validCount > 0 ? sum / validCount : -1000000.0;
}

void modifyPositionByZ(double targetZ, std::array<double, 3>& position,
                      double theta, double phi) {
    if (std::abs(theta) > M_PI || std::abs(phi) > 2 * M_PI) {
        throw std::invalid_argument("Invalid angle input (theta or phi)");
    }

    double deltaZ = targetZ - position[2];
    position[0] += deltaZ * std::tan(theta) * std::cos(phi);
    position[1] += deltaZ * std::tan(theta) * std::sin(phi);
    position[2] = targetZ;
}

double betaToKineticEnergy(double beta) {
    if (beta <= 0.0 || beta >= 1.0) return -9.0;
    
    double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    return (gamma - 1.0) * MASS_UNIT;
}

double kineticEnergyToBeta(double kineticEnergy) {
    if (kineticEnergy < 0.0) return -9.0;
    
    double gamma = kineticEnergy / MASS_UNIT + 1.0;
    return std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

double rigidityToBeta(double rigidity, int charge, int mass, bool isElectron) {
    if (!isElectron && mass < charge) {
        throw std::invalid_argument("Invalid charge/mass combination");
    }

    if (isElectron) {
        constexpr double ELECTRON_MASS = 0.000511;
        double beta = rigidity * std::sqrt(1.0 / (ELECTRON_MASS * ELECTRON_MASS + rigidity * rigidity));
        return (charge == -1) ? beta : -beta;
    } else {
        double particleMass = mass * MASS_UNIT;
        return rigidity * charge * std::sqrt(1.0 / (particleMass * particleMass + rigidity * rigidity * charge * charge));
    }
}

double betaToRigidity(double beta, int charge, int mass, bool isElectron) {
    if (beta <= 0.0 || beta >= 1.0) return -100000.0;

    if (isElectron) {
        constexpr double ELECTRON_MASS = 0.000511;
        double rigidity = beta * ELECTRON_MASS / std::sqrt(1.0 - beta * beta);
        return (charge == -1) ? rigidity : -rigidity;
    } else {
        double particleMass = mass * MASS_UNIT;
        return beta * particleMass / (std::sqrt(1.0 - beta * beta) * charge);
    }
}

double rigidityToKineticEnergy(double rig_gv, int z, int a) {
    if (rig_gv <= 0.0 || z == 0) return -9.0;
    
    double factor = (a * MASS_UNIT) / z;
    double term = std::pow(rig_gv / factor, 2);
    return MASS_UNIT * (std::sqrt(1 + term) - 1);
}

double kineticEnergyToRigidity(double ek_per_nucleon, int z, int a) {
    if (ek_per_nucleon < 0.0 || z == 0) return -100000.0;
    
    double factor = (a * MASS_UNIT) / z;
    double term = std::pow(ek_per_nucleon / MASS_UNIT + 1, 2) - 1;
    return factor * std::sqrt(term);
}

double dR_dEk(double ek_per_nucleon, int z, int a) {
    if (ek_per_nucleon < 0.0 || z == 0) return -100000.0;
    
    double factor = (a * MASS_UNIT) / z;
    double ek_term = ek_per_nucleon / MASS_UNIT + 1;
    return factor * ek_term / std::sqrt(ek_term * ek_term - 1);
}

void setCutStatus(std::bitset<32>& cutStatus, bool expectedValue, int bitPosition) {
    if (bitPosition < 0 || bitPosition > 31) {
        throw std::invalid_argument("Invalid bit position");
    }
    cutStatus[bitPosition] = expectedValue;
}

bool isValidBeta(double beta) {
    return beta > 0 && beta < 1;
}

int findBin(std::vector<double> Rbins_beta, double beta) {
    if (!isValidBeta(beta)) return -1;
    if (beta < Rbins_beta.front() || beta >= Rbins_beta.back()) {
        return -1;
    }
        
    for (int i = 0; i < Rbins_beta.size(); ++i) {
        if (beta >= Rbins_beta[i] && beta < Rbins_beta[i + 1]) {
            return i;
        }
    }
    return -1;
}
bool isBeyondCutoff(double beta_low, double cutoffRig, double safetyFactor, int charge, int UseMass,  bool isMC) {
    if (isMC) return true;
    if (!isValidBeta(beta_low)) return false;
        
    double cutoffBeta = rigidityToBeta(cutoffRig, charge, UseMass, false);
    if (!isValidBeta(cutoffBeta)) return false;
        
    return beta_low > safetyFactor * cutoffBeta;
}

double CorrectCalibrationBias(
                                         double beta,        // Reconstructed calibrated beta                                                                
                                         bool naf_rad   // True if the radiator is naf
                                         ){

    if(naf_rad) return beta;             // So far NaF is not corrected

    double H=47;                          // Reference expansion height (cm)                                                         
    double n=naf_rad?1.332:1.05;          // Reference refractive index                                                              
    double dbeta=naf_rad?0:1.61485e-05;   // Bias in beta due to the finitie rigidity used during calibration
    double dH=0.239239;                   // Bias in the expansion length in cm
    double dn=0;                          // Bias in the average refractive index: for agl it is fully degenerated
    double dh=0;                          // Bias in the average emission point in cm. For agl it is fullly degenerated
    
    double nb=n*beta;
    double b=beta;

    double n2=n*n;
    double b2=b*b;
    double nb2=nb*nb;

   
    double estimate= (-n2*(b2-1)*dn/n
                    -((nb2-1)/b2 *(1-(n2-1)*b2)-(n2-1)*(2-n2)) *dH/H
                    -((nb2-1)/b2 *sqrt(1-(n2-1)*b2)-(n2-1)*sqrt(2-n2))*dh/H)*b-dbeta;


    return beta-estimate;
}

double CorrectCalibrationBiasInData(
    double beta,        // Reconstructed and corrected beta
    bool naf_rad   // True if the radiator is naf
    ){

    double H=47;                                    // Reference expansion height (cm)
    double n=naf_rad?1.332:1.05;          // Reference refractive index
    double dbeta=naf_rad?1.12123e-04:1.61485e-05;   // Bias in beta due to the finitie rigidity used during calibration
    double dH=0.239239;                       // Bias in the expansion length in cm
    double dn=naf_rad?2.56350e-03:0;  // Bias in the average refractive index: for agl it is fully degenerated
    double dh=0;                                     // Bias in the average emission point in cm. For agl it is fullly degenerated

    double nb=n*beta;
    double b=beta;

    double n2=n*n;
    double b2=b*b;
    double nb2=nb*nb;

    double estimate= (-n2*(b2-1)*dn/n
    -((nb2-1)/b2 *(1-(n2-1)*b2)-(n2-1)*(2-n2)) *dH/H
    -((nb2-1)/b2 *sqrt(1-(n2-1)*b2)-(n2-1)*sqrt(2-n2))*dh/H)*b-dbeta;

    return beta+estimate;  // The plus sign is for correcting ISS data, where the shift in expansion length takes place. To correct MC the sign should be positive
}

std::optional<Point2D> calculateXYAtZ(const Float_t positions[9][3],
                                    const Float_t directions[9][3],
                                    double zpl,
                                    int trackIndex = 0) {
    // 检查轨迹索引是否有效
    if (trackIndex < 0) {
        return std::nullopt;
    }

    // 初始化搜索范围
    int lmin = 4, lmax = 4;
    double ww = 1;

    // 确定插值区间
    if (zpl >= positions[0][2]) {
        lmin = lmax = 0;
        ww = 1;
    }
    else if (zpl <= positions[8][2]) {
        lmin = lmax = 8;
        ww = 1;
    }
    else {
        for (int ilay = 0; ilay < 8; ilay++) {
            if (zpl <= positions[ilay][2] && zpl >= positions[ilay+1][2]) {
                ww = (zpl - positions[ilay+1][2]) / 
                     (positions[ilay][2] - positions[ilay+1][2]);
                lmin = ilay;
                lmax = ilay + 1;
                break;
            }
        }
    }

    // 权重检查
    if (!(ww >= 0 && ww <= 1)) {
        std::cerr << "Error_Interpolate=" << ww << std::endl;
        return std::nullopt;
    }

    Point2D result;
    double point[2];
    
    // 计算 x 和 y 坐标
    for (int ixy = 0; ixy < 2; ixy++) {
        double xymin = positions[lmin][ixy] + 
                      (zpl - positions[lmin][2]) * 
                      directions[lmin][ixy] / directions[lmin][2];
        
        double xymax = positions[lmax][ixy] + 
                      (zpl - positions[lmax][2]) * 
                      directions[lmax][ixy] / directions[lmax][2];
        
        point[ixy] = ww * xymin + (1 - ww) * xymax;
    }

    result.x = point[0];
    result.y = point[1];
    
    return result;
}

std::optional<Point2D> calculateXYAtNaFZ(const std::array<double, 3>& pos1,
                                        const std::array<double, 3>& pos2,
                                        double nafZ) {
    double dx = pos2[0] - pos1[0];
    double dy = pos2[1] - pos1[1];
    double dz = pos2[2] - pos1[2];

    if (std::abs(dz) < std::numeric_limits<double>::epsilon()) return std::nullopt;

    double t = (nafZ - pos1[2]) / dz;

    return Point2D{
        pos1[0] + t * dx,
        pos1[1] + t * dy
    };
}

} // namespace Tools
} // namespace AMS_Iso