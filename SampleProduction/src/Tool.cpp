/***********************************************************
 * File: Tool.cpp
 *
 * Modern C++ implementation file for AMS analysis tools.
 *
 * History:
 * 20241029 - created by ZX.Yan
 ***********************************************************/
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <bitset>
#include <optional>
#include <array>
#include <limits>
#include <TRandom3.h>
#include <iostream>
#include "TFile.h"
#include "Tool.h"
#include <TFile.h>

namespace AMS_Iso {
namespace Tools {

const double geneRig_low = 0.9;
const double geneRig_up  = 3300.0;

TF1 f_MC("f_MC", "1/x", geneRig_low, geneRig_up);
TF1 f_Reweight("f_Reweight", "pow(x, -2.7)", geneRig_low, geneRig_up);
const double MC_norm = f_MC.Integral(geneRig_low, geneRig_up);
const double Reweight_norm = f_Reweight.Integral(geneRig_low, geneRig_up);

// 使用智能指针管理 TF1
std::map<std::string, std::shared_ptr<TF1>> fluxMap;
std::map<std::string, double> fluxNorm;

// 线程安全标志
std::once_flag fluxInitFlag;

std::map<std::string, std::shared_ptr<TF1>>& getFluxMap() { return fluxMap; }
std::map<std::string, double>& getFluxNorm() { return fluxNorm; }

void initFluxFunctions(const std::string& filename) {
    std::call_once(fluxInitFlag, [&](){
        TFile fin(filename.c_str(), "READ");
        if (!fin.IsOpen()) {
            std::cerr << "[ERROR] cannot open flux file: " << filename << std::endl;
            return;
        }

        // 先处理 Be 和 B 的总通量
        TF1* f_Be = dynamic_cast<TF1*>(fin.Get("Be"));
        TF1* f_B = dynamic_cast<TF1*>(fin.Get("B"));

        double integral_Be = 0;
        if (f_Be) {
            integral_Be = f_Be->Integral(geneRig_low, geneRig_up);
            std::cout<<integral_Be<<std::endl;
            // 保存总通量函数
            fluxMap["Be"] = std::shared_ptr<TF1>(static_cast<TF1*>(f_Be->Clone("Be_clone")));
            fluxNorm["Be"] = integral_Be;
        } else {
            std::cerr << "[WARN] flux TF1 Be not found in file!" << std::endl;
        }

        double integral_B = 0;
        if (f_B) {
            integral_B = f_B->Integral(geneRig_low, geneRig_up);
            // 保存总通量函数
            fluxMap["B"] = std::shared_ptr<TF1>(static_cast<TF1*>(f_B->Clone("B_clone")));
            fluxNorm["B"] = integral_B;
        } else {
            std::cerr << "[WARN] flux TF1 B not found in file!" << std::endl;
        }

        // 再处理所有同位素
        std::vector<std::string> names = {"Be7","Be9","Be10", "B10","B11", "C","N","O"};
        for (auto& n : names) {
            TF1* f = dynamic_cast<TF1*>(fin.Get(n.c_str()));
            if (!f) {
                std::cerr << "[WARN] flux TF1 " << n << " not found in file!" << std::endl;
                continue;
            }
            auto f_clone = std::shared_ptr<TF1>(static_cast<TF1*>(f->Clone((n+"_clone").c_str())));
            fluxMap[n] = f_clone;

            if (n == "Be7" || n == "Be9" || n == "Be10") {
                // 按比例计算 Be 同位素的 norm
                if (n == "Be7") fluxNorm[n] = integral_Be * 0.7;
                if (n == "Be9") fluxNorm[n] = integral_Be * 0.2;
                if (n == "Be10") fluxNorm[n] = integral_Be * 0.1;
            } else if (n == "B10" || n == "B11") {
                // 按比例计算 B 同位素的 norm
                if (n == "B10") fluxNorm[n] = integral_B * 0.3;
                if (n == "B11") fluxNorm[n] = integral_B * 0.7;
            } else {
                // 其他元素按原始方式计算 norm
                fluxNorm[n] = f_clone->Integral(geneRig_low, geneRig_up);
            }
            std::cout << "[INIT] Loaded flux TF1: " << n 
                      << "  norm=" << fluxNorm[n] << std::endl;
        }

        fin.Close();
    });
}

void cleanupFluxFunctions() {
    fluxMap.clear();
    fluxNorm.clear();
    // 注意：std::once_flag 无法重置，如果需要重新 init，需要换方案
}

std::string selectFluxName(int charge, int mass) {
    if (charge==4) {
        if (mass==7) return "Be7";
        if (mass==9) return "Be9";
        if (mass==10) return "Be10";
        return "Be";
    }
    if (charge==5) {
        if (mass==10) return "B10";
        if (mass==11) return "B11";
        return "B";
    }
    if (charge==6) return "C";
    if (charge==7) return "N";
    if (charge==8) return "O";
    return "";
}

double calculateWeight(double mmom, int charge, int mass, bool isISS) {
    if (isISS) return 1.0;
    if (charge == 0) {
        std::cerr << "[ERROR] Charge is zero, set weight=0." << std::endl;
        return 0.0;
    }

    double geneRig = mmom / charge;
    if (geneRig < geneRig_low || geneRig > geneRig_up) {
        std::cerr << "[WARN] rigidity out of range: " << geneRig << std::endl;
    }

    std::string name = selectFluxName(charge, mass);
    if (name.empty() || fluxMap.find(name) == fluxMap.end()) {
        std::cerr << "[ERROR] flux TF1 not found for (Z=" << charge << ", A=" << mass << ")" << std::endl;
        return 0.0;
    }

    auto f_flux = fluxMap[name];
    double flux_norm = fluxNorm[name];

    double mc_val   = f_MC.Eval(geneRig) / MC_norm;
    double flux_val = f_flux->Eval(geneRig) / flux_norm;

    return flux_val / mc_val;
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
    double ek_term = ek_per_nucleon / MASS_UNIT + 1;
    return factor * std::sqrt(ek_term * ek_term - 1);
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
        
    for (size_t i = 0; i < Rbins_beta.size(); ++i) {
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
    double beta,       
    bool naf_rad  
    ){

    if(naf_rad) return beta;

    double H=47;
    double n=naf_rad?1.332:1.05;
    double dbeta=naf_rad?0:1.61485e-05;
    double dH=0.239239;
    double dn=0;
    double dh=0;
    
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
    double beta,       
    bool naf_rad  
    ){

    double H=47;
    double n=naf_rad?1.332:1.05;
    double dbeta=naf_rad?1.12123e-04:1.61485e-05;
    double dH=0.239239;
    double dn=naf_rad?2.56350e-03:0;
    double dh=0;

    double nb=n*beta;
    double b=beta;

    double n2=n*n;
    double b2=b*b;
    double nb2=nb*nb;

    double estimate= (-n2*(b2-1)*dn/n
    -((nb2-1)/b2 *(1-(n2-1)*b2)-(n2-1)*(2-n2)) *dH/H
    -((nb2-1)/b2 *sqrt(1-(n2-1)*b2)-(n2-1)*sqrt(2-n2))*dh/H)*b-dbeta;

    return beta+estimate;
}

std::optional<Point2D> calculateXYAtZ(const Float_t positions[9][3],
                                      const Float_t directions[9][3],
                                      double zpl,
                                      int trackIndex) {
    if (trackIndex < 0) {
        return std::nullopt;
    }

    int lmin = 4, lmax = 4;
    double ww = 1;

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

    if (!(ww >= 0 && ww <= 1)) {
        std::cerr << "Error_Interpolate=" << ww << std::endl;
        return std::nullopt;
    }

    Point2D result;
    double point[2];
    
    for (int ixy = 0; ixy < 2; ixy++) {
        double xymin = positions[lmin][ixy] + 
                       (zpl - positions[lmin][2]) * directions[lmin][ixy] / directions[lmin][2];
        
        double xymax = positions[lmax][ixy] + 
                       (zpl - positions[lmax][2]) * directions[lmax][ixy] / directions[lmax][2];
        
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

double GetRichWidth(int iz, bool isNaF) {
    const int nch = 7;
    int zch[nch] = {2, 3, 4, 5, 6, 7, 8};
    
    const int nrad = 2;
    double datasig[nrad][nch]={{0.0021171,  0.00168615, 0.00129161, 0.00108304, 0.00100818,0.000897614,0.000859017},
                                {0.000679357,0.000509869,0.000421296,0.000379986,0.00034206,0.000328994,0.000317768}};

    double datamcdiff[nrad][nch]={{1.12,1.14,1.14,1.14,1.21,1.21,1.21},
                                   {1.08,1.13,1.13,1.13,1.13,1.13,1.13}};
    
    int irad = (isNaF) ? 0 : 1;

    int ich = 0;
    if (iz < zch[0]) {
        ich = 0;
    } else if (iz > zch[nch - 1]) {
        ich = nch - 1;
    } else {
        for (int i = 0; i < nch; i++) {
            if (zch[i] == iz) {
                ich = i;
                break;
            }
        }
    }
    
    double dsig = datasig[irad][ich];
    double msig = datasig[irad][ich] / datamcdiff[irad][ich];

    double datasigBer[nrad]={0.001425, 0.000455};
    double datamcdiffBer[nrad][3]={{1.1709120789,1.1632653061, 1.1623164763},
                                    {1.1375,      1.1346633416, 1.131840796}};

    double wch = 0;
    if (dsig > msig) {
        wch = std::sqrt(std::fabs(std::pow(dsig, 2) - std::pow(msig, 2)));
    }
    
    return wch;
}

double GetSmearRichBeta(int iz, double beta, bool isNaF) {
    if (beta <= 0) return beta;
    
    beta = 1 / beta;
    
    double wid = GetRichWidth(iz, isNaF);
    
    std::random_device rd;
    TRandom3 rand(rd());
    
    double smear = rand.Gaus();
    double newbeta = beta + wid * smear;
    newbeta = 1 / newbeta;
    
    return newbeta;
}


} // namespace Tools
} // namespace AMS_Iso
