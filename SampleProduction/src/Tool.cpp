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
#include "Tool.h"
#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"

namespace AMS_Iso {
namespace Tools {

const double geneRig_low = 1.0;     
const double geneRig_up = 2000.0;  

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

        const std::vector<std::string> tags = {
            "Be7", "Be9", "Be10", "B10", "B11", "C12", "N14", "N15", "O16"
        };

        for (const auto& t : tags) {
            const std::string key = t + "_spline_R";
            TF1* f = dynamic_cast<TF1*>(fin.Get(key.c_str()));
            if (!f) {
                std::cerr << "[WARN] flux TF1 not found in file: " << key << std::endl;
                continue;
            }

            auto f_clone = std::shared_ptr<TF1>(static_cast<TF1*>(f->Clone((t + "_clone").c_str())));
            f_clone->SetNpx(2000);

            fluxMap[t]  = f_clone;
            const double integ = f_clone->Integral(geneRig_low, geneRig_up);
            fluxNorm[t] = integ;

            if (!(integ > 0)) {
                std::cerr << "[WARN] norm<=0 for " << t
                          << " over [" << geneRig_low << "," << geneRig_up << "]" << std::endl;
            }
            std::cout << "[INIT] Loaded TF1 " << key << "  norm=" << integ << std::endl;
        }

        fin.Close();
    });
}

void cleanupFluxFunctions() {
    fluxMap.clear();
    fluxNorm.clear();
    // 注意：std::once_flag 无法重置，如果需要重新 init，需要换方案
}

std::string selectFluxName(int charge, double mass) {
    if (charge == 4) {
        if (mass == 7)  return "Be7";
        if (mass == 9)  return "Be9";
        if (mass == 10) return "Be10";
        return ""; // 不再回退到 "Be"
    }
    if (charge == 5) {
        if (mass == 10) return "B10";
        if (mass == 11) return "B11";
        return ""; // 不再回退到 "B"
    }
    if (charge == 6) return "C12";
    if (charge == 7) {
        if (mass == 14) return "N14";
        if (mass == 15) return "N15";
        return "";
    }
    if (charge == 8) return "O16";
    return "";
}

double calculateWeight(double mmom, int charge, double mass, bool isISS) {
    if (isISS) return 1.0;
    if (charge == 0) return 0.0;

    const double geneRig = std::abs(mmom) / std::abs(charge);
    if (geneRig < geneRig_low || geneRig > geneRig_up) {
        return 0.0;
    }

    const std::string name = selectFluxName(charge, mass);
    if (name.empty()) {
        std::cerr << "[ERROR] flux name not found for (Z=" << charge << ", A=" << mass << ")" << std::endl;
        return 0.0;
    }

    auto itF = fluxMap.find(name);
    auto itN = fluxNorm.find(name);
    if (itF == fluxMap.end() || itN == fluxNorm.end()) {
        std::cerr << "[ERROR] flux TF1 or norm missing for " << name << std::endl;
        return 0.0;
    }

    const double mc_den = MC_norm;
    const double fl_den = itN->second;
    if (!(mc_den > 0) || !(fl_den > 0)) return 0.0;

    const double mc_val   = f_MC.Eval(geneRig) / mc_den;
    const double flux_val = itF->second->Eval(geneRig) / fl_den;

    if (!(mc_val > 0)) return 0.0;
    return flux_val / mc_val;
}

MassResult calculateMass(double beta, double alpha, double innerRig, int charge) {
    MassResult result{0.0, 0.0, 0.0, 0.0};
    if (!isValidBeta(beta) || innerRig < 0.0) return result;
    
    result.beta = alpha * beta / std::sqrt(1 - beta * beta + std::pow(alpha * beta, 2));
    if (!isValidBeta(result.beta)) return {0.0, 0.0, 0.0, 0.0};
    
    result.gamma = 1.0 / std::sqrt(1.0 - result.beta * result.beta);
    result.ek = (result.gamma - 1) * MASS_UNIT;
    result.invMass = (result.beta * result.gamma) / (charge * innerRig);
    
    return result;
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

double rigidityToBeta(double rigidity, int charge, double mass, bool isElectron) {
    if (!isElectron && mass < charge) {
        throw std::invalid_argument("Invalid charge/mass combination");
    }

    if (isElectron) {
        constexpr double ELECTRON_MASS = 0.000511;
        double beta = std::abs(rigidity) * std::sqrt(1.0 / (ELECTRON_MASS * ELECTRON_MASS + rigidity * rigidity));
        return (charge == -1) ? beta : -beta;
    } else {
        double particleMass = mass * MASS_UNIT;
        return std::abs(rigidity) * charge * std::sqrt(1.0 / (particleMass * particleMass + rigidity * rigidity * charge * charge));
    }
}

double betaToRigidity(double beta, int charge, double mass, bool isElectron) {
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

double rigidityToKineticEnergy(double rig_gv, int z, double a) {
    if (rig_gv <= 0.0 || z == 0) return -9.0;
    
    double factor = (a * MASS_UNIT) / z;
    double term = std::pow(rig_gv / factor, 2);
    return MASS_UNIT * (std::sqrt(1 + term) - 1);
}

double kineticEnergyToRigidity(double ek_per_nucleon, int z, double a) {
    if (ek_per_nucleon < 0.0 || z == 0) return -100000.0;
    
    double factor = (a * MASS_UNIT) / z;
    double ek_term = ek_per_nucleon / MASS_UNIT + 1;
    return factor * std::sqrt(ek_term * ek_term - 1);
}

double dR_dEk(double ek_per_nucleon, int z, double a) {
    if (ek_per_nucleon < 0.0 || z == 0) return -100000.0;
    
    double factor = 1.0*a / z;
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

bool isBeyondCutoff(double beta_low, double cutoffRig, double safetyFactor, int charge, double UseMass,  bool isMC) {
    if (isMC) return true;
    if (beta_low > 1) return true;
        
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
    double dH=0.239239*1.5; //2025Nov5,new correction factor 1.4 // 11.12 factor 1.5
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

// Internal helper class to manage a single CDF lookup table.
class LookupTable {
public:
    bool isValid = false;

    void loadFromVectors(TVectorD* q_vec, TVectorD* cdf1_vec, TVectorD* cdf2_vec) {
        isValid = false;
        if (q_vec && cdf1_vec && cdf2_vec && q_vec->GetNrows() > 1) {
            q_values.assign(q_vec->GetMatrixArray(), q_vec->GetMatrixArray() + q_vec->GetNrows());
            cdf_l1.assign(cdf1_vec->GetMatrixArray(), cdf1_vec->GetMatrixArray() + cdf1_vec->GetNrows());
            cdf_l2.assign(cdf2_vec->GetMatrixArray(), cdf2_vec->GetMatrixArray() + cdf2_vec->GetNrows());
            q_min = q_values.front();
            q_max = q_values.back();
            dq = (q_values.size() > 1) ? (q_max - q_min) / (q_values.size() - 1) : 0.0;
            if (std::abs(dq) > 1e-9) isValid = true;
        }
    }
    
    double getCDF_L2(double q) const {
        if (!isValid || q < q_min || q > q_max) return -1.0;
        double fidx = (q - q_min) / dq;
        int idx = static_cast<int>(fidx);
        if (idx >= static_cast<int>(cdf_l2.size()) - 1) return cdf_l2.back();
        if (idx < 0) return cdf_l2.front();
        double frac = fidx - idx;
        return cdf_l2[idx] + frac * (cdf_l2[idx+1] - cdf_l2[idx]);
    }
    
    double getInvCDF_L1(double cdf_target) const {
        if (!isValid || cdf_target < 0.0 || cdf_target > 1.0) return -999.0;
        auto it = std::lower_bound(cdf_l1.begin(), cdf_l1.end(), cdf_target);
        if (it == cdf_l1.end()) return q_max;
        if (it == cdf_l1.begin()) return q_min;
        int idx = std::distance(cdf_l1.begin(), it);
        double den = cdf_l1[idx] - cdf_l1[idx-1];
        if (std::abs(den) < 1e-9) return q_values[idx-1];
        double frac = (cdf_target - cdf_l1[idx-1]) / den;
        return q_values[idx-1] + frac * (q_values[idx] - q_values[idx-1]);
    }

private:
    std::vector<double> q_values, cdf_l1, cdf_l2;
    double q_min = 0.0, q_max = 0.0, dq = 0.0;
};

// --- Global variables for pre-loaded tuning data ---
namespace {
    // OPTIMIZATION: Use an integer tuple as the map key for high performance.
    // Key: <chain_idx, nucleus_idx, detector_idx, ekBin>
    std::map<std::tuple<int, int, int, int>, LookupTable> g_chargeTuningData_EGE;
    std::once_flag g_chargeTuningInitFlag;

    // OPTIMIZATION: Define name-to-index mappings for fast string-to-int conversion.
    const std::vector<std::string> CHAIN_NAMES = {"UnbiasedL1Inner", "L1Inner"};
    const std::vector<std::string> NUCLEUS_NAMES = {"Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
    const std::vector<std::string> DETECTOR_NAMES = {"TOF", "NaF", "AGL"};
    
    // Reverse maps for fast lookup, filled during initialization.
    std::map<std::string, int> g_chain_map, g_nucleus_map, g_detector_map;
}

// --- Function Implementations ---

void initChargeTuning(const std::string& filename) {
    std::call_once(g_chargeTuningInitFlag, [&](){
        std::cout << "[INFO] Initializing Charge Tuning System (EGE model, optimized)..." << std::endl;
        
        // Populate the reverse maps for string-to-index conversion.
        for(size_t i=0; i<CHAIN_NAMES.size(); ++i) g_chain_map[CHAIN_NAMES[i]] = i;
        for(size_t i=0; i<NUCLEUS_NAMES.size(); ++i) g_nucleus_map[NUCLEUS_NAMES[i]] = i;
        for(size_t i=0; i<DETECTOR_NAMES.size(); ++i) g_detector_map[DETECTOR_NAMES[i]] = i;

        auto lookupFile = std::unique_ptr<TFile>(TFile::Open(filename.c_str()));
        if (!lookupFile || lookupFile->IsZombie()) {
            std::cerr << "[CRITICAL ERROR] Cannot open charge tuning lookup file: " << filename << std::endl;
            return;
        }

        // Temporarily store pointers to related TVectorD objects.
        std::map<std::string, std::array<TVectorD*, 3>> temp_vectors;
        
        // Iterate through all keys in the ROOT file.
        TIter next(lookupFile->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)next())) {
            // Skip objects that are not TVectorD.
            if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TVectorD")) continue;
            
            TString name = key->GetName();
            // Filter for EGE model vectors only.
            if (!name.EndsWith("_EGE_q") && !name.EndsWith("_EGE_cdf_l1") && !name.EndsWith("_EGE_cdf_l2")) continue;
            
            // Read the object from the file.
            TVectorD* vec = (TVectorD*)key->ReadObj();
            
            // Extract the base name and determine the vector type (q, cdf1, or cdf2).
            TString baseName = name;
            int index = -1;
            if (name.EndsWith("_EGE_q")) { baseName.ReplaceAll("_EGE_q", ""); index = 0; } 
            else if (name.EndsWith("_EGE_cdf_l1")) { baseName.ReplaceAll("_EGE_cdf_l1", ""); index = 1; } 
            else if (name.EndsWith("_EGE_cdf_l2")) { baseName.ReplaceAll("_EGE_cdf_l2", ""); index = 2; }
            
            // Group the three related vectors by their base name.
            if (index != -1) {
                if (temp_vectors.find(baseName.Data()) == temp_vectors.end()) temp_vectors[baseName.Data()] = {nullptr, nullptr, nullptr};
                temp_vectors[baseName.Data()][index] = vec;
            }
        }

        size_t loaded_count = 0;
        // Process the grouped vectors.
        for (auto const& [baseNameStr, vecs] : temp_vectors) {
            TString baseName(baseNameStr);
            TObjArray* tokens = baseName.Tokenize("_");
            if (tokens->GetEntries() >= 4) {
                // Parse the name to get chain, nucleus, detector, and bin index.
                std::string chain_str = ((TObjString*)tokens->At(0))->GetString().Data();
                std::string nucleus_str = ((TObjString*)tokens->At(1))->GetString().Data();
                std::string detector_str = ((TObjString*)tokens->At(2))->GetString().Data();
                TString bin_str = ((TObjString*)tokens->At(3))->GetString();
                bin_str.ReplaceAll("bin", "");
                int ekBin = bin_str.Atoi();

                // Convert names to integer indices using the pre-filled maps.
                if (g_chain_map.count(chain_str) && g_nucleus_map.count(nucleus_str) && g_detector_map.count(detector_str)) {
                    int chain_idx = g_chain_map[chain_str];
                    int nucleus_idx = g_nucleus_map[nucleus_str];
                    int detector_idx = g_detector_map[detector_str];
                    std::tuple<int, int, int, int> key = {chain_idx, nucleus_idx, detector_idx, ekBin};
                    
                    // Create and load the LookupTable object.
                    LookupTable pdata;
                    pdata.loadFromVectors(vecs[0], vecs[1], vecs[2]);
                    if (pdata.isValid) {
                        // Move the loaded data into the global map.
                        g_chargeTuningData_EGE[key] = std::move(pdata);
                        loaded_count++;
                    }
                }
            }
            delete tokens;
            // Free the memory of the TVectorD objects read from the file.
            delete vecs[0]; delete vecs[1]; delete vecs[2];
        }
        
        std::cout << "[INFO] Charge Tuning System initialized. Loaded " << loaded_count << " valid EGE lookup tables using integer keys." << std::endl;
    });
}

double tuneL2Charge(
    const std::string& chain, 
    const std::string& nucleusName, 
    const std::string& detectorName, 
    int ekBin, 
    double q_l2
) {
    if (ekBin < 0) {
        // std::cerr << "[DBG] L2Q Fail: Invalid ekBin=" << ekBin << std::endl;
        return q_l2;
    }

    // OPTIMIZATION: Convert input strings to integer indices for fast lookup.
    auto it_chain = g_chain_map.find(chain);
    auto it_nuc = g_nucleus_map.find(nucleusName);
    auto it_det = g_detector_map.find(detectorName);

    // If any name is not found, it's an invalid call. Return original value.
    if (it_chain == g_chain_map.end() || it_nuc == g_nucleus_map.end() || it_det == g_detector_map.end()) {
         std::cerr << "[DBG] L2Q Fail: Invalid name. chain=" << chain << " nuc=" << nucleusName << " det=" << detectorName << std::endl;
        return q_l2;
    }

    // Construct the integer-based key.
    std::tuple<int, int, int, int> key = {it_chain->second, it_nuc->second, it_det->second, ekBin};

    // Perform the fast lookup in the map.
    auto it = g_chargeTuningData_EGE.find(key);
    if (it == g_chargeTuningData_EGE.end()) {
         //std::cerr << "[DBG] L2Q Fail: Table not found. chain=" << chain << " nuc=" << nucleusName << " det=" << detectorName << " ekBin=" << ekBin << std::endl;
        return q_l2; // No table found for this combination.
    }
    
    // Use the found lookup table to perform the tuning.
    const LookupTable& pdata = it->second;
    double cdf_l2 = pdata.getCDF_L2(q_l2);
    if (cdf_l2 < 0.0) {
        //std::cerr << "[DBG] L2Q Fail: q_l2 out of range. q_l2=" << q_l2 << std::endl;
        return q_l2; // q_l2 is out of the table's range.
    }
    
    double q_tuned = pdata.getInvCDF_L1(cdf_l2);
    if (q_tuned < 0) {
         std::cerr << "[DBG] L2Q Fail: InvCDF failed. cdf_l2=" << cdf_l2 << std::endl;
        return q_l2; // Inverse CDF calculation failed.
    }
    //std::cout<<q_l2<<" tune:"<<q_tuned<<std::endl;

    return q_tuned;
}

} // namespace Tools
} // namespace AMS_Iso
