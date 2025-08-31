#pragma once

#include <string>
#include <vector>
#include <array>
#include <optional>
#include <numeric>
#include <algorithm>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>

namespace IsoToolbox {
namespace Tools {

    // ===================================================================
    // 1. Core Data Structures
    // ===================================================================
    struct Point2D { double x = 0.0; double y = 0.0; };

    template<size_t N>
    struct CutResult {
        bool pass_all = false;
        std::array<bool, N> details{};
        CutResult() = default;
        explicit CutResult(const std::array<bool, N>& cuts) : details(cuts) {
            pass_all = std::all_of(details.begin(), details.end(), [](bool b){ return b; });
        }
    };

    // ===================================================================
    // 2. Core Physics Conversions
    // ===================================================================
    double RigidityToMomentum(double rigidity, int charge);
    double MomentumToBeta(double momentum, int mass);
    double BetaToMomentum(double beta, int mass);
    
    // Replaces all rigidity/beta <-> Ek functions for consistency
    double MomentumToEkPerNucleon(double momentum, int mass);
    double EkPerNucleonToMomentum(double ek_per_nucleon, int mass);

    // Convenience wrappers
    double RigidityToEkPerNucleon(double rigidity, int mass, int charge);
    double EkPerNucleonToRigidity(double ek_per_nucleon, int mass, int charge);
    double BetaToEkPerNucleon(double beta, int mass);
    double EkPerNucleonToBeta(double ek_per_nucleon, int mass);
    double RigidityToBeta(double rigidity, int mass, int charge);
    double BetaToRigidity(double beta, int mass, int charge);
    
    double dR_dEk(double ek_per_nucleon, int mass, int charge); // Derivative
    double dEk_dR(double rigidity, int mass, int charge);     // Derivative

    double CalculateMass(double rigidity, double beta, int charge);

    // ===================================================================
    // 3. MC & Weighting Functions
    // ===================================================================
    TF1* GetMCNormFunction();
    TF1* GetReweightFunction();
    double CalculateWeight(double rigidity, bool is_iss, int charge);

    // ===================================================================
    // 4. Analysis Helper Functions
    // ===================================================================
    std::string ConvertElementName(const std::string& fullName, bool isISS);
    bool IsValidBeta(double beta);
    int FindBin(const std::vector<double>& bins, double value);
    bool IsBeyondCutoff(double beta, double theta_m);
    bool IsBeyondCutoff(double beta, double rigidity, bool is_rich_naf, bool is_agl); // Overloaded version
    std::optional<Point2D> CalculateXYAtZ(double z_target, double z_pos, const double pos[3], const double dir[3]);
    
    // From old code, kept for compatibility
    void CorrectCalibrationBias(double& beta, int charge, const TF1& func);
    void CorrectCalibrationBiasInData(double& beta, const TF1& func);
    void setCutStatus(std::vector<bool>& status, int index, bool value);
    double calculateAverage(const std::vector<double>& vec);
    
    // ===================================================================
    // 5. Plotting & Diagnostic Functions
    // ===================================================================
    TGraphErrors* calculatePull(const TH1* dataHist, const TF1* fitFunc);
    void setupPullPlot(TGraphErrors* pullGraph, double fit_min, double fit_max);
    void DrawTheoryLines(int Z);
    double langaufun(double* x, double* par);
    double funcExpGausExp(double* x, double* par);
    long getCurrentRSS_MB();
    TH1F* getEnergySlice(TH2F* h, int bin);
    TH1* extendHistogram(const TH1* h, int nbins, double xmax);
    TLegend* createLegend(double x1, double y1, double x2, double y2);

} // namespace Tools
} // namespace IsoToolbox