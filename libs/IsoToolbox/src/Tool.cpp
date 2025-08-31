#include "IsoToolbox/Tools.h"
#include "IsoToolbox/PhysicsConstants.h"
#include <stdexcept>
#include <iostream>
#include <sys/resource.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLine.h>

namespace IsoToolbox {
namespace Tools {

    // ===================================================================
    // 2. Core Physics Conversions
    // ===================================================================
    double RigidityToMomentum(double rigidity, int charge) {
        return std::abs(charge * rigidity);
    }

    double MomentumToBeta(double momentum, int mass) {
        if (mass <= 0 || momentum < 0) return 0.0;
        double M_particle = Constants::AMU_GEV * mass;
        return momentum / TMath::Sqrt(momentum * momentum + M_particle * M_particle);
    }
    
    double BetaToMomentum(double beta, int mass) {
        if (mass <= 0 || beta <= 0 || beta >= 1) return 0.0;
        double M_particle = Constants::AMU_GEV * mass;
        return M_particle * beta / TMath::Sqrt(1 - beta * beta);
    }
    
    double MomentumToEkPerNucleon(double momentum, int mass) {
        if (mass <= 0 || momentum < 0) return 0.0;
        double M_particle = Constants::AMU_GEV * mass;
        double E_total = TMath::Sqrt(momentum * momentum + M_particle * M_particle);
        return (E_total - M_particle) / mass;
    }
    
    double EkPerNucleonToMomentum(double ek_per_nucleon, int mass) {
        if (mass <= 0 || ek_per_nucleon <= 0) return 0.0;
        double M_particle = Constants::AMU_GEV * mass;
        double E_total = ek_per_nucleon * mass + M_particle;
        return TMath::Sqrt(E_total * E_total - M_particle * M_particle);
    }

    double RigidityToEkPerNucleon(double rigidity, int mass, int charge) {
        return MomentumToEkPerNucleon(RigidityToMomentum(rigidity, charge), mass);
    }

    double EkPerNucleonToRigidity(double ek_per_nucleon, int mass, int charge) {
        if (charge == 0) return 0.0;
        return EkPerNucleonToMomentum(ek_per_nucleon, mass) / std::abs(charge);
    }

    double BetaToEkPerNucleon(double beta, int mass) {
        return MomentumToEkPerNucleon(BetaToMomentum(beta, mass), mass);
    }

    double EkPerNucleonToBeta(double ek_per_nucleon, int mass) {
        return MomentumToBeta(EkPerNucleonToMomentum(ek_per_nucleon, mass), mass);
    }
    
    double RigidityToBeta(double rigidity, int mass, int charge) {
        return MomentumToBeta(RigidityToMomentum(rigidity, charge), mass);
    }

    double BetaToRigidity(double beta, int mass, int charge) {
        if (charge == 0) return 0.0;
        return BetaToMomentum(beta, mass) / std::abs(charge);
    }

    double dR_dEk(double ek_per_nucleon, int mass, int charge) {
        if (charge == 0 || mass == 0 || ek_per_nucleon <= 0) return 0.0;
        double M_particle = Constants::AMU_GEV * mass;
        double E_kinetic = ek_per_nucleon * mass;
        double E_total = E_kinetic + M_particle;
        return (E_total * E_total) / (std::abs(charge) * mass * TMath::Sqrt(E_total * E_total - M_particle * M_particle));
    }

    double dEk_dR(double rigidity, int mass, int charge) {
        double R_dEk = dR_dEk(RigidityToEkPerNucleon(rigidity, mass, charge), mass, charge);
        return (R_dEk == 0) ? 0.0 : 1.0 / R_dEk;
    }

    double CalculateMass(double rigidity, double beta, int charge) {
        if (beta <= 0 || beta >= 1) return 0.0;
        return (std::abs(charge) * rigidity / beta) * TMath::Sqrt(1 - beta * beta);
    }

    // ===================================================================
    // 3. MC & Weighting Functions
    // ===================================================================
    namespace { // Anonymous namespace for internal linkage
        TF1* gMCNormFunc = nullptr;
        TF1* gReweightFunc = nullptr;
    }

    TF1* GetMCNormFunction() {
        if (!gMCNormFunc) {
            gMCNormFunc = new TF1("f_MC_norm", "1/x", 1.0, 3300.0);
        }
        return gMCNormFunc;
    }
    
    TF1* GetReweightFunction() {
        if (!gReweightFunc) {
            gReweightFunc = new TF1("f_Reweight_norm", "pow(x, -2.7)", 1.0, 3300.0);
        }
        return gReweightFunc;
    }

    double CalculateWeight(double rigidity, bool is_iss, int charge) {
        if (is_iss || rigidity <= 0 || charge == 0) return 1.0;
        
        static const double mc_norm = GetMCNormFunction()->Integral(1.0, 3300.0);
        static const double reweight_norm = GetReweightFunction()->Integral(1.0, 3300.0);
        
        if (mc_norm == 0 || reweight_norm == 0) return 1.0;

        double weight = (GetReweightFunction()->Eval(rigidity) / reweight_norm) / (GetMCNormFunction()->Eval(rigidity) / mc_norm);
        return weight;
    }

    // ===================================================================
    // 4. Analysis Helper Functions
    // ===================================================================
    std::string ConvertElementName(const std::string& fullName, bool isISS) {
        if (isISS) return fullName.substr(0, 3);
        static const std::map<std::string, std::string> mcNameMap = {
            {"Lithium", "Li"}, {"Beryllium", "Be"}, {"Boron", "B"},
            {"Carbon", "C"}, {"Nitrogen", "N"}, {"Oxygen", "O"}
        };
        auto it = mcNameMap.find(fullName);
        return (it != mcNameMap.end()) ? it->second : fullName;
    }
    
    bool IsValidBeta(double beta) {
        return (beta > 0 && beta < 1);
    }
    
    int FindBin(const std::vector<double>& bins, double value) {
        if (value < bins.front() || value >= bins.back()) return -1;
        auto it = std::upper_bound(bins.begin(), bins.end(), value);
        return std::distance(bins.begin(), it) - 1;
    }

    bool IsBeyondCutoff(double beta, double theta_m) {
        if (theta_m < 0.2) return beta > 0.953;
        if (theta_m < 0.4) return beta > 0.96;
        return beta > 0.97;
    }
    
    bool IsBeyondCutoff(double beta, double rigidity, bool is_rich_naf, bool is_agl) {
        if (beta < 0) return false;
        if (rigidity < 2.15) return true;
        if (is_rich_naf && beta > 0.96) return true;
        if (is_agl && beta > 0.99) return true;
        return false;
    }

    std::optional<Point2D> CalculateXYAtZ(double z_target, double z_pos, const double pos[3], const double dir[3]) {
        if (std::abs(dir[2]) < 1e-9) return std::nullopt;
        double t = (z_target - pos[2]) / dir[2];
        Point2D result;
        result.x = pos[0] + t * dir[0];
        result.y = pos[1] + t * dir[1];
        return result;
    }

    void CorrectCalibrationBias(double& beta, int charge, const TF1& func) {
        beta -= func.Eval(charge);
    }
    
    void CorrectCalibrationBiasInData(double& beta, const TF1& func) {
        beta -= func.GetParameter(0);
    }

    void setCutStatus(std::vector<bool>& status, int index, bool value) {
        if (index < status.size()) {
            status[index] = value;
        }
    }
    
    double calculateAverage(const std::vector<double>& vec) {
        if (vec.empty()) return 0.0;
        return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
    }
    
    // ===================================================================
    // 5. Plotting & Diagnostic Functions
    // ===================================================================
    TGraphErrors* calculatePull(const TH1* dataHist, const TF1* fitFunc) {
        if (!dataHist || !fitFunc) return nullptr;
        auto pullGraph = new TGraphErrors();
        int pointIndex = 0;
        for (int i = 1; i <= dataHist->GetNbinsX(); ++i) {
            if (dataHist->GetBinContent(i) > 0) {
                double x = dataHist->GetBinCenter(i);
                double y_data = dataHist->GetBinContent(i);
                double y_err = dataHist->GetBinError(i);
                double y_model = fitFunc->Eval(x);
                double pull = (y_data - y_model) / y_err;
                if (std::isfinite(pull)) {
                    pullGraph->SetPoint(pointIndex, x, pull);
                    pullGraph->SetPointError(pointIndex, 0, 0);
                    pointIndex++;
                }
            }
        }
        return pullGraph;
    }

    void setupPullPlot(TGraphErrors* pullGraph, double fit_min, double fit_max) {
        if (!pullGraph) return;
        pullGraph->SetTitle(";L1 Charge;Pull");
        pullGraph->SetMarkerColor(kBlack);
        pullGraph->SetMarkerStyle(20);
        pullGraph->GetXaxis()->SetRangeUser(fit_min, fit_max);
        pullGraph->GetYaxis()->SetRangeUser(-5, 5);
        // Additional styling...
    }
    
    void DrawTheoryLines(int Z) { /* Implementation needed */ }
    double langaufun(double* x, double* par) { /* Implementation needed */ return 0.0;}
    double funcExpGausExp(double* x, double* par) { /* Implementation needed */ return 0.0;}
    
    long getCurrentRSS_MB() {
        struct rusage rusage;
        getrusage(RUSAGE_SELF, &rusage);
        return rusage.ru_maxrss / 1024; // Convert from KB to MB
    }

    TH1F* getEnergySlice(TH2F* h, int bin) {
        if (!h) return nullptr;
        std::string name = std::string(h->GetName()) + "_slice" + std::to_string(bin);
        return (TH1F*)h->ProjectionY(name.c_str(), bin, bin);
    }
    
    TH1* extendHistogram(const TH1* h, int nbins, double xmax) {
        if (!h) return nullptr;
        std::string newName = std::string(h->GetName()) + "_ext";
        TH1* h_ext = new TH1F(newName.c_str(), h->GetTitle(), nbins, h->GetXaxis()->GetXmin(), xmax);
        for(int i = 1; i <= h->GetNbinsX(); ++i) {
            h_ext->SetBinContent(i, h->GetBinContent(i));
        }
        return h_ext;
    }
    
    TLegend* createLegend(double x1, double y1, double x2, double y2) {
        return new TLegend(x1, y1, x2, y2);
    }

} // namespace Tools
} // namespace IsoToolbox