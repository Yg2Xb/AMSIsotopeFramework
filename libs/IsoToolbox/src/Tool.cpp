/***********************************************************
 * File: Tools.cpp
 *
 * Modern C++ implementation file for AMS analysis tools.
 *
 * History:
 * 20250831 - Ported and refactored by Gemini
 ***********************************************************/

#include <IsoToolbox/Tools.h>
#include <IsoToolbox/PhysicsConstants.h>
#include <map>
#include <numeric>

namespace IsoToolbox {
namespace Tools {

// MC REWEIGHTING HELPERS
namespace {
    static const double geneRig_low = 1.0;
    static const double geneRig_up = 2000.0;
    
    TF1& getMCFunction() {
        static TF1 f("f_MC", "1/x", geneRig_low, geneRig_up);
        return f;
    }

    TF1& getReweightFunction() {
        static TF1 f("f_Reweight", "pow(x, -2.7)", geneRig_low, geneRig_up);
        return f;
    }

    double getMCNorm() {
        static const double norm = getMCFunction().Integral(geneRig_low, geneRig_up);
        return norm;
    }

    double getReweightNorm() {
        static const double norm = getReweightFunction().Integral(geneRig_low, geneRig_up);
        return norm;
    }
}

double calculateWeight(double mmom, int charge, bool isISS) {
    if (isISS) return 1.0;
    if (charge == 0) {
        std::cerr << "Charge is zero, setting weight to 0." << std::endl;
        return 0.0;
    }
    double geneRig = mmom / static_cast<double>(charge);
    if (geneRig < geneRig_low || geneRig > geneRig_up) return 0.0;

    double mc_eval = getMCFunction().Eval(geneRig);
    double reweight_eval = getReweightFunction().Eval(geneRig);

    if (mc_eval == 0) {
        std::cerr << "MC function evaluates to zero at rigidity " << geneRig << ", setting weight to 0." << std::endl;
        return 0.0;
    }

    return (reweight_eval / getReweightNorm()) / (mc_eval / getMCNorm());
}


// PHYSICS & KINEMATICS IMPLEMENTATION
double betaToKineticEnergy(double beta) {
    if (beta <= 0.0 || beta >= 1.0) return -9.0;
    double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    return (gamma - 1.0) * PhysicsConstants::MASS_UNIT;
}

double kineticEnergyToBeta(double kineticEnergy) {
    if (kineticEnergy < 0.0) return -9.0;
    double gamma = kineticEnergy / PhysicsConstants::MASS_UNIT + 1.0;
    return std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

double rigidityToBeta(double rigidity, int charge, int mass, bool isElectron) {
    if (!isElectron && mass < charge) {
        throw std::invalid_argument("Invalid charge/mass combination");
    }

    if (isElectron) {
        constexpr double ELECTRON_MASS = 0.000511; // GeV
        double beta = rigidity / std::sqrt(ELECTRON_MASS * ELECTRON_MASS + rigidity * rigidity);
        return (charge == -1) ? beta : -beta;
    } else {
        double particleMass = static_cast<double>(mass) * PhysicsConstants::MASS_UNIT;
        double charge_d = static_cast<double>(charge);
        return (rigidity * charge_d) / std::sqrt(particleMass * particleMass + rigidity * rigidity * charge_d * charge_d);
    }
}

double betaToRigidity(double beta, int charge, int mass, bool isElectron) {
    if (beta <= 0.0 || beta >= 1.0) return -100000.0;

    if (isElectron) {
        constexpr double ELECTRON_MASS = 0.000511; // GeV
        double rigidity = beta * ELECTRON_MASS / std::sqrt(1.0 - beta * beta);
        return (charge == -1) ? rigidity : -rigidity;
    } else {
        double particleMass = static_cast<double>(mass) * PhysicsConstants::MASS_UNIT;
        return (beta * particleMass) / (std::sqrt(1.0 - beta * beta) * static_cast<double>(charge));
    }
}

double kineticEnergyToRigidity(double ek_per_nucleon, int z, int a) {
    if (ek_per_nucleon < 0.0 || z == 0) return -100000.0;
    double factor = (static_cast<double>(a) * PhysicsConstants::MASS_UNIT) / static_cast<double>(z);
    double term = std::pow(ek_per_nucleon / PhysicsConstants::MASS_UNIT + 1.0, 2) - 1.0;
    if (term < 0) return -100000.0;
    return factor * std::sqrt(term);
}

double rigidityToKineticEnergy(double rig_gv, int z, int a) {
    if (rig_gv <= 0.0 || z == 0) return -9.0;
    double factor = (static_cast<double>(a) * PhysicsConstants::MASS_UNIT) / static_cast<double>(z);
    double term = std::pow(rig_gv / factor, 2);
    return PhysicsConstants::MASS_UNIT * (std::sqrt(1.0 + term) - 1.0);
}

double dR_dEk(double ek_per_nucleon, int z, int a) {
    if (ek_per_nucleon < 0.0 || z == 0) return -100000.0;
    double factor = (static_cast<double>(a) * PhysicsConstants::MASS_UNIT) / static_cast<double>(z);
    double ek_term = ek_per_nucleon / PhysicsConstants::MASS_UNIT + 1.0;
    double denominator = std::sqrt(ek_term * ek_term - 1.0);
    if (denominator == 0) return std::numeric_limits<double>::infinity();
    return factor * ek_term / denominator;
}

double dEk_dR(double rig_gv, int z, int a) {
    if (rig_gv <= 0.0 || z == 0) return -9.0;
    double factor = (static_cast<double>(a) * PhysicsConstants::MASS_UNIT) / static_cast<double>(z);
    double term = std::pow(rig_gv / factor, 2);
    return PhysicsConstants::MASS_UNIT * rig_gv / (factor * factor * std::sqrt(1.0 + term));
}

MassResult calculateMass(double beta, double alpha, double innerRig, int charge) {
    MassResult result{0.0, 0.0, 0.0, 0.0};
    if (!isValidBeta(beta) || innerRig < 0.0 || charge == 0) return result;
    
    result.beta = alpha * beta / std::sqrt(1.0 - beta * beta + std::pow(alpha * beta, 2));
    if (!isValidBeta(result.beta)) return {0.0, 0.0, 0.0, 0.0};
    
    result.gamma = 1.0 / std::sqrt(1.0 - result.beta * result.beta);
    result.ek = (result.gamma - 1.0) * PhysicsConstants::MASS_UNIT;
    result.invMass = (result.beta * result.gamma) / (static_cast<double>(charge) * innerRig);
    
    return result;
}

// EVENT SELECTION & BINNING IMPLEMENTATION
bool isValidBeta(double beta) {
    return beta > 0.0 && beta < 1.0;
}

int findBin(const std::vector<double>& bins, double value) {
    if (bins.empty() || value < bins.front() || value >= bins.back()) return -1;
    auto it = std::upper_bound(bins.begin(), bins.end(), value);
    return std::distance(bins.begin(), it) - 1;
}

int findRigidityLowerBound(double rigidity, double *xAxis, int range) {
    int pos = -1;
    for (int i = 0; i < range; ++i) {
        if (rigidity >= xAxis[i] && rigidity < xAxis[i + 1]) {
            pos = i;
            break;
        }
    }
    return pos;
}

// GEOMETRY & UTILITIES IMPLEMENTATION
std::string convertElementName(const std::string& fullName, bool isISS) {
    if (isISS) {
        return fullName.substr(0, 3);
    } else {
        static const std::map<std::string, std::string> mcNameMap = {
            {"Lithium", "Li"}, {"Beryllium", "Be"}, {"Boron", "B"},
            {"Carbon", "C"}, {"Nitrogen", "N"}, {"Oxygen", "O"}
        };
        auto it = mcNameMap.find(fullName);
        return (it != mcNameMap.end()) ? it->second : fullName;
    }
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
    return validCount > 0 ? sum / static_cast<double>(validCount) : -1000000.0;
}

void modifyPositionByZ(double targetZ, std::array<double, 3>& position, double theta, double phi) {
    if (std::abs(theta) > M_PI || std::abs(phi) > 2.0 * M_PI) {
        throw std::invalid_argument("Invalid angle input (theta or phi)");
    }
    double deltaZ = targetZ - position[2];
    position[0] += deltaZ * std::tan(theta) * std::cos(phi);
    position[1] += deltaZ * std::tan(theta) * std::sin(phi);
    position[2] = targetZ;
}

std::optional<Point2D> calculateXYAtZ(const Float_t positions[9][3], const Float_t directions[9][3], double zpl, int trackIndex) {
    if (trackIndex < 0) return std::nullopt;

    int lmin = 4, lmax = 4;
    double ww = 1.0;

    if (zpl >= positions[0][2]) {
        lmin = lmax = 0;
    } else if (zpl <= positions[8][2]) {
        lmin = lmax = 8;
    } else {
        for (int ilay = 0; ilay < 8; ++ilay) {
            if (zpl <= positions[ilay][2] && zpl >= positions[ilay + 1][2]) {
                double dz_layer = positions[ilay][2] - positions[ilay + 1][2];
                if (std::abs(dz_layer) < std::numeric_limits<double>::epsilon()) continue;
                ww = (zpl - positions[ilay + 1][2]) / dz_layer;
                lmin = ilay;
                lmax = ilay + 1;
                break;
            }
        }
    }

    if (!(ww >= 0.0 && ww <= 1.0)) {
        std::cerr << "Error_Interpolate_Weight=" << ww << std::endl;
        return std::nullopt;
    }

    Point2D result;
    double point[2];
    for (int ixy = 0; ixy < 2; ++ixy) {
        double xymin = positions[lmin][ixy] + (zpl - positions[lmin][2]) * directions[lmin][ixy] / directions[lmin][2];
        double xymax = positions[lmax][ixy] + (zpl - positions[lmax][2]) * directions[lmax][ixy] / directions[lmax][2];
        point[ixy] = ww * xymin + (1.0 - ww) * xymax;
    }
    result.x = point[0];
    result.y = point[1];
    return result;
}

std::optional<Point2D> calculateXYAtNaFZ(const std::array<double, 3>& pos1, const std::array<double, 3>& pos2, double nafZ) {
    double dz = pos2[2] - pos1[2];
    if (std::abs(dz) < std::numeric_limits<double>::epsilon()) return std::nullopt;
    double t = (nafZ - pos1[2]) / dz;
    return Point2D{pos1[0] + t * (pos2[0] - pos1[0]), pos1[1] + t * (pos2[1] - pos1[1])};
}

double getCurrentRSS_MB() {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    #ifdef __APPLE__
        return static_cast<double>(rusage.ru_maxrss) / (1024.0 * 1024.0);
    #else
        return static_cast<double>(rusage.ru_maxrss) / 1024.0;
    #endif
}

std::unique_ptr<TH1D> getEnergySlice(int plusIbin, int Nrebin, TH2D* h2d, int binIndex, const std::string& name) {
    // This function references Binning::NarrowBins, which should be defined in PhysicsConstants.h
    // For now, this is a placeholder for compilation.
    // Example: const std::vector<double>& narrowBins = IsoToolbox::PhysicsConstants::Binning::NarrowBins;
    // if (!h2d || binIndex < 0 || binIndex >= narrowBins.size()-1) return nullptr;
    // double eMin = narrowBins[binIndex];
    // double eMax = narrowBins[binIndex + plusIbin];
    // Since the constant is not yet defined, this function cannot be fully implemented.
    // Returning a nullptr to indicate it's not ready.
    std::cerr << "Warning: getEnergySlice cannot be used until Binning constants are defined in PhysicsConstants.h" << std::endl;
    return nullptr;
}


// FITTING FUNCTIONS IMPLEMENTATION
Double_t landauGaussianFit(Double_t *x, Double_t *par) {
    Double_t invsq2pi = 0.3989422804014;
    Double_t mpshift = -0.22278298;
    Double_t np = 1500.0;
    Double_t sc = 5.0;
    Double_t xx, mpc, fland, sum = 0.0, xlow, xupp, step;
    mpc = par[1] - mpshift * par[0];
    xlow = x[0] * par[4] - sc * par[3];
    xupp = x[0] * par[4] + sc * par[3];
    step = (xupp - xlow) / np;
    for (double i = 1.0; i <= np / 2.0; i++) {
        xx = xlow + (i - 0.5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0] * par[4], xx, par[3]);
        xx = xupp - (i - 0.5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0] * par[4], xx, par[3]);
    }
    return (par[2] * step * sum * invsq2pi / par[3]);
}

double extendedExponentialGaussianFit(double *x, double *p) {
    constexpr double sqrtPiOver2 = 1.2533141373;
    constexpr double sqrt2 = 1.4142135624;
    const double &x0 = p[0], &sigmaL = p[1], &alphaL = p[2], &sigmaR = p[3], &alphaR = p[4], &norm = p[5], &xmin = p[6], &xmax = p[7];
    double t = (x[0] - x0) / (x[0] < x0 ? sigmaL : sigmaR);
    double v = 0;
    if (t < -alphaL) {
        v = std::exp(0.5 * alphaL * alphaL + alphaL * t);
    } else if (t <= alphaR) {
        v = std::exp(-0.5 * t * t);
    } else {
        v = std::exp(0.5 * alphaR * alphaR - alphaR * t);
    }
    double tmin = (xmin - x0) / (xmin < x0 ? sigmaL : sigmaR);
    double tmax = (xmax - x0) / (xmax < x0 ? sigmaL : sigmaR);
    double sum = 0;
    if (tmin < -alphaL) {
        double a = 0.5 * alphaL * alphaL;
        double lv = tmin, uv = std::min(tmax, -alphaL);
        sum += (sigmaL / alphaL) * (std::exp(a + alphaL * uv) - std::exp(a + alphaL * lv));
    }
    if (tmax > alphaR) {
        double a = 0.5 * alphaR * alphaR;
        double lv = std::max(tmin, alphaR), uv = tmax;
        sum += (sigmaR / alphaR) * (std::exp(a - alphaR * lv) - std::exp(a - alphaR * uv));
    }
    if (tmin < alphaR && tmax > -alphaL) {
        double sigmaMin = (tmin < 0.0) ? sigmaL : sigmaR;
        double sigmaMax = (tmax < 0.0) ? sigmaL : sigmaR;
        sum += sqrtPiOver2 * (sigmaMax * std::erf(std::min(tmax, alphaR) / sqrt2) - sigmaMin * std::erf(std::max(tmin, -alphaL) / sqrt2));
    }
    if (sum == 0) return 0;
    return norm * (v / sum);
}

// PLOTTING & ANALYSIS IMPLEMENTATION
void calculatePull(RooPlot* frame, TGraphErrors* pullGraph, double fit_min, double fit_max) {
    if (!frame || !pullGraph) return;
    RooHist* dataHist = (RooHist*)frame->findObject("data_hist");
    RooCurve* totalCurve = (RooCurve*)frame->findObject("total_pdf");
    if (!dataHist || !totalCurve) {
        std::cout << "Warning: Cannot find data or model curve for pull calculation" << std::endl;
        return;
    }
    int pointIndex = 0;
    for (int i = 0; i < dataHist->GetN(); ++i) {
        double x_val, y_data;
        dataHist->GetPoint(i, x_val, y_data);
        if (x_val >= fit_min && x_val <= fit_max && y_data > 0) {
            double y_model = totalCurve->Eval(x_val);
            if (y_model > 0) {
                double pull_val = (y_data - y_model) / sqrt(y_data);
                if (std::isfinite(pull_val) && std::abs(pull_val) < 10) {
                    pullGraph->SetPoint(pointIndex, x_val, pull_val);
                    pullGraph->SetPointError(pointIndex, 0, 0);
                    pointIndex++;
                }
            }
        }
    }
}

void setupPullPlot(TGraphErrors* pullGraph, double fit_min, double fit_max) {
    if (!pullGraph) return;
    pullGraph->SetTitle(";L1 Charge;Pull");
    pullGraph->SetMarkerColor(kBlack);
    pullGraph->SetMarkerStyle(20);
    pullGraph->SetMarkerSize(0.8);
    pullGraph->SetLineColor(kBlack);
    pullGraph->GetXaxis()->SetLabelSize(0.12);
    pullGraph->GetYaxis()->SetLabelSize(0.10);
    pullGraph->GetXaxis()->SetTitleSize(0.14);
    pullGraph->GetYaxis()->SetTitleSize(0.14);
    pullGraph->GetYaxis()->SetTitleOffset(0.4);
    pullGraph->GetXaxis()->SetTitleOffset(0.9);
    pullGraph->GetXaxis()->SetRangeUser(fit_min, fit_max);
    pullGraph->GetYaxis()->SetRangeUser(-6, 6);
    pullGraph->GetYaxis()->SetNdivisions(505);
}

std::unique_ptr<TLegend> createLegend() {
    auto legend = std::make_unique<TLegend>(0.72, 0.72, 0.88, 0.86);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);
    return legend;
}

std::unique_ptr<TH1D> extendHistogram(const TH1D* hist, double new_min, double new_max) {
    if (!hist) return nullptr;
    
    double curr_min = hist->GetXaxis()->GetXmin();
    double curr_max = hist->GetXaxis()->GetXmax();
    int curr_nbins = hist->GetNbinsX();
    double bin_width = (curr_max - curr_min) / static_cast<double>(curr_nbins);
    
    if (curr_min <= new_min && curr_max >= new_max) {
        return std::unique_ptr<TH1D>(static_cast<TH1D*>(hist->Clone((std::string(hist->GetName()) + "_clone").c_str())));
    }
    
    int new_nbins = static_cast<int>((new_max - new_min) / bin_width);
    
    auto new_hist = std::make_unique<TH1D>(
        (std::string(hist->GetName()) + "_ext").c_str(),
        hist->GetTitle(),
        new_nbins, new_min, new_max
    );
    
    new_hist->SetLineColor(hist->GetLineColor());
    new_hist->SetLineWidth(hist->GetLineWidth());
    new_hist->SetLineStyle(hist->GetLineStyle());
    new_hist->SetFillColor(hist->GetFillColor());
    new_hist->SetFillStyle(hist->GetFillStyle());
    new_hist->SetMarkerColor(hist->GetMarkerColor());
    new_hist->SetMarkerStyle(hist->GetMarkerStyle());
    new_hist->SetMarkerSize(hist->GetMarkerSize());
    
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double x = hist->GetBinCenter(i);
        double content = hist->GetBinContent(i);
        if (content > 0) {
            int new_bin = new_hist->FindBin(x);
            new_hist->SetBinContent(new_bin, content);
            new_hist->SetBinError(new_bin, hist->GetBinError(i));
        }
    }
    
    new_hist->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
    new_hist->GetXaxis()->SetRangeUser(new_min, new_max);
    new_hist->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
    new_hist->SetDirectory(nullptr);
    
    return new_hist;
}

} // namespace Tools
} // namespace IsoToolbox
