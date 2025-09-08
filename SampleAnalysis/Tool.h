/***********************************************************
 *  File: Tool.cpp
 *
 *  Modern C++ implementation file for AMS analysis tools.
 *
 *  History:
 *    20250131 - created by ZX.Yan
 ***********************************************************/
#include <cmath>
#include <stdexcept>
#include <bitset>
#include <optional>
#include <array>
#include <limits>
#include <algorithm>  // for std::all_of
#include <numeric>  // 为 std::accumulate
#include <iostream>  // 为了使用 std::cerr
#include <ctime>
#include <sys/resource.h>
#include <TF1.h>  // 为了使用 TF1
#include "./basic_var.h"


namespace AMS_Iso {

std::string ConvertElementName(const std::string& fullName, bool isISS) {
    if (isISS) {
        // ISS数据使用前三个字母
        return fullName.substr(0, 3);
    } else {
        // MC数据使用标准缩写
        static const std::map<std::string, std::string> mcNameMap = {
            {"Lithium", "Li"},
            {"Beryllium", "Be"},
            {"Boron", "B"},
            {"Carbon", "C"},
            {"Nitrogen", "N"},
            {"Oxygen", "O"}
        };

        auto it = mcNameMap.find(fullName);
        return (it != mcNameMap.end()) ? it->second : fullName;
    }
}

// 辅助函数：根据核素和质量数获取bin边界
inline const std::array<double, AMS_Iso::Constants::RIGIDITY_BINS + 1>& getKineticEnergyBins(int input_charge, int input_mass) {
    // 在IsotopeData中查找对应电荷的元素
    for (const auto& isotope : AMS_Iso::IsotopeData) {
        if (isotope.charge_ == input_charge) {  // 使用charge_
            // 在该元素的质量数数组中查找对应质量数的位置
            for (size_t i = 0; i < Constants::MAX_ISOTOPES; ++i) {  // 使用固定大小
                if (isotope.mass_[i] == input_mass) {  // 使用mass_
                    return AMS_Iso::Binning::KineticEnergyBins[input_charge-1][i];  // 使用input_charge
                }
            }
            break;  // 如果找到了正确的元素但没有找到对应的质量数，退出循环
        }
    }
    // 如果没有找到对应的bins，返回空数组
    static const std::array<double, AMS_Iso::Constants::RIGIDITY_BINS + 1> empty{};
    return empty;
}

    namespace {
        static const double geneRig_low = 1.0;
        static const double geneRig_up = 2000.0;
    }

    // 将TF1对象改为函数
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

    // 修改calculateWeight函数
    double calculateWeight(double mmom, int charge, bool isISS) {
        if (isISS) return 1.0;
        if (charge == 0) {
            std::cerr << "Charge is zero, setting weight to 0." << std::endl;
            return 0.0;
        }
        double geneRig = mmom / charge;
        return (getReweightFunction().Eval(geneRig) / getReweightNorm()) / 
               (getMCFunction().Eval(geneRig) / getMCNorm());
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
    return (gamma - 1.0) * Mass_Unit;
}

double kineticEnergyToBeta(double kineticEnergy) {
    if (kineticEnergy < 0.0) return -9.0;
    
    double gamma = kineticEnergy / Mass_Unit + 1.0;
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
        double particleMass = mass * Mass_Unit;
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
        double particleMass = mass * Mass_Unit;
        return beta * particleMass / (std::sqrt(1.0 - beta * beta) * charge);
    }
}

double kineticEnergyToRigidity(double ek_per_nucleon, int z, int a) {
    if (ek_per_nucleon < 0.0 || z == 0) return -100000.0;
    
    double factor = (a * Mass_Unit) / z;
    double term = std::pow(ek_per_nucleon / Mass_Unit + 1, 2) - 1;
    return factor * std::sqrt(term);
}

double dR_dEk(double ek_per_nucleon, int z, int a) {
    if (ek_per_nucleon < 0.0 || z == 0) return -100000.0;
    
    double factor = (a * Mass_Unit) / z;
    double ek_term = ek_per_nucleon / Mass_Unit + 1;
    return factor * ek_term / std::sqrt(ek_term * ek_term - 1);
}

double rigidityToKineticEnergy(double rig_gv, int z, int a) {
    if (rig_gv <= 0.0 || z == 0) return -9.0;
    
    double factor = (a * Mass_Unit) / z;
    double term = std::pow(rig_gv / factor, 2);
    return Mass_Unit * (std::sqrt(1 + term) - 1);
}

double dEk_dR(double rig_gv, int z, int a) {
    if (rig_gv <= 0.0 || z == 0) return -9.0;
    
    double factor = (a * Mass_Unit) / z;
    double term = std::pow(rig_gv / factor, 2);
    return Mass_Unit * rig_gv / (factor * factor * std::sqrt(1 + term));
}

bool isValidBeta(double beta) {
    return beta > 0 && beta < 1;
}

// Optimized findBin function
int findBin(const std::vector<double>& Rbins_beta, double beta) {
    if (!isValidBeta(beta)) return -1;
    if (beta >= Rbins_beta.back()) {
        return -1;
    }

    auto it = std::lower_bound(Rbins_beta.begin(), Rbins_beta.end(), beta);

    if (it == Rbins_beta.begin()) {
        return -1;
    }
    
    return std::distance(Rbins_beta.begin(), it) - 1;
}

int findBin2(std::vector<double>& bins, double value) {
    if (value < bins.front() || value >= bins.back()) return -1;
    for (int i = 0; i < bins.size()-1; ++i) {
        if (value >= bins[i] && value < bins[i+1]) return i;
    }
    return -1;
}

int findEkBin(double beta, int charge, int mass) {
    if (!isValidBeta(beta)) return -1;
    
    double ek = betaToKineticEnergy(beta);
    const auto& bins = Binning::KineticEnergyBins[charge-1][mass-1];
    
    if (ek < bins[0] || ek >= bins[Constants::RIGIDITY_BINS]) {
        return -1;
    }
    
    for (int i = 0; i < Constants::RIGIDITY_BINS; ++i) {
        if (ek >= bins[i] && ek < bins[i + 1]) {
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

bool isBeyondCutoff2(double beta_low, double cutoffRig, double safetyFactor, int charge, int UseMass,  bool isMC) {
    if (isMC) return true;
    if (!isValidBeta(beta_low)) return false;
        
    double rig_beta = betaToRigidity(beta_low, charge, UseMass, false);
        
    return rig_beta > 1.2 * cutoffRig;
}

bool selectPure(double beta, double cutoffRig, double safetyFactor, int charge, double InnerRig,  bool isMC) {
    if (isMC) return true;
    if (!isValidBeta(beta)) return false;
        
    double cutoffBeta10 = rigidityToBeta(cutoffRig, charge, 10, false);
    double cutoffBeta9 = rigidityToBeta(cutoffRig, charge, 9, false);
    if (!isValidBeta(cutoffBeta10)) return false;
    if (!isValidBeta(cutoffBeta9)) return false;

    return  (beta > (safetyFactor * cutoffBeta10)) && (beta < (1.0*cutoffBeta9/safetyFactor));// && InnerRig > 1.*cutoffRig;
}

bool selectPure2(double beta1, double beta2, double cutoffRig, double safetyFactor, int charge, double InnerRig,  bool isMC) {
    if (isMC) return true;
    if (!isValidBeta(beta1)) return false;
        
    double cutoffBeta10 = rigidityToBeta(cutoffRig, charge, 10, false);
    double cutoffBeta9 = rigidityToBeta(cutoffRig, charge, 9, false);
    if (!isValidBeta(cutoffBeta10)) return false;
    if (!isValidBeta(cutoffBeta9)) return false;
    
    return  (beta1 > (safetyFactor * cutoffBeta10)) && (beta2 < (cutoffBeta9/safetyFactor));// && InnerRig > 1.*cutoffRig;
}

bool selectPure3(double beta, double cutoffRig, double cutoffRig25, double cutoffRig40, double safetyFactor, int charge, double InnerRig,  bool isMC) {
    if (isMC) return true;
    if (!isValidBeta(beta)) return false;
        
    double cutoffBeta10 = rigidityToBeta(cutoffRig, charge, 10, false);
    double cutoffBeta9 = rigidityToBeta(cutoffRig, charge, 9, false);
    if (!isValidBeta(cutoffBeta10)) return false;
    if (!isValidBeta(cutoffBeta9)) return false;
   
    return  (beta > (safetyFactor * cutoffBeta10)) && (beta < (cutoffBeta9/safetyFactor));// && InnerRig > 1.*cutoffRig;
}

MassResult calculateMass(double beta, double alpha, double innerRig, int charge) {
    MassResult result{0.0, 0.0, 0.0, 0.0};
    if (!isValidBeta(beta) || innerRig < 0.0) return result;
    
    result.beta = alpha * beta / std::sqrt(1 - beta * beta + std::pow(alpha * beta, 2));
    if (!isValidBeta(result.beta)) return {0.0, 0.0, 0.0, 0.0};
    
    result.gamma = 1.0 / std::sqrt(1.0 - result.beta * result.beta);
    result.ek = (result.gamma - 1) * Mass_Unit;
    result.invMass = (result.beta * result.gamma) / (charge * innerRig);
    
    return result;
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

double calculateChi2(RooPlot* frame, const char *DataName, const char* ModelName, double fitRangeLow, double fitRangeUp) {
    RooHist* dataHist_frame = frame->getHist(Form("%s", DataName));
    RooCurve* modelCurve = frame->getCurve(Form("%s", ModelName));
    
    double chi2 = 0;
    int ndf = 0;
    if (dataHist_frame && modelCurve) {
        double x_val, y_data;
        for (int i = 0; i < dataHist_frame->GetN(); ++i) {
            dataHist_frame->GetPoint(i, x_val, y_data);
            
            if (x_val >= fitRangeLow && x_val <= fitRangeUp) {
                if(y_data <= 0) continue;
                ndf++;
                double y_model = modelCurve->Eval(x_val);
                //double dataError = dataHist_frame->GetErrorY(i);
                double dataError = sqrt(y_data);
                
                double pull = (y_data - y_model) / dataError;
                chi2 += pull * pull;
                /* 
                cout << "x: " << x_val 
                     << " data: " << y_data 
                     << " model: " << y_model 
                     << " error: " << dataError 
                     << " pull: " << pull 
                     << " chi2: " << pull * pull << endl
                     << " ndf: " << ndf << endl;
                     */
            }
        }
    }
    ndf -= 2;
    cout << "Total chi2/ndf: " << chi2/ndf << endl;
    return chi2;
}

// 创建理论曲线函数
void DrawTheoryLines(TCanvas* c, double ymin, double ymax) {
    // 创建理论曲线
    const int nPoints = 1000;
    double x[nPoints], y7[nPoints], y9[nPoints], y10[nPoints];
    
    for(int i = 0; i < nPoints; i++) {
        x[i] = i * 30.0 / nPoints;
        y7[i] = rigidityToBeta(x[i], 4, 7, false);
        y9[i] = rigidityToBeta(x[i], 4, 9, false);
        y10[i] = rigidityToBeta(x[i], 4, 10, false);
        if(x[i] <= 10 && x[i] >=9.99)
        {
            //cout<<y7[i]<<endl;
            //cout<<y9[i]<<endl;
            //cout<<y10[i]<<endl;
        }
    }

    // 创建图形并设置样式
    TGraph* g7 = new TGraph(nPoints, x, y7);
    TGraph* g9 = new TGraph(nPoints, x, y9);
    TGraph* g10 = new TGraph(nPoints, x, y10);

    // 设置样式
    g7->SetLineColor(kBlue);
    g9->SetLineColor(kGreen+2);
    g10->SetLineColor(kOrange+1);

    for(auto g : {g7, g9, g10}) {
        g->SetLineWidth(3);
        g->SetLineStyle(2);
        g->Draw("L SAME");
    }

    // 添加图例
    TLegend* leg = new TLegend(0.63, 0.22, 0.85, 0.45);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(g7, "Be7 cutoff beta", "l");
    leg->AddEntry(g9, "Be9 cutoff beta", "l");
    leg->AddEntry(g10, "Be10 cutoff beta", "l");
    leg->Draw();
}

// Landau-Gauss convolution
Double_t langaufun(Double_t *x, Double_t *par) {
    Double_t invsq2pi = 0.3989422804014;
    Double_t mpshift = -0.22278298;
    Double_t np = 1500.0;
    Double_t sc = 5.0;
    Double_t xx, mpc, fland, sum = 0.0, xlow, xupp, step, i;
    mpc = par[1] - mpshift * par[0];
    xlow = x[0]*par[4] - sc*par[3];
    xupp = x[0]*par[4] + sc*par[3];
    step = (xupp - xlow) / np;
    for (i = 1.0; i <= np/2; i++) {
        xx = xlow + (i-.5)*step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0]*par[4], xx, par[3]);
        xx = xupp - (i-.5)*step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0]*par[4], xx, par[3]);
    }
    return (par[2] * step * sum * invsq2pi / par[3]);
}

// ExpGausExp
double funcExpGausExp(double *x, double *p) {
    constexpr double sqrtPiOver2 = 1.2533141373;
    constexpr double sqrt2 = 1.4142135624;
    const double &x0 = p[0], &sigmaL = p[1], &alphaL = p[2], &sigmaR = p[3], &alphaR = p[4], &norm = p[5], &xmin = p[6], &xmax = p[7];
    double t = (x[0] - x0) / (x[0] < x0 ? sigmaL : sigmaR);
    double v = 0;
    if (t < -alphaL) {
        double a = 0.5 * alphaL * alphaL, b = alphaL * t;
        v = std::exp(a + b);
    } else if (t <= alphaR) {
        v = std::exp(-0.5 * t * t);
    } else {
        double a = 0.5 * alphaR * alphaR, b = alphaR * (-t);
        v = std::exp(a + b);
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
        double sigmaMin = (tmin < double(0)) ? sigmaL : sigmaR;
        double sigmaMax = (tmax < double(0)) ? sigmaL : sigmaR;
        sum += sqrtPiOver2 * (sigmaMax * std::erf(std::min(tmax, alphaR) / sqrt2) - sigmaMin * std::erf(std::max(tmin, -alphaL) / sqrt2));
    }
    return norm * (v / sum);
}

int findStartBin() {
    for (size_t i = 0; i < Binning::NarrowBins.size(); ++i) {
        if (Binning::NarrowBins[i] >= 0.51) {
            return i;
        }
    }
    return 0;
}

double getCurrentRSS_MB() {
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    #ifdef __APPLE__
        return rusage.ru_maxrss / 1024. / 1024.;
    #else
        return rusage.ru_maxrss / 1024.;
    #endif
}

// Get 1D projection from 2D histogram for specific energy bin
std::unique_ptr<TH1D> getEnergySlice(int plusIbin, int Nrebin, TH2D* h2d, int binIndex, const std::string& name) {
    if (!h2d || binIndex < 0 || binIndex >= Binning::NarrowBins.size()-1) {
        return nullptr;
    }
    
    double eMin = Binning::NarrowBins[binIndex];
    double eMax = Binning::NarrowBins[binIndex + plusIbin];
    
    TAxis* yAxis = h2d->GetYaxis();
    int yBinLow = yAxis->FindBin(eMin + 1e-6);
    int yBinHigh = yAxis->FindBin(eMax - 1e-6);
    
    if (yBinLow > yBinHigh || yBinLow < 1 || yBinHigh > yAxis->GetNbins()) {
        std::cout << "Warning: Invalid energy bin range for " << name << std::endl;
        return nullptr;
    }
    
    std::unique_ptr<TH1D> slice(h2d->ProjectionX(name.c_str(), yBinLow, yBinHigh));
    if (slice) {
        slice->SetTitle(Form("E_{k}: %.2f-%.2f GeV/n", eMin, eMax));
        slice->SetDirectory(nullptr);
    }

    slice->Rebin(Nrebin);
    if(binIndex > 0) {
        //slice->Rebin(2);
    }
    else  {
        //slice->Rebin(2);
    }
    
    return slice;
}

// Extend histogram to specified range
std::unique_ptr<TH1D> extendHistogram(const TH1D* hist, double new_min, double new_max) {
    if (!hist) return nullptr;
    
    double curr_min = hist->GetXaxis()->GetXmin();
    double curr_max = hist->GetXaxis()->GetXmax();
    int curr_nbins = hist->GetNbinsX();
    double bin_width = (curr_max - curr_min) / curr_nbins;
    
    if (curr_min <= new_min && curr_max >= new_max) {
        return std::unique_ptr<TH1D>(static_cast<TH1D*>(hist->Clone((std::string(hist->GetName()) + "_clone").c_str())));
    }
    
    int new_nbins = static_cast<int>((new_max - new_min) / bin_width);
    
    std::unique_ptr<TH1D> new_hist(new TH1D(
        (std::string(hist->GetName()) + "_ext").c_str(),
        hist->GetTitle(),
        new_nbins, new_min, new_max
    ));
    
    // Copy style
    new_hist->SetLineColor(hist->GetLineColor());
    new_hist->SetLineWidth(hist->GetLineWidth());
    new_hist->SetLineStyle(hist->GetLineStyle());
    new_hist->SetFillColor(hist->GetFillColor());
    new_hist->SetFillStyle(hist->GetFillStyle());
    new_hist->SetMarkerColor(hist->GetMarkerColor());
    new_hist->SetMarkerStyle(hist->GetMarkerStyle());
    new_hist->SetMarkerSize(hist->GetMarkerSize());
    
    // Fill new histogram
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double x = hist->GetBinCenter(i);
        if (hist->GetBinContent(i) > 0) {
            int new_bin = new_hist->FindBin(x);
            new_hist->SetBinContent(new_bin, hist->GetBinContent(i));
            new_hist->SetBinError(new_bin, hist->GetBinError(i));
        }
        else{
            new_hist->SetBinContent(i, 0);
            new_hist->SetBinError(i, 0);
        }
    }
    
    new_hist->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
    new_hist->GetXaxis()->SetRangeUser(new_min, new_max);
    new_hist->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
    new_hist->SetDirectory(nullptr);
    
    return new_hist;
}

// Calculate Pull distribution
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
                double pull_err = fabs(pull_val * (sqrt(y_data) / y_data));
                if (std::isfinite(pull_val) && abs(pull_val) < 10) {
                    pullGraph->SetPoint(pointIndex, x_val, pull_val);
                    pullGraph->SetPointError(pointIndex, 0, 0);
                    pointIndex++;
                }
            }
        }
    }
}

// Set up Pull plot style
void setupPullPlot(TGraphErrors* pullGraph, double fit_min, double fit_max) {
    if (!pullGraph) return;
    
    pullGraph->SetTitle(";L1 Charge;Pull");
    pullGraph->SetMarkerColor(kBlack);
    pullGraph->SetMarkerStyle(20);
    pullGraph->SetMarkerSize(.8);
    pullGraph->SetLineColor(kBlack);
    
    pullGraph->GetXaxis()->SetLabelSize(0.12);
    pullGraph->GetYaxis()->SetLabelSize(0.10);
    pullGraph->GetXaxis()->SetTitleSize(0.14);
    pullGraph->GetYaxis()->SetTitleSize(0.14);
    pullGraph->GetYaxis()->SetTitleOffset(0.4);
    pullGraph->GetXaxis()->SetTitleOffset(.9);  // Set x title offset to 1.0
    pullGraph->GetXaxis()->SetRangeUser(fit_min, fit_max);
    pullGraph->GetYaxis()->SetRangeUser(-6, 6);  // Change y range to -6 to 6
    pullGraph->GetYaxis()->SetNdivisions(505);
}

// Create legend
std::unique_ptr<TLegend> createLegend() {
    auto legend = std::make_unique<TLegend>(0.72, 0.72, 0.88, 0.86);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);
    return legend;
}


} // namespace AMS_Iso