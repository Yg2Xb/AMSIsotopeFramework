#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <cstdarg>
#include <iomanip>

#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TLatex.h>

#include "../Tool.h" // must contain definitions of langaufun and funcExpGausExp

using namespace AMS_Iso;
using std::string;
using std::vector;
using std::map;
using std::unique_ptr;
using std::cout;
using std::endl;

// Config (keep simple)
struct DebugCfg {
    bool infoOpen = true;
    bool fitSummary = true;
    bool rangeLog = false;
    bool boundaryWarn = true;
} DBG;

static inline void printOpen(const string& tag, const string& path, TFile* f) {
    if (DBG.infoOpen) cout << "[open] " << tag << " '" << path << "' -> " << ((f && !f->IsZombie()) ? "OK" : "FAIL") << endl;
}

static inline void drawFitRangeLines(double fitLow, double fitHigh, int color, TH1D* hist) {
    if (!gPad || !hist) return;
    gPad->Modified(); gPad->Update();

    double yMin, yMax;

    if (gPad->GetLogy()) {
        yMin = TMath::Power(10, gPad->GetUymin());
        yMax = TMath::Power(10, gPad->GetUymax());
    } else {
        yMin = gPad->GetUymin();
        yMax = gPad->GetUymax();
    }

    TLine* L1 = new TLine(fitLow,  yMin, fitLow,  yMax);
    TLine* L2 = new TLine(fitHigh, yMin, fitHigh, yMax);
    L1->SetLineColor(color); L2->SetLineColor(color);
    L1->SetLineStyle(3);      L2->SetLineStyle(3);
    L1->SetLineWidth(2);      L2->SetLineWidth(2);
    L1->Draw("same");         L2->Draw("same");
}

static double findContentLevel(double charge, TH1D* hist, double ratio, bool searchLeft) {
    if (!hist) return 0;
    int binMin = hist->GetXaxis()->FindBin(charge - 0.12);
    int binMax = hist->GetXaxis()->FindBin(charge + 0.12);
    if (binMin < 1) binMin = 1;
    if (binMax > hist->GetNbinsX()) binMax = hist->GetNbinsX();

    int maxBin = binMin;
    double maxContent = hist->GetBinContent(binMin);
    for (int bin = binMin + 1; bin <= binMax; ++bin) {
        double c = hist->GetBinContent(bin);
        if (c > maxContent) { maxContent = c; maxBin = bin; }
    }
    double target = maxContent * ratio;
    double targetX = hist->GetXaxis()->GetBinCenter(maxBin);

    if (searchLeft) {
        for (int bin = maxBin; bin >= 1; --bin) {
            if (hist->GetBinContent(bin) <= target) { targetX = hist->GetXaxis()->GetBinCenter(bin); break; }
        }
    } else {
        for (int bin = maxBin; bin <= hist->GetNbinsX(); ++bin) {
            if (hist->GetBinContent(bin) <= target) { targetX = hist->GetXaxis()->GetBinCenter(bin); break; }
        }
    }
    return targetX;
}

static void findFitRange(TH1D* hist, double chargeValue, double leftRatio, double rightRatio, double& lowEdge, double& highEdge) {
    if (!hist) { lowEdge = highEdge = 0; return; }
    lowEdge = findContentLevel(chargeValue, hist, leftRatio, true);
    highEdge = findContentLevel(chargeValue, hist, rightRatio, false);
    if (lowEdge > chargeValue - 0.24) lowEdge = chargeValue - 0.24;
    if (lowEdge < chargeValue - 0.6) lowEdge = chargeValue - 0.6;
    if (highEdge < chargeValue + 0.35) highEdge = chargeValue + 0.35;
    if (highEdge > chargeValue + 0.7) highEdge = chargeValue + 0.7;
    if (DBG.rangeLog) {
        cout << std::fixed << std::setprecision(4)
             << " [range] Z=" << chargeValue << " -> [" << lowEdge << ", " << highEdge << "]" << endl;
    }
}

static inline bool passEnergyWindow(const string& det, double ekCen) {
    if (det == "TOF") return ekCen > 0.4 && ekCen <= 1.28;
    if (det == "NaF") return ekCen > 0.71 && ekCen <= 5.1;
    if (det == "AGL") return ekCen > 2.8 && ekCen <= 20.00;
    return false;
}

// Parameter initializer
class FitParameterManager {
public:
    FitParameterManager(bool firstFit, const string& splineFilePath, const string& histOriFilePath)
        : firstFit_(firstFit), splineFile_(nullptr), histFile_(nullptr)
    {
        if (!firstFit_) {
            splineFile_ = TFile::Open(splineFilePath.c_str());
            printOpen("spline", splineFilePath, splineFile_);
            histFile_ = TFile::Open(histOriFilePath.c_str());
            printOpen("hist_ori", histOriFilePath, histFile_);
        }
    }

    ~FitParameterManager() {
        if (splineFile_) { splineFile_->Close(); delete splineFile_; }
        if (histFile_) { histFile_->Close(); delete histFile_; }
    }

    void setParameters(TF1& f, const string& fitType,
                       const string& chain, const string& elem, const string& det, const string& type,
                       double ekCen, double charge, double histMax,
                       double fitLow, double fitHigh)
    {
        if (firstFit_) {
            setParameters_firstFit(f, fitType, charge, histMax, fitLow, fitHigh);
        } else {
            setParameters_fromFile(f, fitType, chain, elem, det, type, ekCen, charge, histMax, fitLow, fitHigh);
        }
    }

    double getInitialValue(const string& chain, const string& elem, const string& det, const string& type, 
                           const string& fitType, const string& parName, double ekCen) {
        if (firstFit_ || !splineFile_) return 0.0;
        string tf1Name = buildTF1Name(chain, elem, det, type, fitType, parName);
        TF1* f_param = dynamic_cast<TF1*>(splineFile_->Get(tf1Name.c_str()));
        if (f_param) {
            return f_param->Eval(ekCen);
        }
        return 0.0;
    }

    vector<int> checkBoundaries(const TF1& func, const string& fitName,
                                const string& chain, const string& elem, const string& det, const string& type,
                                double ekCen) const
    {
        vector<int> boundaryParams;
        if (!DBG.boundaryWarn) return boundaryParams;

        for (int i = 0; i < func.GetNpar(); ++i) {
            double v = func.GetParameter(i), lo, hi;
            func.GetParLimits(i, lo, hi);
            if (lo < hi) {
                const double eps = 1e-8;
                if (TMath::Abs(v - lo) < eps || TMath::Abs(v - hi) < eps) {
                    cout << "  [warn] boundary: " << fitName
                         << " par#" << i << " " << func.GetParName(i)
                         << " = " << v << " in [" << lo << "," << hi << "] "
                         << chain << "/" << elem << "/" << det << "/" << type << " EkCen=" << ekCen << endl;
                    boundaryParams.push_back(i);
                }
            }
        }
        return boundaryParams;
    }

private:
    bool firstFit_;
    TFile* splineFile_;
    TFile* histFile_;
    map<string, double> alphaLAverageCache_;

    string buildTF1Name(const string& chain, const string& elem, const string& det, const string& type, const string& fitType, const string& parName) const {
        return chain + "_" + elem + "_" + det + "_" + type + "_" + fitType + "_" + parName + "_spline";
    }

    string buildHistName(const string& chain, const string& elem, const string& det, const string& type, const string& fitType, const string& parName) const {
        return chain + "_" + elem + "_" + det + "_" + type + "_" + fitType + "_" + parName;
    }

    void setParameters_firstFit(TF1& f, const string& fitType, double charge, double histMax, double fitLow, double fitHigh) {
        if (fitType == "LG") {
            f.SetParameters(0.05, charge, histMax, 0.20, 1.0);
            f.FixParameter(4, 1.0);
            f.SetParLimits(0, 0.005, 0.11);
            f.SetParLimits(1, charge - 0.15, charge + 0.18);
            f.SetParLimits(3, 0.06, 0.5);
        } else { // EGE
            f.SetParameters(charge, 0.20, 2.0, 0.20, 2.0, histMax, fitLow, fitHigh);
            f.SetParLimits(1, 0.01, 0.5);
            f.SetParLimits(2, 0.5, 4.0);
            f.SetParLimits(3, 0.01, 0.5);
            f.SetParLimits(4, 0.5, 4.0);
            f.FixParameter(6, fitLow);
            f.FixParameter(7, fitHigh);
        }
    }

    void setParameters_fromFile(TF1& f, const string& fitType,
                                const string& chain, const string& elem, const string& det, const string& type,
                                double ekCen, double charge, double histMax,
                                double fitLow, double fitHigh)
    {
        // Set base parameters first, in case some splines are missing
        setParameters_firstFit(f, fitType, charge, histMax, fitLow, fitHigh);

        if (!splineFile_ || !histFile_) {
            cout << "  [warn] Spline/Hist file not available. Using first fit parameters." << endl;
            return;
        }

        for (int i = 0; i < f.GetNpar(); ++i) {
            string parName = f.GetParName(i);
            if ((fitType == "LG" && (parName == "Area" || i == 4 /* p4 */)) ||
                (fitType == "EGE" && (parName == "Norm" || parName == "xmin" || parName == "xmax"))) {
                continue;
            }
            
            if (fitType == "EGE" && parName == "AlphaL") {
                double avgAlphaL = getAlphaLAverage(chain, elem, det, type);
                if (avgAlphaL > 0) {
                    f.FixParameter(i, avgAlphaL);
                } else {
                    cout << "  [warn] Could not get AlphaL average for " << chain << "/" << elem << "/" << det << "/" << type << ". Using default fixed value 2.0" << endl;
                    f.FixParameter(i, 2.0);
                }
                continue;
            }

            string tf1Name = buildTF1Name(chain, elem, det, type, fitType, parName);
            TF1* f_param = dynamic_cast<TF1*>(splineFile_->Get(tf1Name.c_str()));

            if (f_param) {
                double initialVal = f_param->Eval(ekCen);
                //cout<<" [param] " << tf1Name << " -> " << initialVal << endl;
                f.SetParameter(i, initialVal);
                
                double lowerLimit = 0.75 * initialVal;
                double upperLimit = 1.25 * initialVal;

                if (parName == "AlphaR" || parName == "Width") {
                    lowerLimit = 0.5 * initialVal;
                    upperLimit = 1.5 * initialVal;
                }
                
                if (fitType == "LG" && parName == "Width") {
                    if (upperLimit > 0.110) upperLimit = 0.110;
                }
                
                // Ensure limits are sane
                if (lowerLimit < upperLimit) {
                    f.SetParLimits(i, lowerLimit, upperLimit);
                } else {
                    cout<<" !!!!!!!"<<ekCen<<endl;
                    cout << "  [warn] Insane limits for " << parName << ": [" << lowerLimit << ", " << upperLimit << "]. Using defaults." << endl;
                }

            } else {
                cout << "  [warn] Could not find TF1: " << tf1Name << ". Using default limits." << endl;
                // Defaults are already set by setParameters_firstFit, so nothing to do here.
            }
        }
    }

    double getAlphaLAverage(const string& chain, const string& elem, const string& det, const string& type) {
        string cacheKey = chain + "_" + elem + "_" + det + "_" + type;
        if (alphaLAverageCache_.count(cacheKey)) return alphaLAverageCache_[cacheKey];
        if (!histFile_) return 0.0;

        string histName = buildHistName(chain, elem, det, type, "EGE", "AlphaL");
        if (histName.empty()) return 0.0;

        TH1D* h = dynamic_cast<TH1D*>(histFile_->Get(histName.c_str()));
        if (!h) {
            cout << "  [warn] Could not find AlphaL histogram: " << histName << " in hist_ori file." << endl;
            return 0.0;
        }

        double rangeLow = 0, rangeHigh = 0;
        if (det == "TOF") { rangeLow = 0.42; rangeHigh = 1.55; }
        else if (det == "NaF") { rangeLow = 0.86; rangeHigh = 4.91; }
        else if (det == "AGL") { rangeLow = 2.88; rangeHigh = 16.3; }
        else return 0.0;
        
        int binLow = h->GetXaxis()->FindBin(rangeLow + 0.01);
        int binHigh = h->GetXaxis()->FindBin(rangeHigh - 0.01);
        double sum = 0; int count = 0;
        for (int bin = binLow; bin <= binHigh; ++bin) {
            double content = h->GetBinContent(bin);
            if (content != 0) { sum += content; count++; }
        }
        double average = (count > 0) ? sum / count : 0.0;
        alphaLAverageCache_[cacheKey] = average;
        return average;
    }
};

// Build TF1 (no lambdas)
static TF1* BuildLG(double lo, double hi) { return new TF1("fLG", langaufun, lo, hi, 5); }
static TF1* BuildEGE(double lo, double hi){ return new TF1("fEGE", funcExpGausExp, lo, hi, 8); }

// One-pass fitter with enhanced retry logic
static TF1* doFit(TH1D* h, TF1* (*builder)(double,double),
                  const vector<std::pair<string,int>>& parList,
                  FitParameterManager& pm, const string& fitType,
                  const string& chain, const string& elem, const string& det, const string& type,
                  double ekCen, double chargeZ, double histMax,
                  double& fitLow, double& fitHigh, int maxIter = 4, const char* fitOpt = "RQ0")
{
    TF1* f = nullptr;
    bool fitSucceeded = false; // NEW: Flag to track success

    for (int iter = 0; iter < maxIter; ++iter) {
        if (f) { delete f; }
        f = builder(fitLow, fitHigh);
        for (size_t i = 0; i < parList.size(); ++i)
            f->SetParName(parList[i].second, parList[i].first.c_str());

        pm.setParameters(*f, fitType, chain, elem, det, type, ekCen, chargeZ, histMax, fitLow, fitHigh);
        
        // --- First Fit Attempt ---
        h->Fit(f, fitOpt, "", fitLow, fitHigh);
        vector<int> boundary_params = pm.checkBoundaries(*f, fitType + "-try1", chain, elem, det, type, ekCen);

        // --- If at boundary, adjust parameter limits and refit ---
        if (!boundary_params.empty()) {
            cout << "  [info] Params at boundary. Adjusting limits and refitting..." << endl;
            for (int parIdx : boundary_params) {
                string parName = f->GetParName(parIdx);
                double initialVal = pm.getInitialValue(chain, elem, det, type, fitType, parName, ekCen);
                if (initialVal > 0) {
                    double min_lim, max_lim;
                    f->GetParLimits(parIdx, min_lim, max_lim);
                    double adjustment = 0.2 * initialVal;
                    double new_min = min_lim - adjustment;
                    double new_max = max_lim + adjustment;

                    if (fitType == "LG" && parName == "Width" && new_max > 0.110) {
                        new_max = 0.110;
                    }
                    if (new_min < new_max) {
                        cout << "    -> Adjusting " << parName << " limits to [" << new_min << ", " << new_max << "]" << endl;
                        f->SetParLimits(parIdx, new_min, new_max);
                    }
                }
            }
            // Refit with adjusted parameter limits
            h->Fit(f, fitOpt, "", fitLow, fitHigh);
            boundary_params = pm.checkBoundaries(*f, fitType + "-try2", chain, elem, det, type, ekCen);
        }

        // --- Check Chi2, if good, break. Otherwise, adjust fit range. ---
        double ndf = f->GetNDF();
        double chi2ndf = (ndf > 0) ? f->GetChisquare()/ndf : 1e9;
        if (chi2ndf < 3.0 && boundary_params.empty()) {
            fitSucceeded = true; // NEW: Mark as succeeded
            break;
        }

        // If still failing, adjust fit range for next iteration
        if (iter < maxIter - 1) {
            cout << "  [info] Fit not optimal (Chi2/NDF=" << chi2ndf << ", Boundary Hits=" << boundary_params.size() 
                 << "). Adjusting fit range for iter " << iter + 2 << "..." << endl;
            if (fitType == "LG") {
                if (iter == 0) { findFitRange(h, chargeZ, 0.4, 0.2, fitLow, fitHigh); }
                else { fitLow += 0.05; fitHigh -= 0.05; }
            } else { // EGE
                fitLow += 0.1; fitHigh -= 0.1;
            }
            if(fitLow >= chargeZ - 0.28) fitLow = chargeZ - 0.28;
        }
    }

    // NEW: Check the flag after the loop and print a final failure message if needed
    if (!fitSucceeded) {
        double final_chi2ndf = (f && f->GetNDF() > 0) ? f->GetChisquare() / f->GetNDF() : -1.0;
        if(final_chi2ndf > 4.){
            cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            cout << "  [FIT FAILED] After all iterations for " << fitType << " on "
                 << chain << "/" << elem << "/" << det << "/" << type << " at EkCen=" << ekCen << endl;
            cout << "  Final Chi2/NDF = " << final_chi2ndf << endl;
            cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        }
    }

    return f;
}

void chargeHistFit(
    const string& histFile = "/eos/user/z/zixuan/Isotope/Add/Be_frag4.root",
    const string& pdfOut = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFits_BeToOxy_0.5_iter1.pdf",
    const string& histOut = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFitParams_BeToOxy_0.5_iter1.root",
    int rebin = 2,
    bool firstFit = false,
    bool listKeysOnce = false // kept for compatibility, unused to stay simple
) {
    const vector<string> chains = {"L1Inner", "UnbiasedL1Inner"};
    const vector<string> nuclei = {"Beryllium","Boron","Carbon","Nitrogen","Oxygen"};
    //const vector<string> nuclei = {"Nitrogen"};
    const map<string,double> chargeZ = {{"Beryllium",4.0}, {"Boron",5.0}, {"Carbon",6.0}, {"Nitrogen",7.0}, {"Oxygen",8.0}};
    const vector<string> types = {"L1QTemplate", "L2QTemplate"};
    const vector<string> dets = {"TOF","NaF","AGL"};
    const map<string,int> detColor = {{"TOF", kRed}, {"NaF", kBlue}, {"AGL", kGreen+2}};
    const map<string,int> typeMarker = {{"L1QTemplate", 20}, {"L2QTemplate", 21}};

    const vector<std::pair<string,int>> LG_params = {{"Width",0},{"MPV",1},{"Area",2},{"Sigma",3}};
    const vector<std::pair<string,int>> EGE_params = {{"Peak",0},{"SigmaL",1},{"AlphaL",2},{"SigmaR",3},{"AlphaR",4},{"Norm",5},{"xmin",6},{"xmax",7}};

    // --- FILE PATHS FOR SMART FITTING ---
    const string splineFilePath = "/eos/user/z/zixuan/Isotope/ChargeFit/comparison_plots/allFitHistSplineSmooth_0.8_orig.root";
    const string histOriFilePath = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFitParams_BeToOxy_0.8_orig.root";

    FitParameterManager pm(firstFit, splineFilePath, histOriFilePath);

    unique_ptr<TFile> fin(TFile::Open(histFile.c_str()));
    printOpen("input", histFile, fin.get());
    if (!fin || fin->IsZombie()) return;

    unique_ptr<TFile> fout(TFile::Open(histOut.c_str(), "RECREATE"));
    printOpen("output", histOut, fout.get());
    if (!fout || fout->IsZombie()) return;

    unique_ptr<TCanvas> c(new TCanvas("c","",900,700));
    c->SetLogy(true);
    c->Print((pdfOut + "[").c_str(), "pdf");

    struct ParSet {
        bool init = false;
        int nb = 0;
        std::vector<double> edges;
        map<string, TH1D*> h;
    };
    map<string, ParSet> store;

    auto keyEDTChain = [](const string& chain, const string& e, const string& d, const string& t){ 
        return chain + "_" + e + "_" + d + "_" + t; 
    };
    auto namePar = [](const string& chain, const string& e, const string& d, const string& t, const string& fit, const string& par){
        return chain + "_" + e + "_" + d + "_" + t + "_" + fit + "_" + par;
    };

    for (size_t ich = 0; ich < chains.size(); ++ich)
    for (size_t ie = 0; ie < nuclei.size(); ++ie)
    for (size_t id = 0; id < dets.size(); ++id)
    for (size_t it = 0; it < types.size(); ++it)
    {
        const string& chain = chains[ich];
        const string& elem = nuclei[ie];
        const string& det  = dets[id];
        const string& type = types[it];

        string key = chain + "_ISS_BKG_H2_" + elem + "_" + type + "_" + det;
        TH2* h2 = dynamic_cast<TH2*>(fin->Get(key.c_str()));
        if (!h2) { if (DBG.infoOpen) cout << "[miss] " << key << endl; continue; }

        string edtc = keyEDTChain(chain, elem, det, type);
        ParSet& PS = store[edtc];
        if (!PS.init) {
            int nY = h2->GetYaxis()->GetNbins();
            PS.nb = nY;
            PS.edges.resize(nY+1);
            for (int i=1;i<=nY;++i) PS.edges[i-1] = h2->GetYaxis()->GetBinLowEdge(i);
            PS.edges[nY] = h2->GetYaxis()->GetBinUpEdge(nY);
            
            auto mk = [&](const string& fit, const string& par){
                string nm = namePar(chain, elem, det, type, fit, par);
                string tt = nm + ";E_{k}/n [GeV/n];" + par;
                TH1D* h = new TH1D(nm.c_str(), tt.c_str(), PS.nb, PS.edges.data());
                h->SetDirectory(nullptr);
                h->SetLineColor(detColor.at(det));
                h->SetMarkerColor(detColor.at(det));
                h->SetMarkerStyle(typeMarker.at(type));
                h->SetLineWidth(2);
                PS.h[nm] = h;
            };

            for (size_t iP=0;iP<LG_params.size();++iP) mk("LG", LG_params[iP].first);
            mk("LG","Chi2NDF");
            for (size_t iP=0;iP<EGE_params.size();++iP) mk("EGE", EGE_params[iP].first);
            mk("EGE","Chi2NDF");
            PS.init = true;
            if (DBG.infoOpen) cout << "[init] " << edtc << " with " << PS.nb << " Ek/n bins" << endl;
        }

        int nY = h2->GetYaxis()->GetNbins();
        for (int ybin=1; ybin<=nY; ++ybin) {
            double ekLow  = h2->GetYaxis()->GetBinLowEdge(ybin);
            double ekHigh = h2->GetYaxis()->GetBinUpEdge(ybin);
            double ekCen = 0.5*(ekLow+ekHigh);
            if (!passEnergyWindow(det, ekCen)) continue;

            string projName = key + Form("_projY%d", ybin);
            TH1D* h1 = dynamic_cast<TH1D*>(h2->ProjectionX(projName.c_str(), ybin, ybin));
            if (!h1) continue;
            if (h1->GetMaximum() < 16) { delete h1; continue; }
            if (rebin > 1) h1->Rebin(rebin);
            h1->Sumw2();

            double z = chargeZ.at(elem);

            // LG
            double fitLowLG=0, fitHighLG=0;
            findFitRange(h1, z, 0.10, 0.10, fitLowLG, fitHighLG);
            TF1* fLG = doFit(h1, BuildLG, LG_params, pm, "LG", chain, elem, det, type, ekCen, z, h1->GetMaximum(), fitLowLG, fitHighLG);

            // EGE
            double fitLowEGE=0, fitHighEGE=0;
            findFitRange(h1, z, 0.10, 0.10, fitLowEGE, fitHighEGE);
            TF1* fEGE = doFit(h1, BuildEGE, EGE_params, pm, "EGE", chain, elem, det, type, ekCen, z, h1->GetMaximum(), fitLowEGE, fitHighEGE);

            int b = PS.h[namePar(chain, elem, det, type, "LG", "Width")]->FindBin(ekCen);
            for (size_t iP=0; iP<LG_params.size(); ++iP) {
                const string& par = LG_params[iP].first; int ip = LG_params[iP].second;
                TH1D* hh = PS.h[namePar(chain, elem, det, type, "LG", par)];
                hh->SetBinContent(b, fLG->GetParameter(ip));
                hh->SetBinError(b, fLG->GetParError(ip));
            }
            {
                TH1D* hchi = PS.h[namePar(chain, elem, det, type, "LG", "Chi2NDF")];
                hchi->SetBinContent(b, (fLG->GetNDF()>0)? fLG->GetChisquare()/fLG->GetNDF() : 0);
                hchi->SetBinError(b, 0);
            }

            for (size_t iP=0; iP<EGE_params.size(); ++iP) {
                const string& par = EGE_params[iP].first; int ip = EGE_params[iP].second;
                TH1D* hh = PS.h[namePar(chain, elem, det, type, "EGE", par)];
                hh->SetBinContent(b, fEGE->GetParameter(ip));
                hh->SetBinError(b, fEGE->GetParError(ip));
            }
            {
                TH1D* hchi = PS.h[namePar(chain, elem, det, type, "EGE", "Chi2NDF")];
                hchi->SetBinContent(b, (fEGE->GetNDF()>0)? fEGE->GetChisquare()/fEGE->GetNDF() : 0);
                hchi->SetBinError(b, 0);
            }

            h1->SetStats(0);
            h1->SetTitle(Form("%s | %s | %s | %s | Ek/n: [%.3f, %.3f] GeV/n", chain.c_str(), type.c_str(), elem.c_str(), det.c_str(), ekLow, ekHigh));
            h1->GetXaxis()->SetTitle("Tracker Layer Q");
            h1->GetYaxis()->SetTitle("Counts");
            h1->GetYaxis()->SetTitleOffset(1.4);
            h1->SetMarkerStyle(20); h1->SetMarkerSize(0.9); h1->SetMarkerColor(kBlack); h1->SetLineColor(kBlack);
            h1->GetXaxis()->SetRangeUser(z - 2.0, z + 1.5);

            c->cd(); c->SetLogy(true);
            h1->Draw("E");

            fLG->SetRange(h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
            fLG->SetLineColor(kRed); fLG->SetLineWidth(3); fLG->SetLineStyle(1); fLG->Draw("same");
            fEGE->SetRange(h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
            fEGE->SetLineColor(kGreen+2); fEGE->SetLineWidth(3); fEGE->SetLineStyle(1); fEGE->Draw("same");

            drawFitRangeLines(fitLowLG,  fitHighLG,  kRed, h1);
            drawFitRangeLines(fitLowEGE, fitHighEGE, kGreen+2, h1);

            TLegend leg(0.66, 0.66, 0.93, 0.88);
            leg.SetTextSize(0.03); leg.SetBorderSize(0); leg.SetFillStyle(0);
            leg.AddEntry(h1, "Data", "ep");
            leg.AddEntry(fLG, "Landau-Gauss", "l");
            leg.AddEntry(fEGE, "ExpGausExp", "l");
            leg.Draw();

            TLatex lt; 
            lt.SetNDC(); 
            lt.SetTextFont(62); 
            lt.SetTextSize(0.028);
            lt.SetTextAlign(12);

            double x = 0.16;
            double y = 0.86;
            double dy = 0.032;

            lt.SetTextColor(kRed);
            lt.DrawLatex(x, y, "LandauGauss:"); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("MPV = %.3f #pm %.3f", fLG->GetParameter(1), fLG->GetParError(1))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("Width = %.3f #pm %.3f", fLG->GetParameter(0), fLG->GetParError(0))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("Sigma = %.3f #pm %.3f", fLG->GetParameter(3), fLG->GetParError(3))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("#chi^{2}/ndf = %.0f/%d = %.2f", 
                fLG->GetChisquare(), (int)fLG->GetNDF(), 
                (fLG->GetNDF()>0? fLG->GetChisquare()/fLG->GetNDF():0))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("fit range: %.2f-%.2f", fitLowLG, fitHighLG)); y -= dy*1.2;

            lt.SetTextColor(kGreen+2);
            lt.DrawLatex(x, y, "ExpGausExp:"); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("Peak = %.4f #pm %.4f", fEGE->GetParameter(0), fEGE->GetParError(0))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("#sigma_{L} = %.4f #pm %.4f", fEGE->GetParameter(1), fEGE->GetParError(1))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("#sigma_{R} = %.4f #pm %.4f", fEGE->GetParameter(3), fEGE->GetParError(3))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("#alpha_{L} = %.4f #pm %.4f", fEGE->GetParameter(2), fEGE->GetParError(2))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("#alpha_{R} = %.4f #pm %.4f", fEGE->GetParameter(4), fEGE->GetParError(4))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("#chi^{2}/ndf = %.0f/%d = %.2f", 
                fEGE->GetChisquare(), (int)fEGE->GetNDF(), 
                (fEGE->GetNDF()>0? fEGE->GetChisquare()/fEGE->GetNDF():0))); y -= dy;
            lt.DrawLatex(x + 0.02, y, Form("fit range: %.2f-%.2f", fitLowEGE, fitHighEGE));

            c->Print(pdfOut.c_str(), "pdf");

            if (DBG.fitSummary) {
                cout << std::fixed << std::setprecision(3)
                     << "[ok] " << chain << " " << elem << " " << type << " " << det
                     << " EkCen=" << ekCen
                     << " LG chi2/ndf=" << (fLG->GetNDF()>0? fLG->GetChisquare()/fLG->GetNDF():0)
                     << " EGE chi2/ndf=" << (fEGE->GetNDF()>0? fEGE->GetChisquare()/fEGE->GetNDF():0)
                     << endl;
            }

            delete fLG;
            delete fEGE;
            delete h1;
        }
    }

    c->Print((pdfOut + "]").c_str(), "pdf");

    fout->cd();
    for (auto& kv : store) {
        for (auto& hk : kv.second.h) hk.second->Write();
    }
    fout->Close();
    cout << "All parameter histograms saved to " << histOut << endl;
}