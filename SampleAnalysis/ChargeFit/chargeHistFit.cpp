#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <iomanip>
#include <cmath>
#include <algorithm>

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

#include "../Tool.h" 

using namespace AMS_Iso;
using std::string;
using std::vector;
using std::map;
using std::unique_ptr;
using std::cout;
using std::endl;

// --- Configuration ---
struct DebugCfg {
    bool infoOpen = true;
    bool fitSummary = true;
    bool boundaryWarn = true;
} DBG;

// --- Helpers ---
static inline void printOpen(const string& tag, const string& path, TFile* f) {
    if (DBG.infoOpen) cout << "[open] " << tag << " '" << path << "' -> " << ((f && !f->IsZombie()) ? "OK" : "FAIL") << endl;
}

static inline void drawFitRangeLines(double fitLow, double fitHigh, int color) {
    if (!gPad) return;
    gPad->Update();
    double yMin = gPad->GetLogy() ? TMath::Power(10, gPad->GetUymin()) : gPad->GetUymin();
    double yMax = gPad->GetLogy() ? TMath::Power(10, gPad->GetUymax()) : gPad->GetUymax();

    TLine* L1 = new TLine(fitLow, yMin, fitLow, yMax);
    TLine* L2 = new TLine(fitHigh, yMin, fitHigh, yMax);
    L1->SetLineColor(color); L2->SetLineColor(color);
    L1->SetLineStyle(7); L2->SetLineStyle(7);  
    L1->SetLineWidth(2); L2->SetLineWidth(2);
    L1->Draw("same"); L2->Draw("same");
}

static double findContentLevel_Local(double charge, TH1D* hist, double ratio, bool searchLeft) {
    if (!hist) return charge;
    int maxBin = hist->GetXaxis()->FindBin(charge);
    double target = hist->GetBinContent(maxBin) * ratio;
    
    if (searchLeft) {
        for (int bin = maxBin; bin >= 1; --bin) 
            if (hist->GetBinContent(bin) < target) return hist->GetBinCenter(bin);
        return hist->GetBinLowEdge(1);
    } else {
        for (int bin = maxBin; bin <= hist->GetNbinsX(); ++bin) 
            if (hist->GetBinContent(bin) < target) return hist->GetBinCenter(bin);
        return hist->GetBinLowEdge(hist->GetNbinsX()+1);
    }
}

static void findFitRange(TH1D* hist, double Z, double& lowEdge, double& highEdge) {
    if (!hist) { lowEdge = highEdge = 0; return; }
    
    // He needs wider tail scan (Fat tail)
    bool isHe = (Z < 2.5);
    double rLeft  = isHe ? 0.20 : 0.10;
    double rRight = isHe ? 0.20 : 0.10; 
    
    lowEdge  = findContentLevel_Local(Z, hist, rLeft, true);
    highEdge = findContentLevel_Local(Z, hist, rRight, false);

    double minW_L = isHe ? 0.05 : 0.24;
    double minW_R = isHe ? 0.05 : 0.35;
    double maxD_L = isHe ? 0.80 : 0.60;
    double maxD_R = isHe ? 2.00 : 0.70; 

    if (lowEdge > Z - minW_L) lowEdge = Z - minW_L;
    if (lowEdge < Z - maxD_L) lowEdge = Z - maxD_L;
    if (highEdge < Z + minW_R) highEdge = Z + minW_R;
    if (highEdge > Z + maxD_R) highEdge = Z + maxD_R;
}

static inline bool passEnergyWindow(const string& det, double ekCen) {
    if (det == "TOF") return ekCen > 0.27 && ekCen <= 1.5;
    if (det == "NaF") return ekCen > 0.71 && ekCen <= 6.1;
    if (det == "AGL") return ekCen > 2.7 && ekCen <= 22.0;
    return false;
}

// --- Parameter Manager ---
class FitParameterManager {
public:
    FitParameterManager(bool firstFit, const string& splinePath, const string& histPath)
        : firstFit_(firstFit) {
        if (!firstFit_) {
            splineFile_ = std::unique_ptr<TFile>(TFile::Open(splinePath.c_str()));
            printOpen("spline", splinePath, splineFile_.get());
            histFile_   = std::unique_ptr<TFile>(TFile::Open(histPath.c_str()));
            printOpen("hist_ori", histPath, histFile_.get());
        }
    }

    void setParameters(TF1& f, const string& keyBase, double ekCen, double Z, double maxVal, double fitLo, double fitHi) {
        if (firstFit_) setParameters_Init(f, Z, maxVal, fitLo, fitHi);
        else setParameters_File(f, keyBase, ekCen, Z, maxVal, fitLo, fitHi);
    }

    double getInitialValue(const string& keyBase, const string& parName, double ekCen) {
        if (firstFit_ || !splineFile_) return 0.0;
        string name = keyBase + "_EGE_" + parName + "_spline";
        if (auto fn = dynamic_cast<TF1*>(splineFile_->Get(name.c_str()))) return fn->Eval(ekCen);
        return 0.0;
    }

    vector<int> checkBoundaries(const TF1& f, const string& tag, double ekCen) const {
        vector<int> badParams;
        if (!DBG.boundaryWarn) return badParams;
        for (int i = 0; i < f.GetNpar(); ++i) {
            double v = f.GetParameter(i), lo, hi;
            f.GetParLimits(i, lo, hi);
            if (lo < hi && (std::abs(v - lo) < 1e-6 || std::abs(v - hi) < 1e-6)) {
                //cout << "  [warn] Bound: " << tag << " p" << i << "(" << f.GetParName(i) << ")=" 
                //     << v << " [" << lo << "," << hi << "] Ek=" << ekCen << endl;
                badParams.push_back(i);
            }
        }
        return badParams;
    }

private:
    bool firstFit_;
    std::unique_ptr<TFile> splineFile_, histFile_;

    void setParameters_Init(TF1& f, double Z, double maxVal, double lo, double hi) {
        bool isHe = (Z < 2.5);
        f.SetParameters(Z, isHe?0.25:0.20, 2.0, isHe?0.30:0.20, isHe?1.5:2.0, maxVal, lo, hi);
        
        double peakWin = isHe ? 0.6 : 0.2;
        f.SetParLimits(0, Z - peakWin, Z + peakWin); 
        f.SetParLimits(1, 0.06, isHe ? 1.0 : 0.6);   
        f.SetParLimits(2, 0.1,  8.0);                
        f.SetParLimits(3, 0.06, isHe ? 1.0 : 0.6);   
        f.SetParLimits(4, 0.1,  8.0);                
        
        f.FixParameter(6, lo); f.FixParameter(7, hi);
    }

    void setParameters_File(TF1& f, const string& keyBase, double ek, double Z, double max, double lo, double hi) {
        setParameters_Init(f, Z, max, lo, hi); // Defaults
        if (!splineFile_) return;
        bool isHe = (Z < 2.5);

        for (int i = 0; i < f.GetNpar(); ++i) {
            string pName = f.GetParName(i);
            if (pName == "Norm" || pName == "xmin" || pName == "xmax") continue;
            
            double init = getInitialValue(keyBase, pName, ek);
            if (init > 0) {
                f.SetParameter(i, init);
                if (pName == "Peak") {
                    double w = isHe ? 0.5 : 0.3;
                    f.SetParLimits(i, init - w, init + w);
                } else {
                    double range = (isHe || pName.find("Alpha")!=string::npos) ? 0.6 : 0.3;
                    f.SetParLimits(i, init*(1.0 - range), init*(1.0 + range));
                }
            }
        }
    }
};

static TF1* BuildEGE(double l, double h){ return new TF1("fEGE", funcExpGausExp, l, h, 8); }

// --- UPGRADED: Smart Range Finding Logic ---
struct FitStep { 
    int L, R; 
    int sum; 
    int diff; 
};

static TF1* doFit(TH1D* h, const vector<std::pair<string,int>>& pars, FitParameterManager& pm,
                  const string& keyBase, double ek, double Z, double maxH, double& lo, double& hi) 
{
    TF1* f = nullptr;
    TF1* best_f = nullptr;
    double best_chi2 = 1e9;
    double best_lo = lo, best_hi = hi;

    bool isHe = (Z < 2.5);
    double stepL = isHe ? 0.01 : 0.05;
    double stepR = isHe ? 0.01 : 0.05;
    
    // 1. Generate Scan Plan
    // Constraints: Try Left fixed/Right varying, but don't let them be too asymmetric
    vector<FitStep> steps;
    for (int iL = 0; iL <= 25; ++iL) {
        for (int iR = 0; iR <= 25; ++iR) {
            if (std::abs(iL - iR) > 21 ) continue; // Enforce Symmetry: Difference <= 3 steps
            steps.push_back({iL, iR, iL + iR, std::abs(iL - iR)});
        }
    }

    // 2. Sort Strategy
    // Priority 1: Minimum shrinkage (L+R smallest) -> keeps most data
    // Priority 2: Symmetry (abs(L-R) smallest)
    // Priority 3: Prefer Right Shrink for He (since Right is problematic) -> if equal, prefer larger R
    std::sort(steps.begin(), steps.end(), [&](const FitStep& a, const FitStep& b){
        if (a.sum != b.sum) return a.sum < b.sum;       // Widest range first
        if (a.diff != b.diff) return a.diff < b.diff;   // Most symmetric first
        return a.R > b.R;                               // If same, try shrinking Right side first (for He)
    });

    double limitLo = Z - 0.1;
    double limitHi = Z + 0.05;

    for (const auto& s : steps) {
        double curr_lo = lo + s.L * stepL;
        double curr_hi = hi - s.R * stepR;

        if (curr_lo > limitLo || curr_hi < limitHi) continue;

        if (f) delete f;
        f = BuildEGE(1.5, 9.5);
        for (auto& p : pars) f->SetParName(p.second, p.first.c_str());
        
        pm.setParameters(*f, keyBase, ek, Z, maxH, curr_lo, curr_hi);
        h->Fit(f, "RQ0", "", curr_lo, curr_hi);

        auto badPars = pm.checkBoundaries(*f, Form("Step_L%d_R%d", s.L, s.R), ek);
        
        // Retry boundary logic
        if (!badPars.empty()) {
            for (int idx : badPars) {
                string pName = f->GetParName(idx);
                double init = pm.getInitialValue(keyBase, pName, ek);
                if (init > 0) {
                    double l, u; f->GetParLimits(idx, l, u);
                    double relax = (pName=="Peak") ? 0.4 : (init * 0.3);
                    f->SetParLimits(idx, l - relax, u + relax);
                }
            }
            h->Fit(f, "RQ0", "", curr_lo, curr_hi);
            badPars = pm.checkBoundaries(*f, "Retry", ek);
        }

        double chi2 = (f->GetNDF() > 0) ? f->GetChisquare() / f->GetNDF() : 999;
        
        // Threshold
        double thres = (isHe && keyBase.find("L1Template") != string::npos) ? 3.0 : 3.0;

        // Success: Found wide, symmetric range with good Chi2 -> STOP
        if (badPars.empty() && chi2 < thres) {
            lo = curr_lo; hi = curr_hi;
            if (best_f) delete best_f;
            return f; 
        }

        // Keep best fallback
        if (badPars.empty() && chi2 < best_chi2) {
            best_chi2 = chi2;
            best_lo = curr_lo; best_hi = curr_hi;
            if (best_f) delete best_f;
            best_f = (TF1*)f->Clone();
        }
    }

    if (f) delete f;
    if (best_f) {
        lo = best_lo; hi = best_hi;
        cout << "  [BestFallback] " << keyBase << " Chi2=" << best_chi2 
             << " Range=[" << lo << "," << hi << "]" << endl;
        return best_f;
    }
    
    return BuildEGE(lo, hi); 
}

void chargeHistFit(
    //const string& histFile = "/eos/ams/group/ihep/zixuan/filter/basic_L1Q2p5to8p8.root",
    const string& histFile = "/eos/user/z/zixuan/Isotope/Add/Be_frag4_withBkg_NoTune_full.root",
    const string& pdfOut = "/eos/user/z/zixuan/Isotope/ChargeFit/withBkg_ChargeFits_HeToOxy_NoTune_iter1.pdf",
    const string& histOut = "/eos/user/z/zixuan/Isotope/ChargeFit/withBkg_ChargeFitParams_HeToOxy_NoTune_iter1.root",
    /*
    const string& histFile = "/eos/user/z/zixuan/Isotope/PureChargeTemp/PureChargeTemplates_UnbiasedL1Inner.root",
    const string& pdfOut = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFits_PureL1.pdf",
    const string& histOut = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFits_PureL1.root",
    */
    int rebin = 2, bool firstFit = false, bool = false) 
{
    const vector<string> chains = {"UnbiasedL1Inner", "L1Inner"};
    const vector<string> nuclei = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen","Oxygen"}; 
    //const vector<string> nuclei = {"Carbon"}; 
    const map<string,double> chargeZ = {{"Helium",2.0}, {"Lithium",3.0}, {"Beryllium",4.0}, {"Boron",5.0}, {"Carbon",6.0}, {"Nitrogen",7.0}, {"Oxygen",8.0}};
    
    const vector<string> types = {"L1Template", "L2Template"}; 
    const vector<string> dets = {"TOF", "NaF", "AGL"};
    const map<string,int> detColor = {{"TOF", kRed}, {"NaF", kBlue}, {"AGL", kGreen+2}};
    const map<string,int> typeMarker = {{"L1Template", 20}, {"L2Template", 21}};

    const vector<std::pair<string,int>> EGE_p = {{"Peak",0},{"SigmaL",1},{"AlphaL",2},{"SigmaR",3},{"AlphaR",4},{"Norm",5},{"xmin",6},{"xmax",7}};

    const string splinePath = "/eos/user/z/zixuan/Isotope/ChargeFit/smooth/withBkg_ChargeFitParamsSmooth_HeToOxy_NoTune_iter0.root";
    const string oriPath = "/eos/user/z/zixuan/Isotope/ChargeFit/withBkg_ChargeFitParams_HeToOxy_NoTune_iter0.root";
    FitParameterManager pm(firstFit, splinePath, oriPath);

    unique_ptr<TFile> fin(TFile::Open(histFile.c_str()));
    unique_ptr<TFile> fout(TFile::Open(histOut.c_str(), "RECREATE"));
    if (!fin || fin->IsZombie() || !fout || fout->IsZombie()) return;

    unique_ptr<TCanvas> c(new TCanvas("c","",900,700));
    c->Print((pdfOut + "[").c_str(), "pdf");

    struct ParSet { bool init=false; int n=0; vector<double> x; map<string,TH1D*> h; };
    map<string, ParSet> store;

    for (auto& chain : chains)
    for (auto& elem : nuclei)
    for (auto& det : dets)
    for (auto& type : types) {
        string key = chain + "_BKG_H4_" + elem + "_" + type + "_" + det;
        //string key = "h2d_PureQTemp_" + elem + "_" + det;
        TH2* h2 = dynamic_cast<TH2*>(fin->Get(key.c_str()));
        if (!h2) { if(DBG.infoOpen) cout << "[miss] " << key << endl; continue; }

        string baseKey = chain + "_" + elem + "_" + det + "_" + type;
        ParSet& PS = store[baseKey];
        if (!PS.init) {
            PS.n = h2->GetNbinsY();
            PS.x.resize(PS.n+1);
            for(int i=1;i<=PS.n;++i) PS.x[i-1] = h2->GetYaxis()->GetBinLowEdge(i);
            PS.x[PS.n] = h2->GetYaxis()->GetBinUpEdge(PS.n);
            
            auto mk = [&](string p){
                string n = baseKey + "_EGE_" + p;
                TH1D* h = new TH1D(n.c_str(), (n+";E_{k}/n;"+p).c_str(), PS.n, PS.x.data());
                h->SetLineColor(detColor.at(det)); 
                h->SetMarkerStyle(typeMarker.at(type)); 
                h->SetLineWidth(2);
                PS.h[n] = h;
            };
            for(auto& p : EGE_p) mk(p.first); mk("Chi2NDF");
            PS.init = true;
        }

        for (int y=1; y<=PS.n; ++y) {
            double ek = h2->GetYaxis()->GetBinCenter(y);
            double eklow = h2->GetYaxis()->GetBinLowEdge(y);
            double ekhigh = h2->GetYaxis()->GetBinLowEdge(y+1);
            if (!passEnergyWindow(det, ek)) continue;

            double Z = chargeZ.at(elem);
            
            TH1D* h1 = dynamic_cast<TH1D*>(h2->ProjectionX(Form("%s_py%d",key.c_str(),y), y, y));
            if (!h1 || h1->GetMaximum() < 16) { delete h1; continue; }
            if (h1->GetMaximum()<100) h1->Rebin(rebin);
            h1->Sumw2();

            double lo, hi;
            
            findFitRange(h1, Z, lo, hi);
            TF1* fEGE = doFit(h1, EGE_p, pm, baseKey, ek, Z, h1->GetMaximum(), lo, hi);

            int b = PS.h[baseKey+"_EGE_Peak"]->FindBin(ek);
            for(auto& p : EGE_p) {
                TH1D* hh = PS.h[baseKey+"_EGE_"+p.first];
                hh->SetBinContent(b, fEGE->GetParameter(p.second));
                hh->SetBinError(b, p.second == 0 ? 1.*fEGE->GetParError(p.second) : 1.*fEGE->GetParError(p.second));
            }
            PS.h[baseKey+"_EGE_Chi2NDF"]->SetBinContent(b, (fEGE->GetNDF()>0 ? fEGE->GetChisquare()/fEGE->GetNDF() : 0));

            // Drawing
            h1->SetTitle(Form("%s %s | %s %s | [%.2f-%.2f]GeV", chain.c_str(), type.c_str(), elem.c_str(), det.c_str(), eklow, ekhigh));
            h1->SetMarkerStyle(20); h1->SetMarkerSize(0.8); h1->SetLineColor(kBlack);
            h1->GetXaxis()->SetRangeUser(Z == 8 ? Z - 1.5 : Z-2.0, Z > 2 ? Z+2.5 : Z+0.6); 
            
            c->cd(); c->SetLogy(Z==2?0:1);
            h1->Draw("E");
            fEGE->SetNpx(2000);
            fEGE->SetLineColor(kGreen+2); fEGE->SetLineWidth(3); fEGE->Draw("same");
            drawFitRangeLines(lo, hi, kGreen+2);

            TLegend leg(0.7, 0.7, 0.88, 0.88); leg.SetBorderSize(0);leg.SetFillStyle(0);
            leg.AddEntry(h1, "Data", "ep");
            leg.AddEntry(fEGE, Form("EGE #chi^{2}=%.1f", fEGE->GetChisquare()/fEGE->GetNDF()), "l");
            leg.Draw();
            
            TLatex lt; lt.SetNDC(); lt.SetTextFont(62); lt.SetTextSize(0.028); lt.SetTextAlign(12);
            double tx = 0.16, ty = 0.86, tdy = 0.032;
            lt.SetTextColor(kGreen+2);
            lt.DrawLatex(tx, ty, "ExpGausExp:"); ty -= tdy;
            lt.DrawLatex(tx+0.02, ty, Form("Peak = %.4f #pm %.4f", fEGE->GetParameter(0), fEGE->GetParError(0))); ty -= tdy;
            lt.DrawLatex(tx+0.02, ty, Form("#sigma_{L} = %.4f #pm %.4f", fEGE->GetParameter(1), fEGE->GetParError(1))); ty -= tdy;
            lt.DrawLatex(tx+0.02, ty, Form("#sigma_{R} = %.4f #pm %.4f", fEGE->GetParameter(3), fEGE->GetParError(3))); ty -= tdy;
            lt.DrawLatex(tx+0.02, ty, Form("#alpha_{L} = %.4f #pm %.4f", fEGE->GetParameter(2), fEGE->GetParError(2))); ty -= tdy;
            lt.DrawLatex(tx+0.02, ty, Form("#alpha_{R} = %.4f #pm %.4f", fEGE->GetParameter(4), fEGE->GetParError(4))); ty -= tdy;
            lt.DrawLatex(tx+0.02, ty, Form("#chi^{2}/ndf = %.2f", fEGE->GetNDF()>0?fEGE->GetChisquare()/fEGE->GetNDF():0)); ty -= tdy;
            lt.DrawLatex(tx+0.02, ty, Form("Range: %.2f-%.2f", lo, hi));

            c->Print(pdfOut.c_str(), "pdf");
            
            if (DBG.fitSummary) 
                cout << "[Fit] " << chain << " " << elem << " " << det << " Ek=" << ek 
                     << " EGE:" << fEGE->GetChisquare()/fEGE->GetNDF() << endl;

            delete fEGE; delete h1;
        }
    }

    c->Print((pdfOut + "]").c_str(), "pdf");
    fout->cd();
    for (auto& kv : store) for (auto& hk : kv.second.h) hk.second->Write();
    fout->Close();
    cout << "Done. Saved to " << histOut << endl;
}