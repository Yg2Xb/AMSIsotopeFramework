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

static inline void drawFitRangeLines(TCanvas* c, double fitLow, double fitHigh, int color) {
    if (!c) return;
    double yMin = c->GetUymin(), yMax = c->GetUymax();
    TLine* L1 = new TLine(fitLow, yMin, fitLow, yMax);  L1->SetLineColor(color);  L1->SetLineStyle(3);  L1->SetLineWidth(2);  L1->Draw("same");
    TLine* L2 = new TLine(fitHigh, yMin, fitHigh, yMax);L2->SetLineColor(color);  L2->SetLineStyle(3);  L2->SetLineWidth(2);  L2->Draw("same");
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
    if (highEdge < chargeValue + 0.35) highEdge = chargeValue + 0.35;
    if (DBG.rangeLog) {
        cout << std::fixed << std::setprecision(4)
             << "  [range] Z=" << chargeValue << " -> [" << lowEdge << ", " << highEdge << "]" << endl;
    }
}

static inline bool passEnergyWindow(const string& det, double ekLow) {
    if (det == "TOF") return ekLow >= 0.26 && ekLow < 1.55;
    if (det == "NaF") return ekLow >= 0.50 && ekLow <= 5.00;
    if (det == "AGL") return ekLow >= 2.00 && ekLow <= 18.00;
    return false;
}

static inline void addTextFmt(TPaveText& pt, Color_t color, const char* fmt, ...) {
    va_list ap; va_start(ap, fmt);
    char buf[512]; vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    auto* t = pt.AddText(buf);
    t->SetTextColor(color);
}

// Parameter initializer
class FitParameterManager {
public:
    explicit FitParameterManager(bool firstFit) : firstFit_(firstFit) {}
    void setParameters(TF1& f, const string& fitType,
                       const string&, const string&, const string&,
                       double, double charge, double histMax,
                       double fitLow, double fitHigh)
    {
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
    void checkBoundaries(const TF1& func, const string& fitName,
                         const string& elem, const string& det, const string& type,
                         double ekCen) const
    {
        if (!DBG.boundaryWarn) return;
        for (int i = 0; i < func.GetNpar(); ++i) {
            double v = func.GetParameter(i), lo, hi;
            func.GetParLimits(i, lo, hi);
            if (lo < hi) {
                const double eps = 1e-8;
                if (TMath::Abs(v - lo) < eps || TMath::Abs(v - hi) < eps) {
                    cout << "  [warn] boundary: " << fitName
                         << " par#" << i << " " << func.GetParName(i)
                         << " = " << v << " in [" << lo << "," << hi << "] "
                         << elem << "/" << det << "/" << type << " EkCen=" << ekCen << endl;
                }
            }
        }
    }
private:
    bool firstFit_;
};

// Build TF1 (no lambdas)
static TF1* BuildLG(double lo, double hi) { return new TF1("fLG",  langaufun,       lo, hi, 5); }
static TF1* BuildEGE(double lo, double hi){ return new TF1("fEGE", funcExpGausExp, lo, hi, 8); }

// One-pass fitter without templates/lambdas
static TF1* doFit(TH1D* h, TF1* (*builder)(double,double),
                  const vector<std::pair<string,int>>& parList,
                  FitParameterManager& pm, const string& fitType,
                  const string& chain, const string& elem, const string& det, const string& type,
                  double ekCen, double chargeZ, double histMax,
                  double& fitLow, double& fitHigh, int maxIter = 4, const char* fitOpt = "RQ0")
{
    TF1* f = nullptr;
    for (int iter = 0; iter < maxIter; ++iter) {
        if (f) { delete f; f = nullptr; }
        f = builder(fitLow, fitHigh);
        for (size_t i = 0; i < parList.size(); ++i)
            f->SetParName(parList[i].second, parList[i].first.c_str());

        pm.setParameters(*f, fitType, elem, det, type, ekCen, chargeZ, histMax, fitLow, fitHigh);
        h->Fit(f, fitOpt, "", fitLow, fitHigh);
        pm.checkBoundaries(*f, fitType + (iter==0 ? "" : ("-" + std::to_string(iter+1))), elem, det, type, ekCen);

        double ndf = f->GetNDF();
        double chi2ndf = (ndf > 0) ? f->GetChisquare()/ndf : 1e9;
        if (chi2ndf < 3.0) break;

        // simple range adjust
        if (fitType == "LG") {
            if (iter == 0) { fitLow -= 0.05; }
            else { fitLow += 0.05; fitHigh -= 0.05; }
        } else { // EGE
            if (iter > 0) { fitLow += 0.10; fitHigh -= 0.10; }
        }
    }
    return f;
}

void chargeHistFit(
    const string& histFile = "/eos/user/z/zixuan/Isotope/Add/Be_all.root",
    const string& pdfOut  = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFits_BeToOxy.pdf",
    const string& histOut = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFitParams_BeToOxy.root",
    int rebin = 2,
    bool firstFit = true,
    bool listKeysOnce = false // kept for compatibility, unused to stay simple
) {
    const vector<string> chains = {"L1Inner", "UnbiasedL1Inner"};
    const vector<string> nuclei = {"Beryllium","Boron","Carbon","Nitrogen","Oxygen"};
    const map<string,double> chargeZ = {{"Beryllium",4.0}, {"Boron",5.0}, {"Carbon",6.0}, {"Nitrogen",7.0}, {"Oxygen",8.0}};
    const vector<string> types = {"L1QTemplate", "L2QTemplate"};
    const vector<string> dets = {"TOF","NaF","AGL"};
    const map<string,int> detColor = {{"TOF", kRed}, {"NaF", kBlue}, {"AGL", kGreen+2}};
    const map<string,int> typeMarker = {{"L1QTemplate", 20}, {"L2QTemplate", 21}};

    const vector<std::pair<string,int>> LG_params = {{"Width",0},{"MPV",1},{"Area",2},{"Sigma",3}};
    const vector<std::pair<string,int>> EGE_params = {{"Peak",0},{"SigmaL",1},{"AlphaL",2},{"SigmaR",3},{"AlphaR",4},{"Norm",5},{"xmin",6},{"xmax",7}};

    FitParameterManager pm(firstFit);

    unique_ptr<TFile> fin(TFile::Open(histFile.c_str()));
    printOpen("input", histFile, fin.get());
    if (!fin || fin->IsZombie()) return;

    unique_ptr<TFile> fout(TFile::Open(histOut.c_str(), "RECREATE"));
    printOpen("output", histOut, fout.get());
    if (!fout || fout->IsZombie()) return;

    unique_ptr<TCanvas> c(new TCanvas("c","",900,700));
    c->SetLogy(true);
    c->Print((pdfOut + "[").c_str(), "pdf");

    // Parameter hist store struct
    struct ParSet {
        bool init = false;
        int nb = 0;
        std::vector<double> edges;
        map<string, TH1D*> h;
    };
    map<string, ParSet> store;

    auto keyEDT = [](const string& e, const string& d, const string& t){ return e + "_" + d + "_" + t; };
    auto namePar = [](const string& e,const string& d,const string& t,const string& fit,const string& par){
        return e + "_" + d + "_" + t + "_" + fit + "_" + par;
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

        // init param hists on first time for this elem/det/type using Y-axis edges
        string edt = keyEDT(elem,det,type);
        ParSet& PS = store[edt];
        if (!PS.init) {
            int nY = h2->GetYaxis()->GetNbins();
            PS.nb = nY;
            PS.edges.resize(nY+1);
            for (int i=1;i<=nY;++i) PS.edges[i-1] = h2->GetYaxis()->GetBinLowEdge(i);
            PS.edges[nY] = h2->GetYaxis()->GetBinUpEdge(nY);

            auto mk = [&](const string& fit, const string& par){
                string nm = namePar(elem,det,type,fit,par);
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
            if (DBG.infoOpen) cout << "[init] " << edt << " with " << PS.nb << " Ek/n bins" << endl;
        }

        int nY = h2->GetYaxis()->GetNbins();
        for (int ybin=1; ybin<=nY; ++ybin) {
            double ekLow  = h2->GetYaxis()->GetBinLowEdge(ybin);
            double ekHigh = h2->GetYaxis()->GetBinUpEdge(ybin);
            if (!passEnergyWindow(det, ekLow)) continue;

            double ekCen = 0.5*(ekLow+ekHigh);
            string projName = key + Form("_projY%d", ybin);
            TH1D* h1 = dynamic_cast<TH1D*>(h2->ProjectionX(projName.c_str(), ybin, ybin));
            if (!h1) continue;
            if (h1->GetMaximum() < 50) { delete h1; continue; }
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

            // Fill parameter hists
            int b = store[edt].h[namePar(elem,det,type,"LG","Width")]->FindBin(ekCen);
            for (size_t iP=0;iP<LG_params.size();++iP) {
                const string& par = LG_params[iP].first; int ip = LG_params[iP].second;
                TH1D* hh = store[edt].h[namePar(elem,det,type,"LG",par)];
                hh->SetBinContent(b, fLG->GetParameter(ip));
                hh->SetBinError(b, fLG->GetParError(ip));
            }
            {
                TH1D* hchi = store[edt].h[namePar(elem,det,type,"LG","Chi2NDF")];
                hchi->SetBinContent(b, (fLG->GetNDF()>0)? fLG->GetChisquare()/fLG->GetNDF() : 0);
                hchi->SetBinError(b, 0);
            }

            for (size_t iP=0;iP<EGE_params.size();++iP) {
                const string& par = EGE_params[iP].first; int ip = EGE_params[iP].second;
                TH1D* hh = store[edt].h[namePar(elem,det,type,"EGE",par)];
                hh->SetBinContent(b, fEGE->GetParameter(ip));
                hh->SetBinError(b, fEGE->GetParError(ip));
            }
            {
                TH1D* hchi = store[edt].h[namePar(elem,det,type,"EGE","Chi2NDF")];
                hchi->SetBinContent(b, (fEGE->GetNDF()>0)? fEGE->GetChisquare()/fEGE->GetNDF() : 0);
                hchi->SetBinError(b, 0);
            }

            // Draw
            h1->SetStats(0);
            h1->SetTitle(Form("%s | %s | %s | %s | Ek/n: [%.3f, %.3f] GeV/n", chain.c_str(), type.c_str(), elem.c_str(), det.c_str(), ekLow, ekHigh));
            h1->GetXaxis()->SetTitle("Tracker Layer Q");
            h1->GetYaxis()->SetTitle("Counts");
            h1->GetYaxis()->SetTitleOffset(1.4);
            h1->SetMarkerStyle(20); h1->SetMarkerSize(0.9); h1->SetMarkerColor(kBlack); h1->SetLineColor(kBlack);
            h1->GetXaxis()->SetRangeUser(z - 1.2, z + 1.0);

            c->cd(); c->SetLogy(true);
            h1->Draw("E");

            fLG->SetRange(h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
            fLG->SetLineColor(kRed); fLG->SetLineWidth(3); fLG->SetLineStyle(1); fLG->Draw("same");
            fEGE->SetRange(h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
            fEGE->SetLineColor(kGreen+2); fEGE->SetLineWidth(3); fEGE->SetLineStyle(2); fEGE->Draw("same");

            drawFitRangeLines(c.get(), fitLowLG,  fitHighLG,  kRed);
            drawFitRangeLines(c.get(), fitLowEGE, fitHighEGE, kGreen+2);

            TLegend leg(0.66, 0.66, 0.93, 0.88);
            leg.SetTextSize(0.03); leg.SetBorderSize(0); leg.SetFillStyle(0);
            leg.AddEntry(h1, "Data", "ep");
            leg.AddEntry(fLG, "Landau-Gauss", "l");
            leg.AddEntry(fEGE, "ExpGausExp", "l");
            leg.Draw();

            TPaveText pt(0.15, 0.35, 0.47, 0.88, "NDC");
            pt.SetTextSize(0.033); pt.SetBorderSize(0); pt.SetFillStyle(0); pt.SetTextAlign(13);
            addTextFmt(pt, kRed,      "LandauGauss:");
            addTextFmt(pt, kRed,      "  MPV = %.3f #pm %.3f",  fLG->GetParameter(1), fLG->GetParError(1));
            addTextFmt(pt, kRed,      "  Width = %.3f #pm %.3f",fLG->GetParameter(0), fLG->GetParError(0));
            addTextFmt(pt, kRed,      "  Sigma = %.3f #pm %.3f",fLG->GetParameter(3), fLG->GetParError(3));
            addTextFmt(pt, kRed,      "  #chi^{2}/ndf = %.0f/%d = %.2f", fLG->GetChisquare(), (int)fLG->GetNDF(), (fLG->GetNDF()>0? fLG->GetChisquare()/fLG->GetNDF():0));
            addTextFmt(pt, kRed,      "  fit range: %.2f-%.2f", fitLowLG, fitHighLG);
            addTextFmt(pt, kGreen+2,  "ExpGausExp:");
            addTextFmt(pt, kGreen+2,  "  Peak = %.4f #pm %.4f",   fEGE->GetParameter(0), fEGE->GetParError(0));
            addTextFmt(pt, kGreen+2,  "  #sigma_{L} = %.4f #pm %.4f", fEGE->GetParameter(1), fEGE->GetParError(1));
            addTextFmt(pt, kGreen+2,  "  #sigma_{R} = %.4f #pm %.4f", fEGE->GetParameter(3), fEGE->GetParError(3));
            addTextFmt(pt, kGreen+2,  "  #alpha_{L} = %.4f #pm %.4f", fEGE->GetParameter(2), fEGE->GetParError(2));
            addTextFmt(pt, kGreen+2,  "  #alpha_{R} = %.4f #pm %.4f", fEGE->GetParameter(4), fEGE->GetParError(4));
            addTextFmt(pt, kGreen+2,  "  #chi^{2}/ndf = %.0f/%d = %.2f", fEGE->GetChisquare(), (int)fEGE->GetNDF(), (fEGE->GetNDF()>0? fEGE->GetChisquare()/fEGE->GetNDF():0));
            addTextFmt(pt, kGreen+2,  "  fit range: %.2f-%.2f", fitLowEGE, fitHighEGE);
            pt.Draw();

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