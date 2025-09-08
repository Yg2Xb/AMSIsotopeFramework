#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLine.h>
#include "../Tool.h"

using namespace AMS_Iso;
using namespace std;

void drawFitRangeLines(TCanvas* canvas, double yMax, double fitLow, double fitHigh,
                       int color = kGray+2, int style = 3, int width = 2) {
    if (!canvas) return;
    canvas->cd();
    double yMin = canvas->GetUymin();
    TLine* lowLine = new TLine(fitLow, yMin, fitLow, yMax);
    lowLine->SetLineColor(color);
    lowLine->SetLineStyle(style);
    lowLine->SetLineWidth(width);
    lowLine->Draw("same");
    TLine* highLine = new TLine(fitHigh, yMin, fitHigh, yMax);
    highLine->SetLineColor(color);
    highLine->SetLineStyle(style);
    highLine->SetLineWidth(width);
    highLine->Draw("same");
    canvas->Modified();
    canvas->Update();
}

double findContentLevel(double charge, TH1D* hist, double ratio, bool searchLeft) {
    if (!hist) return 0;
    int binMin = hist->GetXaxis()->FindBin(charge - 0.12);
    int binMax = hist->GetXaxis()->FindBin(charge + 0.12);
    int maxBin = binMin;
    double maxContent = hist->GetBinContent(binMin);
    for (int bin = binMin + 1; bin <= binMax; ++bin) {
        double content = hist->GetBinContent(bin);
        if (content > maxContent) {
            maxContent = content;
            maxBin = bin;
        }
    }
    double targetContent = maxContent * ratio;
    double targetX = hist->GetXaxis()->GetBinCenter(maxBin);
    if (searchLeft) {
        for (int bin = maxBin; bin >= 1; --bin) {
            if (hist->GetBinContent(bin) <= targetContent) {
                targetX = hist->GetXaxis()->GetBinCenter(bin);
                break;
            }
        }
    } else {
        for (int bin = maxBin; bin <= hist->GetNbinsX(); ++bin) {
            if (hist->GetBinContent(bin) <= targetContent) {
                targetX = hist->GetXaxis()->GetBinCenter(bin);
                break;
            }
        }
    }
    return targetX;
}

void findFitRange(TH1D* hist, double chargeValue, double leftratio, double rightratio, double& lowEdge, double& highEdge) {
    if (!hist) return;
    lowEdge = findContentLevel(chargeValue, hist, leftratio, true);
    highEdge = findContentLevel(chargeValue, hist, rightratio, false);
    if (lowEdge > chargeValue-0.24) lowEdge = chargeValue-0.24;
    if (highEdge < chargeValue+0.35) highEdge = chargeValue+0.35;
}

class FitParameterManager {
public:
    FitParameterManager(bool firstFit, const string& splineFilePath, const string& histOriFilePath)
        : m_firstFit(firstFit), m_splineFile(nullptr), m_histFile(nullptr) {
        if (!m_firstFit) {
            m_splineFile = TFile::Open(splineFilePath.c_str());
            if (!m_splineFile || m_splineFile->IsZombie()) {
                cout << "ERROR: Could not open spline file: " << splineFilePath << endl;
                m_splineFile = nullptr;
            }
            m_histFile = TFile::Open(histOriFilePath.c_str());
            if (!m_histFile || m_histFile->IsZombie()) {
                cout << "ERROR: Could not open hist_ori file: " << histOriFilePath << endl;
                m_histFile = nullptr;
            }
        }
    }

    ~FitParameterManager() {
        if (m_splineFile) { m_splineFile->Close(); delete m_splineFile; }
        if (m_histFile) { m_histFile->Close(); delete m_histFile; }
    }

    void setParameters(TF1& func, const string& fitType, const string& elem, const string& det, const string& type, double ekCen, double charge, double histMax, double fitLow = 0, double fitHigh = 0) {
        if (m_firstFit) {
            setParameters_firstFit(func, fitType, charge, histMax, det, fitLow, fitHigh);
        } else {
            setParameters_fromFile(func, fitType, elem, det, type, ekCen, charge, histMax, fitLow, fitHigh);
        }
    }
private:
    bool m_firstFit;
    TFile* m_splineFile;
    TFile* m_histFile;
    map<string, double> m_alphaLAverageCache;

    void setParameters_firstFit(TF1& func, const string& fitType, double charge, double histMax, const string& det, double fitLow, double fitHigh) {
        if (fitType == "LG") {
            func.SetParameters(0.05, charge, histMax, 0.2, 1.0);
            func.FixParameter(4, 1.0);
            func.SetParLimits(0, 0.005, 0.11);
            func.SetParLimits(1, charge - 0.15, charge + 0.18);
            func.SetParLimits(3, 0.06, 0.5);
        } else if (fitType == "EGE") {
            func.SetParameters(charge, 0.20, 2.0, 0.20, 2.0, histMax, fitLow, fitHigh);
            func.SetParLimits(1, 0.01, 0.5);
            func.SetParLimits(2, 0.5, 4.0);
            func.SetParLimits(3, 0.01, 0.5);
            func.SetParLimits(4, 0.5, 4.0);
            func.FixParameter(6, fitLow);
            func.FixParameter(7, fitHigh);
        }
    }
    void setParameters_fromFile(TF1& func, const string& fitType, const string& elem, const string& det, const string& type,double ekCen, double charge, double histMax, double fitLow, double fitHigh) {
        if (fitType == "LG") {
            func.SetParameters(0.05, charge, histMax, 0.2, 1.0);
            func.FixParameter(4, 1.0);
        } else if (fitType == "EGE") {
            func.SetParameters(charge, 0.20, 2.0, 0.20, 2.0, histMax, fitLow, fitHigh);
            func.FixParameter(6, fitLow);
            func.FixParameter(7, fitHigh);
        }
        if (!m_splineFile || !m_histFile) {
            setParameters_firstFit(func, fitType, charge, histMax, det, fitLow, fitHigh);
            return;
        }
        for (int i = 0; i < func.GetNpar(); ++i) {
            string parName = func.GetParName(i);
            if (fitType == "EGE" && parName == "AlphaL") { //elem != "Beryllium" && type != "L1Temp"
                double avgAlphaL = getAlphaLAverage(elem, det, type);
                if (avgAlphaL > 0) {
                    func.FixParameter(i, avgAlphaL);
                    //func.SetParLimits(i, 0.9 * avgAlphaL, 1.1 * avgAlphaL);
                } else {
                    func.FixParameter(i, 3.0);
                }
                continue;
            }
            string tf1Name = elem + "_Combined_" + type + "_" + fitType + "_" + parName + "_SplineSmooth";
            TF1* f_param = (TF1*)m_splineFile->Get(tf1Name.c_str());
            if (f_param) {
                double initialVal = f_param->Eval(ekCen);
                func.SetParameter(i, initialVal);
                func.SetParLimits(i, 0.75 * initialVal, 1.25 * initialVal);
                if(parName == "AlphaR" || parName == "Width")
                {
                    func.SetParameter(i, initialVal);
                    func.SetParLimits(i, 0.5 * initialVal, 1.5 * initialVal);
                }
            } else {
                double min_lim, max_lim;
                func.GetParLimits(i, min_lim, max_lim);
                if (min_lim >= max_lim) continue;
                if (fitType == "LG") {
                    if (parName == "Width") func.SetParLimits(i, 0.005, 0.11);
                    else if (parName == "MPV") func.SetParLimits(i, charge - 0.15, charge + 0.18);
                    else if (parName == "Sigma") func.SetParLimits(i, 0.06, 0.5);
                } else if (fitType == "EGE") {
                    if (parName == "SigmaL") func.SetParLimits(i, 0.01, 0.5);
                    else if (parName == "SigmaR") func.SetParLimits(i, 0.01, 0.5);
                    else if (parName == "AlphaR") func.SetParLimits(i, 0.5, 4.0);
                }
            }
        }
    }
    double getAlphaLAverage(const string& elem, const string& det, const string& type) {
        string cacheKey = elem + "_" + det + "_" + type;
        if (m_alphaLAverageCache.count(cacheKey)) return m_alphaLAverageCache[cacheKey];
        if (!m_histFile) return 0.0;
        string histName = elem + "_" + det + "_" + type + "_EGE_AlphaL";
        TH1D* h = (TH1D*)m_histFile->Get(histName.c_str());
        if (!h) return 0.0;
        double rangeLow = 0, rangeHigh = 0;
        if (det == "TOF") { rangeLow = 0.42; rangeHigh = 1.55; }
        else if (det == "NaF") { rangeLow = 0.86; rangeHigh = 4.91; }
        else if (det == "AGL") { rangeLow = 2.88; rangeHigh = 16.3; }
        else return 0.0;
        int binLow = h->GetXaxis()->FindBin(rangeLow+0.01);
        int binHigh = h->GetXaxis()->FindBin(rangeHigh-0.01);
        double sum = 0; int count = 0;
        for (int bin = binLow; bin <= binHigh; ++bin) {
            double content = h->GetBinContent(bin);
            if (content != 0) { sum += content; count++; }
        }
        double average = (count > 0) ? sum / count : 0.0;
        m_alphaLAverageCache[cacheKey] = average;
        return average;
    }
};

void checkFitBoundaries(const TF1& func, const string& fitName, const string& elem, const string& det, const string& type, double ekCen) {
    for (int i = 0; i < func.GetNpar(); ++i) {
        double val = func.GetParameter(i);
        double min_lim, max_lim;
        func.GetParLimits(i, min_lim, max_lim);
        if (min_lim < max_lim) {
            const double epsilon = 1e-6;
            if (TMath::Abs(val - min_lim) < epsilon || TMath::Abs(val - max_lim) < epsilon) {
                cout << "WARNING: Fit parameter at boundary! -> "
                     << "Fit: " << fitName << ", "
                     << "Param: " << func.GetParName(i) << " (" << i << "), "
                     << "Value: " << val << ", "
                     << "Limits: [" << min_lim << ", " << max_lim << "], "
                     << "for " << elem << "/" << det << "/" << type << " @ Ek~" << ekCen << " GeV/n"
                     << endl;
            }
        }
    }
}

template <typename TF1Builder>
unique_ptr<TF1> multiFit(
    TH1D* h1,
    TF1Builder makeTF1,
    vector<pair<string, int>> parList,
    FitParameterManager& paramManager,
    const string& fitType,
    const string& elem, const string& det, const string& type,
    double ekCen, double charge, double histMax,
    double& fitLow, double& fitHigh,
    function<void(int, double&, double&)> adjustRange,
    int maxTries = 4,
    const char* fitOpt = "RQ0"
) {
    unique_ptr<TF1> func;
    for (int iter = 0; iter < maxTries; ++iter) {
        func.reset(makeTF1(fitLow, fitHigh)); // 每次全新建
        for (const auto& param : parList)
            func->SetParName(param.second, param.first.c_str());
        paramManager.setParameters(*func, fitType, elem, det, type, ekCen, charge, histMax, fitLow, fitHigh);
        h1->Fit(func.get(), fitOpt, "", fitLow, fitHigh);
        checkFitBoundaries(*func, fitType + (iter==0?"":("-"+to_string(iter+1))), elem, det, type, ekCen);

        double chi2ndf = func->GetNDF() > 0 ? func->GetChisquare()/func->GetNDF() : 1e9;
        if (chi2ndf < 3.0) break;
        adjustRange(iter, fitLow, fitHigh); // 调整区间
    }
    return func;
}

// ----- 主程序 -----
void chargeHistFit(
    const string& histFile = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/ChargeTemp_Hist_0p8.root",
    const string& pdfOut = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/allFits_0p8_r_.pdf",
    const string& histOut = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/allFitHist_0p8_r_.root",
    int rebin = 2,
    bool firstFit = false)
{
    vector<string> nuclei = {"Beryllium", "Boron", "Nitrogen"};
    map<string, double> charge = {{"Beryllium",4.0}, {"Boron",5.0}, {"Nitrogen",7.0}};
    vector<string> types = {"L1Temp", "L2Temp", "L1Temp_unb", "L2Temp_unb"};
    vector<string> dets = {"TOF", "NaF", "AGL"};
    const int nbins = Binning::NarrowBins.size()-1;
    const double* bins = Binning::NarrowBins.data();
    map<string, TH1D*> hists;
    map<string, int> detColor = {{"TOF", kRed}, {"NaF", kBlue}, {"AGL", kGreen+2}};
    map<string, int> typeMarker = {{"L1Temp", 20}, {"L2Temp", 21}, {"L1Temp_unb",22}, {"L2Temp_unb",23}};
    vector<pair<string, int>> LG_paramList = {{"Width", 0}, {"MPV", 1}, {"Area", 2}, {"Sigma", 3}};
    vector<pair<string, int>> EGE_paramList = {
        {"Peak", 0}, {"SigmaL", 1}, {"AlphaL", 2},
        {"SigmaR", 3}, {"AlphaR", 4}, {"Norm", 5}, {"xmin", 6}, {"xmax", 7}
    };
    const string splineFilePath = "/eos/user/z/zixuan/Isotope/chargeTempFit/allFitHistSplineSmooth_0p8.root";
    const string histOriFilePath = "/eos/user/z/zixuan/Isotope/chargeTempFit/allFitHist_0p8_ori.root";
    FitParameterManager paramManager(firstFit, splineFilePath, histOriFilePath);

    auto regHist = [&](const string& fit, const string& par, const string& elem, const string& det, const string& type) {
        string hname = elem+"_"+det+"_"+type+"_"+fit+"_"+par;
        if(!hists.count(hname)) {
            string htitle = hname+";E_{k}/n [GeV/n];"+par;
            auto h = new TH1D(hname.c_str(), htitle.c_str(), nbins, bins);
            h->SetLineColor(detColor[det]);
            h->SetMarkerColor(detColor[det]);
            h->SetMarkerStyle(typeMarker[type]);
            h->SetLineWidth(2);
            hists[hname] = h;
        }
        return hists[hname];
    };

    TFile* file = new TFile(histFile.c_str());
    if (!file || file->IsZombie()) { cout << "Error opening file: " << histFile << endl; return; }
    TFile* ohist = new TFile(histOut.c_str(), "RECREATE");
    TCanvas* c = new TCanvas("c", "", 800, 600);
    c->Print((pdfOut+"[").c_str(),"pdf");

    for(const auto& elem : nuclei)
    for(const auto& det : dets)
    for(const auto& type : types)
    {
        string hname2D = Form("%s_%s_%s_Bkg", type.c_str(), elem.c_str(), det.c_str());
        TH2D* h2 = (TH2D*)file->Get(hname2D.c_str());
        if(!h2) continue;

        int nY = h2->GetYaxis()->GetNbins();
        for(int ybin=1; ybin<=nY; ++ybin) {
            double ekLow = h2->GetYaxis()->GetBinLowEdge(ybin);
            double ekHigh= h2->GetYaxis()->GetBinUpEdge(ybin);
            double ekCen = 0.5*(ekLow+ekHigh);

            if ((det=="TOF" && ekLow>1.55) || (det=="NaF" && (ekLow<0.86 || ekLow>4.91)) || (det=="AGL" && ekLow<2.88)) continue;

            TH1D* h1 = h2->ProjectionX(Form("%s_proj%d", hname2D.c_str(), ybin), ybin, ybin);
            if(!h1||h1->GetMaximum()<50) { delete h1; continue; }
            if(rebin>1) h1->Rebin(rebin);
            h1->Sumw2();

            double fitLow_LG, fitHigh_LG;
            findFitRange(h1, charge[elem], 0.2, 0.2, fitLow_LG, fitHigh_LG);
            auto LG_adjust = [&](int iter, double& low, double& high) {
                if(iter==0) { /* nothing */ }
                else if(iter==1) {
                    findFitRange(h1, charge[elem], 0.4, 0.2, low, high);
                } else {
                    low += 0.05; high -= 0.05;
                }
            };
            auto LGBuilder = [](double low, double high) {
                return new TF1("fLG", langaufun, low, high, 5);
            };
            auto fLG = multiFit(
                h1, LGBuilder, LG_paramList, paramManager, "LG",
                elem, det, type, ekCen, charge[elem], h1->GetMaximum(),
                fitLow_LG, fitHigh_LG, LG_adjust
            );

            for (auto& kv : LG_paramList) {
                int idx = kv.second;
                TH1D* h = regHist("LG", kv.first, elem, det, type);
                int bin = h->FindBin(ekCen);
                h->SetBinContent(bin, fLG->GetParameter(idx));
                h->SetBinError(bin, fLG->GetParError(idx));
            }
            TH1D* hchi = regHist("LG", "Chi2NDF", elem, det, type);
            int chibin = hchi->FindBin(ekCen);
            hchi->SetBinContent(chibin, fLG->GetChisquare()/fLG->GetNDF());
            hchi->SetBinError(chibin, 0);

            double fitLow_EGE, fitHigh_EGE;
            double fitr = ((det=="NaF" && ekLow>1.17 && ekLow<2.01) || (det=="TOF" && (ekLow==0.61 || ekLow==1.00))) ? 0.2 : 0.1;
            findFitRange(h1, charge[elem], fitr, fitr, fitLow_EGE, fitHigh_EGE);
            auto EGE_adjust = [](int iter, double& low, double& high) {
                if(iter>0) { low += 0.1; high -= 0.1; }
            };
            auto EGEBuilder = [](double low, double high) {
                return new TF1("fEGE", funcExpGausExp, low, high, 8);
            };
            auto fEGE = multiFit(
                h1, EGEBuilder, EGE_paramList, paramManager, "EGE",
                elem, det, type, ekCen, charge[elem], h1->GetMaximum(),
                fitLow_EGE, fitHigh_EGE, EGE_adjust
            );

            for (auto& kv : EGE_paramList) {
                int idx = kv.second;
                TH1D* h = regHist("EGE", kv.first, elem, det, type);
                int bin = h->FindBin(ekCen);
                h->SetBinContent(bin, fEGE->GetParameter(idx));
                h->SetBinError(bin, fEGE->GetParError(idx));
            }
            TH1D* hchi2 = regHist("EGE", "Chi2NDF", elem, det, type);
            int chibin2 = hchi2->FindBin(ekCen);
            hchi2->SetBinContent(chibin2, fEGE->GetChisquare()/fEGE->GetNDF());
            hchi2->SetBinError(chibin2, 0);

            h1->SetStats(0);
            h1->SetTitle(Form("%s, %s, %s, E_{k}/n: [%.2f, %.2f] GeV/n", type.c_str(), elem.c_str(), det.c_str(), ekLow, ekHigh));
            h1->GetXaxis()->SetTitle("Tracker Layer Q");
            h1->GetYaxis()->SetTitle("Counts");
            h1->GetYaxis()->SetTitleOffset(1.4);
            h1->SetMarkerStyle(20); h1->SetMarkerSize(0.9); h1->SetMarkerColor(kBlack); h1->SetLineColor(kBlack);
            float xlow = charge[elem] > 4 ? 1.2 : 1;
            h1->GetXaxis()->SetRangeUser(charge[elem]-xlow, charge[elem]+1.);
			h1->Draw("E");
            fLG->SetRange(h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
            fLG->SetLineColor(kRed); fLG->SetLineWidth(3); fLG->SetLineStyle(1); fLG->Draw("same");
            fEGE->SetRange(h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
            fEGE->SetLineColor(kGreen+2); fEGE->SetLineWidth(3); fEGE->SetLineStyle(2); fEGE->Draw("same");
            c->Update();
            double yMax = c->GetUymax();
            drawFitRangeLines(c, yMax, fitLow_LG, fitHigh_LG, kRed);
            drawFitRangeLines(c, yMax, fitLow_EGE, fitHigh_EGE, kGreen+2);

            TLegend leg(0.68, 0.68, 0.93, 0.88);
            leg.SetTextSize(0.03); leg.SetBorderSize(0); leg.SetFillStyle(0);
            leg.AddEntry(h1, "Data", "ep");
            leg.AddEntry(fLG.get(), "Landau-Gauss", "l");
            leg.AddEntry(fEGE.get(), "ExpGausExp", "l");
            leg.Draw();

            TPaveText pt(0.15, 0.39, 0.47, 0.88, "NDC");
            pt.SetTextSize(0.033); pt.SetBorderSize(0); pt.SetFillStyle(0); pt.SetTextAlign(13);
            TText* t;
            t = pt.AddText(Form("LandauGauss:")); t->SetTextColor(kRed);
            t = pt.AddText(Form("  MPV = %.3f #pm %.3f", fLG->GetParameter(1), fLG->GetParError(1))); t->SetTextColor(kRed);
            t = pt.AddText(Form("  Width = %.3f #pm %.3f", fLG->GetParameter(0), fLG->GetParError(0))); t->SetTextColor(kRed);
            t = pt.AddText(Form("  Sigma = %.3f #pm %.3f", fLG->GetParameter(3), fLG->GetParError(3))); t->SetTextColor(kRed);
            t = pt.AddText(Form("  #chi^{2}/ndf = %.0f/%d = %.2f", fLG->GetChisquare(), (int)fLG->GetNDF(), fLG->GetChisquare()/fLG->GetNDF())); t->SetTextColor(kRed);
            t = pt.AddText(Form("  fit range: %.2f-%.2f", fitLow_LG, fitHigh_LG)); t->SetTextColor(kRed);
            t = pt.AddText(Form("ExpGausExp:")); t->SetTextColor(kGreen+2);
            t = pt.AddText(Form("  Peak = %.4f #pm %.4f", fEGE->GetParameter(0), fEGE->GetParError(0))); t->SetTextColor(kGreen+2);
            t = pt.AddText(Form("  #sigma_{L} = %.4f #pm %.4f", fEGE->GetParameter(1), fEGE->GetParError(1))); t->SetTextColor(kGreen+2);
            t = pt.AddText(Form("  #sigma_{R} = %.4f #pm %.4f", fEGE->GetParameter(3), fEGE->GetParError(3))); t->SetTextColor(kGreen+2);
            t = pt.AddText(Form("  #alpha_{L} = %.4f #pm %.4f", fEGE->GetParameter(2), fEGE->GetParError(2))); t->SetTextColor(kGreen+2);
            t = pt.AddText(Form("  #alpha_{R} = %.4f #pm %.4f", fEGE->GetParameter(4), fEGE->GetParError(4))); t->SetTextColor(kGreen+2);
            t = pt.AddText(Form("  #chi^{2}/ndf = %.0f/%d = %.2f", fEGE->GetChisquare(), (int)fEGE->GetNDF(), fEGE->GetChisquare()/fEGE->GetNDF())); t->SetTextColor(kGreen+2);
            t = pt.AddText(Form("  fit range: %.2f-%.2f", fitLow_EGE, fitHigh_EGE)); t->SetTextColor(kGreen+2);
            pt.Draw();

            c->Print(pdfOut.c_str(), "pdf");
            delete h1;
        }
    }
    c->Print((pdfOut+"]").c_str(),"pdf");
    ohist->cd();
    for(auto& kv : hists) kv.second->Write();
    ohist->Close();
    file->Close();
    delete file;
    delete ohist;
    delete c;
    cout << "All parameter histograms saved to " << histOut << endl;
}