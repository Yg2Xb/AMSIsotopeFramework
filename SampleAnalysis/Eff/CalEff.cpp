#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <tuple>
#include <numeric>
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TF1.h"
#include "../Tool.h" 

using namespace std;

void SetStyle(TH1* h, int color, int marker) {
    if (!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.8);
    h->SetLineWidth(2);
    h->SetStats(0); 
    h->GetYaxis()->SetTitle("Efficiency");
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.055);
    h->SetTitle(""); 
}

void GetGlobalMinMax(const vector<TH1D*>& hists, double xMin, double xMax, double& yMin, double& yMax) {
    yMax = -1e9; yMin = 1e9;
    bool foundAny = false;
    for (auto h : hists) {
        if (!h) continue;
        int bStart = h->FindBin(xMin);
        int bEnd = h->FindBin(xMax);
        for (int b = bStart; b <= bEnd; ++b) {
            double c = h->GetBinContent(b);
            double e = h->GetBinError(b);
            if (c + e > yMax) yMax = c + e; 
            if (c > 1e-9) {
                if (c - e < yMin) yMin = c - e; 
                foundAny = true;
            }
        }
    }
    if (yMax == -1e9) yMax = 1.0;
    if (!foundAny) yMin = 0.0;
}

void GetRatioMinMax(TH1* h, double xMin, double xMax, double& rMin, double& rMax) {
    rMax = -1e9; rMin = 1e9;
    bool found = false;
    if (!h) { rMin = 0.95; rMax = 1.05; return; }
    int bStart = h->FindBin(xMin);
    int bEnd = h->FindBin(xMax);
    for (int b = bStart; b <= bEnd; ++b) {
        double c = h->GetBinContent(b);
        if (c > 1e-9) { 
            if (c > rMax) rMax = c;
            if (c < rMin) rMin = c;
            found = true;
        }
    }
    if (!found) { rMin = 0.8; rMax = 1.2; }
}

TGraphAsymmErrors* HistToGraph(TH1* h) {
    if (!h) return nullptr;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(h);
    for (int i = g->GetN() - 1; i >= 0; --i) {
        if (g->GetY()[i] <= 0.001 || g->GetEYhigh()[i] > 0.5) g->RemovePoint(i);
    }
    return g;
}

vector<double> GetQuantileNodes(TH1* h, double xmin, double xmax, int n_seg) {
    vector<double> validX;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double center = h->GetBinCenter(i);
        if (center >= xmin && center <= xmax && h->GetBinContent(i) > 1e-9) validX.push_back(center);
    }
    vector<double> nodes;
    if (validX.empty()) { nodes.push_back(xmin); nodes.push_back(xmax); return nodes; }
    int nPoints = validX.size();
    if (nPoints <= n_seg) nodes = validX; 
    else {
        nodes.push_back(xmin); 
        double step = (double)nPoints / (double)n_seg;
        for (int i = 1; i < n_seg; ++i) {
            int idx = (int)(i * step);
            if (validX[idx] > nodes.back() + 0.001) nodes.push_back(validX[idx]);
        }
        if (xmax > nodes.back()) nodes.push_back(xmax);
    }
    return nodes;
}

TF1* RunAdaptiveSplineFit(TH1* hRatio, string name, double xmin, double xmax, string opts) {
    if (!hRatio) return nullptr;
    TGraphAsymmErrors* g = HistToGraph(hRatio);
    if (!g || g->GetN() < 4) return nullptr; 
    TF1* bestFit = nullptr;
    for (int n_seg = 3; n_seg <= 12; ++n_seg) {
        vector<double> nodes = GetQuantileNodes(hRatio, xmin, xmax, n_seg);
        if (nodes.size() < 2) continue;
        TF1* f = SplineFit(g, nodes.data(), nodes.size(), 0x38, opts.c_str(), (name + to_string(n_seg)).c_str(), xmin, xmax);
        if (!f) continue;
        bestFit = f;
        if (f->GetNDF() > 0 && f->GetChisquare()/f->GetNDF() < 2.5) break;
    }
    if (bestFit) { bestFit->SetLineColor(kRed); bestFit->SetLineWidth(2); }
    return bestFit;
}

void CalEff() {
    string dirIn = "/eos/user/z/zixuan/Isotope/Add/";
    string dirOut = "/eos/user/z/zixuan/Isotope/Eff/";
    string outRatioFile = dirOut + "ISSMCRatioFit2.root";
    string fileISS = dirIn + "Be_frag4_withBkg_Tune_full.root";

    gSystem->mkdir(dirOut.c_str(), true);
    gStyle->SetOptStat(0); gStyle->SetPadGridX(true); gStyle->SetPadGridY(true);

    map<int, vector<tuple<string, string, double>>> isoConfigs;
    isoConfigs[2] = {{"He4", "He4_rew_frag4_NoBkg_full.root", 1.0}};
    isoConfigs[3] = {{"Li6", "Li6_rew_frag4_NoBkg_full.root", 0.9}, {"Li7", "Li7_rew_frag4_NoBkg_full.root", 0.1}};
    isoConfigs[4] = {{"Be7", "Be7_rew_frag4_NoBkg_full.root", 0.6}, {"Be9", "Be9_rew_frag4_NoBkg_full.root", 0.3}, {"Be10", "Be10_rew_frag4_NoBkg_full.root", 0.1}};
    isoConfigs[5] = {{"B10", "B10_rew_frag4_NoBkg_full.root", 0.3},{"B11", "B11_rew_frag4_NoBkg_full.root", 0.7}};
    isoConfigs[6] = {{"C12", "C12_rew_frag4_NoBkg_full.root", 1.0}};
    isoConfigs[7] = {{"N14", "N14_rew_frag4_NoBkg_full.root", 0.5}, {"N15", "N15_rew_frag4_NoBkg_full.root", 0.5}};
    isoConfigs[8] = {{"O16", "O16_rew_frag4_NoBkg_full.root", 1.0}};

    int compColors[] = {kBlue, kMagenta, kOrange+7, kCyan+2, kGreen+3};

    TFile* fISS = TFile::Open(fileISS.c_str());
    TFile* fOutRatio = new TFile(outRatioFile.c_str(), "RECREATE");

    vector<string> trackerCuts = {"Trigger", "L1QLowLimit", "L1PickUp", "InnerTracking", "InnerTrackerQ", "UpperTOFQ", "BkgReduction"};
    vector<string> betaDets = {"TOF", "NaF", "AGL"};

    for (int z = 2; z <= 8; ++z) {
        string zSuffix = "Z" + to_string(z);
        vector<pair<string, string>> tasks;
        for (auto& c : trackerCuts) tasks.push_back({c, "Tracker"});
        for (auto& d : betaDets) tasks.push_back({"BetaRecQuality", d});

        for (auto& task : tasks) {
            string cut = task.first; string det = task.second;
            double drawXMin, drawXMax;
            if (cut == "BetaRecQuality") {
                if (det == "TOF") { drawXMin = 0.3; drawXMax = 2.0; }
                else if (det == "NaF") { drawXMin = 0.7; drawXMax = 30.0; }
                else { drawXMin = 2.0; drawXMax = 60.0; }
            } else {
                drawXMin = 2.0; drawXMax = 1000.0;
            }

            double totalIntegral = 0;
            vector<double> isoIntegrals;
            for (auto& config : isoConfigs[z]) {
                TFile* fTemp = TFile::Open((dirIn + get<1>(config)).c_str());
                TH1D* hGen = (TH1D*)fTemp->Get("MC_FLUX_H3");
                double integral = (hGen) ? hGen->Integral() : 0;
                isoIntegrals.push_back(integral);
                totalIntegral += integral;
                fTemp->Close();
            }

            TH1D *hNumTot = nullptr, *hDenTot = nullptr;
            vector<TH1D*> hEffSingles;
            vector<string> isoNames;

            for (size_t i = 0; i < isoConfigs[z].size(); ++i) {
                if (isoIntegrals[i] <= 0) continue;
                TFile* fMC = TFile::Open((dirIn + get<1>(isoConfigs[z][i])).c_str());
                TH1D* hN = (TH1D*)fMC->Get(("Eff_FLUXH1_" + cut + "_Num_" + det + "_" + zSuffix).c_str());
                TH1D* hD = (TH1D*)fMC->Get(("Eff_FLUXH1_" + cut + "_Den_" + det + "_" + zSuffix).c_str());
                if (!hN || !hD) { fMC->Close(); continue; }

                TH1D* hEffSingle = (TH1D*)hN->Clone(Form("hEff_%s_%s_%s", get<0>(isoConfigs[z][i]).c_str(), cut.c_str(), det.c_str()));
                hEffSingle->SetDirectory(0);
                if (cut == "Trigger") {
                    TH1D* hTrigS = (TH1D*)hN->Clone(); hTrigS->Add(hD);
                    hEffSingle->Divide(hN, hTrigS, 1, 1, "B"); delete hTrigS;
                } else hEffSingle->Divide(hN, hD, 1, 1, "B");
                hEffSingles.push_back(hEffSingle);
                isoNames.push_back(get<0>(isoConfigs[z][i]));

                double weight = get<2>(isoConfigs[z][i]) * (totalIntegral / isoIntegrals[i]);
                if (!hNumTot) {
                    hNumTot = (TH1D*)hN->Clone(); hNumTot->Reset(); hNumTot->SetDirectory(0);
                    hDenTot = (TH1D*)hD->Clone(); hDenTot->Reset(); hDenTot->SetDirectory(0);
                }
                hNumTot->Add(hN, weight);
                hDenTot->Add(hD, weight);
                fMC->Close();
            }

            if (!hNumTot) continue;
            TH1D* hEffMC = (TH1D*)hNumTot->Clone();
            if (cut == "Trigger") {
                TH1D *hTrigM = (TH1D*)hNumTot->Clone(); hTrigM->Add(hDenTot);
                hEffMC->Divide(hNumTot, hTrigM, 1, 1, "B"); delete hTrigM;
            } else hEffMC->Divide(hNumTot, hDenTot, 1, 1, "B");

            TH1D *hEffISS = nullptr;
            TH1D *hRatio = nullptr;

            if (cut == "BkgReduction" && (z == 3 || z == 4 || z == 5 || z == 7)) {
                int zL = (z < 6) ? 2 : 6; int zH = (z < 6) ? 6 : 8;
                double alpha = (double)(z - zL) / (zH - zL);
                auto GetRefRatio = [&](int refZ) {
                    double tInt = 0; vector<double> sInts;
                    for (auto& cfg : isoConfigs[refZ]) {
                        TFile* f = TFile::Open((dirIn + get<1>(cfg)).c_str());
                        TH1D* hG = (TH1D*)f->Get("MC_FLUX_H3");
                        double it = hG ? hG->Integral() : 0; sInts.push_back(it); tInt += it; f->Close();
                    }
                    TH1D *hNM = nullptr, *hDM = nullptr;
                    for (size_t i = 0; i < isoConfigs[refZ].size(); ++i) {
                        if (sInts[i] <= 0) continue;
                        TFile* fM = TFile::Open((dirIn + get<1>(isoConfigs[refZ][i])).c_str());
                        TH1D* n = (TH1D*)fM->Get(("Eff_FLUXH1_BkgReduction_Num_" + det + "_Z" + to_string(refZ)).c_str());
                        TH1D* d = (TH1D*)fM->Get(("Eff_FLUXH1_BkgReduction_Den_" + det + "_Z" + to_string(refZ)).c_str());
                        if(!n || !d) { fM->Close(); continue; }
                        double w = get<2>(isoConfigs[refZ][i]) * (tInt / sInts[i]);
                        if (!hNM) { hNM=(TH1D*)n->Clone(); hNM->Reset(); hNM->SetDirectory(0); hDM=(TH1D*)d->Clone(); hDM->Reset(); hDM->SetDirectory(0); }
                        hNM->Add(n, w); hDM->Add(d, w); fM->Close();
                    }
                    if (!hNM) return (TH1D*)nullptr;
                    TH1D* eM = (TH1D*)hNM->Clone(); eM->Divide(hNM, hDM, 1, 1, "B");
                    TH1D* hNI = (TH1D*)fISS->Get(("Eff_FLUXH1_BkgReduction_Num_" + det + "_Z" + to_string(refZ)).c_str());
                    TH1D* hDI = (TH1D*)fISS->Get(("Eff_FLUXH1_BkgReduction_Den_" + det + "_Z" + to_string(refZ)).c_str());
                    if (!hNI || !hDI) { delete hNM; delete hDM; delete eM; return (TH1D*)nullptr; }
                    TH1D* eI = (TH1D*)hNI->Clone(); eI->Divide(hNI, hDI, 1, 1, "B");
                    TH1D* r = (TH1D*)eI->Clone(); r->Divide(eM);
                    delete hNM; delete hDM; delete eM; delete eI; return r;
                };
                TH1D* rL = GetRefRatio(zL); TH1D* rH = GetRefRatio(zH);
                if (rL && rH) {
                    hRatio = (TH1D*)rL->Clone(Form("Ratio_Z%d_%s_%s", z, cut.c_str(), det.c_str()));
                    hRatio->Reset();
                    for(int i=1; i<=rL->GetNbinsX(); ++i) {
                        double vL = rL->GetBinContent(i), vH = rH->GetBinContent(i);
                        double eL = rL->GetBinError(i), eH = rH->GetBinError(i);
                        if(vL > 0 && vH > 0) {
                            hRatio->SetBinContent(i, vL + alpha * (vH - vL));
                            hRatio->SetBinError(i, eL + alpha * (eH - eL));
                        }
                    }
                    hEffISS = (TH1D*)hEffMC->Clone();
                    hEffISS->Multiply(hRatio);
                }
                delete rL; delete rH;
            } else {
                TH1D *hN = (TH1D*)fISS->Get(("Eff_FLUXH1_" + cut + "_Num_" + det + "_" + zSuffix).c_str());
                TH1D *hD = (TH1D*)fISS->Get(("Eff_FLUXH1_" + cut + "_Den_" + det + "_" + zSuffix).c_str());
                if (hN && hD) {
                    hEffISS = (TH1D*)hN->Clone();
                    if (cut == "Trigger") {
                        TH1D *hTrig = (TH1D*)hN->Clone(); hTrig->Add(hD);
                        hEffISS->Divide(hN, hTrig, 1, 1, "B"); delete hTrig;
                    } else hEffISS->Divide(hN, hD, 1, 1, "B");
                    hRatio = (TH1D*)hEffISS->Clone(Form("Ratio_Z%d_%s_%s", z, cut.c_str(), det.c_str()));
                    hRatio->Divide(hEffMC);
                }
            }

            if (!hEffISS || !hRatio) { delete hEffMC; continue; }

            TCanvas* c = new TCanvas("c", "", 800, 800);
            TPad *p1 = new TPad("p1","",0,0.35,1,1); p1->SetBottomMargin(0.02); p1->SetGrid(); p1->SetLogx(); p1->Draw();
            p1->SetRightMargin(0.05); p1->SetLeftMargin(0.15); p1->SetTopMargin(0.05);
            TPad *p2 = new TPad("p2","",0,0,1,0.35); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.25); p2->SetGrid(); p2->SetLogx(); p2->Draw();
            p2->SetRightMargin(0.05); p2->SetLeftMargin(0.15);

            p1->cd();
            SetStyle(hEffISS, kBlack, 20); hEffISS->GetXaxis()->SetLabelSize(0);
            SetStyle(hEffMC, kRed, 21);
            vector<TH1D*> allHists = {hEffISS, hEffMC};
            for(size_t i=0; i<hEffSingles.size(); ++i) {
                SetStyle(hEffSingles[i], compColors[i%5], 24);
                allHists.push_back(hEffSingles[i]);
            }
            double gMin, gMax; 
            GetGlobalMinMax(allHists, drawXMin, drawXMax, gMin, gMax);
            hEffISS->GetYaxis()->SetRangeUser(gMin*0.95, gMax*1.05);
            hEffISS->GetXaxis()->SetRangeUser(drawXMin, drawXMax);
            hEffISS->Draw("PE"); 
            hEffMC->Draw("PE SAME");
            for(auto hS : hEffSingles) hS->Draw("PE SAME");

            TLegend* leg = new TLegend(0.12, 0.15, 0.45, 0.45); leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
            leg->AddEntry(hEffISS, Form("ISS Z=%d", z), "lp");
            leg->AddEntry(hEffMC, "MC Mixed", "lp");
            for(size_t i=0; i<isoNames.size(); ++i) leg->AddEntry(hEffSingles[i], isoNames[i].c_str(), "lp");
            leg->Draw();

            p2->cd();
            SetStyle(hRatio, kBlack, 20); hRatio->GetYaxis()->SetTitle("ISS/MC Mix");
            double rMin, rMax; 
            GetRatioMinMax(hRatio, drawXMin, drawXMax, rMin, rMax);
            hRatio->GetYaxis()->SetRangeUser(rMin*0.95, rMax*1.05);
            hRatio->GetXaxis()->SetRangeUser(drawXMin, drawXMax);
            hRatio->GetYaxis()->SetNdivisions(505);
            hRatio->GetYaxis()->SetTitleOffset(0.7);
            hRatio->GetYaxis()->SetTitleSize(0.1); 
            hRatio->GetYaxis()->SetLabelSize(0.1); 
            hRatio->GetXaxis()->SetTitleSize(0.1);
            hRatio->GetXaxis()->SetLabelSize(0.1);
            hRatio->Draw("PE");

            double fMin = drawXMin, fMax = (cut == "InnerTracking") ? 18 : drawXMax;
            TF1* fit = RunAdaptiveSplineFit(hRatio, Form("fit_Z%d_%s_%s", z, cut.c_str(), det.c_str()), fMin, fMax, "b2e1");
            if (fit) {
                if (cut == "InnerTracking") {
                    double yVal = fit->Eval(18);
                    TF1* comp = new TF1(Form("fit_Z%d_%s_%s", z, cut.c_str(), det.c_str()), [fit, yVal](double* x, double* p){ return (x[0] < 18) ? fit->Eval(x[0]) : yVal; }, fMin, drawXMax, 0);
                    comp->SetLineColor(kRed); comp->Draw("SAME");
                    fOutRatio->cd(); comp->Write();
                } else {
                    fit->Draw("SAME"); fOutRatio->cd(); fit->Write();
                }
            }
            c->SaveAs((dirOut + "Eff_Comp_Z" + to_string(z) + "_" + cut + "_" + det + ".png").c_str());
            for(auto h : hEffSingles) delete h;
            delete c; delete hNumTot; delete hDenTot; delete hEffMC; delete hEffISS; delete hRatio;
        }
    }
    fOutRatio->Close(); fISS->Close();
}