#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLine.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <algorithm>

// --- Basic Structures ---
struct ValueWithError {
    double val = 0.0; double err = 0.0;
    ValueWithError() = default;
    ValueWithError(double v, double e) : val(v), err(e) {}
    ValueWithError operator+(const ValueWithError& o) const { return {val + o.val, std::hypot(err, o.err)}; }
    ValueWithError operator-(const ValueWithError& o) const { return {val - o.val, std::hypot(err, o.err)}; }
    ValueWithError operator*(double s) const { return {val * s, std::abs(err * s)}; }
    ValueWithError operator*(const ValueWithError& o) const { return {val * o.val, std::hypot(err * o.val, val * o.err)}; }
    ValueWithError operator/(const ValueWithError& o) const {
        if (o.val <= 0) return {0, 0};
        double q = val / o.val;
        double e = std::hypot(err / o.val, (val * o.err) / (o.val * o.val));
        return {q, e};
    }
};

struct AnalysisConfig {
    std::string source; std::string fragment; std::string fragFileID; std::string fragAbbr;
    int useMass; int sourceZ;
    std::map<std::string, double> srcAbundance;
    std::map<std::string, std::string> mcFiles;
    std::map<std::string, int> fragIsoMap;
    std::vector<std::string> targetIsotopes;
};

struct DetResults { 
    std::map<std::string, TH1D*> hists; 
};

// --- Helper: Get Min/Max for Y-axis ---
void getMinMax(TH1* h, double& minVal, double& maxVal) {
    if(!h) return;
    for(int i=1; i<=h->GetNbinsX(); ++i) {
        double c = h->GetBinCenter(i);
        if(c < 0.21 || c > 23.0) continue; 
        
        double v = h->GetBinContent(i);
        double e = 0.4*h->GetBinError(i);
        
        if(v <= 0) continue; 

        double top = v + e;
        double bottom = v - e;

        if(top > maxVal) maxVal = top;
        if(bottom > 0) {
            if(bottom < minVal) minVal = bottom;
        } else {
            if(v < minVal) minVal = v;
        }
    }
}

TH1D* GetHistD(TFile* f, TString name) {
    if (!f || f->IsZombie()) return nullptr;
    TH1D* h = (TH1D*)f->Get(name);
    if (h) h->SetDirectory(0);
    return h;
}

TH1D* CreateHist(const char* name, TH1D* ref) {
    TH1D* h = new TH1D(name, "", ref->GetNbinsX(), ref->GetXaxis()->GetXbins()->GetArray());
    h->SetDirectory(0); h->Sumw2();
    return h;
}

ValueWithError GetBinVal(TH1D* h, int bin) {
    if (!h || bin <= 0 || bin > h->GetNbinsX()) return {0, 0};
    return {h->GetBinContent(bin), h->GetBinError(bin) > 0 ? h->GetBinError(bin) : 0};
}

// --- Plotting Function: Overlap ---
void CreateOverlapPlot(DetResults& issRes, DetResults& mcRes, std::string iso, std::string equName, TString outPath, const AnalysisConfig& cfg) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    TCanvas* c1 = new TCanvas(Form("c_%s_%s", equName.c_str(), iso.c_str()), "", 900, 800);
    TPad* p1 = new TPad("p1", "Yields", 0.0, 0.5, 1.0, 1.0); 
    TPad* p2 = new TPad("p2", "Ratios", 0.0, 0.0, 1.0, 0.5); 
    
    p1->SetBottomMargin(0.04); p1->SetTopMargin(0.08); p1->SetLeftMargin(0.15); p1->SetRightMargin(0.05);
    p2->SetTopMargin(0.02); p2->SetBottomMargin(0.2); p2->SetLeftMargin(0.15); p2->SetRightMargin(0.05);
    p1->Draw(); p2->Draw();

    struct DetAttr { std::string name; int color; };
    std::vector<DetAttr> dets = { {"TOF", kRed}, {"NaF", kBlue}, {"AGL", kGreen+2} };

    // --- Pad 1: Distributions ---
    p1->cd(); p1->SetGridy();  p1->SetLogx();

    TLegend* leg = new TLegend(0.65, 0.65, 0.92, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.06);

    double yMin = 1e9, yMax = -1e9;
    std::vector<TH1D*> histsTop;

    for (auto& d : dets) {
        std::string key = iso + "_" + d.name;
        if (issRes.hists.count(key)) {
            TH1D* h = issRes.hists[key];
            h->SetLineColor(d.color); h->SetMarkerColor(d.color); h->SetMarkerStyle(20); h->SetMarkerSize(1.0);
            getMinMax(h, yMin, yMax); 
            histsTop.push_back(h);
            leg->AddEntry(h, Form("ISS %s", d.name.c_str()), "p");
        }
        if (mcRes.hists.count(key)) {
            TH1D* h = mcRes.hists[key];
            h->SetLineColor(d.color); h->SetMarkerColor(d.color); h->SetMarkerStyle(24); h->SetMarkerSize(1.0);
            getMinMax(h, yMin, yMax);
            histsTop.push_back(h);
        }
    }
    leg->AddEntry((TObject*)0, "Open Circle: MC", "");

    if(yMin >= yMax) { yMin = 0.8 * yMin; yMax = 1.2 * yMin; }
    double diff = yMax - yMin;
    double plotMin = yMin - diff * 0.6;
    plotMin = (plotMin <= 0) ? 0.0 : plotMin;
    double plotMax = yMax + diff * 0.6;

    bool first = true;
    for (auto h : histsTop) {
        h->GetXaxis()->SetRangeUser(0.28, 23.0);
        if (first) {
            TString sourceName = cfg.fragAbbr; // Be
            TString motherName = cfg.source.substr(0,1); // 获取元素首字母, 如 Boron -> B

            TString formattedIso = iso;
            if (iso == "Total") formattedIso = sourceName; 

            TString yTitle = "";
            if (equName == "Equ1") {
                yTitle = Form("N_{%s #rightarrow %s} / N_{%s}", motherName.Data(), formattedIso.Data(), motherName.Data());
            } else if (equName == "Equ2") {
                yTitle = Form("N_{%s #rightarrow %s} / N_{%s, survival}", motherName.Data(), formattedIso.Data(), motherName.Data());
            } else if (equName == "Equ3") {
                yTitle = Form("N_{%s #rightarrow %s} / N_{%s, survival}", motherName.Data(), formattedIso.Data(), motherName.Data());
            } else if (equName == "Equ4") {
                yTitle = Form("N_{%s #rightarrow %s} / N_{%s #rightarrow AllFrag}", motherName.Data(), formattedIso.Data(), motherName.Data());
            }
            h->GetYaxis()->SetTitle(yTitle);
            h->GetYaxis()->SetTitleSize(0.08); h->GetYaxis()->SetTitleOffset(0.9);
            h->GetYaxis()->SetLabelSize(0.08);; h->GetYaxis()->SetNdivisions(505);
            h->GetYaxis()->SetRangeUser(plotMin, plotMax);
            h->GetXaxis()->SetLabelSize(0); 
            h->Draw("PZ");
            first = false;
        } else {
            h->Draw("PZ SAME");
        }
    }
    leg->Draw();

    // --- Pad 2: Ratios ---
    p2->cd(); p2->SetGridy(); p2->SetLogx();

    std::vector<TH1D*> histsBot;
    double rMin = 1e9, rMax = -1e9;

    for (auto& d : dets) {
        std::string key = iso + "_" + d.name;
        if (issRes.hists.count(key) && mcRes.hists.count(key)) {
            TH1D* h_iss = issRes.hists[key];
            TH1D* h_mc = mcRes.hists[key];
            TH1D* h_ratio = (TH1D*)h_iss->Clone(Form("ratio_%s", key.c_str()));
            h_ratio->Divide(h_mc);
            h_ratio->SetLineColor(d.color); h_ratio->SetMarkerColor(d.color); h_ratio->SetMarkerStyle(20);
            getMinMax(h_ratio, rMin, rMax);
            histsBot.push_back(h_ratio);
        }
    }

    if(rMin >= rMax) { rMin = 0.5; rMax = 1.5; }
    double rDiff = rMax - rMin;
    double rPlotMin = rMin - rDiff * 0.25; 
    rPlotMin = (rPlotMin < 0) ? 0 : rPlotMin;
    double rPlotMax = rMax + rDiff * 0.25;

    first = true;
    for (auto h : histsBot) {
        h->GetXaxis()->SetRangeUser(0.28, 23.0);
        if (first) {
            h->GetYaxis()->SetTitle("ISS / MC");
            h->GetXaxis()->SetTitle("E_{k}/n [GeV]");
            h->GetYaxis()->SetTitleSize(0.08); h->GetYaxis()->SetTitleOffset(0.9);
            h->GetYaxis()->SetLabelSize(0.08); h->GetYaxis()->SetNdivisions(505);
            h->GetXaxis()->SetTitleSize(0.08); h->GetXaxis()->SetTitleOffset(1.1);
            h->GetXaxis()->SetLabelSize(0.08);
            h->GetYaxis()->SetRangeUser(rPlotMin, rPlotMax);
            h->Draw("PZ");
            first = false;
        } else {
            h->Draw("PZ SAME");
        }
    }
    
    TLine* l = new TLine(0.28, 1.0, 23.0, 1.0);
    l->SetLineStyle(2); l->SetLineColor(kBlack); l->Draw();
    c1->SaveAs(outPath + ".png");
    delete c1;
}

// --- MC Calculation Logic ---
DetResults CalculateMCRatio(const AnalysisConfig& cfg, TH1D* refHist, TString chain, std::string denTag, std::string numTag) {
    DetResults res;
    std::vector<std::string> dets = {"TOF", "NaF", "AGL"};
    std::string mcPath = "/eos/user/z/zixuan/Isotope/Add/";
    
    std::map<std::string, TH1D*> sum_den;
    std::map<std::string, TH1D*> sum_num_tot; 
    std::map<std::string, std::map<std::string, TH1D*>> sum_num; 

    // Initialize
    for(auto d : dets) {
        sum_den[d] = CreateHist(("mc_den_"+d).c_str(), refHist);
        sum_num_tot[d] = CreateHist(("mc_num_tot_"+d).c_str(), refHist);
        for (auto const& [iso, m] : cfg.fragIsoMap) {
            sum_num[iso][d] = CreateHist(("mc_num_"+iso+"_"+d).c_str(), refHist);
        }
    }

    // Load Weights
    double n_gen_all_sources = 0;
    std::map<std::string, double> n_gen_individual;
    for (auto const& [srcIso, srcFileBase] : cfg.mcFiles) {
        TString fullPath = mcPath + srcFileBase.c_str() + "_rew_frag4_NoBkg_full.root";
        TFile* f = TFile::Open(fullPath);
        if(!f || f->IsZombie()) { std::cout << "MC File Missing: " << fullPath << std::endl; continue; }
        TH1D* h_gen = GetHistD(f, "MC_FLUX_H3");
        if(h_gen) {
            double n_gen = h_gen->Integral();
            n_gen_individual[srcIso] = n_gen; n_gen_all_sources += n_gen;
        } else { std::cout << "MC_FLUX_H3 missing in " << srcFileBase << std::endl; }
        f->Close(); delete f;
    }

    std::cout << "\n>>>>>> MC Processing (" << cfg.source << ") n_gen_tot=" << n_gen_all_sources << " <<<<<<" << std::endl;
    
    for (auto const& [srcIso, srcFileBase] : cfg.mcFiles) {
        TString fullPath = mcPath + srcFileBase.c_str() + "_rew_frag4_NoBkg_full.root";
        TFile* f = TFile::Open(fullPath);
        if(!f || f->IsZombie()) continue;
        
        double global_w = (n_gen_individual[srcIso] > 0) ? (cfg.srcAbundance.at(srcIso) * n_gen_all_sources / n_gen_individual[srcIso]) : 0;
        
        for(auto d : dets) {
            double eMin = 0, eMax = 1000;
            if(d == "TOF") { eMin = 0.21; eMax = 1.6; }
            else if(d == "NaF") { eMin = 0.61; eMax = 6.1; }
            else if(d == "AGL") { eMin = 2.7; eMax = 30.0; }

            // Denominator Handling
            ValueWithError vD(0,0);
            if (denTag == "H1a_minus_H1b") { // Equ4 Special
                TString dn1 = TString::Format("%s_BKG_%s_%s_%s", chain.Data(), "H1a", cfg.source.c_str(), d.c_str());
                TString dn2 = TString::Format("%s_BKG_%s_%s_%s", chain.Data(), "H1b", cfg.source.c_str(), d.c_str());
                TH1D* h1 = GetHistD(f, dn1);
                TH1D* h2 = GetHistD(f, dn2);
                if (h1 && h2) {
                    // Fill loop for subtraction
                     h1->Rebin(2); h2->Rebin(2);
                     for(int b=1; b<=sum_den[d]->GetNbinsX(); ++b) {
                         double cent = sum_den[d]->GetBinCenter(b);
                         if(cent >= eMin && cent <= eMax) {
                             ValueWithError v = (GetBinVal(h1, h1->FindBin(cent)) - GetBinVal(h2, h2->FindBin(cent))) * global_w;
                             sum_den[d]->SetBinContent(b, sum_den[d]->GetBinContent(b) + v.val);
                             sum_den[d]->SetBinError(b, std::hypot(sum_den[d]->GetBinError(b), v.err));
                         }
                     }
                }
            } else { // Standard Case
                TString dn = TString::Format("%s_BKG_%s_%s_%s", chain.Data(), denTag.c_str(), cfg.source.c_str(), d.c_str());
                TH1D* hD = GetHistD(f, dn);
                if(hD) {
                    hD->Rebin(2);
                    for(int b=1; b<=sum_den[d]->GetNbinsX(); ++b) {
                        double cent = sum_den[d]->GetBinCenter(b);
                        if(cent >= eMin && cent <= eMax) {
                             ValueWithError v = GetBinVal(hD, hD->FindBin(cent)) * global_w;
                             sum_den[d]->SetBinContent(b, sum_den[d]->GetBinContent(b) + v.val);
                             sum_den[d]->SetBinError(b, std::hypot(sum_den[d]->GetBinError(b), v.err));
                        }
                    }
                } else { std::cout << " [Warn] MC Denom Missing: " << dn << std::endl; }
            }

            // Numerator Isotopes
            for (auto const& [iso, mass] : cfg.fragIsoMap) {
                TString nn = TString::Format("%s_BKG_%s_%s_Z4_Mass%d", chain.Data(), numTag.c_str(), d.c_str(), mass);
                TH1D* hN = GetHistD(f, nn);
                if(hN) {
                    hN->Rebin(2);
                    for(int b=1; b<=sum_num[iso][d]->GetNbinsX(); ++b) {
                         double cent = sum_num[iso][d]->GetBinCenter(b);
                         if(cent >= eMin && cent <= eMax) {
                             ValueWithError v = GetBinVal(hN, hN->FindBin(cent)) * global_w;
                             sum_num[iso][d]->SetBinContent(b, sum_num[iso][d]->GetBinContent(b) + v.val);
                             sum_num[iso][d]->SetBinError(b, std::hypot(sum_num[iso][d]->GetBinError(b), v.err));
                         }
                    }
                } else { std::cout << " [Warn] MC Num Missing: " << nn << std::endl; }
            }
        }
        f->Close(); delete f;
    }

    // Final Calculation
    for(auto d : dets) {
        sum_num_tot[d]->Reset();
        for (auto const& [iso, m] : cfg.fragIsoMap) sum_num_tot[d]->Add(sum_num[iso][d]);

        for(auto const& [iso, mass] : cfg.fragIsoMap) {
            TH1D* r = CreateHist(("mc_r_"+iso+"_"+d).c_str(), refHist);
            r->Divide(sum_num[iso][d], sum_den[d]);
            res.hists[iso + "_" + d] = r;
        }
        TH1D* rt = CreateHist(("mc_r_Total_"+d).c_str(), refHist);
        rt->Divide(sum_num_tot[d], sum_den[d]);
        res.hists["Total_" + d] = rt;
    }
    return res;
}

void RunAnalysis(TString chain, const AnalysisConfig& cfg) {
    std::cout << "\n>>>>>> Processing Source: " << cfg.source.c_str() << " <<<<<<" << std::endl;
    
    TString pMeas = TString::Format("/eos/user/z/zixuan/Isotope/ChargeTemp/NoBkg_PureQFit_InnerBeryllium_%s.root", chain.Data());
    TString pSource = TString::Format("/eos/user/z/zixuan/Isotope/ChargeTemp/NoBkg_PureQFit_Inner%s_%s.root", cfg.source.c_str(), chain.Data());
    
    TFile* fMeas = TFile::Open(pMeas);
    TFile* fSource = TFile::Open(pSource);
    if(!fMeas || !fSource) { std::cout << "Error opening files!" << std::endl; return; }

    TFile* fOut = new TFile(TString::Format("/eos/user/z/zixuan/Isotope/BkgValid/%s_to_%s_%s_Validation.root", cfg.source.c_str(), cfg.fragment.c_str(), chain.Data()), "RECREATE");

    struct EquInfo { std::string name; std::string l1Sig; std::string mcDen; std::string mcNum; std::string massSuffix; std::string denChain; std::string sourceSigType; };
    std::vector<EquInfo> equs = {
        {"Equ1", "L1Sig_Any",       "H1a", "H2a", "_NoBkg", "UnbiasedL1Inner", "L1Sig_Any"},
        {"Equ2", "L1Sig_PassLoose", "H1b", "H2a", "_NoBkg", "UnbiasedL1Inner", "L1Sig_PassLoose"},
        {"Equ3", "L1Sig_Pass",      "H1b", "H2a2", "_NoBkg", "UnbiasedL1Inner", "L1Sig_PassLoose"},
        {"Equ4", "L1Sig_PassLoose", "H1a_minus_H1b", "H2a", "_NoBkg", "UnbiasedL1Inner", "Special"} // New Equ4
    };

    TH1D* ref_raw = GetHistD(fMeas, "h_yield_Beryllium_in_Beryllium_window_TOF_L1Sig_Any_UnbiasedL1Inner");
    if(!ref_raw) return;
    TH1D* ref = (TH1D*)ref_raw->Clone("ref_binned"); ref->Rebin(2);

    TFile* fMF_Contam = TFile::Open("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_Be_UnbiasedL1Inner_H1_UseMass7_NoBkg.root");
    if(!fMF_Contam) { std::cout << "Contamination MassFit file missing!" << std::endl; return; }

    for(const auto& eq : equs) {
        std::cout << "  - Calculation for " << eq.name << std::endl;
        
        DetResults mcRes = CalculateMCRatio(cfg, ref, chain, eq.mcDen, eq.mcNum);
        
        TString pMFN = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H1_UseMass%d_FragFrom%s%s.root", cfg.fragAbbr.c_str(), chain.Data(), cfg.useMass, cfg.source.c_str(), eq.massSuffix.c_str());
        TFile *fMF = TFile::Open(pMFN);
        if(!fMF) { std::cout << "Signal Mass Fit file missing: " << pMFN << std::endl; continue; }

        DetResults issRes;
        std::vector<std::string> detNames = {"TOF", "NaF", "AGL"};

        for(auto iso : cfg.targetIsotopes) {
            for(auto d : detNames) issRes.hists[iso + "_" + d] = CreateHist(("h_iss_"+eq.name+"_"+iso+"_"+d).c_str(), ref);
        }

        for (int i = 1; i <= ref->GetNbinsX(); ++i) {
            double cent = ref->GetBinCenter(i);
            
            for(auto det : detNames) {
                double dMin=0, dMax=0;
                if(det=="TOF") { dMin=0.21; dMax=1.6; }
                else if(det=="NaF") { dMin=0.61; dMax=6.1; }
                else if(det=="AGL") { dMin=2.7; dMax=30.0; }
                if(cent < dMin || cent > dMax) continue; 

                // A. vRaw
                TString measSig = (eq.name == "Equ1") ? "L1Sig_PassLoose" : eq.l1Sig;
                TString nRawName = TString::Format("h_raw_count_in_%s_window_%s_%s_%s", cfg.source.c_str(), det.c_str(), measSig.Data(), chain.Data());
                TH1D* hRaw = GetHistD(fMeas, nRawName);
                ValueWithError vRaw(0,0);
                if(hRaw) { hRaw->Rebin(2); vRaw = GetBinVal(hRaw, hRaw->FindBin(cent)); }

                // B. vSigTotal (Num)
                TString nSigName = TString::Format("h_yield_%s_in_%s_window_%s_%s_%s", cfg.source.c_str(), cfg.source.c_str(), det.c_str(), measSig.Data(), chain.Data());
                TH1D* hSigTotal = GetHistD(fMeas, nSigName);
                ValueWithError vSigTotal(0,0);
                if(hSigTotal) { hSigTotal->Rebin(2); vSigTotal = GetBinVal(hSigTotal, hSigTotal->FindBin(cent)); }

                // C. vContam
                TString nCntName = TString::Format("h_yield_Beryllium_in_%s_window_%s_%s_%s", cfg.source.c_str(), det.c_str(), eq.l1Sig.c_str(), chain.Data());
                TH1D* hCnt = GetHistD(fMeas, nCntName);
                ValueWithError vContam(0,0);
                if(hCnt) { hCnt->Rebin(2); vContam = GetBinVal(hCnt, hCnt->FindBin(cent)); }

                // D. vDen (Denominator)
                ValueWithError vDen(0,0);
                if (eq.name == "Equ4") { // Equ4 Special Denom: Any - PassLoose
                    TString nDen1 = TString::Format("h_yield_%s_in_%s_window_%s_%s_%s", cfg.source.c_str(), cfg.source.c_str(), det.c_str(), "L1Sig_Any", eq.denChain.c_str());
                    TString nDen2 = TString::Format("h_yield_%s_in_%s_window_%s_%s_%s", cfg.source.c_str(), cfg.source.c_str(), det.c_str(), "L1Sig_PassLoose", eq.denChain.c_str());
                    TH1D* hDen1 = GetHistD(fMeas, nDen1);
                    TH1D* hDen2 = GetHistD(fSource, nDen2);
                    if(hDen1) hDen1->Rebin(2); if(hDen2) hDen2->Rebin(2);
                    cout<< "Ek"<<cent<<" Det:"<<det<<" Den1:"<< (hDen1 ? GetBinVal(hDen1, hDen1->FindBin(cent)).val : 0) <<" Den2:"<< (hDen2 ? GetBinVal(hDen2, hDen2->FindBin(cent)).val : 0) << endl;
                    vDen = (hDen1 ? GetBinVal(hDen1, hDen1->FindBin(cent)) : ValueWithError(0,0)) - 
                           (hDen2 ? GetBinVal(hDen2, hDen2->FindBin(cent)) : ValueWithError(0,0));
                } else { // Standard
                    TString nDenName = TString::Format("h_yield_%s_in_%s_window_%s_%s_%s", cfg.source.c_str(), cfg.source.c_str(), det.c_str(), eq.sourceSigType.c_str(), eq.denChain.c_str());
                    TH1D* hDen = GetHistD(eq.name == "Equ1" ? fMeas : fSource, nDenName);
                    if(hDen) hDen->Rebin(2); 
                    vDen = GetBinVal(hDen, hDen ? hDen->FindBin(cent) : 0);
                }

                if(vDen.val <= 1e-4) continue;

                // E. Mass Fractions
                std::map<std::string, ValueWithError> F_Meas, F_Contam;
                for(auto fitIso : {"Be7", "Be9"}) {
                    TString kF = TString::Format("h_best_%s_frac_%s", fitIso, det.c_str());
                    TH1D* hf1 = GetHistD(fMF, kF);
                    F_Meas[fitIso] = hf1 ? GetBinVal(hf1, hf1->FindBin(cent)) : ValueWithError(0,0);
                    TH1D* hf2 = GetHistD(fMF_Contam, kF);
                    F_Contam[fitIso] = hf2 ? GetBinVal(hf2, hf2->FindBin(cent)) : ValueWithError(0,0);
                }
                F_Meas["Be10"] = ValueWithError(1.0, 0.0) - F_Meas["Be7"] - F_Meas["Be9"];
                F_Contam["Be10"] = ValueWithError(1.0, 0.0) - F_Contam["Be7"] - F_Contam["Be9"];

                // F. Fill Histograms (Only if value != 0)
                std::cout << Form("[%s %s Bin%d E=%.2f] Raw:%.1f+/-%.1f SigTot:%.1f+/-%.1f Cont:%.1f+/-%.1f Den:%.1f+/-%.1f", 
                                  eq.name.c_str(), det.c_str(), i, cent, vRaw.val, vRaw.err, vSigTotal.val, vSigTotal.err, vContam.val, vContam.err, vDen.val, vDen.err) << std::endl;

                // 1. Total
                ValueWithError vResTot = vSigTotal/vDen;
                if (std::abs(vResTot.val) > 1e-12) {
                    TH1D* hTot = issRes.hists["Total_" + det];
                    if(hTot) { hTot->SetBinContent(i, vResTot.val); hTot->SetBinError(i, vResTot.err); }
                }

                // 2. Isotopes
                for (auto iso : cfg.targetIsotopes) {
                    if (iso == "Total") continue;
                    // --- 原有误差计算保留并注释 ---
                    ValueWithError vTrue = (vRaw * F_Meas[iso]) - (vContam * F_Contam[iso]);
                    ValueWithError vResIso = vTrue/vDen; 
                    /*
                    // 1. 计算信号项 (Signal Component) 及其统计误差
                    double nSigVal = vRaw.val * F_Meas[iso].val;
                    double nSigErr = (nSigVal > 0) ? std::sqrt(nSigVal) : 0.0;
                    ValueWithError vSigComp(nSigVal, nSigErr);

                    // 2. 计算污染项 (Background Component) 及其统计误差
                    double nBkgVal = vContam.val * F_Contam[iso].val;
                    double nBkgErr = (nBkgVal > 0) ? std::sqrt(nBkgVal) : 0.0;
                    ValueWithError vBkgComp(nBkgVal, nBkgErr);

                    // 3. 两项相减，ValueWithError 内部会自动执行 std::hypot(nSigErr, nBkgErr)
                    ValueWithError vTrue = vSigComp - vBkgComp;

                    // 4. 计算最终比值：vTrue / vDen (分母误差按原逻辑传递)
                    ValueWithError vResIso = vTrue / vDen;
                    */

                    if (std::abs(vResIso.val) > 1e-12) {
                        TH1D* hDest = issRes.hists[iso + "_" + det];
                        if(hDest) { hDest->SetBinContent(i, vResIso.val); hDest->SetBinError(i, vResIso.err); }
                    }
                }
            }
        }

        fOut->mkdir(eq.name.c_str()); fOut->cd(eq.name.c_str());
        for(auto const& [key, h] : issRes.hists) { h->Write(); }
        for(auto const& [key, h] : mcRes.hists) { h->Write(); }

        for(auto iso : cfg.targetIsotopes) {
            TString pBase = "/eos/user/z/zixuan/Isotope/BkgValid/Plots/" + TString::Format("%s_to_%s_%s_%s_%s_NoBkg", cfg.source.c_str(), cfg.fragment.c_str(), chain.Data(), eq.name.c_str(), iso.c_str());
            CreateOverlapPlot(issRes, mcRes, iso, eq.name, pBase, cfg);
        }
        fMF->Close();
    }
    fMF_Contam->Close();
    fMeas->Close(); fSource->Close(); fOut->Close();
}

void CalFrag_new() {
    gROOT->SetBatch(kTRUE);
    
    // Boron -> Be
    AnalysisConfig cB; cB.source="Boron"; cB.fragment="Beryllium"; cB.fragAbbr="Be"; cB.useMass=7; cB.sourceZ=5;
    cB.srcAbundance={{"B10",0.3},{"B11",0.7}}; cB.mcFiles={{"B10","B10"},{"B11","B11"}};
    cB.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cB.targetIsotopes={"Be7","Be9","Be10","Total"};
    RunAnalysis("UnbiasedL1Inner", cB);

    // Carbon -> Be
    AnalysisConfig cC; cC.source="Carbon"; cC.fragment="Beryllium"; cC.fragAbbr="Be"; cC.useMass=7; cC.sourceZ=6;
    cC.srcAbundance={{"C12",1.0}}; cC.mcFiles={{"C12","C12"}};
    cC.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cC.targetIsotopes={"Be7","Be9","Be10","Total"};
    RunAnalysis("UnbiasedL1Inner", cC);

    // Nitrogen -> Be
    AnalysisConfig cN; cN.source="Nitrogen"; cN.fragment="Beryllium"; cN.fragAbbr="Be"; cN.useMass=7; cN.sourceZ=7;
    cN.srcAbundance={{"N14",0.5},{"N15",0.5}}; cN.mcFiles={{"N14","N14"},{"N15","N15"}};
    cN.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cN.targetIsotopes={"Be7","Be9","Be10","Total"};
    RunAnalysis("UnbiasedL1Inner", cN);

    // Oxygen -> Be
    AnalysisConfig cO; cO.source="Oxygen"; cO.fragment="Beryllium"; cO.fragAbbr="Be"; cO.useMass=7; cO.sourceZ=8;
    cO.srcAbundance={{"O16",1.0}}; cO.mcFiles={{"O16","O16"}};
    cO.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cO.targetIsotopes={"Be7","Be9","Be10","Total"};
    RunAnalysis("UnbiasedL1Inner", cO);
}