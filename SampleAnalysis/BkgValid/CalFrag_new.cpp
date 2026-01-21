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
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <iomanip>
#include <algorithm>

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

struct DetResults { std::map<std::string, TH1D*> combine; };

// --- 绘图辅助函数：自动拟合逻辑 ---
void DrawFit(TH1* h, TLegend* leg, double min, double max, int color, TString label) {
    bool hasData = false;
    for(int b=h->GetXaxis()->FindBin(min); b<=h->GetXaxis()->FindBin(max); ++b) {
        if(h->GetBinContent(b) > 0) hasData = true;
    }
    if(!hasData) return;

    TF1* f0 = new TF1(Form("fit0_%s_%f", h->GetName(), min), "pol0", min, max);
    f0->SetLineColor(color); f0->SetLineWidth(2);
    TFitResultPtr r0 = h->Fit(f0, "SQN+"); 

    bool useConst = false;
    if (r0->IsValid()) {
        double ndf = r0->Ndf();
        double chi2 = r0->Chi2();
        if (ndf > 0 && (chi2 / ndf) < 2.0) useConst = true;
        else if (ndf == 0 && chi2 < 1.0) useConst = true;
    }

    if (useConst) {
        f0->Draw("SAME");
        TString txt = Form("%s: %.2f #pm %.2f (Const)", label.Data(), f0->GetParameter(0), f0->GetParError(0));
        leg->AddEntry(f0, txt, "l");
    } else {
        delete f0;
        TF1* f1 = new TF1(Form("fit1_%s_%f", h->GetName(), min), "pol1", min, max);
        f1->SetLineColor(color); f1->SetLineWidth(2);
        TFitResultPtr r1 = h->Fit(f1, "SRQ+"); 
        if (r1->IsValid()) {
            f1->Draw("SAME");
            TString txt = Form("%s: %.2f + %.2f*E", label.Data(), f1->GetParameter(0), f1->GetParameter(1));
            leg->AddEntry(f1, txt, "l");
        }
    }
}

// --- 绘图函数 ---
void CreateComparisonPlots(TH1D* h_iss, TH1D* h_mc, TString title, TString outPath, int detMode) {
    if (!h_iss) return;
    auto getMinMax = [](TH1* h, double& minVal, double& maxVal) {
        if(!h) return;
        for(int i=1; i<=h->GetNbinsX(); ++i) {
            double c = h->GetBinCenter(i);
            if(c < 0.4 || c > 21.5) continue;
            double v = h->GetBinContent(i);
            if(v <= 0) continue;
            if(v < minVal) minVal = v;
            if(v > maxVal) maxVal = v;
        }
    };

    TCanvas* c1 = new TCanvas(Form("c1_%p", (void*)h_iss), title, 800, 600);
    c1->SetGrid();
    h_iss->SetMarkerStyle(20); h_iss->SetLineColor(kBlack);
    h_iss->SetTitle(title); h_iss->GetXaxis()->SetTitle("E_{k}/n [GeV]"); h_iss->GetYaxis()->SetTitle("Ratio");
    
    double yMin = 1e9, yMax = -1e9;
    getMinMax(h_iss, yMin, yMax);
    if(h_mc) getMinMax(h_mc, yMin, yMax);
    if(yMin >= yMax) { yMin = 0.8 * yMin; yMax = 1.2 * yMin; }
    double diff = yMax - yMin;
    h_iss->GetYaxis()->SetRangeUser(yMin - diff * 0.3, yMax + diff * 0.3);
    
    double xMin = 0.4, xMax = 21.5;
    if (detMode == 1) { xMin = 0.3; xMax = 1.5; }
    else if (detMode == 2) { xMin = 0.8; xMax = 6.1; }
    else if (detMode == 3) { xMin = 2.5; xMax = 21.5; }
    h_iss->GetXaxis()->SetRangeUser(xMin, xMax);
    h_iss->Draw("PZ");
    
    TLegend* leg = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg->AddEntry(h_iss, "ISS Data", "ep");
    if (h_mc) {
        h_mc->SetMarkerStyle(24); h_mc->SetMarkerColor(kRed); h_mc->SetLineColor(kRed);
        h_mc->Draw("PZ SAME"); leg->AddEntry(h_mc, "MC Simulation", "ep");
    }
    leg->Draw();
    c1->SaveAs(outPath + "_comp.png"); delete c1;

    if (h_mc) {
        TCanvas* c2 = new TCanvas(Form("c2_%p", (void*)h_iss), title + " Ratio", 800, 600);
        c2->SetGrid();
        TH1D* h_ratio = (TH1D*)h_iss->Clone(Form("ratio_%p", (void*)h_iss));
        h_ratio->Divide(h_mc);
        h_ratio->SetTitle(title + " (ISS/MC)"); h_ratio->GetYaxis()->SetTitle("ISS / MC");
        h_ratio->GetXaxis()->SetRangeUser(xMin, xMax); 
        double rMin = 1e9, rMax = -1e9;
        getMinMax(h_ratio, rMin, rMax);
        if(rMin >= rMax) { rMin = 0.5; rMax = 1.5; }
        double rDiff = rMax - rMin;
        h_ratio->GetYaxis()->SetRangeUser(rMin - rDiff * 1.0, rMax + rDiff * 1.0);
        h_ratio->Draw("PZ");
        TLegend* fitLeg = new TLegend(0.15, 0.78, 0.85, 0.88);
        fitLeg->SetBorderSize(0); fitLeg->SetFillStyle(0); fitLeg->SetTextSize(0.035);
        if (detMode == 0) {
            DrawFit(h_ratio, fitLeg, 0.3, 7.0, kMagenta, "Low"); 
            DrawFit(h_ratio, fitLeg, 7.0, 21.5, kOrange+1, "High"); 
        } else if (detMode == 1) DrawFit(h_ratio, fitLeg, 0.3, 1.5, kMagenta, "TOF");
        else if (detMode == 2) DrawFit(h_ratio, fitLeg, 0.8, 6.1, kMagenta, "NaF");
        else if (detMode == 3) {
            DrawFit(h_ratio, fitLeg, 2.5, 7.0, kMagenta, "AGL-Low");
            DrawFit(h_ratio, fitLeg, 7.0, 21.5, kOrange+1, "AGL-High");
        }
        fitLeg->Draw(); c2->SaveAs(outPath + "_ratio.png"); delete c2; delete h_ratio;
    }
}

TH1D* GetHistD(TFile* f, TString name) {
    if (!f || f->IsZombie()) return nullptr;
    TH1D* h = (TH1D*)f->Get(name);
    if (!h) return nullptr;
    h->SetDirectory(0);
    return h;
}

TH1D* CreateHist(const char* name, TH1D* ref) {
    TH1D* h = new TH1D(name, "", ref->GetNbinsX(), ref->GetXaxis()->GetXbins()->GetArray());
    h->SetDirectory(0); h->Sumw2();
    return h;
}

ValueWithError GetBinVal(TH1D* h, int bin) {
    if (!h || bin <= 0 || bin > h->GetNbinsX()) return {0, 0};
    return {h->GetBinContent(bin), h->GetBinError(bin)};
}

DetResults CalculateMCRatio(const AnalysisConfig& cfg, TH1D* refHist, TString chain, std::string denTag, std::string numTag) {
    DetResults res;
    TString mcPath = "/eos/user/z/zixuan/Isotope/Add/";
    std::vector<std::string> dets = {"TOF", "NaF", "AGL"};
    TH1D* sum_den = CreateHist("mc_sum_den", refHist);
    TH1D* sum_num_tot = CreateHist("mc_sum_num_tot", refHist);
    std::map<std::string, TH1D*> sum_num;
    for (auto const& [iso, m] : cfg.fragIsoMap) sum_num[iso] = CreateHist(("mc_sum_num_" + iso).c_str(), refHist);

    struct BinDebugInfo {
        double cent; std::string det;
        std::map<std::string, ValueWithError> srcDenRaw;
        std::map<std::string, std::map<std::string, ValueWithError>> srcNumRaw;
        std::map<std::string, double> srcWeight;
    };
    std::map<int, BinDebugInfo> binLog;

    double n_gen_all_sources = 0;
    std::map<std::string, double> n_gen_individual;
    for (auto const& [srcIso, srcFileBase] : cfg.mcFiles) {
        TString fullPath = mcPath + srcFileBase.c_str() + "_rew_frag4_withBkg_bkgest_eff.root";
        TFile* f = TFile::Open(fullPath);
        if(!f || f->IsZombie()) continue;
        TH1D* h_gen = GetHistD(f, "MC_FLUX_H3");
        if(h_gen) {
            double n_gen = h_gen->Integral();
            n_gen_individual[srcIso] = n_gen; n_gen_all_sources += n_gen;
        }
        f->Close(); delete f;
    }

    std::cout << "\n>>>>>> MC MIXING DEBUG (" << cfg.source << ") <<<<<<" << std::endl;
    for (auto const& [srcIso, srcFileBase] : cfg.mcFiles) {
        TString fullPath = mcPath + srcFileBase.c_str() + "_rew_frag4_withBkg_bkgest_eff.root";
        TFile* f = TFile::Open(fullPath);
        if(!f || f->IsZombie()) continue;
        double n_this = n_gen_individual[srcIso];
        double abund = cfg.srcAbundance.at(srcIso);
        double global_w = (n_this > 0) ? (abund * n_gen_all_sources / n_this) : 0;
        for(auto d : dets) {
            TString dn = TString::Format("%s_BKG_%s_%s_%s", chain.Data(), denTag.c_str(), cfg.source.c_str(), d.c_str());
            TH1D* h = GetHistD(f, dn);
            if(h) {
                h->Rebin(2);
                for(int b=1; b<=sum_den->GetNbinsX(); ++b) {
                    double cent = sum_den->GetBinCenter(b);
                    std::string curDet = (cent < 1.2) ? "TOF" : (cent < 3.5) ? "NaF" : "AGL";
                    if(curDet == d) {
                        ValueWithError vRaw = GetBinVal(h, h->FindBin(cent));
                        ValueWithError vScaled = vRaw * global_w;
                        sum_den->SetBinContent(b, sum_den->GetBinContent(b) + vScaled.val);
                        sum_den->SetBinError(b, std::hypot(sum_den->GetBinError(b), vScaled.err));
                        binLog[b].cent = cent; binLog[b].det = curDet;
                        binLog[b].srcDenRaw[srcIso] = vRaw; binLog[b].srcWeight[srcIso] = global_w;
                    }
                }
            }
        }
        for (auto const& [iso, mass] : cfg.fragIsoMap) {
            for(auto d : dets) {
                TString nn = TString::Format("%s_BKG_%s_%s_Z4_Mass%d", chain.Data(), numTag.c_str(), d.c_str(), mass);
                TH1D* h = GetHistD(f, nn);
                if(h) {
                    h->Rebin(2);
                    for(int b=1; b<=sum_num[iso]->GetNbinsX(); ++b) {
                        double cent = sum_num[iso]->GetBinCenter(b);
                        std::string curDet = (cent < 1.2) ? "TOF" : (cent < 3.5) ? "NaF" : "AGL";
                        if(curDet == d) {
                            ValueWithError vRaw = GetBinVal(h, h->FindBin(cent));
                            ValueWithError vScaled = vRaw * global_w;
                            sum_num[iso]->SetBinContent(b, sum_num[iso]->GetBinContent(b) + vScaled.val);
                            sum_num[iso]->SetBinError(b, std::hypot(sum_num[iso]->GetBinError(b), vScaled.err));
                            sum_num_tot->SetBinContent(b, sum_num_tot->GetBinContent(b) + vScaled.val);
                            sum_num_tot->SetBinError(b, std::hypot(sum_num_tot->GetBinError(b), vScaled.err));
                            binLog[b].srcNumRaw[srcIso][iso] = vRaw;
                        }
                    }
                }
            }
        }
        f->Close(); delete f;
    }
    for(auto const& [iso, h] : sum_num) {
        TH1D* r = CreateHist(("mc_ratio_" + iso).c_str(), refHist); r->Divide(h, sum_den); res.combine[iso] = r;
    }
    TH1D* rt = CreateHist("mc_ratio_Total", refHist); rt->Divide(sum_num_tot, sum_den); res.combine["Total"] = rt;
    return res;
}

void RunAnalysis(TString chain, const AnalysisConfig& cfg) {
    std::cout << "\n>>>>>> Processing Source: " << cfg.source.c_str() << " <<<<<<" << std::endl;
    TString pFIn = "/eos/user/z/zixuan/Isotope/ChargeTemp/FitResults_L1Inner.root";
    TString pL1N = TString::Format("/eos/user/z/zixuan/Isotope/ChargeTemp/withBkg_PureQFit_%sTo%s_%s.root", cfg.source.c_str(), cfg.fragment.c_str(), chain.Data());
    TString pL1F = TString::Format("/eos/user/z/zixuan/Isotope/ChargeTemp/withBkg_PureQFit_%sTo%s_%s.frag.root", cfg.source.c_str(), cfg.fragment.c_str(), chain.Data());
    
    TFile* fInner = TFile::Open(pFIn);
    TFile* fL1_noFrag = TFile::Open(pL1N);
    TFile* fL1_Frag = TFile::Open(pL1F);
    if(!fL1_Frag || !fL1_noFrag) return;

    TFile* fOut = new TFile(TString::Format("/eos/user/z/zixuan/Isotope/BkgValid/%s_to_%s_%s_Validation.root", cfg.source.c_str(), cfg.fragment.c_str(), chain.Data()), "RECREATE");

    struct EquInfo { std::string name; std::string l1Sig; std::string mcDen; std::string mcNum; std::string massSuffix; std::string denChain; };
    std::vector<EquInfo> equs = {
        {"Equ1", "L1Sig_Any", "H1a", "H2a", "_withBkg", "L1Inner"},
        {"Equ2", "L1Sig_Pass", "H1b", "H2b", "_withBkg", "L1Inner"},
        {"Equ3", "L1Sig_Frag", "H1c", "H2c", "_withBkg", "L1Inner"}
    };

    TH1D* ref_raw = GetHistD(fL1_Frag, TString::Format("h_yield_in_%s_window_from_%s_TOF_L1Sig_Any_L1Inner", cfg.source.c_str(), cfg.source.c_str()));
    if(!ref_raw) return;
    TH1D* ref = (TH1D*)ref_raw->Clone("ref_binned"); ref->Rebin(2);

    for(const auto& eq : equs) {
        std::cout << "  - Calculation for " << eq.name << std::endl;
        DetResults mcRes = CalculateMCRatio(cfg, ref, chain, eq.mcDen, eq.mcNum);
        TString pMFN = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H1_UseMass%d_FragFrom%s%s.root", cfg.fragAbbr.c_str(), chain.Data(), cfg.useMass, cfg.source.c_str(), eq.massSuffix.c_str());
        TString pMNN = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H1_UseMass%d%s.root", cfg.fragAbbr.c_str(), chain.Data(), cfg.useMass, eq.massSuffix.c_str());
        TFile *fMF = TFile::Open(pMFN); TFile *fMN = TFile::Open(pMNN);
        if(!fMF || !fMN) continue;

        DetResults issRes;
        for(auto iso : cfg.targetIsotopes) issRes.combine[iso] = CreateHist(("h_iss_"+eq.name+"_"+iso).c_str(), ref);

        for (int i = 1; i <= ref->GetNbinsX(); ++i) {
            double cent = ref->GetBinCenter(i);
            if(cent < 0.1 || cent > 22.0) continue;
            std::string det = (cent < 1.2) ? "TOF" : (cent < 3.5) ? "NaF" : "AGL";

            TString nVNum = TString::Format("h_fit_Be_Z%d_%s", cfg.sourceZ, det.c_str());
            TH1D* hV = GetHistD(fInner, nVNum);
            ValueWithError vNum(0,0);
            if(hV) { hV->Rebin(2); vNum = GetBinVal(hV, hV->FindBin(cent)); }

            TString nDen = TString::Format("h_yield_in_%s_window_from_%s_%s_%s_%s", cfg.source.c_str(), cfg.source.c_str(), det.c_str(), eq.l1Sig.c_str(), eq.denChain.c_str());
            TH1D* hDen = GetHistD(fL1_noFrag, nDen);
            if(hDen) hDen->Rebin(2); ValueWithError vDen = GetBinVal(hDen, hDen ? hDen->FindBin(cent) : 0);

            TString nCnt = TString::Format("h_yield_in_%s_window_from_%s_%s_%s_%s", cfg.source.c_str(), cfg.fragment.c_str(), det.c_str(), eq.l1Sig.c_str(), eq.denChain.c_str());
            TH1D* hC = GetHistD(fL1_Frag, nCnt);
            if(hC) hC->Rebin(2); ValueWithError vCnt = GetBinVal(hC, hC ? hC->FindBin(cent) : 0);

            if(vDen.val <= 1e-4) continue;

            std::map<std::string, ValueWithError> fracs_F, fracs_N;
            for(auto fitIso : {"Be7", "Be9"}) {
                TString kF = TString::Format("h_best_%s_frac_%s", fitIso, det.c_str());
                TH1D* hf = GetHistD(fMF, kF);
                fracs_F[fitIso] = hf ? GetBinVal(hf, hf->FindBin(cent)) : ValueWithError(0,0);
                TString kN = TString::Format("h_best_%s_frac_%s", fitIso, det.c_str());
                TH1D* hn = GetHistD(fMN, kN);
                fracs_N[fitIso] = hn ? GetBinVal(hn, hn->FindBin(cent)) : ValueWithError(0,0);
            }
            fracs_F["Be10"] = ValueWithError(1.0, 0.0) - fracs_F["Be7"] - fracs_F["Be9"];
            fracs_N["Be10"] = ValueWithError(1.0, 0.0) - fracs_N["Be7"] - fracs_N["Be9"];

            std::cout << Form("[%s Bin %d] vNum=%.1f | vDen=%.1f | vCnt=%.1f | fF10=%.3f", eq.name.c_str(), i, vNum.val, vDen.val, vCnt.val, fracs_F["Be10"].val) << std::endl;

            ValueWithError sumTrue(0,0);
            for (auto iso : cfg.targetIsotopes) {
                if (iso == "Total") continue;
                ValueWithError vTrue = (vNum * fracs_F[iso]) - (vCnt * fracs_N[iso]);
                issRes.combine[iso]->SetBinContent(i, (vTrue/vDen).val); 
                issRes.combine[iso]->SetBinError(i, (vTrue/vDen).err);
                sumTrue = sumTrue + vTrue;
            }
            issRes.combine["Total"]->SetBinContent(i, (sumTrue/vDen).val); 
            issRes.combine["Total"]->SetBinError(i, (sumTrue/vDen).err);
        }

        fOut->mkdir(eq.name.c_str()); fOut->cd(eq.name.c_str());
        for(auto const& [iso, h] : issRes.combine) {
            h->Write();
            TH1D* hm = mcRes.combine.count(iso) ? mcRes.combine[iso] : nullptr;
            if(hm) hm->Write();
            TString pBase = "/eos/user/z/zixuan/Isotope/BkgValid/Plots/" + TString::Format("%s_to_%s_%s_%s_%s", cfg.source.c_str(), cfg.fragment.c_str(), chain.Data(), eq.name.c_str(), iso.c_str());
            CreateComparisonPlots(h, hm, Form("%s %s %s", cfg.source.c_str(), eq.name.c_str(), iso.c_str()), pBase, 0);
        }
        fMF->Close(); fMN->Close();
    }
    fInner->Close(); fL1_noFrag->Close(); fL1_Frag->Close(); fOut->Close();
}

void CalFrag_new() {
    gROOT->SetBatch(kTRUE); gStyle->SetOptStat(0);
    
    // Boron -> Be
    AnalysisConfig cB; cB.source="Boron"; cB.fragment="Beryllium"; cB.fragAbbr="Be"; cB.useMass=7; cB.sourceZ=5;
    cB.srcAbundance={{"B10",0.3},{"B11",0.7}}; cB.mcFiles={{"B10","B10"},{"B11","B11"}};
    cB.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cB.targetIsotopes={"Be7","Be9","Be10","Total"};
    RunAnalysis("L1Inner", cB);

    // Carbon -> Be
    AnalysisConfig cC; cC.source="Carbon"; cC.fragment="Beryllium"; cC.fragAbbr="Be"; cC.useMass=7; cC.sourceZ=6;
    cC.srcAbundance={{"C12",1.0}}; cC.mcFiles={{"C12","C12"}};
    cC.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cC.targetIsotopes={"Be7","Be9","Be10","Total"};
    RunAnalysis("L1Inner", cC);

    // Nitrogen -> Be
    AnalysisConfig cN; cN.source="Nitrogen"; cN.fragment="Beryllium"; cN.fragAbbr="Be"; cN.useMass=7; cN.sourceZ=7;
    cN.srcAbundance={{"N14",0.5},{"N15",0.5}}; cN.mcFiles={{"N14","N14"},{"N15","N15"}};
    cN.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cN.targetIsotopes={"Be7","Be9","Be10","Total"};
    RunAnalysis("L1Inner", cN);

    // Oxygen -> Be
    AnalysisConfig cO; cO.source="Oxygen"; cO.fragment="Beryllium"; cO.fragAbbr="Be"; cO.useMass=7; cO.sourceZ=8;
    cO.srcAbundance={{"O16",1.0}}; cO.mcFiles={{"O16","O16"}};
    cO.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cO.targetIsotopes={"Be7","Be9","Be10","Total"};
    RunAnalysis("L1Inner", cO);
}