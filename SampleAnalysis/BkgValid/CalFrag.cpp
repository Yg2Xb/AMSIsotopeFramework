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
    double val = 0.0;
    double err = 0.0;
    ValueWithError() = default;
    ValueWithError(double v, double e) : val(v), err(e) {}
    ValueWithError operator+(const ValueWithError& o) const { return {val + o.val, std::hypot(err, o.err)}; }
    ValueWithError operator-(const ValueWithError& o) const { return {val - o.val, std::hypot(err, o.err)}; }
    ValueWithError operator*(double s) const { return {val * s, std::abs(err * s)}; }
    ValueWithError operator*(const ValueWithError& o) const { return {val * o.val, std::hypot(err * o.val, val * o.err)}; }
    ValueWithError operator/(const ValueWithError& o) const {
        if (o.val == 0) return {0, 0};
        double q = val / o.val;
        double e = std::hypot(err / o.val, (val * o.err) / (o.val * o.val));
        return {q, e};
    }
};

struct AnalysisConfig {
    std::string source;
    std::string fragment;
    std::string fragFileID;
    std::string fragAbbr;
    int useMass;
    std::map<std::string, double> srcAbundance;
    std::map<std::string, std::string> mcFiles;
    std::map<std::string, int> fragIsoMap;
    std::vector<std::string> fitIsotopes;
    std::vector<std::string> targetIsotopes;
};

struct DetResults {
    std::map<std::string, TH1D*> combine;
    std::map<std::string, TH1D*> tof;
    std::map<std::string, TH1D*> naf;
    std::map<std::string, TH1D*> agl;
};

struct MCSourceData {
    std::string isoName;
    double abundance;
    TH1* h_gen;
    TH1D* h_den_stitched;
    std::map<std::string, TH1D*> h_num_stitched; 
    std::map<std::string, TH1*> h_den_parts; 
    std::map<std::string, std::map<std::string, TH1*>> h_num_parts;
};

TH1* GetHist(TFile* f, TString name) {
    if (!f) return nullptr;
    TH1* h = (TH1*)f->Get(name);
    if (!h) {
        std::cout << "[ERROR] Hist Missing: " << name << " in " << f->GetName() << std::endl;
        return nullptr;
    }
    h->SetDirectory(0);
    return h;
}

TH1D* CreateHist(const char* name, TH1* ref) {
    TH1D* h = new TH1D(name, "", ref->GetNbinsX(), ref->GetXaxis()->GetXbins()->GetArray());
    h->SetDirectory(0);
    h->Sumw2();
    return h;
}

TH1D* StitchHists(const std::map<std::string, TH1*>& hists, const TH1* ref) {
    if(hists.empty() || !hists.begin()->second) return nullptr;
    TH1D* out = CreateHist(TString::Format("%s_stitched", hists.begin()->second->GetName()), (TH1*)ref);
    for (int i = 1; i <= out->GetNbinsX(); ++i) {
        double center = out->GetBinCenter(i);
        std::string det = (center < 1.2) ? "TOF" : (center < 3.2) ? "NaF" : "AGL";
        if (hists.count(det) && hists.at(det)) {
            int srcBin = hists.at(det)->FindBin(center);
            out->SetBinContent(i, hists.at(det)->GetBinContent(srcBin));
            out->SetBinError(i, hists.at(det)->GetBinError(srcBin));
        }
    }
    return out;
}

ValueWithError GetBinVal(TH1* h, int bin) {
    if (!h) return {0, 0};
    return {h->GetBinContent(bin), h->GetBinError(bin)};
}

ValueWithError GetValDebug(TH1* h, double energy, TString fName, TString desc) {
    if (!h) {
        std::cout << "    [FAIL] " << desc << ": Hist is NULL in " << fName << "\n";
        return {0,0};
    }
    int bin = h->FindBin(energy);
    double low = h->GetXaxis()->GetBinLowEdge(bin);
    double up = h->GetXaxis()->GetBinUpEdge(bin);
    ValueWithError v(h->GetBinContent(bin), h->GetBinError(bin));
    
    std::cout << "    " << std::left << std::setw(18) << desc 
              << " | File: " << std::setw(25) << fName.Data()
              << " | Hist: " << std::setw(40) << h->GetName()
              << " | Bin: " << std::setw(3) << bin 
              << " [" << std::fixed << std::setprecision(3) << low << "," << up << "]"
              << " | Val: " << std::scientific << std::setprecision(4) << v.val << "\n";
    return v;
}

// -------------------------------------------------------------------
// Modified DrawFit Function
// -------------------------------------------------------------------
void DrawFit(TH1* h, TLegend* leg, double min, double max, int color, TString label) {
    bool hasData = false;
    for(int b=h->GetXaxis()->FindBin(min); b<=h->GetXaxis()->FindBin(max); ++b) {
        if(h->GetBinContent(b) > 0) hasData = true;
    }
    if(!hasData) return;

    // --- Step 1: Attempt Constant Fit (pol0) ---
    TF1* f0 = new TF1(TString::Format("fit0_%s_%f", h->GetName(), min), "pol0", min, max);
    f0->SetLineColor(color);
    f0->SetLineWidth(2);
    // Use "SQN" first: Store result, Quiet, No draw (yet)
    TFitResultPtr r0 = h->Fit(f0, "SQN+"); 

    bool useConst = false;
    if (r0->IsValid()) {
        double ndf = r0->Ndf();
        double chi2 = r0->Chi2();
        // Check condition: Success AND Chi2/NDF < 1.5
        if (ndf > 0 && (chi2 / ndf) < 2) {
            useConst = true;
        } else if (ndf == 0 && chi2 < 1.0) {
            // Edge case: if NDF is 0 (e.g. 1 point), treat as good fit if chi2 is small
            useConst = true;
        }
    }

    if (useConst) {
        // --- Success with Constant Fit ---
        f0->Draw("SAME");
        double p0 = f0->GetParameter(0); double e0 = f0->GetParError(0);
        
        TString txt = TString::Format("%s: %.2f #pm %.2f (Const)", label.Data(), p0, e0);
        leg->AddEntry(f0, txt, "l");
        return; // Done, exit function
    }

    // --- Step 2: Fallback to Linear Fit (pol1) ---
    // Clean up the unused constant fit function object to avoid clutter
    delete f0; 

    TF1* f1 = new TF1(TString::Format("fit_%s_%f", h->GetName(), min), "pol1", min, max);
    f1->SetLineColor(color);
    f1->SetLineWidth(2);
    TFitResultPtr r1 = h->Fit(f1, "SRQ+"); // Draw this one ("R" for range, "Q" for quiet, "+" to add to list)
    
    if (r1->IsValid()) {
        double p0 = f1->GetParameter(0); double e0 = f1->GetParError(0);
        double p1 = f1->GetParameter(1); double e1 = f1->GetParError(1);
        f1->Draw("SAME");
        
        TString txt = TString::Format("%s: %.2f #pm %.2f + (%.2f #pm %.2f) #times E_{k}/n [GeV]", 
                                      label.Data(), p0, e0, p1, e1);
        
        leg->AddEntry(f1, txt, "l");
    }
}
// -------------------------------------------------------------------

void CreateComparisonPlots(TH1* h_iss, TH1* h_mc, TString title, TString outPath, int detMode) {
    if (!h_iss) return;

    auto getMinMax = [](TH1* h, double& minVal, double& maxVal) {
        if(!h) return;
        for(int i=1; i<=h->GetNbinsX(); ++i) {
            double c = h->GetBinCenter(i);
            if(c < 0.4 || c > 21.5) continue;
            double v = h->GetBinContent(i);
            if(v == 0) continue;
            if(v < minVal) minVal = v;
            if(v > maxVal) maxVal = v;
        }
    };

    TCanvas* c1 = new TCanvas("c1", title, 800, 600);
    c1->SetGrid();
    h_iss->SetMarkerStyle(20); h_iss->SetLineColor(kBlack);
    h_iss->SetTitle(title); h_iss->GetXaxis()->SetTitle("E_{k}/n [GeV]"); h_iss->GetYaxis()->SetTitle("Ratio");
    
    double yMin = 1e9, yMax = -1e9;
    getMinMax(h_iss, yMin, yMax);
    if(h_mc) getMinMax(h_mc, yMin, yMax);
    if(yMin > yMax) { yMin = 0; yMax = 1; }
    double diff = yMax - yMin; if(diff == 0) diff = 1.0;
    h_iss->GetYaxis()->SetRangeUser(yMin - diff * 0.2, yMax + diff * 0.2);
    
    double xMin = 0.4, xMax = 21.5;
    if (detMode == 1) { xMin = 0.3; xMax = 1.5; }
    else if (detMode == 2) { xMin = 0.8; xMax = 6.1; }
    else if (detMode == 3) { xMin = 2.5; xMax = 21.5; }
    
    h_iss->GetXaxis()->SetRangeUser(xMin, xMax);
    h_iss->Draw("PZ");
    
    TLegend* leg = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg->AddEntry(h_iss, "ISS Data", "ep");
    if (h_mc) {
        h_mc->SetMarkerStyle(20); h_mc->SetMarkerColor(kRed); h_mc->SetLineColor(kRed);
        h_mc->Draw("PZ SAME"); leg->AddEntry(h_mc, "MC Simulation", "ep");
    }
    leg->Draw();
    c1->SaveAs(outPath + "_comp.png"); delete c1;
    
    if (h_mc) {
        TCanvas* c2 = new TCanvas("c2", title + " Ratio", 800, 600); c2->SetGrid();
        TH1* h_ratio = (TH1*)h_iss->Clone("ratio"); h_ratio->Divide(h_mc);
        h_ratio->SetTitle(title); h_ratio->GetYaxis()->SetTitle("ISS/MC");
        h_ratio->GetXaxis()->SetRangeUser(xMin, xMax); 
        
        double rMin = 1e9, rMax = -1e9;
        getMinMax(h_ratio, rMin, rMax);
        if(rMin > rMax) { rMin = 0; rMax = 1; }
        double rDiff = rMax - rMin; if(rDiff == 0) rDiff = 1.0;
        h_ratio->GetYaxis()->SetRangeUser(rMin - rDiff * 1.2, rMax + rDiff * 1.2);

        h_ratio->Draw("PZ");
        TLegend* fitLeg = new TLegend(0.15, 0.78, 0.85, 0.88);
        fitLeg->SetBorderSize(0); fitLeg->SetFillStyle(0); fitLeg->SetTextSize(0.035);
        
        if (detMode == 0) { 
            DrawFit(h_ratio, fitLeg, 0.3, 7.0, kMagenta, "Low"); 
            DrawFit(h_ratio, fitLeg, 7.0, 21.5, kOrange+1, "High"); 
        } else if (detMode == 1) {
            DrawFit(h_ratio, fitLeg, 0.3, 1.5, kMagenta, "TOF");
        } else if (detMode == 2) {
            DrawFit(h_ratio, fitLeg, 0.8, 6.1, kMagenta, "NaF");
        } else if (detMode == 3) {
            DrawFit(h_ratio, fitLeg, 2.5, 7.0, kMagenta, "AGL-Low");
            DrawFit(h_ratio, fitLeg, 7.0, 21.5, kOrange+1, "AGL-High");
        }
        fitLeg->Draw(); c2->SaveAs(outPath + "_ratio.png"); delete c2; delete h_ratio;
    }
    delete leg;
}

DetResults CalculateMCRatio(const AnalysisConfig& cfg, TH1* refHist, TString chain, TString denSuffix) {
    std::cout << "\n>>> [MC] Calculation (" << denSuffix << "): " << cfg.source << " -> " << cfg.fragment << std::endl;
    DetResults res;
    TString mcPath = "/eos/user/z/zixuan/Isotope/Add/";
    std::vector<std::string> dets = {"TOF", "NaF", "AGL"};
    
    TH1D* sum_den = CreateHist("mc_sum_den", refHist);
    std::map<std::string, TH1D*> sum_num;
    for (auto const& [iso, m] : cfg.fragIsoMap) sum_num[iso] = CreateHist(("mc_sum_num_" + iso).c_str(), refHist);
    TH1D* sum_num_tot = CreateHist("mc_sum_num_tot", refHist);

    std::map<std::string, TH1D*> sum_den_d, sum_num_tot_d;
    std::map<std::string, std::map<std::string, TH1D*>> sum_num_d;
    for(auto d : dets) {
        sum_den_d[d] = CreateHist(("mc_sum_den_"+d).c_str(), refHist);
        sum_num_tot_d[d] = CreateHist(("mc_sum_num_tot_"+d).c_str(), refHist);
        for(auto const& [iso, m] : cfg.fragIsoMap) sum_num_d[d][iso] = CreateHist(("mc_sum_num_"+d+"_"+iso).c_str(), refHist);
    }

    std::vector<MCSourceData> mcDataList;
    for (auto const& [srcIso, srcFile] : cfg.mcFiles) {
        TString fullPath = mcPath + srcFile.c_str();
        TFile* f = TFile::Open(fullPath);
        if(!f || f->IsZombie()) continue;
        
        MCSourceData mcData;
        mcData.isoName = srcIso;
        mcData.abundance = cfg.srcAbundance.at(srcIso);
        
        TH1* hg = GetHist(f, "MC_FLUX_H3");
        if(hg) { hg->Rebin(2); mcData.h_gen = hg; } 
        else mcData.h_gen = nullptr;

        for(auto det : dets) {
            TString n = TString::Format("%s_BKG_%s_%s_%s", chain.Data(), denSuffix.Data(), cfg.source.c_str(), det.c_str());
            TH1* h = GetHist(f, n);
            if(h) h->Rebin(2);
            mcData.h_den_parts[det] = h;
        }
        mcData.h_den_stitched = StitchHists(mcData.h_den_parts, refHist);

        int fZ = (cfg.fragment == "Beryllium") ? 4 : 5;
        for (auto const& [iso, mass] : cfg.fragIsoMap) {
            for(auto det : dets) {
                TString n = TString::Format("%s_BKG_H2a2_%s_Z%d_Mass%d", chain.Data(), det.c_str(), fZ, mass);
                TH1* h = GetHist(f, n);
                if(h) h->Rebin(2);
                mcData.h_num_parts[iso][det] = h;
            }
            mcData.h_num_stitched[iso] = StitchHists(mcData.h_num_parts[iso], refHist);
        }
        mcDataList.push_back(mcData);
        f->Close(); delete f; 
    }

    std::cout << "\n---------------------------------------------------------------------------------------------------\n";
    std::cout << "MC DETAILED DEBUG OUTPUT (" << denSuffix << ")\n";
    std::cout << "---------------------------------------------------------------------------------------------------\n";

    for (int i = 1; i <= refHist->GetNbinsX(); ++i) {
        double low = refHist->GetXaxis()->GetBinLowEdge(i);
        double up = refHist->GetXaxis()->GetBinUpEdge(i);
        double cent = refHist->GetBinCenter(i);
        
        if (cent < 0.3 || cent > 21.5) continue;
        
        double totGen = 0;
        for(auto& d : mcDataList) if(d.h_gen) totGen += d.h_gen->GetBinContent(i);
        if (totGen <= 0) continue;

        std::cout << "\nBin " << i << " [" << std::fixed << std::setprecision(3) << low << ", " << up << "] TotalGen: " << totGen << "\n";

        for(auto& d : mcDataList) {
            double curGen = d.h_gen ? d.h_gen->GetBinContent(i) : 0;
            if(curGen <= 0) continue;
            double w = (totGen / curGen) * d.abundance;
            
            ValueWithError vDenRaw = GetBinVal(d.h_den_stitched, i);
            ValueWithError vDenW = vDenRaw * w;
            
            std::cout << "  Src " << d.isoName << " (Gen=" << curGen << ", Abd=" << d.abundance << ", W=" << w << "):\n";
            std::cout << "    DenRaw: " << vDenRaw.val << " -> WDen: " << vDenW.val << "\n";
            
            sum_den->SetBinContent(i, sum_den->GetBinContent(i) + vDenW.val);
            sum_den->SetBinError(i, std::hypot(sum_den->GetBinError(i), vDenW.err));

            for(auto det : dets) {
                double minE = (det=="TOF")?0.3:(det=="NaF")?0.8:2.5;
                double maxE = (det=="TOF")?1.5:(det=="NaF")?6.1:21.5;
                if(cent >= minE && cent <= maxE) {
                    ValueWithError vDenPart = GetBinVal(d.h_den_parts[det], d.h_den_parts[det]?d.h_den_parts[det]->FindBin(cent):0) * w;
                    sum_den_d[det]->SetBinContent(i, sum_den_d[det]->GetBinContent(i) + vDenPart.val);
                    sum_den_d[det]->SetBinError(i, std::hypot(sum_den_d[det]->GetBinError(i), vDenPart.err));
                }
            }

            for (auto const& [iso, h] : d.h_num_stitched) {
                ValueWithError vNumRaw = GetBinVal(h, i);
                ValueWithError vNumW = vNumRaw * w;
                std::cout << "    Frag " << iso << ": NumRaw: " << vNumRaw.val << " -> WNum: " << vNumW.val << "\n";

                sum_num[iso]->SetBinContent(i, sum_num[iso]->GetBinContent(i) + vNumW.val);
                sum_num[iso]->SetBinError(i, std::hypot(sum_num[iso]->GetBinError(i), vNumW.err));
                sum_num_tot->SetBinContent(i, sum_num_tot->GetBinContent(i) + vNumW.val);
                sum_num_tot->SetBinError(i, std::hypot(sum_num_tot->GetBinError(i), vNumW.err));

                for(auto det : dets) {
                    double minE = (det=="TOF")?0.3:(det=="NaF")?0.8:2.5;
                    double maxE = (det=="TOF")?1.5:(det=="NaF")?6.1:21.5;
                    if(cent >= minE && cent <= maxE) {
                        ValueWithError vNumPart = GetBinVal(d.h_num_parts[iso][det], d.h_num_parts[iso][det]?d.h_num_parts[iso][det]->FindBin(cent):0) * w;
                        sum_num_d[det][iso]->SetBinContent(i, sum_num_d[det][iso]->GetBinContent(i) + vNumPart.val);
                        sum_num_d[det][iso]->SetBinError(i, std::hypot(sum_num_d[det][iso]->GetBinError(i), vNumPart.err));
                        sum_num_tot_d[det]->SetBinContent(i, sum_num_tot_d[det]->GetBinContent(i) + vNumPart.val);
                        sum_num_tot_d[det]->SetBinError(i, std::hypot(sum_num_tot_d[det]->GetBinError(i), vNumPart.err));
                    }
                }
            }
        }
        
        double finalDen = sum_den->GetBinContent(i);
        double finalNumTot = sum_num_tot->GetBinContent(i);
        double finalRatio = (finalDen > 0) ? finalNumTot / finalDen : 0;
        std::cout << "  => MIXED TOTAL: Den=" << finalDen << " Num=" << finalNumTot << " Ratio=" << finalRatio << "\n";
    }

    for(auto& d : mcDataList) {
        if(d.h_gen) delete d.h_gen;
        if(d.h_den_stitched) delete d.h_den_stitched;
        for(auto& p : d.h_num_stitched) delete p.second;
        for(auto& p : d.h_den_parts) if(p.second) delete p.second;
        for(auto& px : d.h_num_parts) for(auto& py : px.second) if(py.second) delete py.second;
    }

    for(auto const& [iso, h] : sum_num) {
        TH1D* r = CreateHist(("h_mc_ratio_"+iso).c_str(), refHist); r->Divide(h, sum_den);
        res.combine[iso] = r; delete h;
    }
    TH1D* rt = CreateHist("h_mc_ratio_Total", refHist); rt->Divide(sum_num_tot, sum_den);
    res.combine["Total"] = rt; delete sum_num_tot; delete sum_den;

    for(auto d : dets) {
        for(auto const& [iso, h] : sum_num_d[d]) {
            TH1D* r = CreateHist(("h_mc_ratio_"+d+"_"+iso).c_str(), refHist); r->Divide(h, sum_den_d[d]);
            if(d=="TOF") res.tof[iso] = r; else if(d=="NaF") res.naf[iso] = r; else res.agl[iso] = r;
            delete h;
        }
        TH1D* rtd = CreateHist(("h_mc_ratio_"+d+"_Total").c_str(), refHist); rtd->Divide(sum_num_tot_d[d], sum_den_d[d]);
        if(d=="TOF") res.tof["Total"] = rtd; else if(d=="NaF") res.naf["Total"] = rtd; else res.agl["Total"] = rtd;
        delete sum_num_tot_d[d]; delete sum_den_d[d];
    }
    
    return res;
}

void RunAnalysis(TString chain, const AnalysisConfig& cfg) {
    std::cout << "\n======================================================\n";
    std::cout << "[INFO] Starting Analysis for " << cfg.source << " -> " << cfg.fragment << std::endl;

    TString pL1 = TString::Format("/eos/user/z/zixuan/Isotope/ChargeTemp/PureQFit_%sTo%s_UnbiasedL1Inner.root", cfg.source.c_str(), cfg.fragment.c_str());
    TString pCnt = TString::Format("/eos/user/z/zixuan/Isotope/Add/%s.root", cfg.fragFileID.c_str());
    TString pOut = TString::Format("/eos/user/z/zixuan/Isotope/BkgValid/%s_to_%s_%s_Validation.root", cfg.source.c_str(), cfg.fragment.c_str(), chain.Data());
    
    TFile* fL1 = TFile::Open(pL1);
    if(!fL1 || fL1->IsZombie()) std::cout << "[ERROR] Cannot open L1 file: " << pL1 << std::endl;

    TFile* fCnt = TFile::Open(pCnt);
    if(!fCnt || fCnt->IsZombie()) std::cout << "[ERROR] Cannot open Add file: " << pCnt << std::endl;

    TFile* fOut = TFile::Open(pOut, "RECREATE");
    
    TString pFitF = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d_FragFrom%s.root", cfg.fragAbbr.c_str(), chain.Data(), cfg.useMass, cfg.source.c_str());
    TString pFitN = TString::Format("/eos/user/z/zixuan/Isotope/MassTempFit/wide_MassTF_%s_%s_H2_UseMass%d.root", cfg.fragAbbr.c_str(), chain.Data(), cfg.useMass);
    
    TFile* fMF = TFile::Open(pFitF);
    if(!fMF || fMF->IsZombie()) std::cout << "[ERROR] Cannot open Frag MassFit: " << pFitF << std::endl;
    
    TFile* fMN = TFile::Open(pFitN);
    if(!fMN || fMN->IsZombie()) std::cout << "[ERROR] Cannot open Norm MassFit: " << pFitN << std::endl;
    
    TH1* ref = GetHist(fCnt, TString::Format("%s_BKG_H1a_%s_TOF", chain.Data(), cfg.source.c_str()));
    if(!ref) return;
    ref->Rebin(2);
    
    // Loop for Equ1, Equ2, Equ3
    struct EquInfo { std::string name; std::string mcDen; std::string l1Sig; };
    std::vector<EquInfo> equs = {
        {"Equ1", "H1a", "L1Sig_Any"},
        {"Equ2", "H1b", "L1Sig_Pass"},
        {"Equ3", "H1c", "L1Sig_Frag"}
    };

    for(const auto& eq : equs) {
        std::cout << "\n>>> Processing " << eq.name << " <<<\n";
        DetResults mcRes = CalculateMCRatio(cfg, ref, chain, eq.mcDen.c_str());
        DetResults issRes;
        
        for(auto iso : cfg.targetIsotopes) {
            issRes.combine[iso] = CreateHist(("h_iss_ratio_"+iso).c_str(), ref);
            issRes.tof[iso] = CreateHist(("h_iss_ratio_TOF_"+iso).c_str(), ref);
            issRes.naf[iso] = CreateHist(("h_iss_ratio_NaF_"+iso).c_str(), ref);
            issRes.agl[iso] = CreateHist(("h_iss_ratio_AGL_"+iso).c_str(), ref);
        }

        std::cout << "\n--- ISS DEBUG (" << eq.name << ") ---\n";
        std::vector<std::string> dets = {"TOF", "NaF", "AGL"};
        
        for (int i = 1; i <= ref->GetNbinsX(); ++i) {
            double cent = ref->GetBinCenter(i);
            if (cent < 0.3 || cent > 21.5) continue;
            
            // Combine Logic
            std::string detC = (cent < 1.2) ? "TOF" : (cent < 3.2) ? "NaF" : "AGL";
            
            TString nNum = TString::Format("%s_BKG_H2a_%s_%s", chain.Data(), cfg.source.c_str(), detC.c_str());
            TH1* hNum = GetHist(fCnt, nNum); if(hNum) hNum->Rebin(2); 
            ValueWithError vNum = GetValDebug(hNum, cent, "Add(H2a)", "RawNum"); if(hNum) delete hNum;

            TString nDen = TString::Format("h_yield_in_%s_from_%s_%s_%s_%s", cfg.source.c_str(), cfg.source.c_str(), detC.c_str(), eq.l1Sig.c_str(), chain.Data());
            TH1* hDen = GetHist(fL1, nDen); if(hDen) hDen->Rebin(2);
            ValueWithError vDen = GetValDebug(hDen, cent, "PureQFit", "L1Den"); if(hDen) delete hDen;

            TString nCt = TString::Format("h_yield_in_%s_from_%s_%s_%s_%s", cfg.source.c_str(), cfg.fragment.c_str(), detC.c_str(), eq.l1Sig.c_str(), chain.Data());
            TH1* hC = GetHist(fL1, nCt); if(hC) hC->Rebin(2);
            ValueWithError vCnt = GetValDebug(hC, cent, "PureQFit", "L1Contam"); if(hC) delete hC;

            ValueWithError sumNum(0,0);
            ValueWithError sumFracF(0,0), sumFracN(0,0);

            if (vDen.val > 0) {
                for (auto iso : cfg.targetIsotopes) {
                    if (iso == "Total") continue;
                    ValueWithError fF(0,0), fN(0,0);
                    bool isFit = false; for(auto fitIso : cfg.fitIsotopes) if(iso == fitIso) isFit = true;
                    
                    auto GetF = [&](TFile* f, TString nIso, TString detName) {
                        if(!f) return ValueWithError(0,0);
                        TString n = TString::Format("h_best_%s_frac_%s", nIso.Data(), detName.Data());
                        TH1* h = GetHist(f, n); 
                        ValueWithError v = h ? ValueWithError(h->GetBinContent(h->FindBin(cent)), h->GetBinError(h->FindBin(cent))) : ValueWithError(0,0);
                        if(h) delete h; return v;
                    };

                    if(isFit) {
                        fF = GetF(fMF, iso.c_str(), detC.c_str()); sumFracF = sumFracF + fF;
                        fN = GetF(fMN, iso.c_str(), detC.c_str()); sumFracN = sumFracN + fN;
                    } else {
                        fF = ValueWithError(1.0, 0.0) - sumFracF;
                        fN = ValueWithError(1.0, 0.0) - sumFracN;
                    }
                    ValueWithError vTrue = (vNum * fF) - (vCnt * fN);
                    ValueWithError rat = vTrue / vDen;
                    issRes.combine[iso]->SetBinContent(i, rat.val); issRes.combine[iso]->SetBinError(i, rat.err);
                    sumNum = sumNum + vTrue;
                }
                ValueWithError rTot = sumNum / vDen;
                issRes.combine["Total"]->SetBinContent(i, rTot.val); issRes.combine["Total"]->SetBinError(i, rTot.err);
            }

            // Individual Logic
            for(auto d : dets) {
                double minE = (d=="TOF")?0.3:(d=="NaF")?0.8:2.5;
                double maxE = (d=="TOF")?1.5:(d=="NaF")?6.1:21.5;
                if(cent < minE || cent > maxE) continue;

                TString nNumD = TString::Format("%s_BKG_H2a_%s_%s", chain.Data(), cfg.source.c_str(), d.c_str());
                TH1* hNumD = GetHist(fCnt, nNumD); if(hNumD) hNumD->Rebin(2);
                ValueWithError vNumD = GetBinVal(hNumD, hNumD?hNumD->FindBin(cent):0); if(hNumD) delete hNumD;

                TString nDenD = TString::Format("h_yield_in_%s_from_%s_%s_%s_%s", cfg.source.c_str(), cfg.source.c_str(), d.c_str(), eq.l1Sig.c_str(), chain.Data());
                TH1* hDenD = GetHist(fL1, nDenD); if(hDenD) hDenD->Rebin(2);
                ValueWithError vDenD = GetBinVal(hDenD, hDenD?hDenD->FindBin(cent):0); if(hDenD) delete hDenD;

                TString nCtD = TString::Format("h_yield_in_%s_from_%s_%s_%s_%s", cfg.source.c_str(), cfg.fragment.c_str(), d.c_str(), eq.l1Sig.c_str(), chain.Data());
                TH1* hCD = GetHist(fL1, nCtD); if(hCD) hCD->Rebin(2);
                ValueWithError vCntD = GetBinVal(hCD, hCD?hCD->FindBin(cent):0); if(hCD) delete hCD;

                if(vDenD.val <= 0) continue;
                ValueWithError sumNumD(0,0);
                ValueWithError sFF(0,0), sFN(0,0);

                for (auto iso : cfg.targetIsotopes) {
                    if (iso == "Total") continue;
                    ValueWithError fF(0,0), fN(0,0);
                    bool isFit = false; for(auto fitIso : cfg.fitIsotopes) if(iso == fitIso) isFit = true;
                    
                    auto GetFD = [&](TFile* f, TString nIso) {
                        if(!f) return ValueWithError(0,0);
                        TString n = TString::Format("h_best_%s_frac_%s", nIso.Data(), d.c_str());
                        TH1* h = GetHist(f, n); 
                        ValueWithError v = h ? ValueWithError(h->GetBinContent(h->FindBin(cent)), h->GetBinError(h->FindBin(cent))) : ValueWithError(0,0);
                        if(h) delete h; return v;
                    };

                    if(isFit) {
                        fF = GetFD(fMF, iso.c_str()); sFF = sFF + fF;
                        fN = GetFD(fMN, iso.c_str()); sFN = sFN + fN;
                    } else {
                        fF = ValueWithError(1.0, 0.0) - sFF;
                        fN = ValueWithError(1.0, 0.0) - sFN;
                    }
                    ValueWithError vTrue = (vNumD * fF) - (vCntD * fN);
                    ValueWithError rat = vTrue / vDenD;
                    
                    if(d=="TOF") { issRes.tof[iso]->SetBinContent(i, rat.val); issRes.tof[iso]->SetBinError(i, rat.err); }
                    else if(d=="NaF") { issRes.naf[iso]->SetBinContent(i, rat.val); issRes.naf[iso]->SetBinError(i, rat.err); }
                    else { issRes.agl[iso]->SetBinContent(i, rat.val); issRes.agl[iso]->SetBinError(i, rat.err); }
                    sumNumD = sumNumD + vTrue;
                }
                ValueWithError rTotD = sumNumD / vDenD;
                if(d=="TOF") { issRes.tof["Total"]->SetBinContent(i, rTotD.val); issRes.tof["Total"]->SetBinError(i, rTotD.err); }
                else if(d=="NaF") { issRes.naf["Total"]->SetBinContent(i, rTotD.val); issRes.naf["Total"]->SetBinError(i, rTotD.err); }
                else { issRes.agl["Total"]->SetBinContent(i, rTotD.val); issRes.agl["Total"]->SetBinError(i, rTotD.err); }
            }
        }
        
// 修改后的 SaveSet：分级创建目录，确保保存成功
        auto SaveSet = [&](std::map<std::string, TH1D*>& iss, std::map<std::string, TH1D*>& mc, TString dir, int m) {
            fOut->cd(); // 1. 先回到文件顶层
            
            // 2. 检查并创建第一级目录 (例如 Combine, TOF, NaF...)
            if (!fOut->GetDirectory(dir)) {
                fOut->mkdir(dir);
            }
            fOut->cd(dir); // 进入第一级

            // 3. 检查并创建第二级目录 (例如 Equ1, Equ2...)
            // 注意：eq.name 是 std::string，需要 .c_str()
            if (!gDirectory->GetDirectory(eq.name.c_str())) {
                gDirectory->mkdir(eq.name.c_str());
            }
            gDirectory->cd(eq.name.c_str()); // 进入第二级

            // 4. 保存直方图 (逻辑不变)
            for(auto const& [iso, h] : iss) {
                h->Write(); // 写入 ISS 结果
                
                TH1* hm = mc.count(iso) ? mc[iso] : nullptr;
                if(hm) hm->Write(); // 写入 MC 结果
                
                // 画图逻辑 (逻辑不变)
                TString tit = TString::Format("%s->%s %s %s %s", cfg.source.c_str(), cfg.fragment.c_str(), dir.Data(), eq.name.c_str(), iso.c_str());
                TString p = "/eos/user/z/zixuan/Isotope/BkgValid/Plots/" + TString::Format("%s_to_%s_%s_%s_%s_%s", cfg.source.c_str(), cfg.fragment.c_str(), chain.Data(), dir.Data(), eq.name.c_str(), iso.c_str());
                CreateComparisonPlots(h, hm, tit, p, m);
            }
        };
        
        SaveSet(issRes.combine, mcRes.combine, "Combine", 0);
        SaveSet(issRes.tof, mcRes.tof, "TOF", 1);
        SaveSet(issRes.naf, mcRes.naf, "NaF", 2);
        SaveSet(issRes.agl, mcRes.agl, "AGL", 3);
        
        for(auto p : issRes.combine) delete p.second; for(auto p : issRes.tof) delete p.second;
        for(auto p : issRes.naf) delete p.second; for(auto p : issRes.agl) delete p.second;
        for(auto p : mcRes.combine) delete p.second; for(auto p : mcRes.tof) delete p.second;
        for(auto p : mcRes.naf) delete p.second; for(auto p : mcRes.agl) delete p.second;
    }
    
    if(fL1) fL1->Close(); if(fCnt) fCnt->Close(); if(fOut) fOut->Close(); if(fMF) fMF->Close(); if(fMN) fMN->Close();
}

void CalFrag() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPadTopMargin(0.1);
    
    AnalysisConfig cB; cB.source="Boron"; cB.fragment="Beryllium"; cB.fragAbbr="Be"; cB.fragFileID="Be_frag4"; cB.useMass=7;
    cB.srcAbundance={{"B10",0.3},{"B11",0.7}}; cB.mcFiles={{"B10","B10_rew_frag4.root"},{"B11","B11_rew_frag4.root"}};
    cB.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cB.fitIsotopes={"Be7","Be9"}; cB.targetIsotopes={"Be7","Be9","Be10","Total"};
    
    AnalysisConfig cC; cC.source="Carbon"; cC.fragment="Beryllium"; cC.fragAbbr="Be"; cC.fragFileID="Be_frag4"; cC.useMass=7;
    cC.srcAbundance={{"C12",1.0}}; cC.mcFiles={{"C12","C12_rew_frag4.root"}};
    cC.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cC.fitIsotopes={"Be7","Be9"}; cC.targetIsotopes={"Be7","Be9","Be10","Total"};
    
    AnalysisConfig cN; cN.source="Nitrogen"; cN.fragment="Beryllium"; cN.fragAbbr="Be"; cN.fragFileID="Be_frag4"; cN.useMass=7;
    cN.srcAbundance={{"N14",0.5},{"N15",0.5}}; cN.mcFiles={{"N14","N14_rew_frag4.root"},{"N15","N15_rew_frag4.root"}};
    cN.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cN.fitIsotopes={"Be7","Be9"}; cN.targetIsotopes={"Be7","Be9","Be10","Total"};
    
    AnalysisConfig cO; cO.source="Oxygen"; cO.fragment="Beryllium"; cO.fragAbbr="Be"; cO.fragFileID="Be_frag4"; cO.useMass=7;
    cO.srcAbundance={{"O16",1.0}}; cO.mcFiles={{"O16","O16_rew_frag4.root"}};
    cO.fragIsoMap={{"Be7",7},{"Be9",9},{"Be10",10}}; cO.fitIsotopes={"Be7","Be9"}; cO.targetIsotopes={"Be7","Be9","Be10","Total"};

    RunAnalysis("UnbiasedL1Inner", cB);
    RunAnalysis("UnbiasedL1Inner", cC);
    RunAnalysis("UnbiasedL1Inner", cN);
    RunAnalysis("UnbiasedL1Inner", cO);
}