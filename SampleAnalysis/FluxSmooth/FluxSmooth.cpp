#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <numeric> // For std::accumulate
#include <cmath>   // For std::fabs

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"

// Assuming Tool.h is in the parent directory and its functions now accept double for mass number.
#include "../Tool.h"  

using namespace std;
using namespace AMS_Iso;

// ---------------- Configuration ----------------
static const string indir   = "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg";
static const string outPdf  = "/eos/user/z/zixuan/Isotope/FluxSmooth/FluxSmooth.pdf";
static const string outRoot = "/eos/user/z/zixuan/Isotope/FluxSmooth/FluxSmooth.root";

// ---------------- Data sources -----------------
struct Species {
    string elem;
    int    Z;
    string path;
};

static const vector<Species> speciesList = {
    {"Be", 4, indir + "/AMS2011to2018PhysReport_BerylliumFlux.root"},
    {"B",  5, indir + "/AMS2011to2018PhysReport_BoronFlux.root"},
    {"C",  6, indir + "/AMS2011to2018PhysReport_CarbonFlux.root"},
    {"N",  7, indir + "/AMS2011to2018PhysReport_NitrogenFlux.root"},
    {"O",  8, indir + "/AMS2011to2018PhysReport_OxygenFlux.root"}
};

// ---------------- Isotope fractions ------------
static map<string, vector<pair<int, double>>> getIsotopeFractions() {
    return {
        {"Be", {{7, 0.6}, {9, 0.3}, {10, 0.1}}},
        {"B",  {{10, 0.3}, {11, 0.7}}},
        {"C",  {{12, 1.0}}},
        {"N",  {{14, 0.5}, {15, 0.5}}},
        {"O",  {{16, 1.0}}}
    };
}

// ---------------- Spline nodes -----------------
static vector<double> generateXPoints(const TGraphAsymmErrors* g) {
    vector<double> pts;
    if (!g) return pts;
    const int n = g->GetN();
    if (n <= 0) return pts;
    const double* xs = g->GetX();
    int i = 0;
    while (i < n - 3) {
        pts.push_back(xs[i]);
        if (i < 10) i += 6;
        else if (i < 50) i += 7;
        else if (i < 200) i += 6;
        else i += 3;
    }
    if (n >= 2) {
        pts.push_back(xs[n - 2]);
        pts.push_back(xs[n - 1]);
    } else if (n > 0) {
        pts.push_back(xs[n - 1]);
    }
    vector<double> uniq;
    if (!pts.empty()) {
        uniq.reserve(pts.size());
        uniq.push_back(pts[0]);
        for (size_t k = 1; k < pts.size(); ++k) {
            if (pts[k] > pts[k - 1] + 1e-9) uniq.push_back(pts[k]);
        }
    }
    return uniq;
}

// ---------------- Convert R->Ek/n --------------
static unique_ptr<TGraphAsymmErrors> convertRGraphToEk(const TGraphAsymmErrors* gR,
                                                       int Z, double A,
                                                       const string& name) {
    if (!gR) return nullptr;
    const int n = gR->GetN();
    auto g = make_unique<TGraphAsymmErrors>(n);
    g->SetName(name.c_str());

    for (int i = 0; i < n; ++i) {
        const double R  = gR->GetX()[i];
        const double yR = gR->GetY()[i];

        const double Ek  = rigidityToKineticEnergy(R, Z, A);
        const double jac = dR_dEk(Ek, Z, A);
        const double yEk = yR * jac;

        g->SetPoint(i, Ek, yEk);
        g->SetPointError(i, 0.0, 0.0,
                         gR->GetErrorYlow(i)  * fabs(jac),
                         gR->GetErrorYhigh(i) * fabs(jac));
    }
    return g;
}

// ---------------- Spline wrapper ---------------
static TF1* doSplineFit(TGraphAsymmErrors* g, const string& name, double xmin, double xmax, int color) {
    if (!g || g->GetN() < 2) { cerr << "[WARN] Not enough points to fit spline for " << name << endl; return nullptr; }
    auto nodes = generateXPoints(g);
    cout << "=== " << name << " ===" << endl;
    cout << "Data points: " << g->GetN() << endl;
    cout << "Spline nodes: " << nodes.size() << endl;
    if (nodes.size() < 3) {
        nodes.clear();
        const int n = g->GetN();
        const double* xs = g->GetX();
        for (int i = 0; i < n; ++i) nodes.push_back(xs[i]);
        if (nodes.size() < 3) { cerr << "[WARN] Still not enough points for spline " << name << ". Skipping fit." << endl; return nullptr; }
    }
    TF1* f = SplineFit(g, nodes.data(), nodes.size(), 0x38, "b2e2", name.c_str(), xmin, xmax);
    if (f) { f->SetLineWidth(2); f->SetLineColor(color); f->SetNpx(4000); }
    return f;
}

// ---------------- Y-range finder ---------------
static pair<double, double> findYRangePos(const vector<TGraphAsymmErrors*>& gs) {
    double ymin = numeric_limits<double>::infinity();
    double ymax = 0.0;
    for (auto* g : gs) {
        if (!g) continue;
        const int n = g->GetN();
        for (int i = 0; i < n; ++i) {
            const double y = g->GetY()[i];
            if (y > 0.0) { ymin = min(ymin, y); ymax = max(ymax, y); }
        }
    }
    if (!isfinite(ymin) || ymax <= 0.0) { ymin = 1e-12; ymax = 10.0; }
    ymin = max(ymin, 1e-30);
    const double lo = 0.5 * ymin;
    double hi = 10.0 * ymax;
    if (hi <= lo) hi = lo * 10.0;
    return {lo, hi};
}

// ---------------- Main -------------------------
void FluxSmooth() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    TCanvas c_dummy("c_dummy", "dummy", 800, 600);
    c_dummy.Print((outPdf + "[").c_str());

    auto fracMap = getIsotopeFractions();
    
    // Central storage for all graphs and functions to be saved
    map<string, unique_ptr<TGraphAsymmErrors>> graphs_to_save;
    map<string, TF1*> functions_to_save;

    const int colors[] = {kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 1, kAzure + 2, kViolet + 1};

    for (const auto& sp : speciesList) {
        cout << "\n\n======================================================\n";
        cout << "Processing Element: " << sp.elem << " (Z=" << sp.Z << ")\n";
        cout << "======================================================\n";
        
        TFile fin(sp.path.c_str(), "READ");
        if (fin.IsZombie()) { cerr << "[ERROR] Cannot open " << sp.path << endl; continue; }

        auto gTot_R = dynamic_cast<TGraphAsymmErrors*>(fin.Get("graph1"));
        if (!gTot_R) { cerr << "[ERROR] Missing graph1 in " << sp.path << endl; continue; }
        
        // --- Store total raw flux for saving ---
        string gTot_R_name = sp.elem + "_total_flux_R";
        gTot_R->SetName(gTot_R_name.c_str());
        graphs_to_save[gTot_R_name] = unique_ptr<TGraphAsymmErrors>((TGraphAsymmErrors*)gTot_R->Clone());

        auto fit_params = fracMap.find(sp.elem);
        if (fit_params == fracMap.end() || fit_params->second.empty()) { cerr << "[WARN] No isotope fractions for " << sp.elem << endl; continue; }
        
        const auto& isotopes = fit_params->second;
        const int n_isotopes = isotopes.size();
        
        // --- STEP 1: Convert R-Flux to Ek-Flux using a WEIGHTED AVERAGE MASS ---
        double weighted_A_sum = 0.0;
        double sum_of_ratios = 0.0;
        for (const auto& iso : isotopes) {
            weighted_A_sum += iso.first * iso.second;
            sum_of_ratios += iso.second;
        }
        const double A_avg = weighted_A_sum / sum_of_ratios;

        cout << "Step 1: Converting total R flux to an approximate total Ek/n flux..." << endl;
        auto gTot_Ek_approx = convertRGraphToEk(gTot_R, sp.Z, A_avg, sp.elem + "_total_Ek_approx");

        // --- STEP 2: Split Approximate Total Ek-Flux using Initial Constant Ratios ---
        cout << "Step 2: Splitting approximate total Ek/n flux with constant ratios..." << endl;
        vector<unique_ptr<TGraphAsymmErrors>> initial_graphs_Ek(n_isotopes);
        for (int i = 0; i < n_isotopes; ++i) {
            const double frac = isotopes[i].second / sum_of_ratios;
            const string tag = sp.elem + to_string(isotopes[i].first);
            
            initial_graphs_Ek[i] = make_unique<TGraphAsymmErrors>(gTot_Ek_approx->GetN());
            initial_graphs_Ek[i]->SetName((tag + "_initial_Ek").c_str());
            for(int p = 0; p < gTot_Ek_approx->GetN(); ++p) {
                initial_graphs_Ek[i]->SetPoint(p, gTot_Ek_approx->GetX()[p], gTot_Ek_approx->GetY()[p] * frac);
                initial_graphs_Ek[i]->SetPointError(p, 0, 0, gTot_Ek_approx->GetErrorYlow(p) * frac, gTot_Ek_approx->GetErrorYhigh(p) * frac);
            }
        }

        // --- STEP 3: Smooth Isotope Ek-Fluxes to Create Continuous Models ---
        cout << "Step 3: Smoothing initial isotope Ek/n fluxes..." << endl;
        vector<TF1*> splines_Ek(n_isotopes);
        bool step3_failed = false;
        for (int i = 0; i < n_isotopes; ++i) {
            const string tag = sp.elem + to_string(isotopes[i].first);
            splines_Ek[i] = doSplineFit(initial_graphs_Ek[i].get(), tag + "_spline_Ek", 0.08, 1800, colors[i % 6]);
            if (!splines_Ek[i]) {
                cerr << "[FATAL] Failed to create spline for " << tag << ". Aborting this element." << endl;
                step3_failed = true;
                break;
            }
        }
        if (step3_failed) continue;

        // --- NEW: PLOT Ek/n FLUXES (DIAGNOSTIC PLOT) ---
        {
            auto cEk = new TCanvas(("cEk_" + sp.elem).c_str(), (sp.elem + " Ek/n").c_str(), 900, 700);
            cEk->SetLogx(); cEk->SetLogy();

            vector<TGraphAsymmErrors*> gs_for_ek_range;
            gs_for_ek_range.push_back(gTot_Ek_approx.get());
            for(const auto& g : initial_graphs_Ek) gs_for_ek_range.push_back(g.get());
            auto [yminEk, ymaxEk] = findYRangePos(gs_for_ek_range);

            auto frameEk = (TGraphAsymmErrors*)gTot_Ek_approx->Clone((sp.elem + "_frame_Ek").c_str());
            frameEk->SetMarkerSize(0); frameEk->SetLineColor(0);
            frameEk->GetXaxis()->SetTitle("Kinetic Energy per Nucleon E_{k}/n [GeV/n]");
            frameEk->GetYaxis()->SetTitle("Flux [(m^{2} sr s GeV/n)^{-1}]");
            frameEk->GetYaxis()->SetRangeUser(yminEk, ymaxEk);
            frameEk->GetXaxis()->SetRangeUser(0.08, 2000);
            frameEk->Draw("AP");

            auto legEk = new TLegend(0.64, 0.65, 0.84, 0.84);
            legEk->SetBorderSize(0); legEk->SetFillStyle(0);
            
            gTot_Ek_approx->SetMarkerStyle(20);
            gTot_Ek_approx->SetMarkerColor(kBlack);
            gTot_Ek_approx->SetLineColor(kBlack);
            gTot_Ek_approx->Draw("P SAME");
            legEk->AddEntry(gTot_Ek_approx.get(), "Total Flux (Approx. Ek)", "p");

            for (int i = 0; i < n_isotopes; ++i) {
                initial_graphs_Ek[i]->SetMarkerStyle(24);
                initial_graphs_Ek[i]->SetMarkerColor(colors[i % 6]);
                initial_graphs_Ek[i]->SetLineColor(colors[i % 6]);
                initial_graphs_Ek[i]->Draw("P SAME");
                if (splines_Ek[i]) splines_Ek[i]->Draw("SAME");
                legEk->AddEntry(initial_graphs_Ek[i].get(), (sp.elem + to_string(isotopes[i].first) + " (Initial Split)").c_str(), "p");
            }
            legEk->Draw();
            cEk->Print(outPdf.c_str());
            delete legEk; delete frameEk; delete cEk;
        }
        
        // --- Store Ek/n graphs and splines for saving ---
        graphs_to_save[gTot_Ek_approx->GetName()] = std::move(gTot_Ek_approx);
        for(auto& g : initial_graphs_Ek) { graphs_to_save[g->GetName()] = std::move(g); }
        for(auto* f : splines_Ek) { if(f) functions_to_save[f->GetName()] = f; }

        // --- STEP 4: Calculate Final Isotope R-Fluxes with Dynamic Ratios ---
        cout << "Step 4: Re-splitting ORIGINAL total R-flux with dynamic ratios..." << endl;
        vector<unique_ptr<TGraphAsymmErrors>> temp_final_graphs_R(n_isotopes);
        for(int i=0; i<n_isotopes; ++i) {
            temp_final_graphs_R[i] = make_unique<TGraphAsymmErrors>(gTot_R->GetN());
            temp_final_graphs_R[i]->SetName((sp.elem + to_string(isotopes[i].first) + "_final_R").c_str());
        }

        for (int p = 0; p < gTot_R->GetN(); ++p) {
            const double R = gTot_R->GetX()[p];
            const double flux_R_total = gTot_R->GetY()[p];
            vector<double> dynamic_flux_values(n_isotopes);
            double sum_dynamic_flux = 0.0;
            for (int i = 0; i < n_isotopes; ++i) {
                const int A_i = isotopes[i].first;
                const double Ek_i = rigidityToKineticEnergy(R, sp.Z, A_i);
                if (splines_Ek[i] && Ek_i > splines_Ek[i]->GetXmin() && Ek_i < splines_Ek[i]->GetXmax()) {
                    dynamic_flux_values[i] = splines_Ek[i]->Eval(Ek_i);
                    if (dynamic_flux_values[i] < 0) dynamic_flux_values[i] = 0;
                    sum_dynamic_flux += dynamic_flux_values[i];
                } else {
                    dynamic_flux_values[i] = 0;
                }
            }
            if (sum_dynamic_flux > 1e-30) {
                for (int i = 0; i < n_isotopes; ++i) {
                    const double dynamic_frac = dynamic_flux_values[i] / sum_dynamic_flux;
                    temp_final_graphs_R[i]->SetPoint(p, R, flux_R_total * dynamic_frac);
                    temp_final_graphs_R[i]->SetPointError(p, gTot_R->GetErrorXlow(p), gTot_R->GetErrorXhigh(p),
                                                          gTot_R->GetErrorYlow(p) * dynamic_frac, gTot_R->GetErrorYhigh(p) * dynamic_frac);
                }
            } else {
                 for (int i = 0; i < n_isotopes; ++i) {
                    temp_final_graphs_R[i]->SetPoint(p, R, 0);
                    temp_final_graphs_R[i]->SetPointError(p, gTot_R->GetErrorXlow(p), gTot_R->GetErrorXhigh(p), 0, 0);
                }
            }
        }
        
        for(int i=0; i<n_isotopes; ++i) {
            const string tag = sp.elem + to_string(isotopes[i].first);
            graphs_to_save[tag + "_final_R"] = std::move(temp_final_graphs_R[i]);
        }

        // --- STEP 5 & PLOTTING ---
        cout << "Step 5: Smoothing final isotope R fluxes..." << endl;
        for (int i = 0; i < n_isotopes; ++i) {
            const string tag = sp.elem + to_string(isotopes[i].first);
            string graph_name = tag + "_final_R";
            string spline_name = tag + "_spline_R";
            auto* final_graph = graphs_to_save[graph_name].get();
            functions_to_save[spline_name] = doSplineFit(final_graph, spline_name, 0.8, 3300, colors[i % 6]);
        }

        auto cR = new TCanvas(("cR_" + sp.elem).c_str(), (sp.elem + " R").c_str(), 900, 700);
        cR->SetLogx(); cR->SetLogy();
        vector<TGraphAsymmErrors*> graphs_for_r_range;
        graphs_for_r_range.push_back(gTot_R);
        for(const auto& iso : isotopes) { graphs_for_r_range.push_back(graphs_to_save[sp.elem + to_string(iso.first) + "_final_R"].get()); }
        auto [yminR, ymaxR] = findYRangePos(graphs_for_r_range);

        auto frameR = (TGraphAsymmErrors*)gTot_R->Clone((sp.elem + "_frame_R").c_str());
        frameR->SetMarkerSize(0); frameR->SetLineColor(0);
        frameR->GetXaxis()->SetTitle("Rigidity R [GV]");
        frameR->GetYaxis()->SetTitle("Flux [m^{-2} sr^{-1} s^{-1} GV^{-1}]");
        frameR->GetYaxis()->SetRangeUser(yminR, ymaxR);
        frameR->GetXaxis()->SetRangeUser(0.8, 3500);
        frameR->Draw("AP");

        auto legR = new TLegend(0.64, 0.60, 0.84, 0.88);
        legR->SetBorderSize(0); legR->SetFillStyle(0);
        
        gTot_R->SetMarkerColor(kBlack);
        gTot_R->SetLineColor(kBlack);
        gTot_R->SetMarkerStyle(20);
        gTot_R->Draw("P SAME");
        legR->AddEntry(gTot_R, "Total Flux (Data)", "p");

        for (int i = 0; i < n_isotopes; ++i) {
            const string tag = sp.elem + to_string(isotopes[i].first);
            auto* g = graphs_to_save[tag + "_final_R"].get();
            auto* f = functions_to_save.count(tag + "_spline_R") ? functions_to_save[tag + "_spline_R"] : nullptr;
            g->SetMarkerStyle(24);
            g->SetMarkerColor(colors[i % 6]);
            g->SetLineColor(colors[i % 6]);
            g->Draw("P SAME");
            if (f) f->Draw("SAME");
            legR->AddEntry(g, (tag).c_str(), "p");
        }
        legR->Draw();
        cR->Print(outPdf.c_str());
        delete legR; delete frameR; delete cR;
        fin.Close();
    }

    c_dummy.Print((outPdf + "]").c_str());

    TFile fout(outRoot.c_str(), "RECREATE");
    if (fout.IsZombie()) { cerr << "[ERROR] Cannot create output file: " << outRoot << endl; return; }
    
    cout << "\n[INFO] Saving all graphs and functions to " << outRoot << "..." << endl;
    for (auto& kv : graphs_to_save) {
        if(kv.second) kv.second->Write();
    }
    for (auto& kv : functions_to_save) {
        if (kv.second) kv.second->Write();
    }

    fout.Close();
    cout << "[INFO] Done." << endl;
}