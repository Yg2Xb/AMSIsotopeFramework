#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TError.h>
#include <TLine.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <limits>

// 假设 SplineFit 在此定义
#include "../Tool.h"

using namespace std;
using namespace AMS_Iso;

// =================================================================================
// 简单结构体处理误差传播
// =================================================================================
struct ValueWithError {
    double val = 0.0;
    double err = 0.0;

    ValueWithError() = default;
    ValueWithError(double v, double e) : val(v), err(e) {}

    ValueWithError operator+(const ValueWithError& o) const {
        return {val + o.val, std::hypot(err, o.err)};
    }
    ValueWithError operator-(const ValueWithError& o) const {
        return {val - o.val, std::hypot(err, o.err)};
    }
    ValueWithError operator*(const ValueWithError& o) const {
        double new_val = val * o.val;
        if (new_val == 0.0) return {0.0, 0.0};
        double rel_err_sq = (val != 0 ? (err/val)*(err/val) : 0) + (o.val != 0 ? (o.err/o.val)*(o.err/o.val) : 0);
        return {new_val, std::abs(new_val) * std::sqrt(rel_err_sq)};
    }
    ValueWithError operator/(const ValueWithError& o) const {
        if (o.val == 0.0) return {0.0, 0.0};
        double new_val = val / o.val;
        if (new_val == 0.0) return {0.0, 0.0};
        double rel_err_sq = (val != 0 ? (err/val)*(err/val) : 0) + (o.val != 0 ? (o.err/o.val)*(o.err/o.val) : 0);
        return {new_val, std::abs(new_val) * std::sqrt(rel_err_sq)};
    }
};

// =================================================================================
// 全局配置 & 节点定义
// =================================================================================
const vector<string> CHAINS = {"UnbiasedL1Inner", "L1Inner"};
const string FLUX_FILE_PATH = "/eos/user/z/zixuan/Isotope/FluxSmooth/FluxSmooth.root";
const string ACC_FILE_DIR = "/eos/user/z/zixuan/Isotope/Add/";
const string OUTPUT_DIR = "/eos/user/z/zixuan/Isotope/BkgEst/";

double xpoints2_arr[] = {0.2, 1.5, 5.62, 10.00, 27};
const int n_xpoints2 = sizeof(xpoints2_arr)/sizeof(double);

double xpoints_global_arr[] = {0.2, 0.77, 1.5, 3.00, 5.62, 10.00, 27};
const int n_xpoints_global = sizeof(xpoints_global_arr)/sizeof(double);

double xpoints_tof_arr[] = {0.2, 0.77, 1.5};
const int n_xpoints_tof = sizeof(xpoints_tof_arr)/sizeof(double);

double xpoints_naf_arr[] = {0.61, 2.50, 6.10};
const int n_xpoints_naf = sizeof(xpoints_naf_arr)/sizeof(double);

double xpoints_agl_arr[] = {2.7, 5.62, 10.00, 27};
const int n_xpoints_agl = sizeof(xpoints_agl_arr)/sizeof(double);

// 计算 Ratio 和 Fit 的能区范围 (BinCenter)
const map<string, pair<double, double>> DET_CALC_RANGES = {
    {"TOF", {0.27, 1.5}},
    {"NaF", {0.61, 6.1}},
    {"AGL", {2.7, 23.0}}
};

// 最终拼接能区 (Epsilon Stitching)
const double STITCH_TOF_NAF = 1.2;
const double STITCH_NAF_AGL = 3.2;

const bool DO_REBIN = true;
const vector<string> ALL_PARTICLES = {"Be7", "Be9", "Be10", "B10", "B11", "C12", "N14", "N15","O16"};
const vector<string> ALL_SECONDARY_ISOTOPES = {"Be7", "Be9", "Be10"};
const vector<string> DETECTORS = {"TOF", "NaF", "AGL"};

// 绘图颜色
const int COLOR_TOF = kRed;
const int COLOR_NAF = kBlue;
const int COLOR_AGL = kGreen+2;
const int COLOR_COMBINED = kBlack;

// =================================================================================
// 辅助函数
// =================================================================================

int getMassNumber(const string& particleName) {
    size_t first_digit = particleName.find_first_of("0123456789");
    if (first_digit == string::npos) return 0;
    return stoi(particleName.substr(first_digit));
}

// 恢复正常的 Style 设置
void setHistStyle(TH1* h, int color, int markerStyle=20) {
    if(!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerSize(1.0);
    h->SetStats(0);
    h->GetYaxis()->SetTitleOffset(1.0);
}

// 设置 Legend 样式
void setLegendStyle(TLegend* leg) {
    if(!leg) return;
    leg->SetFillStyle(0);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.04); 
}

// 计算 Acceptance
map<int, ValueWithError> calculateAcceptance(const string& accFilePath, const string& eventHistName) {
    map<int, ValueWithError> acceptance_map;

    TFile* file = TFile::Open(accFilePath.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open " << accFilePath << endl;
        if(file) delete file;
        return acceptance_map;
    }

    TH1* h_events_raw = (TH1*)file->Get(eventHistName.c_str());
    TH1* h_mc_flux_raw = (TH1*)file->Get("MC_FLUX_H3");

    if (!h_events_raw || !h_mc_flux_raw) {
        cerr << "Error: Hists not found in " << accFilePath << " (" << eventHistName << ")" << endl;
        return acceptance_map;
    }

    TH1* h_events = (TH1*)h_events_raw->Clone(); 
    h_events->SetDirectory(0);
    
    TH1* h_mc_flux = (TH1*)h_mc_flux_raw->Clone(); 
    h_mc_flux->SetDirectory(0);


    if (DO_REBIN) {
        h_events->Rebin(2);
        h_mc_flux->Rebin(2);
    }

    double scale = TMath::Power(3.9, 2) * TMath::Pi();

    for (int i = 1; i <= h_events->GetNbinsX(); ++i) {
        double k = h_events->GetBinContent(i);
        double N = h_mc_flux->GetBinContent(i);

        if (N > 0) {
            double eff = k / N;
            double eff_err = (k > 0 && k < N) ? sqrt(eff * (1.0 - eff) / N) : (k==0 ? 1.0/N : 0.0);
            acceptance_map[i] = ValueWithError(eff * scale, eff_err * scale);
        } else {
            acceptance_map[i] = ValueWithError(0.0, 0.0);
        }
    }

    delete h_events;
    delete h_mc_flux;

    return acceptance_map;
}

// 封装 SplineFit
TF1* smoothRatio(TH1* h_ratio, const string& fitName, bool isGlobal) {
    if(!h_ratio) return nullptr;
    
    double* knots = xpoints_global_arr;
    int n_knots = n_xpoints_global;

    if (isGlobal) {
        knots = xpoints_global_arr;
        n_knots = n_xpoints_global;
        if (fitName.find("C12_to_Be7") != string::npos || fitName.find("C12_to_B11") != string::npos || 
            fitName.find("C12_to_Be10") != string::npos || fitName.find("N15_to_Be7") != string::npos || 
            fitName.find("O16_to_Be10") != string::npos || fitName.find("O16_to_B11") != string::npos) {
            knots = xpoints2_arr;
            n_knots = n_xpoints2;
        }
    } else {
        knots = xpoints_agl_arr;
        n_knots = n_xpoints_agl;
        if (fitName.find("TOF") != string::npos) {
            knots = xpoints_tof_arr;
            n_knots = n_xpoints_tof;
            if(fitName.find("O16_to_Be7") != string::npos)
            {
                knots = new double[5]{0.2, 0.4, 0.7, 1.1, 1.5};
                n_knots = 5;
            }
        }
        else if (fitName.find("NaF") != string::npos) {
            knots = xpoints_naf_arr;
            n_knots = n_xpoints_naf;
        }
    }

    TF1* fit = nullptr;
    try {
        fit = SplineFit(h_ratio, knots, n_knots, 0x38, "b1e1", fitName.c_str(), 0.27, 23.5); // 调整 fitting range
    } catch (...) {
        return nullptr;
    }
    return fit;
}

// =================================================================================
// 主程序
// =================================================================================
void est_frag() {
    auto GetNonZeroMin = [](TH1* h) -> double {
        double min_val = std::numeric_limits<double>::max();
        bool found = false;
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double c = h->GetBinContent(i);
            if (c > 0 && c < min_val) {
                min_val = c;
                found = true;
            }
        }
        return found ? min_val : 1e-5;
    };
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0); 
    gErrorIgnoreLevel = kWarning;

    TFile* fluxFile = TFile::Open(FLUX_FILE_PATH.c_str());
    if (!fluxFile || fluxFile->IsZombie()) {
        cerr << "FATAL: Flux file error." << endl;
        return;
    }

    map<string, TF1*> map_fluxes;
    for (const auto& p : ALL_PARTICLES) {
        TF1* f = (TF1*)fluxFile->Get((p + "_spline_Ek").c_str());
        if(!f) {
            cerr << "FATAL: Flux for " << p << " not found." << endl;
            return;
        }
        map_fluxes[p] = f;
    }

    // 用于最后 Summary PNG (Chain -> Sec -> Hists)
    map<string, map<string, vector<TH1D*>>> all_summary_hists;
    map<string, map<string, vector<string>>> all_summary_labels;

    for (const auto& secondary : ALL_SECONDARY_ISOTOPES) {
        int sec_Z = (secondary.find("Be") == 0) ? 4 : 5;
        int sec_A = getMassNumber(secondary);

        vector<string> primaries;
        for (const auto& p : ALL_PARTICLES) {
            if (getMassNumber(p) > sec_A) {
                if (!(p.find("Be") == 0 && secondary.find("Be") == 0)) {
                    primaries.push_back(p);
                }
            }
        }

        for (const auto& primary : primaries) {
            if (primary == "Be7") continue; 

            cout << "\n>>> Processing " << primary << " -> " << secondary << "..." << endl;
            
            // Output Root File (One per Primary->Secondary, contains all chains)
            string root_out = OUTPUT_DIR + "Epsilon_" + primary + "_to_" + secondary + ".root";
            TFile* fout = new TFile(root_out.c_str(), "RECREATE");

            // Output PDF (One per Primary->Secondary, contains all chains pages)
            string pdf_path = OUTPUT_DIR + "rew_results_" + primary + "_to_" + secondary + ".pdf";
            TCanvas* c1 = new TCanvas("c1", "Canvas", 900, 700);
            c1->Print((pdf_path + "[").c_str());

            // ---------------- CHAIN LOOP ----------------
            for (const auto& chain : CHAINS) {
                cout << "   > Chain: " << chain << endl;

                string channel_label = primary + "->" + secondary;
                
                TF1* flux_X = map_fluxes[primary];
                TF1* flux_Y = map_fluxes[secondary];

                // 1. Plot Fluxes
                c1->Clear(); c1->SetLogy(1); c1->SetLogx(1);
                flux_X->SetLineColor(kRed); 
                flux_X->SetTitle(Form("Fluxes %s", channel_label.c_str()));
                flux_X->GetYaxis()->SetTitle("Flux");
                flux_X->GetXaxis()->SetRangeUser(0.27, 100.0);
                flux_X->Draw(); 
                flux_Y->SetLineColor(kBlue); 
                flux_Y->Draw("SAME");
                TLegend* leg1 = new TLegend(0.78, 0.78, 0.88, 0.88);
                setLegendStyle(leg1);
                leg1->AddEntry(flux_X, primary.c_str(), "l");
                leg1->AddEntry(flux_Y, secondary.c_str(), "l");
                leg1->Draw();
                c1->Print(pdf_path.c_str());
                delete leg1;

                // 2. Flux Ratio
                c1->Clear(); c1->SetLogy(0); c1->SetLogx(1);
                string hname_fr = Form("h_fr_%s_%s_%s", chain.c_str(), primary.c_str(), secondary.c_str());
                TH1D* h_fr = new TH1D(hname_fr.c_str(), Form("Flux Ratio %s;E_{k}/n;Flux Ratio", channel_label.c_str()), 200, 0.27, 100);
                h_fr->SetStats(0);
                for(int b=1; b<=h_fr->GetNbinsX(); ++b) {
                    double x = h_fr->GetBinCenter(b);
                    double vy = flux_Y->Eval(x);
                    if(vy > 0) h_fr->SetBinContent(b, flux_X->Eval(x)/vy);
                }
                setHistStyle(h_fr, kBlack, 1);
                h_fr->Draw("L");
                c1->Print(pdf_path.c_str());
                
                // SAVE FLUX RATIO
                fout->cd(); h_fr->Write();

                auto calc_flux_ratio = [&](double ek) {
                    double vx = flux_X->Eval(ek);
                    double vy = flux_Y->Eval(ek);
                    if(vy == 0) return ValueWithError(0,0);
                    ValueWithError wx(vx, vx*0.01);
                    ValueWithError wy(vy, vy*0.01);
                    return wx/wy;
                };

                string suffix = (secondary.find("Be") == 0) ? "_rew_frag4_withBkg_bkgest_eff.root" : "_rew_frag5.root";
                string acc_file_Y = ACC_FILE_DIR + secondary + suffix;
                string acc_file_XtoY = ACC_FILE_DIR + primary + suffix;

                // Lambda for hist name with chain
                auto getHName = [&](const string& d) {
                    return string(Form("%s_BKG_H3b_%s_Z%d_Mass%d", chain.c_str(), d.c_str(), sec_Z, sec_A));
                };

                // Template histogram
                TFile* ftmp = TFile::Open(acc_file_Y.c_str());
                TH1* h_tpl_raw = (TH1*)ftmp->Get(getHName("TOF").c_str()); 
                TH1D* h_template = (TH1D*)h_tpl_raw->Clone("h_template");
                h_template->SetDirectory(0);
                if(DO_REBIN) h_template->Rebin(2);

                vector<TH1D*> vec_det_eps; 
                vector<string> vec_det_names;

                // ---------------- DETECTOR LOOP ----------------
                // 处理独立探测器：TOF, NaF, AGL
                for(const auto& det : DETECTORS) {
                    auto calc_range = DET_CALC_RANGES.at(det);
                    int color = (det=="TOF")?COLOR_TOF : (det=="NaF"?COLOR_NAF : COLOR_AGL);

                    map<int, ValueWithError> ay = calculateAcceptance(acc_file_Y, getHName(det));
                    map<int, ValueWithError> axy = calculateAcceptance(acc_file_XtoY, getHName(det));

                    // Acc Hists
                    string suffix_name = "_" + chain + "_" + det + "_" + primary + "_" + secondary;
                    TH1D* h_ay_det = (TH1D*)h_template->Clone(("h_ay" + suffix_name).c_str());
                    h_ay_det->Reset(); h_ay_det->SetTitle(Form("%s Acc (%s)", det.c_str(), chain.c_str()));
                    
                    TH1D* h_axy_det = (TH1D*)h_template->Clone(("h_axy" + suffix_name).c_str());
                    h_axy_det->Reset(); h_axy_det->SetTitle(Form("%s Frag Acc (%s)", det.c_str(), chain.c_str()));

                    TH1D* h_ar_det = (TH1D*)h_template->Clone(("h_ar" + suffix_name).c_str());
                    h_ar_det->Reset(); h_ar_det->SetTitle(Form("%s Acc Ratio (%s)", det.c_str(), chain.c_str()));

                    // Fill Data
                    for(int i=1; i<=h_ar_det->GetNbinsX(); ++i) {
                        double ek = h_ar_det->GetBinCenter(i);
                        // 只在指定计算能区内填充
                        if(ek >= calc_range.first && ek <= calc_range.second) {
                            h_ay_det->SetBinContent(i, ay[i].val); h_ay_det->SetBinError(i, ay[i].err);
                            h_axy_det->SetBinContent(i, axy[i].val); h_axy_det->SetBinError(i, axy[i].err);
                            ValueWithError r = axy[i] / ay[i];
                            h_ar_det->SetBinContent(i, r.val); h_ar_det->SetBinError(i, r.err);
                        }
                    }

                    // Save Acc and Ratio
                    fout->cd();
                    h_ay_det->Write(); h_axy_det->Write(); h_ar_det->Write();

                    // Draw Single Acc
                    c1->Clear(); c1->SetLogx(1); c1->SetLogy(1);
                    h_ay_det->GetXaxis()->SetRangeUser(calc_range.first*0.9, calc_range.second*1.1);
                    h_ay_det->GetYaxis()->SetTitle("Acceptance [m^{2} sr]");
                    setHistStyle(h_ay_det, kBlue);
                    h_ay_det->GetYaxis()->SetRangeUser(GetNonZeroMin(h_axy_det)*0.1, h_ay_det->GetMaximum()*10);
                    h_ay_det->Draw("P");
                    setHistStyle(h_axy_det, kRed);
                    h_axy_det->Draw("SAME P");
                    TLegend* legAcc = new TLegend(0.75, 0.78, 0.9, 0.88);
                    setLegendStyle(legAcc);
                    legAcc->AddEntry(h_ay_det, "Acc", "p");
                    legAcc->AddEntry(h_axy_det, "Frag Acc", "p");
                    legAcc->Draw();
                    c1->Print(pdf_path.c_str());
                    delete legAcc; delete h_ay_det; delete h_axy_det;

                    // Fit
                    string fit_name = "fit" + suffix_name;
                    TF1* fit_det = smoothRatio(h_ar_det, fit_name, false);
                    
                    // Save Fit
                    if(fit_det) { fout->cd(); fit_det->Write(); }

                    // Draw Ratio & Fit
                    c1->Clear(); c1->SetLogx(1); c1->SetLogy(0);
                    h_ar_det->GetXaxis()->SetRangeUser(calc_range.first*0.9, calc_range.second*1.1);
                    h_ar_det->GetYaxis()->SetTitle("Acceptance Ratio");
                    setHistStyle(h_ar_det, color);
                    h_ar_det->Draw("P");
                    if(fit_det) {
                        fit_det->SetLineColor(kBlack);
                        fit_det->Draw("SAME");
                    }
                    c1->Print(pdf_path.c_str());

                    // Calc Epsilon (Det)
                    TH1D* h_eps_det = (TH1D*)h_template->Clone(("h_eps" + suffix_name).c_str());
                    h_eps_det->Reset();
                    h_eps_det->SetTitle(Form("%s Bkg Frac (%s)", det.c_str(), chain.c_str()));
                    
                    if(fit_det) {
                        for(int i=1; i<=h_eps_det->GetNbinsX(); ++i) {
                            double ek = h_eps_det->GetBinCenter(i);
                            if(ek >= calc_range.first && ek <= calc_range.second) {
                                ValueWithError fr = calc_flux_ratio(ek);
                                ValueWithError ar(fit_det->Eval(ek), h_ar_det->GetBinError(i));
                                ValueWithError eps = fr * ar;
                                h_eps_det->SetBinContent(i, eps.val);
                                h_eps_det->SetBinError(i, eps.err);
                            }
                        }
                    }
                    setHistStyle(h_eps_det, color);
                    h_eps_det->GetYaxis()->SetTitle("Background Fraction");
                    
                    // Save Epsilon Det
                    fout->cd(); h_eps_det->Write();
                    
                    vec_det_eps.push_back(h_eps_det);
                    vec_det_names.push_back(det);
                    delete h_ar_det;
                } // End Detector Loop

                // ---------------- COMBINED (STITCHING) ----------------
                // 拼接逻辑：TOF < 1.2, NaF [1.2, 3.2), AGL >= 3.2
                TH1D* h_eps_comb = (TH1D*)h_template->Clone(Form("h_eps_comb_%s_%s_%s", chain.c_str(), primary.c_str(), secondary.c_str()));
                h_eps_comb->Reset();
                h_eps_comb->SetTitle(channel_label.c_str());
                h_eps_comb->GetYaxis()->SetTitle("Background Fraction");

                TH1D* h_tof = vec_det_eps[0]; // TOF
                TH1D* h_naf = vec_det_eps[1]; // NaF
                TH1D* h_agl = vec_det_eps[2]; // AGL

                for(int i=1; i<=h_eps_comb->GetNbinsX(); ++i) {
                    double ek = h_eps_comb->GetBinCenter(i);
                    double val = 0, err = 0;
                    
                    if (ek < STITCH_TOF_NAF) {
                        val = h_tof->GetBinContent(i); err = h_tof->GetBinError(i);
                    } else if (ek >= STITCH_TOF_NAF && ek < STITCH_NAF_AGL) {
                        val = h_naf->GetBinContent(i); err = h_naf->GetBinError(i);
                    } else if (ek >= STITCH_NAF_AGL) {
                        val = h_agl->GetBinContent(i); err = h_agl->GetBinError(i);
                    }

                    if(val > 0) {
                        h_eps_comb->SetBinContent(i, val);
                        h_eps_comb->SetBinError(i, err);
                    }
                }

                // Save Combined
                fout->cd(); h_eps_comb->Write();
                
                // Add to Summary Collection
                all_summary_hists[chain][secondary].push_back((TH1D*)h_eps_comb->Clone());
                all_summary_labels[chain][secondary].push_back(channel_label);

                // Draw Summary for this Chain
                c1->Clear(); c1->SetLogx(1);
                TH1D* h_draw = (TH1D*)h_eps_comb->Clone();
                h_draw->GetXaxis()->SetRangeUser(0.27, 21.5);
                setHistStyle(h_draw, kBlack, 20);
                h_draw->Draw("P");
                TLegend* legSum = new TLegend(0.7, 0.68, 0.9, 0.88);
                setLegendStyle(legSum);
                legSum->AddEntry(h_draw, "Combined", "p");
                for(size_t i=0; i<vec_det_eps.size(); ++i) {
                    vec_det_eps[i]->Draw("P SAME");
                    legSum->AddEntry(vec_det_eps[i], vec_det_names[i].c_str(), "p");
                }
                legSum->Draw();
                c1->Print(pdf_path.c_str());

                // Cleanup Chain Loop
                delete h_fr; delete h_template; delete h_draw; delete legSum;
                for(auto h : vec_det_eps) delete h;

            } // End Chain Loop

            c1->Print((pdf_path + "]").c_str());
            delete c1;

        } // End Primary Loop
    } // End Secondary Loop

    cout << "\n--- Generating Summary PNGs ---\n" << endl;
    
    TFile* resFile = new TFile((OUTPUT_DIR + "Epsilon_Results_Summary.root").c_str(), "RECREATE");

    // Loop Chains to draw summary separately
    for (const auto& chain : CHAINS) {
        for(auto const& [sec, hists] : all_summary_hists[chain]) {
            if(hists.empty()) continue;

            TH1D* h_sum = (TH1D*)hists[0]->Clone(Form("h_sum_%s_%s", chain.c_str(), sec.c_str()));
            h_sum->Reset();
            h_sum->SetTitle(Form("Total %s to %s", chain.c_str(), sec.c_str()));
            
            for(auto h : hists) h_sum->Add(h);

            TCanvas* c2 = new TCanvas("c2", "Sum", 1200, 600);
            c2->SetGrid();
            
            double ymax = h_sum->GetMaximum() * 1.25;
            if(ymax <= 0) ymax = 1.0;

            TH1* frame = c2->DrawFrame(0.27, 0, 21.5, ymax);
            frame->SetTitle(Form("Fragmentation to %s (%s);E_{k}/n;Background Fraction", sec.c_str(), chain.c_str()));

            TLegend* leg = new TLegend(0.68, 0.5, 0.88, 0.88);
            setLegendStyle(leg);

            const vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta, kOrange-3, kCyan+1, kPink+7, kSpring-5};

            for(size_t i=0; i<hists.size(); ++i) {
                setHistStyle(hists[i], colors[i % colors.size()]);
                hists[i]->GetXaxis()->SetRangeUser(0.27, 21.5);
                hists[i]->Draw("P SAME");
                resFile->cd(); hists[i]->Write();
                leg->AddEntry(hists[i], all_summary_labels[chain][sec][i].c_str(), "p");
            }

            setHistStyle(h_sum, kBlack);
            h_sum->Draw("P SAME");
            resFile->cd(); h_sum->Write();
            leg->AddEntry(h_sum, "Total", "p");
            leg->Draw();
            c2->SetLeftMargin(0.15); c2->SetRightMargin(0.05);c2->SetBottomMargin(0.15);
            c2->SaveAs((OUTPUT_DIR + "rew_frag_epsilon_" + chain + "_" + sec + ".png").c_str());
            
            delete c2; delete leg; delete h_sum;
            for(auto h : hists) delete h;
        }
    }

    resFile->Close(); delete resFile;
    fluxFile->Close(); delete fluxFile;

    cout << "Done." << endl;
}