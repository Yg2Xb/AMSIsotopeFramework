#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>

#include "../Tool.h" 

using namespace std;
using namespace AMS_Iso;

struct FitResult {
    double mean = 0.0;
    double mean_err = 0.0;
    double sigma = 0.0;
    double sigma_err = 0.0;
    double chi2 = 0.0;
    double ndf = 0.0;
    double LR = 0.0;
    double LR_err = 0.0;
    double RR = 0.0;
    double RR_err = 0.0;
    bool IsValid() const { return mean_err > 0 || sigma_err > 0; }
};

struct HistInfo {
    string suffix;
    string x_label;
    string out_suffix;
    string desc;
};

struct Config {
    const vector<double> rig_edges = {30.0, 50.0, 80.0, 120.0, 160.0, 240.0};
    const vector<string> rig_labels = {"30-50 GV", "50-80 GV", "80-120 GV", "120-160 GV", "160-240 GV"};
    const vector<double> charge_edges = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
    const vector<string> charge_labels = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
    const vector<string> suffixes_h5 = {"ID_H5a", "ID_H5b"};
    const vector<string> suffixes_h4 = {"ID_H4a", "ID_H4b"};
    
    const vector<int> colors_iss = {kYellow+2, kMagenta, kBlack, kRed, kBlue, kGreen + 2, kOrange + 1};
    const vector<int> colors_mc = {kYellow+2, kMagenta, kBlack, kRed, kBlue, kGreen + 2, kOrange + 1, kViolet, kCyan};
    const vector<int> colors_rig = {kBlack, kRed, kBlue, kGreen+2, kMagenta};

    string path_base = "/eos/user/z/zixuan/Isotope/Add/";
    string file_iss = "Be_frag4.root";
    string dir_out_iss = "/eos/user/z/zixuan/Isotope/Beta/ISS/";
    string dir_out_mc = "/eos/user/z/zixuan/Isotope/Beta/MC/";
    
    vector<string> files_mc = {
        "B10_rew_frag4.root", "B11_rew_frag4.root", "Be10_rew_frag4.root", "Be7_rew_frag4.root",
        "Be9_rew_frag4.root", "C12_rew_frag4.root", "N15_rew_frag4.root", "O16_rew_frag4.root"
    };

    int n_rig() const { return rig_edges.size() - 1; }
    int n_chg() const { return charge_edges.size() - 1; }
};

class BetaAnalysis {
public:
    BetaAnalysis() {
        gStyle->SetOptFit(0);
        gStyle->SetErrorX(0);
    }

    void Run() {
        ProcessISS();
        ProcessMC();
    }

private:
    Config cfg;
    TFile* f_out = nullptr;
    string curr_dir;

    void SetupOutput(string dir, string fname) {
        curr_dir = dir;
        if (gSystem->AccessPathName(dir.c_str())) gSystem->mkdir(dir.c_str(), kTRUE);
        if (f_out) { f_out->Close(); delete f_out; }
        f_out = TFile::Open((dir + fname).c_str(), "RECREATE");
    }

    HistInfo GetInfo(string s) {
        if (s == "ID_H5a") return {s, "Rigidity [GV]", "Rigidity", "NaF-Tracker #Delta(1/#beta)"};
        if (s == "ID_H5b") return {s, "Rigidity [GV]", "Rigidity", "AGL-Tracker #Delta(1/#beta)"};
        if (s == "ID_H4a") return {s, "Rigidity [GV]", "Rigidity", "NaF 1/#beta (Rig > 80GV)"};
        if (s == "ID_H4b") return {s, "Rigidity [GV]", "Rigidity", "AGL 1/#beta (Rig > 150GV)"};
        return {s, "Y", "Unk", "Unknown"};
    }

    string ParseNuclide(string s) {
        size_t i = s.find_last_of('/');
        if (i != string::npos) s = s.substr(i + 1);
        size_t j = s.find("_rew");
        if (j != string::npos) return s.substr(0, j);
        return s.substr(0, s.find_last_of('.'));
    }

    double GetZ(string n) {
        if (n.empty()) return 0;
        if (n[0] == 'H') return 2.0;
        if (n[0] == 'L') return 3.0;
        if (n[0] == 'B') return (n.size() > 1 && (n[1] == 'e' || n[1] == 'E')) ? 4.0 : 5.0;
        if (n[0] == 'C') return 6.0;
        if (n[0] == 'N') return 7.0;
        if (n[0] == 'O') return 8.0;
        return 0.0;
    }

    FitResult Fit(TH1* h, TCanvas* c, string title, double center, string rig_lbl) {
        FitResult r;
        if (!h || h->GetEntries() < 10) return r;
        h->Sumw2();

        int bin = h->GetXaxis()->FindFixBin(center);
        if (bin < 1 || bin > h->GetNbinsX()) bin = h->GetMaximumBin();
        
        double integ_tgt = h->Integral() * 0.8;
        double integ_curr = h->GetBinContent(bin);
        int l = bin - 1, u = bin + 1;
        double x1 = h->GetXaxis()->GetBinLowEdge(bin);
        double x2 = h->GetXaxis()->GetBinUpEdge(bin);

        while (integ_curr < integ_tgt) {
            bool exp = false;
            double cl = (l >= 1) ? h->GetBinContent(l) : 0;
            double cu = (u <= h->GetNbinsX()) ? h->GetBinContent(u) : 0;
            if (l >= 1 && u <= h->GetNbinsX()) { integ_curr += cl + cu; x1 = h->GetXaxis()->GetBinLowEdge(l--); x2 = h->GetXaxis()->GetBinUpEdge(u++); exp = true; }
            else if (l >= 1) { integ_curr += cl; x1 = h->GetXaxis()->GetBinLowEdge(l--); exp = true; }
            else if (u <= h->GetNbinsX()) { integ_curr += cu; x2 = h->GetXaxis()->GetBinUpEdge(u++); exp = true; }
            if (!exp) break;
        }
        double fc = (x1 + x2) / 2.0;
        double hw = max(abs(x2 - fc), abs(x1 - fc));
        x1 = fc - hw; x2 = fc + hw;

        if (x1 >= x2) return r;

        TF1* f1 = new TF1("f1", "gaus", x1, x2);
        f1->SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
        if (h->Fit(f1, "QRS") != 0) { delete f1; return r; }
        double m0 = f1->GetParameter(1), s0 = f1->GetParameter(2);
        delete f1;

        double nsig = (rig_lbl == "30-50 GV" && title.find("Be") != string::npos) ? 2.5 : 4.0;
        double fit_min = m0 - nsig * abs(s0);
        double fit_max = m0 + nsig * abs(s0);

        h->SetTitle(title.c_str());
        h->GetYaxis()->SetTitle("Events");
        double buff = 0.6 * (fit_max - fit_min);
        h->GetXaxis()->SetRangeUser(fit_min - buff, fit_max + buff);

        vector<vector<double>> res = DoGausPlusAsymGausFit(h, fit_min, fit_max, c, true);
        if (res.size() == 2 && res[1][4] > 0) {
            r.mean = res[0][0]; r.mean_err = res[1][0];
            r.sigma = res[0][1]; r.sigma_err = res[1][1];
            r.LR = res[0][2]; r.LR_err = res[1][2];
            r.RR = res[0][3]; r.RR_err = res[1][3];
            r.chi2 = res[0][4]; r.ndf = res[1][4];
            c->cd();
            TLine l1(fit_min, 0, fit_min, h->GetMaximum()); l1.SetLineStyle(2); l1.SetLineColor(kRed); l1.DrawClone();
            TLine l2(fit_max, 0, fit_max, h->GetMaximum()); l2.SetLineStyle(2); l2.SetLineColor(kRed); l2.DrawClone();
            c->Update();
        } else {
            c->cd(); h->Draw("hist");
            TLatex t; t.SetNDC(); t.SetTextSize(0.04); t.DrawLatex(0.5, 0.8, "Fit Failed"); c->Update();
        }
        return r;
    }

    void Draw(const vector<TGraphErrors*>& gs, const vector<string>& lg, string xt, string yt, string t, string suf) {
        if (gs.empty()) return;
        string cn = "c_" + t + suf;
        replace(cn.begin(), cn.end(), ' ', '_'); replace(cn.begin(), cn.end(), '#', '_');
        TCanvas* c = new TCanvas(cn.c_str(), t.c_str(), 700, 500);
        c->SetGrid();

        double xmin = gs[0]->GetXaxis()->GetXmin(), xmax = gs[0]->GetXaxis()->GetXmax();
        if (suf.find("Charge") != string::npos) { xmin = 1.0; xmax = 9.0; }
        else { xmin = cfg.rig_edges.front(); xmax = cfg.rig_edges.back() * 1.5; }

        double ymin = 1e9, ymax = -1e9;
        bool f = false;
        for (auto g : gs) {
            for (int i = 0; i < g->GetN(); ++i) {
                double x, y; g->GetPoint(i, x, y);
                if (x >= xmin && x <= xmax) {
                    ymin = min(ymin, y - g->GetErrorY(i));
                    ymax = max(ymax, y + g->GetErrorY(i));
                    f = true;
                }
            }
        }
        if (!f) { ymin = -0.1; ymax = 0.1; }
        else { double r = ymax - ymin; if (r == 0) r = 0.01; ymin -= r * 0.1; ymax += r * 0.1; }

        TLegend* leg = new TLegend(0.8, 0.8, 0.99, 0.99);
        for (size_t i = 0; i < gs.size(); ++i) {
            gs[i]->Draw(i == 0 ? "APZ" : "PZ same");
            if (i == 0) {
                gs[i]->SetTitle((t + " vs " + xt).c_str());
                gs[i]->GetXaxis()->SetTitle(xt.c_str());
                gs[i]->GetYaxis()->SetTitle(yt.c_str());
                gs[i]->GetXaxis()->SetRangeUser(xmin, xmax);
                gs[i]->GetYaxis()->SetRangeUser(ymin, ymax);
            }
            if (i < lg.size()) leg->AddEntry(gs[i], lg[i].c_str(), "p");
            f_out->cd();
            string gn = "g_" + t + "_" + to_string(i);
            replace(gn.begin(), gn.end(), ' ', '_');
            gs[i]->SetName(gn.c_str()); gs[i]->Write();
        }
        leg->Draw();
        c->SaveAs((curr_dir + t + suf + ".png").c_str());
        delete c; delete leg;
    }

    void ProcessISS() {
        SetupOutput(cfg.dir_out_iss, "ISS_RICHBetaStudy.root");
        TFile* fin = TFile::Open((cfg.path_base + cfg.file_iss).c_str());
        if (!fin || fin->IsZombie()) return;
        cout << "Analyzing ISS..." << endl;

        for (string suf : cfg.suffixes_h5) {
            HistInfo info = GetInfo(suf);
            string base = info.out_suffix + suf.substr(3);
            string pdf = curr_dir + "FitResults_" + base + ".pdf";
            TCanvas* cf = new TCanvas("cf", "f", 800, 600); cf->Print((pdf + "[").c_str());

            f_out->cd();
            TH2F* h2m = new TH2F(("h2m_" + suf).c_str(), "Mean", cfg.n_rig(), &cfg.rig_edges[0], cfg.n_chg(), &cfg.charge_edges[0]);
            TH2F* h2s = new TH2F(("h2s_" + suf).c_str(), "Sigma", cfg.n_rig(), &cfg.rig_edges[0], cfg.n_chg(), &cfg.charge_edges[0]);

            for (size_t p = 0; p < cfg.charge_labels.size(); ++p) {
                string hn = "UnbiasedL1Inner_" + cfg.charge_labels[p] + "_" + suf;
                TH2F* h2 = (TH2F*)fin->Get(hn.c_str());
                if (!h2) continue;
                for (int r = 0; r < cfg.n_rig(); ++r) {
                    TH1D* h1 = h2->ProjectionX(Form("p_%s_%d", hn.c_str(), r), h2->GetYaxis()->FindFixBin(cfg.rig_edges[r]), h2->GetYaxis()->FindFixBin(cfg.rig_edges[r+1] - 1e-6));
                    if (h1->GetEntries() > 0) {
                        while (h1->GetMaximum() < 80 && h1->GetBinWidth(1) < 1e-6 && h1->GetNbinsX() > 10) h1->Rebin(2);
                        FitResult res = Fit(h1, cf, cfg.charge_labels[p] + " " + info.out_suffix + " " + cfg.rig_labels[r], 0.0, cfg.rig_labels[r]);
                        if (res.IsValid()) {
                            h2m->SetBinContent(r + 1, p + 1, res.mean); h2m->SetBinError(r + 1, p + 1, res.mean_err);
                            h2s->SetBinContent(r + 1, p + 1, res.sigma); h2s->SetBinError(r + 1, p + 1, res.sigma_err);
                            cf->Print(pdf.c_str());
                        }
                    }
                    delete h1;
                }
            }
            cf->Print((pdf + "]").c_str()); delete cf; h2m->Write(); h2s->Write();

            vector<TGraphErrors*> gmr, gsr;
            for (int p = 1; p <= cfg.n_chg(); ++p) {
                TGraphErrors *gm = new TGraphErrors(), *gs = new TGraphErrors();
                int c = cfg.colors_iss[p - 1];
                gm->SetMarkerStyle(20); gm->SetMarkerColor(c); gm->SetLineColor(c);
                gs->SetMarkerStyle(20); gs->SetMarkerColor(c); gs->SetLineColor(c);
                for (int r = 1; r <= cfg.n_rig(); ++r) {
                    if (h2m->GetBinError(r, p) > 0) {
                        double x = (cfg.rig_edges[r - 1] + cfg.rig_edges[r]) / 2.0;
                        gm->SetPoint(gm->GetN(), x, h2m->GetBinContent(r, p)); gm->SetPointError(gm->GetN() - 1, 0, h2m->GetBinError(r, p));
                        gs->SetPoint(gs->GetN(), x, h2s->GetBinContent(r, p)); gs->SetPointError(gs->GetN() - 1, 0, h2s->GetBinError(r, p));
                    }
                }
                gmr.push_back(gm); gsr.push_back(gs);
            }
            Draw(gmr, cfg.charge_labels, "Rigidity [GV]", "#mu", base + "_Mean", "_vs_Rig");
            Draw(gsr, cfg.charge_labels, "Rigidity [GV]", "#sigma", base + "_Sigma", "_vs_Rig");

            vector<TGraphErrors*> gmz, gsz;
            for (int r = 1; r <= cfg.n_rig(); ++r) {
                TGraphErrors *gm = new TGraphErrors(), *gs = new TGraphErrors();
                int c = cfg.colors_rig[(r - 1) % cfg.colors_rig.size()];
                gm->SetMarkerStyle(20); gm->SetMarkerColor(c); gm->SetLineColor(c);
                gs->SetMarkerStyle(20); gs->SetMarkerColor(c); gs->SetLineColor(c);
                for (int p = 1; p <= cfg.n_chg(); ++p) {
                    if (h2m->GetBinError(r, p) > 0) {
                        double z = (cfg.charge_edges[p - 1] + cfg.charge_edges[p]) / 2.0;
                        gm->SetPoint(gm->GetN(), z, h2m->GetBinContent(r, p)); gm->SetPointError(gm->GetN() - 1, 0, h2m->GetBinError(r, p));
                        gs->SetPoint(gs->GetN(), z, h2s->GetBinContent(r, p)); gs->SetPointError(gs->GetN() - 1, 0, h2s->GetBinError(r, p));
                    }
                }
                gmz.push_back(gm); gsz.push_back(gs);
            }
            Draw(gmz, cfg.rig_labels, "Charge (Z)", "#mu", base + "_Mean", "_vs_Charge");
            Draw(gsz, cfg.rig_labels, "Charge (Z)", "#sigma", base + "_Sigma", "_vs_Charge");
        }

        for (string suf : cfg.suffixes_h4) {
            HistInfo info = GetInfo(suf);
            string base = info.out_suffix + suf.substr(3);
            string pdf = curr_dir + "FitResults_" + base + ".pdf";
            TCanvas* cf = new TCanvas("cf", "f", 800, 600); cf->Print((pdf + "[").c_str());

            TGraphErrors *gm = new TGraphErrors(), *gs = new TGraphErrors();
            gm->SetMarkerStyle(20); gm->SetLineColor(kBlack); gs->SetMarkerStyle(20); gs->SetLineColor(kBlack);

            for (size_t p = 0; p < cfg.charge_labels.size(); ++p) {
                string hn = "UnbiasedL1Inner_" + cfg.charge_labels[p] + "_" + suf;
                TH1F* h = (TH1F*)fin->Get(hn.c_str());
                if (!h) continue;
                TH1F* hc = (TH1F*)h->Clone();
                while (hc->GetMaximum() < 80 && hc->GetNbinsX() > 10) hc->Rebin(2);
                FitResult res = Fit(hc, cf, cfg.charge_labels[p] + " " + info.desc, 1.0, "");
                if (res.IsValid()) {
                    double z = (cfg.charge_edges[p] + cfg.charge_edges[p+1]) / 2.0;
                    gm->SetPoint(gm->GetN(), z, res.mean); gm->SetPointError(gm->GetN() - 1, 0, res.mean_err);
                    gs->SetPoint(gs->GetN(), z, res.sigma); gs->SetPointError(gs->GetN() - 1, 0, res.sigma_err);
                    cf->Print(pdf.c_str());
                }
                delete hc;
            }
            cf->Print((pdf + "]").c_str()); delete cf;
            Draw({gm}, {"ISS"}, "Charge (Z)", "#mu", base + "_Mean", "_vs_Charge");
            Draw({gs}, {"ISS"}, "Charge (Z)", "#sigma", base + "_Sigma", "_vs_Charge");
        }
        fin->Close();
    }

    void ProcessMC() {
        SetupOutput(cfg.dir_out_mc, "MC_RICHBetaStudy.root");
        cout << "Analyzing MC..." << endl;
        
        map<string, TH2F*> mh5, sh5;
        map<string, FitResult> rh4;
        vector<string> nucs;

        for (string fn : cfg.files_mc) {
            string nuc = ParseNuclide(fn);
            TFile* fin = TFile::Open((cfg.path_base + fn).c_str());
            if (!fin || fin->IsZombie()) continue;
            nucs.push_back(nuc);

            for (string suf : cfg.suffixes_h5) {
                HistInfo info = GetInfo(suf);
                string base = info.out_suffix + suf.substr(3) + "_" + nuc;
                string pdf = curr_dir + "FitResults_" + base + ".pdf";
                TCanvas* cf = new TCanvas("cf", "f", 800, 600); cf->Print((pdf + "[").c_str());
                
                string hn = "UnbiasedL1Inner_Helium_" + suf;
                TH2F* h2 = (TH2F*)fin->Get(hn.c_str());
                if (h2) {
                    f_out->cd();
                    TH2F* h2m = new TH2F(("h2m_" + nuc + "_" + suf).c_str(), "Mean", cfg.n_rig(), &cfg.rig_edges[0], 1, 1.5, 2.5);
                    TH2F* h2s = new TH2F(("h2s_" + nuc + "_" + suf).c_str(), "Sigma", cfg.n_rig(), &cfg.rig_edges[0], 1, 1.5, 2.5);
                    for (int r = 0; r < cfg.n_rig(); ++r) {
                        TH1D* h1 = h2->ProjectionX(Form("p_%s_%s_%d", nuc.c_str(), suf.c_str(), r), h2->GetYaxis()->FindFixBin(cfg.rig_edges[r]), h2->GetYaxis()->FindFixBin(cfg.rig_edges[r+1] - 1e-6));
                        if (h1->GetEntries() > 0) {
                            while (h1->GetMaximum() < 60 && h1->GetNbinsX() > 10) h1->Rebin(2);
                            FitResult res = Fit(h1, cf, nuc + " " + cfg.rig_labels[r], 0.0, cfg.rig_labels[r]);
                            if (res.IsValid()) {
                                h2m->SetBinContent(r + 1, 1, res.mean); h2m->SetBinError(r + 1, 1, res.mean_err);
                                h2s->SetBinContent(r + 1, 1, res.sigma); h2s->SetBinError(r + 1, 1, res.sigma_err);
                                cf->Print(pdf.c_str());
                            }
                        }
                        delete h1;
                    }
                    mh5[nuc + "_" + suf] = h2m; sh5[nuc + "_" + suf] = h2s;
                }
                cf->Print((pdf + "]").c_str()); delete cf;
            }

            for (string suf : cfg.suffixes_h4) {
                HistInfo info = GetInfo(suf);
                string base = info.out_suffix + suf.substr(3) + "_" + nuc;
                string pdf = curr_dir + "FitResults_" + base + ".pdf";
                TCanvas* cf = new TCanvas("cf", "f", 800, 600); cf->Print((pdf + "[").c_str());

                string hn = "UnbiasedL1Inner_Helium_" + suf;
                TH1F* h = (TH1F*)fin->Get(hn.c_str());
                if (h) {
                    TH1F* hc = (TH1F*)h->Clone();
                    while (hc->GetMaximum() < 60 && hc->GetNbinsX() > 10) hc->Rebin(2);
                    FitResult res = Fit(hc, cf, nuc + " " + info.desc, 1.0, "");
                    if (res.IsValid()) { rh4[nuc + "_" + suf] = res; cf->Print(pdf.c_str()); }
                    delete hc;
                }
                cf->Print((pdf + "]").c_str()); delete cf;
            }
            fin->Close();
        }

        for (string suf : cfg.suffixes_h5) {
            string base = GetInfo(suf).out_suffix + suf.substr(3) + "_MC_Combined";
            vector<TGraphErrors*> gm, gs; vector<string> lg;
            int ci = 0;
            for (string nuc : nucs) {
                string k = nuc + "_" + suf;
                if (mh5.count(k)) {
                    f_out->cd(); mh5[k]->Write(); sh5[k]->Write();
                    TGraphErrors *m = new TGraphErrors(), *s = new TGraphErrors();
                    int c = cfg.colors_mc[ci++ % cfg.colors_mc.size()];
                    m->SetMarkerStyle(20 + (ci%4)); m->SetMarkerColor(c); m->SetLineColor(c);
                    s->SetMarkerStyle(20 + (ci%4)); s->SetMarkerColor(c); s->SetLineColor(c);
                    for (int r = 1; r <= cfg.n_rig(); ++r) {
                        if (mh5[k]->GetBinError(r, 1) > 0) {
                            double x = (cfg.rig_edges[r - 1] + cfg.rig_edges[r]) / 2.0;
                            m->SetPoint(m->GetN(), x, mh5[k]->GetBinContent(r, 1)); m->SetPointError(m->GetN() - 1, 0, mh5[k]->GetBinError(r, 1));
                            s->SetPoint(s->GetN(), x, sh5[k]->GetBinContent(r, 1)); s->SetPointError(s->GetN() - 1, 0, sh5[k]->GetBinError(r, 1));
                        }
                    }
                    gm.push_back(m); gs.push_back(s); lg.push_back(nuc);
                }
            }
            Draw(gm, lg, "Rigidity [GV]", "#mu", base + "_Mean", "_vs_Rig");
            Draw(gs, lg, "Rigidity [GV]", "#sigma", base + "_Sigma", "_vs_Rig");
        }

        for (string suf : cfg.suffixes_h4) {
            string base = GetInfo(suf).out_suffix + suf.substr(3) + "_MC_Combined";
            TGraphErrors *gm = new TGraphErrors(), *gs = new TGraphErrors();
            gm->SetMarkerStyle(20); gm->SetMarkerColor(kRed); gm->SetLineColor(kRed);
            gs->SetMarkerStyle(20); gs->SetMarkerColor(kBlue); gs->SetLineColor(kBlue);
            for (string nuc : nucs) {
                string k = nuc + "_" + suf;
                if (rh4.count(k)) {
                    double z = GetZ(nuc);
                    if (z > 0) {
                        gm->SetPoint(gm->GetN(), z, rh4[k].mean); gm->SetPointError(gm->GetN() - 1, 0, rh4[k].mean_err);
                        gs->SetPoint(gs->GetN(), z, rh4[k].sigma); gs->SetPointError(gs->GetN() - 1, 0, rh4[k].sigma_err);
                    }
                }
            }
            Draw({gm}, {"MC"}, "Charge (Z)", "#mu", base + "_Mean", "_vs_Charge");
            Draw({gs}, {"MC"}, "Charge (Z)", "#sigma", base + "_Sigma", "_vs_Charge");
        }
    }
};

void BetaStudy() {
    BetaAnalysis app;
    app.Run();
}