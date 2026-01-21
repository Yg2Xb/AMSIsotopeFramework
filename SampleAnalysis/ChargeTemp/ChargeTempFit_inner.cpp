#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <map>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooArgList.h>
#include <RooMsgService.h>
#include <RooFormulaVar.h>
#include <RooCmdArg.h>
#include "../Tool.h" 

using namespace AMS_Iso; 
using namespace RooFit;
using namespace std;

const std::string outputDir = "/eos/user/z/zixuan/Isotope/ChargeTemp/";
const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
const std::map<std::string, std::pair<double, double>> detector_ek_ranges = {
  {"TOF", {0.25, 1.50}}, {"NaF", {0.61, 6.10}}, {"AGL", {2.50, 23.0}}
};

struct UnifiedFitResult {
    double ekpernuc_low = 0.0, ekpernuc_up = 0.0;
    double chi2ndf = 0.0;
    int ndf = 0;
    std::map<std::string, double> yield_val;
    std::map<std::string, double> yield_err;
};

// Fitter类保持不变
class UnifiedChargeFitter {
public:
    UnifiedChargeFitter(const std::string& det, int bin, TH1D* h_sig, const std::map<std::string, TH1D*>& h_temps,
                        double fMin, double fMax, const std::vector<std::string>& elms)
        : detName_(det), energyBin_(bin), fitMin_(fMin), fitMax_(fMax), templateElements_(elms) {
        h_signal_raw_ = std::unique_ptr<TH1D>((TH1D*)h_sig->Clone());
        for(auto const& [el, h] : h_temps) h_templates_raw_[el] = std::unique_ptr<TH1D>((TH1D*)h->Clone());
    }

    bool initialize() {
        if (!h_signal_raw_) return false;
        double sigEntries = h_signal_raw_->GetEntries();
        if (sigEntries < 5) return false;

        charge_ = std::make_unique<RooRealVar>("charge", "InnerQ", fitMin_, fitMax_);
        h_signal_ext_ = extendHistogram(h_signal_raw_.get(), fitMin_, fitMax_);
        data_hist_ = std::make_unique<RooDataHist>("data_hist", "Data", *charge_, h_signal_ext_.get());
        
        double totalEntries = h_signal_ext_->Integral();
        RooArgList pdfList;
        RooArgList yieldList;

        for (const auto& el : templateElements_) {
            if (h_templates_raw_.find(el) == h_templates_raw_.end()) continue;
            h_templates_ext_[el] = extendHistogram(h_templates_raw_[el].get(), fitMin_, fitMax_);
            for(int i=1; i<=h_templates_ext_[el]->GetNbinsX(); ++i) {
                if(h_templates_ext_[el]->GetBinContent(i) <= 0) h_templates_ext_[el]->SetBinContent(i, 1e-9);
            }
            tmp_dists_[el] = std::make_unique<RooDataHist>(Form("dh_%s", el.c_str()), "", *charge_, h_templates_ext_[el].get());
            tmp_pdfs_[el] = std::make_unique<RooHistPdf>(Form("pdf_%s", el.c_str()), "", *charge_, *tmp_dists_[el]);
            auto y = std::make_unique<RooRealVar>(Form("yield_%s", el.c_str()), el.c_str(), totalEntries/templateElements_.size(), 0, totalEntries * 1.5);
            yield_ptr_map_[el] = y.get();
            yieldParams_.push_back(std::move(y));
            pdfList.add(*tmp_pdfs_[el]);
            yieldList.add(*yield_ptr_map_[el]);
        }
        total_pdf_ = std::make_unique<RooAddPdf>("total_pdf", "total_pdf", pdfList, yieldList);
        return true;
    }

    bool runFit() {
        if (!total_pdf_) return false;
        fitResult_ = std::unique_ptr<RooFitResult>(total_pdf_->fitTo(*data_hist_, Extended(kTRUE), Save(), PrintLevel(-1), Strategy(1)));
        return (fitResult_ && (fitResult_->status() == 0 || fitResult_->status() == 1));
    }

    UnifiedFitResult getResult() {
        UnifiedFitResult res;
        for (const auto& el : templateElements_) {
            if (yield_ptr_map_.find(el) != yield_ptr_map_.end()) {
                res.yield_val[el] = yield_ptr_map_[el]->getVal();
                res.yield_err[el] = yield_ptr_map_[el]->getError();
            }
        }
        return res;
    }

    void draw(TCanvas* c, const std::string& pdfName, double ek_l, double ek_u, UnifiedFitResult& res) {
        c->Clear();
        c->Divide(1, 2);
        TPad* p1 = (TPad*)c->cd(1);
        p1->SetPad(0, 0.3, 1, 1);
        p1->SetLogy(1);
        p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.12); p1->SetRightMargin(0.05);

        auto frame = std::unique_ptr<RooPlot>(charge_->frame());
        data_hist_->plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
        total_pdf_->plotOn(frame.get(), Name("total_pdf"), LineColor(kRed), LineWidth(2));

        std::vector<int> colors = {kCyan-3, kBlue+1, kGreen+2, kOrange+7, kMagenta+1, kAzure-2};
        for(size_t i=0; i<templateElements_.size(); ++i) {
            std::string compName = Form("pdf_%s", templateElements_[i].c_str());
            total_pdf_->plotOn(frame.get(), Components(compName.c_str()), Name(Form("plot_comp_%s", templateElements_[i].c_str())),
                               LineColor(colors[i % colors.size()]), LineStyle(1), LineWidth(2));
        }

        frame->SetTitle(""); frame->GetYaxis()->SetTitle("Events");
        frame->GetYaxis()->SetTitleSize(0.055); frame->GetYaxis()->SetTitleOffset(0.85);
        frame->GetXaxis()->SetLabelSize(0); frame->SetMinimum(0.5);
        frame->Draw();

        auto legend = new TLegend(0.65, 0.65, 0.88, 0.88);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.035);
        legend->AddEntry("data_hist", "Data", "pe");
        legend->AddEntry("total_pdf", "Total Fit", "l");
        for(auto& e : templateElements_) legend->AddEntry(Form("plot_comp_%s", e.c_str()), e.substr(0, 3).c_str(), "l");
        legend->Draw();

        auto info = new TPaveText(0.14, 0.45, 0.44, 0.88, "NDC");
        info->SetFillStyle(0); info->SetBorderSize(0); info->SetTextAlign(12); info->SetTextSize(0.035);
        info->AddText(Form("%s: %.2f - %.2f GeV/n", detName_.c_str(), ek_l, ek_u));
        
        if (fitResult_) {
            int nFreeParams = fitResult_->floatParsFinal().getSize();
            int nBinsInRange = 0;
            for(int i = 1; i <= h_signal_ext_->GetNbinsX(); ++i) {
                if (h_signal_ext_->GetBinCenter(i) >= fitMin_ && h_signal_ext_->GetBinCenter(i) <= fitMax_) nBinsInRange++;
            }
            res.ndf = nBinsInRange - nFreeParams;
            double chi2 = calculateChi2(frame.get(), "data_hist", "total_pdf", fitMin_, fitMax_);
            if (res.ndf > 0) {
                res.chi2ndf = chi2 / res.ndf;
                info->AddText(Form("#chi^{2}/ndf = %.1f / %d = %.2f", chi2, res.ndf, res.chi2ndf));
            }
            info->AddLine(0, 0, 1, 0);
            for (const auto& el : templateElements_) info->AddText(Form("%s: %.1f #pm %.1f", el.substr(0, 3).c_str(), res.yield_val[el], res.yield_err[el]));
        }
        info->Draw();

        TPad* p2 = (TPad*)c->cd(2);
        p2->SetPad(0, 0, 1, 0.3);
        p2->SetTopMargin(0.02); p2->SetBottomMargin(0.4); p2->SetLeftMargin(0.12); p2->SetRightMargin(0.05); p2->SetGridy();

        TGraphErrors* pullG = new TGraphErrors();
        calculatePull(frame.get(), pullG, fitMin_, fitMax_);
        pullG->SetTitle(""); pullG->GetXaxis()->SetTitle("InnerQ"); pullG->GetXaxis()->SetTitleSize(0.12);
        pullG->GetXaxis()->SetLabelSize(0.12); pullG->GetYaxis()->SetTitle("Pull"); pullG->GetYaxis()->SetTitleSize(0.12);
        pullG->GetYaxis()->SetLabelSize(0.12); pullG->GetXaxis()->SetRangeUser(fitMin_, fitMax_);
        pullG->GetYaxis()->SetRangeUser(-5.5, 5.5); pullG->SetMarkerStyle(20); pullG->SetMarkerSize(0.7);
        pullG->Draw("AP");
        c->Print(pdfName.c_str());
    }

private:
    std::string detName_; int energyBin_; double fitMin_, fitMax_;
    std::vector<std::string> templateElements_;
    std::unique_ptr<TH1D> h_signal_raw_, h_signal_ext_;
    std::map<std::string, std::unique_ptr<TH1D>> h_templates_raw_, h_templates_ext_;
    std::unique_ptr<RooRealVar> charge_;
    std::unique_ptr<RooAddPdf> total_pdf_;
    std::unique_ptr<RooFitResult> fitResult_;
    std::unique_ptr<RooDataHist> data_hist_;
    std::map<std::string, std::unique_ptr<RooDataHist>> tmp_dists_;
    std::map<std::string, std::unique_ptr<RooHistPdf>> tmp_pdfs_;
    std::vector<std::unique_ptr<RooRealVar>> yieldParams_;
    std::map<std::string, RooRealVar*> yield_ptr_map_;
};

void runUnifiedChargeAnalysis(const std::string& chain, const std::string& mode) {
    auto inFile = std::unique_ptr<TFile>(TFile::Open("/eos/user/z/zixuan/Isotope/Add/Be_frag4_withBkg_NoTune_full.root"));
    if (!inFile || inFile->IsZombie()) return;

    // 创建输出文件保存直方图
    TFile* fOut = TFile::Open(Form("%sFitResults_%s.root", outputDir.c_str(), chain.c_str()), "RECREATE");

    std::string pdf_fn = outputDir + mode + "QFit_YieldBased_" + chain + ".pdf";
    auto c = std::make_unique<TCanvas>("c", "", 800, 1000);
    c->Print((pdf_fn + "[").c_str());

    for (const auto& det : detectors) {
        TH3F* h3 = (TH3F*)inFile->Get(Form("%s_BKG_H5_%s_type0", chain.c_str(), det.c_str()));
        if (!h3) continue;

        int nBinsEk = h3->GetZaxis()->GetNbins();

        // 核心修改：增加 Z 循环 (5到8)
        for (int Z = 5; Z <= 8; ++Z) {
            // 初始化保存 Be 拟合数和 Cut 数的直方图，名字包含 Z
            TH1D* h_fit_Be_Z = (h3->GetZaxis()->GetXbins()->GetSize() > 0) ? 
                new TH1D(Form("h_fit_Be_Z%d_%s", Z, det.c_str()), "", nBinsEk, h3->GetZaxis()->GetXbins()->GetArray()) :
                new TH1D(Form("h_fit_Be_Z%d_%s", Z, det.c_str()), "", nBinsEk, h3->GetZaxis()->GetXmin(), h3->GetZaxis()->GetXmax());
            TH1D* h_cut_Be_Z = (TH1D*)h_fit_Be_Z->Clone(Form("h_cut_Be_Z%d_%s", Z, det.c_str()));

            std::vector<std::string> elms = {"Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
            std::map<std::string, TH2F*> h2_tmps;
            for(auto& e : elms) h2_tmps[e] = (TH2F*)inFile->Get(Form("%s_BKG_H4_%s_InnerQTemplate_%s", chain.c_str(), e.c_str(), det.c_str()));

            for (int i=1; i<=nBinsEk; ++i) {
                double ek = h3->GetZaxis()->GetBinCenter(i);
                if (ek < detector_ek_ranges.at(det).first || ek > detector_ek_ranges.at(det).second) continue;

                h3->GetZaxis()->SetRange(i, i);
                TH2D* h2_yx = (TH2D*)h3->Project3D("yx");
                
                // 根据 Z 选择 y 轴范围 Z-0.25 到 Z+0.4
                int binY_low = h2_yx->GetYaxis()->FindBin(Z - 0.2);
                int binY_up  = h2_yx->GetYaxis()->FindBin(Z + 0.4);
                TH1D* h_sig = h2_yx->ProjectionX(Form("hsig_%s_%d_Z%d", det.c_str(), i, Z), binY_low, binY_up);
                h_sig->Rebin(2);

                // 保存 X 轴 3.5-4.5 内的事例数
                double n_cut = h_sig->Integral(h_sig->FindBin(3.5), h_sig->FindBin(4.5));
                h_cut_Be_Z->SetBinContent(i, n_cut); h_cut_Be_Z->SetBinError(i, sqrt(n_cut));

                std::map<std::string, TH1D*> h_tmps_1d;
                for(auto& e : elms) {
                    if(h2_tmps[e]) { h_tmps_1d[e] = h2_tmps[e]->ProjectionX(Form("ht_%s_%d_Z%d", e.c_str(), i, Z), i, i); h_tmps_1d[e]->Rebin(2); }
                }

                UnifiedChargeFitter fitter(det, i, h_sig, h_tmps_1d, 2.5, 8.5, elms);
                if (fitter.initialize() && fitter.runFit()) {
                    UnifiedFitResult res = fitter.getResult();
                    // 只有在 Z 循环内部进行 draw 才能看到不同 Z 范围的拟合 PDF
                    fitter.draw(c.get(), pdf_fn, h3->GetZaxis()->GetBinLowEdge(i), h3->GetZaxis()->GetBinUpEdge(i), res);
                    
                    // 保存拟合得到的 Be 数量
                    h_fit_Be_Z->SetBinContent(i, res.yield_val["Beryllium"]);
                    h_fit_Be_Z->SetBinError(i, res.yield_err["Beryllium"]);
                }
                delete h2_yx; delete h_sig;
                for(auto& [el, h] : h_tmps_1d) if(h) delete h;
            }
            // 写入文件
            fOut->cd();
            h_fit_Be_Z->Write();
            h_cut_Be_Z->Write();
        }
    }
    c->Print((pdf_fn + "]").c_str());
    fOut->Close();
}

void ChargeTempFit_inner() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    runUnifiedChargeAnalysis("L1Inner", "Pure");
}