#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <cmath>

#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2D.h> // Added for TH2D
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TGraphErrors.h>

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooArgList.h>
#include <RooMsgService.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <RooFormulaVar.h>
#include <RooConstVar.h>

#include "../Tool.h"
using namespace AMS_Iso;
using namespace RooFit;
using namespace std;

const std::string inputFileName = "/eos/ams/group/ihep/zixuan/filter/basic_L1Q2p5to8p8.root";
const std::string outputDir = "/eos/user/z/zixuan/Isotope/PureChargeTemp/";
const std::string chainName = "UnbiasedL1Inner";
const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
const double Q_GLOBAL_MIN = 1.0;
const double Q_GLOBAL_MAX = 9.0;

const std::vector<std::string> requiredElements = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
const std::map<std::string, int> elementZ = {{"Helium", 2}, {"Lithium", 3}, {"Beryllium", 4}, {"Boron", 5}, {"Carbon", 6}, {"Nitrogen", 7}, {"Oxygen", 8}};

const std::map<std::string, std::pair<double, double>> detector_ek_ranges = {
  {"TOF", {0.25, 1.50}},
  {"NaF", {0.80, 6.10}},
  {"AGL", {2.50, 22.0}}
};

std::string getTemplateHistName(const std::string& elName, const std::string& detector, bool isL2) {
  std::string tag = isL2 ? "L2Template" : "L1Template";
  return chainName + "_BKG_H4_" + elName + "_" + tag + "_" + detector;
}

struct PurgeResult {
  int fitStatus = -1;
  double chi2ndf = 0.0;
  std::map<std::string, double> finalFractions;
  std::unique_ptr<TH1D> purifiedHist;
};

class TemplatePurger {
public:
  TemplatePurger(const string& primary, const vector<string>& contaminants, double min, double max)
    : primaryElement_(primary), contaminantElements_(contaminants), fitMin_(min), fitMax_(max) {}

  PurgeResult runFitAndPurge(
    const string& fitName,
    TH1D* h_target,
    const map<string, TH1D*>& components,
    const map<string, bool>& isFixed,
    TH1D* h_pure_temp_input,
    TCanvas* c_pdf,
    const std::string& pdfFileName,
    int fitStatus
  );

private:
  std::unique_ptr<RooPlot> setupPlot(
    const string& fitName,
    TH1D* h_target,
    RooRealVar& charge,
    RooDataHist& dataHist,
    RooAddPdf& totalPdf,
    const map<string, RooHistPdf*>& compPdfs,
    double& chi2ndf,
    TH1D* h_pure_temp_input,
    TGraphErrors*& pullGraph,
    int fitStatus
  );

  string primaryElement_;
  vector<string> contaminantElements_;
  double fitMin_, fitMax_;
};

std::unique_ptr<RooPlot> TemplatePurger::setupPlot(const string& fitName, TH1D* h_target, RooRealVar& charge, RooDataHist& dataHist, RooAddPdf& totalPdf, const map<string, RooHistPdf*>& compPdfs, double& chi2ndf, TH1D* h_pure_temp_input, TGraphErrors*& pullGraph, int fitStatus) {
  
  auto frame = std::unique_ptr<RooPlot>(charge.frame(Range(fitMin_, fitMax_), Title(Form("%s Fit", fitName.c_str()))));
  frame->GetXaxis()->SetRangeUser(fitMin_, fitMax_);
  
  double targetMax = h_target->GetMaximum();
  double targetMin = h_target->GetMinimum(0.0);
  double yMinPlot = std::max(1., targetMin * 0.1);
  
  frame->GetYaxis()->SetRangeUser(yMinPlot, 10*targetMax);
  frame->SetMinimum(yMinPlot);
  
  dataHist.plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
  
  if (fitStatus <= 1) {
    totalPdf.plotOn(frame.get(), Name("total_pdf"), LineColor(kRed), LineWidth(2));

    map<string, int> elementColors = {{"Helium", kGray}, {"Lithium", kBlue + 2}, {"Beryllium", kAzure + 7}, {"Boron", kOrange - 3}, {"Carbon", kGreen + 2}, {"Nitrogen", kMagenta - 3}, {"Oxygen", kCyan + 2}, {"Temp1", kBlue}};
    
    string primaryCompName = h_pure_temp_input ? "Temp1" : primaryElement_;

    for (const auto& el : contaminantElements_) {
      if (compPdfs.count(el)) {
        totalPdf.plotOn(frame.get(), Components(*compPdfs.at(el)), Name(Form("comp_%s", el.c_str())), LineColor(elementColors.at(el)), LineStyle(1), LineWidth(3));
      }
    }
    
    if (compPdfs.count(primaryElement_) && !h_pure_temp_input) {
      totalPdf.plotOn(frame.get(), Components(*compPdfs.at(primaryElement_)), Name(Form("comp_%s", primaryElement_.c_str())), LineColor(elementColors.at(primaryElement_)), LineWidth(2));
    } else if (h_pure_temp_input && compPdfs.count("Temp1")) {
      totalPdf.plotOn(frame.get(), Components(*compPdfs.at("Temp1")), Name(Form("comp_%s", primaryCompName.c_str())), LineColor(elementColors.at(primaryElement_)), LineWidth(2));
    }
    
    double chi2 = calculateChi2(frame.get(), "data_hist", "total_pdf", fitMin_, fitMax_);

    RooHist* h_data = frame->getHist("data_hist");
    int nBinsInRange = (h_data) ? h_data->GetN() : 0;
    int nFreeParams = totalPdf.getParameters(dataHist)->selectByAttrib("Constant", false)->getSize();
    int ndf_approx = nBinsInRange - nFreeParams;
    
    chi2ndf = (ndf_approx > 0) ? chi2 / ndf_approx : 0.0;
    
    pullGraph = new TGraphErrors();
    calculatePull(frame.get(), pullGraph, fitMin_, fitMax_);
  } else {
    chi2ndf = 9999.0;
    pullGraph = nullptr;
  }

  return frame;
}

PurgeResult TemplatePurger::runFitAndPurge(
  const string& fitName,
  TH1D* h_target,
  const map<string, TH1D*>& components,
  const map<string, bool>& isFixed,
  TH1D* h_pure_temp_input,
  TCanvas* c_pdf,
  const std::string& pdfFileName,
  int fitStatus
) {
  PurgeResult result;

  if (!c_pdf) {
    cerr << " [ERROR] c_pdf pointer is null." << endl;
    return result;
  }

  if (h_target->GetEntries() < 50) {
    cerr << " [ERROR] Target histogram too few entries: " << h_target->GetEntries() << " for " << fitName << endl;
    return result;
  }
  RooRealVar charge("charge", "Charge", Q_GLOBAL_MIN, Q_GLOBAL_MAX);
  auto dataHist = std::make_unique<RooDataHist>("data_target", "Target Data", charge, h_target);
  
  map<string, unique_ptr<RooHistPdf>> pdfs;
  map<string, unique_ptr<RooDataHist>> dHists;
  vector<unique_ptr<RooRealVar>> fracVars;
  map<string, RooAbsReal*> fracMap;
  RooArgList pdfList, fracList;

  vector<string> allComps = contaminantElements_;
  string primaryElName = primaryElement_;

  if (h_pure_temp_input) {
    primaryElName = "Temp1";
    allComps.push_back(primaryElName);
    dHists[primaryElName] = std::make_unique<RooDataHist>("dhist_temp1", "", charge, h_pure_temp_input);
    pdfs[primaryElName] = std::make_unique<RooHistPdf>("pdf_temp1", "", charge, *dHists[primaryElName]);
  } else {
    allComps.push_back(primaryElement_);
  }
  
  for (const auto& el : allComps) {
    if (el == "Temp1") continue;
    TH1D* h = components.at(el);
    dHists[el] = std::make_unique<RooDataHist>(Form("dhist_%s", el.c_str()), "", charge, h);
    pdfs[el] = std::make_unique<RooHistPdf>(Form("pdf_%s", el.c_str()), "", charge, *dHists[el]);
    pdfList.add(*pdfs.at(el));
  }
  
  if (pdfs.count(primaryElName) && !pdfList.contains(*pdfs.at(primaryElName))) {
    pdfList.add(*pdfs.at(primaryElName));
  }

  for (const auto& el : contaminantElements_) {
    double guess = 0.01;
    auto fracVar = std::make_unique<RooRealVar>(Form("frac_%s", el.c_str()), "", guess, 0.0, 1.0);
    if (isFixed.count(el) && isFixed.at(el)) fracVar->setConstant(true);
    fracMap[el] = fracVar.get();
    fracList.add(*fracMap.at(el));
    fracVars.push_back(std::move(fracVar));
  }
  
  string formula = "1.0";
  RooArgList formulaArgs;
  for (const auto& param : fracVars) {
    formula += " - @" + std::to_string(formulaArgs.getSize());
    formulaArgs.add(*param.get());
  }
  auto lastFraction = std::make_unique<RooFormulaVar>(Form("frac_%s", primaryElName.c_str()), "", formula.c_str(), formulaArgs);
  fracMap[primaryElName] = lastFraction.get();
  fracList.add(*fracMap.at(primaryElName));  

  auto totalPdf = std::make_unique<RooAddPdf>("total_pdf", "", pdfList, fracList, false);

  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  unique_ptr<RooFitResult> fitResult(
    totalPdf->fitTo(*dataHist, Save(true), PrintLevel(-1), Range(fitMin_, fitMax_), Strategy(2), Minimizer("Minuit2", "migrad"))
  );
  
  result.fitStatus = fitResult ? fitResult->status() : -1;
  
  for (const auto& el : allComps) {
    if (el == primaryElName) {
      result.finalFractions[el] = lastFraction->getVal();
    } else {
      RooAbsReal* frac_param = (RooAbsReal*)totalPdf->getParameters(*dataHist)->find(Form("frac_%s", el.c_str()));
      result.finalFractions[el] = frac_param ? frac_param->getVal() : 0.0;
    }
  }

  result.purifiedHist = unique_ptr<TH1D>((TH1D*)h_target->Clone(Form("h_pure_%s", h_target->GetName())));
  
  if (result.fitStatus <= 1) {
    double totalFitEvents = dataHist->sumEntries();

    for (const auto& el : contaminantElements_) {
      if (result.finalFractions.count(el) && components.count(el)) {
        
        double F_A = result.finalFractions.at(el);
        
        TH1D* h_template_A = components.at(el);  
        double N_template_A = h_template_A->GetSumOfWeights();  
        
        if (N_template_A > 0 && F_A > 0) {
          double N_target_A = F_A * totalFitEvents;
          double scaleFactor = N_target_A / N_template_A;

          auto h_contam_scaled = unique_ptr<TH1D>((TH1D*)h_template_A->Clone(Form("h_contam_scaled_%s_%s", el.c_str(), fitName.c_str())));
          h_contam_scaled->Scale(scaleFactor);
          
          result.purifiedHist->Add(h_contam_scaled.get(), -1.0);
          
          cout << "   [Purge] Element: " << el  
            << ", Fraction: " << F_A  
            << ", N_template: " << N_template_A
            << ", N_purged: " << N_target_A  
            << ", Scale: " << scaleFactor << endl;
        }
      }
    }
    
    for (int i = 1; i <= result.purifiedHist->GetNbinsX(); ++i) {
      double content = result.purifiedHist->GetBinContent(i);
      if (content < 0) result.purifiedHist->SetBinContent(i, 0.0);
    }

  } else {
    result.purifiedHist = unique_ptr<TH1D>((TH1D*)h_target->Clone(Form("h_failed_pure_%s", h_target->GetName())));
  }

  map<string, RooHistPdf*> compPdfs_raw_ptrs;
  for (const auto& el : allComps) {
    if (pdfs.count(el)) {
      compPdfs_raw_ptrs[el] = pdfs.at(el).get();
    }
  }

  double chi2ndf_val = 0.0;
  TGraphErrors* pullGraphPtr = nullptr;
  
  auto frame = setupPlot(fitName, h_target, charge, *dataHist, *totalPdf,
               compPdfs_raw_ptrs,
               chi2ndf_val, h_pure_temp_input,
               pullGraphPtr, result.fitStatus);
  
  unique_ptr<TGraphErrors> pullGraph(pullGraphPtr);
  result.chi2ndf = chi2ndf_val;

  c_pdf->Clear();
  c_pdf->Divide(1, 2);

  c_pdf->cd(1);
  TPad* pad1 = (TPad*)gPad;
  pad1->SetPad(0, 0.3, 1, 1);
  h_target->SetMinimum(1);
  h_target->GetXaxis()->SetRangeUser(fitMin_, fitMax_);
  h_target->Draw("PZ");
  frame->SetMinimum(1);
  
  pad1->SetLogy();
  pad1->SetBottomMargin(0.02);

  pad1->cd();

  frame->GetYaxis()->SetTitle("Events");
  frame->GetXaxis()->SetLabelSize(0);
  frame->Draw();
  
  auto legend = std::make_unique<TLegend>(0.7, 0.55, 0.88, 0.88);
  legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.03);
  legend->AddEntry("data_hist", Form("%s L1 Raw", primaryElement_.c_str()), "pe");
  if (result.fitStatus <= 1 && frame->findObject("total_pdf")) legend->AddEntry("total_pdf", "Total Fit", "l");
  
  string primaryLegendName = primaryElement_;
  if (h_pure_temp_input) primaryLegendName += " Pure Temp1"; else primaryLegendName += " L2";

  string primaryCompDrawName = h_pure_temp_input ? "comp_Temp1" : Form("comp_%s", primaryElement_.c_str());
  if(result.fitStatus <= 1 && frame->findObject(primaryCompDrawName.c_str())) legend->AddEntry(primaryCompDrawName.c_str(), primaryLegendName.c_str(), "l");

  for (const auto& el : contaminantElements_) {
    if (result.fitStatus <= 1 && frame->findObject(Form("comp_%s", el.c_str()))) {
      double frac = result.finalFractions.count(el) ? result.finalFractions.at(el) : 0.0;
      legend->AddEntry(Form("comp_%s", el.c_str()), Form("%s Contam. (F=%.3f)", el.c_str(), frac), "l");
    }
  }
  legend->Draw("same");

  auto info = std::make_unique<TPaveText>(0.15, 0.65, 0.65, 0.88, "NDC");
  info->SetFillStyle(0); info->SetBorderSize(0); info->SetTextAlign(12);
  info->SetTextSize(0.03);
  info->AddText(fitName.c_str());
  info->AddText(Form("#chi^{2}/NDF = %.2f (Status: %d)", result.chi2ndf, result.fitStatus));
  if (result.finalFractions.count(primaryElName)) {
    string fracs_line = Form("F_{%s}=%.3f", primaryElement_.c_str(), result.finalFractions.at(primaryElName));
    info->AddText(fracs_line.c_str());
  }
  info->Draw("same");
  
  c_pdf->cd(2);
  TPad* pad2 = (TPad*)gPad;
  pad2->SetPad(0, 0, 1, 0.3); pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.3); pad2->SetGridy();
  
  if (pullGraph) {
    pullGraph->SetTitle("");
    pullGraph->SetMarkerStyle(20);
    pullGraph->SetMarkerSize(0.8);
    pullGraph->Draw("AP");
    
    pullGraph->GetYaxis()->SetTitle("Pull");
    pullGraph->GetYaxis()->SetTitleSize(0.1);
    pullGraph->GetYaxis()->SetTitleOffset(0.4);
    pullGraph->GetYaxis()->SetLabelSize(0.1);
    pullGraph->GetYaxis()->SetRangeUser(-6, 6);
    pullGraph->GetYaxis()->SetNdivisions(505);

    pullGraph->GetXaxis()->SetTitle("Charge");
    pullGraph->GetXaxis()->SetTitleSize(0.12);
    pullGraph->GetXaxis()->SetTitleOffset(1.0);
    pullGraph->GetXaxis()->SetLabelSize(0.1);
    
    TLine line;
    line.SetLineColor(kRed);
    line.SetLineWidth(1);
    line.DrawLine(fitMin_, 0, fitMax_, 0);

  }

  c_pdf->Update();
  c_pdf->Print(pdfFileName.c_str());

  return result;
}

void CalPureQTemp() {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  TH1::AddDirectory(kFALSE);
  
  std::string pdfFileName = outputDir + "CalPureQTemp_" + chainName + ".pdf";
  std::string rootFileName = outputDir + "PureChargeTemplates_" + chainName + ".root";
  system(Form("mkdir -p %s", outputDir.c_str()));

  auto inputFile = std::unique_ptr<TFile>(TFile::Open(inputFileName.c_str()));
  if (!inputFile || inputFile->IsZombie()) {
    cerr << "CRITICAL: Could not open input file: " << inputFileName << endl;
    return;
  }
  auto outputFile = std::make_unique<TFile>(rootFileName.c_str(), "RECREATE");
  
  TCanvas* c_pdf = new TCanvas("c_pdf", "PDF Canvas", 800, 600);
  c_pdf->Print((pdfFileName + "[").c_str());

  const std::vector<std::string> orderedElements = {"Oxygen", "Nitrogen", "Carbon", "Boron", "Beryllium", "Lithium", "Helium"};

  const std::map<std::string, std::vector<std::string>> contaminationMap = {
    {"Oxygen", {}},
    {"Nitrogen", {"Oxygen"}},
    {"Carbon", {"Nitrogen", "Oxygen"}},
    {"Boron", {"Carbon", "Nitrogen", "Oxygen"}},
    {"Beryllium", {"Boron", "Carbon", "Nitrogen", "Oxygen"}},
    {"Lithium", {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"}},
    {"Helium", {}}
  };
  
  for (const auto& detector : detectors) {
    cout << "\n--- Processing Detector: " << detector << " ---" << endl;

    if (detector_ek_ranges.find(detector) == detector_ek_ranges.end()) {
      cerr << "WARNING: No Ek range defined for detector " << detector << ". Skipping." << endl;
      continue;
    }
    double ek_min = detector_ek_ranges.at(detector).first;
    double ek_max = detector_ek_ranges.at(detector).second;

    map<string, map<int, unique_ptr<TH1D>>> finalPureTemplates;
    
    map<string, unique_ptr<TH2F>> h2d_raw_L1, h2d_raw_L2;
    for (const auto& el : requiredElements) {
      string nameL1 = getTemplateHistName(el, detector, false);
      string nameL2 = getTemplateHistName(el, detector, true);

      TH2F* h2dL1_raw = (TH2F*)inputFile->Get(nameL1.c_str());
      TH2F* h2dL2_raw = (TH2F*)inputFile->Get(nameL2.c_str());

      if (!h2dL1_raw) { cerr << "L1 template not found for " << el << " in " << detector << endl; continue; }
      
      h2d_raw_L1[el] = unique_ptr<TH2F>((TH2F*)h2dL1_raw->Clone(Form("L1_RAW_%s_%s", el.c_str(), detector.c_str())));
      h2d_raw_L1[el]->RebinX(2);
      
      if (h2dL2_raw) {
        h2d_raw_L2[el] = unique_ptr<TH2F>((TH2F*)h2dL2_raw->Clone(Form("L2_RAW_%s_%s", el.c_str(), detector.c_str())));
        h2d_raw_L2[el]->RebinX(2);
      } else {
        cerr << "L2 template not found for " << el << " in " << detector << " (Using L1 raw as L2 fallback)" << endl;
        h2d_raw_L2[el] = unique_ptr<TH2F>((TH2F*)h2dL1_raw->Clone(Form("L2_FALLBACK_%s_%s", el.c_str(), detector.c_str())));
        h2d_raw_L2[el]->RebinX(2);
      }
    }
    
    if (h2d_raw_L1.empty()) continue;

    const TAxis* y_axis = h2d_raw_L1.begin()->second->GetYaxis();
    int n_bins_y = y_axis->GetNbins();

    for (int y_bin = 1; y_bin <= n_bins_y; ++y_bin) {
      double ek_center = y_axis->GetBinCenter(y_bin);
      
      if (ek_center < ek_min || ek_center > ek_max) {
        continue;
      }
      
      cout << " [Bin " << y_bin << "] Ek=" << ek_center << " GeV/n. Analyzing..." << endl;

      map<string, unique_ptr<TH1D>> currentBinPureTemplates;

      for (const auto& elPrimary : orderedElements) {
        if (!h2d_raw_L1.count(elPrimary) || !h2d_raw_L2.count(elPrimary)) {
          cerr << "WARNING: Missing raw templates for " << elPrimary << ". Skipping purge for this element." << endl;
          continue;
        }
        
        const auto& contaminants = contaminationMap.at(elPrimary);
        map<string, TH1D*> comps1;

        if (elPrimary == "Oxygen" || elPrimary == "Helium") {
          auto h_raw = unique_ptr<TH1D>(h2d_raw_L1.at(elPrimary)->ProjectionX(Form("h_slice_L1_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
          h_raw->Smooth(1,"G");
          currentBinPureTemplates[elPrimary] = unique_ptr<TH1D>((TH1D*)h_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin)));
          finalPureTemplates[elPrimary][y_bin] = unique_ptr<TH1D>((TH1D*)currentBinPureTemplates.at(elPrimary)->Clone());
          continue;
        }
        
        bool allContaminantsReady = true;
        for (const auto& contam : contaminants) {
          if (!currentBinPureTemplates.count(contam)) {
            allContaminantsReady = false;
            cerr << "ERROR: Pure template for contaminant " << contam << " (needed for " << elPrimary << ") not found. Skipping." << endl;
            break;
          }
          comps1[contam] = currentBinPureTemplates.at(contam).get();
        }
        if (!allContaminantsReady) continue;

        auto h_Primary_raw = unique_ptr<TH1D>(h2d_raw_L1.at(elPrimary)->ProjectionX(Form("h_slice_L1_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
        h_Primary_raw->Smooth(1,"G");
        auto h_Primary_L2 = unique_ptr<TH1D>(h2d_raw_L2.at(elPrimary)->ProjectionX(Form("h_slice_L2_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
        h_Primary_L2->Smooth(1,"G");

        comps1[elPrimary] = h_Primary_L2.get();
        
        double Q_min = elementZ.at(elPrimary) - 0.4;
        double Q_max = Q_GLOBAL_MAX-0.6;
        
        TemplatePurger purger1(elPrimary, contaminants, Q_min, Q_max);
        PurgeResult res1 = purger1.runFitAndPurge(
          Form("%s Purge 1 (L2+PureContams) E%d %s", elPrimary.c_str(), y_bin, detector.c_str()), 
          h_Primary_raw.get(), comps1, {}, nullptr, c_pdf, pdfFileName, -1
        );
        
        if (res1.fitStatus <= 1) {
          TemplatePurger purger2(elPrimary, contaminants, Q_min, Q_max);
          map<string, TH1D*> comps2 = comps1;
          comps2[elPrimary] = res1.purifiedHist.get();

          PurgeResult res2 = purger2.runFitAndPurge(
            Form("%s Purge 2 (Temp1+PureContams) E%d %s", elPrimary.c_str(), y_bin, detector.c_str()), 
            h_Primary_raw.get(), comps2, {}, res1.purifiedHist.get(), c_pdf, pdfFileName, res1.fitStatus
          );
          
          if (res2.fitStatus <= 1) {
            currentBinPureTemplates[elPrimary] = unique_ptr<TH1D>((TH1D*)res2.purifiedHist->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin)));
            finalPureTemplates[elPrimary][y_bin] = unique_ptr<TH1D>((TH1D*)currentBinPureTemplates.at(elPrimary)->Clone());
          } else {
            currentBinPureTemplates[elPrimary] = unique_ptr<TH1D>((TH1D*)h_Primary_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin)));
          }
        } else {
          currentBinPureTemplates[elPrimary] = unique_ptr<TH1D>((TH1D*)h_Primary_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin)));
        }
      }
    }
    
    // Save Pure L1 Templates
    for (const auto& el : requiredElements) {
      outputFile->cd();
      if (!finalPureTemplates.count(el) || !h2d_raw_L1.count(el)) continue;

      const auto& h_base_2d = h2d_raw_L1.at(el);
      std::string pureHistName = Form("h2d_PureQTemp_%s_%s", el.c_str(), detector.c_str());
      auto h_pure_2d = unique_ptr<TH2F>((TH2F*)h_base_2d->Clone(pureHistName.c_str()));
            
      h_pure_2d->Reset();

      for (int y_bin = 1; y_bin <= n_bins_y; ++y_bin) {
        if (finalPureTemplates.at(el).count(y_bin)) {
          TH1D* h_slice = finalPureTemplates.at(el).at(y_bin).get();
          for (int x_bin = 1; x_bin <= h_slice->GetNbinsX(); ++x_bin) {
            h_pure_2d->SetBinContent(x_bin, y_bin, h_slice->GetBinContent(x_bin));
          }
        }
      }
      h_pure_2d->Write();
      cout << "Saved " << pureHistName << endl;
    }

    // Save Tune L1 Templates (He & Oxy)
    for (const auto& el : {"Helium", "Oxygen"}) {
      if (h2d_raw_L1.count(el)) {
        outputFile->cd();
        // Retrieve as TH2F
        const auto& h_base_l1 = h2d_raw_L1.at(el);
        std::string tuneHistNameL1 = Form("h2d_TuneQTemp_%s_%s", el, detector.c_str());
        // Convert to TH2D for saving
        auto h_tune_2d_l1 = std::make_unique<TH2D>(tuneHistNameL1.c_str(), h_base_l1->GetTitle(), 
                                                   h_base_l1->GetNbinsX(), h_base_l1->GetXaxis()->GetXmin(), h_base_l1->GetXaxis()->GetXmax(),
                                                   h_base_l1->GetNbinsY(), h_base_l1->GetYaxis()->GetXmin(), h_base_l1->GetYaxis()->GetXmax());
        h_tune_2d_l1->GetYaxis()->Set(h_base_l1->GetNbinsY(), h_base_l1->GetYaxis()->GetXbins()->GetArray());
        
        for(int x=1; x<=h_base_l1->GetNbinsX(); ++x) {
            for(int y=1; y<=h_base_l1->GetNbinsY(); ++y) {
                h_tune_2d_l1->SetBinContent(x, y, h_base_l1->GetBinContent(x, y));
                h_tune_2d_l1->SetBinError(x, y, h_base_l1->GetBinError(x, y));
            }
        }
        h_tune_2d_l1->Write();
        cout << "Saved " << tuneHistNameL1 << endl;
      }
    }

    // Save Tune L2 Templates (Li to Nit)
    for (const auto& el : {"Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen"}) {
      if (h2d_raw_L2.count(el)) {
        outputFile->cd();
        // Retrieve as TH2F
        const auto& h_base_l2 = h2d_raw_L2.at(el);
        std::string tuneHistNameL2 = Form("h2d_TuneQTemp_%s_%s", el, detector.c_str());
        // Convert to TH2D for saving
        auto h_tune_2d_l2 = std::make_unique<TH2D>(tuneHistNameL2.c_str(), h_base_l2->GetTitle(), 
                                                   h_base_l2->GetNbinsX(), h_base_l2->GetXaxis()->GetXmin(), h_base_l2->GetXaxis()->GetXmax(),
                                                   h_base_l2->GetNbinsY(), h_base_l2->GetYaxis()->GetXmin(), h_base_l2->GetYaxis()->GetXmax());
        h_tune_2d_l2->GetYaxis()->Set(h_base_l2->GetNbinsY(), h_base_l2->GetYaxis()->GetXbins()->GetArray());

        for(int x=1; x<=h_base_l2->GetNbinsX(); ++x) {
            for(int y=1; y<=h_base_l2->GetNbinsY(); ++y) {
                h_tune_2d_l2->SetBinContent(x, y, h_base_l2->GetBinContent(x, y));
                h_tune_2d_l2->SetBinError(x, y, h_base_l2->GetBinError(x, y));
            }
        }
        h_tune_2d_l2->Write();
        cout << "Saved " << tuneHistNameL2 << endl;
      }
    }
  }

  c_pdf->Print((pdfFileName + "]").c_str());
}