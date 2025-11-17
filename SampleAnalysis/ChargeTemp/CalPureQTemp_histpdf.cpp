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
#include <RooConstVar.h> // 引入 RooConstVar

// ** 引入用户自定义的工具函数 **
#include "../Tool.h"
using namespace AMS_Iso;
using namespace RooFit;
using namespace std;

// --- 全局配置 ---
// 请确保这些路径和名称正确
const std::string inputFileName = "/eos/user/z/zixuan/Isotope/Add/Be_frag4.root";
const std::string outputDir = "/eos/user/z/zixuan/Isotope/PureChargeTemp/";
const std::string chainName = "UnbiasedL1Inner";
const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
const double Q_GLOBAL_MIN = 3.0; // 统一的Q范围
const double Q_GLOBAL_MAX = 9.0;

const std::vector<std::string> requiredElements = {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
const std::map<std::string, int> elementZ = {{"Beryllium", 4}, {"Boron", 5}, {"Carbon", 6}, {"Nitrogen", 7}, {"Oxygen", 8}};

// ** 新增：探测器能量范围配置 **
const std::map<std::string, std::pair<double, double>> detector_ek_ranges = {
  {"TOF", {0.33, 1.29}},
  {"NaF", {0.90, 5.10}},
  {"AGL", {2.90, 21.0}}
};

// 模板名称辅助函数
std::string getTemplateHistName(const std::string& elName, const std::string& detector, bool isL2) {
  std::string tag = isL2 ? "L2QTemplate" : "L1QTemplate";
  return chainName + "_ISS_BKG_H2_" + elName + "_" + tag + "_" + detector;
}

// 统一的拟合-扣除辅助类
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
    TH1D* h_pure_temp_input, // 用于Fit 2的临时纯净模板 (Temp1)
    TCanvas* c_pdf, // 传入用于绘图的 Canvas
    const std::string& pdfFileName, // 传入 PDF 文件名
    int fitStatus // 新增：传入拟合状态
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
    int fitStatus // 新增：传入拟合状态
  );

  string primaryElement_;
  vector<string> contaminantElements_;
  double fitMin_, fitMax_;
};

// --- PurgeTool 实现 ---

std::unique_ptr<RooPlot> TemplatePurger::setupPlot(const string& fitName, TH1D* h_target, RooRealVar& charge, RooDataHist& dataHist, RooAddPdf& totalPdf, const map<string, RooHistPdf*>& compPdfs, double& chi2ndf, TH1D* h_pure_temp_input, TGraphErrors*& pullGraph, int fitStatus) {
 
  auto frame = std::unique_ptr<RooPlot>(charge.frame(Range(fitMin_, fitMax_), Title(Form("%s Fit", fitName.c_str()))));
  frame->GetXaxis()->SetRangeUser(fitMin_, fitMax_);
 
  // 调整 Y 轴下限以适应 Log 轴
  double targetMax = h_target->GetMaximum();
  double targetMin = h_target->GetMinimum(0.0);
  double yMinPlot = std::max(1., targetMin * 0.1);
 
  frame->GetYaxis()->SetRangeUser(yMinPlot, 10*targetMax);
  frame->SetMinimum(yMinPlot);
 
  // Plot Data (Raw L1)
  dataHist.plotOn(frame.get(), Name("data_hist"), MarkerStyle(20), MarkerSize(0.8));
 
  // Plot Total PDF (仅在拟合成功时才绘制Total PDF和分量)
  if (fitStatus <= 1) {
    // ** 修正点 1：plotOn 返回 RooPlot*，我们不需要捕获它 **
    totalPdf.plotOn(frame.get(), Name("total_pdf"), LineColor(kRed), LineWidth(2));

    map<string, int> elementColors = {{"Beryllium", kAzure + 7}, {"Boron", kOrange - 3}, {"Carbon", kGreen + 2}, {"Nitrogen", kMagenta - 3}, {"Oxygen", kCyan + 2}, {"Temp1", kBlue}};
   
    string primaryCompName = h_pure_temp_input ? "Temp1" : primaryElement_;

    // 绘制污染成分
    for (const auto& el : contaminantElements_) {
      if (compPdfs.count(el)) {
        totalPdf.plotOn(frame.get(), Components(*compPdfs.at(el)), Name(Form("comp_%s", el.c_str())), LineColor(elementColors.at(el)), LineStyle(2), LineWidth(2));
      }
    }
   
    // 绘制主成分 (L2 或 Temp1)
    if (compPdfs.count(primaryElement_) && !h_pure_temp_input) {
      totalPdf.plotOn(frame.get(), Components(*compPdfs.at(primaryElement_)), Name(Form("comp_%s", primaryElement_.c_str())), LineColor(elementColors.at(primaryElement_)), LineWidth(2));
    } else if (h_pure_temp_input && compPdfs.count("Temp1")) {
      totalPdf.plotOn(frame.get(), Components(*compPdfs.at("Temp1")), Name(Form("comp_%s", primaryCompName.c_str())), LineColor(elementColors.at(primaryElement_)), LineWidth(2));
    }
   
    // ** 只有在成功拟合的情况下才计算 Chi2 和 Pull **
    double chi2 = calculateChi2(frame.get(), "data_hist", "total_pdf", fitMin_, fitMax_);

    RooHist* h_data = frame->getHist("data_hist");
    int nBinsInRange = (h_data) ? h_data->GetN() : 0;
    int nFreeParams = totalPdf.getParameters(dataHist)->selectByAttrib("Constant", false)->getSize();
    int ndf_approx = nBinsInRange - nFreeParams;
   
    chi2ndf = (ndf_approx > 0) ? chi2 / ndf_approx : 0.0;
   
    pullGraph = new TGraphErrors();
    calculatePull(frame.get(), pullGraph, fitMin_, fitMax_);
  } else {
    // 拟合失败，只绘制数据点，不绘制曲线
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
  int fitStatus // 传入拟合状态，用于绘图判断
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
    pdfList.add(*pdfs.at(el)); // ** 将所有 PDF 加入列表 **
  }
 
  // ** 确保 primary 模板的 PDF 已经被加入列表 **
  if (pdfs.count(primaryElName) && !pdfList.contains(*pdfs.at(primaryElName))) {
    pdfList.add(*pdfs.at(primaryElName));
  }


  // 创建 Fraction Variables (Contaminants)
  for (const auto& el : contaminantElements_) {
    double guess = 0.01;
    auto fracVar = std::make_unique<RooRealVar>(Form("frac_%s", el.c_str()), "", guess, 0.0, 1.0);
    if (isFixed.count(el) && isFixed.at(el)) fracVar->setConstant(true);
    fracMap[el] = fracVar.get();
    fracList.add(*fracMap.at(el));
    fracVars.push_back(std::move(fracVar));
  }
 
  // Primary Fraction (1 - sum(contaminant fractions))
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

  // Run fit
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  unique_ptr<RooFitResult> fitResult(
    totalPdf->fitTo(*dataHist, Save(true), PrintLevel(-1), Range(fitMin_, fitMax_), Strategy(2), Minimizer("Minuit2", "migrad"))
  );
 
  result.fitStatus = fitResult ? fitResult->status() : -1;
 
  // 拟合后更新 Fraction 映射
  for (const auto& el : allComps) {
    // 尝试从拟合结果中获取变量，如果是最后一个 FormulaVar 则直接取值
    if (el == primaryElName) {
      result.finalFractions[el] = lastFraction->getVal();
    } else {
      RooAbsReal* frac_param = (RooAbsReal*)totalPdf->getParameters(*dataHist)->find(Form("frac_%s", el.c_str()));
      result.finalFractions[el] = frac_param ? frac_param->getVal() : 0.0;
    }
  }

    // ***************************************************************
    // ** 扣除逻辑修改：使用模板缩放法 (精确全范围扣除) **
    // ***************************************************************
    
  // Purge: H_pure = H_target - SUM(F_contaminant * h_contam_template * N_total / N_template)
    // 首先克隆目标直方图作为纯净直方图的起点
  result.purifiedHist = unique_ptr<TH1D>((TH1D*)h_target->Clone(Form("h_pure_%s", h_target->GetName())));
  
  if (result.fitStatus <= 1) { // 拟合成功
    double totalFitEvents = dataHist->sumEntries(); // 拟合区间内的总事件数 (N_total)

    // 遍历所有污染物，进行缩放和扣除
    for (const auto& el : contaminantElements_) {
      // 确保拟合结果包含该分数，并且 components 映射中包含原始模板
      if (result.finalFractions.count(el) && components.count(el)) {
        
        // 1. 获取拟合的比例 F_A
        double F_A = result.finalFractions.at(el);
        
        // 2. 获取原始模板 h_template_A (这里是 L2/Pure N/O 等)
        TH1D* h_template_A = components.at(el); 
        // 获取模板的总事件数，用于计算缩放因子
        double N_template_A = h_template_A->GetSumOfWeights(); 
        
        if (N_template_A > 0 && F_A > 0) {
          // 3. 计算目标事件数 N_A = F_A * N_total (在拟合范围内的总事件数)
          double N_target_A = F_A * totalFitEvents;
          
          // 4. 计算缩放因子：S = N_target_A / N_template_A
          // 这个缩放因子会将 h_template_A 的总事件数缩放到 N_target_A
          double scaleFactor = N_target_A / N_template_A;

          // 5. 克隆模板并缩放。克隆是为了不影响原始模板 h_template_A
          auto h_contam_scaled = unique_ptr<TH1D>((TH1D*)h_template_A->Clone(Form("h_contam_scaled_%s_%s", el.c_str(), fitName.c_str())));
          h_contam_scaled->Scale(scaleFactor);
          
          // 6. 执行全范围扣除：h_pure = h_pure - h_contam_scaled
          result.purifiedHist->Add(h_contam_scaled.get(), -1.0);
          
          cout << "   [Purge] Element: " << el 
            << ", Fraction: " << F_A 
            << ", N_template: " << N_template_A
            << ", N_purged: " << N_target_A 
            << ", Scale: " << scaleFactor << endl;
        }
      }
    }
    
    // 7. 清理负值 (可选，但推荐)
    for (int i = 1; i <= result.purifiedHist->GetNbinsX(); ++i) {
      double content = result.purifiedHist->GetBinContent(i);
      if (content < 0) result.purifiedHist->SetBinContent(i, 0.0);
    }

  } else { // 拟合失败，保留原始直方图
    // 确保 purifiedHist 在失败时仍然是 h_target 的克隆
    result.purifiedHist = unique_ptr<TH1D>((TH1D*)h_target->Clone(Form("h_failed_pure_%s", h_target->GetName())));
  }
    // ***************************************************************
    // ** 扣除逻辑修改结束 **
    // ***************************************************************


  // Plotting
  map<string, RooHistPdf*> compPdfs_raw_ptrs;
  for (const auto& el : allComps) {
    if (pdfs.count(el)) {
      compPdfs_raw_ptrs[el] = pdfs.at(el).get();
    }
  }

  double chi2ndf_val = 0.0;
  TGraphErrors* pullGraphPtr = nullptr;
 
  // ** 传递拟合状态到 setupPlot **
  auto frame = setupPlot(fitName, h_target, charge, *dataHist, *totalPdf,
              compPdfs_raw_ptrs,
              chi2ndf_val, h_pure_temp_input,
              pullGraphPtr, result.fitStatus);
 
  unique_ptr<TGraphErrors> pullGraph(pullGraphPtr);
  result.chi2ndf = chi2ndf_val;

  // --- 绘图到传入的 Canvas ---
  c_pdf->Clear();
  c_pdf->Divide(1, 2);

  // pad1
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
 
  // Legend and Info
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
 

  // Pad 2: Pull plot
  c_pdf->cd(2);
  TPad* pad2 = (TPad*)gPad;
  pad2->SetPad(0, 0, 1, 0.3); pad2->SetTopMargin(0.02); pad2->SetBottomMargin(0.3); pad2->SetGridy();
 
  if (pullGraph) {
    // ... (Pull Plot 绘制逻辑保持不变)
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

// --- 主分析函数 CalPureQTemp ---

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

  for (const auto& detector : detectors) {
    cout << "\n--- Processing Detector: " << detector << " ---" << endl;

    // 检查探测器是否有预定义的能量范围
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
      h2d_raw_L1[el]->RebinX(4);
      h2d_raw_L1[el]->Smooth(1);
     
      if (h2dL2_raw) {
        h2d_raw_L2[el] = unique_ptr<TH2F>((TH2F*)h2dL2_raw->Clone(Form("L2_RAW_%s_%s", el.c_str(), detector.c_str())));
        h2d_raw_L2[el]->RebinX(4);
        h2d_raw_L2[el]->Smooth(1);
      } else {
        cerr << "L2 template not found for " << el << " in " << detector << " (Using L1 raw as L2 fallback)" << endl;
        h2d_raw_L2[el] = unique_ptr<TH2F>((TH2F*)h2dL1_raw->Clone(Form("L2_FALLBACK_%s_%s", el.c_str(), detector.c_str())));
        h2d_raw_L2[el]->RebinX(4);
      }
    }
   
    if (h2d_raw_L1.empty()) continue;

    const TAxis* y_axis = h2d_raw_L1.begin()->second->GetYaxis();
    int n_bins_y = y_axis->GetNbins();

    // 按 EkBin 循环
    for (int y_bin = 1; y_bin <= n_bins_y; ++y_bin) {
      double ek_center = y_axis->GetBinCenter(y_bin);
     
      // ** 能量范围检查 **
      if (ek_center < ek_min || ek_center > ek_max) {
        cout << " [Bin " << y_bin << "] Ek=" << ek_center << " GeV/n is OUTSIDE range [" << ek_min << ", " << ek_max << "]. Skipping." << endl;
        continue;
      }
     
      cout << " [Bin " << y_bin << "] Ek=" << ek_center << " GeV/n. Analyzing..." << endl;

      map<string, unique_ptr<TH1D>> currentBinPureTemplates;

      // --- 0. Oxygen (O): 假设纯净 ---
      string elO = "Oxygen";
      if (h2d_raw_L1.count(elO)) {
        auto h_O_raw = unique_ptr<TH1D>(h2d_raw_L1.at(elO)->ProjectionX(Form("h_slice_L1_%s_%s_E%d", elO.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
        currentBinPureTemplates[elO] = unique_ptr<TH1D>((TH1D*)h_O_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elO.c_str(), detector.c_str(), y_bin)));
        finalPureTemplates[elO][y_bin] = unique_ptr<TH1D>((TH1D*)currentBinPureTemplates.at(elO)->Clone());
      }

      // --- 1. Nitrogen (N): 污染 O (Pure) ---
      string elN = "Nitrogen";
      if (h2d_raw_L1.count(elN) && h2d_raw_L2.count(elN) && currentBinPureTemplates.count(elO)) {
        auto h_N_raw = unique_ptr<TH1D>(h2d_raw_L1.at(elN)->ProjectionX(Form("h_slice_L1_%s_%s_E%d", elN.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
        auto h_N_L2 = unique_ptr<TH1D>(h2d_raw_L2.at(elN)->ProjectionX(Form("h_slice_L2_%s_%s_E%d", elN.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
       
        TemplatePurger purger(elN, {elO}, elementZ.at(elN) - 0.6, elementZ.at(elO) + 0.6);
       
        map<string, TH1D*> comps1 = {{elN, h_N_L2.get()}, {elO, currentBinPureTemplates.at(elO).get()}};
        PurgeResult res1 = purger.runFitAndPurge(Form("%s Purge 1 (L2+%s) E%d %s", elN.c_str(), elO.c_str(), y_bin, detector.c_str()), h_N_raw.get(), comps1, {}, nullptr, c_pdf, pdfFileName, -1);
       
        if (res1.fitStatus <= 1) {
          TemplatePurger purger2(elN, {elO}, elementZ.at(elN) - 0.8, elementZ.at(elO) + 0.8);
          map<string, TH1D*> comps2 = {{elN, res1.purifiedHist.get()}, {elO, currentBinPureTemplates.at(elO).get()}};
         
          PurgeResult res2 = purger2.runFitAndPurge(Form("%s Purge 2 (Temp1+%s) E%d %s", elN.c_str(), elO.c_str(), y_bin, detector.c_str()), h_N_raw.get(), comps2, {}, res1.purifiedHist.get(), c_pdf, pdfFileName, res1.fitStatus);
         
          if (res2.fitStatus <= 1) {
            currentBinPureTemplates[elN] = unique_ptr<TH1D>((TH1D*)res2.purifiedHist->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elN.c_str(), detector.c_str(), y_bin)));
            finalPureTemplates[elN][y_bin] = unique_ptr<TH1D>((TH1D*)currentBinPureTemplates.at(elN)->Clone());
          } else {
            currentBinPureTemplates[elN] = unique_ptr<TH1D>((TH1D*)h_N_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elN.c_str(), detector.c_str(), y_bin)));
          }
        } else {
          currentBinPureTemplates[elN] = unique_ptr<TH1D>((TH1D*)h_N_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elN.c_str(), detector.c_str(), y_bin)));
        }
      }

      // --- 2. Carbon (C): 污染 N (Pure), O (Pure)
      string elC = "Carbon";
      if (h2d_raw_L1.count(elC) && h2d_raw_L2.count(elC) && currentBinPureTemplates.count(elN) && currentBinPureTemplates.count(elO)) {
        auto h_C_raw = unique_ptr<TH1D>(h2d_raw_L1.at(elC)->ProjectionX(Form("h_slice_L1_%s_%s_E%d", elC.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
        auto h_C_L2 = unique_ptr<TH1D>(h2d_raw_L2.at(elC)->ProjectionX(Form("h_slice_L2_%s_%s_E%d", elC.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
       
        TemplatePurger purger(elC, {elN, elO}, elementZ.at(elC) - 0.6, elementZ.at(elO) + 0.6);
       
        map<string, TH1D*> comps1 = {{elC, h_C_L2.get()}, {elN, currentBinPureTemplates.at(elN).get()}, {elO, currentBinPureTemplates.at(elO).get()}};
        PurgeResult res1 = purger.runFitAndPurge(Form("%s Purge 1 (L2+PureN+PureO) E%d %s", elC.c_str(), y_bin, detector.c_str()), h_C_raw.get(), comps1, {}, nullptr, c_pdf, pdfFileName, -1);
       
        if (res1.fitStatus <= 1) {
          TemplatePurger purger2(elC, {elN, elO}, elementZ.at(elC) - 0.8, elementZ.at(elO) + 0.8);
          map<string, TH1D*> comps2 = {{elC, res1.purifiedHist.get()}, {elN, currentBinPureTemplates.at(elN).get()}, {elO, currentBinPureTemplates.at(elO).get()}};
          PurgeResult res2 = purger2.runFitAndPurge(Form("%s Purge 2 (Temp1+PureN+PureO) E%d %s", elC.c_str(), y_bin, detector.c_str()), h_C_raw.get(), comps2, {}, res1.purifiedHist.get(), c_pdf, pdfFileName, res1.fitStatus);
         
          if (res2.fitStatus <= 1) {
            currentBinPureTemplates[elC] = unique_ptr<TH1D>((TH1D*)res2.purifiedHist->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elC.c_str(), detector.c_str(), y_bin)));
            finalPureTemplates[elC][y_bin] = unique_ptr<TH1D>((TH1D*)currentBinPureTemplates.at(elC)->Clone());
          } else {
            currentBinPureTemplates[elC] = unique_ptr<TH1D>((TH1D*)h_C_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elC.c_str(), detector.c_str(), y_bin)));
          }
        } else {
          currentBinPureTemplates[elC] = unique_ptr<TH1D>((TH1D*)h_C_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elC.c_str(), detector.c_str(), y_bin)));
        }
      }
     
      // --- 3. Boron (B): 污染 C (Pure), N (Pure)
      string elB = "Boron";
      if (h2d_raw_L1.count(elB) && h2d_raw_L2.count(elB) && currentBinPureTemplates.count(elC) && currentBinPureTemplates.count(elN)) {
        auto h_B_raw = unique_ptr<TH1D>(h2d_raw_L1.at(elB)->ProjectionX(Form("h_slice_L1_%s_%s_E%d", elB.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
        auto h_B_L2 = unique_ptr<TH1D>(h2d_raw_L2.at(elB)->ProjectionX(Form("h_slice_L2_%s_%s_E%d", elB.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
       
        TemplatePurger purger(elB, {elC, elN}, elementZ.at(elB) - 0.6, elementZ.at(elC) + 0.6);
       
        map<string, TH1D*> comps1 = {{elB, h_B_L2.get()}, {elC, currentBinPureTemplates.at(elC).get()}, {elN, currentBinPureTemplates.at(elN).get()}};
        PurgeResult res1 = purger.runFitAndPurge(Form("%s Purge 1 (L2+PureC+PureN) E%d %s", elB.c_str(), y_bin, detector.c_str()), h_B_raw.get(), comps1, {}, nullptr, c_pdf, pdfFileName, -1);
       
        if (res1.fitStatus <= 1) {
          TemplatePurger purger2(elB, {elC, elN}, elementZ.at(elB) - 0.8, elementZ.at(elC) + 0.8);
          map<string, TH1D*> comps2 = {{elB, res1.purifiedHist.get()}, {elC, currentBinPureTemplates.at(elC).get()}, {elN, currentBinPureTemplates.at(elN).get()}};
          PurgeResult res2 = purger2.runFitAndPurge(Form("%s Purge 2 (Temp1+PureC+PureN) E%d %s", elB.c_str(), y_bin, detector.c_str()), h_B_raw.get(), comps2, {}, res1.purifiedHist.get(), c_pdf, pdfFileName, res1.fitStatus);
         
          if (res2.fitStatus <= 1) {
            currentBinPureTemplates[elB] = unique_ptr<TH1D>((TH1D*)res2.purifiedHist->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elB.c_str(), detector.c_str(), y_bin)));
            finalPureTemplates[elB][y_bin] = unique_ptr<TH1D>((TH1D*)currentBinPureTemplates.at(elB)->Clone());
          } else {
            currentBinPureTemplates[elB] = unique_ptr<TH1D>((TH1D*)h_B_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elB.c_str(), detector.c_str(), y_bin)));
          }
        } else {
          currentBinPureTemplates[elB] = unique_ptr<TH1D>((TH1D*)h_B_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elB.c_str(), detector.c_str(), y_bin)));
        }
      }

      // --- 4. Beryllium (Be): 污染 B (Pure), C (Pure)
      string elBe = "Beryllium";
      if (h2d_raw_L1.count(elBe) && h2d_raw_L2.count(elBe) && currentBinPureTemplates.count(elB) && currentBinPureTemplates.count(elC)) {
        auto h_Be_raw = unique_ptr<TH1D>(h2d_raw_L1.at(elBe)->ProjectionX(Form("h_slice_L1_%s_%s_E%d", elBe.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
        auto h_Be_L2 = unique_ptr<TH1D>(h2d_raw_L2.at(elBe)->ProjectionX(Form("h_slice_L2_%s_%s_E%d", elBe.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e"));
       
        TemplatePurger purger(elBe, {elB, elC}, elementZ.at(elBe) - 0.6, elementZ.at(elC) + 0.6);
       
        map<string, TH1D*> comps1 = {{elBe, h_Be_L2.get()}, {elB, currentBinPureTemplates.at(elB).get()}, {elC, currentBinPureTemplates.at(elC).get()}};
        PurgeResult res1 = purger.runFitAndPurge(Form("%s Purge 1 (L2+PureB+PureC) E%d %s", elBe.c_str(), y_bin, detector.c_str()), h_Be_raw.get(), comps1, {}, nullptr, c_pdf, pdfFileName, -1);
       
        if (res1.fitStatus <= 1) {
          TemplatePurger purger2(elBe, {elB, elC}, elementZ.at(elBe) - 0.8, elementZ.at(elC) + 0.8);
          map<string, TH1D*> comps2 = {{elBe, res1.purifiedHist.get()}, {elB, currentBinPureTemplates.at(elB).get()}, {elC, currentBinPureTemplates.at(elC).get()}};
          PurgeResult res2 = purger2.runFitAndPurge(Form("%s Purge 2 (Temp1+PureB+PureC) E%d %s", elBe.c_str(), y_bin, detector.c_str()), h_Be_raw.get(), comps2, {}, res1.purifiedHist.get(), c_pdf, pdfFileName, res1.fitStatus);
         
          if (res2.fitStatus <= 1) {
            currentBinPureTemplates[elBe] = unique_ptr<TH1D>((TH1D*)res2.purifiedHist->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elBe.c_str(), detector.c_str(), y_bin)));
            finalPureTemplates[elBe][y_bin] = unique_ptr<TH1D>((TH1D*)currentBinPureTemplates.at(elBe)->Clone());
          } else {
            currentBinPureTemplates[elBe] = unique_ptr<TH1D>((TH1D*)h_Be_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elBe.c_str(), detector.c_str(), y_bin)));
          }
        } else {
          // ** 修正点 3：使用 currentBinPureTemplates **
          currentBinPureTemplates[elBe] = unique_ptr<TH1D>((TH1D*)h_Be_raw->Clone(Form("L1QTemp_Pure_%s_%s_E%d", elBe.c_str(), detector.c_str(), y_bin)));
        }
      }
    }
   
    // 5. 写入 ROOT 文件
    for (const auto& el : requiredElements) {
      outputFile->cd();
      if (!finalPureTemplates.count(el) || !h2d_raw_L1.count(el)) continue;

      const auto& h_base_2d = h2d_raw_L1.at(el);
        std::string pureHistName = Form("h2d_PureQTemp_%s_%s", el.c_str(), detector.c_str());
        auto h_pure_2d = unique_ptr<TH2F>((TH2F*)h_base_2d->Clone(pureHistName.c_str()));
            
        h_pure_2d->Reset();

      for (int y_bin = 1; y_bin <= n_bins_y; ++y_bin) {
        // 确保只有在分析过的 bin 中才写入数据
        if (finalPureTemplates.at(el).count(y_bin)) {
          TH1D* h_slice = finalPureTemplates.at(el).at(y_bin).get();
          for (int x_bin = 1; x_bin <= h_slice->GetNbinsX(); ++x_bin) {
            h_pure_2d->SetBinContent(x_bin, y_bin, h_slice->GetBinContent(x_bin));
          }
        }
      }
      h_pure_2d->Write();
    }
  }

  c_pdf->Print((pdfFileName + "]").c_str());
}