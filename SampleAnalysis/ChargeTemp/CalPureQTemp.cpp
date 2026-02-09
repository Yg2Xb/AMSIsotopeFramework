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
#include <TH2D.h>
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

const std::string inputFileName = "/eos/user/z/zixuan/Isotope/Add/Be_frag4_NoBkg_Tune_full.root";
const std::string outputDir = "/eos/user/z/zixuan/Isotope/PureChargeTemp/";
const std::string chainName = "UnbiasedL1Inner";
const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
const double Q_GLOBAL_MIN = 1.0;
const double Q_GLOBAL_MAX = 9.0;

const std::vector<std::string> requiredElements = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
const std::map<std::string, int> elementZ = {{"Helium", 2}, {"Lithium", 3}, {"Beryllium", 4}, {"Boron", 5}, {"Carbon", 6}, {"Nitrogen", 7}, {"Oxygen", 8}};

const std::map<std::string, std::pair<double, double>> detector_ek_ranges = {
	{"TOF", {0.25, 1.5}},
	{"NaF", {0.61, 6.10}},
	{"AGL", {2.50, 25.0}}
};

std::string getTemplateHistName(const std::string& elName, const std::string& detector, std::string tag) {
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
	double yMinPlot = 0.01;

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
			totalPdf->fitTo(*dataHist, 
				Save(true), 
				PrintLevel(-1), 
				Range(fitMin_, fitMax_), 
				Strategy(0), 
				Minimizer("Minuit2", "migrad"),
				Offset(kTRUE),
				SumW2Error(kTRUE),
				NumCPU(4))
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
		double binWidth = h_target->GetBinWidth(1);
		RooArgSet normSet(charge);

		// 2. Smooth PDF Subtraction (Initial)
		for (int i = 1; i <= result.purifiedHist->GetNbinsX(); ++i) {
			double x = result.purifiedHist->GetBinCenter(i);
			if (x < fitMin_ || x > fitMax_) continue;

			double originalContent = h_target->GetBinContent(i);
			double totalContamExpected = 0;

			charge.setVal(x);
			for (const auto& el : contaminantElements_) {
				if (pdfs.count(el)) {
					double pdfVal = pdfs[el]->getVal(normSet);
					totalContamExpected += result.finalFractions[el] * totalFitEvents * pdfVal * binWidth;
				}
			}
			result.purifiedHist->SetBinContent(i, originalContent - totalContamExpected);
			result.purifiedHist->SetBinError(i, h_target->GetBinError(i));
		}

		// 3. Monotonic Back-tracing Correction using L2 Relative Ratio Anomaly Detection
		TH1D* h_ref = h_pure_temp_input ? h_pure_temp_input : components.at(primaryElement_);
		double maxRef = h_ref->GetMaximum();
		double maxL1 = result.purifiedHist->GetMaximum();
		double startQ = elementZ.at(primaryElement_) + 0.6;
		int startBin = result.purifiedHist->FindBin(startQ);
		int nBins = result.purifiedHist->GetNbinsX();

		for (int i = startBin; i <= nBins; ++i) {
			double cCurr = result.purifiedHist->GetBinContent(i);
			double cPrev = result.purifiedHist->GetBinContent(i-1);
			
			double rL1 = (maxL1 > 0) ? (cCurr / maxL1) : 0;
			double rRef = (maxRef > 0) ? (h_ref->GetBinContent(i) / maxRef) : 0;

			bool ratioAnomaly = (rRef > 0 && (rL1 > 2.0 * rRef || rL1 < 0.5 * rRef));
			bool monotonicAnomaly = (cCurr > 0 && (cCurr - sqrt(cCurr)) > cPrev);
			bool negativeAnomaly = (cCurr <= 0.0);

			if (ratioAnomaly || monotonicAnomaly || negativeAnomaly) {
				double refPrev2 = h_ref->GetBinContent(i-2);
				double refPrev = h_ref->GetBinContent(i-1);
				if (refPrev2 > 0) {
					double ratio1 = refPrev / refPrev2;
					result.purifiedHist->SetBinContent(i-1, result.purifiedHist->GetBinContent(i-2) * ratio1);
				}
				double refCurr = h_ref->GetBinContent(i);
				double refPrevNow = h_ref->GetBinContent(i-1);
				if (refPrevNow > 0) {
					double ratio2 = refCurr / refPrevNow;
					result.purifiedHist->SetBinContent(i, result.purifiedHist->GetBinContent(i-1) * ratio2);
				} else {
					result.purifiedHist->SetBinContent(i, 0.0);
				}
			}
		}

		int maxBinRef = h_ref->GetMaximumBin();
		for (int i = maxBinRef; i <= nBins; ++i) {
			if (h_ref->GetBinContent(i) <= 0 && (i == nBins || h_ref->GetBinContent(i+1) <= 0)) {
				for (int j = i; j <= nBins; ++j) {
					result.purifiedHist->SetBinContent(j, 0.0);
					result.purifiedHist->SetBinError(j, 0.0);
				}
				break;
			}
		}
		
		for (int i = 1; i <= nBins; ++i) if (result.purifiedHist->GetBinContent(i) < 0) result.purifiedHist->SetBinContent(i, 0);

	} else {
		result.purifiedHist = unique_ptr<TH1D>((TH1D*)h_target->Clone(Form("h_failed_pure_%s", h_target->GetName())));
	}

	map<string, RooHistPdf*> compPdfs_raw_ptrs;
	for (const auto& el : allComps) if (pdfs.count(el)) compPdfs_raw_ptrs[el] = pdfs.at(el).get();

	double chi2ndf_val = 0.0;
	TGraphErrors* pullGraphPtr = nullptr;
	auto frame = setupPlot(fitName, h_target, charge, *dataHist, *totalPdf, compPdfs_raw_ptrs, chi2ndf_val, h_pure_temp_input, pullGraphPtr, result.fitStatus);
	unique_ptr<TGraphErrors> pullGraph(pullGraphPtr);
	result.chi2ndf = chi2ndf_val;

	c_pdf->Clear();
	c_pdf->Divide(1, 2);

	c_pdf->cd(1);
	TPad* pad1 = (TPad*)gPad;
	pad1->SetPad(0, 0.3, 1, 1);
	h_target->SetMinimum(0.01);
	h_target->GetXaxis()->SetRangeUser(fitMin_, fitMax_);
	h_target->Draw("PZ");
	frame->SetMinimum(0.01);

	pad1->SetLogy();
	pad1->SetBottomMargin(0.02);

	pad1->cd();

	frame->GetYaxis()->SetTitle("Events");
	frame->GetXaxis()->SetLabelSize(0);
	frame->Draw();

	auto legend = std::make_unique<TLegend>(0.7, 0.55, 0.88, 0.88);
	legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.03);
	legend->AddEntry("data_hist", Form("%s L1 Raw", primaryElement_.c_str()), "pe");
	if (result.fitStatus <= 1) legend->AddEntry("total_pdf", "Total Fit", "l");
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
	pad2->SetPad(0, 0, 1, 0.3); pad2->SetBottomMargin(0.3); pad2->SetGridy();
	if (pullGraph) {
		pullGraph->SetMarkerStyle(20);
		pullGraph->Draw("AP");
		pullGraph->GetYaxis()->SetRangeUser(-6, 6);
		TLine line;
		line.SetLineColor(kRed);
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

	std::string pdfFileName = outputDir + "NoBkg_CalPureQTemp_" + chainName + ".pdf";
	std::string rootFileName = outputDir + "NoBkg_PureChargeTemplates_" + chainName + ".root";
	system(Form("mkdir -p %s", outputDir.c_str()));

	auto inputFile = std::unique_ptr<TFile>(TFile::Open(inputFileName.c_str()));
	if (!inputFile || inputFile->IsZombie()) return;
	auto outputFile = std::make_unique<TFile>(rootFileName.c_str(), "RECREATE");

	TCanvas* c_pdf = new TCanvas("c_pdf", "PDF Canvas", 800, 600);
	c_pdf->Print((pdfFileName + "[").c_str());

	const std::vector<std::string> orderedElements = {"Oxygen", "Nitrogen", "Carbon", "Boron", "Beryllium", "Lithium", "Helium"};
	const std::map<std::string, std::vector<std::string>> contaminationMap = { {"Oxygen", {}}, {"Nitrogen", {"Oxygen"}}, {"Carbon", {"Nitrogen", "Oxygen"}}, {"Boron", {"Carbon", "Nitrogen", "Oxygen"}}, {"Beryllium", {"Boron", "Carbon", "Nitrogen", "Oxygen"}}, {"Lithium", {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"}}, {"Helium", {}} };

	for (const auto& detector : detectors) {
		if (detector_ek_ranges.find(detector) == detector_ek_ranges.end()) continue;
		double ek_min = detector_ek_ranges.at(detector).first;
		double ek_max = detector_ek_ranges.at(detector).second;

		map<string, map<int, unique_ptr<TH1D>>> finalPureTemplates;
		map<string, unique_ptr<TH2F>> h2d_raw_L1, h2d_raw_L2, h2d_raw_Inner;
		for (const auto& el : requiredElements) {
			string nameL1 = getTemplateHistName(el, detector, "L1Template");
			string nameL2 = getTemplateHistName(el, detector, "L2Template");
			string nameInner = getTemplateHistName(el, detector, "InnerQTemplate"); 
			TH2F* h2dL1_raw = (TH2F*)inputFile->Get(nameL1.c_str());
			TH2F* h2dL2_raw = (TH2F*)inputFile->Get(nameL2.c_str());
			TH2F* h2dInner_raw = (TH2F*)inputFile->Get(nameInner.c_str()); 
			if (h2dL1_raw) { h2d_raw_L1[el] = unique_ptr<TH2F>((TH2F*)h2dL1_raw->Clone(Form("L1_RAW_%s_%s", el.c_str(), detector.c_str()))); h2d_raw_L1[el]->RebinX(2); }
			if (h2dL2_raw) { h2d_raw_L2[el] = unique_ptr<TH2F>((TH2F*)h2dL2_raw->Clone(Form("L2_RAW_%s_%s", el.c_str(), detector.c_str()))); h2d_raw_L2[el]->RebinX(2); }
			if (h2dInner_raw) { h2d_raw_Inner[el] = unique_ptr<TH2F>((TH2F*)h2dInner_raw->Clone(Form("Inner_RAW_%s_%s", el.c_str(), detector.c_str()))); h2d_raw_Inner[el]->RebinX(2); }
		}

		if (h2d_raw_L1.empty()) continue;
		const TAxis* y_axis = h2d_raw_L1.begin()->second->GetYaxis();

		for (int y_bin = 1; y_bin <= y_axis->GetNbins(); ++y_bin) {
			if (y_axis->GetBinCenter(y_bin) < ek_min || y_axis->GetBinCenter(y_bin) > ek_max) continue;
			map<string, unique_ptr<TH1D>> currentBinPureTemplates;
			for (const auto& elPrimary : orderedElements) {
				if (!h2d_raw_L1.count(elPrimary)) continue;
				if (elPrimary == "Oxygen" || elPrimary == "Helium") {
					auto h = unique_ptr<TH1D>(h2d_raw_L1.at(elPrimary)->ProjectionX(Form("h_slice_L1_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e")); h->Smooth(1,"G");
					finalPureTemplates[elPrimary][y_bin] = unique_ptr<TH1D>((TH1D*)h->Clone()); currentBinPureTemplates[elPrimary] = unique_ptr<TH1D>((TH1D*)h->Clone()); continue;
				}
				map<string, TH1D*> comps1; for (const auto& contam : contaminationMap.at(elPrimary)) if (currentBinPureTemplates.count(contam)) comps1[contam] = currentBinPureTemplates.at(contam).get();
				auto h_Primary_raw = unique_ptr<TH1D>(h2d_raw_L1.at(elPrimary)->ProjectionX(Form("h_slice_L1_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e")); h_Primary_raw->Smooth(1,"G");
				auto h_Primary_L2 = unique_ptr<TH1D>(h2d_raw_L2.at(elPrimary)->ProjectionX(Form("h_slice_L2_%s_%s_E%d", elPrimary.c_str(), detector.c_str(), y_bin), y_bin, y_bin, "e")); h_Primary_L2->Smooth(1,"G");
				comps1[elPrimary] = h_Primary_L2.get();
				TemplatePurger purger1(elPrimary, contaminationMap.at(elPrimary), elementZ.at(elPrimary)-0.4, Q_GLOBAL_MAX-0.6);
				PurgeResult res1 = purger1.runFitAndPurge(Form("%s Purge 1 Ek/n=%.2f %s", elPrimary.c_str(), y_axis->GetBinCenter(y_bin), detector.c_str()), h_Primary_raw.get(), comps1, {}, nullptr, c_pdf, pdfFileName, -1);
				if (res1.fitStatus <= 1) {
					TemplatePurger purger2(elPrimary, contaminationMap.at(elPrimary), elementZ.at(elPrimary)-0.4, Q_GLOBAL_MAX-0.6);
					map<string, TH1D*> comps2 = comps1; comps2[elPrimary] = res1.purifiedHist.get();
					PurgeResult res2 = purger2.runFitAndPurge(Form("%s Purge 2 Ek/n=%.2f %s", elPrimary.c_str(), y_axis->GetBinCenter(y_bin), detector.c_str()), h_Primary_raw.get(), comps2, {}, res1.purifiedHist.get(), c_pdf, pdfFileName, res1.fitStatus);
					if (res2.fitStatus <= 1) { currentBinPureTemplates[elPrimary] = unique_ptr<TH1D>((TH1D*)res2.purifiedHist->Clone()); finalPureTemplates[elPrimary][y_bin] = unique_ptr<TH1D>((TH1D*)currentBinPureTemplates.at(elPrimary)->Clone()); }
				}
			}
		}

		for (const auto& el : requiredElements) {
			outputFile->cd(); 
			if (finalPureTemplates[el].empty()) continue;

			TH2F* h_pure_2d = (TH2F*)h2d_raw_L1[el]->Clone(Form("h2d_PureQTemp_%s_%s", el.c_str(), detector.c_str())); 
			h_pure_2d->Reset();
			for (int y_bin = 1; y_bin <= y_axis->GetNbins(); ++y_bin) {
				if (finalPureTemplates.at(el).count(y_bin)) { 
					TH1D* h_slice = finalPureTemplates.at(el).at(y_bin).get(); 
					for (int x_bin = 1; x_bin <= h_slice->GetNbinsX(); ++x_bin) h_pure_2d->SetBinContent(x_bin, y_bin, h_slice->GetBinContent(x_bin)); 
				}
			}
			h_pure_2d->Write();

			if (h2d_raw_L2.count(el)) {
				TH2D h_tune(Form("h2d_TuneQTemp_%s_%s", el.c_str(), detector.c_str()), "", 
					h2d_raw_L2[el]->GetNbinsX(), h2d_raw_L2[el]->GetXaxis()->GetXmin(), h2d_raw_L2[el]->GetXaxis()->GetXmax(), 
					h2d_raw_L2[el]->GetNbinsY(), h2d_raw_L2[el]->GetYaxis()->GetXbins()->GetArray());
				for(int x=1; x<=h2d_raw_L2[el]->GetNbinsX(); ++x) {
					for(int y=1; y<=h2d_raw_L2[el]->GetNbinsY(); ++y) {
						h_tune.SetBinContent(x, y, h2d_raw_L2[el]->GetBinContent(x, y));
						h_tune.SetBinError(x, y, h2d_raw_L2[el]->GetBinError(x, y));
					}
				}
				h_tune.Write();
			}

			if (h2d_raw_Inner.count(el)) { 
				TH2D h_inner(Form("h2d_InnerQTemp_%s_%s", el.c_str(), detector.c_str()), "", 
					h2d_raw_Inner[el]->GetNbinsX(), h2d_raw_Inner[el]->GetXaxis()->GetXmin(), h2d_raw_Inner[el]->GetXaxis()->GetXmax(), 
					h2d_raw_Inner[el]->GetNbinsY(), h2d_raw_Inner[el]->GetYaxis()->GetXbins()->GetArray()); 
				for(int x=1; x<=h2d_raw_Inner[el]->GetNbinsX();++x) {
					for(int y=1; y<=h2d_raw_Inner[el]->GetNbinsY();++y) { 
						h_inner.SetBinContent(x,y,h2d_raw_Inner[el]->GetBinContent(x,y)); 
						h_inner.SetBinError(x,y,h2d_raw_Inner[el]->GetBinError(x,y)); 
					} 
				}
				h_inner.Write(); 
			}
		}
	}
	c_pdf->Print((pdfFileName + "]").c_str()); 
	outputFile->Close();
}