#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <array>
#include <map>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TPaveText.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooArgList.h>
#include <RooMsgService.h>
#include <RooFormulaVar.h>

using namespace RooFit;

// Configuration
const std::string inputFileName = "/eos/user/z/zixuan/Isotope/Add/Be_all.root";
const std::string outputDir = "/eos/user/z/zixuan/Isotope/chargeTempFit/";

const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
const std::vector<std::string> elements = {"Carbon", "Nitrogen", "Oxygen"};
const std::vector<std::string> chargeTypes = {"UnbiasedL1Inner", "L1Inner"};

// Template configuration for each element
struct TemplateConfig {
	std::vector<std::string> templates;
	double fit_min, fit_max;
	double frac_min, frac_max;
};

const std::map<std::string, TemplateConfig> elementConfig = {
	{"Carbon", {{"L2QTemplate", "L1QTemplate", "L2QTemplate"}, 5.3, 6.7, 5.6, 6.4}},  // C: BL2, CL1, NL2
	{"Nitrogen", {{"L1QTemplate", "L2QTemplate", "L1QTemplate"}, 6.3, 7.7, 6.8, 7.4}}, // N: CL1, NL2, OL1  
	{"Oxygen", {{"L2QTemplate", "L1QTemplate"}, 7.3, 8.7, 7.6, 8.4}}                   // O: NL2, OL1
};

// Fraction calculation ranges
const std::map<std::string, std::pair<double, double>> fractionRanges = {
	{"Carbon", {5.6, 6.4}},    // Z+-0.4
	{"Nitrogen", {6.8, 7.4}},  // Z-0.2 to Z+0.4  
	{"Oxygen", {7.6, 8.4}}     // Z+-0.4
};

struct FractionResult {
	std::vector<double> fractions;
	std::vector<double> errors;
	FractionResult(int n_elements) : fractions(n_elements, 0.0), errors(n_elements, 0.0) {}
};

struct EnergyBinResult {
	double energy;
	double energy_width;
	FractionResult fit_fractions;
	FractionResult narrow_fractions;
	double chi2_ndf;
	int ndf;
	int fit_status;
	double n_signal;
	bool no_fit_needed;  // 新增：标记是否需要拟合

	EnergyBinResult(int n_elements) : fit_fractions(n_elements), narrow_fractions(n_elements), no_fit_needed(false) {}
};

// Helper function to create legend
std::unique_ptr<TLegend> createLegend() {
	auto legend = std::make_unique<TLegend>(0.7, 0.7, 0.95, 0.9);
	legend->SetTextSize(0.035);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	return legend;
}

// Helper function to create info box
std::unique_ptr<TPaveText> createInfoBox(
		const std::string& element,
		double energy,
		const EnergyBinResult& result
		) {
	auto info = std::make_unique<TPaveText>(0.15, 0.65, 0.65, 0.9, "NDC");
	info->SetFillStyle(0);
	info->SetBorderSize(0);
	info->SetTextSize(0.03);
	info->SetTextAlign(12);

	info->AddText(Form("Element: %s", element.c_str()));
	info->AddText(Form("Energy: %.2f GeV/n", energy));

	if (result.no_fit_needed) {
		info->AddText("No fit needed - insufficient template coverage");
		info->AddText("Default fractions applied");
	} else if (result.fit_status == 0) {
		info->AddText(Form("#chi^{2}/ndf = %.2f/%d = %.2f", result.chi2_ndf*result.ndf, result.ndf, result.chi2_ndf));

		// Display fit fractions
		std::string fit_text = "Fit fractions: ";
		for (size_t i = 0; i < result.fit_fractions.fractions.size(); ++i) {
			if (i > 0) fit_text += ", ";
			fit_text += Form("%.3f±%.3f", result.fit_fractions.fractions[i], result.fit_fractions.errors[i]);
		}
		info->AddText(fit_text.c_str());

		// Display narrow range fractions
		if (result.n_signal > 0) {
			std::string narrow_text = "Narrow fractions: ";
			for (size_t i = 0; i < result.narrow_fractions.fractions.size(); ++i) {
				if (i > 0) narrow_text += ", ";
				narrow_text += Form("%.4f±%.4f", result.narrow_fractions.fractions[i], result.narrow_fractions.errors[i]);
			}
			info->AddText(narrow_text.c_str());
		}
	} else {
		info->AddText("Fit Failed");
	}
	return info;
}

// Chi2 calculation function
double calculateChi2(RooPlot* frame, const char* dataName, const char* pdfName, double xmin, double xmax) {
	double chi2 = 0.0;
	int npoints = 0;

	RooCurve* pdfCurve = frame->getCurve(pdfName);
	RooHist* dataHist = frame->getHist(dataName);

	if (!pdfCurve || !dataHist) return 0.0;

	for (int i = 0; i < dataHist->GetN(); ++i) {
		double x, y, ex, ey;
		dataHist->GetPoint(i, x, y);
		ex = dataHist->GetErrorXhigh(i);
		ey = dataHist->GetErrorYhigh(i);

		if (x < xmin || x > xmax) continue;
		if (ey <= 0) continue;

		double pdf_val = pdfCurve->Eval(x);
		chi2 += pow((y - pdf_val) / ey, 2);
		npoints++;
	}

	return chi2;
}

// Pull calculation
void calculatePull(RooPlot* frame, TGraphErrors* pullGraph, double xmin, double xmax) {
	RooCurve* pdfCurve = frame->getCurve("total_pdf");
	RooHist* dataHist = frame->getHist("data_hist");

	if (!pdfCurve || !dataHist) return;

	int point = 0;
	for (int i = 0; i < dataHist->GetN(); ++i) {
		double x, y, ex, ey;
		dataHist->GetPoint(i, x, y);
		ex = dataHist->GetErrorXhigh(i);
		ey = dataHist->GetErrorYhigh(i);

		if (x < xmin || x > xmax || ey <= 0) continue;

		double pdf_val = pdfCurve->Eval(x);
		double pull = (y - pdf_val) / ey;

		pullGraph->SetPoint(point, x, pull);
		pullGraph->SetPointError(point, ex, 1.0);
		point++;
	}
}

// Setup pull plot
void setupPullPlot(TGraphErrors* pullGraph, double xmin, double xmax) {
	pullGraph->SetTitle("");
	pullGraph->GetXaxis()->SetTitle("Charge");
	pullGraph->GetYaxis()->SetTitle("Pull");
	pullGraph->GetXaxis()->SetRangeUser(xmin, xmax);
	pullGraph->GetYaxis()->SetRangeUser(-4, 4);
	pullGraph->SetMarkerStyle(20);
	pullGraph->SetMarkerSize(0.6);
	pullGraph->GetXaxis()->SetTitleSize(0.12);
	pullGraph->GetYaxis()->SetTitleSize(0.12);
	pullGraph->GetXaxis()->SetLabelSize(0.1);
	pullGraph->GetYaxis()->SetLabelSize(0.1);
	pullGraph->GetXaxis()->SetTitleOffset(1.2);
	pullGraph->GetYaxis()->SetTitleOffset(0.4);
}

// Get 1D charge slice from 2D histogram for specific energy bin
std::unique_ptr<TH1D> getChargeSlice(TH2F* h2d, int energy_bin, const char* name) {
	if (!h2d) return nullptr;

	auto slice = std::make_unique<TH1D>(name, name, 
			h2d->GetNbinsX(), h2d->GetXaxis()->GetXmin(), h2d->GetXaxis()->GetXmax());

	for (int i = 1; i <= h2d->GetNbinsX(); ++i) {
		double content = h2d->GetBinContent(i, energy_bin);
		double error = h2d->GetBinError(i, energy_bin);
		slice->SetBinContent(i, content);
		slice->SetBinError(i, error);
	}

	return slice;
}

// Check if template has sufficient coverage in fit range
bool checkTemplateCoverage(TH1D* hist, double fit_min, double fit_max) {
	if (!hist) return false;
	
	int first_bin = hist->FindBin(fit_min);
	int last_bin = hist->FindBin(fit_max);
	double total_events = hist->Integral(first_bin, last_bin);
	
	return total_events > 1.0;  // At least 1 event in fit range
}

// Main fit processor class for single energy bin
class EnergyBinFitProcessor {
	public:
		EnergyBinFitProcessor(const std::string& element, const std::string& detector, 
				const std::string& chargeType, TFile* inputFile, int energyBin, 
				TH2F* rebinnedSignal);
		~EnergyBinFitProcessor() = default;

		bool initialize();
		std::pair<RooPlot*, EnergyBinResult> runFit();

	private:
		std::string elementName, detectorName, chargeTypeName;
		TFile* inputFilePtr;
		int energyBinIndex;
		TH2F* h2d_signal_rebinned; 

		std::unique_ptr<TH1D> h_signal;
		std::vector<std::unique_ptr<TH1D>> template_hists;
		std::vector<std::string> templateElements;
		std::vector<std::unique_ptr<TH1D>> extended_template_hists;

		RooRealVar charge;
		std::vector<std::unique_ptr<RooRealVar>> fractionParams;
		std::vector<RooRealVar*> fractionVars;
		std::unique_ptr<RooFormulaVar> lastFraction;

		std::vector<std::unique_ptr<RooDataHist>> template_data_hists;
		std::vector<std::unique_ptr<RooHistPdf>> template_pdfs;

		std::unique_ptr<RooAddPdf> total_pdf;
		std::unique_ptr<RooDataHist> data_hist;

		bool has_sufficient_templates;  // 新增：标记是否有足够的模板覆盖

		void setupTemplateElements();
		FractionResult calculateNarrowRangeFractions();
		EnergyBinResult createDefaultResult();  // 新增：创建默认结果
};

EnergyBinFitProcessor::EnergyBinFitProcessor(const std::string& element, const std::string& detector, 
                      const std::string& chargeType, TFile* inputFile, int energyBin, 
                      TH2F* rebinnedSignal)
	: elementName(element), detectorName(detector), chargeTypeName(chargeType), 
	inputFilePtr(inputFile), energyBinIndex(energyBin), h2d_signal_rebinned(rebinnedSignal),
	charge("charge", "Charge", elementConfig.at(element).fit_min, elementConfig.at(element).fit_max),
	has_sufficient_templates(false)
{
	setupTemplateElements();
}

void EnergyBinFitProcessor::setupTemplateElements() {
	if (elementName == "Carbon") {
		templateElements = {"Beryllium", "Carbon", "Nitrogen"};  // BL2, CL1, NL2
	} else if (elementName == "Nitrogen") {
		templateElements = {"Carbon", "Nitrogen", "Oxygen"};     // CL1, NL2, OL1
	} else if (elementName == "Oxygen") {
		templateElements = {"Nitrogen", "Oxygen"};               // NL2, OL1
	}

	// Create fraction variables (n-1 free parameters)
	for (size_t i = 0; i < templateElements.size() - 1; ++i) {
		auto fracVar = std::make_unique<RooRealVar>(
				Form("frac_%s_bin%d", templateElements[i].c_str(), energyBinIndex),
				Form("%s fraction", templateElements[i].c_str()),
				1.0 / templateElements.size(), 0.0, 1.0
				);
		fractionVars.push_back(fracVar.get());
		fractionParams.push_back(std::move(fracVar));
	}

	// Last fraction constrained to sum to 1
	if (templateElements.size() > 1) {
		std::string formula = "1.0";
		RooArgList argList;
		for (auto& fvar : fractionVars) {
			formula += " - " + std::string(fvar->GetName());
			argList.add(*fvar);
		}

		lastFraction = std::make_unique<RooFormulaVar>(
				Form("frac_%s_bin%d", templateElements.back().c_str(), energyBinIndex),
				Form("%s fraction", templateElements.back().c_str()),
				formula.c_str(), argList
				);
	}
}

bool EnergyBinFitProcessor::initialize() {
	// Get charge slice for this energy bin from rebinned signal
	h_signal = getChargeSlice(h2d_signal_rebinned, energyBinIndex, 
		Form("signal_energy_bin_%d", energyBinIndex));
	if (!h_signal || h_signal->GetEntries() < 20) {
		std::cerr << "Insufficient signal events in energy bin " << energyBinIndex << std::endl;
		return false;
	}

	// Load template histograms
	auto config = elementConfig.at(elementName);
	int valid_templates = 0;
	
	for (size_t i = 0; i < templateElements.size(); ++i) {
		std::string templateName = chargeTypeName + "_ISS_BKG_H2_" + templateElements[i] + "_" + 
			config.templates[i] + "_" + detectorName;

		TH2F* h2d_template = (TH2F*)inputFilePtr->Get(templateName.c_str());
		if (!h2d_template) {
			std::cerr << "Cannot find template histogram: " << templateName << std::endl;
			return false;
		}

		TH2F* h2d_template_rebinned = (TH2F*)h2d_template->Clone(Form("%s_rebinned", h2d_template->GetName()));
		h2d_template_rebinned->RebinY(2);

		auto h_template = getChargeSlice(h2d_template_rebinned, energyBinIndex, 
				Form("%s_template_energy_bin_%d", templateElements[i].c_str(), energyBinIndex));
		
		if (!h_template || h_template->GetEntries() < 5) {
			std::cerr << "Insufficient template events for " << templateElements[i] 
				<< " in energy bin " << energyBinIndex << std::endl;
			return false;
		}

		// Check if template has coverage in fit range
		if (checkTemplateCoverage(h_template.get(), charge.getMin(), charge.getMax())) {
			valid_templates++;
			std::cout << "        " << templateElements[i] << " has coverage in fit range" << std::endl;
		} else {
			std::cout << "        " << templateElements[i] << " has NO coverage in fit range" << std::endl;
		}

		template_hists.push_back(std::move(h_template));
		delete h2d_template_rebinned;
	}

	// Check if we have sufficient template coverage
	has_sufficient_templates = (valid_templates >= 2);  // Need at least 2 templates with coverage
	
	if (!has_sufficient_templates) {
		std::cout << "        Insufficient template coverage for fitting" << std::endl;
		return true;  // Still return true to generate default result
	}

	// Create extended histograms and RooFit objects only if we have sufficient coverage
	for (size_t i = 0; i < template_hists.size(); ++i) {
		// Create extended histogram with wider range to ensure coverage
		double fit_min = charge.getMin();
		double fit_max = charge.getMax();
		int n_bins = 100;
		
		auto extended_hist = std::make_unique<TH1D>(
			Form("%s_extended_bin%d", templateElements[i].c_str(), energyBinIndex),
			template_hists[i]->GetTitle(),
			n_bins, fit_min, fit_max
		);
		
		// Fill extended histogram
		TH1D* original = template_hists[i].get();
		for (int j = 1; j <= original->GetNbinsX(); ++j) {
			double x = original->GetBinCenter(j);
			double content = original->GetBinContent(j);
			if (content > 0 && x >= fit_min && x <= fit_max) {
				extended_hist->Fill(x, content);
			}
		}
		
		// If no coverage in fit range, add minimal uniform distribution
		if (extended_hist->GetEntries() == 0) {
			for (int bin = 1; bin <= n_bins; ++bin) {
				extended_hist->SetBinContent(bin, 0.01);
			}
			std::cout << "        " << templateElements[i] << " filled with minimal uniform distribution" << std::endl;
		}
		
		std::cout << "        " << templateElements[i] << " extended histogram entries: " 
				  << extended_hist->GetEntries() << std::endl;
		
		auto rooData = std::make_unique<RooDataHist>(
				Form("%s_data_bin%d", templateElements[i].c_str(), energyBinIndex), 
				"", charge, extended_hist.get()
				);
		template_pdfs.push_back(std::make_unique<RooHistPdf>(
					Form("%s_pdf_bin%d", templateElements[i].c_str(), energyBinIndex), 
					"", charge, *rooData
					));
		template_data_hists.push_back(std::move(rooData));
		extended_template_hists.push_back(std::move(extended_hist));
	}

	// Create total PDF only if we have sufficient templates
	if (has_sufficient_templates) {
		RooArgList pdfList, fracList;
		for (auto& pdf : template_pdfs) {
			pdfList.add(*pdf);
		}
		for (auto& fvar : fractionVars) {
			fracList.add(*fvar);
		}
		if (lastFraction) {
			fracList.add(*lastFraction);
		}

		total_pdf = std::make_unique<RooAddPdf>(
				Form("total_pdf_bin%d", energyBinIndex), "Total PDF", pdfList, fracList, false);
	}

	return true;
}

EnergyBinResult EnergyBinFitProcessor::createDefaultResult() {
	EnergyBinResult result(templateElements.size());
	result.no_fit_needed = true;
	result.fit_status = 0;  // Mark as successful (no fit needed)
	result.chi2_ndf = 0;
	result.ndf = 1;
	result.n_signal = h_signal ? h_signal->GetEntries() : 0;
	
	// Set target element fraction to 1.0, others to 0.0
	std::string target_element;
	if (elementName == "Carbon") target_element = "Carbon";
	else if (elementName == "Nitrogen") target_element = "Nitrogen"; 
	else if (elementName == "Oxygen") target_element = "Oxygen";
	
	for (size_t i = 0; i < templateElements.size(); ++i) {
		if (templateElements[i] == target_element) {
			result.fit_fractions.fractions[i] = 1.0;
			result.fit_fractions.errors[i] = 0.1;  // Assign some uncertainty
			result.narrow_fractions.fractions[i] = 1.0;
			result.narrow_fractions.errors[i] = 0.1;
		} else {
			result.fit_fractions.fractions[i] = 0.0;
			result.fit_fractions.errors[i] = 0.0;
			result.narrow_fractions.fractions[i] = 0.0;
			result.narrow_fractions.errors[i] = 0.0;
		}
	}
	
	return result;
}

std::pair<RooPlot*, EnergyBinResult> EnergyBinFitProcessor::runFit() {
	// If insufficient template coverage, return default result
	if (!has_sufficient_templates) {
		auto result = createDefaultResult();
		
		// Create simple plot showing signal only
		RooPlot* frame = charge.frame(Title(Form("%s %s %s Energy Bin %d (No Fit)", 
				elementName.c_str(), detectorName.c_str(), 
				chargeTypeName.c_str(), energyBinIndex)));
		
		data_hist = std::make_unique<RooDataHist>(
				Form("data_bin%d", energyBinIndex), "data", charge, h_signal.get());
		
		data_hist->plotOn(frame, Name("data_hist"), MarkerStyle(20), MarkerSize(0.8), 
				MarkerColor(kBlack), LineColor(kBlack), XErrorSize(0), DrawOption("PZ"));
		
		return {frame, result};
	}

	// Proceed with normal fitting
	EnergyBinResult result(templateElements.size());

	data_hist = std::make_unique<RooDataHist>(
			Form("data_bin%d", energyBinIndex), "data", charge, h_signal.get());

	// Set initial values: target element = 0.99, others = small values
	std::string target_element;
	if (elementName == "Carbon") target_element = "Carbon";
	else if (elementName == "Nitrogen") target_element = "Nitrogen"; 
	else if (elementName == "Oxygen") target_element = "Oxygen";
	
	for (size_t i = 0; i < fractionVars.size(); ++i) {
		if (templateElements[i] == target_element) {
			fractionVars[i]->setVal(0.99);
		} else {
			fractionVars[i]->setVal(0.005);
		}
	}

	auto fitResult = std::unique_ptr<RooFitResult>(
			total_pdf->fitTo(*data_hist, Save(true), PrintLevel(-1), Strategy(2), Minimizer("Minuit2"))
			);

	result.fit_status = fitResult->status();

	if (result.fit_status != 0) {
		std::cout << "        Fit failed with status " << result.fit_status << ", using default result" << std::endl;
		return {nullptr, createDefaultResult()};
	}

	// Store fit results
	for (size_t i = 0; i < fractionVars.size(); ++i) {
		result.fit_fractions.fractions[i] = fractionVars[i]->getVal();
		result.fit_fractions.errors[i] = fractionVars[i]->getError();
	}
	if (lastFraction) {
		result.fit_fractions.fractions.back() = lastFraction->getVal();
		// Error propagation for last fraction
		double last_error = 0;
		for (auto& fvar : fractionVars) {
			last_error += pow(fvar->getError(), 2);
		}
		result.fit_fractions.errors.back() = sqrt(last_error);
	}

	// Create plot
	RooPlot* frame = charge.frame(Title(Form("%s %s %s Energy Bin %d", 
					elementName.c_str(), detectorName.c_str(), 
					chargeTypeName.c_str(), energyBinIndex)));
	data_hist->plotOn(frame, Name("data_hist"), MarkerStyle(20), MarkerSize(0.8), 
			MarkerColor(kBlack), LineColor(kBlack), XErrorSize(0), DrawOption("PZ"));
	total_pdf->plotOn(frame, Name("total_pdf"), LineColor(kRed), LineWidth(2));

	// Plot components
	std::vector<int> colors = {kBlue, kGreen+2, kMagenta, kOrange, kCyan};
	for (size_t i = 0; i < template_pdfs.size(); ++i) {
		total_pdf->plotOn(frame, Components(*template_pdfs[i]), 
				Name((templateElements[i] + "_pdf").c_str()), 
				LineColor(colors[i % colors.size()]), LineStyle(1), LineWidth(2));
	}

	// Calculate chi2
	double chi2 = calculateChi2(frame, "data_hist", "total_pdf", charge.getMin(), charge.getMax());
	int first_bin = h_signal->FindBin(charge.getMin());
	int last_bin = h_signal->FindBin(charge.getMax());
	result.ndf = (last_bin - first_bin + 1) - fitResult->floatParsFinal().getSize();
	result.chi2_ndf = (result.ndf > 0) ? chi2 / result.ndf : 0;

	// Calculate narrow range fractions
	result.narrow_fractions = calculateNarrowRangeFractions();
	auto range = fractionRanges.at(elementName);
	int bin_low = h_signal->FindBin(range.first);
	int bin_high = h_signal->FindBin(range.second);
	result.n_signal = h_signal->Integral(bin_low, bin_high);

	return {frame, result};
}

FractionResult EnergyBinFitProcessor::calculateNarrowRangeFractions() {
	auto range = fractionRanges.at(elementName);
	charge.setRange("narrow", range.first, range.second);

	FractionResult result(templateElements.size());

	if (!has_sufficient_templates || !total_pdf) {
		// Return default result
		std::string target_element;
		if (elementName == "Carbon") target_element = "Carbon";
		else if (elementName == "Nitrogen") target_element = "Nitrogen"; 
		else if (elementName == "Oxygen") target_element = "Oxygen";
		
		for (size_t i = 0; i < templateElements.size(); ++i) {
			if (templateElements[i] == target_element) {
				result.fractions[i] = 1.0;
				result.errors[i] = 0.1;
			}
		}
		return result;
	}

	auto& pdf_list = total_pdf->pdfList();
	std::vector<double> integrals;
	std::vector<double> fractionValues;

	// Get integrals for each component
	for (int i = 0; i < pdf_list.getSize(); ++i) {
		auto integral_obj = std::unique_ptr<RooAbsReal>(
				static_cast<RooAbsPdf*>(pdf_list.at(i))->createIntegral(charge, NormSet(charge), Range("narrow"))
				);
		integrals.push_back(integral_obj->getVal());
	}

	// Get fraction values
	for (auto& fvar : fractionVars) {
		fractionValues.push_back(fvar->getVal());
	}
	if (lastFraction) {
		fractionValues.push_back(lastFraction->getVal());
	}

	// Calculate narrow range signal events
	auto range_pair = fractionRanges.at(elementName);
	int bin_low = h_signal->FindBin(range_pair.first);
	int bin_high = h_signal->FindBin(range_pair.second);
	double n_signal_narrow = h_signal->Integral(bin_low, bin_high);

	// Calculate normalized fractions in narrow range
	double total_integral = 0;
	for (size_t i = 0; i < integrals.size(); ++i) {
		total_integral += fractionValues[i] * integrals[i];
	}

	if (total_integral > 1e-9) {
		for (size_t i = 0; i < result.fractions.size(); ++i) {
			result.fractions[i] = fractionValues[i] * integrals[i] / total_integral;
			// Simple error estimation
			if (n_signal_narrow > 0) {
				result.errors[i] = sqrt(result.fractions[i] * (1 - result.fractions[i]) / n_signal_narrow);
			}
		}
	}

	return result;
}

// Main analysis function
void runChargeTempFit() {
	gROOT->SetBatch(kTRUE);
	gStyle->SetOptStat(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	auto inputFile = std::unique_ptr<TFile>(TFile::Open(inputFileName.c_str()));
	if (!inputFile || inputFile->IsZombie()) {
		std::cerr << "Error: Cannot open input file: " << inputFileName << std::endl;
		return;
	}

	// Process each charge type
	for (const auto& chargeType : chargeTypes) {
		std::string outputPDF = outputDir + "chargeTempFit_" + chargeType + ".pdf";
		std::string outputROOT = outputDir + "chargeTempFit_" + chargeType + ".root";

		TCanvas* c0 = new TCanvas("c0", "pdf_canvas");
		c0->Print((outputPDF + "[").c_str());

		auto outputFile = std::unique_ptr<TFile>(TFile::Open(outputROOT.c_str(), "RECREATE"));

		for (const auto& element : elements) {
			for (const auto& detector : detectors) {
				std::cout << "\nProcessing " << chargeType << " " << element << " " << detector << std::endl;

				// Get signal histogram to determine number of energy bins
				std::string signalName = chargeType + "_ISS_BKG_H2_" + element + "_L1QSignal_" + detector;
				TH2F* h2d_signal = (TH2F*)inputFile->Get(signalName.c_str());
				if (!h2d_signal) {
					std::cerr << "Cannot find signal histogram: " << signalName << std::endl;
					continue;
				}

				// Rebin Y axis to see how many bins we have after rebinning
				TH2F* h2d_signal_rebinned = (TH2F*)h2d_signal->Clone("temp_rebinned");
				h2d_signal_rebinned->RebinY(2);
				int n_energy_bins = h2d_signal_rebinned->GetNbinsY();

				std::cout << "  Total energy bins after rebin(2): " << n_energy_bins << std::endl;

				// Storage for results vs energy
				std::vector<EnergyBinResult> energy_results;
				std::vector<double> energy_centers, energy_widths;

				// Loop over all energy bins
				for (int energy_bin = 1; energy_bin <= n_energy_bins; ++energy_bin) {
					double energy = h2d_signal_rebinned->GetYaxis()->GetBinCenter(energy_bin);
					double energy_width = h2d_signal_rebinned->GetYaxis()->GetBinWidth(energy_bin);

					std::cout << "    Processing energy bin " << energy_bin 
						<< " (E = " << energy << " GeV/n)" << std::endl;

					auto processor = std::make_unique<EnergyBinFitProcessor>(
							element, detector, chargeType, inputFile.get(), energy_bin, h2d_signal_rebinned);

					if (!processor->initialize()) {
						std::cerr << "      -> Initialization failed, skipping." << std::endl;
						continue;
					}

					auto [frame, result] = processor->runFit();
					if (!frame) {
						std::cerr << "      -> No plot generated, skipping." << std::endl;
						continue;
					}

					result.energy = energy;
					result.energy_width = energy_width;
					energy_results.push_back(result);
					energy_centers.push_back(energy);
					energy_widths.push_back(energy_width);

					// Create canvas with fit and pull plots
					TCanvas* canvas = new TCanvas(
							Form("c_%s_%s_%s_bin%d", chargeType.c_str(), element.c_str(), detector.c_str(), energy_bin), 
							"Fit", 800, 600);
					canvas->Divide(1, 2);

					// Main fit plot
					TPad* pad1 = (TPad*)canvas->cd(1);
					pad1->SetPad(0, 0.26, 1, 1);
					pad1->SetBottomMargin(0.02);
					pad1->SetLogy();

					double fit_min = elementConfig.at(element).fit_min;
					double fit_max = elementConfig.at(element).fit_max;
					frame->GetYaxis()->SetTitle("Events");
					frame->GetYaxis()->SetRangeUser(0.5, 2.5 * frame->GetMaximum());
					frame->GetXaxis()->SetRangeUser(fit_min, fit_max);
					frame->GetXaxis()->SetLabelSize(0);
					frame->Draw();

					auto legend = createLegend();
					legend->AddEntry("data_hist", "Signal", "pze");
					if (!result.no_fit_needed) {
						legend->AddEntry("total_pdf", "Fit", "l");
					}
					legend->Draw();

					auto info = createInfoBox(element, energy, result);
					info->Draw();

					// Pull plot (only if fit was performed)
					TPad* pad2 = (TPad*)canvas->cd(2);
					pad2->SetPad(0, 0, 1, 0.26);
					pad2->SetBottomMargin(0.3);
					pad2->SetGridy();
					pad2->SetTopMargin(0.02);

					if (!result.no_fit_needed) {
						auto pullGraph = std::make_unique<TGraphErrors>();
						calculatePull(frame, pullGraph.get(), fit_min, fit_max);
						setupPullPlot(pullGraph.get(), fit_min, fit_max);

						if (pullGraph->GetN() > 0) {
							pullGraph->Draw("AP");
							TLine zeroLine(fit_min, 0, fit_max, 0);
							zeroLine.SetLineStyle(2);
							zeroLine.SetLineColor(kGray + 1);
							zeroLine.Draw("SAME");
						}
					} else {
						// Empty pad for no fit case
						pad2->Clear();
						TPaveText* no_pull = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
						no_pull->AddText("No pull plot - no fit performed");
						no_pull->Draw();
					}

					canvas->Print(outputPDF.c_str());

					delete frame;
					delete canvas;

					std::cout << "      -> Completed successfully" << std::endl;
				}

				// Create summary plots of fractions vs energy
				if (!energy_results.empty()) {
					TCanvas* summary_canvas = new TCanvas(
							Form("summary_%s_%s_%s", chargeType.c_str(), element.c_str(), detector.c_str()), 
							"Fraction vs Energy", 1200, 800);

					// Create TGraphErrors for each template element
					std::vector<int> colors = {kBlue, kRed, kGreen+2, kMagenta, kOrange};
					auto legend_summary = createLegend();

					double y_min = 0, y_max = 1;
					bool first_plot = true;

					for (size_t i = 0; i < energy_results[0].fit_fractions.fractions.size(); ++i) {
						auto graph = std::make_unique<TGraphErrors>();

						for (size_t j = 0; j < energy_results.size(); ++j) {
							graph->SetPoint(j, energy_results[j].energy, energy_results[j].fit_fractions.fractions[i]);
							graph->SetPointError(j, energy_results[j].energy_width/2, energy_results[j].fit_fractions.errors[i]);
						}

						graph->SetMarkerStyle(20 + i);
						graph->SetMarkerColor(colors[i % colors.size()]);
						graph->SetLineColor(colors[i % colors.size()]);
						graph->SetMarkerSize(1.2);

						if (first_plot) {
							graph->SetTitle(Form("%s %s %s Fractions vs Energy;E_{k}/n [GeV/n];Fraction", 
										element.c_str(), detector.c_str(), chargeType.c_str()));
							graph->GetYaxis()->SetRangeUser(0, 1.1);
							graph->Draw("AP");
							first_plot = false;
						} else {
							graph->Draw("P SAME");
						}

						// Determine template element name
						std::string template_name;
						if (element == "Carbon") {
							std::vector<std::string> names = {"Beryllium", "Carbon", "Nitrogen"};
							template_name = names[i];
						} else if (element == "Nitrogen") {
							std::vector<std::string> names = {"Carbon", "Nitrogen", "Oxygen"};
							template_name = names[i];
						} else if (element == "Oxygen") {
							std::vector<std::string> names = {"Nitrogen", "Oxygen"};
							template_name = names[i];
						}

						legend_summary->AddEntry(graph.get(), template_name.c_str(), "p");

						// Save to ROOT file
						outputFile->cd();
						graph->Write(Form("fraction_%s_%s_%s_%s_vs_energy", 
									chargeType.c_str(), element.c_str(), detector.c_str(), template_name.c_str()));
					}

					legend_summary->Draw();
					summary_canvas->SetGrid();
					summary_canvas->Print(outputPDF.c_str());

					delete summary_canvas;
				}

				delete h2d_signal_rebinned;
			}
		}

		c0->Print((outputPDF + "]").c_str());
		delete c0;
		outputFile->Close();

		std::cout << "\nCompleted " << chargeType << " analysis. Results saved to:" << std::endl;
		std::cout << "  PDF: " << outputPDF << std::endl;
		std::cout << "  ROOT: " << outputROOT << std::endl;
	}

	std::cout << "\nAll analyses complete!" << std::endl;
}

// Entry point
void ChargeTempFit() {
	runChargeTempFit();
}