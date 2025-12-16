#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLine.h>
#include <TError.h>
#include <TROOT.h>
#include <TString.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <memory>
#include "../Tool.h"

using namespace std;
using namespace AMS_Iso;

// =================================================================================
// Data Structures
// =================================================================================

struct FitResult {
	TF1* fit = nullptr;
	double chi2 = 0;
	int ndf = 0;
};

struct Config {
	string main_analysis_file_path = "/eos/user/z/zixuan/Isotope/ChargeFit/ChargeFitParams_HeToOxy_iter0.root";
	string output_dir = "/eos/user/z/zixuan/Isotope/ChargeFit/comparison_plots/";
	//vector<string> elements = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
	vector<string> elements = {"Helium"};
	vector<string> detectors = {"TOF", "NaF", "AGL"};
	vector<string> chains = {"L1Inner", "UnbiasedL1Inner"};
	vector<string> templates = {"L1Template", "L2Template"};
	vector<pair<string, string>> paramList = {
		{"LG", "Width"}, {"LG", "MPV"}, {"LG", "Sigma"}, {"LG", "Chi2NDF"},
		{"EGE", "Peak"}, {"EGE", "SigmaL"}, {"EGE", "SigmaR"},
		{"EGE", "AlphaL"}, {"EGE", "AlphaR"}, {"EGE", "Chi2NDF"}
	};
	map<string, int> detColors = {{"TOF", kRed}, {"NaF", kBlue}, {"AGL", kGreen+2}};
	map<string, int> chainColors = {{"L1Inner", kBlack}, {"UnbiasedL1Inner", kRed}};
	map<string, int> templateColors = {{"L1Template", kBlack}, {"L2Template", kRed}};
	map<string, int> modelColors = {{"LG", kRed}, {"EGE", kBlue}};
	map<string, pair<double, double>> detRanges = {
		{"TOF", {0.3, 1.3}}, {"NaF", {0.71, 5.1}}, {"AGL", {2.8, 20.0}}
	};
};

Config gConfig;

// --- START OF MODIFICATION 1: Replace global vector with a global file pointer ---
unique_ptr<TFile> gSplineOutFile = nullptr;
// --- END OF MODIFICATION 1 ---


// =================================================================================
// Helper Functions
// =================================================================================
TH1* getHist(TFile* file, const string& name) {
	if (!file || file->IsZombie()) return nullptr;
	auto* h = dynamic_cast<TH1*>(file->Get(name.c_str()));
	if (!h) return nullptr;
	h->SetDirectory(nullptr); h->SetStats(0);
	return h;
}

bool isValidHist(TH1* h) {
	if (!h) return false;
	for (int b = 1; b <= h->GetNbinsX(); b++) if (h->GetBinContent(b) > 0) return true;
	return false;
}

pair<double, double> getYRange(const vector<TH1*>& hists) {
	double minY = 1e30, maxY = -1e30; bool hasData = false;
	for (auto h : hists) {
		if (!isValidHist(h)) continue; hasData = true;
		for (int b = 1; b <= h->GetNbinsX(); b++) {
			double v = h->GetBinContent(b), e = h->GetBinError(b);
			if (v > 0 && isfinite(v)) { minY = min(minY, v - e); maxY = max(maxY, v + e); }
		}
	}
	if (!hasData || minY >= maxY) return {0.1, 1.0};
	double margin = (maxY - minY) * 0.4;
	return {max(1e-9, minY - margin), maxY + margin};
}

vector<double> generateSplineNodes(int segments, const string& detector, const vector<double>& valid_bin_centers, double fitStart, double fitEnd) {
	vector<double> xpoints;
	xpoints.push_back(fitStart);

	if (false) {
		if (valid_bin_centers.size() > 1) {
			xpoints.push_back(valid_bin_centers[1] + 0.01);
		}
		if (segments > 1 && valid_bin_centers.size() > 2) {
			int num_remaining_points = valid_bin_centers.size() - 2;
			for (int i = 1; i < segments; ++i) {
				int point_idx = 2 + (num_remaining_points * i) / segments;
				if (point_idx < valid_bin_centers.size()) {
					xpoints.push_back(valid_bin_centers[point_idx]);
				}
			}
		}
	} else {
		for (int i = 1; i < segments; ++i) {
			int point_idx = (valid_bin_centers.size() * i) / segments;
			if (point_idx < valid_bin_centers.size()) {
				xpoints.push_back(valid_bin_centers[point_idx]);
			}
		}
	}

	xpoints.push_back(fitEnd);
	sort(xpoints.begin(), xpoints.end());
	xpoints.erase(unique(xpoints.begin(), xpoints.end()), xpoints.end());
	return xpoints;
}

// =================================================================================
// Spline Fitting - Returns a FitResult struct
// =================================================================================
FitResult performSplineFit(TH1* hist, const string& detector, const string& baseName) {
	if (!isValidHist(hist)) {
		return {};
	}

	double fitStart = gConfig.detRanges[detector].first;
	double fitEnd = gConfig.detRanges[detector].second;

	vector<double> valid_bin_centers;
	for (int b = 1; b <= hist->GetNbinsX(); ++b) {
		double bin_center = hist->GetBinCenter(b);
		if (hist->GetBinContent(b) > 0 && bin_center >= fitStart && bin_center <= fitEnd) {
			valid_bin_centers.push_back(bin_center);
		}
	}

	if (valid_bin_centers.size() < 2) {
		cout << "[WARN] Not enough valid points for " << baseName << endl;
		return {};
	}

	cout << "\n[INFO] Starting smart spline fit for: " << baseName << " (" << valid_bin_centers.size() << " points)" << endl;

	const map<int, double> chi2ndf_thresholds = {
		{1, 1.8}, {2, 2.5}, {3, 3.0}, {4, 4.0}, {5, 5.0}, {6, 4.0}, {7, 4.0}
	};
	const int max_segments = (detector == "TOF") ? 5 : 7;

	map<int, pair<double, int>> fit_quality;
	int best_seg_idx = -1;

	for (int segments = 1; segments <= max_segments; ++segments) {
		vector<double> xpoints = generateSplineNodes(segments, detector, valid_bin_centers, fitStart, fitEnd);
		if (xpoints.size() < 2) continue;

		TF1* temp_fit = nullptr;
		try {
			temp_fit = SplineFit(hist, xpoints.data(), xpoints.size(), 0x38, "b1e1", "temp_fit", 0.25, 22);
		} catch (...) { continue; }

		if (!temp_fit || temp_fit->GetNDF() < 1) {
			if(temp_fit) delete temp_fit;
			continue;
		}

		double current_chi2 = temp_fit->GetChisquare();
		int current_ndf = temp_fit->GetNDF();
		double current_chi2ndf = current_ndf > 0 ? current_chi2 / current_ndf : 1e9;

		fit_quality[segments] = {current_chi2, current_ndf};
		delete temp_fit;

		if (current_chi2ndf < chi2ndf_thresholds.at(segments) && current_ndf > 1) {
			if (segments == 1) {
				best_seg_idx = 1;
			} else {
				double prev_chi2ndf = fit_quality[segments - 1].first / fit_quality[segments - 1].second;
				if (prev_chi2ndf <= current_chi2ndf + 0.5) {
					best_seg_idx = segments - 1;
				} else {
					best_seg_idx = segments;
				}
			}
			cout << "  [FIT] Found good fit candidate. Best segments so far: " << best_seg_idx << ". Stopping search." << endl;
			break;
		}
	}

	if (best_seg_idx == -1 && !fit_quality.empty()) {
		cout << "  [FIT] No fit met threshold. Falling back to best chi2/ndf." << endl;
		double min_chi2ndf = 1e9;
		for (auto const& [seg, quality] : fit_quality) {
			if (quality.second > 0) {
				double chi2ndf = quality.first / quality.second;
				if (chi2ndf < min_chi2ndf) {
					min_chi2ndf = chi2ndf;
					best_seg_idx = seg;
				}
			}
		}
	}

	if (best_seg_idx != -1) {
		cout << "  [FIT] Finalizing with loop variable 'segments' = " << best_seg_idx << endl;
		vector<double> xpoints = generateSplineNodes(best_seg_idx, detector, valid_bin_centers, fitStart, fitEnd);
		string final_fit_name = baseName + "_spline";
		TF1* final_fit = SplineFit(hist, xpoints.data(), xpoints.size(), 0x38, "b2e2", final_fit_name.c_str(), 0.25, 22);

		if (!final_fit) return {};

		double final_chi2 = fit_quality[best_seg_idx].first;
		int final_ndf = fit_quality[best_seg_idx].second;

		cout << "[SUCCESS] Created final fit for " << final_fit->GetName() << " with chi2/ndf = " << (final_ndf > 0 ? final_chi2 / final_ndf : 0.0) << endl;

		// --- START OF MODIFICATION 2: Save the fit immediately ---
		if (gSplineOutFile && gSplineOutFile->IsOpen()) {
			gSplineOutFile->cd(); // Switch to the output file's directory context
			final_fit->Write();
		} else {
			cout << "[ERROR] Global output file is not open. Cannot save fit." << endl;
		}
		// --- END OF MODIFICATION 2 ---

		return {final_fit, final_chi2, final_ndf};
	}

	cout << "[FAIL] Could not produce any valid spline fit for " << baseName << endl;
	return {};
}

// =================================================================================
// Plotting and Main Functions
// =================================================================================

void drawPage(TCanvas* c, const vector<TH1*>& hists, const vector<int>& colors,
		const vector<string>& labels, const string& title, const string& yTitle,
		double xMin, double xMax, const string& pdfName, bool drawBoundaries = false,
		const vector<FitResult>& fitResults = {}) {
	c->Clear(); c->SetLogx();
	cout<<"begin drawing "<<title<<endl;
	if (!any_of(hists.begin(), hists.end(), isValidHist)) return;

	auto yRange = getYRange(hists);
	TH1F* frame = c->DrawFrame(xMin, yRange.first, xMax, yRange.second);
	frame->SetTitle(""); frame->GetXaxis()->SetTitle("E_{k}/n [GeV/n]"); frame->GetYaxis()->SetTitle(yTitle.c_str());

	for (size_t i = 0; i < hists.size(); i++) {
		if (isValidHist(hists[i])) {
			hists[i]->SetLineColor(colors[i]); hists[i]->SetMarkerColor(colors[i]);
			hists[i]->SetLineWidth(2); hists[i]->SetMarkerStyle(20);
			hists[i]->Draw("E SAME");
		}
	}

	TLatex chi2_latex;
	chi2_latex.SetNDC();
	chi2_latex.SetTextSize(0.045);

	for (size_t i = 0; i < fitResults.size(); i++) {
		const auto& result = fitResults[i];
		if (result.fit) {
			result.fit->SetLineColor(colors[i]); 
			result.fit->SetLineWidth(2); 
			result.fit->Draw("SAME");

			if (result.ndf > 0) {
				chi2_latex.SetTextColor(colors[i]);
				TString text_content = TString::Format("#chi^{2}/ndf = %.2f/%d = %.2f", result.chi2, result.ndf, result.chi2 / result.ndf);
				chi2_latex.DrawLatex(0.18, 0.83 - i * 0.06, text_content.Data());
			}
		}
	}

	if (drawBoundaries) {
		TLine l(1.1, yRange.first, 1.1, yRange.second);
		l.SetLineColor(kGray+2); l.SetLineStyle(2); l.SetLineWidth(2); l.DrawClone();
		l.SetX1(5.0); l.SetX2(5.0); l.DrawClone();
	}

	TLatex main_latex; main_latex.SetNDC(); main_latex.SetTextSize(0.045);
	main_latex.DrawLatex(0.12, 0.93, title.c_str());

	double y_pos = 0.84;
	for (size_t i = 0; i < labels.size(); i++) {
		main_latex.SetTextColor(colors[i]);
		main_latex.DrawLatex(0.79, y_pos, labels[i].c_str());
		y_pos -= 0.06;
	}

	c->Update(); 
	c->Print(pdfName.c_str());
}

void plotTemplateComparison(TFile* fin, TCanvas* c, const string& pdfName) {
	cout << "\n[PROCESS] Starting Template Comparison plots..." << endl;
	for (const auto& element : gConfig.elements) for (const auto& chain : gConfig.chains) for (const auto& det : gConfig.detectors)
		for (const auto& param : gConfig.paramList) {
			vector<TH1*> hists; 
			vector<FitResult> fit_results; 
			vector<string> labels; 
			vector<int> colors;

			for (const auto& temp : gConfig.templates) {
				string name = chain+"_"+element+"_"+det+"_"+temp+"_"+param.first+"_"+param.second;
				TH1* h = getHist(fin, name);
				hists.push_back(h);

				if (param.second == "AlphaL") {
					// --- 修改开始：计算探测器范围内的平均值 ---

					double sumAlpha = 0.0;
					int countBins = 0;

					// 获取当前探测器的有效能量范围
					double rangeMin = gConfig.detRanges[det].first;
					double rangeMax = gConfig.detRanges[det].second;

					for (int b = 1; b <= h->GetNbinsX(); ++b) {
						double center = h->GetBinCenter(b);
						double content = h->GetBinContent(b);

						// 仅累加：1. 内容大于0 (有效值) 且 2. 在探测器定义的能量范围内
						if (content > 0 && center >= rangeMin && center <= rangeMax) {
							sumAlpha += content;
							countBins++;
						}
					}

					// 防止除以0，如果没有有效bin，给一个默认值(例如 h->GetMaximum() 或 0)
					double alpha_avg_val = (countBins > 0) ? (sumAlpha / countBins) : h->GetMaximum();

					// --- 修改结束 ---

					string fit_name = name + "_spline";
					// 使用计算出的平均值创建常数函数
					TF1* const_fit = new TF1(fit_name.c_str(), Form("%f", alpha_avg_val), 0.25, 22);

					if (gSplineOutFile && gSplineOutFile->IsOpen()) {
						gSplineOutFile->cd();
						const_fit->Write();
						cout << "  [SUCCESS] Wrote constant TF1 '" << fit_name << "' (Avg=" << alpha_avg_val << ") to file." << endl;
					}
					fit_results.push_back({const_fit, 0, 0});

				} else if (param.second != "Chi2NDF") {
					fit_results.push_back(performSplineFit(h, det, name));
				} else {
					fit_results.push_back({});
				}

				labels.push_back(temp == "L1Template" ? "L1Temp" : "L2Temp");
				colors.push_back(gConfig.templateColors[temp]);
			}

			string title = element+" "+chain+" "+det+" "+param.first+"_"+param.second;
			drawPage(c, hists, colors, labels, title, param.second, 
					gConfig.detRanges[det].first, gConfig.detRanges[det].second, pdfName, false, fit_results);

			// --- START OF MODIFICATION 3: Clean up memory after drawing ---
			// The TF1 objects have been saved and drawn, now we delete them from memory.
			for (auto& result : fit_results) {
				if (result.fit) {
					delete result.fit;
				}
			}
			for(auto h : hists) {
				if (h) delete h;
			}
			// --- END OF MODIFICATION 3 ---
		}
}

void plotDetectorComparison(TFile* fin, TCanvas* c, const string& pdfName) {
	cout << "\n[PROCESS] Starting Detector Comparison plots..." << endl;
	for (const auto& element : gConfig.elements) for (const auto& chain : gConfig.chains) for (const auto& temp : gConfig.templates)
		for (const auto& param : gConfig.paramList) {
			vector<TH1*> hists;
			for (const auto& det : gConfig.detectors) {
				string name = chain+"_"+element+"_"+det+"_"+temp+"_"+param.first+"_"+param.second;
				hists.push_back(getHist(fin, name));
			}
			vector<int> colors;
			for(const auto& det : gConfig.detectors) colors.push_back(gConfig.detColors[det]);
			string title = element+" "+chain+" "+temp+" "+param.first+"_"+param.second;
			drawPage(c, hists, colors, gConfig.detectors, title, param.second, 0.4, 20.0, pdfName, true);
			for(auto h : hists) if (h) delete h;
		}
}

void plotChainComparison(TFile* fin, TCanvas* c, const string& pdfName) {
	cout << "\n[PROCESS] Starting Chain Comparison plots..." << endl;
	for (const auto& element : gConfig.elements) for (const auto& det : gConfig.detectors) for (const auto& temp : gConfig.templates)
		for (const auto& param : gConfig.paramList) {
			vector<TH1*> hists;
			vector<int> colors;
			for (const auto& chain : gConfig.chains) {
				string name = chain+"_"+element+"_"+det+"_"+temp+"_"+param.first+"_"+param.second;
				hists.push_back(getHist(fin, name));
				colors.push_back(gConfig.chainColors[chain]);
			}
			string title = element+" "+det+" "+temp+" "+param.first+"_"+param.second;
			auto range = gConfig.detRanges[det];
			drawPage(c, hists, colors, gConfig.chains, title, param.second, range.first, range.second, pdfName);
			for(auto h : hists) if (h) delete h;
		}
}

void plotModelChi2Comparison(TFile* fin, TCanvas* c, const string& pdfName) {
	cout << "\n[PROCESS] Starting Model Chi2 Comparison plots..." << endl;
	vector<string> models = {"LG", "EGE"};
	vector<string> modelLabels = {"Landau-Gauss", "ExpGausExp"};
	for (const auto& element : gConfig.elements) for (const auto& chain : gConfig.chains) for (const auto& det : gConfig.detectors)
		for (const auto& temp : gConfig.templates) {
			vector<TH1*> hists;
			vector<int> colors;
			for (const auto& model : models) {
				string name = chain+"_"+element+"_"+det+"_"+temp+"_"+model+"_Chi2NDF";
				hists.push_back(getHist(fin, name));
				colors.push_back(gConfig.modelColors[model]);
			}
			string title = element+" "+chain+" "+det+" "+temp+" Chi2/NDF";
			auto range = gConfig.detRanges[det];
			drawPage(c, hists, colors, modelLabels, title, "#chi^{2}/ndf", range.first, range.second, pdfName);
			for(auto h : hists) if(h) delete h;
		}
}

// --- START OF MODIFICATION 4: The saveSplineFits function is no longer needed ---
// void saveSplineFits(const string& outputFile) { ... } // This function is removed.
// --- END OF MODIFICATION 4 ---

void compareFitResults() {
	cout << "[INFO] Starting compareFitResults..." << endl;
	gStyle->SetOptStat(0); gErrorIgnoreLevel = kWarning;
	unique_ptr<TFile> fin_main(TFile::Open(gConfig.main_analysis_file_path.c_str()));
	TCanvas* c = new TCanvas("c", "c", 800, 400); c->SetGrid();

	// --- START OF MODIFICATION 5: Open the output file at the beginning ---
	string splineOutFile = gConfig.output_dir + "allFitHistSplineSmooth_iter0.root";
	gSplineOutFile = make_unique<TFile>(splineOutFile.c_str(), "RECREATE");
	if (!gSplineOutFile || gSplineOutFile->IsZombie()) {
		cout << "[ERROR] Cannot create output spline file: " << splineOutFile << endl;
		return;
	}
	cout << "[INFO] Opened output file for splines: " << splineOutFile << endl;

	string pdfName_template = gConfig.output_dir + "FitParam_TemplateCompare_iter0.pdf";
	c->Print((pdfName_template + "[").c_str());
	plotTemplateComparison(fin_main.get(), c, pdfName_template);
	c->Print((pdfName_template + "]").c_str());

	string pdfName_detector = gConfig.output_dir + "FitParam_DetectorCompare_iter0.pdf";
	c->Print((pdfName_detector + "[").c_str());
	plotDetectorComparison(fin_main.get(), c, pdfName_detector);
	c->Print((pdfName_detector + "]").c_str());

	string pdfName_chain = gConfig.output_dir + "FitParam_ChainCompare_iter0.pdf";
	c->Print((pdfName_chain + "[").c_str());
	plotChainComparison(fin_main.get(), c, pdfName_chain);
	c->Print((pdfName_chain + "]").c_str());

	string pdfName_chi2 = gConfig.output_dir + "FitParam_Chi2ModelCompare_iter0.pdf";
	c->Print((pdfName_chi2 + "[").c_str());
	plotModelChi2Comparison(fin_main.get(), c, pdfName_chi2);
	c->Print((pdfName_chi2 + "]").c_str());

	// --- START OF MODIFICATION 6: Close the file at the end ---
	gSplineOutFile->Close();
	cout << "\n[INFO] All spline fits have been saved to " << splineOutFile << endl;
	// --- END OF MODIFICATION 6 ---

	cout << "\n[INFO] All done!" << endl;
	delete c;
}
