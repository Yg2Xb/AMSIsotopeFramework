#include <TSystem.h>
#include <TAxis.h>
#include <iostream>
#include <memory>
#include <algorithm>
#include "HistTempFit_hist.h"
#include "./helper_func.cpp"
#include "../Tool.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::unique_ptr;

using namespace AMS_Iso;

// Main fitting function
// ===== MODIFICATION START: Added fitFragMass and sourceName parameters =====
void HistTempFit(const string& isotype, int UseMass, bool useUnbiasedChain = true, bool usePureTemplates = false, int rebinX, int ProNbin, bool fitFragMass = false, const string& sourceName = "") {
// ===== MODIFICATION END =====
    
    // 1. ============================ Configuration and Initialization ============================
    cout << "====================================================================\n";
    // ===== MODIFICATION START: Adjust log message to reflect fit type =====
    if (fitFragMass) {
        cout << "Starting FRAGMENT Mass Template Fit for: " << isotype << " from source " << sourceName << " (UseMass: " << UseMass << ")\n";
    } else {
        cout << "Starting SIGNAL Mass Template Fit for: " << isotype << " (UseMass: " << UseMass << ")\n";
    }
    // ===== MODIFICATION END =====
    cout << "Chain: " << (useUnbiasedChain ? "UnbiasedL1Inner" : "L1Inner") << "\n";
    cout << "Templates: " << (usePureTemplates ? "Pure (H2)" : "Mixed (H1)") << "\n";
    cout << "====================================================================\n";

    // Get isotope configuration
    if (IsoFitConstants::configs.find(isotype) == IsoFitConstants::configs.end()) {
        cout << "[ERROR] Isotope configuration for '" << isotype << "' not found." << endl;
        return;
    }
    const auto& config = IsoFitConstants::configs.at(isotype);
    const string chainName = useUnbiasedChain ? "UnbiasedL1Inner" : "L1Inner";
    const string templateType = usePureTemplates ? "H2" : "H1";
    const string suffix = usePureTemplates ? "" : "MC_";

    // Construct input file paths (UNCHANGED as requested)
    string dataFilePath = "/eos/user/z/zixuan/Isotope/Add/" + config.name + Form("_frag%d_NoBkg_NoTune_full.root",config.charge);
    vector<unique_ptr<TFile>> f_mc_vec;
    for (int mass : config.masses) {
        string mcFilePath = Form("/eos/user/z/zixuan/Isotope/Add/%s%d_rew_frag%d_NoBkg_full.root", config.name.c_str(), mass, config.charge);
        f_mc_vec.emplace_back(TFile::Open(mcFilePath.c_str()));
        if (!f_mc_vec.back() || f_mc_vec.back()->IsZombie()) {
            cout << "[ERROR] Failed to open MC file: " << mcFilePath << endl;
            return;
        }
    }

    unique_ptr<TFile> f_data(TFile::Open(dataFilePath.c_str()));
    if (!f_data || f_data->IsZombie()) {
        cout << "[ERROR] Failed to open data file: " << dataFilePath << endl;
        return;
    }
    cout << "[INFO] Data file loaded: " << f_data->GetName() << endl;

    // Construct output file paths
    string outputDir = "/eos/user/z/zixuan/Isotope/MassTempFit";
    gSystem->mkdir(outputDir.c_str(), true);

    // --- MODIFICATION TO INCLUDE FRAG/SOURCE INFO ---
    // Create an extra suffix only if it's a fragment fit
    string fitTypeSuffix = fitFragMass ? Form("_FragFrom%s", sourceName.c_str()) : "";
    // Append this extra suffix to the original file suffix
    string fileSuffix = Form("_%s_%s_UseMass%d%s", chainName.c_str(), templateType.c_str(), UseMass, fitTypeSuffix.c_str());
    // --- END OF MODIFICATION ---

    string outputPdfPath = outputDir + "/wide_MassTF_" + config.name + fileSuffix + "_NoBkg.pdf";
    string outputRootPath = outputDir + "/wide_MassTF_" + config.name + fileSuffix + "_NoBkg.root";

    unique_ptr<TFile> output_file(TFile::Open(outputRootPath.c_str(), "RECREATE"));
    cout << "[INFO] Output ROOT file: " << outputRootPath << endl;

    // Create canvas for plotting
    TCanvas* canvas = new TCanvas("canvas", "Template Fit", 800, 600);
    canvas->Divide(1, 2);
    TPad* pad1 = (TPad*)canvas->cd(1);
    pad1->SetPad(0, 0.25, 1, 1);
    TPad* pad2 = (TPad*)canvas->cd(2);
    pad2->SetPad(0, 0, 1, 0.27);
    pad2->SetBottomMargin(0.3);
    pad2->SetGridy();
    canvas->Print((outputPdfPath + "[").c_str());

    // 2. ============================ Detector and Energy Loop ============================
    for (int idet = 0; idet < 3; ++idet) {
        cout << "\n--- Processing Detector: " << DetName[idet] << " ---\n";

        // ===== MODIFICATION START: Dynamically select the data histogram to be fitted =====
        string dataHistName;
        if (fitFragMass) {
            // New logic for fitting fragment mass from BKG histogram
            dataHistName = Form("%s_BKG_H2b_%s_%s", 
                                chainName.c_str(), 
                                sourceName.c_str(), 
                                DetName[idet], 
                                config.charge, 
                                UseMass);
        } else {
            // Original logic for fitting signal mass from ID histogram
            dataHistName = Form("%s_ID_H2_%s_Mass%dBin", chainName.c_str(), DetName[idet], UseMass);
        }
        // ===== MODIFICATION END =====
        
        TH2F* data_hist_2d = (TH2F*)f_data->Get(dataHistName.c_str());
        if (!data_hist_2d) {
            cout << "[WARN] Data histogram not found: " << dataHistName << ". Skipping detector " << DetName[idet] << "." << endl;
            continue;
        }
        
        // Extract energy bins from the Y-axis
        TAxis* y_axis = data_hist_2d->GetYaxis();
        vector<double> ek_bins;
        for(int i = 1; i <= y_axis->GetNbins() + 1; ++i) {
            ek_bins.push_back(y_axis->GetBinLowEdge(i));
        }

        // Create result histograms for this detector
        RooRealVar inv_mass("inv_mass", "1/mass", config.fitRangeLow[idet], config.fitRangeUp[idet]);
        vector<unique_ptr<TH1F>> h_best_fractions;
        for (int i = 0; i < config.nFitParams; ++i) {
            h_best_fractions.emplace_back(new TH1F(
                Form("h_best_%s%d_frac_%s", config.name.c_str(), config.masses[i], DetName[idet]),
                Form("%s Best %s%d Fraction;E_{k}/n [GeV/n];Fraction", DetName[idet], config.name.c_str(), config.masses[i]),
                ek_bins.size() - 1, ek_bins.data()));
        }
        auto h_best_chi2 = std::make_unique<TH1F>(Form("h_best_chi2_%s", DetName[idet]),
            Form("%s Best #chi^{2}/NDF;E_{k}/n [GeV/n];#chi^{2}/NDF", DetName[idet]), ek_bins.size() - 1, ek_bins.data());
        auto h_best_entries = std::make_unique<TH1F>(Form("h_best_entries_%s", DetName[idet]),
            Form("%s Entries;E_{k}/n [GeV/n];Entries", DetName[idet]), ek_bins.size() - 1, ek_bins.data());

        // Loop over energy bins
        for (int ibin = 1; ibin < ek_bins.size(); ibin = ibin+ProNbin) {
            double ek_center = 0.5 * (ek_bins[ibin-1] + ek_bins[ibin+ProNbin-1]);

            // Check if bin center is within the detector's valid energy range
            if (ek_center < DetRanges[idet][0] || ek_center > DetRanges[idet][1]) {
                continue;
            }

            // Get MC templates for this energy bin (UNCHANGED as requested)
            vector<TH1D*> mc_hists;
            bool templates_ok = true;
            for (size_t i = 0; i < config.masses.size(); ++i) {
                string mcHistName = Form("%s_%sID_%s_%s_Mass%dBin", chainName.c_str(), suffix.c_str(), templateType.c_str(), DetName[idet], config.masses[i]);
                TH2F* mc_hist_2d = (TH2F*)f_mc_vec[i]->Get(mcHistName.c_str());
                if (!mc_hist_2d) {
                    cout << "[WARN] MC template not found: " << mcHistName << endl;
                    templates_ok = false;
                    break;
                }
                TH1D* mc_proj = mc_hist_2d->ProjectionX(Form("mc_proj_%d_%d", ibin, i), ibin, ibin+ProNbin-1);
                cout<<ProNbin<<endl;
                mc_proj->Rebin(rebinX);
                mc_proj->Smooth(1);
                mc_hists.push_back(mc_proj);
            }
            if (!templates_ok) {
                for(auto h : mc_hists) delete h;
                continue;
            }

            // Project data histogram for this energy bin
            TH1D* data_hist = data_hist_2d->ProjectionX(Form("data_proj_%d", ibin), ibin, ibin+ProNbin-1);
            data_hist->Rebin(rebinX);
            data_hist->Sumw2();

            // Check statistics
            if (data_hist->GetEntries() == 0 ||
                std::any_of(mc_hists.begin(), mc_hists.end(), [](TH1D* h) { return h->GetEntries() == 0; })) {
                cout<<"[WARN] Low statistics in bin " << ibin<<" max:"<< data_hist->GetMaximum() << ". Skipping." << endl;
                delete data_hist;
                for (auto* h : mc_hists) delete h;
                continue;
            }
            
            // 3. ============================ RooFit Fitting ============================
            RooDataHist data("data", "Data", RooArgList(inv_mass), data_hist);
            vector<unique_ptr<RooDataHist>> templates;
            vector<unique_ptr<RooHistPdf>> pdfs;
            vector<RooRealVar*> fractions;

            for (size_t i = 0; i < mc_hists.size(); ++i) {
                templates.emplace_back(new RooDataHist(Form("template%d", config.masses[i]), "Template", RooArgList(inv_mass), mc_hists[i]));
                pdfs.emplace_back(new RooHistPdf(Form("pdf%d", config.masses[i]), "PDF", RooArgSet(inv_mass), *templates.back()));
                if (i < config.nFitParams) {
                    fractions.push_back(new RooRealVar(Form("frac_%d", config.masses[i]), "Fraction", config.initFractions[i], 0., 1.));
                }
            }

            string formula = "1.0";
            RooArgList fracList;
            for (auto* frac : fractions) {
                formula += " - " + string(frac->GetName());
                fracList.add(*frac);
            }
            auto last_frac = std::make_unique<RooFormulaVar>(Form("frac_%d", config.masses.back()), "Last Fraction", formula.c_str(), fracList);

            RooArgList pdf_list, frac_list_model;
			for (size_t i = 0; i < pdfs.size(); ++i) {
				pdf_list.add(*pdfs[i]);
				if (i < fractions.size()) {
					frac_list_model.add(*fractions[i]);
				} else {
					frac_list_model.add(*last_frac);
				}
			}
            RooAddPdf model("model", "Combined Model", pdf_list, frac_list_model);

            // Perform the fit
            cout << "Fitting bin " << ibin << " (Ek: " << ek_bins[ibin-1] << "-" << ek_bins[ibin+ProNbin-1] << " GeV/n)..." << endl;
            unique_ptr<RooFitResult> fit_res(model.fitTo(data, RooFit::Save(), RooFit::SumW2Error(kTRUE), RooFit::PrintLevel(-1)));
            
            // 4. ============================ Plotting and Saving Results ============================
            RooPlot* frame = inv_mass.frame();
            data.plotOn(frame, RooFit::XErrorSize(0), RooFit::Name("data"));
            model.plotOn(frame, RooFit::LineColor(kRed), RooFit::LineWidth(3), RooFit::Name("model"));

            // Calculate Chi2/NDF
            int nBinsInRange = data_hist->FindBin(config.fitRangeUp[idet]) - data_hist->FindBin(config.fitRangeLow[idet]);
            int nParams = model.getParameters(data)->getSize();
            int ndf = nBinsInRange - nParams;
            double chi2 = calculateChi2(frame, "data", "model", config.fitRangeLow[idet], config.fitRangeUp[idet]);
            double chi2_ndf = (ndf > 0 && fit_res->status() == 0) ? chi2 / ndf : 5000;
            cout << "  Fit Result: Chi2/NDF = " << chi2 << "/" << ndf << " = " << chi2_ndf << ", Status = " << fit_res->status() << endl;

            // Plot components
            const vector<int> lineColors = {kBlue, kMagenta, kOrange + 1, kViolet};
            for (size_t i = 0; i < pdfs.size(); ++i) {
                model.plotOn(frame, RooFit::Components(*pdfs[i]), RooFit::LineStyle(kSolid), RooFit::LineColor(lineColors[i % lineColors.size()]), RooFit::LineWidth(3), RooFit::Name(Form("comp%d", config.masses[i])));
            }

            // Draw main plot
            canvas->cd(1);
            pad1->SetLogy(0);
            frame->SetXTitle(Form("1/%s Mass", DetName[idet]));
            frame->SetYTitle("Events");
            frame->GetYaxis()->SetTitleOffset(1.2);
            frame->SetTitle(Form("%s: E_{k}/n in [%.2f, %.2f] GeV/n", DetName[idet], ek_bins[ibin-1], ek_bins[ibin+ProNbin-1]));
            frame->Draw();

            // Legend
            TLegend* legend = new TLegend(0.16, 0.61, 0.43, 0.85);
            // ===== MODIFICATION START: Adjust legend entry for data source =====
            if (fitFragMass) {
                legend->AddEntry("data", Form("ISS Data (Frag. from %s)", sourceName.c_str()), "ep");
            } else {
                legend->AddEntry("data", "ISS Data", "ep");
            }
            // ===== MODIFICATION END =====
            legend->AddEntry("model", "Total Fit", "l");
            for (size_t i = 0; i < config.masses.size(); ++i) {
                legend->AddEntry(Form("comp%d", config.masses[i]), Form("%s%d Template", config.name.c_str(), config.masses[i]), "l");
            }
            setLegend(legend);
            legend->Draw("same");

            // Info box
            TPaveText* pt = new TPaveText(0.68, 0.55, 0.92, 0.88, "NDC");
            double n_entries = data_hist->Integral(data_hist->FindBin(config.fitRangeLow[idet]), data_hist->FindBin(config.fitRangeUp[idet]));
            pt->AddText(Form("Entries: %.0f", n_entries));
            for (size_t i = 0; i < fractions.size(); ++i) {
                pt->AddText(Form("%s%d frac: %.3f #pm %.3f", config.name.c_str(), config.masses[i], fractions[i]->getVal(), fractions[i]->getError()));
            }
            if (fit_res->status() == 0) {
                pt->AddText(Form("#chi^{2}/NDF: %.2f", chi2_ndf));
            } else {
                pt->AddText(Form("Fit Status: %d", fit_res->status()));
            }
            setPaveText(pt);
            pt->Draw("same");

            // Pull Plot
            canvas->cd(2);
            TGraphErrors* pullGraph = new TGraphErrors();
            pullGraph->SetTitle(Form(";Pull;1/%s Mass", DetName[idet]));
            setupPullPlot(pullGraph, Form(";Pull;1/%s Mass", DetName[idet]), config.fitRangeLow[idet], config.fitRangeUp[idet]);
            calculatePull(frame, pullGraph, "data", "model", config.fitRangeLow[idet], config.fitRangeUp[idet]);
            pullGraph->GetXaxis()->SetRangeUser(config.fitRangeLow[idet], config.fitRangeUp[idet]);
            pullGraph->GetYaxis()->SetRangeUser(-6, 6);
            pullGraph->Draw("AP");

            canvas->Print(outputPdfPath.c_str());

            // Store results if fit was successful
            if (fit_res->status() == 0) {
                for (size_t i = 0; i < fractions.size(); ++i) {
                    h_best_fractions[i]->SetBinContent(ibin, fractions[i]->getVal());
                    h_best_fractions[i]->SetBinError(ibin, fractions[i]->getError());
                }
                h_best_chi2->SetBinContent(ibin, chi2_ndf);
                h_best_entries->SetBinContent(ibin, n_entries);
            }

            // Cleanup for this bin
            delete legend;
            delete pt;
            delete pullGraph;
            delete frame;
            for (auto* frac : fractions) delete frac;
            delete data_hist;
            for (auto* h : mc_hists) delete h;
        } // End of energy bin loop

        // Write results for this detector to the output file
        output_file->cd();
        for (auto& h : h_best_fractions) 
        {
            h->Rebin(ProNbin);
            h->Write();
        }
        h_best_chi2->Rebin(ProNbin);
        h_best_chi2->Write();
        h_best_entries->Rebin(ProNbin);
        h_best_entries->Write();
    } // End of detector loop

    // 5. ============================ Finalization ============================
    canvas->Print((outputPdfPath + "]").c_str());
    delete canvas;
    output_file->Close();
    cout << "\n[SUCCESS] Fitting process completed. Results saved to:\n" << outputPdfPath << "\n" << outputRootPath << endl;
}

// Entry point to run the analysis
void HistTempFit_hist() {
    HistTempFit("Be", 7, true, false, 1, 1);
    HistTempFit("Be", 7, true, false, 1, 1, true, "Boron");
    HistTempFit("Be", 7, true, false, 1, 1, true, "Carbon");
    HistTempFit("Be", 7, true, false, 1, 1, true, "Nitrogen");
    HistTempFit("Be", 7, true, false, 1, 1, true, "Oxygen");

}