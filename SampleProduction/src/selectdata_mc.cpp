#define selectdata_cxx
#include "selectdata.h"
#include "IsotopeAnalyzer.h"
#include "ModelManager.h"
#include "RTICut.h"
#include "TrackerCut.h"
#include "TOFCut.h"
#include "RICHCut.h"
#include "Tool.h"
#include "ProductRegistry.h"
#include "HistManager.h"
#include "BinningManager.h" // <-- Include the BinningManager

using namespace AMS_Iso;

void selectdata::SetAnalyzer(IsotopeAnalyzer* analyzer) {
	analyzer_ = analyzer;
}

void selectdata::LoopMC() {
	if (!fChain || !analyzer_) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout << "Total entries: " << nentries << std::endl;

	int charge = analyzer_->getCharge();
	int UseMass = analyzer_->getUseMass();
	std::cout << "Using Mass = " << UseMass << " for binning and ek convert" << std::endl;

	const IsotopeVar* isotope = analyzer_->getIsotope();
	int tempo_Mass = isotope->getMass(0);
	std::cout << "Using tempo_Mass = " << tempo_Mass << " for cutoff convert" << std::endl;

	auto histManager = analyzer_->getHistManager();
	if (!histManager) {
		std::cerr << "Failed to get HistManager" << std::endl;
		return;
	}

	std::vector<std::string> active_chains = analyzer_->getActiveChains();
	if (active_chains.empty()) {
		std::cerr << "No active chains found" << std::endl;
		return;
	}

    // --- NEW: Get an instance of the BinningManager ---
    auto& binningManager = BinningManager::GetInstance();

	std::map<std::string, std::map<std::string, TH1F*>> hist1F_chains;
	std::map<std::string, std::map<std::string, TH2F*>> hist2F_chains;
	
	std::vector<std::string> isotope_list;
	isotope_list.push_back(isotope->getName() + std::to_string(UseMass));
	
	const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
	const std::vector<std::string> cut_groups = {"CutGroup1", "CutGroup2", "CutGroup3"};
	const std::vector<std::string> num_den = {"Num", "Den"};
	const std::vector<std::string> gen_rec = {"gen", "rec"};

	auto buildHistName = [&](const std::string& chain, const std::string& template_id, const std::map<std::string, std::string>& params = {}) -> std::string {
		std::string name = chain + "_" + template_id;
		for (const auto& [key, value] : params) {
			name += "_" + value;
		}
		return name;
	};

	auto populateHistMaps = [&]() {
		for (const auto& chain : active_chains) {
			hist1F_chains[chain] = std::map<std::string, TH1F*>();
			hist2F_chains[chain] = std::map<std::string, TH2F*>();

			for (const auto& det : detectors) {
				for (const auto& iso : isotope_list) {
					std::string name = buildHistName(chain, "MC.ID.H1a", {{"detector", det}, {"isotope", iso}});
					hist2F_chains[chain]["mass_sel_" + det + "_" + iso] = histManager->GetHist2F(name);
					
					name = buildHistName(chain, "MC.ID.H1b", {{"detector", det}, {"isotope", iso}});
					hist2F_chains[chain]["mass_sel_notoffrag_" + det + "_" + iso] = histManager->GetHist2F(name);
				}
				
				std::string name = buildHistName(chain, "MC.ID.H2", {{"detector", det}});
				hist2F_chains[chain]["cutoff_" + det] = histManager->GetHist2F(name);
			}
			
			hist1F_chains[chain]["rich_beta_naf"] = histManager->GetHist1F(buildHistName(chain, "MC.ID.H3a"));
			hist1F_chains[chain]["rich_beta_agl"] = histManager->GetHist1F(buildHistName(chain, "MC.ID.H3b"));
			
			hist2F_chains[chain]["delta_beta_naf_trk"] = histManager->GetHist2F(buildHistName(chain, "MC.ID.H4a"));
			hist2F_chains[chain]["delta_beta_agl_trk"] = histManager->GetHist2F(buildHistName(chain, "MC.ID.H4b"));
			hist2F_chains[chain]["delta_beta_tof_naf"] = histManager->GetHist2F(buildHistName(chain, "MC.ID.H5a"));
			hist2F_chains[chain]["delta_beta_tof_agl"] = histManager->GetHist2F(buildHistName(chain, "MC.ID.H5b"));
			hist2F_chains[chain]["delta_beta_tof_naf_bg"] = histManager->GetHist2F(buildHistName(chain, "MC.ID.H6a"));
			hist2F_chains[chain]["delta_beta_tof_agl_bg"] = histManager->GetHist2F(buildHistName(chain, "MC.ID.H6b"));
			
			for (const auto& det : detectors) {
				hist1F_chains[chain]["source_" + det] = histManager->GetHist1F(buildHistName(chain, "MC.BKG.H1", {{"detector", det}}));
				
				for (const auto& iso : isotope_list) {
					std::string name = buildHistName(chain, "MC.BKG.H2", {{"detector", det}, {"isotope", iso}});
					hist1F_chains[chain]["frag_" + det + "_" + iso] = histManager->GetHist1F(name);
					
					name = buildHistName(chain, "MC.BKG.H3a", {{"detector", det}, {"isotope", iso}});
					hist1F_chains[chain]["uptof_frag_" + det + "_" + iso] = histManager->GetHist1F(name);
					
					name = buildHistName(chain, "MC.BKG.H3b", {{"detector", det}, {"isotope", iso}});
					hist1F_chains[chain]["uptof_frag_rich_" + det + "_" + iso] = histManager->GetHist1F(name);
				}
			}
			
			for (const auto& cut_group : cut_groups) {
				for (const auto& gr : gen_rec) {
					std::string name = buildHistName(chain, "MC.FLUX.H1", {{"cut_group", cut_group}, {"gen_rec", gr}});
					hist1F_chains[chain]["passed_" + cut_group + "_" + gr] = histManager->GetHist1F(name);
				}
			}
			
			for (const auto& det : detectors) {
				for (const auto& iso : isotope_list) {
					for (const auto& cut_group : cut_groups) {
						for (const auto& nd : num_den) {
							std::string name = buildHistName(chain, "MC.FLUX.H2", {{"cut_group", cut_group}, {"num_den", nd}, {"detector", det}, {"isotope", iso}});
							hist1F_chains[chain]["eff_" + cut_group + "_" + nd + "_" + det + "_" + iso] = histManager->GetHist1F(name);
						}
					}
				}
			}
			
			hist1F_chains[chain]["mc_total"] = histManager->GetHist1F(buildHistName(chain, "MC.FLUX.H3"));
		}
	};

	populateHistMaps();

	std::string currentIsotope = isotope_list[0];
	
	std::map<std::string, std::vector<TH1F*>> h_cuts_chains[3];
	for (const auto& chain : active_chains) {
		h_cuts_chains[0][chain].resize(11);
		h_cuts_chains[1][chain].resize(11);
		h_cuts_chains[2][chain].resize(11);
		
		for(int i = 0; i <= 10; ++i) {
			int group_idx = (i / 4) % 3;
			std::string cut_group = cut_groups[group_idx];
			std::string nd = (i % 2 == 0) ? "Den" : "Num";
			
			h_cuts_chains[0][chain][i] = hist1F_chains[chain]["eff_" + cut_group + "_" + nd + "_TOF_" + currentIsotope];
			h_cuts_chains[1][chain][i] = hist1F_chains[chain]["eff_" + cut_group + "_" + nd + "_NaF_" + currentIsotope];
			h_cuts_chains[2][chain][i] = hist1F_chains[chain]["eff_" + cut_group + "_" + nd + "_AGL_" + currentIsotope];
		}
	}

	std::cout << "Loading models... " << std::endl;
    // --- MODIFIED: Using relative paths ---
	ModelManager::init("../../model_data.root", "../../model_mc.root");
	std::cout << "model n:" << ModelManager::model[0][0].index_correction.GetEntries() << std::endl;

	UInt_t current_run = 0;
	UInt_t min_event = 0;
	UInt_t max_event = 0;
	int event_count = 1;
	std::vector<double> mc_events;
	
	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		
		if (current_run != run) {
			if (current_run != 0) {
				mc_events.push_back(max_event - min_event + 1 +
						(max_event - min_event + 1) / event_count);
			}
			current_run = run;
			min_event = event;
			max_event = event;
			event_count = 1;
		} else {
			min_event = std::min(min_event, event);
			max_event = std::max(max_event, event);
			event_count++;

			if (jentry == nentries - 1) {
				mc_events.push_back(max_event - min_event + 1 +
						(max_event - min_event + 1) / event_count);
			}
		}

		RTICut rti_cut(this);
		TrackerCut tracker_cut(this);
		TOFCut tof_cut(this);
		RICHCut rich_cut(this);

		double InnerRig = tracker_cut.getRigidity();
		double L1InnerRig = tracker_cut.getRigidity(1,2,2);
		double cutOffRig = rti_cut.getCutoffRigidity();
		double richBeta = rich_cut.getBeta();
		double TOFBeta = tof_cut.getBeta();

		auto modiRichPos = rich_cut.getModifiedPosition(true);
		double modiRichX = modiRichPos[0];
		double modiRichY = modiRichPos[1];
		Rad rad = (rich_NaF) ? NAF : AGL;
		double rich_beta_corr = ModelManager::corrected_beta(
				richBeta, rad, run, charge, modiRichX, modiRichY,
				rich_theta, rich_phi, rich_usedm, rich_hit, 1);

		if (jentry % 1000000 == 0) {
			printf("rich beta before corr=%.6f, after corr=%.6f\n", richBeta, rich_beta_corr);
		}
		richBeta = rich_beta_corr;
		
		if (jentry % 1000000 == 0) {
			std::cout << "Processing entry " << jentry << "/" << nentries << std::endl;
		}

		double NaFBeta = rich_NaF ? richBeta : -1;
		double AGLBeta = !rich_NaF ? richBeta : -1;
		double NaFEk = rich_NaF ? Tools::betaToKineticEnergy(richBeta) : -1;
		double AGLEk = !rich_NaF ? Tools::betaToKineticEnergy(richBeta) : -1;
		double TOFEk = Tools::betaToKineticEnergy(TOFBeta);

		auto TrackerCutResult = tracker_cut.cutTracker(charge, false);
		auto TOFCutResult = tof_cut.cutTOF(charge, false);
		auto RICHCutResult = rich_cut.cutRICH(charge, false, true);

		double MC_weight_27 = Tools::calculateWeight(mmom, mch, false);

		if(TrackerCutResult.total) {
			if(RICHCutResult.total) {
				double invRichBeta = 1. / richBeta;
				
				if(L1InnerRig > 150) {
					for (const auto& chain : active_chains) {
						if(rich_NaF) {
							if (hist1F_chains[chain]["rich_beta_naf"]) hist1F_chains[chain]["rich_beta_naf"]->Fill(invRichBeta);
						} else {
							if (hist1F_chains[chain]["rich_beta_agl"]) hist1F_chains[chain]["rich_beta_agl"]->Fill(invRichBeta);
						}
					}
				}
				
				if(L1InnerRig > 150 && TOFCutResult.total) {
					double deltaBeta = (1. / TOFBeta) - invRichBeta;
					for (const auto& chain : active_chains) {
						if(rich_NaF) {
							if (hist2F_chains[chain]["delta_beta_tof_naf"]) hist2F_chains[chain]["delta_beta_tof_naf"]->Fill(deltaBeta, NaFBeta*L1InnerRig);
						} else {
							if (hist2F_chains[chain]["delta_beta_tof_agl"]) hist2F_chains[chain]["delta_beta_tof_agl"]->Fill(deltaBeta, AGLBeta*L1InnerRig);
						}
					}
				}
			}
		}

		double generatedRig = (mch != 0) ? (mmom/mch) : 0;
		double generatedEk = Tools::rigidityToKineticEnergy(generatedRig, mch, UseMass);

		bool passTracker = true;
		for (size_t i = 0; i < 7; ++i) {
			passTracker &= TrackerCutResult.details[i+1];
			if(passTracker) {
				for (const auto& chain : active_chains) {
					if (h_cuts_chains[0][chain][i]) h_cuts_chains[0][chain][i]->Fill(generatedEk);
					if (h_cuts_chains[1][chain][i]) h_cuts_chains[1][chain][i]->Fill(generatedEk);
					if (h_cuts_chains[2][chain][i]) h_cuts_chains[2][chain][i]->Fill(generatedEk);
				}
			}
		}

		bool passTOF = passTracker;
		for (size_t i = 0; i < 2; ++i) {
			passTOF &= TOFCutResult.details[i];
			if(passTOF) {
				for (const auto& chain : active_chains) {
					if (h_cuts_chains[0][chain][i+7]) h_cuts_chains[0][chain][i+7]->Fill(generatedEk);
				}
			}
		}

		bool passRICH = passTracker;
		for (size_t i = 0; i < 3; ++i) {
			passRICH &= RICHCutResult.details[i];
			if(passRICH && rich_NaF) {
				for (const auto& chain : active_chains) {
					if (h_cuts_chains[1][chain][i+7]) h_cuts_chains[1][chain][i+7]->Fill(generatedEk);
				}
			}
			if(passRICH && !rich_NaF) {
				for (const auto& chain : active_chains) {
					if (h_cuts_chains[2][chain][i+7]) h_cuts_chains[2][chain][i+7]->Fill(generatedEk);
				}
			}
		}
	}

	double total_events = std::accumulate(mc_events.begin(), mc_events.end(), 0.0);
	for (const auto& chain : active_chains) {
		auto h_mc_num = hist1F_chains[chain]["mc_total"];
		if (h_mc_num) {
			h_mc_num->SetBinContent(1, total_events - 2);
			std::cout << "Total MC Number for " << chain << ": " << h_mc_num->GetBinContent(1) << std::endl;
		}
	}

	std::cout << "Event processing completed" << std::endl;
}