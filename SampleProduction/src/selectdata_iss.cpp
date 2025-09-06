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
#include "BinningManager.h"
#include <numeric>
#include <algorithm>

using namespace AMS_Iso;

void selectdata::SetAnalyzer(IsotopeAnalyzer* analyzer) {
	analyzer_ = analyzer;
}

void selectdata::LoopISS() {
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
    
    auto& binningManager = BinningManager::GetInstance();

	std::map<std::string, std::map<std::string, TH1F*>> hist1F_chains;
	std::map<std::string, std::map<std::string, TH2F*>> hist2F_chains;
	
	std::vector<std::string> isotope_list;
	for (int i = 0; i < isotope->getIsotopeCount(); ++i) {
		isotope_list.push_back(isotope->getName() + std::to_string(isotope->getMass(i)));
	}
	
	const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
	const std::vector<std::string> sources = {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
	const std::vector<std::string> cut_groups = {"CutGroup1", "CutGroup2", "CutGroup3"};
	const std::vector<std::string> num_den = {"Num", "Den"};
	const std::vector<std::string> charge_types = {"L1QSignal", "L1QTemplate", "L2QTemplate"};

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
					std::string name = buildHistName(chain, "ISS.ID.H1", {{"detector", det}, {"isotope", iso}});
					hist1F_chains[chain]["event_" + det + "_" + iso] = histManager->GetHist1F(name);
				}
			}
			
			for (const auto& det : detectors) {
				for (const auto& iso : isotope_list) {
					std::string name = buildHistName(chain, "ISS.ID.H2", {{"detector", det}, {"isotope", iso}});
					hist2F_chains[chain]["mass_" + det + "_" + iso] = histManager->GetHist2F(name);
				}
				std::string name = buildHistName(chain, "ISS.ID.H3", {{"detector", det}});
				hist2F_chains[chain]["cutoff_" + det] = histManager->GetHist2F(name);
			}
			
			hist1F_chains[chain]["rich_beta_naf"] = histManager->GetHist1F(buildHistName(chain, "ISS.ID.H4a"));
			hist1F_chains[chain]["rich_beta_agl"] = histManager->GetHist1F(buildHistName(chain, "ISS.ID.H4b"));
			
			hist2F_chains[chain]["delta_beta_naf_trk"] = histManager->GetHist2F(buildHistName(chain, "ISS.ID.H5a"));
			hist2F_chains[chain]["delta_beta_agl_trk"] = histManager->GetHist2F(buildHistName(chain, "ISS.ID.H5b"));
			hist2F_chains[chain]["delta_beta_tof_naf"] = histManager->GetHist2F(buildHistName(chain, "ISS.ID.H6a"));
			hist2F_chains[chain]["delta_beta_tof_agl"] = histManager->GetHist2F(buildHistName(chain, "ISS.ID.H6b"));
			hist2F_chains[chain]["delta_beta_tof_naf_bg"] = histManager->GetHist2F(buildHistName(chain, "ISS.ID.H7a"));
			hist2F_chains[chain]["delta_beta_tof_agl_bg"] = histManager->GetHist2F(buildHistName(chain, "ISS.ID.H7b"));
			
			for (const auto& det : detectors) {
				for (const auto& src : sources) {
					for (const auto& charge_type : charge_types) {
						std::string name = buildHistName(chain, "ISS.BKG.H1", {{"source", src}, {"charge_type", charge_type}, {"detector", det}});
						hist2F_chains[chain]["charge_" + src + "_" + charge_type + "_" + det] = histManager->GetHist2F(name);
					}
					std::string name = buildHistName(chain, "ISS.BKG.H2", {{"source", src}, {"detector", det}});
					hist1F_chains[chain]["source_" + src + "_" + det] = histManager->GetHist1F(name);
					
					name = buildHistName(chain, "ISS.BKG.H3", {{"source", src}, {"detector", det}});
					hist1F_chains[chain]["frag_" + src + "_" + det] = histManager->GetHist1F(name);
					
					name = buildHistName(chain, "ISS.BKG.H4", {{"source", src}, {"detector", det}});
					hist2F_chains[chain]["frag_mass_" + src + "_" + det] = histManager->GetHist2F(name);
				}
			}
			
			hist1F_chains[chain]["exp_rig"] = histManager->GetHist1F(buildHistName(chain, "ISS.FLUX.H1"));
			for (const auto& det : detectors) {
				for (const auto& iso : isotope_list) {
					std::string name = buildHistName(chain, "ISS.FLUX.H2", {{"detector", det}, {"isotope", iso}});
					hist1F_chains[chain]["exp_" + det + "_" + iso] = histManager->GetHist1F(name);
				}
			}
			
			for (const auto& det : detectors) {
				for (const auto& iso : isotope_list) {
					for (const auto& cut_group : cut_groups) {
						for (const auto& nd : num_den) {
							std::string name = buildHistName(chain, "ISS.FLUX.H3", {{"cut_group", cut_group}, {"num_den", nd}, {"detector", det}, {"isotope", iso}});
							hist1F_chains[chain]["eff_" + cut_group + "_" + nd + "_" + det + "_" + iso] = histManager->GetHist1F(name);
						}
					}
				}
			}
		}
	};

	populateHistMaps();

	std::map<std::string, std::vector<std::vector<TH1F*>>> h_beta_exposure_chains;
	for (const auto& chain : active_chains) {
		h_beta_exposure_chains[chain].resize(isotope->getIsotopeCount());
		for(int is = 0; is < isotope->getIsotopeCount(); ++is) {
			std::string isoName = isotope->getName() + std::to_string(isotope->getMass(is));
			h_beta_exposure_chains[chain][is].resize(3);
			h_beta_exposure_chains[chain][is][0] = hist1F_chains[chain]["exp_TOF_" + isoName];
			h_beta_exposure_chains[chain][is][1] = hist1F_chains[chain]["exp_NaF_" + isoName];
			h_beta_exposure_chains[chain][is][2] = hist1F_chains[chain]["exp_AGL_" + isoName];
		}
	}

	std::cout << "Loading models... " << std::endl;
	ModelManager::init("../../model_data.root", "../../model_mc.root");
	std::cout << "model n:" << ModelManager::model[0][0].index_correction.GetEntries() << std::endl;
    
    // MODIFIED: 获取用于不同目的的、特定于同位素的专属分箱
    // 1. Beta分箱 (源于标准刚度分箱)，更通用，用于事件的binning
    const std::vector<double>& Rbins_beta_std_for_UseMass = binningManager.GetBetaBins(charge, UseMass);
    // 2. Ek/n分箱 (源于标准刚度分箱)
    const std::vector<double>& Rbins_ek_for_UseMass = binningManager.GetEkPerNucleonBins(charge, UseMass);
	
	std::vector<unsigned int> timeTag;

	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
        
		RTICut rti_cut(this);
		TrackerCut tracker_cut(this);
		TOFCut tof_cut(this);
		RICHCut rich_cut(this);

		double InnerRig = tracker_cut.getRigidity();
		double L1InnerRig = tracker_cut.getRigidity(1,2,2);
		double cutOffRig = rti_cut.getCutoffRigidity();
		double richBeta = rich_cut.getBeta();
		double TOFBeta = tof_cut.getBeta();

		double rigCut = Constants::SAFE_FACTOR_RIG * cutOffRig;

		if (!rti_cut.cutRTI().total) continue;
		
		if (std::find(timeTag.begin(), timeTag.end(), time[0]) == timeTag.end()) {
			float exposureTime = rti_cut.calculateExposure().value;

			for (const auto& chain : active_chains) {
				auto h_exp_rig = histManager->GetHist1F(chain + "_ISS.FLUX.H1_");
				if (h_exp_rig) {
					for (int ibin = 1; ibin <= Constants::RIGIDITY_BINS; ibin++) {
						if (Binning::RigidityBins[ibin - 1] >= rigCut) {
							for (int ib = ibin; ib <= Constants::RIGIDITY_BINS; ib++) {
								double currentExp = h_exp_rig->GetBinContent(ib);
								h_exp_rig->SetBinContent(ib, currentExp + exposureTime);
							}
							break;
						}
					}
				}
				
				for (int is = 0; is < isotope->getIsotopeCount(); is++) {
                    int current_mass = isotope->getMass(is);
                    const auto& rb_expo_current = binningManager.GetIsotopeBetaBins(current_mass);
					double betaCO = Tools::rigidityToBeta(cutOffRig, charge, current_mass, false);

					for (int ib = 0; ib < Constants::BETA_TYPES; ib++) {
						double betaCut = Detector::BetaTypes[ib].getSafetyFactor() * betaCO;
						for (size_t ibin = 1; ibin < rb_expo_current.size(); ibin++) {
							if (rb_expo_current[ibin - 1] >= betaCut) {
								for (size_t ip = ibin; ip < rb_expo_current.size(); ip++) {
                                    std::string iso_name_str = isotope->getName() + std::to_string(current_mass);
                                    std::string hist_name = chain + "_ISS.FLUX.H2_" + Detector::BetaTypes[ib].getName() + "_" + iso_name_str;
                                    TH1F* h_exp_beta = histManager->GetHist1F(hist_name);
									if (h_exp_beta) {
										double currentExp = h_exp_beta->GetBinContent(ip);
										h_exp_beta->SetBinContent(ip, currentExp + exposureTime);
									}
								}
								break;
							}
						}
					}
				}
			}
			timeTag.push_back(time[0]);
		}

		auto modiRichPos = rich_cut.getModifiedPosition(true);
		double modiRichX = modiRichPos[0];
		double modiRichY = modiRichPos[1];
		Rad rad = (rich_NaF) ? NAF : AGL;
		double rich_beta_corr = ModelManager::corrected_beta(
				richBeta, rad, run, charge, modiRichX, modiRichY,
				rich_theta, rich_phi, rich_usedm, rich_hit, 0);

		double corr_cali_richBeta = Tools::CorrectCalibrationBiasInData(rich_beta_corr, rich_NaF);

		if (jentry % 1000000 == 0) {
			printf("rich beta before corr=%.6f, after corr=%.6f, after cali_corr=%.6f\n", richBeta, rich_beta_corr, corr_cali_richBeta);
		}
		richBeta = corr_cali_richBeta;
		
		if (jentry % 1000000 == 0) {
			std::cout << "Processing entry " << jentry << "/" << nentries << std::endl;
		}

		double NaFBeta = rich_NaF ? richBeta : -1;
		double AGLBeta = !rich_NaF ? richBeta : -1;
		double NaFEk = rich_NaF ? Tools::betaToKineticEnergy(richBeta) : -1;
		double AGLEk = !rich_NaF ? Tools::betaToKineticEnergy(richBeta) : -1;
		double TOFEk = Tools::betaToKineticEnergy(TOFBeta);

		auto h_exp_rig_first = histManager->GetHist1F(active_chains[0] + "_ISS.FLUX.H1_");
		double binLowForCutoffRig = h_exp_rig_first ? h_exp_rig_first->GetBinLowEdge(h_exp_rig_first->FindBin(L1InnerRig)) : 0;
		bool beyondCutoffRig = binLowForCutoffRig > Constants::SAFE_FACTOR_RIG * cutOffRig;

		int tofBin = Tools::findBin(Rbins_beta_std_for_UseMass, TOFBeta);
		int richBin = Tools::findBin(Rbins_beta_std_for_UseMass, richBeta);
		
        bool beyondCutoffTOF = TOFBeta >= 1 || ((tofBin >= 0) &&
				Tools::isBeyondCutoff(Rbins_beta_std_for_UseMass[tofBin], cutOffRig, Detector::BetaTypes[0].getSafetyFactor(), analyzer_->getCharge(), UseMass, false));
		bool beyondCutoffNaF = richBeta >= 1 || ((richBin >= 0) &&
				Tools::isBeyondCutoff(Rbins_beta_std_for_UseMass[richBin], cutOffRig, Detector::BetaTypes[1].getSafetyFactor(), analyzer_->getCharge(), UseMass, false));
		bool beyondCutoffAGL = richBeta >= 1 || ((richBin >= 0) &&
				Tools::isBeyondCutoff(Rbins_beta_std_for_UseMass[richBin], cutOffRig, Detector::BetaTypes[2].getSafetyFactor(), analyzer_->getCharge(), UseMass, false));

		auto TrackerCutResult = tracker_cut.cutTracker(charge, true);
		auto TOFCutResult = tof_cut.cutTOF(charge, true);
		auto RICHCutResult = rich_cut.cutRICH(charge, true, true);

		if(TrackerCutResult.total) {
			if(RICHCutResult.total) {
				double invRichBeta = 1. / richBeta;
				
				if(beyondCutoffRig) {
					for (const auto& chain : active_chains) {
						if(rich_NaF && L1InnerRig > 100) {
                            TH1F* h_rich_beta = histManager->GetHist1F(chain + "_ISS.ID.H4a_");
							if (h_rich_beta) h_rich_beta->Fill(invRichBeta);
						} else if(!rich_NaF && L1InnerRig > 200) {
                            TH1F* h_rich_beta = histManager->GetHist1F(chain + "_ISS.ID.H4b_");
							if (h_rich_beta) h_rich_beta->Fill(invRichBeta);
						}
					}
				}
				
				if(beyondCutoffRig && TOFCutResult.total) {
					double deltaBeta = (1. / TOFBeta) - invRichBeta;
					for (const auto& chain : active_chains) {
						if(rich_NaF) {
                            TH2F* h_delta_beta = histManager->GetHist2F(chain + "_ISS.ID.H6a_");
							if (h_delta_beta) h_delta_beta->Fill(deltaBeta, NaFBeta*L1InnerRig);
						} else if(!rich_NaF) {
                            TH2F* h_delta_beta = histManager->GetHist2F(chain + "_ISS.ID.H6b_");
							if (h_delta_beta) h_delta_beta->Fill(deltaBeta, AGLBeta*L1InnerRig);
						}
					}
				}
			}
		}

		if(tracker_cut.cutTracker(6, true).total && beyondCutoffRig){
            // Potential future logic here
		}
	}

	std::cout << "Event processing completed" << std::endl;
}

