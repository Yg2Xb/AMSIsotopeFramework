#define selectdata_cxx
#include "selectdata.h"
#include "basic_var.h"
#include "IsotopeAnalyzer.h"
#include "ModelManager.h"
#include "RTICut.h"
#include "TrackerCut.h"
#include "TOFCut.h"
#include "RICHCut.h"
#include "Tool.h"
#include "HistManager.h"
#include "BinningManager.h"

#include <numeric>
#include <algorithm>
#include <iostream>
#include <map>
#include <unordered_map>
#include <cmath>

using namespace AMS_Iso;

void selectdata::SetAnalyzer(IsotopeAnalyzer* a){ analyzer_ = a; }

void selectdata::Loop() {
	if (!fChain || !analyzer_) {
		std::cerr << "[FATAL] fChain or analyzer_ is null" << std::endl;
		return;
	}

	const bool isISS = analyzer_->isISS();
	const int charge = analyzer_->getCharge();
	const int UseMass = analyzer_->getUseMass();
	auto* histManager = analyzer_->getHistManager();
	bool forBackground = analyzer_->isNoBkgCut();
	int FragmentZ = analyzer_->getFragmentZ();
	if (!histManager) { std::cerr << "Failed to get HistManager\n"; return; }

	TTree* filteredTree = histManager->GetFilteredTree();

	auto& binMgr = BinningManager::GetInstance();
	auto StdBetaBins = binMgr.GetBetaBins(2, 4);
	auto* bins = histManager->IDH5a[0][0]->GetYaxis()->GetXbins();
	std::vector<double> StdRigBins(bins->GetArray(), bins->GetArray() + bins->GetSize());
	const IsotopeVar* iso = analyzer_->getIsotope();
	if (!iso) { std::cerr << "[FATAL] IsotopeVar is null" << std::endl; return; }
	std::vector<std::string> chains = analyzer_->getActiveChains();
	if (chains.empty()) { std::cerr << "No active chains found\n"; return; }

	const Long64_t nentries = fChain->GetEntries();
	std::cout << "Total entries: " << nentries << " | isISS=" << isISS << " | Z=" << charge << " | A=" << UseMass << std::endl;

	// Initialize models and flux helpers
	ModelManager::init("/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/model_data.root",
			"/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/model_mc.root");
	if(!isISS)
	{
		AMS_Iso::Tools::initFluxFunctions();
	}
	if(isISS)
	{
		//AMS_Iso::Tools::initChargeTuning();
	}

	// Book-keeping for ISS exposure and MC generation accounting
	std::vector<unsigned int> timeTag;
	UInt_t current_run = 0, min_event = 0, max_event = 0; int event_count = 1;
	std::vector<double> mc_events;

	// Dimensions
	const int NchainLoc = static_cast<int>(chains.size()); // chain index
	const int NdetLoc   = 3;                               // 0:TOF, 1:NaF, 2:AGL
	const int NisoLoc   = iso->getIsotopeCount();

	// ISS sources (parents) for background study: Z = 4..8 (Be,B,C,N,O)
	const std::vector<int> source_Z = {2,3,4,5,6,7,8};
	const int NsrcLoc = static_cast<int>(source_Z.size());
	const std::vector<float> minRig_NaF = {120,70,80,100,70};
	const std::vector<float> minRig_AGL = {120,70,80,100,70};

	auto getMinAForZ = [&](int Z)->int {
		int Amin = -1;
		for (int n = 0; n < Constants::N_nuc; ++n) {
			if (Constants::nuclei_Z[n] != Z) continue;
			int A = Constants::nuclei_A[n];
			if (Amin < 0 || A < Amin) Amin = A;
		}
		return Amin;
	};

	// Working buffers
	double weight_NucFlux = 1.;
	double cutOffRig = -1, TOFBeta = -1, richBeta = -1, RecCharge_Int = -1;
	double rig_chain[2] = {-1, -1};
	double beta_det[3] = {-1, -1, -1};
	double ek_det[3] = {-1, -1, -1};

	bool beyondRigCutoff = false;
	// Beyond-cutoff masks per detector and per (Z,A) in Constants lists
	bool beyondBetaCutoff[3][Constants::N_nuc] = {};
	std::map<std::pair<int,int>, int> ZAMap;
	for (int n = 0; n < Constants::N_nuc; ++n)
		ZAMap[{Constants::nuclei_Z[n], Constants::nuclei_A[n]}] = n;

	auto getBeyondBetaCutoffCut = [&](int det, int Z, int A){
		if (det < 0 || det >= NdetLoc) return false;
		auto it = ZAMap.find({Z, A});
		if (it == ZAMap.end()) return false;
		return beyondBetaCutoff[det][it->second];
	};

	// MC: single input-mother GeneID (as discussed)
	const int geneID_MC = analyzer_->getGeneID(charge, UseMass);

	// Global fragZ and its fragment IDs (defines NisoBKG)
	const int fragZ = isISS ? charge : FragmentZ;
	const std::vector<int> fragIDs_global = analyzer_->getBkgFragIDs(fragZ);
	const int NisoBKG = static_cast<int>(fragIDs_global.size());
	std::unordered_map<int,int> fragID_to_index;
	for (int i = 0; i < NisoBKG; ++i) fragID_to_index[fragIDs_global[i]] = i;

	// Fragment isotope slots (align with HistManager booking order)
	const std::vector<int> FragA = [&]() {
		const auto& weights = isotopeWeights.at(fragZ);
		std::vector<int> result(weights.size());
		
		std::transform(weights.begin(), weights.end(), result.begin(), 
			[](const auto& p) { return p.first; });
		
		return result;
	}();


	if (!isISS) {
		std::cout << "[INFO] MC mode: geneID_MC=" << geneID_MC << " | fragZ=" << fragZ << " | NisoBKG=" << NisoBKG << std::endl;
	}

	for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) {
			std::cerr << "[WARN] LoadTree returned <0 at entry " << jentry << std::endl;
			break;
		}
		fChain->GetEntry(jentry);

		if (jentry % 1000000 == 0)
			std::cout << "Processing entry " << jentry << "/" << nentries << std::endl;

		// MC run-event accounting (keep original style)
		if (!isISS) {
			if (current_run != run) {
				if (current_run != 0) {
					mc_events.push_back(max_event - min_event + 1 + (max_event - min_event + 1) / event_count);
				}
				current_run = run; min_event = event; max_event = event; event_count = 1;
			} else {
				min_event = std::min(min_event, event);
				max_event = std::max(max_event, event);
				event_count++;
				if (jentry == nentries - 1) {
					mc_events.push_back(max_event - min_event + 1 + (max_event - min_event + 1) / event_count);
				}
			}
		}

		// Build cut helpers
		RTICut rti_cut(this);
		TrackerCut tracker_cut(this);
		TOFCut tof_cut(this);
		RICHCut rich_cut(this);

		// RTI cut for ISS
		if (isISS && !rti_cut.cutRTI().total) continue;

		// Basic kinematics
		RecCharge_Int = tk_qin[0][2] > 0 ? static_cast<int>(tk_qin[0][2] + 0.5f) : -1; //use rec charge for both data and MC, never use mc info if not necessary
		rig_chain[0] = tracker_cut.getRigidity();            // Inner
		rig_chain[1] = tracker_cut.getRigidity(1,2,2);       // L1Inner
		if(rig_chain[0] > 0){
			rig_chain[0] = 1./(1./rig_chain[0] - 1./35000.); // Apply rigidity bias correction
		}
		if(rig_chain[1] > 0){
			rig_chain[1] = 1./(1./rig_chain[1] - 1./35000.); // Apply rigidity bias correction
		}
		cutOffRig = rti_cut.getCutoffRigidity();
		richBeta = yanzx_dst ? rich_cut.getBeta(1) : rich_cut.getBeta(0);
		TOFBeta = tof_cut.getBeta();

		// Exposure filling for ISS (only the two flux hists)
		//bool isExpoTcount = false;
		if (isISS && (std::find(timeTag.begin(), timeTag.end(), time[0]) == timeTag.end())) {
			//isExpoTcount = true;
			float exposureTime = rti_cut.calculateExposure().value;

			// [FILL] ISS.FLUX.H2 (exposure time vs rigidity threshold)
			if (auto* h_exp_rig = histManager->ISS_FLUXH2[0].get()) {
				double rigCut = Constants::SAFE_FACTOR_RIG * cutOffRig;
				for (int ibin = 1; ibin <= h_exp_rig->GetNbinsX(); ++ibin) {
					if (h_exp_rig->GetBinLowEdge(ibin) >= rigCut) {
						for (int ib = ibin; ib <= h_exp_rig->GetNbinsX(); ++ib)
							h_exp_rig->SetBinContent(ib, h_exp_rig->GetBinContent(ib) + exposureTime);
						break;
					}
				}
			} else {
				std::cerr << "[ERROR] ISS_FLUXH2[0] is null" << std::endl;
			}

			// [FILL] ISS.FLUX.H3 (exposure time vs E_k/n) for each detector and isotope
			for (int d = 0; d < NdetLoc; ++d) {
				for (int i = 0; i < NisoLoc; ++i) {
					int mass = iso->getMass(i);
					double betaCO = Tools::rigidityToBeta(cutOffRig, charge, mass, false);
					double betaCut = Detector::BetaTypes[d].getSafetyFactor() * betaCO;
					if (auto* h = histManager->ISS_FLUXH3[d][i].get()) {
						for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
							if (h->GetBinLowEdge(ibin) >= betaCut) {
								for (int ip = ibin; ip <= h->GetNbinsX(); ++ip)
									h->SetBinContent(ip, h->GetBinContent(ip) + exposureTime);
								break;
							}
						}
					} else {
						std::cerr << "[ERROR] ISS_FLUXH3["<<d<<"]["<<i<<"] is null" << std::endl;
					}
				}
			}
			timeTag.push_back(time[0]);
		}

		// RICH beta correction (data) and smearing (MC)
		// 计算用于校正的电荷（优先级：Inner > TOF > L1 > charge）
		int Charge_forCorr = -1;

		if (RecCharge_Int > 0) {
			Charge_forCorr = RecCharge_Int;
		} else {
			// 计算 TOF 下两层的平均电荷
			double lowtofQmean = (tof_ql[2] + tof_ql[3]) / 2.0;
			
			// 如果其中一层失效（相乘为0），将均值乘2恢复有效层的值
			if (tof_ql[2] * tof_ql[3] == 0) {
				lowtofQmean *= 2.0;
			}
			
			// 按优先级尝试取整
			if (lowtofQmean > 0) {
				Charge_forCorr = static_cast<int>(lowtofQmean + 0.5f);
			} else if (tk_qln[0][7][2] > 0) {
				Charge_forCorr = static_cast<int>(tk_qln[0][7][2] + 0.5f);
			} else {
				Charge_forCorr = charge;
			}
		}
		double rich_beta_corr = richBeta;
		if (!yanzx_dst && richBeta > 0) {
			auto modiRichPos = rich_cut.getModifiedPosition(true);
			double modiRichX = modiRichPos[0];
			double modiRichY = modiRichPos[1];
			Rad rad = (rich_NaF) ? NAF : AGL;
			int is_mc = (isISS) ? 0 : 1;
			
			rich_beta_corr = ModelManager::corrected_beta(
					richBeta, rad, run, Charge_forCorr, modiRichX, modiRichY,
					rich_theta, rich_phi, rich_usedm, rich_hit, is_mc);
		}
		double corr_cali_richBeta = isISS ? Tools::CorrectCalibrationBiasInData(rich_beta_corr, rich_NaF) : rich_beta_corr;
		richBeta = corr_cali_richBeta; 
		if(!isISS){
			richBeta = Tools::GetSmearRichBeta(Charge_forCorr, richBeta, rich_NaF);
		}
		double generatedRig = mmom / mch;
		double generatedEk_n = Tools::rigidityToKineticEnergy(generatedRig, charge, UseMass);
		double generatedBeta = Tools::rigidityToBeta(generatedRig, charge, UseMass, false);
		double L2TruthRig = mtrmom[1] / mch;
		double L2TruthEk_n = Tools::rigidityToKineticEnergy(L2TruthRig, charge, UseMass);

		// Per-detector beta and Ek/n
		beta_det[0] = TOFBeta;
		beta_det[1] = rich_NaF ? richBeta : -9;
		beta_det[2] = !rich_NaF ? richBeta : -9;

		ek_det[0] = Tools::betaToKineticEnergy(TOFBeta);
		ek_det[1] = rich_NaF ? Tools::betaToKineticEnergy(richBeta) : -9;
		ek_det[2] = (!rich_NaF) ? Tools::betaToKineticEnergy(richBeta) : -9;


		// Event weight
		weight_NucFlux = isISS ? 1.0 : Tools::calculateWeight(mmom, charge, UseMass, isISS);
		if(!isISS && jentry % 100000 == 0) {
			std::cout << "Momentum: " << mmom << " GeV, Charge: " << mch << std::endl;
			std::cout << "Event weight (NucFlux) = " << weight_NucFlux << std::endl;
		}


		// Tracker two-acc selection per chain
		auto twoAcc = tracker_cut.TwoAccTrackerCut(charge, isISS, false);
		bool TwoAccTrackerCutResult[2] = { twoAcc.details[0], twoAcc.details[1] };

		// Detector beta quality selection for ID/BKG usage
		bool BetaDetectorCutResult[3] = {//charge-dependent
			tof_cut.cutTOF(charge, isISS).total,
			rich_NaF && rich_cut.cutRICH(charge, isISS, true).total,
			!rich_NaF && rich_cut.cutRICH(charge, isISS, true).total
		};
		// Detector validity for BKG (quality-only, full version) !!!2025.Nov.3 check how beta chi2 cut influent frag 
		bool detValidBkg[3] = { //charge-independent
			tof_cut.cutTOFforBkg(charge, isISS).details[1], //no chi2 cut, only geo cut
			rich_NaF && rich_cut.cutRICHforBkg(charge, isISS, true).total,
			!rich_NaF && rich_cut.cutRICHforBkg(charge, isISS, true).total
		};

		if (rig_chain[1] < 0.8 || rig_chain[1] > 3300) {
			beyondRigCutoff = false;
		} else {
			int rigBin = Tools::findBin(StdRigBins, rig_chain[1]);
			float rigbinlow = StdRigBins[rigBin];
			beyondRigCutoff = isISS ? (rigbinlow > Constants::SAFE_FACTOR_RIG * cutOffRig) : true;
		}
		//std::cout<<"rig:"<<rig_chain[0]<<" rigbin:"<<rigBin<<" rigbinlow:"<<rigbinlow<<" cutoffrig:"<<cutOffRig<<" cutoff:"<<beyondRigCutoff<<std::endl;
		// directly beyond strictest beta cutoff, Be(4) to O(8)
		/*
		bool beyondBetaCutoff_direct[NdetLoc][5] = {};
		int Amin_forZ[5];
		for (int z = 4; z <= 8; ++z) {
			Amin_forZ[z - 4] = getMinAForZ(z);
		}
		for (int d = 0; d < NdetLoc; ++d) { // beta < 0 means false
			if (!(beta_det[d] > 0) || !std::isfinite(beta_det[d])) {
				for (int k = 0; k < 5; ++k) {
					beyondBetaCutoff_direct[d][k] = false;
				}
				continue;
			}
			for (int z = 4; z <= 8; ++z) {
				const int idx = z - 4; // 0..4 对应 Z=4..8
				const int Amin_src = Amin_forZ[idx];
				if (Amin_src <= 0) {
					beyondBetaCutoff_direct[d][idx] = false;
					continue;
				}
				const double betaCO  = Tools::rigidityToBeta(cutOffRig, z, Amin_src, false);
				const double betaCut = Detector::BetaTypes[d].getSafetyFactor() * betaCO;
				beyondBetaCutoff_direct[d][idx] = (beta_det[d] > betaCut); //attention! using normal beta
			}
		}
		*/

		// Precompute beyondBetaCutoff mask for all (Z,A) in Constants lists, attention! using normal beta!
		for (int d = 0; d < NdetLoc; ++d) {
			int betaBin = Tools::findBin(StdBetaBins, beta_det[d]);
			for (int n = 0; n < Constants::N_nuc; ++n) {
				int Z = Constants::nuclei_Z[n]; 
				int A = Constants::nuclei_A[n];
				//auto StdBetaBins = binMgr.GetBetaBins(Z, A);
				//union ek bin == union std beta bin
				if (!isISS || beta_det[d] >= 1) { beyondBetaCutoff[d][n] = true; continue; }
				if (betaBin >= 0) {
					double betaLow = StdBetaBins[betaBin];
					beyondBetaCutoff[d][n] = Tools::isBeyondCutoff(
							betaLow, cutOffRig, Detector::BetaTypes[d].getSafetyFactor(), Z, A, !isISS
							);
				} else {
					beyondBetaCutoff[d][n] = false;
				}
			}
		}

		double tk_ql1_unbiased = tk_exqln[Tracker::ChargeReco::DEFAULT][0][Tracker::Direction::DEFAULT];
		double tk_ql1 = tk_qln[Tracker::ChargeReco::DEFAULT][0][Tracker::Direction::DEFAULT];
		double tk_ql2 = tk_qln[Tracker::ChargeReco::DEFAULT][1][Tracker::Direction::DEFAULT];
		
		//-------------------filter tree begin--------------------
		auto l1n = tracker_cut.cutL1Norm(charge, isISS);
		auto l1u = tracker_cut.cutL1Unbiased(charge, isISS);
		//fill tree
		if (filteredTree) {
			bool L1Norm_Need = l1n.details[2] && l1n.details[3] && l1n.details[4] && tk_ql1 > 2.3 && tk_ql1 < 8.8; 
    		bool L1Unb_Need  = l1u.details[2] && l1u.details[3] && tk_ql1_unbiased > 2.3 && tk_ql1_unbiased < 8.8; 
            bool passFilter = tracker_cut.Q_L1_BkgIndependCut(charge, isISS) && (L1Norm_Need || L1Unb_Need);
            if (passFilter && (detValidBkg[0] || detValidBkg[1] || detValidBkg[2])) {
                filteredTree->Fill();
            }
        }
		//-------------------filter tree end--------------------

		// =========================
		// ID histograms 
		// =========================
		for (int c = 0; c < NchainLoc; ++c) {
			if (c >= 2) continue;
			if (!TwoAccTrackerCutResult[c]) continue;

			for (int d = 0; d < NdetLoc; ++d) {
				if (!BetaDetectorCutResult[d]) continue;

				auto mres = Tools::calculateMass(beta_det[d], 1.0, rig_chain[c], RecCharge_Int);
				if (mres.invMass <= 0) continue;

				// [FILL] ISS.ID.H1
				if (isISS) {
					for (int i = 0; i < NisoLoc; ++i) {
						const int A = iso->getMass(i);
						if (!getBeyondBetaCutoffCut(d, charge, A)) continue;
						histManager->ISS_IDH1[c][d][i].get()->Fill(ek_det[d], weight_NucFlux);
					}
				}

				for (int i = 0; i < NisoLoc; ++i) {
					const int A = iso->getMass(i);
					if (!getBeyondBetaCutoffCut(d, charge, A)) continue;

					// [FILL] ID.H2
					histManager->IDH2[c][d][i].get()->Fill(mres.invMass, ek_det[d], weight_NucFlux);

					// [FILL] MC.ID.H1
					if (!isISS && geneID_MC != -1 && mtrpar[1] == geneID_MC) {
						histManager->MC_IDH1[c][d][i].get()->Fill(mres.invMass, ek_det[d], weight_NucFlux);
					}
				}
			}
		}

		// =========================beta study histograms =========================
		int nZ = isISS ? 8 - 2 + 1 : 1;
		double AveMass = isISS ? 0 : UseMass;
		if (isISS) {
			auto it = isotopeWeights.find(RecCharge_Int);
			if (it != isotopeWeights.end()) {
				double totalWeight = 0.0;
				for (const auto& [mass, weight] : it->second) {
					AveMass += mass * weight;
					totalWeight += weight;
				}
				if (totalWeight > 0) AveMass /= totalWeight;
			}
		}

		//sample select
		for (int iz = 0; iz < nZ; ++iz) {
			const int z_use = isISS ? source_Z[iz] : charge;
			const int A_use = isISS ? getMinAForZ(z_use) : UseMass;
			if (A_use < 0) continue;
			
			auto twoAcc_use = tracker_cut.TwoAccTrackerCut(z_use, isISS, false);
			bool TwoAccTrackerCutResult_use[2] = {twoAcc_use.details[0], twoAcc_use.details[1] };	
			
			bool BetaDetectorCutResult_use[3] = {
				tof_cut.cutTOF(z_use, isISS).total,
				rich_NaF && rich_cut.cutRICH(z_use, isISS, true).total,
				!rich_NaF && rich_cut.cutRICH(z_use, isISS, true).total
			};

			for (int c = 0; c < std::min(2, NchainLoc); ++c) {
				if (!TwoAccTrackerCutResult_use[c]) continue;
				
				double beta_tracker = 0.0;
				//if frag, we don't know A, need further study
				if (rig_chain[c] > 0 && AveMass > 0 && RecCharge_Int > 0) {
					beta_tracker = isISS ? Tools::rigidityToBeta(rig_chain[c], z_use, AveMass, false) : Tools::rigidityToBeta(rig_chain[c], charge, UseMass, false);
				}

				//[FILL] ID.H4 1/richbeta
				if (BetaDetectorCutResult_use[1] && rig_chain[c] > 80) {
					histManager->IDH4a[iz][c].get()->Fill(1.0 / beta_det[1], weight_NucFlux);
				}
				if (BetaDetectorCutResult_use[2] && rig_chain[c] > 150) {
					histManager->IDH4b[iz][c].get()->Fill(1.0 / beta_det[2], weight_NucFlux);
				}
		
				// [FILL] ID.H5a/H5b:  (RICH - Tracker) Δ(1/β) vs 测量刚度
				// [FILL] ID.H5a2/H5b2:  (RICH - Tracker) Δ(1/β) vs 产生刚度
				// [FILL] ID.H5a3/H5b3:  (mea beta - gene beta) Δ(1/β) vs 产生刚度
					if (richBeta > 0 && beyondRigCutoff) {
						const double dx = beta_tracker > 0 ? 1.0 / richBeta - 1.0 / beta_tracker : -9;
						const double dx_true = 1.0 / richBeta - 1.0 / generatedBeta;
							if (BetaDetectorCutResult_use[1]) {
								histManager->IDH5a[iz][c].get()->Fill(dx, rig_chain[c], weight_NucFlux);
								if(!isISS){
									histManager->IDH5a2[iz][c].get()->Fill(dx, generatedRig, weight_NucFlux);	
									histManager->IDH5a3[iz][c].get()->Fill(dx_true, generatedRig, weight_NucFlux);	
								}
							}
							if(BetaDetectorCutResult_use[2]) {
								histManager->IDH5b[iz][c].get()->Fill(dx, rig_chain[c], weight_NucFlux);
								if(!isISS){
									histManager->IDH5b2[iz][c].get()->Fill(dx, generatedRig, weight_NucFlux);	
									histManager->IDH5b3[iz][c].get()->Fill(dx_true, generatedRig, weight_NucFlux);	
								}
							}
					}
				
				// [FILL] ID.H6a/H6b: Δ(1/β) vs EK (TOF-RICH)
				// [FILL] ID.H7a/H7b (per chain): Δ(1/β) vs RICH beta*gamma
				if (BetaDetectorCutResult_use[0]) {
					if (beta_det[0] > 0 && richBeta > 0) {
						const double dx = 1.0 / beta_det[0] - 1.0 / richBeta;
							if (BetaDetectorCutResult_use[1] && getBeyondBetaCutoffCut(1, z_use, A_use)) {
								histManager->IDH6a[iz][c].get()->Fill(dx, ek_det[1], weight_NucFlux);
								histManager->IDH7a[iz][c].get()->Fill(dx, beta_det[1]*rig_chain[c], weight_NucFlux);
							} 
							if (BetaDetectorCutResult_use[2] && getBeyondBetaCutoffCut(2, z_use, A_use)) {
								histManager->IDH6b[iz][c].get()->Fill(dx, ek_det[2], weight_NucFlux);
								histManager->IDH7b[iz][c].get()->Fill(dx, beta_det[2]*rig_chain[c], weight_NucFlux);
							}
					}
				}
				
			}
		}
		// =========================
		// BKG histograms 
		// =========================

		if (isISS) {
			for (int s = 0; s < NsrcLoc; ++s) {
				const int zsrc = source_Z[s];
				const int Amin_src = getMinAForZ(zsrc);
				if (Amin_src < 0) continue;
				bool BetaDetectorCutResult_zsrc[3] = {//charge-dependent
					tof_cut.cutTOF(zsrc, isISS).total,
					rich_NaF && rich_cut.cutRICH(zsrc, isISS, true).total,
					!rich_NaF && rich_cut.cutRICH(zsrc, isISS, true).total
				};

				// Legacy decisions for BKGH1/BKGH3
				auto l1Pass = tracker_cut.BkgSourceOrFragCut(zsrc, /*isISS=*/true, fragZ, /*isL2Frag=*/false, forBackground);
				auto l2Pass = tracker_cut.BkgSourceOrFragCut(zsrc, /*isISS=*/true, fragZ, /*isL2Frag=*/true, forBackground);

				// New: 6-bit charge template decisions for BKGH2
				auto cuts6 = tracker_cut.chargeTempCut(zsrc, /*fragZ=*/fragZ, /*isISS=*/true, forBackground);

				for (int c = 0; c < std::min(2, NchainLoc); ++c) {
					// [FILL] ISS.BKG.H1 (legacy, unchanged)
					if (l1Pass[c]) {
						for (int d = 0; d < NdetLoc; ++d) {
							if (!detValidBkg[d]) continue;
							//using strictest cutoff(max Z/A) for all BKG hists
							if (!getBeyondBetaCutoffCut(d, zsrc, Amin_src)) continue; //L1source, source domi
							histManager->ISS_BKGH1[c][s][d].get()->Fill(ek_det[d], 1.0);
						}
					}

					// NEW: [FILL] ISS.BKG.H2 (charge vs Ek/n) without charge_types dependency
					for (int d = 0; d < NdetLoc; ++d) {
						// direct cutoff beta cut uses (Z=zsrc, A=Amin_src)
						if (!getBeyondBetaCutoffCut(d, zsrc, Amin_src)) continue;
						//if(beyondRigCutoff == false) continue;

						// t = 0 assumed to be L1QSignal
						if(s==2){
							bool pass = (c == 1) ? cuts6.details[0] : cuts6.details[1];
							if (pass && detValidBkg[d]) { // L1QSignal uses forBkg detector mask
								double x_charge = (c == 1) ? tk_ql1 : tk_ql1_unbiased;
								histManager->ISS_BKGH2[c][s][d][0].get()->Fill(x_charge, ek_det[d], 1.0);
							}
						}

						// t = 1 assumed to be L1QTemplate
						{
							bool pass = (c == 1) ? cuts6.details[2] : cuts6.details[3];
							if (pass && BetaDetectorCutResult_zsrc[d]) { // L1QTemplate uses FULL detector quality
								double x_charge = (c == 1) ? tk_ql1 : tk_ql1_unbiased;
								histManager->ISS_BKGH2[c][s][d][1].get()->Fill(x_charge, ek_det[d], 1.0);
							}
						}

						// t = 2 assumed to be L2QTemplate with L2Q Tuning
						{
							bool pass = (c == 1) ? cuts6.details[4] : cuts6.details[5];
							if (pass && BetaDetectorCutResult_zsrc[d]) { // L2QTemplate uses FULL detector quality

								// --- Start of L2Q Tuning ---
								/*
								// 1. Determine the energy bin for the current event using the BinningManager.
								int ekBin = Tools::findBin(StdBetaBins, beta_det[d]);
								// 2. Define the string parameters for the tuning function.
								const std::string current_chain = (c == 1) ? "L1Inner" : "UnbiasedL1Inner";
								const std::string& current_nucleus = sources[s];
								const std::string& current_detector = detectors[d];
								// 3. Call the tuning function to get the corrected charge.
								double used_ql2 = s < 2 ? tk_ql2 : Tools::tuneL2Charge(
									current_chain, current_nucleus, current_detector, ekBin, tk_ql2
								);
								// 4. Fill the histogram with the TUNED L2 charge.
								*/
								histManager->ISS_BKGH2[c][s][d][2].get()->Fill(tk_ql2, ek_det[d], 1.0);
								// --- End of L2Q Tuning ---
							}
						}
					}

					// [FILL] ISS.BKG.H3
					if (l2Pass[c]) {
						for (int d = 0; d < NdetLoc; ++d) {
							if (!detValidBkg[d]) continue;
							if (!getBeyondBetaCutoffCut(d, zsrc, Amin_src)) continue;
							histManager->ISS_BKGH3[c][s][d].get()->Fill(ek_det[d], weight_NucFlux);
						}
					}
					// [FILL] ISS.BKG.H4 (frag isotope resolved, 1/Mrec vs Ek/n)
					if (l2Pass[c]) {
						for (int d = 0; d < NdetLoc; ++d) {
							if (!BetaDetectorCutResult[d]) continue; //use full quality cut for mass rec, just get isotope fraction, so maybe is ok(contain cut influence fragment)

							// cutoff beta cut uses (Z=zsrc, A=Amin_src) -- same as H1/H3
							if (!getBeyondBetaCutoffCut(d, zsrc, Amin_src)) continue;

							// Reconstruct mass at L2 with fragment charge hypothesis (fragZ)
							// Use per-chain rigidity consistent with your H2/H3 logic
							auto mres_frag = Tools::calculateMass(beta_det[d], 1.0, rig_chain[c], RecCharge_Int);
							const double invMass_rec = mres_frag.invMass;

							// Fill each isotope bin's H4 with the same reconstructed 1/M and Ek/n
							for (size_t bi = 0; bi < FragA.size(); ++bi) {
								histManager->ISS_BKGH4[c][s][d][static_cast<int>(bi)].get()->Fill(invMass_rec, ek_det[d], weight_NucFlux);
							}
						}
					}
				}
			}
		} else {
			// MC-specific (no source dimension)

			// L1: input-mother (single source in MC)
			auto l1Pass_MC = tracker_cut.BkgSourceOrFragCut(charge, /*isISS=*/false, fragZ, /*isL2Frag=*/false, forBackground);
			// L2: fragmentation
			auto l2Pass_MC = tracker_cut.BkgSourceOrFragCut(charge, /*isISS=*/false, fragZ, /*isL2Frag=*/true, forBackground);

			//only basic cut
			auto QandL1andBkgIndepCuts = tracker_cut.Q_L1_BkgIndependCut(charge, false);

			// [FILL] MC.BKG.H1
			if (geneID_MC != -1 && mtrpar[0] == geneID_MC) {
				for (int c = 0; c < std::min(2, NchainLoc); ++c) {
					if (!l1Pass_MC[c]) continue;
					for (int d = 0; d < NdetLoc; ++d) {
						if (!detValidBkg[d]) continue;
						histManager->MC_BKGH1[c][d].get()->Fill(ek_det[d], weight_NucFlux);
					}
				}
			}
			// [FILL] MC.BKG.H1b and H2b
			for (int c = 0; c < std::min(2, NchainLoc); ++c) {
				if (!QandL1andBkgIndepCuts) continue;
				for (int d = 0; d < NdetLoc; ++d) {
					if (!detValidBkg[d]) continue;
					//H1b, aboveL2 total source mc events, only basic selection
					histManager->MC_BKGH1b[c][d].get()->Fill(ek_det[d], weight_NucFlux);
					//H2b, aboveL2 fragmentation mc events
					for (int bi = 0; bi < NisoBKG; ++bi) {
						const int targetFragID = fragIDs_global[bi];
						if (mtrpar[1] == targetFragID && tk_qin[0][2] > fragZ - 0.55 && tk_qin[0][2] < fragZ + 0.45 && tk_qrmn[0][2] < 0.55){
							histManager->MC_BKGH2b[c][d][bi].get()->Fill(ek_det[d], weight_NucFlux);
						}
					}
				}
			}

			// [FILL] MC.BKG.H2,H2c
			if (geneID_MC != -1 && mtrpar[0] == geneID_MC && NisoBKG > 0) {
				for (int c = 0; c < std::min(2, NchainLoc); ++c) {
					if (!l2Pass_MC[c]) continue;
					for (int d = 0; d < NdetLoc; ++d) {
						if (!detValidBkg[d]) continue;
						for (int bi = 0; bi < NisoBKG; ++bi) {
							const int targetFragID = fragIDs_global[bi];
							if(mtrpar[1] == targetFragID || mtrpar[2] == targetFragID || mtrpar[3] == targetFragID)
							{
								histManager->MC_BKGH2c[c][d][bi].get()->Fill(ek_det[d], weight_NucFlux);
							}
							if (mtrpar[1] == targetFragID){
								histManager->MC_BKGH2[c][d][bi].get()->Fill(ek_det[d], weight_NucFlux);
							}
						}
					}
				}
			}

			// Note: MC.BKG.H3a/H3b
			if (NisoBKG > 0) {
				auto twoAcc_frag = tracker_cut.TwoAccTrackerCut(fragZ, isISS, forBackground);
				bool TwoAccTrackerCutResult_frag[2] = {twoAcc_frag.details[0], twoAcc_frag.details[1]};

				bool BetaDetectorCutResult_frag[3] = {
					tof_cut.cutTOF(fragZ, isISS).total,
					rich_NaF && rich_cut.cutRICH(fragZ, isISS, true).total,
					!rich_NaF && rich_cut.cutRICH(fragZ, isISS, true).total
				};

				for (int c = 0; c < std::min(2, NchainLoc); ++c) {
					if (!TwoAccTrackerCutResult_frag[c]) continue;

					for (int d = 0; d < NdetLoc; ++d) {
						if (!BetaDetectorCutResult_frag[d]) continue;

						for (int bi = 0; bi < NisoBKG; ++bi) {
							const int targetFragID = fragIDs_global[bi];

							//int MCID = analyzer_->getGeneID(fragZ, FragA[bi]); 
							bool IsFragOrig = (fragZ == charge) && (FragA[bi] == UseMass);

							// H3a: upTOF fragmentation truth L2 Ek
							// H3b: upTOF fragmentation measured Ek
							if (IsFragOrig || (mtrpar[1] == targetFragID)) {
								histManager->MC_BKGH3a[c][d][bi].get()->Fill(L2TruthEk_n, weight_NucFlux);
								histManager->MC_BKGH3b[c][d][bi].get()->Fill(ek_det[d], weight_NucFlux);
							}
							if(IsFragOrig || (mtrpar[1] == targetFragID || mtrpar[2] == targetFragID || mtrpar[3] == targetFragID))
							{
								histManager->MC_BKGH3c[c][d][bi].get()->Fill(ek_det[d], weight_NucFlux);
							}
						
						}
					}
				}
			}
		}
	} // end main loop

	// [FILL] MC.FLUX.H3 (generated counts vs E_k/n)
	if (!isISS) {
		const std::string fluxName = AMS_Iso::Tools::selectFluxName(charge, UseMass);
		const auto& fmap  = AMS_Iso::Tools::getFluxMap();
		const auto& fnorm = AMS_Iso::Tools::getFluxNorm();

		auto itF = fmap.find(fluxName);
		auto itN = fnorm.find(fluxName);
		if (fluxName.empty() || itF == fmap.end() || itN == fnorm.end()) {
			std::cerr << "[ERROR] No flux TF1 or norm for (Z=" << mch << ", A=" << UseMass << ")\n";
			AMS_Iso::Tools::cleanupFluxFunctions();
			return;
		}

		TF1* f_flux = itF->second.get();
		//TF1 f_flux = AMS_Iso::Tools::f_MC;
		const double base_fluxIntegral = itN->second; // 全范围积分（init 时计算好的 fluxNorm）
		//const double base_fluxIntegral = AMS_Iso::Tools::MC_norm;	
		if (!(base_fluxIntegral > 0)) {
			std::cerr << "[ERROR] base_fluxIntegral <= 0 for " << fluxName << std::endl;
			AMS_Iso::Tools::cleanupFluxFunctions();
			return;
		}

		TH1D* h_flux = histManager->MC_FLUXH3[0].get();
		if (!h_flux) {
			std::cerr << "[ERROR] MC_FLUXH3 histogram not found.\n";
			AMS_Iso::Tools::cleanupFluxFunctions();
			return;
		}

		const TAxis* xAxis = h_flux->GetXaxis();
		if (!xAxis || !xAxis->GetXbins() || !xAxis->GetXbins()->GetArray()) {
			std::cerr << "[FATAL] MC_FLUXH3 x-axis bins invalid." << std::endl;
			AMS_Iso::Tools::cleanupFluxFunctions();
			return;
		}
		const double* ekBins = xAxis->GetXbins()->GetArray();

		const double N_gen = std::accumulate(mc_events.begin(), mc_events.end(), 0.0) - 2.0;

		double intcheck = 0;
		for (int j = 0; j < xAxis->GetNbins(); ++j) {
			const double EkLow = ekBins[j];
			const double EkUp  = ekBins[j + 1];
			//convert standard Ek/n bin to specific rigidity range corresponding to (Z,A) 
			double Rlow = AMS_Iso::Tools::kineticEnergyToRigidity(EkLow, charge, UseMass);
			double Rup  = AMS_Iso::Tools::kineticEnergyToRigidity(EkUp, charge, UseMass);


			std::cout << "Bin " << j << ": Ek [" << EkLow << ", " << EkUp << "] GeV/n <=> R [" << Rlow << ", " << Rup << "] GV" << std::endl;

			if (Rlow < AMS_Iso::Tools::geneRig_low) Rlow = AMS_Iso::Tools::geneRig_low;
			if (Rup  > AMS_Iso::Tools::geneRig_up)  Rup  = AMS_Iso::Tools::geneRig_up;

			double N_gen_bin = 0.0;
			if (Rup > Rlow) {
				const double fluxIntegral_bin = f_flux->Integral(Rlow, Rup);
				std::cout << "  Flux integral in bin = " << fluxIntegral_bin << std::endl;
				intcheck += fluxIntegral_bin;
				if (fluxIntegral_bin > 0) {
					N_gen_bin = N_gen * (fluxIntegral_bin / base_fluxIntegral);
				}
			}

			h_flux->SetBinContent(j + 1, N_gen_bin);
		}
		h_flux->SetBinContent(0, N_gen);
		std::cout<<"Total flux integral check: " << intcheck << " (base: " << base_fluxIntegral << ")\n";
	}

	AMS_Iso::Tools::cleanupFluxFunctions();
    if (filteredTree) {
        std::cout << "[selectdata] Filled " << filteredTree->GetEntries() 
                  << " events into filtered tree (out of " << nentries << " total)." << std::endl;
    }
	std::cout << "Event processing completed" << std::endl;
}
