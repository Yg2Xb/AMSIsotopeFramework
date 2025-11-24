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
    if (!fChain || !analyzer_) { std::cerr << "[FATAL] fChain or analyzer_ is null" << std::endl; return; }

    const bool isISS = analyzer_->isISS();
    const int charge = analyzer_->getCharge(); // Primary X
    const int UseMass = analyzer_->getUseMass();
    auto* histManager = analyzer_->getHistManager();
    bool forBackground = analyzer_->isNoBkgCut();
    int FragmentZ = analyzer_->getFragmentZ(); // Target Y
    
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

    // Init Models
    ModelManager::init("/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/model_data.root",
            "/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/model_mc.root");
    if(!isISS) AMS_Iso::Tools::initFluxFunctions();
    if(isISS)  AMS_Iso::Tools::initChargeTuning();

    std::vector<unsigned int> timeTag;
    UInt_t current_run = 0, min_event = 0, max_event = 0; int event_count = 1;
    std::vector<double> mc_events;

    const int NchainLoc = static_cast<int>(chains.size());
    const int NdetLoc   = 3; 
    const int NisoLoc   = iso->getIsotopeCount();
    const std::vector<int> source_Z = {2,3,4,5,6,7,8};
    const int NsrcLoc = static_cast<int>(source_Z.size());
	const std::vector<float> minRig_NaF = {120,70,80,100,70};
	const std::vector<float> minRig_AGL = {120,70,80,100,70};

    auto getMinAForZ = [&](int Z)->int {
        int Amin = -1;
        for (int n = 0; n < Constants::N_nuc; ++n) {
            if (Constants::nuclei_Z[n] != Z) continue;
            if (Amin < 0 || Constants::nuclei_A[n] < Amin) Amin = Constants::nuclei_A[n];
        }
        return Amin;
    };

    bool beyondBetaCutoff[3][Constants::N_nuc] = {};
    std::map<std::pair<int,int>, int> ZAMap;
    for (int n = 0; n < Constants::N_nuc; ++n) ZAMap[{Constants::nuclei_Z[n], Constants::nuclei_A[n]}] = n;

    auto getBeyondBetaCutoffCut = [&](int det, int Z, int A){
        if (det < 0 || det >= NdetLoc) return false;
        auto it = ZAMap.find({Z, A});
        return (it != ZAMap.end()) ? beyondBetaCutoff[det][it->second] : false;
    };

    // MC Config
    const int geneID_MC = analyzer_->getGeneID(charge, UseMass);
    const int fragZ = isISS ? charge : FragmentZ;
    
    const std::vector<int> FragA = [&]() {
        const auto& weights = isotopeWeights.at(fragZ);
        std::vector<int> result(weights.size());
        std::transform(weights.begin(), weights.end(), result.begin(), [](const auto& p) { return p.first; });
        return result;
    }();
    const int NisoBKG = static_cast<int>(FragA.size());
    
    const std::vector<int> fragIDs_global = analyzer_->getBkgFragIDs(fragZ);
	
	if (!isISS) {
		std::cout << "[INFO] MC mode: geneID_MC=" << geneID_MC << " | fragZ=" << fragZ << " | NisoBKG=" << NisoBKG << std::endl;
	}

    // Buffers
    double weight_NucFlux = 1.;
    double cutOffRig = -1, TOFBeta = -1, richBeta = -1;
    int RecCharge_Int = -1;
    double rig_chain[2] = {-1, -1}; 
    double beta_det[3] = {-1, -1, -1};
    double ek_det[3] = {-1, -1, -1};

    // --- Event Loop ---
    for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        fChain->GetEntry(jentry);

        if (jentry % 1000000 == 0) std::cout << "Processing entry " << jentry << "/" << nentries << std::endl;

        // MC Accounting
        if (!isISS) {
            if (current_run != run) {
                if (current_run != 0) mc_events.push_back(max_event - min_event + 1 + (max_event - min_event + 1) / event_count);
                current_run = run; min_event = event; max_event = event; event_count = 1;
            } else {
                min_event = std::min(min_event, event); max_event = std::max(max_event, event); event_count++;
                if (jentry == nentries - 1) mc_events.push_back(max_event - min_event + 1 + (max_event - min_event + 1) / event_count);
            }
        }

        RTICut rti_cut(this);
        TrackerCut tracker_cut(this);
        TOFCut tof_cut(this);
        RICHCut rich_cut(this);

        if (isISS && !rti_cut.cutRTI().total) continue;

        // Vars from Cache
        double tk_ql1 = tracker_cut.getL1Q_Normal();
        double tk_ql1_unb = tracker_cut.getL1Q_Unbiased();
        double tk_ql2 = tracker_cut.getL2Q(); 
        double tk_qinner = tracker_cut.getInnerQ();
        RecCharge_Int = tk_qinner > 0 ? static_cast<int>(tk_qinner + 0.5f) : -1;
        
        rig_chain[0] = tracker_cut.getRigidity();      
        rig_chain[1] = tracker_cut.getRigidity(1,2,2); 
        if(rig_chain[0] > 0) rig_chain[0] = 1./(1./rig_chain[0] - 1./35000.);
        if(rig_chain[1] > 0) rig_chain[1] = 1./(1./rig_chain[1] - 1./35000.);

        cutOffRig = rti_cut.getCutoffRigidity();
        richBeta = yanzx_dst ? rich_cut.getBeta(1) : rich_cut.getBeta(0);
        TOFBeta = tof_cut.getBeta();

        // Exposure
        if (isISS && (std::find(timeTag.begin(), timeTag.end(), time[0]) == timeTag.end())) {
            float exposureTime = rti_cut.calculateExposure().value;
            if (auto* h = histManager->ISS_FLUXH2[0].get()) {
                double rigCut = Constants::SAFE_FACTOR_RIG * cutOffRig;
                for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
                    if (h->GetBinLowEdge(ib) >= rigCut) h->SetBinContent(ib, h->GetBinContent(ib) + exposureTime);
                }
            }
            for (int d = 0; d < NdetLoc; ++d) {
                for (int i = 0; i < NisoLoc; ++i) {
                    double betaCut = Detector::BetaTypes[d].getSafetyFactor() * Tools::rigidityToBeta(cutOffRig, charge, iso->getMass(i), false);
                    if (auto* h = histManager->ISS_FLUXH3[d][i].get()) {
                        for (int ip = 1; ip <= h->GetNbinsX(); ++ip) {
                            if (h->GetBinLowEdge(ip) >= betaCut) h->SetBinContent(ip, h->GetBinContent(ip) + exposureTime);
                        }
                    }
                }
            }
            timeTag.push_back(time[0]);
        }

        // Rich Beta Corr
        int Charge_forCorr = -1;
        if (RecCharge_Int > 0) Charge_forCorr = RecCharge_Int;
        else {
            double lowQ = (tof_ql[2] + tof_ql[3]) * 0.5;
            if (tof_ql[2] * tof_ql[3] == 0) lowQ *= 2.0;
            if (lowQ > 0) Charge_forCorr = static_cast<int>(lowQ + 0.5f);
            else if (tk_qln[0][7][2] > 0) Charge_forCorr = static_cast<int>(tk_qln[0][7][2] + 0.5f);
            else Charge_forCorr = charge;
        }

        double rBetaCorr = richBeta;
        if (!yanzx_dst && richBeta > 0) {
            auto pos = rich_cut.getModifiedPosition(true);
            Rad rad = rich_NaF ? NAF : AGL;
            rBetaCorr = ModelManager::corrected_beta(richBeta, rad, run, Charge_forCorr, pos[0], pos[1], rich_theta, rich_phi, rich_usedm, rich_hit, isISS ? 0 : 1);
        }
        richBeta = isISS ? Tools::CorrectCalibrationBiasInData(rBetaCorr, rich_NaF) : Tools::GetSmearRichBeta(Charge_forCorr, rBetaCorr, rich_NaF);

        double L2TruthRig = double(mtrmom[1]) / (mtrz[1] & 0x3F);
        double L2TruthEk_n = Tools::rigidityToKineticEnergy(L2TruthRig, charge, UseMass);
		double generatedRig = mmom / mch;
		double generatedEk_n = Tools::rigidityToKineticEnergy(generatedRig, charge, UseMass);
		double generatedBeta = Tools::rigidityToBeta(generatedRig, charge, UseMass, false);

        beta_det[0] = TOFBeta;
        beta_det[1] = rich_NaF ? richBeta : -9;
        beta_det[2] = !rich_NaF ? richBeta : -9;
        ek_det[0] = Tools::betaToKineticEnergy(TOFBeta);
        ek_det[1] = (rich_NaF && richBeta > 0) ? Tools::betaToKineticEnergy(richBeta) : -9;
        ek_det[2] = (!rich_NaF && richBeta > 0) ? Tools::betaToKineticEnergy(richBeta) : -9;

        weight_NucFlux = isISS ? 1.0 : Tools::calculateWeight(mmom, mch, UseMass, isISS);

        // Base Status (ID/BKG/Truth use strict quality if specified)
        bool BetaDetQual[3] = { tof_cut.cutTOF(charge, isISS).total, rich_NaF && rich_cut.cutRICH(charge, isISS, true).total, !rich_NaF && rich_cut.cutRICH(charge, isISS, true).total };
        bool BetaDetQual_frag[3] = { tof_cut.cutTOF(fragZ, isISS).total, rich_NaF && rich_cut.cutRICH(fragZ, isISS, true).total, !rich_NaF && rich_cut.cutRICH(fragZ, isISS, true).total };
        // For BKG Reco (H1/H2/H4) use forBkg validity
        bool DetValidBkg[3] = { tof_cut.cutTOF(charge, isISS).total, rich_NaF && rich_cut.cutRICHforBkg(charge, isISS, true).total, !rich_NaF && rich_cut.cutRICHforBkg(charge, isISS, true).total };

        bool beyondRigCutoff = false;
        if (rig_chain[1] >= 0.8 && rig_chain[1] <= 3300) {
            int bin = Tools::findBin(StdRigBins, rig_chain[1]);
            if (bin >= 0) beyondRigCutoff = isISS ? (StdRigBins[bin] > Constants::SAFE_FACTOR_RIG * cutOffRig) : true;
        }
        for (int d = 0; d < NdetLoc; ++d) {
            int bBin = Tools::findBin(StdBetaBins, beta_det[d]);
            for (int n = 0; n < Constants::N_nuc; ++n) {
                if (!isISS || beta_det[d] >= 1) { beyondBetaCutoff[d][n] = true; continue; }
                beyondBetaCutoff[d][n] = (bBin >= 0) ? Tools::isBeyondCutoff(StdBetaBins[bBin], cutOffRig, Detector::BetaTypes[d].getSafetyFactor(), Constants::nuclei_Z[n], Constants::nuclei_A[n], !isISS) : false;
            }
        }

        // ------------------- Filter Tree --------------------
        if (filteredTree) {
            auto l1n = tracker_cut.cutL1Norm(charge, isISS);
            auto l1u = tracker_cut.cutL1Unbiased(charge, isISS);
            if (tracker_cut.Q_L1_BkgIndependCut(charge, isISS) && 
               ((l1n.details[2] && l1n.details[3] && l1n.details[4] && tk_ql1 > 2.5 && tk_ql1 < 8.8) || 
                (l1u.details[2] && l1u.details[3] && tk_ql1_unb > 2.5 && tk_ql1_unb < 8.8))) {
                filteredTree->Fill();
            }
        }

        auto twoAcc = tracker_cut.TwoAccTrackerCut(charge, isISS, false);
        bool PassTwoAcc[2] = { twoAcc.details[0], twoAcc.details[1] }; // 0:Unb, 1:Norm
        
        auto twoAcc_frag = tracker_cut.TwoAccTrackerCut(fragZ, isISS, false);
        bool PassTwoAcc_frag[2] = { twoAcc_frag.details[0], twoAcc_frag.details[1] }; // 0:Unb, 1:Norm

        // ========================= ID Histograms =========================
        for (int c = 0; c < std::min(2, NchainLoc); ++c) { 
            if (!PassTwoAcc[c]) continue;
            for (int d = 0; d < NdetLoc; ++d) {
                if (!BetaDetQual[d]) continue;
                auto mres = Tools::calculateMass(beta_det[d], 1.0, rig_chain[c], RecCharge_Int);
                for (int i = 0; i < NisoLoc; ++i) {
                    int A = iso->getMass(i);
                    if (!getBeyondBetaCutoffCut(d, charge, A)) continue;
                    if (isISS) histManager->ISS_IDH1[c][d][i]->Fill(ek_det[d], weight_NucFlux);
                    histManager->IDH2[c][d][i]->Fill(mres.invMass, ek_det[d], weight_NucFlux);
                    if (!isISS && geneID_MC != -1 && mtrpar[1] == geneID_MC)
                        histManager->MC_IDH1[c][d][i]->Fill(mres.invMass, ek_det[d], weight_NucFlux);
                }
            }
        }

        // ========================= Beta Study =========================
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

        int nZ_Beta = isISS ? 7 : 1; 
        for (int iz = 0; iz < nZ_Beta; ++iz) {
            int z_use = isISS ? source_Z[iz] : charge;
			if(z_use != RecCharge_Int) continue;
            
			int A_use = isISS ? getMinAForZ(z_use) : UseMass;
            if (A_use < 0) continue;

            auto acc = tracker_cut.TwoAccTrackerCut(z_use, isISS, false);
            bool pAcc[2] = {acc.details[0], acc.details[1]};
            bool bQual[3] = { tof_cut.cutTOF(z_use, isISS).total, rich_NaF && rich_cut.cutRICH(z_use, isISS, true).total, !rich_NaF && rich_cut.cutRICH(z_use, isISS, true).total };

            for (int c = 0; c < std::min(2, NchainLoc); ++c) {
                if (!pAcc[c]) continue;

				double beta_tracker = 0.0;
				if (rig_chain[1] > 0 && AveMass > 0) {
					beta_tracker = isISS ? Tools::rigidityToBeta(rig_chain[1], z_use, AveMass, false) : Tools::rigidityToBeta(rig_chain[1], z_use, UseMass, false);
				}

                if (bQual[1] && rig_chain[c] > 80)  histManager->IDH4a[iz][c]->Fill(1.0 / beta_det[1], weight_NucFlux);
                if (bQual[2] && rig_chain[c] > 150) histManager->IDH4b[iz][c]->Fill(1.0 / beta_det[2], weight_NucFlux);

                if (richBeta > 0 && beyondRigCutoff) {
                    double dx = 1.0/richBeta - 1.0/beta_tracker;
                    double dx_truth = 1.0/richBeta - 1.0/generatedBeta;
                    if (bQual[1]) {
                        histManager->IDH5a[iz][c]->Fill(dx, rig_chain[c], weight_NucFlux);
                        if(!isISS) { 
							histManager->IDH5a2[iz][c]->Fill(dx, generatedRig, weight_NucFlux); 
							histManager->IDH5a3[iz][c]->Fill(dx_truth, generatedRig, weight_NucFlux); 
						}
                    }
                    if (bQual[2]) {
                        histManager->IDH5b[iz][c]->Fill(dx, rig_chain[c], weight_NucFlux);
                        if(!isISS) { 
							histManager->IDH5b2[iz][c]->Fill(dx, generatedRig, weight_NucFlux); 
							histManager->IDH5b3[iz][c]->Fill(dx_truth, generatedRig, weight_NucFlux); 
						}
                    }
                }
                if (bQual[0] && beta_det[0]>0 && richBeta>0) {
                    double dx = 1.0/beta_det[0] - 1.0/richBeta;
                    if (bQual[1] && getBeyondBetaCutoffCut(1, z_use, A_use)) histManager->IDH6a[iz][c]->Fill(dx, beta_det[1], weight_NucFlux);
                    if (bQual[2] && getBeyondBetaCutoffCut(2, z_use, A_use)) histManager->IDH6b[iz][c]->Fill(dx, beta_det[2], weight_NucFlux);
                }
            }
        }

        // ========================= BKG =========================

        // 1. Fragmentation Samples
        int startSrc = isISS ? 0 : (charge-2);
        int endSrc   = isISS ? NsrcLoc : (charge-2+1);

        for (int s = startSrc; s < endSrc; ++s) {
            int zsrc = source_Z[s];
            int Amin = getMinAForZ(zsrc);
            if (Amin < 0) continue;
            
            //don't use mc info here, consistent with data 
            auto charge_cuts = tracker_cut.chargeTempCut(zsrc, fragZ, isISS, forBackground);
            auto eq_num = tracker_cut.FragSampleSel(zsrc, fragZ, 0, isISS, forBackground); 
            auto eq_d1  = tracker_cut.FragSampleSel(zsrc, fragZ, 1, isISS, forBackground); 
            auto eq_d2  = tracker_cut.FragSampleSel(zsrc, fragZ, 2, isISS, forBackground); 
            auto eq_d3  = tracker_cut.FragSampleSel(zsrc, fragZ, 3, isISS, forBackground); 
            
            //DetValidBkg is charge independent actually, but BetaDetQual is 
            bool BetaDetQual_zsrc[3] = {
				tof_cut.cutTOF(zsrc, isISS).total,
				rich_NaF && rich_cut.cutRICH(zsrc, isISS, true).total,
				!rich_NaF && rich_cut.cutRICH(zsrc, isISS, true).total
			};

            for (int c = 0; c < std::min(2, NchainLoc); ++c) { 
                for (int d = 0; d < NdetLoc; ++d) {
                    if (!DetValidBkg[d] || !getBeyondBetaCutoffCut(d, fragZ, getMinAForZ(fragZ))) continue; // using min A for fragZ

                    double ek = ek_det[d];
                    int srcIdx = isISS ? s : 0;
                                
                    if (isISS || mtrpar[0] == geneID_MC){//for mc, use mc info to count den sample, den without mc info is L1signal Sample
                        if (eq_d1[c]) histManager->BKG_H1a[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                        if (eq_d2[c]) histManager->BKG_H1b[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                        if (eq_d3[c]) histManager->BKG_H1c[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                    }

                    if (eq_num[c]) {
                        histManager->BKG_H2a[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                        if(BetaDetQual_frag[d]){//this 1/mass for fit,so must use good beta, need all beta quality cut for Y sel.
                            auto mres = Tools::calculateMass(beta_det[d], 1.0, rig_chain[c], fragZ);
                            histManager->BKG_H2b[c][srcIdx][d]->Fill(mres.invMass, ek, weight_NucFlux);
                        }
                    }

                    // [FILL] H4 Charge Templates
                    // Indices: 0-2(3type L1Sig), 3(L1Temp), 4(L2Temp), 5(InnerSig), 6(InnerTemp)
                    double fillQ = (c==0) ? tk_ql1_unb : tk_ql1; 
                    for (int t = 0; t < 7; ++t) {
                        // c=0(Unb)->idx=2t, c=1(Norm)->idx=2t+1
                        if (!charge_cuts.details[2 * t + c]) continue;
                        // Templates (3, 4, 6) require Strict Beta Quality, Signals (0, 1, 2, 5) only require Basic Bkg Validity (Already checked by DetValidBkg[d] outside)
                        bool isTemplate = (t == 3 || t == 4 || t == 6);
                        if (isTemplate && !BetaDetQual_zsrc[d]) continue;
                        double q_val = 0;
                        if (t <= 3) { q_val = fillQ;} 
                        else if (t == 4) { 
                            q_val = tk_ql2; 
                            if (isISS && s >= 1) { // Tune L2
                                std::string chStr = (c==0) ? "UnbiasedL1Inner" : "L1Inner";
                                int ebin = Tools::findBin(StdBetaBins, beta_det[d]);
                                q_val = Tools::tuneL2Charge(chStr, sources[s], detectors[d], ebin, tk_ql2);
                            }
                        } 
                        else {q_val = tracker_cut.getInnerQ();}
                        histManager->BKG_H4[c][srcIdx][d][t]->Fill(q_val, ek, weight_NucFlux);
                    }
                }
            }
        } 

        // 2. MC Truth (H2a2, H3)
        if (!isISS && NisoBKG > 0) {

            auto numPass = tracker_cut.FragSampleSel(charge, fragZ, 0, isISS, forBackground);
            
            for (int c = 0; c < std::min(2, NchainLoc); ++c) {
                for (int d = 0; d < NdetLoc; ++d) {
                    for (int bi = 0; bi < NisoBKG; ++bi) {
                        int fid = fragIDs_global[bi];
                        bool IsFragOrig = (fragZ == charge) && (FragA[bi] == UseMass);
                        bool isL1SourceTruth = (mtrpar[0] == geneID_MC); 
                        bool isL2FragTruth = (mtrpar[1] == fid); 
                        
                        //same num cut, and bkgbeta cut, no full beta cuts
                        if (numPass[c] && DetValidBkg[d] && isL1SourceTruth && isL2FragTruth) {//sel num sample by cut and mc info,L1X->L2Yiso
                            histManager->BKG_H2a2[c][d][bi]->Fill(ek_det[d], weight_NucFlux);
                        }
                        // H3 Bkg Est, full cut, as same as IDH1 ISS Events sel.
                        if (PassTwoAcc_frag[c] && BetaDetQual_frag[d] && (IsFragOrig || isL2FragTruth)) {//fragZ cut and mc info
                            histManager->BKG_H3a[c][d][bi]->Fill(L2TruthEk_n, weight_NucFlux);
                            histManager->BKG_H3b[c][d][bi]->Fill(ek_det[d], weight_NucFlux);
                        }
                    }
                }
            }
        }

        // ========================= Flux(eff) =========================
        // FLUX_H1 (Efficiency) - Optimized Loop
        auto eff_num = tracker_cut.TwoAccTrackerCut(charge, isISS, false);
        auto eff_den = tracker_cut.TwoAccTrackerCut(charge, isISS, true);
        bool pass_d = eff_den.details[1];
        bool pass_n = eff_num.details[1];

        for (int d = 0; d < NdetLoc; ++d) {
            // Full Beta Cut
            if (!BetaDetQual[d]) continue;
            for (int c = 0; c < std::min(2, NchainLoc); ++c) {
                for(int i=0; i<NisoLoc; ++i) {
                    if(isISS && !getBeyondBetaCutoffCut(d, charge, iso->getMass(i))) continue;
                    if (pass_d) histManager->FLUXH1[c][0][1][d][i]->Fill(ek_det[d], weight_NucFlux);
                    if (pass_n) histManager->FLUXH1[c][0][0][d][i]->Fill(ek_det[d], weight_NucFlux);
                }
            }
        }

    } // End Event Loop

    // MC Flux H3
    if (!isISS) {
        TH1F* h_flux = histManager->MC_FLUXH3[0].get();
        const std::string fluxName = AMS_Iso::Tools::selectFluxName(charge, UseMass);
        const auto& fmap = AMS_Iso::Tools::getFluxMap();
        const auto& fnorm = AMS_Iso::Tools::getFluxNorm();
        if (fmap.count(fluxName) && fnorm.count(fluxName)) {
            TF1* f = fmap.at(fluxName).get();
            double norm = fnorm.at(fluxName);
            double Ngen = std::accumulate(mc_events.begin(), mc_events.end(), 0.0) - 2.0;
            const TAxis* xAx = h_flux->GetXaxis();
            for (int j = 1; j <= xAx->GetNbins(); ++j) {
                double e1 = xAx->GetBinLowEdge(j);
                double e2 = xAx->GetBinUpEdge(j);
                double r1 = Tools::kineticEnergyToRigidity(e1, charge, UseMass);
                double r2 = Tools::kineticEnergyToRigidity(e2, charge, UseMass);
                if (r1 < AMS_Iso::Tools::geneRig_low) r1 = AMS_Iso::Tools::geneRig_low;
                if (r2 > AMS_Iso::Tools::geneRig_up)  r2 = AMS_Iso::Tools::geneRig_up;
                if (r2 > r1) h_flux->SetBinContent(j, Ngen * (f->Integral(r1, r2) / norm));
            }
            h_flux->SetBinContent(0, Ngen);
        }
        else{
            std::cerr << "[ERROR] No flux TF1 or norm for (Z=" << mch << ", A=" << UseMass << ")\n";
			AMS_Iso::Tools::cleanupFluxFunctions();
        }
    }

    AMS_Iso::Tools::cleanupFluxFunctions();
    if (filteredTree) {
        std::cout << "[selectdata] Filled " << filteredTree->GetEntries() 
                  << " events into filtered tree (out of " << nentries << " total)." << std::endl;
    }
    std::cout << "Event processing completed" << std::endl;
}