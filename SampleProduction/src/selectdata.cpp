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

    // --- Configuration & Initialization ---
    const bool isISS = analyzer_->isISS();
    const int charge = analyzer_->getCharge(); 
    const int UseMass = analyzer_->getUseMass();
    const int FragmentZ = analyzer_->getFragmentZ();
    const int fragZ = isISS ? charge : FragmentZ;
    const bool forBackground = analyzer_->isNoBkgCut();
    
    auto* histManager = analyzer_->getHistManager();
    if (!histManager) { std::cerr << "Failed to get HistManager\n"; return; }
    TTree* filteredTree = histManager->GetFilteredTree();

    const IsotopeVar* iso = analyzer_->getIsotope();
    if (!iso) { std::cerr << "[FATAL] IsotopeVar is null" << std::endl; return; }
    
    const std::vector<std::string> chains = analyzer_->getActiveChains();
    if (chains.empty()) { std::cerr << "No active chains found\n"; return; }

    const Long64_t nentries = fChain->GetEntries();
    std::cout << "Total entries: " << nentries << " | isISS=" << isISS << " | Z=" << charge << " | A=" << UseMass << std::endl;

    ModelManager::init("/eos/ams/group/ihep/zixuan/ForSampleProduction/model_data.root",
                       "/eos/ams/group/ihep/zixuan/ForSampleProduction/model_mc.root");
    if(!isISS) AMS_Iso::Tools::initFluxFunctions();
    if(isISS)  AMS_Iso::Tools::initChargeTuning();

    // --- Constants & Helpers ---
    const int NchainLoc = (int)chains.size();
    const int NdetLoc = 3; 
    const int NisoLoc = iso->getIsotopeCount();
    const std::vector<int> source_Z = {2,3,4,5,6,7,8};
    const int NsrcLoc = (int)source_Z.size();
    
    // MC Specifics
    const int geneID_MC = analyzer_->getGeneID(charge, UseMass);
    const std::vector<int> FragA = [&]() {
        if (isotopeWeights.find(fragZ) == isotopeWeights.end()) return std::vector<int>();
        const auto& w = isotopeWeights.at(fragZ);
        std::vector<int> r(w.size());
        std::transform(w.begin(), w.end(), r.begin(), [](auto& p){ return p.first; });
        return r;
    }();
    const int NisoBKG = (int)FragA.size();
    const std::vector<int> fragIDs_global = analyzer_->getBkgFragIDs(fragZ);
    int checkedFragID = (charge == fragZ) ? fragIDs_global.front() : fragIDs_global.back();

    if (!isISS) std::cout << "[INFO] MC: geneID=" << geneID_MC << " | fragZ=" << fragZ << " | NisoBKG=" << NisoBKG << std::endl;

    // Cutoff & Binning Helpers
    auto& binMgr = BinningManager::GetInstance();
    auto StdBetaBins = binMgr.GetBetaBins(2, 4);
    auto* bins = histManager->IDH5a[0][0]->GetYaxis()->GetXbins();
    std::vector<double> StdRigBins(bins->GetArray(), bins->GetArray() + bins->GetSize());

    bool beyondBetaCutoff[3][Constants::N_nuc] = {};
    std::map<std::pair<int,int>, int> ZAMap;
    for (int n = 0; n < Constants::N_nuc; ++n) ZAMap[{Constants::nuclei_Z[n], Constants::nuclei_A[n]}] = n;

    auto getBeyondBetaCutoffCut = [&](int det, int Z, int A) {
        if (det < 0 || det >= NdetLoc) return false;
        auto it = ZAMap.find({Z, A});
        return (it != ZAMap.end()) ? beyondBetaCutoff[det][it->second] : false;
    };
    
    auto getMinAForZ = [&](int Z) {
        int Amin = -1;
        for (int n = 0; n < Constants::N_nuc; ++n) {
            if (Constants::nuclei_Z[n] == Z && (Amin < 0 || Constants::nuclei_A[n] < Amin)) 
                Amin = Constants::nuclei_A[n];
        }
        return Amin;
    };

    // Runtime Buffers
    std::vector<unsigned int> timeTag;
    UInt_t current_run = 0, min_event = 0, max_event = 0; 
    int event_count = 1;
    std::vector<double> mc_events;

    // --- Event Loop ---
    for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
        if (LoadTree(jentry) < 0) break;
        fChain->GetEntry(jentry);
        if (jentry % 1000000 == 0) std::cout << "Processing entry " << jentry << "/" << nentries << std::endl;

        // 1. MC Accounting
        if (!isISS) {
            if (current_run != run) {
                if (current_run != 0) mc_events.push_back(max_event - min_event + 1 + (double)(max_event - min_event + 1) / event_count);
                current_run = run; min_event = event; max_event = event; event_count = 1;
            } else {
                min_event = std::min(min_event, event); max_event = std::max(max_event, event); event_count++;
                if (jentry == nentries - 1) mc_events.push_back(max_event - min_event + 1 + (double)(max_event - min_event + 1) / event_count);
            }
        }

        // 2. Variable Calculation & Cuts
        RTICut rti_cut(this);
        //if (isISS && !rti_cut.cutRTI().total) continue;

        TrackerCut tracker_cut(this);
        TOFCut tof_cut(this);
        RICHCut rich_cut(this);

        // Basic Vars
        double tk_ql1     = tracker_cut.getL1Q_Normal();
        double tk_ql1_unb = tracker_cut.getL1Q_Unbiased();
        double tk_ql2     = tracker_cut.getL2Q();
        double tk_qinner  = tracker_cut.getInnerQ();
        int RecCharge_Int = tk_qinner > 0 ? (int)(tk_qinner + 0.5f) : -1;

        double rig_chain[2] = { tracker_cut.getRigidity(), tracker_cut.getRigidity(1,2,2) };
        for(double& r : rig_chain) if(r > 0) r = 1./(1./r - 1./35000.);

        double cutOffRig = rti_cut.getCutoffRigidity();
        double TOFBeta   = tof_cut.getBeta();
        double richBeta  = yanzx_dst ? rich_cut.getBeta(1) : rich_cut.getBeta(0);

        // Exposure Calculation
        /*
        if (isISS && (std::find(timeTag.begin(), timeTag.end(), time[0]) == timeTag.end())) {
            float expTime = rti_cut.calculateExposure().value;
            // Fill ISS_FLUXH2 and ISS_FLUXH3
            if (auto* h = histManager->ISS_FLUXH2[0].get()) {
                double rigCut = Constants::SAFE_FACTOR_RIG * cutOffRig;
                for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                    if (h->GetBinLowEdge(ib) >= rigCut) h->AddBinContent(ib, expTime);
            }
            for (int d = 0; d < NdetLoc; ++d) {
                for (int i = 0; i < NisoLoc; ++i) {
                    double betaCut = Detector::BetaTypes[d].getSafetyFactor() * Tools::rigidityToBeta(cutOffRig, charge, iso->getMass(i), false);
                    if (auto* h = histManager->ISS_FLUXH3[d][i].get()) {
                        for (int ip = 1; ip <= h->GetNbinsX(); ++ip)
                             if (h->GetBinLowEdge(ip) >= betaCut) h->AddBinContent(ip, expTime);
                    }
                }
            }
            timeTag.push_back(time[0]);
        }
        */

        // Beta Correction & Calibration
        int Charge_forCorr = RecCharge_Int > 0 ? RecCharge_Int : charge; // Fallback logic simplified for readablity
        if(RecCharge_Int <= 0) {
            double lowQ = (tof_ql[2] + tof_ql[3]) * 0.5;
            if (tof_ql[2] * tof_ql[3] == 0) lowQ *= 2.0;
            if (lowQ > 0) Charge_forCorr = (int)(lowQ + 0.5f);
            else if (tk_qln[0][7][2] > 0) Charge_forCorr = (int)(tk_qln[0][7][2] + 0.5f);
        }

        double rBetaCorr = richBeta;
        if (!yanzx_dst && richBeta > 0) {
            //auto pos = rich_cut.getModifiedPosition(true);
            rBetaCorr = ModelManager::corrected_beta(richBeta, rich_NaF ? NAF : AGL, run, Charge_forCorr, 
                                                     rich_pos[0], rich_pos[1], rich_theta, rich_phi, rich_usedm, rich_hit, isISS ? 0 : 1);
            if(isISS){
                float correction=1;
                if(run>1.58013e+09){
                    float delta=(run/1.58013e+09-1);

                    if(!rich_NaF) correction+=-2.4406e-4*delta-5.9031e-3*pow(delta,2)-1.2066e-2*pow(delta,3);  // This is the correction for aerogel
                    if(rich_NaF) correction+=2.5011e-3*delta+3.7659e-2*pow(delta,2)+1.4799e-1*pow(delta,3);         // this is the correction for naf
                }
                rBetaCorr*=correction; // This is the corrected beta
            }
        }
        richBeta = isISS ? Tools::CorrectCalibrationBiasInData(rBetaCorr, rich_NaF) : Tools::GetSmearRichBeta(Charge_forCorr, rBetaCorr, rich_NaF);
        // MC Truth & Reconstructed Vectors
        int L2_Z = mtrz[1] & 0x3F;
        double L2TruthRig  = double(mtrmom[1]) / L2_Z;
        double L2_A        = double(findIsotopeMass(mtrpar[1], L2_Z));
        double L2TruthBeta = L2_A > 0 ? Tools::rigidityToBeta(L2TruthRig, L2_Z, L2_A, false) : -9;
        double L2TruthEk_n = L2_A > 0 ? Tools::rigidityToKineticEnergy(L2TruthRig, L2_Z, L2_A) : -9;
        
        double generatedRig  = mmom / mch;
        double generatedBeta = !isISS ? Tools::rigidityToBeta(generatedRig, mch, UseMass, false) : -9;
        double weight_NucFlux = isISS ? 1.0 : Tools::calculateWeight(mmom, mch, UseMass, isISS);

        double beta_det[3] = { TOFBeta, rich_NaF ? richBeta : -9, !rich_NaF ? richBeta : -9 };
        double ek_det[3]   = { Tools::betaToKineticEnergy(TOFBeta), 
                               (rich_NaF && richBeta > 0) ? Tools::betaToKineticEnergy(richBeta) : -9,
                               (!rich_NaF && richBeta > 0) ? Tools::betaToKineticEnergy(richBeta) : -9 };

        // Pre-compute Selection Results
        auto twoAcc = tracker_cut.TwoAccTrackerCut(charge, isISS, false);
        bool PassTwoAcc[2] = { twoAcc.details[0], twoAcc.details[1] };
        
        auto twoAcc_frag = tracker_cut.TwoAccTrackerCut(fragZ, isISS, false);
        bool PassTwoAcc_frag[2] = { twoAcc_frag.details[0], twoAcc_frag.details[1] };

        // 0:TOF, 1:NaF, 2:Agl
        bool BetaDetQual[3]      = { tof_cut.cutTOF(charge, isISS).total, rich_NaF && rich_cut.cutRICH(charge, isISS, true).total, !rich_NaF && rich_cut.cutRICH(charge, isISS, true).total };
        bool BetaDetQual_frag[3] = { tof_cut.cutTOF(fragZ, isISS).total,  rich_NaF && rich_cut.cutRICH(fragZ, isISS, true).total,  !rich_NaF && rich_cut.cutRICH(fragZ, isISS, true).total };
        bool DetValidBkg[3]      = { tof_cut.cutTOF(charge, isISS).total, rich_NaF && rich_cut.cutRICHforBkg(charge, isISS, true).total, !rich_NaF && rich_cut.cutRICHforBkg(charge, isISS, true).total };

        // Beta Cutoff Logic
        bool beyondRigCutoff = (!isISS);
        if (rig_chain[1] >= 0.8 && rig_chain[1] <= 3300) {
            int bin = Tools::findBin(StdRigBins, rig_chain[1]);
            if (bin >= 0 && StdRigBins[bin] > Constants::SAFE_FACTOR_RIG * cutOffRig) beyondRigCutoff = true;
        }

        for (int d = 0; d < NdetLoc; ++d) {
            int bBin = Tools::findBin(StdBetaBins, beta_det[d]);
            for (int n = 0; n < Constants::N_nuc; ++n) {
                if (!isISS || beta_det[d] >= 1) { beyondBetaCutoff[d][n] = true; continue; }
                beyondBetaCutoff[d][n] = (bBin >= 0) && Tools::isBeyondCutoff(StdBetaBins[bBin], cutOffRig, Detector::BetaTypes[d].getSafetyFactor(), Constants::nuclei_Z[n], Constants::nuclei_A[n], true);
            }
        }

        // ---------------------------------------------------------
        // A. Filtered Tree
        // ---------------------------------------------------------
        if (filteredTree) {
            auto l1n = tracker_cut.cutL1Norm(charge, isISS);
            auto l1u = tracker_cut.cutL1Unbiased(charge, isISS);
            bool l1n_pass = l1n.details[2] && l1n.details[3] && l1n.details[4] && tk_ql1 > 2.5 && tk_ql1 < 8.8;
            bool l1u_pass = l1u.details[2] && l1u.details[3] && tk_ql1_unb > 2.5 && tk_ql1_unb < 8.8;
            
            if (tracker_cut.Q_L1_BkgIndependCut(charge, isISS) && (l1n_pass || l1u_pass)) 
                filteredTree->Fill();
        }

        // ---------------------------------------------------------
        // B. Signal Histograms (ID)
        // ---------------------------------------------------------
        for (int c = 0; c < std::min(2, NchainLoc); ++c) {
            if (!PassTwoAcc[c]) continue;
            for (int d = 0; d < NdetLoc; ++d) {
                if (!BetaDetQual[d]) continue;
                auto mres = Tools::calculateMass(beta_det[d], 1.0, Tools::GetSmearRigidity(rig_chain[c], isISS, d), RecCharge_Int);
                
                for (int i = 0; i < NisoLoc; ++i) {
                    if (!getBeyondBetaCutoffCut(d, charge, iso->getMass(i))) continue;
                    
                    if (isISS) histManager->ISS_IDH1[c][d][i]->Fill(ek_det[d], weight_NucFlux);
                    histManager->IDH2[c][d][i]->Fill(mres.invMass, ek_det[d], weight_NucFlux);
                    
                    if (!isISS && geneID_MC != -1 && mtrpar[1] == geneID_MC)
                        histManager->MC_IDH1[c][d][i]->Fill(mres.invMass, ek_det[d], weight_NucFlux);
                }
            }
        }

        // ---------------------------------------------------------
        // C. Beta & Rigidity Resolution (IDH4, 5, 6, 7)
        // ---------------------------------------------------------
        double AveMass = isISS ? 0 : UseMass;
        if (isISS && isotopeWeights.count(RecCharge_Int)) {
             double sumW = 0, sumM = 0;
             for (auto& p : isotopeWeights.at(RecCharge_Int)) { sumM += p.first * p.second; sumW += p.second; }
             if (sumW > 0) AveMass = sumM / sumW;
        }

        int nZ_Beta = isISS ? 7 : 1;
        for (int iz = 0; iz < nZ_Beta; ++iz) {
            int z_use = isISS ? source_Z[iz] : charge;
            if (z_use != RecCharge_Int && fragZ != RecCharge_Int ) continue;
            int A_use = isISS ? getMinAForZ(z_use) : UseMass;
            if (A_use < 0) continue;

            // Recalculate cuts for specific Z if different from primary charge
            bool pAcc[2], bQual[3];
            if (z_use == charge) {
                std::copy(std::begin(PassTwoAcc), std::end(PassTwoAcc), pAcc);
                std::copy(std::begin(BetaDetQual), std::end(BetaDetQual), bQual);
            } else {
                auto acc = tracker_cut.TwoAccTrackerCut(z_use, isISS, false);
                pAcc[0] = acc.details[0]; pAcc[1] = acc.details[1];
                bQual[0] = tof_cut.cutTOF(z_use, isISS).total;
                bQual[1] = rich_NaF && rich_cut.cutRICH(z_use, isISS, true).total;
                bQual[2] = !rich_NaF && rich_cut.cutRICH(z_use, isISS, true).total;
            }
            for (int c = 0; c < std::min(2, NchainLoc); ++c) {
                // --- Part 1: Standard Resolution (Charge Cuts) ---
                if (pAcc[c]) {
                    double beta_tracker = (rig_chain[1] > 0 && AveMass > 0) ? 
                        Tools::rigidityToBeta(rig_chain[1], z_use, isISS ? AveMass : UseMass, false) : 0.0;

                    if (bQual[1] && rig_chain[c] > 80)  histManager->IDH4a[iz][c]->Fill(1.0/beta_det[1], weight_NucFlux);
                    if (bQual[2] && rig_chain[c] > 150) histManager->IDH4b[iz][c]->Fill(1.0/beta_det[2], weight_NucFlux);

                    if (!isISS) histManager->IDH5d3[iz][c]->Fill(1.0/rig_chain[c] - 1.0/generatedRig, generatedRig, weight_NucFlux);
                    if (!isISS && bQual[0] && beta_det[0]>0) histManager->IDH5c3[iz][c]->Fill(1.0/beta_det[0] - 1.0/generatedBeta, generatedBeta, weight_NucFlux);

                    // RICH Resolution
                    if (richBeta > 0 && beyondRigCutoff) {
                        double inv_rBeta = 1.0/richBeta;
                        double dx       = inv_rBeta - 1.0/beta_tracker;
                        double dx_truth = inv_rBeta - 1.0/generatedBeta;

                        auto fillRichRes = [&](int d_idx, auto* h1, auto* h2, auto* h3) {
                            if (bQual[d_idx]) {
                                h1->Fill(dx, rig_chain[c], weight_NucFlux);
                                if (!isISS) {   
                                    h2->Fill(dx, generatedRig, weight_NucFlux);
                                    h3->Fill(dx_truth, generatedBeta, weight_NucFlux);
                                }
                            }
                        };

                        fillRichRes(1, histManager->IDH5a[iz][c].get(),
                                    isISS ? nullptr : histManager->IDH5a2[iz][c].get(),
                                    isISS ? nullptr : histManager->IDH5a3[iz][c].get());
                        fillRichRes(2, histManager->IDH5b[iz][c].get(),
                                    isISS ? nullptr : histManager->IDH5b2[iz][c].get(),
                                    isISS ? nullptr : histManager->IDH5b3[iz][c].get());                   
                    }

                    // TOF-RICH Consistency
                    if (bQual[0] && beta_det[0] > 0 && richBeta > 0) {
                        double dx = 1.0/beta_det[0] - 1.0/richBeta;
                        if (bQual[1] && getBeyondBetaCutoffCut(1, z_use, A_use)) histManager->IDH6a[iz][c]->Fill(dx, beta_det[1], weight_NucFlux);
                        if (bQual[2] && getBeyondBetaCutoffCut(2, z_use, A_use)) histManager->IDH6b[iz][c]->Fill(dx, beta_det[2], weight_NucFlux);
                    }
                }

                // --- Part 2: Frag Truth Resolution (Frag Cuts) ---
                if (!isISS && L2_Z == fragZ && PassTwoAcc_frag[c]) {
                    histManager->IDH5d4[iz][c]->Fill(1.0/rig_chain[c] - 1.0/generatedRig, generatedRig, weight_NucFlux);
                    if (BetaDetQual_frag[0]) 
                        histManager->IDH5c4[iz][c]->Fill(1.0/beta_det[0] - 1.0/generatedBeta, generatedBeta, weight_NucFlux);
                    
                    if (richBeta > 0) {
                        double dx_truth = 1.0/richBeta - 1.0/generatedBeta;
                        if (BetaDetQual_frag[1]) histManager->IDH5a4[iz][c]->Fill(dx_truth, generatedBeta, weight_NucFlux);
                        if (BetaDetQual_frag[2]) histManager->IDH5b4[iz][c]->Fill(dx_truth, generatedBeta, weight_NucFlux);
                    }
                }
            }
        }

        // IDH7: Frag Rig Change
        for (int c = 0; c < std::min(2, NchainLoc); ++c) {
            if (isISS || !(PassTwoAcc_frag[c] && mtrpar[1] == checkedFragID)) continue;
            for (int d = 0; d < NdetLoc; ++d) {
                if (!BetaDetQual_frag[d]) continue;
                histManager->IDH7[c][d]->Fill((generatedRig - L2TruthRig)/generatedRig, generatedRig, weight_NucFlux);
            }
            histManager->IDH7[c][3]->Fill((generatedRig - L2TruthRig)/generatedRig, generatedRig, weight_NucFlux);
        }

        // ---------------------------------------------------------
        // D. Background Histograms
        // ---------------------------------------------------------
        int startSrc = isISS ? 0 : (charge - 2);
        int endSrc   = isISS ? NsrcLoc : (charge - 1);

        for (int s = startSrc; s < endSrc; ++s) {
            int zsrc = source_Z[s];
            if (getMinAForZ(zsrc) < 0) continue;

            auto charge_cuts = tracker_cut.chargeTempCut(zsrc, fragZ, isISS, forBackground);
            auto eq_num = tracker_cut.FragSampleSel(zsrc, fragZ, 0, isISS, forBackground);
            auto eq_d1  = tracker_cut.FragSampleSel(zsrc, fragZ, 1, isISS, forBackground);
            auto eq_d2  = tracker_cut.FragSampleSel(zsrc, fragZ, 2, isISS, forBackground);
            auto eq_d3  = tracker_cut.FragSampleSel(zsrc, fragZ, 3, isISS, forBackground);

            bool BetaDetQual_zsrc[3] = {
                tof_cut.cutTOF(zsrc, isISS).total,
                rich_NaF && rich_cut.cutRICH(zsrc, isISS, true).total,
                !rich_NaF && rich_cut.cutRICH(zsrc, isISS, true).total
            };

            for (int c = 0; c < std::min(2, NchainLoc); ++c) {
                for (int d = 0; d < NdetLoc; ++d) {
                    if (!DetValidBkg[d] || !getBeyondBetaCutoffCut(d, fragZ, getMinAForZ(fragZ))) continue;

                    double ek = ek_det[d];
                    int srcIdx = isISS ? s : 0;

                    // MC BKG Components
                    if (!isISS && mtrpar[0] == geneID_MC) {
                        if (eq_d1[c]) histManager->BKG_H1a[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                        if (eq_d2[c] && mtrpar[1] == geneID_MC) histManager->BKG_H1b[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                        if (eq_d3[c] && mtrpar[1] != geneID_MC && L2_Z <= mch) histManager->BKG_H1c[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                    }

                    // Data-driven BKG Est
                    if (eq_num[c] && BetaDetQual_frag[d]) {
                        auto mres = Tools::calculateMass(beta_det[d], 1.0, Tools::GetSmearRigidity(rig_chain[c], isISS, d), fragZ);
                        histManager->BKG_H2b[c][srcIdx][d]->Fill(mres.invMass, ek, weight_NucFlux);
                    }

                    // H4: Charge Templates (Nested switch for clarity)
                    double fillQ = (c == 0) ? tk_ql1_unb : tk_ql1;
                    for (int t = 0; t < 5; ++t) {
                        if (t <= 1 && zsrc != charge) continue;
                        if (!charge_cuts.details[2 * t + c]) continue;
                        if (t > 1 && !BetaDetQual_zsrc[d]) continue;

                        double q_val = fillQ;
                        if (t == 3) {
                            q_val = tk_ql2;
                            //if (isISS) q_val = Tools::tuneL2Charge((c==0?"UnbiasedL1Inner":"L1Inner"), sources[s], detectors[d], Tools::findBin(StdBetaBins, beta_det[d]), tk_ql2);
                        } else if (t == 4) {
                            q_val = tk_qinner;
                        }
                        histManager->BKG_H4[c][srcIdx][d][t]->Fill(q_val, ek, weight_NucFlux);
                    }
                }
            }
        }

        // MC Truth for BKG (H2a2, H3)
        if (!isISS && NisoBKG > 0) {
            auto numPass = tracker_cut.FragSampleSel(charge, fragZ, 1, isISS, forBackground);//use den1 cut here because no InnerQ cut
            for (int c = 0; c < std::min(2, NchainLoc); ++c) {
                for (int d = 0; d < NdetLoc; ++d) {
                    for (int bi = 0; bi < NisoBKG; ++bi) {
                        bool IsFragOrig = (fragZ == charge) && (FragA[bi] == UseMass);
                        bool isL2FragTruth = (mtrpar[1] == fragIDs_global[bi]);

                        if (numPass[c] && DetValidBkg[d] && mtrpar[0] == geneID_MC && isL2FragTruth)
                            histManager->BKG_H2a2[c][d][bi]->Fill(ek_det[d], weight_NucFlux);
                        
                        if (PassTwoAcc_frag[c] && BetaDetQual_frag[d] && (IsFragOrig || isL2FragTruth)) {
                            histManager->BKG_H3a[c][d][bi]->Fill(L2TruthEk_n, weight_NucFlux);
                            histManager->BKG_H3b[c][d][bi]->Fill(ek_det[d], weight_NucFlux);
                        }
                    }
                }
            }
        }
        
        // BKG H5 (Correlation)
        for (int c = 0; c < std::min(2, NchainLoc); ++c) {
            if (!tracker_cut.chargeTempCut(8, fragZ, isISS, forBackground).details[c]) continue;
            for (int d = 0; d < NdetLoc; ++d) {
                if (DetValidBkg[d] && getBeyondBetaCutoffCut(d, 2, 4))
                    histManager->BKG_H5[c][d][0]->Fill(tk_qinner, (c == 0 ? tk_ql1_unb : tk_ql1), ek_det[d], weight_NucFlux);
            }
        }

        // ---------------------------------------------------------
        // E. Flux Efficiency
        // ---------------------------------------------------------
        for (int iz = 0; iz < (isISS ? 3 : 1); ++iz) {
            int useZ = isISS ? (iz == 0 ? 2 : (iz == 1 ? 6 : 8)) : charge;
            bool BetaDetQual_Eff[3] = { tof_cut.cutTOF(useZ, isISS).total, rich_NaF && rich_cut.cutRICH(useZ, isISS, true).total, !rich_NaF && rich_cut.cutRICH(useZ, isISS, true).total };
            
            auto eff_num = tracker_cut.TwoAccTrackerCut(useZ, isISS, false, 0.2);
            auto eff_den = tracker_cut.TwoAccTrackerCut(useZ, isISS, true, 0.2);

            for (int d = 0; d < NdetLoc; ++d) {
                if (isISS && !getBeyondBetaCutoffCut(d, 6, 12)) continue;
                if (!BetaDetQual_Eff[d]) continue;
                for (int c = 0; c < std::min(2, NchainLoc); ++c) {
                    if (eff_den.details[c]) histManager->FLUXH1[c][0][1][d][iz]->Fill(ek_det[d], weight_NucFlux);
                    if (eff_num.details[c]) histManager->FLUXH1[c][0][0][d][iz]->Fill(ek_det[d], weight_NucFlux);
                }
            }
        }
            

        if(isISS && rig_chain[0] > 1.2*cutOffRig) histManager->ISS_FLUXH4[0]->Fill(run, btstat);

    } // End Event Loop

    // --- Finalization: MC Flux ---
    if (!isISS) {
        auto* h_flux = histManager->MC_FLUXH3[0].get();
        std::string fluxName = AMS_Iso::Tools::selectFluxName(charge, UseMass);
        const auto& fmap = AMS_Iso::Tools::getFluxMap();
        
        if (fmap.count(fluxName)) {
            double norm = AMS_Iso::Tools::getFluxNorm().at(fluxName);
            double Ngen = std::accumulate(mc_events.begin(), mc_events.end(), 0.0) - 2.0;
            TF1* f = fmap.at(fluxName).get();
            
            for (int j = 1; j <= h_flux->GetNbinsX(); ++j) {
                double r1 = std::max(AMS_Iso::Tools::geneRig_low, Tools::kineticEnergyToRigidity(h_flux->GetBinLowEdge(j), charge, UseMass));
                double r2 = std::min(AMS_Iso::Tools::geneRig_up,  Tools::kineticEnergyToRigidity(h_flux->GetBinLowEdge(j+1),  charge, UseMass));
                if (r2 > r1) h_flux->SetBinContent(j, Ngen * (f->Integral(r1, r2) / norm));
            }
            h_flux->SetBinContent(0, Ngen);
        } else {
            std::cerr << "[ERROR] No flux found for (Z=" << charge << ", A=" << UseMass << ")\n";
        }
    }
        

    AMS_Iso::Tools::cleanupFluxFunctions();
    if (filteredTree) std::cout << "[selectdata] Filtered Tree: " << filteredTree->GetEntries() << " events.\n";
    std::cout << "Event processing completed" << std::endl;
}