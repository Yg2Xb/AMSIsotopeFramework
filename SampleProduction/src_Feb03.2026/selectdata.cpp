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
    double AveMass_charge = isISS ? 0 : UseMass;
    if (isISS && isotopeWeights.count(charge)) {
            double sumW = 0, sumM = 0;
            for (auto& p : isotopeWeights.at(charge)) { sumM += p.first * p.second; sumW += p.second; }
            if (sumW > 0) AveMass_charge = sumM / sumW;
    }
    double AveMass_fragZ = isISS ? 0 : UseMass;
    if (isISS && isotopeWeights.count(fragZ)) {
            double sumW = 0, sumM = 0;
            for (auto& p : isotopeWeights.at(fragZ)) { sumM += p.first * p.second; sumW += p.second; }
            if (sumW > 0) AveMass_fragZ = sumM / sumW;
    }
    
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
    //if(isISS)  AMS_Iso::Tools::initChargeTuning("/eos/ams/group/ihep/zixuan/ForSampleProduction/withBkg_CDFLookupTable_fromSpline.root");

    // --- Constants & Helpers ---
    const int NchainLoc = (int)chains.size();
    const int NdetLoc = 3; 
	const int NcutGroups = (int)cut_groups.size();
    const int NisoLoc = iso->getIsotopeCount();
    const std::vector<int> source_Z = {2,3,4,5,6,7,8};
    const int NsrcLoc = isISS ? (int)source_Z.size() : 1;
    
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
    auto StdBetaBins = binMgr.GetBetaBins(4, 7);
    auto* bins = histManager->FLUXH1[isISS?3:0][0][0][3]->GetXaxis()->GetXbins();
    std::vector<double> StdRigBins(bins->GetArray(), bins->GetArray() + bins->GetSize());

    std::map<std::pair<int,int>, int> ZAMap;
    for (int n = 0; n < Constants::N_nuc; ++n) ZAMap[{Constants::nuclei_Z[n], Constants::nuclei_A[n]}] = n;
    
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
        if (isISS && !rti_cut.cutRTI().total) continue;

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
                        for (int ip = 1; ip <= h->GetNbinsX(); ++ip){
                             if (StdBetaBins[ip-1] >= betaCut) h->AddBinContent(ip, expTime);
                        }
                    }
                }
            }
            timeTag.push_back(time[0]);
        }
        */
        

        // Beta Correction & Calibration
        //make sure must be a valid rec charge
        int Charge_forCorr = RecCharge_Int > 0 ? RecCharge_Int : charge; // Fallback logic simplified for readablity
        if(RecCharge_Int <= 0) {
            double lowQ = (tof_ql[2] + tof_ql[3]) * 0.5;
            if (tof_ql[2] * tof_ql[3] == 0) lowQ *= 2.0;
            if (lowQ > 0) Charge_forCorr = (int)(lowQ + 0.5f);
            else if (tk_qln[0][7][2] > 0) Charge_forCorr = (int)(tk_qln[0][7][2] + 0.5f);
        }
        double AveMass = isISS ? 0 : UseMass;
        if (isISS && isotopeWeights.count(Charge_forCorr)) {
             double sumW = 0, sumM = 0;
             for (auto& p : isotopeWeights.at(Charge_forCorr)) { sumM += p.first * p.second; sumW += p.second; }
             if (sumW > 0) AveMass = sumM / sumW;
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
        int L1_Z = mtrz[0] & 0x3F;
        double L1TruthRig  = double(mtrmom[0]) / L1_Z;
        double L1_A        = double(findIsotopeMass(mtrpar[0], L1_Z));
        double L1TruthBeta = L1_A > 0 ? Tools::rigidityToBeta(L1TruthRig, L1_Z, L1_A, false) : -9;
        double L1TruthEk_n = L1_A > 0 ? Tools::rigidityToKineticEnergy(L1TruthRig, L1_Z, L1_A) : -9;

        int L2_Z = mtrz[1] & 0x3F;
        double L2TruthRig  = double(mtrmom[1]) / L2_Z;
        double L2_A        = double(findIsotopeMass(mtrpar[1], L2_Z));
        double L2TruthBeta = L2_A > 0 ? Tools::rigidityToBeta(L2TruthRig, L2_Z, L2_A, false) : -9;
        double L2TruthEk_n = L2_A > 0 ? Tools::rigidityToKineticEnergy(L2TruthRig, L2_Z, L2_A) : -9;
        
        double generatedRig  = mmom / mch;
        double generatedEk_n  = Tools::rigidityToKineticEnergy(generatedRig, mch, UseMass);
        double generatedBeta = !isISS ? Tools::rigidityToBeta(generatedRig, mch, UseMass, false) : -9;
        double weight_NucFlux = isISS ? 1.0 : Tools::calculateWeight(mmom, mch, UseMass, isISS);

        double beta_det[3] = { TOFBeta, rich_NaF ? richBeta : -9, !rich_NaF ? richBeta : -9 };
        double ek_det[3]   = { Tools::betaToKineticEnergy(TOFBeta), 
                               (rich_NaF && richBeta > 0) ? Tools::betaToKineticEnergy(richBeta) : -9,
                               (!rich_NaF && richBeta > 0) ? Tools::betaToKineticEnergy(richBeta) : -9 };

        // Pre-compute Selection Results
        auto twoAcc = tracker_cut.TwoAccTrackerCut(charge, isISS, forBackground);
        bool PassTwoAcc[2] = { twoAcc.details[0], twoAcc.details[1] };
        
        auto twoAcc_frag = tracker_cut.TwoAccTrackerCut(fragZ, isISS, forBackground);
        bool PassTwoAcc_frag[2] = { twoAcc_frag.details[0], twoAcc_frag.details[1] };

        // 0:TOF, 1:NaF, 2:Agl
        bool BetaDetQual[3]      = { tof_cut.cutTOF(charge, isISS).total, rich_NaF && rich_cut.cutRICH(charge, isISS, true).total, !rich_NaF && rich_cut.cutRICH(charge, isISS, true).total };
        bool BetaDetQual_frag[3] = { tof_cut.cutTOF(fragZ, isISS).total,  rich_NaF && rich_cut.cutRICH(fragZ, isISS, true).total,  !rich_NaF && rich_cut.cutRICH(fragZ, isISS, true).total };
        bool DetValidBkg[3]      = { tof_cut.cutTOF(charge, isISS).total, rich_NaF && rich_cut.cutRICHforBkg(charge, isISS, true).total, !rich_NaF && rich_cut.cutRICHforBkg(charge, isISS, true).total };
        bool BetaDetGeo[3]      = { tof_cut.cutTOF(charge, isISS).details[0], rich_cut.cutGeometry(true, 1).total, rich_cut.cutGeometry(true, 0).total };

        // Beta Cutoff Logic
        
        bool beyondRigCutoff = (!isISS);
        if (isISS && rig_chain[0] >= 0.8 && rig_chain[0] <= 3300) {
            int bin = Tools::findBin(StdRigBins, rig_chain[0]);
            if (bin >= 0 && StdRigBins[bin] > Constants::SAFE_FACTOR_RIG * cutOffRig) beyondRigCutoff = true;
        }
        

        bool beyondBetaCutoff[3][Constants::N_nuc] = {};
        for (int d = 0; d < NdetLoc; ++d) {
            int bBin = Tools::findBin(StdBetaBins, beta_det[d]);
            for (int n = 0; n < Constants::N_nuc; ++n) {
                if (!isISS || beta_det[d] >= StdBetaBins.back()) { beyondBetaCutoff[d][n] = true; continue; }
                beyondBetaCutoff[d][n] = (bBin >= 0) && Tools::isBeyondCutoff(StdBetaBins[bBin], cutOffRig, Detector::BetaTypes[d].getSafetyFactor(), Constants::nuclei_Z[n], Constants::nuclei_A[n], !isISS);
            }
        }
        
        auto getBeyondBetaCutoffCut = [&](int det, int Z, double mass) {
            if (det < 0 || det >= NdetLoc) return false;
            bool isIntegerLike = std::abs(mass - std::round(mass)) < 1e-9;
            if (isIntegerLike) {
                int A = static_cast<int>(std::round(mass)); 
                auto it = ZAMap.find({Z, A});
                if (it != ZAMap.end()) {
                    return beyondBetaCutoff[det][it->second];
                }
            }
            int bBin = Tools::findBin(StdBetaBins, beta_det[det]);
            if (!isISS || beta_det[det] > StdBetaBins.back()) return true; 
            return (bBin >= 0) && Tools::isBeyondCutoff(StdBetaBins[bBin], cutOffRig, 
                                        Detector::BetaTypes[det].getSafetyFactor(), 
                                        Z, mass, !isISS);
        };

        //-------------------------!!!-----------------------
        bool isFilled = false;
        //-------------------------!!!-----------------------

        // ---------------------------------------------------------
        // B. Signal Histograms (ID)
        // ---------------------------------------------------------
        /*
        for (int c = 0; c < std::min(2, NchainLoc); ++c) {
            if (!PassTwoAcc[c]) continue;
            for (int d = 0; d < NdetLoc; ++d) {
                if (!BetaDetQual[d]) continue;
                auto mres = Tools::calculateMass(beta_det[d], 1.0, Tools::GetSmearRigidity(rig_chain[c], isISS, d), RecCharge_Int);
                
                for (int i = 0; i < NisoLoc; ++i) {
                    //if (!getBeyondBetaCutoffCut(d, charge, iso->getMass(i))) continue;
                    if (!getBeyondBetaCutoffCut(d, charge, AveMass_charge)) continue;//use ave mass for cutoff check now, in the future we can try to use isotope mass
                    isFilled = true;
                    
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
        int nZ_Beta = isISS ? 7 : 1;
        for (int iz = 0; iz < nZ_Beta; ++iz) {
            int z_use = isISS ? source_Z[iz] : charge;
            if (z_use != RecCharge_Int && fragZ != RecCharge_Int ) continue;
            int A_use = isISS ? AveMass : UseMass;
            if (A_use < 0) continue;

            // Recalculate cuts for specific Z if different from primary charge
            bool pAcc[2], bQual[3];
            if (z_use == charge) {
                std::copy(std::begin(PassTwoAcc), std::end(PassTwoAcc), pAcc);
                std::copy(std::begin(BetaDetQual), std::end(BetaDetQual), bQual);
            } else {
                auto acc = tracker_cut.TwoAccTrackerCut(z_use, isISS, forBackground);
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

                    if (bQual[1] && rig_chain[c] > 80)  {
                        histManager->IDH4a[iz][c]->Fill(1.0/beta_det[1], weight_NucFlux);
                        isFilled = true;
                    }
                    if (bQual[2] && rig_chain[c] > 150) {
                        histManager->IDH4b[iz][c]->Fill(1.0/beta_det[2], weight_NucFlux);
                        isFilled = true;
                    }

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
                                if(z_use > 2) isFilled = true;
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
                        if (bQual[1] && getBeyondBetaCutoffCut(1, z_use, A_use)) {
                            histManager->IDH6a[iz][c]->Fill(dx, beta_det[1], weight_NucFlux);
                            if(z_use > 2) isFilled = true;
                        }
                        if (bQual[2] && getBeyondBetaCutoffCut(2, z_use, A_use)) {
                            histManager->IDH6b[iz][c]->Fill(dx, beta_det[2], weight_NucFlux);
                            if(z_use > 2) isFilled = true;
                        }
                    }
                }

                // --- Part 2: Frag Truth Resolution (Frag Cuts) ---
                if (!isISS && mtrpar[1] == checkedFragID && PassTwoAcc_frag[c]) {
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
                histManager->IDH7[c][d]->Fill((L1TruthBeta - L2TruthBeta)/L1TruthBeta, L1TruthBeta, weight_NucFlux);
            }
            histManager->IDH7[c][3]->Fill((L1TruthBeta - L2TruthBeta)/L1TruthBeta, L1TruthBeta, weight_NucFlux);
        }
        
        // ---------------------------------------------------------
        // D. Background Histograms
        // ---------------------------------------------------------
        // all sample should be under Z/A = frag/AveMass_fragZ beta cutoff, because we are estimating the background of fragZ, we need make sure each zsrc use same cutoff 
        int startSrc = isISS ? 0 : (charge - 2);
        int endSrc   = isISS ? NsrcLoc : (charge - 1);

        for (int s = startSrc; s < endSrc; ++s) {
            int zsrc = source_Z[s];
            if (getMinAForZ(zsrc) < 0) continue;

            auto charge_cuts = tracker_cut.chargeTempCut(zsrc, fragZ, isISS, forBackground);
            
            bool BetaDetQual_zsrc[3] = {
                tof_cut.cutTOF(zsrc, isISS).total,
                rich_NaF && rich_cut.cutRICH(zsrc, isISS, true).total,
                !rich_NaF && rich_cut.cutRICH(zsrc, isISS, true).total
            };

            for (int c = 0; c < std::min(2, NchainLoc); ++c) {
                bool L1BeamL2Any = tracker_cut.FragSampleSel(zsrc, fragZ, c, 0, isISS, forBackground);
                bool L1BeamL2BeamLoose  = tracker_cut.FragSampleSel(zsrc, fragZ, c, 1, isISS, forBackground);
                bool L1BeamL2BeamFull  = tracker_cut.FragSampleSel(zsrc, fragZ, c, 2, isISS, forBackground);
                bool L1BeamL2FragLoose  = tracker_cut.FragSampleSel(zsrc, fragZ, c, 3, isISS, forBackground);
                bool L1BeamL2FragFull  = tracker_cut.FragSampleSel(zsrc, fragZ, c, 4, isISS, forBackground);
                
                for (int d = 0; d < NdetLoc; ++d) {
                    if (!BetaDetGeo[d] || !getBeyondBetaCutoffCut(d, 6, 12)) continue;

                    double ek = ek_det[d];
                    int srcIdx = isISS ? s : 0;

                    if (!isISS && mtrpar[0] == geneID_MC) {
                        if (L1BeamL2Any) histManager->BKG_H1a[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                        if (L1BeamL2BeamLoose && mtrpar[1] == geneID_MC) histManager->BKG_H1b[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                        if (L1BeamL2BeamFull && mtrpar[1] == geneID_MC) histManager->BKG_H1b2[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                        if (L1BeamL2Any && tk_qinner < charge-0.5 && mtrpar[1] != geneID_MC && L2_Z <= mch) histManager->BKG_H1c[c][srcIdx][d]->Fill(ek, weight_NucFlux);
                    }

                    if (L1BeamL2FragLoose && BetaDetQual_frag[d]) {
                        auto mres = Tools::calculateMass(beta_det[d], 1.0, Tools::GetSmearRigidity(rig_chain[c], isISS, d), fragZ);
                        histManager->BKG_H2b[c][srcIdx][d]->Fill(mres.invMass, ek, weight_NucFlux);
                        if(L1BeamL2FragFull)
                            histManager->BKG_H2b2[c][srcIdx][d]->Fill(mres.invMass, ek, weight_NucFlux);
                        isFilled = true;
                    }

                    // H4: Charge Templates (Nested switch for clarity)
                    double fillQ = (c == 0) ? tk_ql1_unb : tk_ql1;
                    for (int t = 0; t <= 5; ++t) {
                        if (!charge_cuts.details[2 * t + c]) continue;
                        if (t > 2 && !BetaDetQual_zsrc[d]) continue;

                        double q_val = fillQ;
                        if (t == 4) {
                            q_val = tk_ql2;
                            //if (isISS) q_val = Tools::tuneL2Charge((c==0?"UnbiasedL1Inner":"L1Inner"), sources[s], detectors[d], Tools::findBin(StdBetaBins, beta_det[d]), tk_ql2);
                        } else if (t == 5) {
                            q_val = tk_qinner;
                        }
                        histManager->BKG_H4[c][srcIdx][d][t]->Fill(q_val, ek, weight_NucFlux);
                        if (t == 0) {
                            if (tk_ql1_unb > 2.5 && tk_ql1_unb < 8.8) {
                                isFilled = true;
                            }
                        } else {
                            isFilled = true;
                        }
                    }
                }
            }
        }

        // MC Truth for BKG (H2a, H3),beam = charge
        if (!isISS && NisoBKG > 0) {
            for (int c = 0; c < std::min(2, NchainLoc); ++c) {
                bool L1BeamL2FragLoose  = tracker_cut.FragSampleSel(charge, fragZ, c, 3, isISS, forBackground);
                bool L1BeamL2FragFull  = tracker_cut.FragSampleSel(charge, fragZ, c, 4, isISS, forBackground);
                
                for (int d = 0; d < NdetLoc; ++d) {
                    for (int bi = 0; bi < NisoBKG; ++bi) {
                        bool IsFragOrig = (fragZ == charge) && (FragA[bi] == UseMass);
                        bool isL2FragTruth = (mtrpar[1] == fragIDs_global[bi]);

                        if (L1BeamL2FragLoose && BetaDetGeo[d] && mtrpar[0] == geneID_MC && isL2FragTruth)
                            histManager->BKG_H2a[c][d][bi]->Fill(ek_det[d], weight_NucFlux);
                        if (L1BeamL2FragFull && BetaDetGeo[d] && mtrpar[0] == geneID_MC && isL2FragTruth)
                            histManager->BKG_H2a2[c][d][bi]->Fill(ek_det[d], weight_NucFlux);
                        
                        if (PassTwoAcc_frag[c] && BetaDetQual_frag[d] && (IsFragOrig || mtrpar[0] == fragIDs_global[bi])) { //L1 frag truth
                            histManager->BKG_H3a[c][d][bi]->Fill(generatedEk_n, weight_NucFlux);
                            histManager->BKG_H3b[c][d][bi]->Fill(ek_det[d], weight_NucFlux);
                        }
                    }
                }
            }
        }
        
        // BKG H5 (Correlation)
        for (int c = 0; c < 1; ++c) {
            if (!tracker_cut.chargeTempCut(8, fragZ, isISS, forBackground).details[c]) continue;
            for (int d = 0; d < NdetLoc; ++d) {
                if (BetaDetGeo[d] && getBeyondBetaCutoffCut(d, 6, 12)){
                    double Qinner = tk_qinner >= 1.5 ? tk_qinner : tk_q[1];
                    if(Qinner > 0.5 && Qinner < 9 && tk_ql1_unb > 2.5 && tk_ql1_unb < 9){
                        histManager->BKG_H5[c][d][0]->Fill(Qinner, (c == 0 ? tk_ql1_unb : tk_ql1), ek_det[d], weight_NucFlux);
                        isFilled = true;
                    }
                }
            }
        }
        */
        // ---------------------------------------------------------
        // E. Flux Efficiency Filling
        // ---------------------------------------------------------
        
        // Pre-calc common reference variables
        double ref_rig_inner = isISS ? (0.939615 + 0.739231 * cutOffRig) : generatedRig;
        double unbiasedTOFEk = AMS_Iso::Tools::betaToKineticEnergy(betahs);

        // 1. Loop Source Elements (index s: 0=He, 1=Li, 2=Be, 3=B, 4=C, 5=N, 6=O)
        for (int s = 2; s < 3; ++s) {
            // Map index to actual charge: targetZ = 2, 3, 4, 5, 6, 7, 8, for mc only charge
            int tZ = isISS ? s + 2 : charge; 
            // 2. Loop Cut Groups
            for (int cg = 0; cg < NcutGroups; ++cg) {
                const std::string& cut = cut_groups[cg];
                bool isBkgRed = (cut == "BkgReduction"), isBetaRec = (cut == "BetaRecQuality"), isInner = (cut == "InnerTracking");
                
                // Rule: Helium (s=0) only fills BkgReduction
                //if (tZ == 2 && !isBkgRed) continue;

                // 3. Get Efficiency Result for current tZ
                auto res = [&]() -> AMS_Iso::CutResult<2> {
                    if (cut == "Trigger")       return tracker_cut.getEfficiencyTrigger(tZ, isISS, forBackground);
                    if (cut == "L1QLowLimit")   return tracker_cut.getEfficiencyL1QLowLimit(tZ, isISS, forBackground);
                    if (cut == "L1PickUp")      return tracker_cut.getEfficiencyL1PickUp(tZ, isISS, forBackground);
                    if (cut == "InnerTracking") return tracker_cut.getEfficiencyInnerTracking(tZ, isISS);
                    if (cut == "InnerTrackerQ") return tracker_cut.getEfficiencyInnerTrackerQ(tZ, isISS, forBackground);
                    if (cut == "UpperTOFQ")     return tracker_cut.getEfficiencyUTOFQ(tZ, isISS, forBackground);
                    if (cut == "BkgReduction")  return tracker_cut.getEfficiencyBkgReduction(tZ, isISS);
                    return CutResult<2>();
                }();

                // 4. Pre-calc Beta Detector Quality for current tZ
                bool BetaQ_eff[3] = { 
                    tof_cut.cutTOF(tZ, isISS).total, 
                    rich_NaF && rich_cut.cutRICH(tZ, isISS, true).total, 
                    !rich_NaF && rich_cut.cutRICH(tZ, isISS, true).total 
                };

                // 5. Detector Loop (0:TOF, 1:NaF, 2:Agl, 3:Tracker)
                for (int d = 0; d < 4; ++d) { 
                    bool isTrk = (d == 3);

                    if (isTrk == isBetaRec) continue; 

                    // Z/A Ratio Rule: use particle's charge if tZ matches event charge, else use Carbon (6/12.0)
                    int cZ = (tZ == charge) ? charge : 6; 
                    double cA = (tZ == charge) ? AveMass_charge : 12.0;
                    // Cutoff checks
                    if (isTrk) { if (!beyondRigCutoff) continue; } 
                    else { if (!getBeyondBetaCutoffCut(d, cZ, cA)) continue; }

                    // --- Determine Value to fill ---

                    double val = -3;
                    if (isInner) {
                        double unbiasedRig = AMS_Iso::Tools::kineticEnergyToRigidity(unbiasedTOFEk, cZ, cA);
                        val = (unbiasedRig > 0 && unbiasedRig <= 6.3) ? unbiasedRig : ref_rig_inner;
                    } 
                    else {
                        // Tracker fills Rigidity, Detectors fill converted Ek
                        val = isTrk ? rig_chain[0] : AMS_Iso::Tools::rigidityToKineticEnergy(rig_chain[0], cZ, cA);
                    }

                    // --- Pass/Fail Condition ---
                    bool pDen = false, pNum = false;
                    if (isBetaRec) {
                        bool trkAcc = (d == -1) ? tracker_cut.TwoAccTrackerCut_OneTrk(tZ, isISS).details[0] 
                                               : tracker_cut.TwoAccTrackerCut(tZ, isISS, false).details[0];
                        pDen = trkAcc && BetaDetGeo[d];
                        pNum = pDen && BetaQ_eff[d];
                    } else {
                        pDen = res.details[1];
                        pNum = res.details[0];
                    }

                    // 6. Fill Histograms: [source_idx][cut_idx][numden][det]
                    if (pDen && histManager->FLUXH1[s][cg][1][d]) {
                        histManager->FLUXH1[s][cg][1][d]->Fill(val, weight_NucFlux); 
                        if(tZ > 2) isFilled = true;
                    }
                    if (pNum && histManager->FLUXH1[s][cg][0][d]) {
                        histManager->FLUXH1[s][cg][0][d]->Fill(val, weight_NucFlux);
                        if(tZ > 2) isFilled = true;
                    }
                }
            }
        }
        
        /*
        //ISS FLUXH5
        if(isISS && PassTwoAcc[0]) histManager->ISS_FLUXH5[0]->Fill(cutOffRig, rig_chain[0], weight_NucFlux);
        if(isISS && PassTwoAcc[1]) histManager->ISS_FLUXH5[1]->Fill(cutOffRig, rig_chain[1], weight_NucFlux);
            

        if(isISS && rig_chain[0] > 1.2*cutOffRig) histManager->ISS_FLUXH4[0]->Fill(run, btstat);

        // ---------------------------------------------------------
        // A. Filtered Tree
        // ---------------------------------------------------------
        if (filteredTree) {
            auto l1n = tracker_cut.cutL1Norm(charge, isISS);
            auto l1u = tracker_cut.cutL1Unbiased(charge, isISS);
            bool l1n_pass = l1n.details[2] && l1n.details[3] && l1n.details[4] && tk_ql1 > 2.5 && tk_ql1 < 8.8;
            bool l1u_pass = l1u.details[2] && l1u.details[3] && tk_ql1_unb > 2.5 && tk_ql1_unb < 8.8;
            
            if (isFilled || tracker_cut.Q_L1_BkgIndependCut(charge, isISS) && (l1n_pass || l1u_pass)){
                if(tk_ql1_unb < 9) filteredTree->Fill();
            } 
        }
        */
        
    } // End Event Loop

    // --- Finalization: MC Flux ---
    if (!isISS) {
        auto* h_flux = histManager->MC_FLUXH3[0].get();
        std::string fluxName = AMS_Iso::Tools::selectFluxName(charge, UseMass);
        const auto& fmap = AMS_Iso::Tools::getFluxMap();
        
        TF1* f = nullptr;
        double norm = 0.0;

        if (fmap.count(fluxName)) {
            f = fmap.at(fluxName).get();
            norm = AMS_Iso::Tools::getFluxNorm().at(fluxName);
        } else {
            f = &AMS_Iso::Tools::f_Reweight;
            norm = AMS_Iso::Tools::Reweight_norm;
        }

        if (f && norm > 0) {
            double Ngen = std::accumulate(mc_events.begin(), mc_events.end(), 0.0) - 2.0;
            for (int j = 1; j <= h_flux->GetNbinsX(); ++j) {
                double r1 = std::max(AMS_Iso::Tools::geneRig_low, Tools::kineticEnergyToRigidity(h_flux->GetBinLowEdge(j), charge, UseMass));
                double r2 = std::min(AMS_Iso::Tools::geneRig_up,  Tools::kineticEnergyToRigidity(h_flux->GetBinLowEdge(j+1),  charge, UseMass));
                if (r2 > r1) h_flux->SetBinContent(j, Ngen * (f->Integral(r1, r2) / norm));
            }
            h_flux->SetBinContent(0, Ngen);
        }
    }
        

    AMS_Iso::Tools::cleanupFluxFunctions();
    if (filteredTree) std::cout << "[selectdata] Filtered Tree: " << filteredTree->GetEntries() << " events.\n";
    std::cout << "Event processing completed" << std::endl;
}