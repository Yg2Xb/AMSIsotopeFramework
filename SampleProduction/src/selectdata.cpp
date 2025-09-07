#include "selectdata.h"
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

using namespace AMS_Iso;

static inline void EnsureModelLoaded() {
    static bool ok=false;
    if(!ok){
        std::cout<<"Loading models... "<<std::endl;
        ModelManager::init("/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/model_data.root",
                           "/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/model_mc.root");
        std::cout<<"model n:"<<ModelManager::model[0][0].index_correction.GetEntries()<<std::endl;
        ok=true;
    }
}

void selectdata::SetAnalyzer(IsotopeAnalyzer* a){ analyzer_=a; }

void selectdata::Loop() {
    if (!fChain || !analyzer_) return;

    const bool isISS = analyzer_->isISS();
    const int charge = analyzer_->getCharge();
    const int UseMass = analyzer_->getUseMass();
    auto* histManager = analyzer_->getHistManager();
    if(!histManager){ std::cerr<<"Failed to get HistManager\n"; return; }

    auto& binMgr = BinningManager::GetInstance();
    const IsotopeVar* iso = analyzer_->getIsotope();
    std::vector<std::string> chains = analyzer_->getActiveChains();
    if(chains.empty()){ std::cerr<<"No active chains found\n"; return; }

    const Long64_t nentries = fChain->GetEntries();
    std::cout<<"Total entries: "<<nentries<<std::endl;

    EnsureModelLoaded();

    std::vector<unsigned int> timeTag;
    UInt_t current_run=0, min_event=0, max_event=0; int event_count=1;
    std::vector<double> mc_events;

    for (Long64_t jentry=0;jentry<nentries;++jentry){
        Long64_t ientry = LoadTree(jentry); if (ientry<0) break;
        fChain->GetEntry(jentry);

        // --- MC run-event 统计 ---
        if(!isISS){
            if (current_run != run){
                if (current_run != 0){
                    mc_events.push_back(max_event - min_event + 1 + (max_event - min_event + 1)/event_count);
                }
                current_run = run; min_event = event; max_event = event; event_count = 1;
            } else {
                min_event = std::min(min_event, event);
                max_event = std::max(max_event, event);
                event_count++;
                if (jentry == nentries-1){
                    mc_events.push_back(max_event - min_event + 1 + (max_event - min_event + 1)/event_count);
                }
            }
        }

        RTICut rti(this);
        TrackerCut trk(this);
        TOFCut tof(this);
        RICHCut ric(this);

        double InnerRig = trk.getRigidity();
        double L1       = trk.getRigidity(1,2,2);
        double cutOffRig= rti.getCutoffRigidity();
        double richBeta = ric.getBeta();
        double TOFBeta  = tof.getBeta();

        // --- ISS 曝光时间 histogram 填充 ---
        if (isISS){
            if (!rti.cutRTI().total) continue;
            if (std::find(timeTag.begin(), timeTag.end(), time[0]) == timeTag.end()){
                float exposureTime = rti.calculateExposure().value;

                for (size_t c=0;c<chains.size();++c){
                    TH1F* h_exp_rig = histManager->FLUXH2[c];
                    if (h_exp_rig){
                        double rigCut = Constants::SAFE_FACTOR_RIG * cutOffRig;
                        for (int ibin=1; ibin<=h_exp_rig->GetNbinsX(); ++ibin){
                            if (h_exp_rig->GetBinLowEdge(ibin) >= rigCut) {
                                for (int ib=ibin; ib<=h_exp_rig->GetNbinsX(); ++ib){
                                    h_exp_rig->SetBinContent(ib, h_exp_rig->GetBinContent(ib)+exposureTime);
                                }
                                break;
                            }
                        }
                    }
                    for (int d=0; d<3; ++d){
                        for (int i=0;i<iso->getIsotopeCount();++i){
                            int mass=iso->getMass(i);
                            const auto& rb=binMgr.GetIsotopeBetaBins(mass);
                            double betaCO = Tools::rigidityToBeta(cutOffRig, charge, mass, false);
                            double betaCut = (d==0?Detector::BetaTypes[0].getSafetyFactor():d==1?Detector::BetaTypes[1].getSafetyFactor():Detector::BetaTypes[2].getSafetyFactor())*betaCO;
                            TH1F* h=histManager->FLUXH3[c][d][i];
                            if (h){
                                for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin){
                                    if (rb[ibin-1] >= betaCut){
                                        for (int ip=ibin; ip<=h->GetNbinsX(); ++ip){
                                            h->SetBinContent(ip, h->GetBinContent(ip)+exposureTime);
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                timeTag.push_back(time[0]);
            }
        }

        // --- MC RICH beta 修正 ---
        if(!isISS){
            auto mp = ric.getModifiedPosition(true);
            Rad rad = (rich_NaF) ? NAF : AGL;
            int is_mc = 1;
            double corr = ModelManager::corrected_beta(richBeta, rad, run, charge, mp[0], mp[1],
                                                       rich_theta, rich_phi, rich_usedm, rich_hit, is_mc);
            richBeta = corr;
        }

        // --- Monitor ---
        if (jentry % 1000000 == 0)
            std::cout<<"Processing entry "<<jentry<<"/"<<nentries<<std::endl;

        double NaFBeta = rich_NaF ? richBeta : -1;
        double AGLBeta = !rich_NaF ? richBeta : -1;

        bool beyondCutoffRig = true;
        if (isISS){
            TH1F* h_exp_rig = histManager->FLUXH2[0];
            double binLow = h_exp_rig ? h_exp_rig->GetBinLowEdge(h_exp_rig->FindBin(L1)) : 0.0;
            beyondCutoffRig = binLow > Constants::SAFE_FACTOR_RIG * cutOffRig;
        }

        auto TrkRes = trk.cutTracker(charge, isISS);
        auto TOFRes = tof.cutTOF(charge, isISS);
        auto RICHRes= ric.cutRICH(charge, isISS, true);

        // --- ID 区域 histogram 填充 ---
        if (TrkRes.total && RICHRes.total) {
            double invRich = 1.0 / richBeta;
            bool passRigTh = !isISS ? (L1 > 150)
                                    : (beyondCutoffRig && (rich_NaF ? (L1>100) : (L1>200)));
            if (passRigTh) {
                for (size_t c=0;c<chains.size();++c) {
                    if (rich_NaF) histManager->IDH4a[c]->Fill(invRich);
                    else          histManager->IDH4b[c]->Fill(invRich);
                }
            }
            if (passRigTh && TOFRes.total) {
                double dB = (1.0/TOFBeta) - invRich;
                for (size_t c=0;c<chains.size();++c) {
                    if (rich_NaF) histManager->IDH6a[c]->Fill(InnerRig, dB);
                    else          histManager->IDH6b[c]->Fill(InnerRig, dB);
                }
            }
        }

        // --- MC 专用 FLUX 累积 ---
        if (!isISS){
            double generatedRig = (mch!=0)?(mmom/mch):0;
            double generatedEk = Tools::rigidityToKineticEnergy(generatedRig, mch, UseMass);

            for (size_t c=0;c<chains.size();++c){
                // CutEff 填充：用 FLUXH1
                for (size_t cg=0;cg<3;++cg){
                    for (size_t nd=0;nd<2;++nd){
                        for (size_t d=0;d<3;++d){
                            histManager->FLUXH1[c][cg][nd][d][0]->Fill(generatedEk);
                        }
                    }
                }
            }
        }
    }

    // --- MC total events ---
    if(!isISS){
        double total_events = std::accumulate(mc_events.begin(), mc_events.end(), 0.0);
        for (size_t c=0;c<chains.size();++c){
            TH1F* h = histManager->MC_FLUXH3[c];
            if (h){ h->SetBinContent(1, total_events - 2);
                std::cout<<"Total MC Number for "<<chains[c]<<": "<<h->GetBinContent(1)<<std::endl;
            }
        }
    }

    std::cout<<"Event processing completed"<<std::endl;
}