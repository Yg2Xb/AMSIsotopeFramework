#define selectdata_cxx
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

	std::cout<<"Loading models... "<<std::endl;
	ModelManager::init("/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/model_data.root",
			"/afs/cern.ch/user/z/zixuan/public/AMSIsotopeFramework/SampleProduction/model_mc.root");
	std::cout<<"model n:"<<ModelManager::model[0][0].index_correction.GetEntries()<<std::endl;
	AMS_Iso::Tools::initFluxFunctions();

	std::vector<unsigned int> timeTag;
	UInt_t current_run=0, min_event=0, max_event=0; int event_count=1;
	std::vector<double> mc_events;

    const int NchainLoc = static_cast<int>(chains.size()); // 2
    const int NdetLoc   = 3;                               // 0:TOF, 1:NaF, 2:AGL
    const int NisoLoc   = iso->getIsotopeCount();

	double weight_NucFlux = 1.;
	double cutOffRig = -1, TOFBeta = -1, richBeta = -1;
	double rig_chain[2] = {-1, -1};
	double beta_det[3] = {-1, -1, -1};
	double ek_det[3] = {-1, -1, -1};

    bool beyondBetaCutoff[3][Constants::N_nuc];
    std::map<std::pair<int,int>, int> ZAMap;
    for(int n=0; n<Constants::N_nuc; ++n) ZAMap[{Constants::nuclei_Z[n], Constants::nuclei_A[n]}] = n;
    // --- Lambda 查询函数 ---
    auto getBeyondBetaCutoffCut = [&](int det, int Z, int A){
        if(det<0 || det>=3) return false;
        auto it = ZAMap.find({Z,A});
        if(it==ZAMap.end()) return false;
        return beyondBetaCutoff[det][it->second];
    };

	for (Long64_t jentry=0;jentry<nentries;++jentry){
		Long64_t ientry = LoadTree(jentry); 
		if (ientry<0) break;
		fChain->GetEntry(jentry);

		// --- Monitor ---
		if (jentry % 1000000 == 0)
			std::cout<<"Processing entry "<<jentry<<"/"<<nentries<<std::endl;

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

		// Initialize cut objects
		RTICut rti_cut(this);
		TrackerCut tracker_cut(this);
		TOFCut tof_cut(this);
		RICHCut rich_cut(this);

		// Apply RTI cuts for ISS data
		if (isISS && !rti_cut.cutRTI().total) continue;

		// Calculate basic variables
		rig_chain[0] = tracker_cut.getRigidity(); // GBL V6 Inner
		rig_chain[1] = tracker_cut.getRigidity(1,2,2);//GBL V6 L1Inner
		cutOffRig = rti_cut.getCutoffRigidity();
		richBeta = rich_cut.getBeta();
		TOFBeta = tof_cut.getBeta();

		// --- ISS 曝光时间 histogram 填充 ---
		//ok
		if (isISS && (std::find(timeTag.begin(), timeTag.end(), time[0]) == timeTag.end())) {
			float exposureTime = rti_cut.calculateExposure().value;
			TH1F* h_exp_rig = histManager->ISS_FLUXH2[0].get();
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
					double betaCO = Tools::rigidityToBeta(cutOffRig, charge, mass, false);
					double betaCut = Detector::BetaTypes[d].getSafetyFactor()*betaCO;
					TH1F* h=histManager->ISS_FLUXH3[d][i].get();
					if (h){
						for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin){
							if (h->GetBinLowEdge(ibin) >= betaCut){
								for (int ip=ibin; ip<=h->GetNbinsX(); ++ip){
									h->SetBinContent(ip, h->GetBinContent(ip)+exposureTime);
								}
								break;
							}
						}
					}
				}
			}
			timeTag.push_back(time[0]);
		}

		// --- MC RICH beta 修正 ---
		//ok
		bool yanzx_dst = isISS ? true : false;
		if(!yanzx_dst){
			auto modiRichPos = rich_cut.getModifiedPosition(true); 
			double modiRichX = modiRichPos[0];
			double modiRichY = modiRichPos[1];
			Rad rad = (rich_NaF) ? NAF : AGL;
			int is_mc = (isISS) ? 0 : 1;
			double rich_beta_corr = ModelManager::corrected_beta(
					richBeta, rad, run, charge, modiRichX, modiRichY,
					rich_theta, rich_phi, rich_usedm, rich_hit, is_mc);

			double corr_cali_richBeta = isISS ? Tools::CorrectCalibrationBiasInData(rich_beta_corr, rich_NaF) : rich_beta_corr ;

			if (jentry % 1000000 == 0) {
				printf("rich beta before corr=%.6f, after corr=%.6f, after cali_corr=%.6f\n", richBeta, rich_beta_corr, corr_cali_richBeta);
			}
			richBeta = corr_cali_richBeta;
		}
		if(!isISS) richBeta = Tools::GetSmearRichBeta(charge, richBeta, rich_NaF);
		//ok
		beta_det[0] = TOFBeta;
		beta_det[1] = rich_NaF ? richBeta : -9;
		beta_det[2] = !rich_NaF ? richBeta : -9;

		ek_det[0] = Tools::betaToKineticEnergy(TOFBeta);
		ek_det[1] = rich_NaF ? Tools::betaToKineticEnergy(richBeta) : -9;  
		ek_det[2] = (!rich_NaF) ? Tools::betaToKineticEnergy(richBeta) : -9;
		//ok
		weight_NucFlux = isISS ? 1.0 : Tools::calculateWeight(mmom, mch, UseMass, isISS);
		if(isISS && jentry%100000==0) 
		{
			std::cout<<"mmom="<<mmom<<",mch="<<mch<<",UseMass="<<UseMass<<std::endl;
			std::cout<<"weight_NucFlux="<<weight_NucFlux<<std::endl;
		}
		//---------------------------

		//ok--- Cut application
		auto TrackerCutResult = tracker_cut.cutTracker(charge, isISS);
		bool TwoAccTrackerCutResult[2] = {
			tracker_cut.TwoAccTrackerCut(charge, isISS).details[0],
			tracker_cut.TwoAccTrackerCut(charge, isISS).details[1]
		}; 
		bool BetaDetectorCutResult[3] = {
			tof_cut.cutTOF(charge, isISS).total,
			rich_NaF && rich_cut.cutRICH(charge, isISS, true).total,
			!rich_NaF && rich_cut.cutRICH(charge, isISS, true).total
		};

		bool beyondCutoffRig[2] = {true,true}; // UnbiasedL1Inner, L1Inner
		if (isISS){
			for(int c = 0; c < 2; c++){
				if(rig_chain[c]>0){
					TH1F* h_exp_rig = histManager->ISS_FLUXH2[0].get();
					double binLow = h_exp_rig ? h_exp_rig->GetBinLowEdge(h_exp_rig->FindBin(rig_chain[c])) : -1;
					beyondCutoffRig[c] = binLow > Constants::SAFE_FACTOR_RIG * cutOffRig;
				}
			}
		}

        for(int d = 0; d < 3; ++d){ // detector
            for(int n = 0; n < Constants::N_nuc; ++n){ // 核
                int Z = Constants::nuclei_Z[n]; int A = Constants::nuclei_A[n];
                auto betaBins = binMgr.GetBetaBins(Z,A);
                if(beta_det[d] >= 1){
                    beyondBetaCutoff[d][n] = true;
                    continue;
                }
                int betaBin = Tools::findBin(betaBins, beta_det[d]);
                if(betaBin >= 0){
                    double betaLow = betaBins[betaBin];
                    beyondBetaCutoff[d][n] = Tools::isBeyondCutoff(
                        betaLow, cutOffRig, Detector::BetaTypes[d].getSafetyFactor(), Z, A, !isISS
                    );
                } else {
                    beyondBetaCutoff[d][n] = false;
                }
            }
        }

		// --- ID  histogram 填充 ---
        for (int c = 0; c < NchainLoc; ++c) {
            if (!TwoAccTrackerCutResult[c]) continue; // chain gate
            for (int d = 0; d < NdetLoc; ++d) {
                if (!BetaDetectorCutResult[d]) continue; // detector gate
                if (ek_det[d] <= 0 || !Tools::isValidBeta(beta_det[d])) continue;
                // Reconstruct inverse mass; alpha=1.0 nominal
                auto mres = calculateMass(beta_det[d], 1.0, rig_chain[c], charge);
                if (!mres.isValid() || mres.invMass <= 0) continue;
                for (int i = 0; i < NisoLoc; ++i) {
                    const int A = iso->getMass(i);
                    if (!getBeyondBetaCutoffCut(d, charge, A)) continue;
                    if (auto* h = histManager->IDH2[c][d][i].get())
                        h->Fill(mres.invMass, ek_det[d], weight_NucFlux);
                }
            }
        }
		// --- BKG  histogram 填充 ---
		// --- FLUX  histogram 填充 ---
        
        if(isISS){
		    // --- ISS ID histogram 填充 ---
		    //ISS_ID1
            for(int c = 0; c < 2; c++){ // UnbiasedL1Inner, L1Inner
                for(int d = 0; d < 3; d++){ // TOF, NaF, AGL
                    for(int i = 0; i < iso->getIsotopeCount(); i++){
                        if(getBeyondBetaCutoffCut(d, charge, iso->getMass(i)) && TwoAccTrackerCutResult[c] && BetaDetectorCutResult[d]){
                            histManager->ISS_IDH1[c][d][i]->Fill(ek_det[d]);
                        }
                    }
                }
            }
		    // --- ISS BKG histogram 填充 ---

            // --- ISS FLUX histogram 填充 ---
        
        } 
		// --- MC ID histogram 填充 ---
        if(!isISS)
        {
            for (int c = 0; c < Nchain; ++c) {
                for (int d = 0; d < Ndet; ++d) {
                    if(){
                        histManager->MC_IDH1[c][d];
                    }
            }
        }
        }
		// --- MC BKG histogram 填充 ---
		
		// --- MC FLUX histogram 填充 ---
		if (!isISS){
			double generatedRig = (mch!=0)?(mmom/mch):0;
			double generatedEk = Tools::rigidityToKineticEnergy(generatedRig, mch, UseMass);
			for (size_t c=0;c<chains.size();++c){
				for (size_t cg=0;cg<3;++cg){
					for (size_t nd=0;nd<2;++nd){
						for (size_t d=0;d<3;++d){
							//...
						}
					}
				}
			}
		}
	}

	// --- MC total events ---
	if (!isISS) {
		std::string fluxName = Tools::selectFluxName(charge, UseMass);
		TF1* f_flux = nullptr;
		double flux_norm = 1.0;
		if (!fluxName.empty() && Tools::getFluxMap().count(fluxName)) {
			f_flux = Tools::getFluxMap()[fluxName].get(); 
			flux_norm = Tools::getFluxNorm()[fluxName];
		} else {
			std::cerr << "[ERROR] No flux TF1 for (Z="<<charge<<", A="<<UseMass<<")"<<std::endl;
			return;
		}

		TH1F* h_flux = histManager->MC_FLUXH3[0].get();
		if (!h_flux) {
			std::cerr << "[ERROR] MC_FLUXH3 histogram not found." << std::endl;
			return;
		}

		const TAxis* xAxis = h_flux->GetXaxis();
		const double* ekBins = xAxis->GetXbins()->GetArray();

		double N_gen = std::accumulate(mc_events.begin(), mc_events.end(), 0.0) - 2;
		double base_fluxIntegral = f_flux->Integral(Tools::geneRig_low, Tools::geneRig_up);

		for (int j = 0; j < xAxis->GetNbins(); ++j) {
			double EkLow = ekBins[j];
			double EkUp  = ekBins[j+1];

			double Rlow = Tools::betaToRigidity(
					Tools::kineticEnergyToBeta(EkLow),
					charge, UseMass, false);
			double Rup  = Tools::betaToRigidity(
					Tools::kineticEnergyToBeta(EkUp),
					charge, UseMass, false);

			if (Rlow < Tools::geneRig_low) Rlow = Tools::geneRig_low;
			if (Rup  > Tools::geneRig_up)  Rup  = Tools::geneRig_up;

			if (Rup <= Rlow) continue;

			double fluxIntegral_bin = f_flux->Integral(Rlow, Rup);
			double N_gen_bin = N_gen * (fluxIntegral_bin / base_fluxIntegral);

			h_flux->SetBinContent(j+1, N_gen_bin);
		}

		std::cout << "[INFO] Filled MC_FLUXH3 (generated spectrum based on flux "
			<< fluxName << ")" << std::endl;
	}

	AMS_Iso::Tools::cleanupFluxFunctions();
	std::cout<<"Event processing completed"<<std::endl;
}

