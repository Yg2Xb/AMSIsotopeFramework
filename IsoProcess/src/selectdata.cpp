#define selectdata_cxx
/***********************************************************
 *  File: selectdata.cpp
 *  Implementation of selectdata class for AMS analysis
 *
 *  History:
 *    20241029 - created by ZX.Yan
 *    20250416 - main_bkgyyh select Z = 4~7
 *    20250418 - main_be_lu only select Z = 4
 *    20250419 - normal main ISS
 *    20250419 - normal main MC
 *    20250604 - for bkg valid study
 ***********************************************************/
#include "selectdata.h"
#include "IsotopeAnalyzer.h"
#include "ModelManager.h"
#include "RTICut.h"
#include "TrackerCut.h"
#include "TOFCut.h"
#include "RICHCut.h"
#include "Tool.h"

using namespace AMS_Iso;

void selectdata::SetAnalyzer(IsotopeAnalyzer* analyzer) {
	analyzer_ = analyzer;
}

void selectdata::Loop() {
	if (!fChain || !analyzer_) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout << "Total entries: " << nentries << std::endl;

	// Get basic parameters from analyzer
	int charge = analyzer_->getCharge();  
	int UseMass = analyzer_->getUseMass(); 
	std::cout<<"Using Mass = "<<UseMass<<" for binning and ek convert"<<std::endl;
	bool isISS = analyzer_->isISS();  
	auto saveTree = analyzer_->getSaveTree();

	const IsotopeVar* isotope = analyzer_->getIsotope();
	int tempo_Mass = isotope->getMass(0);
	std::cout<<"Using tempo_Mass = "<<tempo_Mass<<" for cutoff convert"<<std::endl;

	if (!saveTree) {
		std::cerr << "Failed to get Tree" << std::endl;
		return;
	}

	// Basic event histograms
	auto h_event_ek = analyzer_->getHist1D(HistType::Event, "EventEk");
	auto h_measured_tof = analyzer_->getHist1D(HistType::Event, "MeasuredEkTOF");
	auto h_measured_naf = analyzer_->getHist1D(HistType::Event, "MeasuredEkNaF");
	auto h_measured_agl = analyzer_->getHist1D(HistType::Event, "MeasuredEkAGL");

	// MC number histogram
	auto h_mc_num = analyzer_->getHist1D(HistType::MC, "MCNum");

	// Exposure histograms
	auto h_exp_rig = analyzer_->getHist1D(HistType::Exposure, "ExpRig");
	std::vector<std::vector<TH1D*>> h_beta_exposure;  // [isotope][beta_type]
	for(int is = 0; is < analyzer_->getIsotope()->getIsotopeCount(); ++is) {
		std::vector<TH1D*> beta_hists;
		for(int ib = 0; ib < Constants::BETA_TYPES; ++ib) {
			beta_hists.push_back(analyzer_->getBetaExposureHist(is, ib));
		}
		h_beta_exposure.push_back(beta_hists);
	}

	// acceptance histograms
	std::vector<TH1D*> h_tof_cuts(11);
	std::vector<TH1D*> h_naf_cuts(11);
	std::vector<TH1D*> h_agl_cuts(11);
	for(int i = 0; i <= 10; ++i) {
		h_tof_cuts[i] = analyzer_->getEventHist(DetType::TOF, i);
		h_naf_cuts[i] = analyzer_->getEventHist(DetType::NaF, i);
		h_agl_cuts[i] = analyzer_->getEventHist(DetType::AGL, i);
	}

	// TOF special region acc histograms
	std::vector<TH1D*> h_tof_special;
	const std::vector<std::string> special_regions = {
		"Exclude4", "NaF", "isNaF", "NaFGeo", "NaFGeo2"
	};
	for(const auto& region : special_regions) {
		h_tof_special.push_back(
				analyzer_->getHist1D(HistType::Event, Form("EventEkTOF_%s", region.c_str()))
				);
	}

	// Resolution histograms
	std::vector<TH1D*> h_rich_beta_rsl;
	std::vector<TH2D*> h_rich_beta_theta;  // 新增
	std::vector<TH2D*> h_rich_beta_phi;    // 新增
	std::vector<TH2D*> h_tof_beta_ek_rsl;
	std::vector<std::vector<TH2D*>> h_tof_geo_rsl;
	for(int i = 0; i < 2; ++i) {
		// RICH Beta resolution
		h_rich_beta_rsl.push_back(
				analyzer_->getHist1D(HistType::Resolution, 
					Form("RSL_RichBeta_%s", Detector::RichDetectorNames[i].c_str()))
				);
		// RICH Beta vs theta and phi
		h_rich_beta_theta.push_back(
				analyzer_->getHist2D(HistType::Resolution,
					Form("RBeta_Theta_%s", Detector::RichDetectorNames[i].c_str()))
				);

		h_rich_beta_phi.push_back(
				analyzer_->getHist2D(HistType::Resolution,
					Form("RBeta_Phi_%s", Detector::RichDetectorNames[i].c_str()))
				);
		// TOF-RICH Beta resolution
		h_tof_beta_ek_rsl.push_back(
				analyzer_->getHist2D(HistType::Resolution,
					Form("RSL_TOFBeta_Ek_%s", Detector::RichDetectorNames[i].c_str()))
				);
		// TOF GeoCut resolution
		std::vector<TH2D*> geo_hists;
		for(const auto& type : {"AllEdgeTOFBeta_Ek", "NaFGeoTOFBeta_Ek"}) {
			geo_hists.push_back(
					analyzer_->getHist2D(HistType::Resolution,
						Form("RSL_%s_%s", type, Detector::RichDetectorNames[i].c_str()))
					);
		}
		h_tof_geo_rsl.push_back(geo_hists);
	}
	auto h_beta_check = analyzer_->getHist1D(HistType::Resolution, "BetaCheck");
	auto h_rig_check = analyzer_->getHist1D(HistType::Resolution, "RigCheck"); 
	auto h_rigbeta_check = analyzer_->getHist1D(HistType::Resolution, "RigBetaCheck");
	auto h_rigbeta_Ekcheck = analyzer_->getHist2D(HistType::Resolution, "RigBetaCheck_Ek");

	// Efficiency histograms
	std::vector<std::pair<TH1D*, TH1D*>> h_effs;  // pairs of (denominator, numerator)
	for(const auto& type : {"L1PickUp", "RICHRec"}) {
		h_effs.emplace_back(
				analyzer_->getHist1D(HistType::Efficiency, Form("%sEff_den_MC_0", type)),
				analyzer_->getHist1D(HistType::Efficiency, Form("%sEff_num_MC_1", type))
				);
	}

	// Physics trigger histograms
	std::vector<TH1D*> h_triggers;
	for(const auto& type : {"PhysTrig", "TOFPhysTrig", "RICHPhysTrig", "UnPhys", "TOFUnPhys", "RICHUnPhys"}) {
		h_triggers.push_back(analyzer_->getHist1D(HistType::Efficiency, type));
	}

	//charge tempFit hist
	const std::vector<std::string> select_types = {"L1Temp", "L2Temp"};
	const std::vector<std::string> l1_types = {"L1Normal", "L1Unbiased"};
	const std::vector<std::string> coe_values = {"0p2", "0p4", "1p0"};
	const std::vector<std::string> elements = {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
	// const std::vector<std::string> bkg_types = {"Bkg"}; 
	std::vector<TH2D*> h_charges; 
	/*
	for (const auto& elem : elements) {
		// L1Temp histograms (only Coe1p0)
		for (const auto& l1_type : l1_types) {
			const std::string& coe = "1p0"; // L1Temp only uses Coe1p0
			const std::string& select = "L1Temp";
			std::string name = Form("%s_%s_Coe%s_%s_Bkg",
					select.c_str(), l1_type.c_str(), coe.c_str(), elem.c_str());
			h_charges.push_back(analyzer_->getHist2D(HistType::ChargeTempFit, name));
			//std::cout << "Defined: " << name << " at index " << h_charges.size() - 1 << std::endl;
		}
	
		// L2Temp histograms (all Coe values)
		for (const auto& l1_type : l1_types) {
			for (const auto& coe : coe_values) {
				const std::string& select = "L2Temp";
				std::string name = Form("%s_%s_Coe%s_%s_Bkg",
						select.c_str(), l1_type.c_str(), coe.c_str(), elem.c_str());
				h_charges.push_back(analyzer_->getHist2D(HistType::ChargeTempFit, name));
				//std::cout << "Defined: " << name << " at index " << h_charges.size() - 1 << std::endl;
			}
		}
	}
	*/
	
	std::vector<TH1D*> h_frag;
	std::vector<TH1D*> h_source;
	const std::vector<std::string> det_types = {"TOF", "NaF", "AGL"}; 
	    for (const auto& det_type : det_types) {
        std::string fragname = Form("FragNucCounts_%s", det_type.c_str());
        h_frag.push_back(analyzer_->getHist1D(HistType::ChargeTempFit, fragname));
        std::string sourcename = Form("CutL1SourceCounts_%s", det_type.c_str());
        h_source.push_back(analyzer_->getHist1D(HistType::ChargeTempFit, sourcename));
    }

	// Create output branches 
	std::bitset<32> cutStatus;
	unsigned int cutStatusInt = 0;
	double InnerRig = -1, L1InnerRig = -1, TOFBeta = -1, richBeta = -1,  NaFBeta = -1, AGLBeta = -1, cutOffRig = -1;

	//for bkg
	double MC_weight_27 = 1.;
	//2: Complete or loose Be selection
	//3: With or without bkg cut or strict bkg cut(ntrack==1)
	//3 det * 3 iso
	int BkgStudyCut[2][3][3] = {{{0}}};
	int MCInterCut[3] = {0};
	double generatedRig = -1, generatedEk = -1;
	double TOFEk = -1, NaFEk = -1, AGLEk = -1;
	
	saveTree->Branch("cutStatus", &cutStatusInt, "cutStatus/i");
	saveTree->Branch("InnerRig", &InnerRig, "InnerRig/D");
	saveTree->Branch("L1InnerRig", &L1InnerRig, "L1InnerRig/D");
	saveTree->Branch("richBeta", &richBeta, "richBeta/D");

	saveTree->Branch("cutOffRig", &cutOffRig, "cutOffRig/D");

	saveTree->Branch("TOFBeta", &TOFBeta, "TOFBeta/D");
	saveTree->Branch("NaFBeta", &NaFBeta, "NaFBeta/D");
	saveTree->Branch("AGLBeta", &AGLBeta, "AGLBeta/D");

	saveTree->Branch("TOFEk", &TOFEk, "TOFEk/D");
	saveTree->Branch("NaFEk", &NaFEk, "NaFEk/D");
	saveTree->Branch("AGLEk", &AGLEk, "AGLEk/D");
	if(!isISS){
		saveTree->Branch("generatedRig", &generatedRig, "generatedRig/D");
		saveTree->Branch("generatedEk", &generatedEk, "generatedEk/D");

		saveTree->Branch("MC_weight_27", &MC_weight_27, "MC_weight_27/D");
		saveTree->Branch("BkgStudyCut", BkgStudyCut, "BkgStudyCut[2][3][3]/I");
		saveTree->Branch("MCInterCut", MCInterCut, "MCInterCut[3]/I");
	}
	
	//for bkg valid study, type 1 = L1 signal, 2 = L1 Temp, 3 = L2 Temp, Z = 4-8
	bool betacutoffcut[5];
	Float_t trk_ql2, trk_ql1, trk_ql1_unbias;
	Float_t QTrkInner[3], QTrkInnerx[3], QTrkInnery[3]; //0=innerQ, 1and2=L3-L8 aveQ
	double usedEk, usedBeta;
	bool ChargeCutsSelectJudge[5][3];
    saveTree->Branch("trk_ql1",&trk_ql1,"trk_ql1/F");
    saveTree->Branch("trk_ql1_unbias",&trk_ql1_unbias,"trk_ql1_unbias/F");
    saveTree->Branch("trk_ql2",&trk_ql2,"trk_ql2/F");
    saveTree->Branch("QTrkInner",QTrkInner,"QTrkInner[3]/F");   
    saveTree->Branch("QTrkInnerx",QTrkInnerx,"QTrkInnerx[3]/F");   
    saveTree->Branch("QTrkInnery",QTrkInnery,"QTrkInnery[3]/F");   
    saveTree->Branch("usedEk",&usedEk,"usedEk/D");
    saveTree->Branch("usedBeta",&usedBeta,"usedBeta/D");
    saveTree->Branch("betacutoffcut", betacutoffcut, "betacutoffcut[5]/O");
    saveTree->Branch("ChargeCutsSelectJudge", ChargeCutsSelectJudge, "ChargeCutsSelectJudge[5][3]/O");


	// Initialize RICH beta correction model
	std::cout << "Loading models... " << std::endl;
	ModelManager::init("/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium/selection/model_data.root", 
			"/afs/cern.ch/work/z/zuhao/public/yanzx/isotpes_code/lithium/selection/model_mc.root");
	std::cout << "model n:" << ModelManager::model[0][0].index_correction.GetEntries() << std::endl;

	// Variables for MC total number calculation
	UInt_t current_run = 0;  // 改为 UInt_t 类型
	UInt_t min_event = 0;
	UInt_t max_event = 0;
	int event_count = 1;
	std::vector<double> mc_events;

	// exposure calculation
	std::vector<unsigned int> timeTag;
	// Beta bins conversion
	std::vector<std::vector<double>> Rbins_beta_Expo;
	Rbins_beta_Expo.resize(isotope->getIsotopeCount());
	for(int is = 0; is < isotope->getIsotopeCount(); is++) {
		for(int ip = 0; ip <= Constants::RIGIDITY_BINS; ip++) {
			double beta = Tools::kineticEnergyToBeta(
				Binning::KineticEnergyBins[analyzer_->getCharge()-1][is][ip]);
			Rbins_beta_Expo[is].push_back(beta);
		}
	}
	
	std::vector<double> Rbins_beta;
	Rbins_beta = Rbins_beta_Expo[findIsotopeIndex(UseMass, analyzer_->getCharge())]; 
	
	std::vector<double> Rbins_betaBC;
	for (int ip = 0; ip <= Constants::RIGIDITY_BINS; ip++) {
		double beta = Tools::kineticEnergyToBeta(Binning::EkWideBin[ip]);
		Rbins_betaBC.push_back(beta);
	}

	std::vector<double> Rbins_beta_bkg[5]; // 0:Be, 1:B, 2:C, 3:N, 4:O
	for (int iz = 0; iz < 5; ++iz) { // iz=0(Be),...,4(O)
		int z = iz + 3; // z=3(Be),4(B),5(C),6(N),7(O)
		Rbins_beta_bkg[iz].clear();
		for (int ip = 0; ip <= Constants::RIGIDITY_BINS; ++ip) {
			double beta = Tools::kineticEnergyToBeta(Binning::EkWideBin[ip]);
			Rbins_beta_bkg[iz].push_back(beta);
		}
	}
    
	// basic mask: bits 7,8,9,10,12,21,22,23
    unsigned int BKG_BASIC_MASK = 0x0E01780;
	// basic mask: bits 7,8,9,10,16,17,18,23, plus L1 Norm cut, no inner q rms and l1unbias cut
    //unsigned int BKG_BASIC_MASK = 0x870780;
    // TOF Cut mask: bits 26,27
    unsigned int BKG_TOF_MASK = 0xC000000; 
    // RICH Cut mask: bits 28,29
    unsigned int BKG_RICH_MASK = 0x30000000;  
	
	int mass_nuc[5] = {7, 10, 12, 14, 16};

	// Main event loop
	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		fChain->GetEntry(jentry);
		// Calculate MC total number
		if (!isISS) {
			if (current_run != run) {
				if (current_run != 0) {  // Skip first run change
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

				// Handle last entry
				if (jentry == nentries - 1) {
					mc_events.push_back(max_event - min_event + 1 + 
							(max_event - min_event + 1) / event_count);
				}
			}
		}

		// Initialize cut objects
		RTICut rti_cut(this);
		TrackerCut tracker_cut(this);
		TOFCut tof_cut(this);
		RICHCut rich_cut(this);

		// Calculate basic variables
		InnerRig = tracker_cut.getRigidity();
		L1InnerRig = tracker_cut.getRigidity(1,2,2);//GBL V6 L1Inner
		cutOffRig = rti_cut.getCutoffRigidity();
		richBeta = rich_cut.getBeta();
		TOFBeta = tof_cut.getBeta();

		double rigCut = Constants::SAFE_FACTOR_RIG * cutOffRig;

		// Apply RTI cuts for ISS data
		if (isISS && !rti_cut.cutRTI().total) continue;
		
		//exposure time
		if (isISS && (std::find(timeTag.begin(), timeTag.end(), time[0]) == timeTag.end())) {
			float exposureTime = rti_cut.calculateExposure().value;

			// Accumulate exposure time for rigidity bins
			for (int ibin = 1; ibin <= Constants::RIGIDITY_BINS; ibin++) {
				if (Binning::RigidityBins[ibin - 1] >= rigCut) {
					for (int ib = ibin; ib <= Constants::RIGIDITY_BINS; ib++) {
						double currentExp = h_exp_rig->GetBinContent(ib);
						h_exp_rig->SetBinContent(ib, currentExp + exposureTime);
					}
					break;
				}
			}
			
			// Accumulate exposure time for beta bins, only for usemass
			for (int is = 0; is < isotope->getIsotopeCount(); is++) {
				double betaCO = Tools::rigidityToBeta(cutOffRig, charge, isotope->getMass(is), false);
				//using most strict beta cutoff (min mass)
				//double betaCO = Tools::rigidityToBeta(cutOffRig, analyzer_->getCharge(), isotope->getMass(0), false);

				for (int ib = 0; ib < Constants::BETA_TYPES; ib++) {
					double betaCut = Detector::BetaTypes[ib].getSafetyFactor() * betaCO;

					for (int ibin = 1; ibin <= Constants::RIGIDITY_BINS; ibin++) {
						if (Rbins_beta_Expo[is][ibin - 1] >= betaCut) {
							for (int ip = ibin; ip <= Constants::RIGIDITY_BINS; ip++) {
								double currentExp = h_beta_exposure[is][ib]->GetBinContent(ip);
								h_beta_exposure[is][ib]->SetBinContent(ip, currentExp + exposureTime);
							}
							break;
						}
					}
				}
			}
			timeTag.push_back(time[0]);
			
		}
		/*dst version!!!!!!!!!!!!!!*/
		bool yanzx_dst = isISS ? true : false;
		if(!yanzx_dst){
			// Apply RICH beta correction if not corrected in yanzx dst
			auto modiRichPos = rich_cut.getModifiedPosition(true); 
			double modiRichX = modiRichPos[0];
			double modiRichY = modiRichPos[1];
			Rad rad = (rich_NaF) ? NAF : AGL;
			int is_mc = (isISS) ? 0 : 1;
			double rich_beta_corr = ModelManager::corrected_beta(
					richBeta, rad, run, charge, modiRichX, modiRichY,
					rich_theta, rich_phi, rich_usedm, rich_hit, is_mc);

			double corr_cali_richBeta = isISS ? Tools::CorrectCalibrationBiasInData(rich_beta_corr, rich_NaF) : rich_beta_corr ;

			// Print progress every million events
			if (jentry % 1000000 == 0) {
				printf("rich beta before corr=%.6f, after corr=%.6f, after cali_corr=%.6f\n", richBeta, rich_beta_corr, corr_cali_richBeta);
			}
			richBeta = corr_cali_richBeta;
		}
		
		// Print progress every million events
		if (jentry % 1000000 == 0) {
			std::cout << "Processing entry " << jentry << "/" << nentries << std::endl;
		}

		// Calculate basic variables
		NaFBeta = rich_NaF ? richBeta : -1;
		AGLBeta = !rich_NaF ? richBeta : -1;
		NaFEk = rich_NaF ? Tools::betaToKineticEnergy(richBeta) : -1;
		AGLEk = !rich_NaF ? Tools::betaToKineticEnergy(richBeta) : -1;
		TOFEk = Tools::betaToKineticEnergy(TOFBeta);

		// Check cutoff conditions
		double binLowForCutoffRig = h_exp_rig->GetBinLowEdge(h_exp_rig->FindBin(L1InnerRig));
		bool beyondCutoffRig = !isISS ? true : binLowForCutoffRig > Constants::SAFE_FACTOR_RIG * cutOffRig; 

		int tofBin = Tools::findBin(Rbins_beta, TOFBeta);
		int richBin = Tools::findBin(Rbins_beta, richBeta);
		bool beyondCutoffTOF = TOFBeta >= 1 || ((tofBin >= 0) && 
				Tools::isBeyondCutoff(Rbins_beta[tofBin], cutOffRig, Detector::BetaTypes[0].getSafetyFactor(), analyzer_->getCharge(), UseMass, !isISS));
		bool beyondCutoffNaF = richBeta >= 1 || ((richBin >= 0) && 
				Tools::isBeyondCutoff(Rbins_beta[richBin], cutOffRig, Detector::BetaTypes[1].getSafetyFactor(), analyzer_->getCharge(), UseMass, !isISS));
		bool beyondCutoffAGL = richBeta >= 1 || ((richBin >= 0) && 
				Tools::isBeyondCutoff(Rbins_beta[richBin], cutOffRig, Detector::BetaTypes[2].getSafetyFactor(), analyzer_->getCharge(), UseMass, !isISS));

		std::array<bool, 3> beyondCutoffBC;
		int tofBinBC = Tools::findBin(Rbins_betaBC, TOFBeta);
		int richBinBC = Tools::findBin(Rbins_betaBC, richBeta);
		beyondCutoffBC[0] = TOFBeta >= 1 || ((tofBinBC >= 0) && 
				Tools::isBeyondCutoff(Rbins_betaBC[tofBinBC], cutOffRig, 
					Detector::BetaTypes[0].getSafetyFactor(), 5, 10, !isISS));
		for(int i = 1; i < 3; i++) {
			beyondCutoffBC[i] = richBeta >= 1 || ((richBinBC >= 0) && 
					Tools::isBeyondCutoff(Rbins_betaBC[richBinBC], cutOffRig,
						Detector::BetaTypes[i].getSafetyFactor(), 5, 10, !isISS));
		}

		if(!isISS)
		{
			beyondCutoffTOF = beyondCutoffNaF = beyondCutoffAGL = beyondCutoffBC[0] = beyondCutoffBC[1] = beyondCutoffBC[2] = true;	
		}


		auto TrackerCutResult = tracker_cut.cutTracker(charge, isISS);
		auto TOFCutResult = tof_cut.cutTOF(charge, isISS);
		auto RICHCutResult = rich_cut.cutRICH(charge, isISS, true); //interpolate

		auto TOFCutGeo = tof_cut.cutTrapezoidEdges();
		auto RICHCutForBkg = rich_cut.cutRICHforBkg(charge, isISS, true); //interpolate 

		// BkgCuts与z无关，只需要定义一次
		auto onlyBkgCut = tracker_cut.cutBackground(charge, isISS); // 使用任意z都可以
		
		MC_weight_27 = isISS ? 1 : Tools::calculateWeight(mmom, mch, isISS);

		//------------------------------------------

		if(TrackerCutResult.total)
		{
			//rsl hist
			if(RICHCutResult.total) {  
				double invRichBeta = 1. / richBeta;
				//rich rsl hist
				if(beyondCutoffRig) {
					if(rich_NaF && L1InnerRig > 100) {
						h_rich_beta_rsl[0]->Fill(invRichBeta);
					} else if(!rich_NaF && L1InnerRig > 200) {
						h_rich_beta_rsl[1]->Fill(invRichBeta);
					}
					if(rich_NaF && L1InnerRig > 50) {
						h_rich_beta_theta[0]->Fill(rich_theta_tkItr, invRichBeta);
						h_rich_beta_phi[0]->Fill(rich_phi_tkItr, invRichBeta);
					} else if(!rich_NaF && L1InnerRig > 100) {
						h_rich_beta_theta[1]->Fill(rich_theta_tkItr, invRichBeta);
						h_rich_beta_phi[1]->Fill(rich_phi_tkItr, invRichBeta);
					}
				}
				//tof beta rsl hist
				if(beyondCutoffRig && TOFCutResult.total) {
					double deltaBeta = (1. / TOFBeta) - invRichBeta;
					if(rich_NaF) {
						h_tof_beta_ek_rsl[0]->Fill(deltaBeta, NaFBeta*L1InnerRig);
					} else if(!rich_NaF) {
						h_tof_beta_ek_rsl[1]->Fill(deltaBeta, AGLBeta*L1InnerRig);
					}
				}
				if(beyondCutoffRig && tof_cut.cutTOFExcludeLayer4(charge, isISS).total) {
					double deltaBeta = (1. / TOFBeta) - invRichBeta;
					if(rich_NaF) {
						h_tof_geo_rsl[0][0]->Fill(deltaBeta, NaFBeta*L1InnerRig);
					} else if(!rich_NaF) {
						h_tof_geo_rsl[0][1]->Fill(deltaBeta, AGLBeta*L1InnerRig);
					}
				}

			}
			//-----------------------------------
			//saveTree, rsl hist, eff hist  done
			//-----------------------------------
		}
		//using carbon estimate trigger eff
		if(tracker_cut.cutTracker(6,isISS).total && beyondCutoffRig){
			h_triggers[0]->Fill(L1InnerRig);
			if(tof_cut.cutTOF(6, isISS).total) {
				h_triggers[1]->Fill(L1InnerRig);
			}
			if(rich_cut.cutRICH(6, isISS, true).total) {
				h_triggers[2]->Fill(L1InnerRig);
			}
		}
		//fill unphyics
		if(tracker_cut.cutUnphysical(6, isISS).total){
			if(beyondCutoffRig) {
				h_triggers[3]->Fill(L1InnerRig);
			}
			if(tof_cut.cutTOF(6, isISS).total) {
				h_triggers[4]->Fill(L1InnerRig);
			}
			if(rich_cut.cutRICH(6, isISS, true).total) {
				h_triggers[5]->Fill(L1InnerRig);
			}
		}

		//rich beta check
		if (TrackerCutResult.total && RICHCutResult.total && L1InnerRig > 0 && beyondCutoffAGL) 
		{
			double Rbeta_rig = 1. / Tools::rigidityToBeta(L1InnerRig, analyzer_->getCharge(), UseMass, false);
			if(L1InnerRig > 150.){
				h_rigbeta_check->Fill(Rbeta_rig - 1./AGLBeta, MC_weight_27);
			}
			h_rigbeta_Ekcheck->Fill(Rbeta_rig - 1./AGLBeta, AGLEk, MC_weight_27);
		}

		//only mc
		if(isISS){
			//if(TrackerCutResult.total) saveTree->Fill(); 
			continue;
		}

		generatedRig = (mmom/mch);
		generatedEk = Tools::rigidityToKineticEnergy(generatedRig, mch, UseMass); 
		// 初始化MCInterCut数组
		for (int iso = 0; iso < 3; iso++) {
			MCInterCut[iso] = 0;
		}
		// 初始化BkgStudyCut数组
		for (int sel = 0; sel < 2; sel++) {
			for (int bkg = 0; bkg < 3; bkg++) {
				for (int det = 0; det < 3; det++) {
					BkgStudyCut[sel][bkg][det] = 0;
				}
			}
		}

		//rich beta check
		if (TrackerCutResult.total && RICHCutResult.total && L1InnerRig > 150) 
		{

			double Rbeta_truth = 1. / Tools::rigidityToBeta(mtrmom[7]/mch, analyzer_->getCharge(), UseMass, false);
			h_beta_check->Fill(Rbeta_truth - 1./AGLBeta, MC_weight_27);

			h_rig_check->Fill(1./(mtrmom[7]/mch) - 1./L1InnerRig, MC_weight_27);
		}


		
		// Tracker cuts (累积式)
		//h_tof_cuts[9]->Fill(generatedEk); 
		//if(L1InnerRig > 0){h_tof_cuts[10]->Fill(generatedEk);}
		bool passTracker = true;// && L1InnerRig > 0;
			for (size_t i = 0; i < 7; ++i) {
			passTracker &= TrackerCutResult.details[i+1];
			if(passTracker) {
				h_tof_cuts[i]->Fill(generatedEk);
				h_naf_cuts[i]->Fill(generatedEk);
				h_agl_cuts[i]->Fill(generatedEk);
			}
		}

		// TOF cuts (累积式)
		bool passTOF = passTracker;  // 继承Tracker的结果
		for (size_t i = 0; i < 2; ++i) {
			passTOF &= TOFCutResult.details[i];
			if(passTOF) {
				h_tof_cuts[i+7]->Fill(generatedEk);
			}
		}

		// RICH cuts (累积式)
		bool passRICH = passTracker;  // 继承Tracker的结果
		for (size_t i = 0; i < 3; ++i) {
			passRICH &= RICHCutResult.details[i];
			if(passRICH && rich_NaF) {
				h_naf_cuts[i+7]->Fill(generatedEk);
			}
			if(passRICH && !rich_NaF) {
				h_agl_cuts[i+7]->Fill(generatedEk);
			}
		}
		

		//--------------------------
		//loop one event end
		//--------------------------
	}
	// Fill MC total number histogram
	if(!isISS) {
		double total_events = std::accumulate(mc_events.begin(), mc_events.end(), 0.0);
		h_mc_num->SetBinContent(1, total_events - 2);
		std::cout << "Total MC Number: " << h_mc_num->GetBinContent(1) << std::endl;
	}

	std::cout << "Event processing completed" << std::endl;

	}
