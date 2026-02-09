#include "HistManager.h"
#include "BinningManager.h"
#include <iostream>
#include <stdexcept>
#include "TString.h"
#include <tuple>
#include <vector>
#include <functional>
#include <algorithm>

using namespace AMS_Iso;

static inline const std::vector<double>& safeBins(const std::string& key) {
	try {
		const auto& bins = BinningManager::GetInstance().Get(key);
		if (!bins.empty()) return bins;
	} catch (const std::exception& e) {
		std::cerr << "Warning: Binning key '" << key << "' not found. Using default bins. Msg: " << e.what() << std::endl;
	} catch (...) {
		std::cerr << "Warning: Unknown error getting bins for key '" << key << "'. Using default bins." << std::endl;
	}
	static const std::vector<double> default_bins = {0, 1};
	return default_bins;
}

std::vector<std::string> HistManager::GetActiveBranches() {
    return {
        // RTI
        "run", "isreal", "event", "mcutoffi", "time",// "zenith", "isbadrun",
        //"issaa","rtilf", "rtinev", "rtinerr", "rtintrig","rtinpar", "rtigood", "rtinexl", "irti",
        
        // Tracker
        "physbpatt2", "itrtrack", "ntrack", "tk_hitb", "tk_qin", "tk_qrmn", "tk_qln", "tk_qls",
        "tk_exqln", "tk_exqls", "tk_rigidity1", "tk_chis1", "tk_pos", "tk_dir",
        "betah2hb", "betah2r","btstat","btstat_pR","btstat_new","cutoffpi",
        
        // TOF
        "ibetah", "tof_betah", "tof_btype", "tof_ql", "tof_barid", "tof_pos",
        "tof_pass", "tof_chisc", "tof_chist", "tof_trapezoidedge", "tof_edge",
        "tof_goodgeo",
        
        // RICH
        "rich_beta", "rich_q", "rich_pb", "rich_pmt", "rich_npe", "rich_good", "rich_clean",
        "rich_NaF", "rich_pos", "rich_theta", "rich_phi", "rich_usedm", "rich_hit",
        "rich_tile", "rich_goodgeo",
        
        // MC
        "mmom", "mch", "mtrpar", "mtrmom", "mtrz"
    };
}

HistManager::HistManager(const std::string& output_filename,
		bool isISS,
		const std::vector<std::string>& chains,
		int charge,
		const IsotopeVar* iso,
		int UseMass, int FragmentZ) {
	m_outputFile = std::make_unique<TFile>(output_filename.c_str(), "RECREATE");
	if (!m_outputFile || m_outputFile->IsZombie()) {
		throw std::runtime_error("Failed to create output ROOT file: " + output_filename);
	}
	std::cout << "HistManager: Output file '" << output_filename << "' opened." << std::endl;

	auto& binMgr = BinningManager::GetInstance();

	const int fragZ = isISS ? charge : FragmentZ;
	const std::vector<int> FragA = [&]() {
		const auto& weights = isotopeWeights.at(fragZ);
		std::vector<int> result(weights.size());
		std::transform(weights.begin(), weights.end(), result.begin(), 
			[](const auto& p) { return p.first; });
		return result;
	}();

	int Nchain = static_cast<int>(chains.size());
	int Ndet = static_cast<int>(detectors.size());
	int Niso = iso->getIsotopeCount();
	int Nsrc = isISS ? static_cast<int>(sources.size()) : 1; // MC source dim is 1
	int Nct = static_cast<int>(charge_types.size());
	int NcutGroups = static_cast<int>(cut_groups.size());
	int NnumDen = static_cast<int>(num_den.size());
	int NgeneRec = static_cast<int>(gene_rec.size());
	
	auto ekBinsStd = safeBins("EkPerNucleon");
	auto rigBins = safeBins("Rigidity");
	auto betaBins = binMgr.GetBetaBins(4, 7);

	std::cout<<"DEBUG: Niso="<<Niso<<" charge="<<charge<<std::endl;
	std::cout<<"DEBUG: iso ptr="<<iso<<std::endl;

	// ============ ID AREA ============
	/*
	std::cout<<"DEBUG: ID Hists"<<std::endl;
	
	auto createMass2D = [&](auto& container, const char* prefix, const char* titleTemplate) {
		container.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			container[c].resize(Ndet);
			for (int d = 0; d < Ndet; ++d) {
				container[c][d].resize(Niso);
				for (int i = 0; i < Niso; ++i) {
					int mass = iso->getMass(i);
					//auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
					container[c][d][i] = createHist<TH2F>(
						Form("%s_%s_%s_Mass%dBin", chains[c].c_str(), prefix, detectors[d].c_str(), mass),
						Form(titleTemplate, chains[c].c_str(), detectors[d].c_str(), mass, detectors[d].c_str(), detectors[d].c_str()),
						200, 0, 0.5, 
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				}
			}
		}
	};
	
	// ISS_IDH1 / MC_IDH1
	if (isISS) {
		ISS_IDH1.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			ISS_IDH1[c].resize(Ndet);
			for (int d = 0; d < Ndet; ++d) {
				ISS_IDH1[c][d].resize(Niso);
				for (int i = 0; i < Niso; ++i) {
					int mass = iso->getMass(i);
					//auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
					ISS_IDH1[c][d][i] = createHist<TH1F>(
						Form("%s_ISS_ID_H1_%s_Mass%dBin", chains[c].c_str(), detectors[d].c_str(), mass),
						Form("%s %s Mass%dBin isotope counts;Counts;%s E_{k}/n [GeV/n]",
							chains[c].c_str(), detectors[d].c_str(), mass, detectors[d].c_str()),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				}
			}
		}
	} else {
		createMass2D(MC_IDH1, "MC_ID_H1", "%s %s Mass%dBin Isotope MC NoFragCut 1/Mass template;%s 1/Mass;%s E_{k}/n [GeV/n]");
	}

	// IDH2 & IDH3
	createMass2D(IDH2, "ID_H2", "%s %s UseMass%dBin Isotope MC 1/Mass vs E_{k}/n;%s 1/Mass;%s E_{k}/n [GeV/n]");
	IDH3.resize(Nchain);
	for (int c = 0; c < Nchain; ++c) {
		IDH3[c].resize(Ndet);
		for (int d = 0; d < Ndet; ++d) {
			IDH3[c][d] = createHist<TH2F>(
				Form("%s_ID_H3_%s", chains[c].c_str(), detectors[d].c_str()),
				Form("%s %s Heaviest iso 1/Mass;%s 1/Mass;%s E_{k}/n [GeV/n]",
					chains[c].c_str(), detectors[d].c_str(), detectors[d].c_str(), detectors[d].c_str()),
				200, 0, 0.5, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
		}
	}
	//--------- Beta Study Histograms ---------
	int zMin = isISS ? 2 : charge;
	int zMax = isISS ? 8 : charge;
	int nZ = zMax - zMin + 1;

	// IDH4a/b
	auto createRICH1D_Z = [&](auto& container, const char* suffix, const char* material, int richIdx) {
		container.resize(nZ);
		for (int iz = 0; iz < nZ; ++iz) {
			int z = zMin + iz;
			int idx = z - 1; 
			const std::string& elemName = isISS ? sources[iz] : sources[charge - 2];
			container[iz].resize(Nchain);
			for (int c = 0; c < Nchain; ++c) {
				container[iz][c] = createHist<TH1F>(
					Form("%s_%s_ID_H4%s", chains[c].c_str(), elemName.c_str(), suffix),
					Form("%s %s %s 1/#beta (beta~1);1/#beta (%s);Counts", 
						chains[c].c_str(), elemName.c_str(), material, material),
					100*Detector::RichBins[richIdx][idx], 
					1.0 - Detector::RichAxis[richIdx], 
					1.0 + Detector::RichAxis[richIdx]);
			}
		}
	};
	createRICH1D_Z(IDH4a, "a", "NaF", 0);
	createRICH1D_Z(IDH4b, "b", "AGL", 1);
	// IDH5/6
	auto createDeltaBeta2D_Z = [&](auto& container, const char* prefix, const char* suffix, 
								const char* titleFmt, const std::vector<double>& yBins) {
		container.resize(nZ);
		for (int iz = 0; iz < nZ; ++iz) {
			int z = zMin + iz;
			const std::string& elemName = isISS ? sources[iz] : sources[charge - 2];
			container[iz].resize(Nchain);
			for (int c = 0; c < Nchain; ++c) {
				container[iz][c] = createHist<TH2F>(
					Form("%s_%s_ID_%s%s", chains[c].c_str(), elemName.c_str(), prefix, suffix),
					Form(titleFmt, chains[c].c_str(), elemName.c_str()), 
					80000, -0.4, 0.4,
					static_cast<int>(yBins.size()) - 1, yBins.data());
			}
		}
	};

	createDeltaBeta2D_Z(IDH5a, "H5", "a", "%s %s NaF-Tracker #Delta(1/#beta);NaF-Tracker #Delta(1/#beta);Rigidity [GV]", rigBins);
	createDeltaBeta2D_Z(IDH5b, "H5", "b", "%s %s AGL-Tracker #Delta(1/#beta);AGL-Tracker #Delta(1/#beta);Rigidity [GV]", rigBins);
	if(!isISS)
	{
		createDeltaBeta2D_Z(IDH5a2, "H5", "a2", "%s %s NaF-Tracker #Delta(1/#beta);NaF-Tracker #Delta(1/#beta);gene Rigidity [GV]", rigBins);
		createDeltaBeta2D_Z(IDH5b2, "H5", "b2", "%s %s AGL-Tracker #Delta(1/#beta);AGL-Tracker #Delta(1/#beta);gene Rigidity [GV]", rigBins);
		
		createDeltaBeta2D_Z(IDH5a3, "H5", "a3", "%s %s NaF-Gene #Delta(1/#beta);NaF-Gene #Delta(1/#beta);gene Beta", betaBins);
		createDeltaBeta2D_Z(IDH5b3, "H5", "b3", "%s %s AGL-Gene #Delta(1/#beta);AGL-Gene #Delta(1/#beta);gene Beta", betaBins);
		createDeltaBeta2D_Z(IDH5c3, "H5", "c3", "%s %s TOF-Gene #Delta(1/#beta);TOF-Gene #Delta(1/#beta);gene Beta", betaBins);
		createDeltaBeta2D_Z(IDH5d3, "H5", "d3", "%s %s Track-Gene #Delta(1/rig);Track-Gene #Delta(1/rig);gene Rigidity [GV]", rigBins);
		
		createDeltaBeta2D_Z(IDH5a4, "H5", "a4", "%s %s frag NaF-Gene #Delta(1/#beta);frag NaF-Gene #Delta(1/#beta);L2True Beta", betaBins);
		createDeltaBeta2D_Z(IDH5b4, "H5", "b4", "%s %s frag AGL-Gene #Delta(1/#beta);frag AGL-Gene #Delta(1/#beta);L2True Beta", betaBins);
		createDeltaBeta2D_Z(IDH5c4, "H5", "c4", "%s %s frag TOF-Gene #Delta(1/#beta);frag TOF-Gene #Delta(1/#beta);L2True Beta", betaBins);
		createDeltaBeta2D_Z(IDH5d4, "H5", "d4", "%s %s frag Track-Gene #Delta(1/rig);frag Track-Gene #Delta(1/rig);L2True Rigidity [GV]", rigBins);
		
		IDH7.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			IDH7[c].resize(Ndet+1);
			for (int d = 0; d < Ndet; ++d) {
				IDH7[c][d] = createHist<TH2F>(
					Form("%s_IDH7_%s", chains[c].c_str(), detectors[d].c_str()),
					Form("%s %s Frag Rig Change; (GeneRig - L2TruthRig)/GeneRig; gene Rigidity [GV]",chains[c].c_str(), detectors[d].c_str()),
					2000, -0.5, 0.5, static_cast<int>(betaBins.size()) - 1, betaBins.data());
			}
			IDH7[c][3] = createHist<TH2F>(
				Form("%s_IDH7_Tracker", chains[c].c_str()),
				Form("%s Tracker Frag Rig Change; (GeneRig - L2TruthRig)/GeneRig; gene Rigidity [GV]",chains[c].c_str()),
				2000, -0.5, 0.5, static_cast<int>(betaBins.size()) - 1, betaBins.data());
		}	
	}

	createDeltaBeta2D_Z(IDH6a, "H6", "a", "%s %s TOF-NaF #Delta(1/#beta);TOF-NaF Delta(1/#beta);NaF #beta [GeV/n]", betaBins);
	createDeltaBeta2D_Z(IDH6b, "H6", "b", "%s %s TOF-AGL #Delta(1/#beta);TOF-AGL #Delta(1/#beta);AGL #beta [GeV/n]", betaBins);
	*/

	/*
	// ============ BKG AREA ============
	std::cout<<"DEBUG: BKG Hists"<<std::endl;

	// Initialize containers
	if(!isISS){
		BKG_H1a.resize(Nchain);
		BKG_H1b.resize(Nchain);
		BKG_H1b2.resize(Nchain);
		BKG_H1c.resize(Nchain);
	}
	BKG_H2b.resize(Nchain);
	BKG_H2b2.resize(Nchain);
	BKG_H4.resize(Nchain);
	for(int c=0; c<Nchain; ++c) {
		if(!isISS){
			BKG_H1a[c].resize(Nsrc);
			BKG_H1b[c].resize(Nsrc);
			BKG_H1b2[c].resize(Nsrc);
			BKG_H1c[c].resize(Nsrc);
		}
		BKG_H2b[c].resize(Nsrc);
		BKG_H2b2[c].resize(Nsrc);
		BKG_H4[c].resize(Nsrc);

		for(int s=0; s<Nsrc; ++s) {
			// Determine source name: MC uses own charge, ISS uses list
			std::string srcName = isISS ? sources[s] : sources[charge-2]; 
			if(!isISS){
				BKG_H1a[c][s].resize(Ndet);
				BKG_H1b[c][s].resize(Ndet);
				BKG_H1b2[c][s].resize(Ndet);
				BKG_H1c[c][s].resize(Ndet);
			}
			BKG_H2b[c][s].resize(Ndet);
			BKG_H2b2[c][s].resize(Ndet);
			BKG_H4[c][s].resize(Ndet);

			for(int d=0; d<Ndet; ++d) {
				// --- BKG_H1 Series (Denominators) --- only mc need it, because iss use L1-InnerQ 2D hist to get everything!
				if(!isISS){
					BKG_H1a[c][s][d] = createHist<TH1F>(
						Form("%s_BKG_H1a_%s_%s", chains[c].c_str(), srcName.c_str(), detectors[d].c_str()),
						Form("%s %s %s Truth L1 X->L2 Any;E_{k}/n [GeV/n];Counts", chains[c].c_str(), detectors[d].c_str(), srcName.c_str()),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

					BKG_H1b[c][s][d] = createHist<TH1F>(
						Form("%s_BKG_H1b_%s_%s", chains[c].c_str(), srcName.c_str(), detectors[d].c_str()),
						Form("%s %s %s Truth L1 X->L2 Still X;E_{k}/n [GeV/n];Counts", chains[c].c_str(), detectors[d].c_str(), srcName.c_str()),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					
					BKG_H1b2[c][s][d] = createHist<TH1F>(
						Form("%s_BKG_H1b2_%s_%s", chains[c].c_str(), srcName.c_str(), detectors[d].c_str()),
						Form("%s %s %s Truth L1 X->L2 Still X fullcut; E_{k}/n [GeV/n]; Counts", chains[c].c_str(), detectors[d].c_str(), srcName.c_str()),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

					BKG_H1c[c][s][d] = createHist<TH1F>(
						Form("%s_BKG_H1c_%s_%s", chains[c].c_str(), srcName.c_str(), detectors[d].c_str()),
						Form("%s %s %s Truth L1 X->L2 Any Frag ;E_{k}/n [GeV/n];Counts", chains[c].c_str(), detectors[d].c_str(), srcName.c_str()),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				}

				BKG_H2b[c][s][d] = createHist<TH2F>(
					Form("%s_BKG_H2b_%s_%s", chains[c].c_str(), srcName.c_str(), detectors[d].c_str()),
					Form("%s %s %s L1 X->L2 Frag Elem Y 1/Mass (L1Q Window Cut and loose selection);1/Mass;E_{k}/n [GeV/n]", chains[c].c_str(), detectors[d].c_str(), srcName.c_str()),
					200, 0, 0.5, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				BKG_H2b2[c][s][d] = createHist<TH2F>(
					Form("%s_BKG_H2b2_%s_%s", chains[c].c_str(), srcName.c_str(), detectors[d].c_str()),
					Form("%s %s %s L1 X->L2 Frag Elem Y 1/Mass (L1Q Window Cut and Standard Selection);1/Mass;E_{k}/n [GeV/n]", chains[c].c_str(), detectors[d].c_str(), srcName.c_str()),
					200, 0, 0.5, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

				// --- BKG_H4 Series (Charge Study) ---
				BKG_H4[c][s][d].resize(Nct);
				for(int t=0; t<Nct; ++t) {
					int x_bins = 800; double x_min = 1.5; double x_max = 9.5;
					BKG_H4[c][s][d][t] = createHist<TH2F>(
						Form("%s_BKG_H4_%s_%s_%s", chains[c].c_str(), srcName.c_str(), charge_types[t].c_str(), detectors[d].c_str()),
						Form("%s %s %s %s Charge vs Ek;Charge;E_{k}/n [GeV/n]", chains[c].c_str(), detectors[d].c_str(), srcName.c_str(), charge_types[t].c_str()),
						x_bins, x_min, x_max, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				}
			}
		}
	}
	
	double binStep = 0.01;
    double valMin = 0.5,  valMax = isISS ? 9. : charge + 1.;
    int nBinsVal = static_cast<int>((valMax - valMin) / binStep + 0.5); // +0.5 用于防止浮点误差
    std::vector<double> binsVal(nBinsVal + 1);
    for (int i = 0; i <= nBinsVal; ++i) { binsVal[i] = valMin + i * binStep;}
    double diffMin = isISS ? -8. : -1.*charge, diffMax = isISS ? 8. : 1.*charge;
    int nBinsDiff = static_cast<int>((diffMax - diffMin) / binStep + 0.5);
    std::vector<double> binsDiff(nBinsDiff + 1);
    for (int i = 0; i <= nBinsDiff; ++i) { binsDiff[i] = diffMin + i * binStep;}
	BKG_H5.resize(Nchain);
    for(int c = 0; c < Nchain; ++c) {
        BKG_H5[c].resize(Ndet);
        for(int d = 0; d < Ndet; ++d) {
            BKG_H5[c][d].resize(2); 
            for(int t = 0; t < 1; ++t) {
                const char* viewType; const char* xTitle; const char* yTitle; 
				const std::vector<double>* xBinsPtr; const std::vector<double>* yBinsPtr;
                if (t == 0) {
                    viewType = "Standard"; xTitle = "InnerQ"; yTitle = "L1Q"; xBinsPtr = &binsVal; yBinsPtr = &binsVal; 
                } else {
                    viewType = "Rotated45d"; xTitle = "(L1Q+InnerQ)/2"; yTitle = "L1Q-InnerQ"; xBinsPtr = &binsVal; yBinsPtr = &binsDiff; 
                }
                BKG_H5[c][d][t] = createHist<TH3F>(
                    Form("%s_BKG_H5_%s_type%d", chains[c].c_str(), detectors[d].c_str(), t),
                    Form("%s %s %s;%s;%s;E_{k}/n", chains[c].c_str(), detectors[d].c_str(), viewType, xTitle, yTitle),
                    static_cast<int>(xBinsPtr->size()) - 1, xBinsPtr->data(), 
                    static_cast<int>(yBinsPtr->size()) - 1, yBinsPtr->data(), 
                    static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data() 
                );
            }
        }
    }
	
	// --- MC Only BKG Histograms ---
	if(!isISS) {
		int NisoBKG = static_cast<int>(FragA.size());
		
		BKG_H2a.resize(Nchain);
		BKG_H2a2.resize(Nchain);
		BKG_H3a.resize(Nchain);
		BKG_H3b.resize(Nchain);

		for(int c=0; c<Nchain; ++c) {
			BKG_H2a[c].resize(Ndet);
			BKG_H2a2[c].resize(Ndet);
			BKG_H3a[c].resize(Ndet);
			BKG_H3b[c].resize(Ndet);

			for(int d=0; d<Ndet; ++d) {
				BKG_H2a[c][d].resize(NisoBKG);
				BKG_H2a2[c][d].resize(NisoBKG);
				BKG_H3a[c][d].resize(NisoBKG);
				BKG_H3b[c][d].resize(NisoBKG);

				for(int i=0; i<NisoBKG; ++i) {
					int massA = FragA[i];
					
					BKG_H2a[c][d][i] = createHist<TH1F>(
						Form("%s_BKG_H2a_%s_Z%d_Mass%d", chains[c].c_str(), detectors[d].c_str(), fragZ, massA),
						Form("%s %s MC L1 X->L2 Truth Frag Isotope Mass%d Counts;E_{k}/n [GeV/n];Counts", chains[c].c_str(), detectors[d].c_str(), massA),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

					BKG_H2a2[c][d][i] = createHist<TH1F>(
						Form("%s_BKG_H2a2_%s_Z%d_Mass%d", chains[c].c_str(), detectors[d].c_str(), fragZ, massA),
						Form("%s %s MC L1 X->L2 Truth Frag Isotope Mass%d Counts full selection;E_{k}/n [GeV/n];Counts", chains[c].c_str(), detectors[d].c_str(), massA),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

					BKG_H3a[c][d][i] = createHist<TH1F>(
						Form("%s_BKG_H3a_%s_Z%d_Mass%d", chains[c].c_str(), detectors[d].c_str(), fragZ, massA),
						Form("%s %s MC L1 Truth Frag Isotop Mass%d Counts;generated E_{k}/n [GeV/n];Counts", chains[c].c_str(), detectors[d].c_str(), massA),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

					BKG_H3b[c][d][i] = createHist<TH1F>(
						Form("%s_BKG_H3b_%s_Z%d_Mass%d", chains[c].c_str(), detectors[d].c_str(), fragZ, massA),
						Form("%s %s MC L1 Truth Frag Isotop Mass%d Counts;%s E_{k}/n [GeV/n];Counts", chains[c].c_str(), detectors[d].c_str(), massA,  detectors[d].c_str()),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				}
			}
		}
	}
	*/
	
	// ============ FLUX AREA ============
	
	std::cout<<"DEBUG: Flux Hists"<<std::endl;

	// Resize the first dimension for 7 sources (He to O)
    FLUXH1.resize(Nsrc);
    for (int s = 0; s < Nsrc; ++s) {
        int targetZ = isISS ? s + 2 : charge; 
        FLUXH1[s].resize(NcutGroups);

        for (int cg = 0; cg < NcutGroups; ++cg) {
            std::string cut = cut_groups[cg];
            FLUXH1[s][cg].resize(NnumDen);

            for (int nd = 0; nd < NnumDen; ++nd) {
                FLUXH1[s][cg][nd].resize(Ndet + 1);

                for (int d = 0; d < Ndet + 1; ++d) {
                    bool isTracker = (d == Ndet);
                    std::string detName = isTracker ? "Tracker" : detectors[d];

                    // --- Refined Boolean Logic to Flatten the Control Flow ---
                    bool isHe = (targetZ == 2);
                    bool isBetaCut = (cut == "BetaRecQuality");
                    bool isBkgRed = (cut == "BkgReduction");

                    bool shouldCreate = false;
                    if (isHe) {
                        // Helium Rule: Only BkgReduction
                        shouldCreate = isBkgRed && isTracker;
                    } else {
                        // Li-O Rule: Tracker gets all except BetaRecQuality; Detectors get ONLY BetaRecQuality
                        shouldCreate = (isTracker != isBetaCut); 
                    }

                    // Exit early if this combination is not needed
                    if (!shouldCreate) {
                        FLUXH1[s][cg][nd][d] = nullptr;
                        continue;
                    }

                    // --- Setup Binning and Axis Titles ---
                    const double* bins = isTracker ? rigBins.data() : ekBinsStd.data();
                    int nBins = isTracker ? (rigBins.size() - 1) : (ekBinsStd.size() - 1);
                    std::string xTitle = isTracker ? "Rigidity [GV]" : "E_{k}/n [GeV/n]";

                    // --- Single, Unified Creation Call ---
                    FLUXH1[s][cg][nd][d] = createHist<TH1F>(
                        Form("Eff_FLUXH1_%s_%s_%s_Z%d", cut.c_str(), num_den[nd].c_str(), detName.c_str(), targetZ),
                        Form("%s %s %s %s Efficiency (Z=%d);%s;Counts", chains[0].c_str(), detName.c_str(), cut.c_str(), num_den[nd].c_str(), targetZ, xTitle.c_str()),
                        nBins, bins
                    );
                }
            }
        }
    }

		/*
	if (isISS) {
		ISS_FLUXH2.resize(1);
		ISS_FLUXH2[0] = createHist<TH1F>(
			"ISS_FLUX_H2",
			"ISS Exposure time;Rigidity [GV];Exposure Time [s]",
			static_cast<int>(rigBins.size()) - 1, rigBins.data());
		
		ISS_FLUXH3.resize(Ndet);
		for (int d = 0; d < Ndet; ++d) {
			ISS_FLUXH3[d].resize(Niso);
			for (int i = 0; i < Niso; ++i) {
				int mass = iso->getMass(i);
				auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
				ISS_FLUXH3[d][i] = createHist<TH1F>(
					Form("ISS_FLUX_H3_%s_Mass%dBin", detectors[d].c_str(), mass),
					Form("%s Mass%dBin Exposure time;E_{k}/n [GeV/n];Exposure Time [s]",
						detectors[d].c_str(), mass),
					static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
			}
		}
		ISS_FLUXH4.resize(1);
		ISS_FLUXH4[0] = createHist<TH2F>(
			"ISS_FLUXH4","BTstatus vs Run;Run;BTstatus",
			70985, 1305853512,1731763512, 3, 0.5,3.5);
		ISS_FLUXH5.resize(2);
		ISS_FLUXH5[0] = createHist<TH2F>("ISS_FLUX_H5_0","30deg Max CutoffRig vs InnerRig;30deg Max CutoffRig;InnerRig",2500,3,28,2500,3,28);
		ISS_FLUXH5[1] = createHist<TH2F>("ISS_FLUX_H5_1","30deg Max CutoffRig vs L1InnerRig;30deg Max CutoffRig;L1InnerRig",2500,3,28,2500,3,28);
		
	} else {
		MC_FLUXH3.resize(1);
		MC_FLUXH3[0] = createHist<TH1F>(
			"MC_FLUX_H3",
			"MC Generated counts;E_{k}^{gen}/n [GeV/n];Counts",
			static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
	}
		*/
	
	

	std::cout<<"DEBUG Finish Hist Defination"<<std::endl;
}

void HistManager::PrepareFilteredTree(TChain* dataChain) {
    if (!dataChain || !m_outputFile || !m_outputFile->IsOpen()) {
        std::cerr << "[HistManager] Cannot prepare tree: invalid dataChain or output file." << std::endl;
        return;
    }

    std::cout << "[HistManager] Preparing filtered tree with selected branches..." << std::endl;

    auto activeBranches = GetActiveBranches();
    
    dataChain->SetBranchStatus("*", 0);
    for (const auto& br : activeBranches) {
        dataChain->SetBranchStatus(br.c_str(), 1);
    }
    
    m_outputFile->cd();
    m_filteredTree = dataChain->CloneTree(0);
    m_filteredTree->SetName("amstreea");
    m_filteredTree->SetTitle("Filtered AMS events");
    m_filteredTree->SetAutoSave(100 * 1024 * 1024);
    
    std::cout << "[HistManager] Filtered tree prepared with " << activeBranches.size() 
              << " branches. Ready for filling." << std::endl;
	dataChain->SetBranchStatus("*", 1);
}

void HistManager::Save(bool saveTree) {
    if (!m_outputFile || !m_outputFile->IsOpen()) {
        std::cerr << "HistManager: Error: Output file is not open or is invalid." << std::endl;
        return;
    }
    
    m_outputFile->cd();

    auto writeHists = [&](const auto& container) {
        auto impl = [](const auto& cont, auto& self) -> void {
            using T = std::decay_t<decltype(cont)>;
            if constexpr (std::is_same_v<T, std::vector<H1Ptr>> || 
                          std::is_same_v<T, std::vector<H2Ptr>> ||
						  std::is_same_v<T, std::vector<H3Ptr>>) {
                for (const auto& h : cont) if (h) h->Write();
            } else {
                for (const auto& item : cont) self(item, self);
            }
        };
        impl(container, impl);
    };

    if (!ISS_IDH1.empty()) writeHists(ISS_IDH1);
    if (!MC_IDH1.empty()) writeHists(MC_IDH1);
    
	writeHists(IDH2); writeHists(IDH3);
    writeHists(IDH4a); writeHists(IDH4b);
    writeHists(IDH5a); writeHists(IDH5a2); writeHists(IDH5a3); writeHists(IDH5a4);
    writeHists(IDH5b); writeHists(IDH5b2); writeHists(IDH5b3); writeHists(IDH5b4);
    writeHists(IDH5c3); writeHists(IDH5c4); writeHists(IDH5d3); writeHists(IDH5d4);
    writeHists(IDH6a); writeHists(IDH6b);
    if (!IDH7.empty()) writeHists(IDH7);
	if (!BKG_H1a.empty()) {
		writeHists(BKG_H1a); writeHists(BKG_H1b); writeHists(BKG_H1c);
		writeHists(BKG_H2a); writeHists(BKG_H1b2);
    }
	writeHists(BKG_H2b);
	writeHists(BKG_H2b2);
	writeHists(BKG_H4);
	if (!BKG_H5.empty()) writeHists(BKG_H5);


	if (!BKG_H2a2.empty()) writeHists(BKG_H2a2);
	if (!BKG_H3a.empty()) {
		writeHists(BKG_H3a); writeHists(BKG_H3b);
	}
    
    writeHists(FLUXH1);

    if (!ISS_FLUXH2.empty()) { writeHists(ISS_FLUXH2); writeHists(ISS_FLUXH3); writeHists(ISS_FLUXH4); writeHists(ISS_FLUXH5);}
    if (!MC_FLUXH3.empty()) { writeHists(MC_FLUXH3);} //writeHists(MC_FLUXH2); }

    if (saveTree && m_filteredTree) {
        m_outputFile->cd();
        m_filteredTree->Write("", TObject::kOverwrite);
    }
    
    m_outputFile->Close();
}