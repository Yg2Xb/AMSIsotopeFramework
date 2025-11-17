#include "HistManager.h"
#include "BinningManager.h"
#include <iostream>
#include <stdexcept>
#include "TString.h"
#include <tuple>
#include <vector>
#include <functional>

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
        "run", "isreal", "event", "time", "mcutoffi","zenith", "isbadrun",
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
        
        // MC (如果需要)
        "mmom", "mch", "mtrpar", "mtrmom"
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
	int Nsrc = static_cast<int>(sources.size());
	int Nct = static_cast<int>(charge_types.size());
	int NcutGroups = static_cast<int>(cut_groups.size());
	int NnumDen = static_cast<int>(num_den.size());
	int NgeneRec = static_cast<int>(gene_rec.size());
	
	auto ekBinsStd = safeBins("EkPerNucleon");
	auto rigBins = safeBins("Rigidity");

	std::cout<<"DEBUG: Niso="<<Niso<<" charge="<<charge<<std::endl;
	std::cout<<"DEBUG: iso ptr="<<iso<<std::endl;

	// ============ ID 区域 ============
	std::cout<<"DEBUG: ID Hists"<<std::endl;
	
	auto createMass2D = [&](auto& container, const char* prefix, const char* titleTemplate) {
		container.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			container[c].resize(Ndet);
			for (int d = 0; d < Ndet; ++d) {
				container[c][d].resize(Niso);
				for (int i = 0; i < Niso; ++i) {
					int mass = iso->getMass(i);
					auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
					container[c][d][i] = createHist<TH2D>(
						Form("%s_%s_%s_Mass%dBin", chains[c].c_str(), prefix, detectors[d].c_str(), mass),
						Form(titleTemplate, chains[c].c_str(), detectors[d].c_str(), mass, detectors[d].c_str(), detectors[d].c_str()),
						200, 0, 0.5, 
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				}
			}
		}
	};
	
	// IDH1: ISS/MC
	if (isISS) {
		ISS_IDH1.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			ISS_IDH1[c].resize(Ndet);
			for (int d = 0; d < Ndet; ++d) {
				ISS_IDH1[c][d].resize(Niso);
				for (int i = 0; i < Niso; ++i) {
					int mass = iso->getMass(i);
					auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
					ISS_IDH1[c][d][i] = createHist<TH1D>(
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

	// IDH2
	createMass2D(IDH2, "ID_H2", "%s %s UseMass%dBin Isotope MC 1/Mass vs E_{k}/n;%s 1/Mass;%s E_{k}/n [GeV/n]");


	// IDH3
	IDH3.resize(Nchain);
	for (int c = 0; c < Nchain; ++c) {
		IDH3[c].resize(Ndet);
		for (int d = 0; d < Ndet; ++d) {
			int mass_heaviest = iso->getMass(Niso-1);
			auto ekBinsH = binMgr.GetEkPerNucleonBins(charge, mass_heaviest);
			IDH3[c][d] = createHist<TH2D>(
				Form("%s_ID_H3_%s", chains[c].c_str(), detectors[d].c_str()),
				Form("%s %s Heaviest iso 1/Mass;%s 1/Mass;%s E_{k}/n [GeV/n]",
					chains[c].c_str(), detectors[d].c_str(), detectors[d].c_str(), detectors[d].c_str()),
				200, 0, 0.5, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
		}
	}

	//---------beta study histograms---------
	int zMin = isISS ? 2 : charge;
	int zMax = isISS ? 8 : charge;
	int nZ = zMax - zMin + 1;

	// IDH4a/b: RICH 1D - 带Z维度和sources命名
	auto createRICH1D_Z = [&](auto& container, const char* suffix, const char* material, int richIdx) {
		container.resize(nZ);
		for (int iz = 0; iz < nZ; ++iz) {
			int z = zMin + iz;
			int idx = z - 1; // RichBins索引
			const std::string& elemName = sources[iz]; // 使用sources名称
			container[iz].resize(Nchain);
			
			for (int c = 0; c < Nchain; ++c) {
				container[iz][c] = createHist<TH1D>(
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

	// IDH5/6/7 系列: Delta Beta 2D - 带Z维度和sources命名
	auto createDeltaBeta2D_Z = [&](auto& container, const char* prefix, const char* suffix, 
								const char* titleFmt, const std::vector<double>& yBins) {
		container.resize(nZ);
		for (int iz = 0; iz < nZ; ++iz) {
			int z = zMin + iz;
			const std::string& elemName = sources[iz];
			container[iz].resize(Nchain);
			
			for (int c = 0; c < Nchain; ++c) {
				container[iz][c] = createHist<TH2D>(
					Form("%s_%s_ID_%s%s", chains[c].c_str(), elemName.c_str(), prefix, suffix),
					Form(titleFmt, chains[c].c_str(), elemName.c_str()), // 标题包含chain和元素名
					80000, -0.2, 0.2,
					static_cast<int>(yBins.size()) - 1, yBins.data());
			}
		}
	};

	// IDH5 系列: RICH-Tracker Delta Beta
	createDeltaBeta2D_Z(IDH5a, "H5", "a", "%s %s NaF-Tracker #Delta(1/#beta);NaF-Tracker #Delta(1/#beta);Rigidity [GV]", rigBins);
	createDeltaBeta2D_Z(IDH5b, "H5", "b", "%s %s AGL-Tracker #Delta(1/#beta);AGL-Tracker #Delta(1/#beta);Rigidity [GV]", rigBins);
	createDeltaBeta2D_Z(IDH5a2, "H5", "a2", "%s %s NaF-Tracker #Delta(1/#beta);NaF-Tracker #Delta(1/#beta);gene Rigidity [GV]", rigBins);
	createDeltaBeta2D_Z(IDH5b2, "H5", "b2", "%s %s AGL-Tracker #Delta(1/#beta);AGL-Tracker #Delta(1/#beta);gene Rigidity [GV]", rigBins);
	createDeltaBeta2D_Z(IDH5a3, "H5", "a3", "%s %s NaF-true #Delta(1/#beta);NaF-True #Delta(1/#beta);gene Rigidity [GV]", rigBins);
	createDeltaBeta2D_Z(IDH5b3, "H5", "b3", "%s %s AGL-True #Delta(1/#beta);AGL-True #Delta(1/#beta);gene Rigidity [GV]", rigBins);

	// IDH6 系列: TOF-RICH Delta Beta
	createDeltaBeta2D_Z(IDH6a, "H6", "a", "%s %s TOF-NaF #Delta(1/#beta);TOF-NaF Delta(1/#beta);NaF E_{k}/n [GeV/n]", ekBinsStd);
	createDeltaBeta2D_Z(IDH6b, "H6", "b", "%s %s TOF-AGL #Delta(1/#beta);TOF-AGL #Delta(1/#beta);AGL E_{k}/n [GeV/n]", ekBinsStd);


	auto betarigBins = binMgr.GetBetaRigBins(2, 4);
	createDeltaBeta2D_Z(IDH7a, "H7", "a", "%s %s TOF-NaF #Delta(1/#beta);TOF-NaF #Delta(1/#beta);NaF #beta Rig[GV]", betarigBins);
	createDeltaBeta2D_Z(IDH7b, "H7", "b", "%s %s TOF-AGL #Delta(1/#beta);TOF-AGL #Delta(1/#beta);AGL #beta Rig[GV]", betarigBins);

	// ============ BKG 区域 ============
	std::cout<<"DEBUG: BKG Hists"<<std::endl;
	
	if (isISS) {
		ISS_BKGH1.resize(Nchain);
		ISS_BKGH2.resize(Nchain);
		ISS_BKGH3.resize(Nchain);
		ISS_BKGH4.resize(Nchain);
		
		for (int c = 0; c < Nchain; ++c) {
			ISS_BKGH1[c].resize(Nsrc);
			ISS_BKGH2[c].resize(Nsrc);
			ISS_BKGH3[c].resize(Nsrc);
			ISS_BKGH4[c].resize(Nsrc);
			
			for (int s = 0; s < Nsrc; ++s) {
				int zsrc = s + 2;
				ISS_BKGH1[c][s].resize(Ndet);
				ISS_BKGH2[c][s].resize(Ndet);
				ISS_BKGH3[c][s].resize(Ndet);
				ISS_BKGH4[c][s].resize(Ndet);
				
				for (int d = 0; d < Ndet; ++d) {
					ISS_BKGH1[c][s][d] = createHist<TH1D>(
						Form("%s_ISS_BKG_H1_%s_%s", chains[c].c_str(), sources[s].c_str(), detectors[d].c_str()),
						Form("%s %s L1 Source %s counts;%s E_{k}/n [GeV/n];Counts",
							chains[c].c_str(), detectors[d].c_str(), sources[s].c_str(), detectors[d].c_str()),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

					ISS_BKGH2[c][s][d].resize(Nct);
					for (int t = 0; t < Nct; ++t) {
						if(s!=2 && t==0) continue;
						const char* x_title = "TrackerLayer Charge";
						int x_bins = 700; double x_min = 2.5; double x_max = 9.5;
						if (charge_types[t] == "L1QSignal") {
							x_bins = 700; x_min = 2.5; x_max = 9.5;
						}
						ISS_BKGH2[c][s][d][t] = createHist<TH2D>(
							Form("%s_ISS_BKG_H2_%s_%s_%s", chains[c].c_str(), sources[s].c_str(), charge_types[t].c_str(), detectors[d].c_str()),
							Form("%s %s %s %s charge vs E_{k}/n;%s;%s E_{k}/n [GeV/n]",
								chains[c].c_str(), detectors[d].c_str(), charge_types[t].c_str(), sources[s].c_str(), x_title, detectors[d].c_str()),
							x_bins, x_min, x_max, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					}

					ISS_BKGH3[c][s][d] = createHist<TH1D>(
						Form("%s_ISS_BKG_H3_%s_%s", chains[c].c_str(), sources[s].c_str(), detectors[d].c_str()),
						Form("%s %s L1%s L2 frag counts;%s E_{k}/n [GeV/n];Counts",
							chains[c].c_str(), detectors[d].c_str(), sources[s].c_str(), detectors[d].c_str()),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

					ISS_BKGH4[c][s][d].resize(static_cast<int>(FragA.size()));
					for (size_t bi = 0; bi < FragA.size(); ++bi) {
						int A = FragA[bi];
						auto ekBinsFrag = binMgr.GetEkPerNucleonBins(fragZ, A);
						ISS_BKGH4[c][s][d][bi] = createHist<TH2D>(
							Form("%s_ISS_BKG_H4_%s_%s_Z%d_Mass%d", chains[c].c_str(), sources[s].c_str(), detectors[d].c_str(), fragZ, A),
							Form("%s %s L1%s L2frag 1/Mass vs E_{k}/n (Z=%d Mass%d);%s 1/Mass;%s E_{k}/n [GeV/n]",
								chains[c].c_str(), detectors[d].c_str(), sources[s].c_str(), fragZ, A, detectors[d].c_str(), detectors[d].c_str()),
							200, 0, 0.5, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					}
				}
			}
		}
	} else {
		int NisoBKG = static_cast<int>(FragA.size());

		MC_BKGH1.resize(Nchain);
		MC_BKGH1b.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			MC_BKGH1[c].resize(Ndet);
			MC_BKGH1b[c].resize(Ndet);
			for (int d = 0; d < Ndet; ++d) {
				MC_BKGH1[c][d] = createHist<TH1D>(
					Form("%s_MC_BKG_H1_%s", chains[c].c_str(), detectors[d].c_str()),
					Form("%s %s MC input L1 counts;%s E_{k}/n [GeV/n];Counts",
						chains[c].c_str(), detectors[d].c_str(), detectors[d].c_str()),
					static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				MC_BKGH1b[c][d] = createHist<TH1D>(
					Form("%s_MC_BKG_H1b_%s", chains[c].c_str(), detectors[d].c_str()),
					Form("%s %s MC input gene counts;%s E_{k}/n [GeV/n];Counts",
						chains[c].c_str(), detectors[d].c_str(), detectors[d].c_str()),
					static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
			}
		}

		auto createMcBkgHists = [&](auto& container, const char* prefix, const char* titleFmt) {
			container.resize(Nchain);
			for (int c = 0; c < Nchain; ++c) {
				container[c].resize(Ndet);
				for (int d = 0; d < Ndet; ++d) {
					container[c][d].resize(NisoBKG);
					for (int i = 0; i < NisoBKG; ++i) {
						int massA = FragA[i];
						auto ekBinsFrag = binMgr.GetEkPerNucleonBins(fragZ, massA);
						container[c][d][i] = createHist<TH1D>(
							Form("%s_MC_BKG_%s_%s_Z%d_Mass%d", chains[c].c_str(), prefix, detectors[d].c_str(), fragZ, massA),
							Form(titleFmt, chains[c].c_str(), detectors[d].c_str(), massA, detectors[d].c_str()),
							static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					}
				}
			}
		};

		createMcBkgHists(MC_BKGH2, "H2", "%s %s MC L1Orig-L2 frag Isotope Mass%d Counts vs E_{k}/n;%s E_{k}/n [GeV/n];Counts");
		createMcBkgHists(MC_BKGH2b, "H2b", "%s %s MC aboveL2(no L1 cut) frag Isotope Mass%d Counts vs E_{k}/n;%s E_{k}/n [GeV/n];Counts");
		createMcBkgHists(MC_BKGH2c, "H2c", "%s %s MC L1Orig-L2toL4 frag Isotope Mass%d Counts vs E_{k}/n;%s E_{k}/n [GeV/n];Counts");
		createMcBkgHists(MC_BKGH3a, "H3a", "%s %s MC L2 frag Isotope Mass%d Counts vs E_{k}/n;L2 Truth E_{k}/n [GeV/n];Counts");
		createMcBkgHists(MC_BKGH3b, "H3b", "%s %s MC L2 frag Isotope Mass%d Counts vs E_{k}/n;%s E_{k}/n [GeV/n];Counts");
		createMcBkgHists(MC_BKGH3c, "H3c", "%s %s MC L2toL4 frag Isotope Mass%d Counts vs E_{k}/n;%s E_{k}/n [GeV/n];Counts");
	}

	// ============ FLUX 区域 ============
	std::cout<<"DEBUG: Flux Hists"<<std::endl;
	
	// FLUXH1 - 这个结构复杂，保留原样
	FLUXH1.resize(Nchain);
	for (int c = 0; c < Nchain; ++c) {
		FLUXH1[c].resize(NcutGroups);
		for (int cg = 0; cg < NcutGroups; ++cg) {
			FLUXH1[c][cg].resize(NnumDen);
			for (int nd = 0; nd < NnumDen; ++nd) {
				FLUXH1[c][cg][nd].resize(Ndet);
				for (int d = 0; d < Ndet; ++d) {
					const int Ni = isISS ? Niso : 1;
					FLUXH1[c][cg][nd][d].resize(Ni);
					for (int i = 0; i < Ni; ++i) {
						int mass = isISS ? iso->getMass(i) : UseMass;
						auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
						FLUXH1[c][cg][nd][d][i] = createHist<TH1D>(
							Form("%s_FLUX_H1_%s_%s_%s_Mass%dBin",
								chains[c].c_str(), cut_groups[cg].c_str(),
								num_den[nd].c_str(), detectors[d].c_str(), mass),
							Form("%s %s %s %s Mass%dBin counts;%s E_{k}/n [GeV/n];Counts",
								chains[c].c_str(), detectors[d].c_str(),
								cut_groups[cg].c_str(), num_den[nd].c_str(),
								mass, detectors[d].c_str()),
							static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					}
				}
			}
		}
	}

	if (isISS) {
		ISS_FLUXH2.resize(1);
		ISS_FLUXH2[0] = createHist<TH1D>(
			"ISS_FLUX_H2",
			"ISS Exposure time;Rigidity [GV];Exposure Time [s]",
			static_cast<int>(rigBins.size()) - 1, rigBins.data());
		
		ISS_FLUXH3.resize(Ndet);
		for (int d = 0; d < Ndet; ++d) {
			ISS_FLUXH3[d].resize(Niso);
			for (int i = 0; i < Niso; ++i) {
				int mass = iso->getMass(i);
				auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
				ISS_FLUXH3[d][i] = createHist<TH1D>(
					Form("ISS_FLUX_H3_%s_Mass%dBin", detectors[d].c_str(), mass),
					Form("%s Mass%dBin Exposure time;E_{k}/n [GeV/n];Exposure Time [s]",
						detectors[d].c_str(), mass),
					static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
			}
		}
	} else {
		MC_FLUXH3.resize(1);
		auto ekBinsUse = binMgr.GetEkPerNucleonBins(charge, UseMass);
		MC_FLUXH3[0] = createHist<TH1D>(
			"MC_FLUX_H3",
			"MC Generated counts;E_{k}^{gen}/n [GeV/n];Counts",
			static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
		
		MC_FLUXH2.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			MC_FLUXH2[c].resize(NcutGroups);
			for (int cg = 0; cg < NcutGroups; ++cg) {
				MC_FLUXH2[c][cg].resize(Ndet);
				for (int d = 0; d < Ndet; ++d) {
					MC_FLUXH2[c][cg][d].resize(NgeneRec);
					for (int gr = 0; gr < NgeneRec; ++gr) {
						MC_FLUXH2[c][cg][d][gr] = createHist<TH1D>(
							Form("%s_MC_FLUX_H2_%s_%s_%s",
								chains[c].c_str(), cut_groups[cg].c_str(),
								detectors[d].c_str(), gene_rec[gr].c_str()),
							Form("%s %s MC %s %s counts;%s E_{k}/n [GeV/n];Counts",
								chains[c].c_str(), detectors[d].c_str(),
								cut_groups[cg].c_str(), gene_rec[gr].c_str(), detectors[d].c_str()),
							static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					}
				}
			}
		}
	}

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
    m_filteredTree->SetAutoSave(0);
    
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

    // Y-combinator 递归写入
    auto writeHists = [&](const auto& container) {
        auto impl = [](const auto& cont, auto& self) -> void {
            using T = std::decay_t<decltype(cont)>;
            if constexpr (std::is_same_v<T, std::vector<H1Ptr>> || 
                          std::is_same_v<T, std::vector<H2Ptr>>) {
                for (const auto& h : cont) if (h) h->Write();
            } else {
                for (const auto& item : cont) self(item, self);
            }
        };
        impl(container, impl);
    };

    // 简洁调用
    if (!ISS_IDH1.empty()) writeHists(ISS_IDH1);
    if (!MC_IDH1.empty()) writeHists(MC_IDH1);
    
	writeHists(IDH2); writeHists(IDH3);
    writeHists(IDH4a); writeHists(IDH4b);
    writeHists(IDH5a); writeHists(IDH5a2); writeHists(IDH5a3);
    writeHists(IDH5b); writeHists(IDH5b2); writeHists(IDH5b3);
    writeHists(IDH6a); writeHists(IDH6b); writeHists(IDH7a); writeHists(IDH7b);
    
    
    if (!ISS_BKGH1.empty()) {
        writeHists(ISS_BKGH1); writeHists(ISS_BKGH2);
        writeHists(ISS_BKGH3); writeHists(ISS_BKGH4);
    }
    if (!MC_BKGH1.empty()) {
        writeHists(MC_BKGH1); writeHists(MC_BKGH1b);writeHists(MC_BKGH2); writeHists(MC_BKGH2b); writeHists(MC_BKGH2c);
        writeHists(MC_BKGH3a); writeHists(MC_BKGH3b); writeHists(MC_BKGH3c);
    }
    
    //writeHists(FLUXH1);
    if (!ISS_FLUXH2.empty()) { writeHists(ISS_FLUXH2); writeHists(ISS_FLUXH3); }
    if (!MC_FLUXH2.empty()) { writeHists(MC_FLUXH3); /*writeHists(MC_FLUXH2)*/; }

    if (saveTree && m_filteredTree) {
        m_outputFile->cd();
        m_filteredTree->Write("", TObject::kOverwrite);
    }
    
    m_outputFile->Close();
}