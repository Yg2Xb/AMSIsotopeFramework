#include "HistManager.h"
#include "BinningManager.h"
#include <iostream>
#include <stdexcept>
#include "TString.h"
#include <tuple>
#include <vector>
#include <functional>

using namespace AMS_Iso;

// 简单常数快速判断碎裂产物元素（与 selectdata 保持一致的简化口径）
// 当 isBeFragment() 返回 true → Z=4(Be)；否则 → Z=5(B)
static constexpr bool FragZFast = false; // true: Be, false: B
static inline bool isBeFragment() { return FragZFast; }

static inline const std::vector<double>& safeBins(const std::string& key) {
	try {
		const auto& bins = BinningManager::GetInstance().Get(key);
		if (!bins.empty()) return bins;
	} catch (const std::exception& e) { // 放宽异常类型，兼容 runtime_error
		std::cerr << "Warning: Binning key '" << key << "' not found. Using default bins. Msg: " << e.what() << std::endl;
	} catch (...) {
		std::cerr << "Warning: Unknown error getting bins for key '" << key << "'. Using default bins." << std::endl;
	}
	static const std::vector<double> default_bins = {0, 1};
	return default_bins;
}

HistManager::HistManager(const std::string& output_filename,
		bool isISS,
		const std::vector<std::string>& chains,
		int charge,
		const IsotopeVar* iso,
		int UseMass) {
	m_outputFile = std::make_unique<TFile>(output_filename.c_str(), "RECREATE");
	if (!m_outputFile || m_outputFile->IsZombie()) {
		throw std::runtime_error("Failed to create output ROOT file: " + output_filename);
	}
	std::cout << "HistManager: Output file '" << output_filename << "' opened." << std::endl;

	auto& binMgr = BinningManager::GetInstance();
	const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
	const std::vector<std::string> cut_groups = {
		"BasicAndFiducial", "Trigger", "InnerTracker", "L1BigZ", "L1PickUp", "L1UpperQ", "UpperTOFQ", "BkgReduction",
		"TOFGeo", "TOFBetaQuality", "NaFGeo", "NaFReconstruction", "AGLGeo", "AGLReconstruction"
	};
	const std::vector<std::string> num_den = {"Num", "Den"};
	const std::vector<std::string> charge_types = {"L1QSignal", "L1QTemplate", "L2QTemplate"};
	const std::vector<std::string> sources = {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
	const std::vector<std::string> gene_rec = {"Gene", "Rec"};

	// 动态的碎裂产物配置：Be 或 B
	const bool beFrag = isBeFragment();
	const int fragZ = beFrag ? 4 : 5;
	const std::vector<int> FragA = beFrag ? std::vector<int>{7, 9, 10} : std::vector<int>{10, 11};

	int Nchain = static_cast<int>(chains.size());
	int Ndet = static_cast<int>(detectors.size());
	int Niso = iso->getIsotopeCount(); // ISS: 输入元素的同位素数；MC: 1（用于 IDH2/3 等保持结构）
	int Nsrc = static_cast<int>(sources.size());
	int Nct = static_cast<int>(charge_types.size());
	int NcutGroups = static_cast<int>(cut_groups.size());
	int NnumDen = static_cast<int>(num_den.size());
	int NgeneRec = static_cast<int>(gene_rec.size());
	
	auto ekBinsStd = safeBins("EkPerNucleon");
	auto rigBins = safeBins("Rigidity");

	std::cout<<"DEBUG: Niso="<<Niso<<" charge="<<charge<<std::endl;
	std::cout<<"DEBUG: iso ptr="<<iso<<std::endl;

	// ---------------- ID 区域 ----------------
	std::cout<<"DEBUG: ID Hists"<<std::endl;
	if (isISS) {
		// ISS_IDH1: 每个链、每个探测器、每个输入元素的同位素，使用该(Z,A)专属的 Ek/n 分箱
		ISS_IDH1.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			ISS_IDH1[c].resize(Ndet);
			for (int d = 0; d < Ndet; ++d) {
				ISS_IDH1[c][d].resize(Niso);
				for (int i = 0; i < Niso; ++i) {
					int mass = iso->getMass(i);
					auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
					ISS_IDH1[c][d][i] = createHist<TH1F>(
							Form("%s_ISS_ID_H1_%s_Mass%dBin", chains[c].c_str(), detectors[d].c_str(), mass),
							Form("%s %s Mass%dBin isotope counts;Counts;%s E_{k}/n [GeV/n]",
								chains[c].c_str(), detectors[d].c_str(), mass, detectors[d].c_str()),
							static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				}
			}
		}
	} else {
		// MC_IDH1: 扩展为 Niso 维度（对输入元素的每个同位素 A 分别建模板，X=1/Mass，Y=Ek/n(A,Z)）
		int NisoMC = iso->getIsotopeCount();
		MC_IDH1.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			MC_IDH1[c].resize(Ndet);
			for (int d = 0; d < Ndet; ++d) {
				MC_IDH1[c][d].resize(NisoMC);
				for (int i = 0; i < NisoMC; ++i) {
					int mass_i = iso->getMass(i);
					auto ekBins_i = binMgr.GetEkPerNucleonBins(charge, mass_i);
					MC_IDH1[c][d][i] = createHist<TH2F>(
							Form("%s_MC_ID_H1_%s_Mass%dBin", chains[c].c_str(), detectors[d].c_str(), mass_i),
							Form("%s %s Mass%dBin Isotope MC NoFragCut 1/Mass template;%s 1/Mass;%s E_{k}/n [GeV/n]",
								chains[c].c_str(), detectors[d].c_str(), mass_i, detectors[d].c_str(), detectors[d].c_str()),
							200, 0, 0.5, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
				}
			}
		}
	}

	IDH2.resize(Nchain);
	IDH3.resize(Nchain);
	IDH4a.resize(Nchain);
	IDH4b.resize(Nchain);
	IDH5a.resize(Nchain);
	IDH5b.resize(Nchain);
	IDH6a.resize(Nchain);
	IDH6b.resize(Nchain);
	IDH7a.resize(Nchain);
	IDH7b.resize(Nchain);
	for (int c = 0; c < Nchain; ++c) {
		IDH2[c].resize(Ndet);
		IDH3[c].resize(Ndet);
		for (int d = 0; d < Ndet; ++d) {
			IDH2[c][d].resize(Niso);
			for (int i = 0; i < Niso; ++i) {
				int mass = iso->getMass(i);
				auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
				IDH2[c][d][i] = createHist<TH2F>(
						Form("%s_ID_H2_%s_Mass%dBin", chains[c].c_str(), detectors[d].c_str(), mass),
						Form("%s %s UseMass%dBin Isotope MC 1/Mass vs E_{k}/n;%s 1/Mass;%s E_{k}/n [GeV/n]",
							chains[c].c_str(), detectors[d].c_str(), mass, detectors[d].c_str(), detectors[d].c_str()),
						200, 0, 0.5, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
			}
			int mass_heaviest = iso->getMass(Niso-1); // Heaviest isotope
			auto ekBinsH = binMgr.GetEkPerNucleonBins(charge, mass_heaviest);
			IDH3[c][d] = createHist<TH2F>(
					Form("%s_ID_H3_%s", chains[c].c_str(), detectors[d].c_str()),
					Form("%s %s Heaviest iso 1/Mass;%s 1/Mass;%s E_{k}/n [GeV/n]",
						chains[c].c_str(), detectors[d].c_str(), detectors[d].c_str(), detectors[d].c_str()),
					200, 0, 0.5, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
		}
		// RICH
		int idx = charge - 1;
		IDH4a[c] = createHist<TH1F>(
				Form("%s_ID_H4a", chains[c].c_str()),
				Form("%s NaF 1/#beta (beta~1);1/#beta (NaF);Counts", chains[c].c_str()),
				Detector::RichBins[0][idx], 1.0 - Detector::RichAxis[0], 1.0 + Detector::RichAxis[0]);
		IDH4b[c] = createHist<TH1F>(
				Form("%s_ID_H4b", chains[c].c_str()),
				Form("%s AGL 1/#beta (beta~1);1/#beta (AGL);Counts", chains[c].c_str()),
				Detector::RichBins[1][idx], 1.0 - Detector::RichAxis[1], 1.0 + Detector::RichAxis[1]);

		auto betarigBins = binMgr.GetBetaRigBins(charge, iso->getMass(Niso-1));

		// Δβ: X=Δβ, Y=rig/ek/bg
		IDH5a[c] = createHist<TH2F>(
				Form("%s_ID_H5a", chains[c].c_str()),
				"NaF-Tracker #Delta(1/#beta);NaF-Tracker #Delta(1/#beta);Rigidity [GV]",
				400, -0.2, 0.2, static_cast<int>(rigBins.size()) - 1, rigBins.data());
		IDH5b[c] = createHist<TH2F>(
				Form("%s_ID_H5b", chains[c].c_str()),
				"AGL-Tracker #Delta(1/#beta);AGL-Tracker #Delta(1/#beta);Rigidity [GV]",
				400, -0.2, 0.2, static_cast<int>(rigBins.size()) - 1, rigBins.data());
		IDH6a[c] = createHist<TH2F>(
				Form("%s_ID_H6a", chains[c].c_str()),
				"TOF-NaF #Delta(1/#beta);#TOF-NaF Delta(1/#beta);NaF E_{k}/n [GeV/n]",
				400, -0.2, 0.2, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
		IDH6b[c] = createHist<TH2F>(
				Form("%s_ID_H6b", chains[c].c_str()),
				"TOF-AGL #Delta(1/#beta);TOF-AGL #Delta(1/#beta);AGL E_{k}/n [GeV/n]",
				400, -0.2, 0.2, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
		IDH7a[c] = createHist<TH2F>(
				Form("%s_ID_H7a", chains[c].c_str()),
				"TOF-NaF #Delta(1/#beta);TOF-NaF #Delta(1/#beta);NaF #beta R[GV]",
				400, -0.2, 0.2, static_cast<int>(betarigBins.size()) - 1, betarigBins.data());
		IDH7b[c] = createHist<TH2F>(
				Form("%s_ID_H7b", chains[c].c_str()),
				"TOF-AGL #Delta(1/#beta);TOF-AGL #Delta(1/#beta);AGL #beta R[GV]",
				400, -0.2, 0.2, static_cast<int>(betarigBins.size()) - 1, betarigBins.data());
	}

	// ---------------- BKG 区域 ----------------
	//Nfragiso dimension for BKG histograms (Be7, Be9, Be10) or (B10, B11), BUT union EKbin!!!
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
				int zsrc = s + 4; // Z=4..8
				ISS_BKGH1[c][s].resize(Ndet);
				ISS_BKGH2[c][s].resize(Ndet);
				ISS_BKGH3[c][s].resize(Ndet);
				ISS_BKGH4[c][s].resize(Ndet);
				for (int d = 0; d < Ndet; ++d) {
					// H1：源计数（仍使用标准 Ek 分箱，作为总体参考）
					ISS_BKGH1[c][s][d] = createHist<TH1F>(
							Form("%s_ISS_BKG_H1_%s_%s", chains[c].c_str(), sources[s].c_str(), detectors[d].c_str()),
							Form("%s %s L1 Source %s counts;%s E_{k}/n [GeV/n];Counts",
								chains[c].c_str(), detectors[d].c_str(), sources[s].c_str(), detectors[d].c_str()),
							static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

					// H2：不同 charge_type 的 charge vs Ek/n（仍为总体参考）
					ISS_BKGH2[c][s][d].resize(Nct);
					for (int t = 0; t < Nct; ++t) {
						const char* x_title = "TrackerLayer Charge";
						int x_bins = 400;
						double x_min = zsrc - 2;
						double x_max = zsrc + 2;
						if (charge_types[t] == "L1QSignal") {
							x_bins = 600;
							x_min = 3;
							x_max = 9;
						}
						ISS_BKGH2[c][s][d][t] = createHist<TH2F>(
								Form("%s_ISS_BKG_H2_%s_%s_%s", chains[c].c_str(), sources[s].c_str(), charge_types[t].c_str(), detectors[d].c_str()),
								Form("%s %s %s %s charge vs E_{k}/n;%s;%s E_{k}/n [GeV/n]",
									chains[c].c_str(), detectors[d].c_str(), charge_types[t].c_str(), sources[s].c_str(), x_title, detectors[d].c_str()),
								x_bins, x_min, x_max, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					}

					// H3：碎裂计数（总体参考）
					ISS_BKGH3[c][s][d] = createHist<TH1F>(
							Form("%s_ISS_BKG_H3_%s_%s", chains[c].c_str(), sources[s].c_str(), detectors[d].c_str()),
							Form("%s %s L1%s L2 frag counts;%s E_{k}/n [GeV/n];Counts",
								chains[c].c_str(), detectors[d].c_str(), sources[s].c_str(), detectors[d].c_str()),
							static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());

					// H4：碎裂 1/Mass vs Ek/n，按碎裂产物（Be 或 B）对应的 A 列表维度展开，Ek 分箱按 (fragZ, A)
					ISS_BKGH4[c][s][d].resize(static_cast<int>(FragA.size()));
					for (size_t bi = 0; bi < FragA.size(); ++bi) {
						int A = FragA[bi];
						auto ekBinsFrag = binMgr.GetEkPerNucleonBins(fragZ, A);
						ISS_BKGH4[c][s][d][bi] = createHist<TH2F>(
								Form("%s_ISS_BKG_H4_%s_%s_Z%d_Mass%d", chains[c].c_str(), sources[s].c_str(), detectors[d].c_str(), fragZ, A),
								Form("%s %s L1%s L2frag 1/Mass vs E_{k}/n (Z=%d Mass%d);%s 1/Mass;%s E_{k}/n [GeV/n]",
									chains[c].c_str(), detectors[d].c_str(), sources[s].c_str(), fragZ, A, detectors[d].c_str(), detectors[d].c_str()),
								200, 0, 0.5, static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					}
				}
			}
		}
	} else { // MC
		int NisoBKG = static_cast<int>(FragA.size());

		MC_BKGH1.resize(Nchain);
		for (int c = 0; c < Nchain; ++c) {
			MC_BKGH1[c].resize(Ndet);
			for (int d = 0; d < Ndet; ++d) {
				MC_BKGH1[c][d] = createHist<TH1F>(
						Form("%s_MC_BKG_H1_%s", chains[c].c_str(), detectors[d].c_str()),
						Form("%s %s MC input counts;%s E_{k}/n [GeV/n];Counts",
							chains[c].c_str(), detectors[d].c_str(), detectors[d].c_str()),
						static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
			}
		}

		auto createMcBkgHists = [&](auto& container, const std::string& prefix, const std::string& titleFmt) {
			container.resize(Nchain);
			for (int c = 0; c < Nchain; ++c) {
				container[c].resize(Ndet);
				for (int d = 0; d < Ndet; ++d) {
					container[c][d].resize(NisoBKG);
					for (int i = 0; i < NisoBKG; ++i) {
						int massA = FragA[i];
						auto ekBinsFrag = binMgr.GetEkPerNucleonBins(fragZ, massA);
						container[c][d][i] = createHist<TH1F>(
								Form("%s_MC_BKG_%s_%s_Z%d_Mass%d", chains[c].c_str(), prefix.c_str(), detectors[d].c_str(), fragZ, massA),
								Form(titleFmt.c_str(), chains[c].c_str(), detectors[d].c_str(), massA, detectors[d].c_str()),
								static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					}
				}
			}
		};

		// 标题模板将 Mass%d 放到第三个占位符（与上方 Form 调用保持一致）
		createMcBkgHists(MC_BKGH2, "H2", "%s %s MC frag Isotope Mass%d Counts vs E_{k}/n;%s E_{k}/n [GeV/n];Counts");
		createMcBkgHists(MC_BKGH3a, "H3a", "%s %s MC upTOF frag Isotope Mass%d Counts vs E_{k}/n;Generated E_{k}/n [GeV/n];Counts");
		createMcBkgHists(MC_BKGH3b, "H3b", "%s %s MC upTOF frag survival in rich Isotope Mass%d Counts vs E_{k}/n;Generated E_{k}/n [GeV/n];Counts");
	}

	// ---------------- FLUX 区域 ----------------
	std::cout<<"DEBUG: Flux Hists"<<std::endl;
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

						std::string hname = Form("%s_FLUX_H1_%s_%s_%s_Mass%dBin",
								(c < (int)chains.size() ? chains[c].c_str() : "<bad>"),
								(cg < (int)cut_groups.size() ? cut_groups[cg].c_str() : "<bad>"),
								(nd < (int)num_den.size() ? num_den[nd].c_str() : "<bad>"),
								(d < (int)detectors.size() ? detectors[d].c_str() : "<bad>"),
								mass);
						std::string htitle = Form("%s %s %s %s Mass%dBin counts;%s E_{k}/n [GeV/n];Counts",
								(c < (int)chains.size() ? chains[c].c_str() : "<bad>"),
								(d < (int)detectors.size() ? detectors[d].c_str() : "<bad>"),
								(cg < (int)cut_groups.size() ? cut_groups[cg].c_str() : "<bad>"),
								(nd < (int)num_den.size() ? num_den[nd].c_str() : "<bad>"),
								mass,
								(d < (int)detectors.size() ? detectors[d].c_str() : "<bad>"));
						// 真正创建
						FLUXH1[c][cg][nd][d][i] = createHist<TH1F>(
								hname.c_str(),
								htitle.c_str(),
								static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
					}
				}
			}
		}
	}

	if (isISS) {
		ISS_FLUXH2.resize(1);
		auto rigBins = safeBins("Rigidity");
		ISS_FLUXH2[0] = createHist<TH1F>(
				"ISS_FLUX_H2", "ISS Exposure time;Rigidity [GV];Exposure Time [s]", static_cast<int>(rigBins.size()) - 1, rigBins.data());
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
	} else {
		MC_FLUXH2.resize(Nchain);
		MC_FLUXH3.resize(1);
		auto ekBinsUse = binMgr.GetEkPerNucleonBins(charge, UseMass);
		MC_FLUXH3[0] = createHist<TH1F>(//std bin for bkg analysis
				"MC_FLUX_H3",
				"MC Generated counts;E_{k}^{gen}/n [GeV/n];Counts",
				static_cast<int>(ekBinsStd.size()) - 1, ekBinsStd.data());
		for (int c = 0; c < Nchain; ++c) {
			MC_FLUXH2[c].resize(NcutGroups);
			for (int cg = 0; cg < NcutGroups; ++cg) {
				MC_FLUXH2[c][cg].resize(Ndet);
				for (int d = 0; d < Ndet; ++d) {
					MC_FLUXH2[c][cg][d].resize(NgeneRec);
					for (int gr = 0; gr < NgeneRec; ++gr) {
						MC_FLUXH2[c][cg][d][gr] = createHist<TH1F>(
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
	std::cout<<"DEBUG Finish Hist Defination="<<std::endl;
}

void HistManager::Save() {
	if (!m_outputFile || !m_outputFile->IsOpen()) {
		std::cerr << "HistManager: Error: Output file is not open or is invalid." << std::endl;
		return;
	}
	m_outputFile->cd();

	// Helper functions for writing different histogram types
	auto writeH1 = [](const H1Ptr& hist) {
		if (hist) hist->Write();
	};

	auto writeH2 = [](const H2Ptr& hist) {
		if (hist) hist->Write();
	};

	// Separate functions for different container depths to avoid recursion issues
	auto write1D = [&](const auto& container1D) {
		for (const auto& item : container1D) {
			if constexpr (std::is_same_v<std::decay_t<decltype(item)>, H1Ptr>) {
				writeH1(item);
			} else if constexpr (std::is_same_v<std::decay_t<decltype(item)>, H2Ptr>) {
				writeH2(item);
			}
		}
	};

	auto write2D = [&](const auto& container2D) {
		for (const auto& inner : container2D) {
			write1D(inner);
		}
	};

	auto write3D = [&](const auto& container3D) {
		for (const auto& inner : container3D) {
			write2D(inner);
		}
	};

	auto write4D = [&](const auto& container4D) {
		for (const auto& inner : container4D) {
			write3D(inner);
		}
	};

	auto write5D = [&](const auto& container5D) {
		for (const auto& inner : container5D) {
			write4D(inner);
		}
	};

	// ID histograms
	write3D(IDH2);
	write2D(IDH3);

    /*
	write1D(IDH4a);
	write1D(IDH4b);
	write1D(IDH5a);
	write1D(IDH5b);
	write1D(IDH6a);
	write1D(IDH6b);
    */
	write1D(IDH7a);
	write1D(IDH7b);

	if (!ISS_IDH1.empty()) write3D(ISS_IDH1);
	if (!MC_IDH1.empty()) write3D(MC_IDH1);

	// BKG histograms
	if (!ISS_BKGH1.empty()) {
		write3D(ISS_BKGH1);
		write4D(ISS_BKGH2);
		write3D(ISS_BKGH3);
		write4D(ISS_BKGH4); // 由原 3D 改为 4D（最后一维为 Frag 同位素）
	}
	if (!MC_BKGH1.empty()) {
		write2D(MC_BKGH1);
		write3D(MC_BKGH2);
		write3D(MC_BKGH3a);
		write3D(MC_BKGH3b);
	}

	// FLUX histograms
	//write5D(FLUXH1);

	if (!ISS_FLUXH2.empty()) {
		write1D(ISS_FLUXH2);
		write2D(ISS_FLUXH3);
	}
	if (!MC_FLUXH2.empty()) {
		//write4D(MC_FLUXH2);
		write1D(MC_FLUXH3);
	}

	std::cout << "HistManager: all histograms saved." << std::endl;
	m_outputFile->Close(); 
}
