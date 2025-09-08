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
    } catch (const std::out_of_range& e) {
        std::cerr << "Warning: Binning key '" << key << "' not found. Using default bins." << std::endl;
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

    int Nchain = chains.size();
    int Ndet = detectors.size();
    int Niso = isISS ? iso->getIsotopeCount() : 1;
    int Nsrc = sources.size();
    int Nct = charge_types.size();
    int NcutGroups = cut_groups.size();
    int NnumDen = num_den.size();
    int NgeneRec = gene_rec.size();
    std::cout<<"DEBUG: Niso="<<Niso<<" charge="<<charge<<std::endl;
    std::cout<<"DEBUG: iso ptr="<<iso<<std::endl;
    
    // ---------------- ID 区域 ----------------
    std::cout<<"DEBUG: ID Hists"<<std::endl;
    if (isISS) {
        ISS_IDH1.resize(Nchain);
        for (int c = 0; c < Nchain; ++c) {
            ISS_IDH1[c].resize(Ndet);
            for (int d = 0; d < Ndet; ++d) {
                ISS_IDH1[c][d].resize(Niso);
                for (int i = 0; i < Niso; ++i) {
                    int mass = iso->getMass(i);
                    auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
                    ISS_IDH1[c][d][i] = createHist<TH1F>(
                        Form("%s_ISS_ID_H1_%s_Mass%d", chains[c].c_str(), detectors[d].c_str(), mass),
                        Form("%s %s Mass%d isotope counts;Counts;%s E_{k}/n [GeV/n]",
                             chains[c].c_str(), detectors[d].c_str(), mass, detectors[d].c_str()),
                        ekBins.size() - 1, ekBins.data());
                }
            }
        }
    } else {
        MC_IDH1.resize(Nchain);
        for (int c = 0; c < Nchain; ++c) {
            MC_IDH1[c].resize(Ndet);
            for (int d = 0; d < Ndet; ++d) {
                auto ekBins = binMgr.GetEkPerNucleonBins(charge, UseMass);
                MC_IDH1[c][d] = createHist<TH2F>(
                    Form("%s_MC_ID_H1_%s", chains[c].c_str(), detectors[d].c_str()),
                    Form("%s %s Mass%dIsotope MC 1/Mass template;%s 1/Mass;%s E_{k}/n [GeV/n]",
                         chains[c].c_str(), detectors[d].c_str(), UseMass, detectors[d].c_str(), detectors[d].c_str()),
                    200, 0, 0.5, ekBins.size() - 1, ekBins.data());
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
            for (int i = 0; i < Niso; ++i) {
                int mass = isISS ? iso->getMass(i) : UseMass;
                auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
                IDH2[c][d].resize(Niso);
                IDH2[c][d][i] = createHist<TH2F>(
                    Form("%s_ID_H2_%s_%d", chains[c].c_str(), detectors[d].c_str(), mass),
                    Form("%s %s UseMass%d 1/Mass vs E_{k}/n;%s 1/Mass;%s E_{k}/n [GeV/n]",
                         chains[c].c_str(), detectors[d].c_str(), mass, detectors[d].c_str(), detectors[d].c_str()),
                    200, 0, 0.5, ekBins.size() - 1, ekBins.data());
            }
            int mass = iso->getMass(Niso-1); // Heaviest isotope
            auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
            IDH3[c][d] = createHist<TH2F>(
                Form("%s_ID_H3_%s", chains[c].c_str(), detectors[d].c_str()),
                Form("%s %s Heaviest iso 1/Mass;%s 1/Mass;%s E_{k}/n [GeV/n]",
                     chains[c].c_str(), detectors[d].c_str(), detectors[d].c_str(), detectors[d].c_str()),
                200, 0, 0.5, ekBins.size() - 1, ekBins.data());
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

        auto rigBins = safeBins("Rigidity");
        auto ekBins = safeBins("EkPerNucleon");
        auto betagammaBins = binMgr.GetBetaGammaBins(charge, UseMass);

        // Δβ: X=Δβ, Y=rig/ek/bg
        IDH5a[c] = createHist<TH2F>(
            Form("%s_ID_H5a", chains[c].c_str()),
            "NaF-Tracker #Delta(1/#beta);NaF-Tracker #Delta(1/#beta);Rigidity [GV]",
            400, -0.2, 0.2, rigBins.size() - 1, rigBins.data());
        IDH5b[c] = createHist<TH2F>(
            Form("%s_ID_H5b", chains[c].c_str()),
            "AGL-Tracker #Delta(1/#beta);AGL-Tracker #Delta(1/#beta);Rigidity [GV]",
            400, -0.2, 0.2, rigBins.size() - 1, rigBins.data());
        IDH6a[c] = createHist<TH2F>(
            Form("%s_ID_H6a", chains[c].c_str()),
            "TOF-NaF #Delta(1/#beta);#TOF-NaF Delta(1/#beta);NaF E_{k}/n [GeV/n]",
            400, -0.2, 0.2, ekBins.size() - 1, ekBins.data());
        IDH6b[c] = createHist<TH2F>(
            Form("%s_ID_H6b", chains[c].c_str()),
            "TOF-AGL #Delta(1/#beta);TOF-AGL #Delta(1/#beta);AGL E_{k}/n [GeV/n]",
            400, -0.2, 0.2, ekBins.size() - 1, ekBins.data());
        IDH7a[c] = createHist<TH2F>(
            Form("%s_ID_H7a", chains[c].c_str()),
            "TOF-NaF #Delta(1/#beta);TOF-NaF #Delta(1/#beta);NaF #beta#gamma",
            400, -0.2, 0.2, betagammaBins.size() - 1, betagammaBins.data());
        IDH7b[c] = createHist<TH2F>(
            Form("%s_ID_H7b", chains[c].c_str()),
            "TOF-AGL #Delta(1/#beta);TOF-AGL #Delta(1/#beta);AGL #beta#gamma",
            400, -0.2, 0.2, betagammaBins.size() - 1, betagammaBins.data());
    }

    // ---------------- BKG 区域 ----------------
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
                ISS_BKGH1[c][s].resize(Ndet);
                ISS_BKGH2[c][s].resize(Ndet);
                ISS_BKGH3[c][s].resize(Ndet);
                ISS_BKGH4[c][s].resize(Ndet);
                for (int d = 0; d < Ndet; ++d) {
                    auto ekBins = safeBins("EkPerNucleon");
                    ISS_BKGH1[c][s][d] = createHist<TH1F>(
                        Form("%s_ISS_BKG_H1_%s_%s", chains[c].c_str(), sources[s].c_str(), detectors[d].c_str()),
                        Form("%s %s L1 Source %s counts;%s E_{k}/n [GeV/n];Counts",
                             chains[c].c_str(), detectors[d].c_str(), sources[s].c_str(), detectors[d].c_str()),
                        ekBins.size() - 1, ekBins.data());
                    ISS_BKGH2[c][s][d].resize(Nct);
                    for (int t = 0; t < Nct; ++t) {
                        const char* x_title = "TrackerLayer Charge";
                        int x_bins = 400;
                        double x_min = charge - 2;
                        double x_max = charge + 2;
                        if (charge_types[t] == "L1QSignal") {
                            x_bins = 600;
                            x_min = 3;
                            x_max = 9;
                        }
                        ISS_BKGH2[c][s][d][t] = createHist<TH2F>(
                            Form("%s_ISS_BKG_H2_%s_%s_%s", chains[c].c_str(), sources[s].c_str(), charge_types[t].c_str(), detectors[d].c_str()),
                            Form("%s %s %s %s charge vs E_{k}/n;%s;%s E_{k}/n [GeV/n]",
                                 chains[c].c_str(), detectors[d].c_str(), charge_types[t].c_str(), sources[s].c_str(), x_title, detectors[d].c_str()),
                            x_bins, x_min, x_max, ekBins.size() - 1, ekBins.data());
                    }
                    ISS_BKGH3[c][s][d] = createHist<TH1F>(
                        Form("%s_ISS_BKG_H3_%s_%s", chains[c].c_str(), sources[s].c_str(), detectors[d].c_str()),
                        Form("%s %s L1%s L2 frag counts;%s E_{k}/n [GeV/n];Counts",
                             chains[c].c_str(), detectors[d].c_str(), sources[s].c_str(), detectors[d].c_str()),
                        ekBins.size() - 1, ekBins.data());
                    ISS_BKGH4[c][s][d] = createHist<TH2F>(
                        Form("%s_ISS_BKG_H4_%s_%s", chains[c].c_str(), sources[s].c_str(), detectors[d].c_str()),
                        Form("%s %s L1%s L2frag 1/Mass vs E_{k}/n;%s 1/Mass;%s E_{k}/n [GeV/n]",
                             chains[c].c_str(), detectors[d].c_str(), sources[s].c_str(), detectors[d].c_str(), detectors[d].c_str()),
                        200, 0, 0.5, ekBins.size() - 1, ekBins.data());
                }
            }
        }
    } else { // MC
        int NisoBKG = iso->getIsotopeCount();
        auto ekBins = safeBins("EkPerNucleon");

        MC_BKGH1.resize(Nchain);
        for (int c = 0; c < Nchain; ++c) {
            MC_BKGH1[c].resize(Ndet);
            for (int d = 0; d < Ndet; ++d) {
                 MC_BKGH1[c][d] = createHist<TH1F>(
                    Form("%s_MC_BKG_H1_%s", chains[c].c_str(), detectors[d].c_str()),
                    Form("%s %s MC input counts;%s E_{k}/n [GeV/n];Counts",
                         chains[c].c_str(), detectors[d].c_str(), detectors[d].c_str()),
                    ekBins.size() - 1, ekBins.data());
            }
        }
        
        auto createMcBkgHists = [&](auto& container, const std::string& prefix, const std::string& titleFmt) {
            container.resize(Nchain);
            for (int c = 0; c < Nchain; ++c) {
                container[c].resize(Ndet);
                for (int d = 0; d < Ndet; ++d) {
                    container[c][d].resize(NisoBKG);
                    for (int i = 0; i < NisoBKG; ++i) {
                        int mass = iso->getMass(i);
                        container[c][d][i] = createHist<TH1F>(
                            Form("%s_MC_BKG_%s_%s_Mass%d", chains[c].c_str(), prefix.c_str(), detectors[d].c_str(), mass),
                            Form(titleFmt.c_str(), chains[c].c_str(), detectors[d].c_str(), mass, detectors[d].c_str()),
                            ekBins.size() - 1, ekBins.data());
                    }
                }
            }
        };
        
        createMcBkgHists(MC_BKGH2, "H2", "%s %s MC frag Isotope Mass%d Counts vs E_{k}/n;%s E_{k}/n [GeV/n];Counts");
        createMcBkgHists(MC_BKGH3a, "H3a", "%s %s MC upTOF frag Isotope Mass%d Counts vs E_{k}/n;%s E_{k}/n [GeV/n];Counts");
        createMcBkgHists(MC_BKGH3b, "H3b", "%s %s MC upTOF frag survival in rich Isotope Mass%d Counts vs E_{k}/n;%s E_{k}/n [GeV/n];Counts");
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
                    FLUXH1[c][cg][nd][d].resize(Niso);
                    for (int i = 0; i < Niso; ++i) {
                        int mass = isISS ? iso->getMass(i) : UseMass;
                        auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
                        FLUXH1[c][cg][nd][d][i] = createHist<TH1F>(
                            Form("%s_FLUX_H1_%s_%s_%s_%d",
                                 chains[c].c_str(), cut_groups[cg].c_str(),
                                 num_den[nd].c_str(), detectors[d].c_str(), mass),
                            Form("%s %s %s %s Mass%d counts;%s E_{k}/n [GeV/n];Counts",
                                 chains[c].c_str(), detectors[d].c_str(),
                                 cut_groups[cg].c_str(), num_den[nd].c_str(), mass, detectors[d].c_str()),
                            ekBins.size() - 1, ekBins.data());
                    }
                }
            }
        }
    }

    if (isISS) {
        ISS_FLUXH2.resize(1);
        auto rigBins = safeBins("Rigidity");
        ISS_FLUXH2[0] = createHist<TH1F>(
            "ISS_FLUX_H2", "ISS Exposure time;Rigidity [GV];Exposure Time [s]", rigBins.size() - 1, rigBins.data());
        ISS_FLUXH3.resize(Nchain);
        for (int c = 0; c < Nchain; ++c) {
            ISS_FLUXH3[c].resize(Ndet);
            for (int d = 0; d < Ndet; ++d) {
                ISS_FLUXH3[c][d].resize(Niso);
                for (int i = 0; i < Niso; ++i) {
                    int mass = iso->getMass(i);
                    auto ekBins = binMgr.GetEkPerNucleonBins(charge, mass);
                    ISS_FLUXH3[c][d][i] = createHist<TH1F>(
                        Form("%s_ISS_FLUX_H3_%s_Mass%d", chains[c].c_str(), detectors[d].c_str(), mass),
                        Form("%s %s Exposure time;E_{k}/n [GeV/n];Exposure Time [s]",
                             chains[c].c_str(), detectors[d].c_str()),
                        ekBins.size() - 1, ekBins.data());
                }
            }
        }
    } else {
        MC_FLUXH2.resize(Nchain);
        MC_FLUXH3.resize(1);
        auto ekBins = binMgr.GetEkPerNucleonBins(charge, UseMass);
        MC_FLUXH3[0] = createHist<TH1F>(
            "MC_FLUX_H3",
            "MC Generated counts;E_{k}^{gen}/n [GeV/n];Counts",
            ekBins.size() - 1, ekBins.data());
        auto ekBins_flux = safeBins("EkPerNucleon");
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
                            ekBins_flux.size() - 1, ekBins_flux.data());
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
    
    write1D(IDH4a);
    write1D(IDH4b);
    write1D(IDH5a);
    write1D(IDH5b);
    write1D(IDH6a);
    write1D(IDH6b);
    write1D(IDH7a);
    write1D(IDH7b);
    
    if (!ISS_IDH1.empty()) write3D(ISS_IDH1);
    if (!MC_IDH1.empty()) write2D(MC_IDH1);

    // BKG histograms
    if (!ISS_BKGH1.empty()) {
        write3D(ISS_BKGH1);
        write4D(ISS_BKGH2);
        write3D(ISS_BKGH3);
        write3D(ISS_BKGH4);
    }
    if (!MC_BKGH1.empty()) {
        write2D(MC_BKGH1);
        write3D(MC_BKGH2);
        write3D(MC_BKGH3a);
        write3D(MC_BKGH3b);
    }

    // FLUX histograms
    write5D(FLUXH1);
    
    if (!ISS_FLUXH2.empty()) {
        write1D(ISS_FLUXH2);
        write3D(ISS_FLUXH3);
    }
    if (!MC_FLUXH2.empty()) {
        write4D(MC_FLUXH2);
        write1D(MC_FLUXH3);
    }

    std::cout << "HistManager: all histograms saved." << std::endl;
}
