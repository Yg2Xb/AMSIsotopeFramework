#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <algorithm>
#include <ctime>
#include <bitset>
#include <sys/resource.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TVectorD.h>
#include <TString.h>
#include "./EventProcessor.hh"

using namespace AMS_Iso;

const std::string lookupFileName = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/CDFLookupTable_fromSpline_0p8.root";
const std::string dataFileName  = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Bor_BeToC.root";
const std::string output_base_path = "/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/";

// Manages a lookup table from a ROOT file, with validity checks.
struct LookupTable {
    std::vector<double> q_values, cdf_l1, cdf_l2;
    double q_min = 0.0, q_max = 0.0, dq = 0.0;
    bool isValid = false;

    void loadFromFile(TFile* f, const std::string& baseName, const std::string& type) {
        isValid = false; 

        auto get_vec = [&](const std::string& name) {
            return dynamic_cast<TVectorD*>(f->Get(name.c_str()));
        };

        TVectorD* q_vec = get_vec(baseName + "_" + type + "_q");
        TVectorD* cdf1_vec = get_vec(baseName + "_" + type + "_cdf_l1");
        TVectorD* cdf2_vec = get_vec(baseName + "_" + type + "_cdf_l2");
        
        if (q_vec && cdf1_vec && cdf2_vec && q_vec->GetNrows() > 1) {
            q_values.assign(q_vec->GetMatrixArray(), q_vec->GetMatrixArray() + q_vec->GetNrows());
            cdf_l1.assign(cdf1_vec->GetMatrixArray(), cdf1_vec->GetMatrixArray() + cdf1_vec->GetNrows());
            cdf_l2.assign(cdf2_vec->GetMatrixArray(), cdf2_vec->GetMatrixArray() + cdf2_vec->GetNrows());
            
            q_min = q_values.front();
            q_max = q_values.back();
            dq = (q_values.size() > 1) ? (q_max - q_min) / (q_values.size() - 1) : 0.0;
            
            if (std::abs(dq) > 1e-9) {
                isValid = true;
            }
        }
    }
    
    double getCDF_L2(double q) const {
        if (!isValid || q < q_min || q > q_max) return -1.0;
        double fidx = (q - q_min) / dq;
        int idx = static_cast<int>(fidx);
        if (idx >= static_cast<int>(cdf_l2.size()) - 1) return cdf_l2.back();
        if (idx < 0) return cdf_l2.front();
        double frac = fidx - idx;
        return cdf_l2[idx] + frac * (cdf_l2[idx+1] - cdf_l2[idx]);
    }
    
    double getInvCDF_L1(double cdf_target) const {
        if (!isValid || cdf_target < 0.0 || cdf_target > 1.0) return -999.0;
        auto it = std::lower_bound(cdf_l1.begin(), cdf_l1.end(), cdf_target);
        if (it == cdf_l1.end()) return q_max;
        if (it == cdf_l1.begin()) return q_min;
        int idx = it - cdf_l1.begin();
        double den = cdf_l1[idx] - cdf_l1[idx-1];
        if (std::abs(den) < 1e-9) return q_values[idx-1];
        double frac = (cdf_target - cdf_l1[idx-1]) / den;
        return q_values[idx-1] + frac * (q_values[idx] - q_values[idx-1]);
    }
};

struct PdfData { 
    LookupTable lg_lookup, ege_lookup; 
    double ege_peak_l1=0, ege_peak_l2=0, ege_sigmaL_l1=0, ege_sigmaR_l1=0, ege_sigmaL_l2=0, ege_sigmaR_l2=0;
    double lg_mpv_l1=0, lg_mpv_l2=0, lg_width_l1=0, lg_width_l2=0, lg_sigma_l1=0, lg_sigma_l2=0;
    bool isValid = false;
};

struct AnalysisHistograms {
    std::vector<std::vector<TH2D*>> hists;
    AnalysisHistograms() : hists(NUM_ELEMENTS, std::vector<TH2D*>(NUM_DETECTORS, nullptr)) {}
};

// --- Main Analysis Function ---
void L2TempTuning_fast() {
    std::cout << "Creating tuned L2 template histograms from saveTree (Refactored)..." << std::endl;

    auto lookupFile = std::unique_ptr<TFile>(TFile::Open(lookupFileName.c_str()));
    auto inFile = std::unique_ptr<TFile>(TFile::Open(dataFileName.c_str()));
    
    if (!lookupFile || lookupFile->IsZombie() || !inFile || inFile->IsZombie()) {
        std::cerr << "Error opening input files!" << std::endl;
        return;
    }

    const std::vector<std::string> analysis_types = {"LG", "EGE", "LGLinear", "EGELinear"};
    std::map<std::string, std::unique_ptr<TFile>> out_files;
    std::map<std::string, std::array<AnalysisHistograms, 2>> all_hists; // Index 0 for normal L1, 1 for unbiased L1

    for (const auto& type : analysis_types) {
        std::string outFileName = output_base_path + "L2TempTuned_" + type + "_detfuncT.root";
        out_files[type] = std::make_unique<TFile>(outFileName.c_str(), "RECREATE");
        
        for (int htype = 0; htype < 2; ++htype) { // 0: biased, 1: unbiased
            std::string suffix = (htype == 1) ? "_unb" : "";
            for (int elem_idx = 0; elem_idx < NUM_ELEMENTS; ++elem_idx) {
                const auto& elem = elements[elem_idx];
                for (int idet = 0; idet < NUM_DETECTORS; ++idet) {
                    std::string det = detNames[idet];
                    TString hname = TString::Format("L2TempTuned_%s_%s_%s%s", type.c_str(), elem.name.c_str(), det.c_str(), suffix.c_str());
                    TString htitle = hname;
                    all_hists[type][htype].hists[elem_idx][idet] = new TH2D(hname, htitle, 400, elem.charge - 2.0, elem.charge + 2.0, Binning::NarrowBins.size() - 1, Binning::NarrowBins.data());
                }
            }
        }
    }

    std::map<std::string, PdfData> lookup_map;
    std::cout << "Loading lookup tables..." << std::endl;
    int startBin = findStartBin();
    std::cout << "Starting from bin " << startBin << " (Ek = " << Binning::NarrowBins[startBin] << " GeV)" << std::endl;

    for (int elem_idx = 0; elem_idx < NUM_ELEMENTS; ++elem_idx) {
        const auto& elem = elements[elem_idx];
        for (int idet = 0; idet < NUM_DETECTORS; ++idet) {
            std::string det = detNames[idet];
            for (size_t iy = startBin; iy < Binning::NarrowBins.size() - 1; ++iy) {
                double ekLow = Binning::NarrowBins[iy];
                if ((det=="TOF" && ekLow>1.55) || (det=="NaF" && (ekLow<0.86 || ekLow>4.91)) || (det=="AGL" && ekLow<2.88)) continue;

                std::string baseName = Form("%s_%s_bin%d", elem.name.c_str(), det.c_str(), iy);
                PdfData pdata;
                pdata.lg_lookup.loadFromFile(lookupFile.get(), baseName, "LG");
                pdata.ege_lookup.loadFromFile(lookupFile.get(), baseName, "EGE");

                auto get_params = [&](const std::string& type, const std::string& layer) {
                    return dynamic_cast<TVectorD*>(lookupFile->Get((baseName + "_" + type + "_params_" + layer).c_str()));
                };
                
                TVectorD* ege_params_l1 = get_params("EGE", "l1");
                TVectorD* ege_params_l2 = get_params("EGE", "l2");
                TVectorD* lg_params_l1 = get_params("LG", "l1");   
                TVectorD* lg_params_l2 = get_params("LG", "l2");  

                bool params_ok = true;
                if (ege_params_l1 && ege_params_l2) {
                    pdata.ege_peak_l1 = (*ege_params_l1)[0]; pdata.ege_sigmaL_l1 = (*ege_params_l1)[1]; pdata.ege_sigmaR_l1 = (*ege_params_l1)[3];
                    pdata.ege_peak_l2 = (*ege_params_l2)[0]; pdata.ege_sigmaL_l2 = (*ege_params_l2)[1]; pdata.ege_sigmaR_l2 = (*ege_params_l2)[3];
                } else { params_ok = false; }

                if (lg_params_l1 && lg_params_l2) {
                    pdata.lg_width_l1 = (*lg_params_l1)[0]; pdata.lg_mpv_l1 = (*lg_params_l1)[1]; pdata.lg_sigma_l1 = (*lg_params_l1)[2];
                    pdata.lg_width_l2 = (*lg_params_l2)[0]; pdata.lg_mpv_l2 = (*lg_params_l2)[1]; pdata.lg_sigma_l2 = (*lg_params_l2)[2];
                } else { params_ok = false; }
                
                pdata.isValid = pdata.lg_lookup.isValid && pdata.ege_lookup.isValid && params_ok;
                lookup_map[baseName] = std::move(pdata);
            }
        }
    }

    TTree* saveTree = dynamic_cast<TTree*>(inFile->Get("saveTree"));
    if (!saveTree) {
        std::cerr << "Error: Cannot find saveTree in the file" << std::endl;
        return;
    }

    EventProcessor ep;
    ep.setBranchAddresses(saveTree);

    std::clock_t start_clock = std::clock();
    time_t start_time = time(nullptr);
    Long64_t nEntries = saveTree->GetEntries();
    std::cout << "Processing " << nEntries << " entries..." << std::endl;
    long long tunedEventCount = 0;

    auto in_cut = [](double val, double z_val, double coe, double low, double high) {
        return val > z_val - coe*low && val < z_val + coe*high;
    };
    double coe = 0.8;

    for (Long64_t i = 0; i < nEntries; i++) {
        if (i > 0 && i % 1000000 == 0) {
            std::clock_t now_clock = std::clock();
            time_t now_time = time(nullptr);
            double elapsed_sec = double(now_clock - start_clock) / CLOCKS_PER_SEC;
            double wall_sec = difftime(now_time, start_time);
            double mem_mb = getCurrentRSS_MB();
            printf("Entry %lld / %lld | CPU(s): %.1f | Wall(s): %.1f | Mem: %.1f MB | Tuned: %lld\n",
                i, nEntries, elapsed_sec, wall_sec, mem_mb, tunedEventCount);
        }

        ep.getEntry(i);

        if (!ep.passBasicCut() && !ep.passUnbiasedCut()) continue;

        double beta_det[NUM_DETECTORS] = {ep.TOFBeta, ep.NaFBeta, ep.AGLBeta};
        double ek_det[NUM_DETECTORS]   = {ep.TOFEk, ep.NaFEk, ep.AGLEk};
        
        std::vector<double> Rbins_beta;
        Rbins_beta.reserve(Binning::NarrowBins.size());
        for (auto ek : Binning::NarrowBins) {
            Rbins_beta.push_back(kineticEnergyToBeta(ek));
        }

        for (int elem_idx = 0; elem_idx < NUM_ELEMENTS-1; elem_idx++) {
            const auto& elem = elements[elem_idx];
            int z = elem.charge;
            int mass = elem.mass;


            for (int idet = 0; idet < NUM_DETECTORS; ++idet) {
                if (!ep.passDetectorCut(idet)) continue;
                if (!ep.passBeyondCutoff(idet, z, mass, Rbins_beta)) continue;
                
                std::string det = detNames[idet];
                bool richQCut = (idet == 0) || ep.passRichQCut(z);

                bool shouldFill_L2Temp = ep.passBasicCut() && ep.getNormalL1XY() && ep.getL2XY() && ep.getQl2StatusCut() && ep.passQTrkInnerCut(2, z, coe) && ep.passTrkQL1Cut_L2Template(z, coe) && ep.passTofQUpCut(z, coe) && richQCut;
               
                bool shouldFill_L2Temp_unb = ep.passUnbiasedCut() && ep.getL2XY() && ep.getQl2StatusCut() && ep.passQTrkInnerCut(2, z, coe) && ep.passTrkQL1UnbiasCut_L2Template(z, coe) && ep.passTofQUpCut(z, coe) && richQCut;

                if (!shouldFill_L2Temp && !shouldFill_L2Temp_unb) continue;
                
                auto it_upper = std::lower_bound(Binning::NarrowBins.begin(), Binning::NarrowBins.end(), ek_det[idet]);
                int ekBin = std::distance(Binning::NarrowBins.begin(), it_upper) - 1;
                if (ekBin < 0) continue;

                std::string lookup_key = Form("%s_%s_bin%d", elem.name.c_str(), det.c_str(), ekBin);
                auto it = lookup_map.find(lookup_key);
                if (it == lookup_map.end() || !it->second.isValid) continue;
                
                const PdfData& pdata = it->second;
                
                // Calculate all tuned q values
                double q_tuned_LG = pdata.lg_lookup.getInvCDF_L1(pdata.lg_lookup.getCDF_L2(ep.trk_ql2));
                double q_tuned_EGE = pdata.ege_lookup.getInvCDF_L1(pdata.ege_lookup.getCDF_L2(ep.trk_ql2));
                double qLinear_lg = pdata.lg_mpv_l1 + (ep.trk_ql2 - pdata.lg_mpv_l2);
                double qLinear_ege = pdata.ege_peak_l1 + (ep.trk_ql2 - pdata.ege_peak_l2);
                
                std::map<std::string, double> tuned_qs = {
                    {"LG", q_tuned_LG}, {"EGE", q_tuned_EGE}, {"LGLinear", qLinear_lg}, {"EGELinear", qLinear_ege}
                };
                
                if (shouldFill_L2Temp) {
                    if (tunedEventCount < 10) { 
                        std::cout << elem.name << " normal " << det << " " << i << " ek:" << ek_det[idet]
                                << " ori q=" << ep.trk_ql2 << " lg q=" << q_tuned_LG << " ege q=" << q_tuned_EGE
                                << " lglinear q=" << qLinear_lg << " egelinear q=" << qLinear_ege<< std::endl;
                    }
                    for (const auto& [type, q_val] : tuned_qs) {
                        if (q_val > -900) all_hists[type][0].hists[elem_idx][idet]->Fill(q_val, ek_det[idet]);
                    }
                }
                if (shouldFill_L2Temp_unb) {
                    if (tunedEventCount < 10) { // Debug print for the first few events
                        std::cout << elem.name << " unbias " << det << " " << i << " ek:" << ek_det[idet]
                                << " ori q=" << ep.trk_ql2 << " lg q=" << q_tuned_LG << " ege q=" << q_tuned_EGE
                                << " lglinear q=" << qLinear_lg << " egelinear q=" << qLinear_ege<< std::endl;
                    }
                    for (const auto& [type, q_val] : tuned_qs) {
                        if (q_val > -900) all_hists[type][1].hists[elem_idx][idet]->Fill(q_val, ek_det[idet]);
                    }
                }

                if (shouldFill_L2Temp || shouldFill_L2Temp_unb) {
                    tunedEventCount++;
                }
            }
        }
    }
    
    std::clock_t end_clock = std::clock();
    time_t end_time = time(nullptr);
    printf("Total CPU time: %.1f s, Wall time: %.1f s, Max memory: %.1f MB\n",
        double(end_clock - start_clock) / CLOCKS_PER_SEC, difftime(end_time, start_time), getCurrentRSS_MB());
    std::cout << "Total tuned events: " << tunedEventCount << std::endl;

    std::cout << "Saving results..." << std::endl;
    for (const auto& type : analysis_types) {
        out_files[type]->cd();
        for (int htype = 0; htype < 2; ++htype) {
            for (int elem_idx = 0; elem_idx < NUM_ELEMENTS; ++elem_idx) {
                for (int idet = 0; idet < NUM_DETECTORS; ++idet) {
                    if(all_hists[type][htype].hists[elem_idx][idet]) {
                        all_hists[type][htype].hists[elem_idx][idet]->Write();
                    }
                }
            }
        }
        std::cout << "File " << out_files[type]->GetName() << " has been written." << std::endl;
        out_files[type]->Close();
    }
    std::cout << "L2 template tuning completed." << std::endl;
}