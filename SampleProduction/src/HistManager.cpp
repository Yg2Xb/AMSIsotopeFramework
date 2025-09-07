#include "HistManager.h"
#include "BinningManager.h"
#include "basic_var.h"
#include <iostream>
#include <stdexcept>

namespace AMS_Iso {

static inline std::vector<double> safeBins(const std::string& key) {
    try {
        const auto& bins = BinningManager::GetInstance().Get(key);
        if (!bins.empty()) return bins;
    } catch (...) {}
    return std::vector<double>{0,1};
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
    const std::vector<std::string> detectors   = {"TOF","NaF","AGL"};
    const std::vector<std::string> cut_groups  = {"CutGroup1","CutGroup2","CutGroup3"};
    const std::vector<std::string> num_den     = {"Num","Den"};
    const std::vector<std::string> charge_types= {"L1QSignal","L1QTemplate","L2QTemplate"};
    const std::vector<std::string> sources     = {"Beryllium","Boron","Carbon","Nitrogen","Oxygen"};

    int Nchain = chains.size();
    int Ndet   = detectors.size();
    int Niso   = isISS ? iso->getIsotopeCount() : 1;
    int Nsrc   = sources.size();
    int Nct    = charge_types.size();
    int C      = cut_groups.size();
    int Nd     = num_den.size();

    // ----------------  ID 区域 ----------------
    if (isISS) {
        ISS_IDH1.resize(Nchain,std::vector<std::vector<TH1F*>>(Ndet,std::vector<TH1F*>(Niso,nullptr)));
    } else {
        MC_IDH1.resize(Nchain,std::vector<TH2F*>(Ndet,nullptr));
    }
    IDH2.resize(Nchain,std::vector<std::vector<TH2F*>>(Ndet,std::vector<TH2F*>(Niso,nullptr)));
    IDH3.resize(Nchain,std::vector<TH2F*>(Ndet,nullptr));
    IDH4a.resize(Nchain,nullptr);
    IDH4b.resize(Nchain,nullptr);
    IDH5a.resize(Nchain,nullptr);
    IDH5b.resize(Nchain,nullptr);
    IDH6a.resize(Nchain,nullptr);
    IDH6b.resize(Nchain,nullptr);
    IDH7a.resize(Nchain,nullptr);
    IDH7b.resize(Nchain,nullptr);

    for (int c=0;c<Nchain;++c){
        for (int d=0;d<Ndet;++d){
            for (int i=0;i<Niso;++i){
                int mass = isISS ? iso->getMass(i) : UseMass;
                auto ekBins  = binMgr.GetEkPerNucleonBins(charge,mass);
                auto invMBins= safeBins("InverseMass");

                if (isISS) {
                    TString hname=Form("%s_IDH1_%s_%d",chains[c].c_str(),detectors[d].c_str(),mass);
                    ISS_IDH1[c][d][i]= new TH1F(hname,hname,ekBins.size()-1,ekBins.data());
                } else {
                    TString hname=Form("%s_MC_IDH1_%s",chains[c].c_str(),detectors[d].c_str());
                    MC_IDH1[c][d]= new TH2F(hname,hname,ekBins.size()-1,ekBins.data(),
                                             invMBins.size()-1,invMBins.data());
                }
                TString hname2=Form("%s_IDH2_%s_%d",chains[c].c_str(),detectors[d].c_str(),mass);
                IDH2[c][d][i]= new TH2F(hname2,hname2,ekBins.size()-1,ekBins.data(),
                                        invMBins.size()-1,invMBins.data());
            }
            auto ekBins  = safeBins("EkPerNucleon");
            auto invMBins= safeBins("InverseMass");
            TString h3=Form("%s_IDH3_%s",chains[c].c_str(),detectors[d].c_str());
            IDH3[c][d]= new TH2F(h3,h3,ekBins.size()-1,ekBins.data(),
                                 invMBins.size()-1,invMBins.data());
        }
        auto betaBins = safeBins("Beta");
        auto dBetaBins= safeBins("DeltaBeta");
        auto rigBins  = safeBins("Rigidity");
        auto bgBins   = safeBins("BetaGamma");
        auto ekBins   = safeBins("EkPerNucleon");
        IDH4a[c]=new TH1F(Form("%s_IDH4a",chains[c].c_str()),"",betaBins.size()-1,betaBins.data());
        IDH4b[c]=new TH1F(Form("%s_IDH4b",chains[c].c_str()),"",betaBins.size()-1,betaBins.data());
        IDH5a[c]=new TH2F(Form("%s_IDH5a",chains[c].c_str()),"",rigBins.size()-1,rigBins.data(),dBetaBins.size()-1,dBetaBins.data());
        IDH5b[c]=new TH2F(Form("%s_IDH5b",chains[c].c_str()),"",rigBins.size()-1,rigBins.data(),dBetaBins.size()-1,dBetaBins.data());
        IDH6a[c]=new TH2F(Form("%s_IDH6a",chains[c].c_str()),"",ekBins.size()-1,ekBins.data(),dBetaBins.size()-1,dBetaBins.data());
        IDH6b[c]=new TH2F(Form("%s_IDH6b",chains[c].c_str()),"",ekBins.size()-1,ekBins.data(),dBetaBins.size()-1,dBetaBins.data());
        IDH7a[c]=new TH2F(Form("%s_IDH7a",chains[c].c_str()),"",bgBins.size()-1,bgBins.data(),dBetaBins.size()-1,dBetaBins.data());
        IDH7b[c]=new TH2F(Form("%s_IDH7b",chains[c].c_str()),"",bgBins.size()-1,bgBins.data(),dBetaBins.size()-1,dBetaBins.data());
    }

    // ----------------  BKG 区域 ----------------
    if (isISS){
        ISS_BKGH1.resize(Nchain,std::vector<std::vector<TH1F*>>(Nsrc,std::vector<TH1F*>(Ndet,nullptr)));
        ISS_BKGH2.resize(Nchain,std::vector<std::vector<std::vector<TH2F*>>>(Nsrc,std::vector<std::vector<TH2F*>>(Nct,std::vector<TH2F*>(Ndet,nullptr))));
        ISS_BKGH3.resize(Nchain,std::vector<std::vector<TH1F*>>(Nsrc,std::vector<TH1F*>(Ndet,nullptr)));
        ISS_BKGH4.resize(Nchain,std::vector<std::vector<TH2F*>>(Nsrc,std::vector<TH2F*>(Ndet,nullptr)));
        for(int c=0;c<Nchain;++c){
            for(int s=0;s<Nsrc;++s){
                for(int d=0;d<Ndet;++d){
                    auto ekBins= safeBins("EkPerNucleon");
                    auto invMBins= safeBins("InverseMass");
                    ISS_BKGH1[c][s][d]= new TH1F(Form("%s_BKGH1_%s_%s",chains[c].c_str(),sources[s].c_str(),detectors[d].c_str()),"",ekBins.size()-1,ekBins.data());
                    for(int t=0;t<Nct;++t){
                        ISS_BKGH2[c][s][t][d]= new TH2F(Form("%s_BKGH2_%s_%s_%s",chains[c].c_str(),sources[s].c_str(),charge_types[t].c_str(),detectors[d].c_str()),"",ekBins.size()-1,ekBins.data(),10,0,10);
                    }
                    ISS_BKGH3[c][s][d]= new TH1F(Form("%s_BKGH3_%s_%s",chains[c].c_str(),sources[s].c_str(),detectors[d].c_str()),"",ekBins.size()-1,ekBins.data());
                    ISS_BKGH4[c][s][d]= new TH2F(Form("%s_BKGH4_%s_%s",chains[c].c_str(),sources[s].c_str(),detectors[d].c_str()),"",ekBins.size()-1,ekBins.data(),invMBins.size()-1,invMBins.data());
                }
            }
        }
    } else {
        MC_BKGH1.resize(Nchain,std::vector<TH1F*>(Ndet,nullptr));
        MC_BKGH2.resize(Nchain,std::vector<std::vector<TH1F*>>(Ndet,std::vector<TH1F*>(Niso,nullptr)));
        MC_BKGH3a=MC_BKGH2; MC_BKGH3b=MC_BKGH2; // same dimensions
        for(int c=0;c<Nchain;++c){
            for(int d=0;d<Ndet;++d){
                auto ekBins= safeBins("EkPerNucleon");
                MC_BKGH1[c][d]=new TH1F(Form("%s_MC_BKGH1_%s",chains[c].c_str(),detectors[d].c_str()),"",ekBins.size()-1,ekBins.data());
                for(int i=0;i<Niso;++i){
                    int mass=UseMass;
                    auto ekBinsIso=binMgr.GetEkPerNucleonBins(charge,mass);
                    MC_BKGH2[c][d][i]=new TH1F(Form("%s_MC_BKGH2_%s_%d",chains[c].c_str(),detectors[d].c_str(),mass),"",ekBinsIso.size()-1,ekBinsIso.data());
                    MC_BKGH3a[c][d][i]=new TH1F(Form("%s_MC_BKGH3a_%s_%d",chains[c].c_str(),detectors[d].c_str(),mass),"",ekBinsIso.size()-1,ekBinsIso.data());
                    MC_BKGH3b[c][d][i]=new TH1F(Form("%s_MC_BKGH3b_%s_%d",chains[c].c_str(),detectors[d].c_str(),mass),"",ekBinsIso.size()-1,ekBinsIso.data());
                }
            }
        }
    }

    // ---------------- FLUX 区域 ----------------
    FLUXH1.resize(Nchain,std::vector<std::vector<std::vector<std::vector<TH1F*>>>>(
                  C,std::vector<std::vector<std::vector<TH1F*>>>(
                  Nd,std::vector<std::vector<TH1F*>>(Ndet,std::vector<TH1F*>(Niso,nullptr)))));
    FLUXH2.resize(Nchain,nullptr);
    FLUXH3.resize(Nchain,std::vector<std::vector<TH1F*>>(Ndet,std::vector<TH1F*>(Niso,nullptr)));

    for(int c=0;c<Nchain;++c){
        for(int cg=0;cg<C;++cg)for(int nd=0;nd<Nd;++nd)for(int d=0;d<Ndet;++d)for(int i=0;i<Niso;++i){
            int mass= isISS ? iso->getMass(i):UseMass;
            auto ekBins=binMgr.GetEkPerNucleonBins(charge,mass);
            FLUXH1[c][cg][nd][d][i]=new TH1F(Form("%s_FLUXH1_%s_%s_%s_%d",chains[c].c_str(),cut_groups[cg].c_str(),num_den[nd].c_str(),detectors[d].c_str(),mass),"",ekBins.size()-1,ekBins.data());
        }
        auto rigBins=safeBins("Rigidity");
        FLUXH2[c]=new TH1F(Form("%s_FLUXH2",chains[c].c_str()),"",rigBins.size()-1,rigBins.data());
        for(int d=0;d<Ndet;++d)for(int i=0;i<Niso;++i){
            int mass=isISS?iso->getMass(i):UseMass;
            auto ekBins=binMgr.GetEkPerNucleonBins(charge,mass);
            FLUXH3[c][d][i]=new TH1F(Form("%s_FLUXH3_%s_%d",chains[c].c_str(),detectors[d].c_str(),mass),"",ekBins.size()-1,ekBins.data());
        }
    }
    MC_FLUXH2.resize(Nchain,nullptr);
    MC_FLUXH3.resize(Nchain,nullptr);
    if(!isISS){
        auto ekBins=safeBins("EkGen");
        for(int c=0;c<Nchain;++c){
            MC_FLUXH2[c]=new TH1F(Form("%s_MC_FLUXH2",chains[c].c_str()),"",ekBins.size()-1,ekBins.data());
            MC_FLUXH3[c]=new TH1F(Form("%s_MC_FLUXH3",chains[c].c_str()),"",10,0,10);
        }
    }

    std::cout<<"HistManager initialized all hist arrays"<<std::endl;
}

HistManager::~HistManager(){
    if (m_outputFile && m_outputFile->IsOpen()){
        m_outputFile->Write();
        m_outputFile->Close();
    }
}

void HistManager::Save(){
    if (m_outputFile && m_outputFile->IsOpen()){
        m_outputFile->cd();
        for(auto& c:ISS_IDH1) for(auto& d:c) for(auto* h:d) if(h)h->Write();
        for(auto& c:MC_IDH1) for(auto* h:c) if(h)h->Write();
        for(auto& c:IDH2) for(auto& d:c) for(auto* h:d) if(h)h->Write();
        for(auto& c:IDH3) for(auto* h:c) if(h)h->Write();
        for(auto* h:IDH4a) if(h)h->Write();
        for(auto* h:IDH4b) if(h)h->Write();
        for(auto* h:IDH5a) if(h)h->Write();
        for(auto* h:IDH5b) if(h)h->Write();
        for(auto* h:IDH6a) if(h)h->Write();
        for(auto* h:IDH6b) if(h)h->Write();
        for(auto* h:IDH7a) if(h)h->Write();
        for(auto* h:IDH7b) if(h)h->Write();
        std::cout<<"HistManager saved histograms"<<std::endl;
    }
}

} // namespace AMS_Iso