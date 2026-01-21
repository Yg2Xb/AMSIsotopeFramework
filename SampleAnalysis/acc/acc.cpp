#include <TFile.h>
#include <TH1.h>
#include <iostream>

void acc() {
    // 常量
    const double scale = 3.9 * 3.9 * 3.141592653589793;
    
    // 创建输出文件
    TFile* fout = new TFile("/eos/user/z/zixuan/Isotope/Acc/Acceptance.root", "RECREATE");
    
    // 粒子、链、探测器列表
    const char* particles[] = {"Be7", "Be9", "Be10"};
    const char* chains[] = {"UnbiasedL1Inner", "L1Inner"};
    const char* detectors[] = {"TOF", "NaF", "AGL"};
    int masses[] = {7, 9, 10};
    
    for (int i = 0; i < 3; i++) {
        // 打开输入文件
        TString fin_name = TString::Format("/eos/user/z/zixuan/Isotope/Add/%s_rew_frag4_withBkg_full.root", particles[i]);
        TFile* fin = TFile::Open(fin_name);
        if (!fin) continue;
        
        // 获取MC_FLUX_H3
        TH1* h_mc = (TH1*)fin->Get("MC_FLUX_H3");
        
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 3; k++) {
                // 构建直方图名
                TString h3a_name = TString::Format("%s_BKG_H3a_%s_Z4_Mass%d", chains[j], detectors[k], masses[i]);
                TH1* h_h3a = (TH1*)fin->Get(h3a_name);
                if (!h_h3a || !h_mc) continue;
                
                // 计算Acceptance
                TH1* h_acc = (TH1*)h_h3a->Clone();
                h_acc->Divide(h_mc);
                h_acc->Scale(scale);
                
                // 重命名保存
                TString acc_name = TString::Format("Acceptance_%s_%s_%s", chains[j], detectors[k], particles[i]);
                h_acc->SetName(acc_name);
                h_acc->SetTitle(TString::Format("Acceptance %s %s %s", chains[j], detectors[k], particles[i]));
                
                // 保存
                fout->cd();
                h_acc->Write();
                
                delete h_acc;
            }
        }
        fin->Close();
        delete fin;
    }
    
    fout->Close();
    delete fout;
    
    std::cout << "Done! Acceptance saved to /eos/user/z/zixuan/Isotope/Acc/Acceptance.root" << std::endl;
}