#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLegend.h>
#include <iostream>

using std::cout;
using std::endl;

static TH1F* GetH1(TFile* f, const char* name) {
    TH1F* h = f ? dynamic_cast<TH1F*>(f->Get(name)) : nullptr;
    if (!h) {
        std::cerr << "[ERROR] Histogram not found: " << (name ? name : "(null)") 
                  << " in file " << (f ? f->GetName() : "(null)") << std::endl;
    }
    return h;
}

// 做一次简单的拷贝，避免修改原对象
static TH1F* CloneAs(const TH1F* src, const char* newname) {
    if (!src) return nullptr;
    TH1F* h = dynamic_cast<TH1F*>(src->Clone(newname));
    if (h) h->SetDirectory(nullptr);
    return h;
}

void cp_plot() {
    gStyle->SetOptStat(0);

    // ========== 第一组对比 ==========
    // 对象1：/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/BeTotal_analysis.root 里的 h_acc_orig_Complete_Bkg_TOF_Be10
    // 对象2：/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/Beryllium10_bkg_analysis_mtrpar1.root 里的 h_acc_orig_Complete_Bkg_TOF_Be10
    // 对象3：/eos/user/z/zixuan/Isotope/Add/Be10_all_w1_frag4.root 中的
    //        UnbiasedL1Inner_MC_BKG_H3a_TOF_Z4_Mass10 divide MC_FLUX_H3 再 scale 3.9×3.9×Pi

    const char* f1_path = "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/BeTotal_analysis.root";
    const char* f2_path = "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/Beryllium10_bkg_analysis_mtrpar1.root";
    const char* f3_path = "/eos/user/z/zixuan/Isotope/Add/Be10_all_w1_frag4.root";

    TFile* f1 = TFile::Open(f1_path, "READ");
    TFile* f2 = TFile::Open(f2_path, "READ");
    TFile* f3 = TFile::Open(f3_path, "READ");

    TH1F* h1_a = GetH1(f1, "h_acc_orig_Complete_Bkg_TOF_Be10"); // 第一组 对象1
    TH1F* h1_b = GetH1(f2, "h_acc_orig_Complete_Bkg_TOF_Be10"); // 第一组 对象2

    // 对象3的构造：h = H3a / MC_FLUX_H3 * (3.9*3.9*pi)
    TH1F* h3_h3a = GetH1(f3, "UnbiasedL1Inner_MC_BKG_H3a_TOF_Z4_Mass10");
    TH1F* h3_flux = GetH1(f3, "MC_FLUX_H3");

    TH1F* h1_c = nullptr;
    if (h3_h3a && h3_flux) {
        h1_c = CloneAs(h3_h3a, "h_first_group_obj3");
        // 不画误差，简单按内容算
        h1_c->Divide(h3_flux);
        h1_c->Scale(3.9 * 3.9 * TMath::Pi());
    }

    // 作图设置与绘制（红、黑、蓝；marker 20；无误差，直线/点均可，这里统一用 P SAME）
    TCanvas* c1 = new TCanvas("c_group1", "Group 1 Comparison", 1000, 700);
    c1->SetGrid();
    c1->SetLogx();

    // 为了合理的y轴范围，拿一个框架：复制对象1（若为空则找下一个）
    TH1F* frame_src = h1_a ? h1_a : (h1_b ? h1_b : h1_c);
    if (!frame_src) {
        std::cerr << "[FATAL] Group 1: no histogram available to draw." << std::endl;
    } else {
        TH1F* frame = CloneAs(frame_src, "frame_group1");
        frame->Reset("ICESM");
        frame->SetTitle("Group 1: Be10 - TOF (acc/bkg/H3a-based);E_{k}/n;value");
        frame->Draw("AXIS"); // 用于坐标轴，随后用 hist 同步上去
        // 先把第一条真正画上来作为范围参考
        if (h1_a) {
            TH1F* h = CloneAs(h1_a, "h1_a_draw");
            h->SetLineColor(kRed);
            h->SetMarkerColor(kRed);
            h->SetMarkerStyle(20);
            h->SetLineWidth(2);
            h->Draw("HIST");
        }
        if (h1_b) {
            TH1F* h = CloneAs(h1_b, "h1_b_draw");
            h->SetLineColor(kBlack);
            h->SetMarkerColor(kBlack);
            h->SetMarkerStyle(20);
            h->SetLineWidth(2);
            h->Draw(h1_a ? "HIST SAME" : "HIST");
        }
        if (h1_c) {
            TH1F* h = CloneAs(h1_c, "h1_c_draw");
            h->SetLineColor(kBlue);
            h->SetMarkerColor(kBlue);
            h->SetMarkerStyle(20);
            h->SetLineWidth(2);
            h->Draw((h1_a || h1_b) ? "HIST SAME" : "HIST");
        }
        c1->Update();
        c1->SaveAs("/eos/user/z/zixuan/Isotope/BkgEst/group1.png");
    }

    // ========== 第二组对比 ==========
    // 对象1：BeTotal_analysis.root 里的 h_bkg_orig_Complete_Bkg_TOF_Be10
    // 对象2：Beryllium10_bkg_analysis_mtrpar1.root 里的 h_bkg_orig_Complete_Bkg_TOF_Be10
    // 对象3：同第一组中的 UnbiasedL1Inner_MC_BKG_H3a_TOF_Z4_Mass10 处理（保持一致对比该“bkg”）
    // 注：你要求“hist name 里的 acc 换为 bkg 其他不变，以及 UnbiasedL1Inner_MC_BKG_H3a_TOF_Z4_Mass10”

    TH1F* h2_a = GetH1(f1, "h_bkg_orig_Complete_Bkg_TOF_Be10"); // 第二组 对象1
    TH1F* h2_b = GetH1(f2, "h_bkg_orig_Complete_Bkg_TOF_Be10"); // 第二组 对象2

    // 对象3沿用第一组的 h1_c（基于 H3a 的构造）。为避免共享同一指针，克隆一份。
    TH1F* h2_c = h1_c ? CloneAs(h3_h3a, "h_second_group_obj3") : nullptr;

    TCanvas* c2 = new TCanvas("c_group2", "Group 2 Comparison", 1000, 700);
    c2->SetGrid();
    c2->SetLogx();

    TH1F* frame_src2 = h2_a ? h2_a : (h2_b ? h2_b : h2_c);
    if (!frame_src2) {
        std::cerr << "[FATAL] Group 2: no histogram available to draw." << std::endl;
    } else {
        TH1F* frame2 = CloneAs(frame_src2, "frame_group2");
        frame2->Reset("ICESM");
        frame2->SetTitle("Group 2: Be10 - TOF (bkg vs H3a-based);E_{k}/n;value");
        frame2->Draw("AXIS");

        if (h2_a) {
            TH1F* h = CloneAs(h2_a, "h2_a_draw");
            h->SetLineColor(kRed);
            h->SetMarkerColor(kRed);
            h->SetMarkerStyle(20);
            h->SetLineWidth(2);
            h->Draw("HIST");
        }
        if (h2_b) {
            TH1F* h = CloneAs(h2_b, "h2_b_draw");
            h->SetLineColor(kBlack);
            h->SetMarkerColor(kBlack);
            h->SetMarkerStyle(20);
            h->SetLineWidth(2);
            h->Draw(h2_a ? "HIST SAME" : "HIST");
        }
        if (h2_c) {
            TH1F* h = CloneAs(h2_c, "h2_c_draw");
            h->SetLineColor(kBlue);
            h->SetMarkerColor(kBlue);
            h->SetMarkerStyle(20);
            h->SetLineWidth(2);
            h->Draw((h2_a || h2_b) ? "HIST SAME" : "HIST");
        }
        c2->Update();
        c2->SaveAs("/eos/user/z/zixuan/Isotope/BkgEst/group2.png");
    }

    // 清理与关闭
    if (f1) f1->Close();
    if (f2) f2->Close();
    if (f3) f3->Close();

    cout << "Done. Saved plots:\n"
         << "  /eos/user/z/zixuan/Isotope/BkgEst/group1.png\n"
         << "  /eos/user/z/zixuan/Isotope/BkgEst/group2.png\n";
}