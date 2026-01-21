#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void drawacc() {
    // 打开两个根文件
    TFile *fileB10 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/BoronIsotopeResults/Acc/B10_acc.root");
    TFile *fileB11 = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/BoronIsotopeResults/Acc/B11_acc.root");

    // 检查文件是否成功打开
    if (!fileB10 || fileB10->IsZombie() || !fileB11 || fileB11->IsZombie()) {
        printf("Error opening one or both of the files!\n");
        return;
    }

    // 获取 B10 的三个直方图
    TH1F *hB10_0 = (TH1F*)fileB10->Get("h_Acc_hist_10_0");
    TH1F *hB10_1 = (TH1F*)fileB10->Get("h_Acc_hist_10_1");
    TH1F *hB10_2 = (TH1F*)fileB10->Get("h_Acc_hist_10_2");

    // 获取 B11 的三个直方图
    TH1F *hB11_0 = (TH1F*)fileB11->Get("h_Acc_hist_10_0");
    TH1F *hB11_1 = (TH1F*)fileB11->Get("h_Acc_hist_10_1");
    TH1F *hB11_2 = (TH1F*)fileB11->Get("h_Acc_hist_10_2");

    // 创建画布
    TCanvas *c1 = new TCanvas("c1", "Boron Isotope Acceptance", 800, 600);
    c1->SetLogx();

    int colori[] = {600, 880, 632, 797, 910, 415, 413, 1, 843, 860, 900, 800, 900};
    
    float max1 = TMath::Max(hB10_0->GetMaximum(), TMath::Max(hB10_1->GetMaximum(), hB10_2->GetMaximum()));
    float max2 = TMath::Max(hB11_0->GetMaximum(), TMath::Max(hB11_1->GetMaximum(), hB11_2->GetMaximum()));
    float maxVal = 1.0 * TMath::Max(max1, max2); // 给最大值留出5%的空间

    hB10_0->SetMaximum(0.07);
    hB10_0->SetMinimum(0.001); // 设置最小值为0

    hB10_0->SetTitle("");  // 去掉标题
    hB10_0->SetMarkerColor(colori[1]);
    hB10_0->SetMarkerStyle(20);  // marker 样式
    hB10_0->SetMarkerSize(1.0);  // marker 大小
    hB10_0->GetXaxis()->SetTitleOffset(1.15);
    hB10_0->GetYaxis()->SetTitleOffset(1.15);
    hB10_0->GetXaxis()->SetLabelSize(0.06);  // 设置 x 轴标签大小
    hB10_0->GetXaxis()->SetTitleSize(0.06);  // 设置 x 轴标题大小
    hB10_0->GetXaxis()->SetLabelSize(0.06);  // 设置 x 轴标签大小
    hB10_0->GetYaxis()->SetTitleSize(0.06);  // 设置 y 轴标题大小
    hB10_0->GetYaxis()->SetLabelSize(0.06);  // 设置 y 轴标签大小
    hB10_0->Draw("P");  // "P" 表示绘制点

    hB10_1->SetMarkerColor(colori[0]);
    hB10_1->SetMarkerStyle(20);
    hB10_1->SetMarkerSize(1.0);
    hB10_1->Draw("SAME P");  // "SAME" 表示在同一张图上绘制

    hB10_2->SetMarkerColor(colori[6]);
    hB10_2->SetMarkerStyle(20);
    hB10_2->SetMarkerSize(1.0);
    hB10_2->Draw("SAME P");

    hB11_0->SetMarkerColor(colori[3]);
    hB11_0->SetMarkerStyle(20);
    hB11_0->SetMarkerSize(1.0);
    hB11_0->Draw("SAME P");

    hB11_1->SetMarkerColor(colori[4]);
    hB11_1->SetMarkerStyle(20);
    hB11_1->SetMarkerSize(1.0);
    hB11_1->Draw("SAME P");

    hB11_2->SetMarkerColor(colori[5]);
    hB11_2->SetMarkerStyle(20);
    hB11_2->SetMarkerSize(1.0);
    hB11_2->Draw("SAME P");

    // 保存为 PNG 文件
    c1->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/BoronIsotopeResults/Acc/acc_iso_cp.png");

    // 关闭 ROOT 文件
    fileB10->Close();
    fileB11->Close();

    // 新增功能：处理 BorMix32_68_acc.root 文件中的直方图
    TFile *fileBorMix = TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/BoronIsotopeResults/Acc/BorMix32_68_acc.root");

    if (!fileBorMix || fileBorMix->IsZombie()) {
        printf("Error opening BorMix32_68_acc.root file!\n");
        return;
    }

    // 获取 BorMix32_68_acc.root 文件中的四个直方图
    TH1F *hBorMix_0 = (TH1F*)fileBorMix->Get("h_Acc_hist_10_0");
    TH1F *hBorMix_1 = (TH1F*)fileBorMix->Get("h_Acc_hist_11_0");
    TH1F *hBorMix_2 = (TH1F*)fileBorMix->Get("h_Acc_hist_13_0");
    TH1F *hBorMix_3 = (TH1F*)fileBorMix->Get("h_Acc_hist_14_0");

    // 创建第二张画布
    TCanvas *c2 = new TCanvas("c2", "Boron Isotope TOF Acceptance", 800, 600);
    c2->SetLogx();  // 设置 x 轴为对数坐标

    // 设置颜色：蓝、红、kGreen+2、Magnet
    hBorMix_0->SetMarkerColor(kBlue);
    hBorMix_1->SetMarkerColor(kRed);
    hBorMix_2->SetMarkerColor(kGreen+2);
    hBorMix_3->SetMarkerColor(kMagenta);

    hBorMix_0->SetTitle("");  // 去掉标题
    hBorMix_0->GetXaxis()->SetTitleOffset(1.15);
    hBorMix_0->GetYaxis()->SetTitleOffset(1.15);
    hBorMix_0->GetXaxis()->SetLabelSize(0.06);  // 设置 x 轴标签大小
    hBorMix_0->GetXaxis()->SetTitleSize(0.06);  // 设置 x 轴标题大小
    hBorMix_0->GetXaxis()->SetLabelSize(0.06);  // 设置 x 轴标签大小
    hBorMix_0->GetYaxis()->SetTitleSize(0.06);  // 设置 y 轴标题大小
    hBorMix_0->GetYaxis()->SetLabelSize(0.06);  // 设置 y 轴标签大小
    // 设置 marker 样式和大小
    hBorMix_0->SetMarkerStyle(20);
    hBorMix_0->SetMarkerSize(1.0);
    hBorMix_1->SetMarkerStyle(20);
    hBorMix_1->SetMarkerSize(1.0);
    hBorMix_2->SetMarkerStyle(20);
    hBorMix_2->SetMarkerSize(1.0);
    hBorMix_3->SetMarkerStyle(20);
    hBorMix_3->SetMarkerSize(1.0);

    // 设置 y 轴范围
    hBorMix_0->SetMaximum(0.07);
    hBorMix_0->SetMinimum(0.001);  // 设置最小值

    // 绘制直方图并叠加
    hBorMix_0->Draw("P");
    hBorMix_1->Draw("SAME P");
    hBorMix_2->Draw("SAME P");
    hBorMix_3->Draw("SAME P");

    // 保存第二张图片为 PNG 文件
    c2->SaveAs("/eos/ams/user/z/zuhao/yanzx/Isotope/BoronIsotopeResults/Acc/acc_tof_cp.png");

    // 关闭 ROOT 文件
    fileBorMix->Close();
}