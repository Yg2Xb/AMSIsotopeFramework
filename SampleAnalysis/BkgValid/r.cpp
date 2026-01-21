#include <iostream>
#include <vector>
#include <cmath>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TAxis.h"

void r() {
    // --- 准备数据 ---
    // 能量分箱 E_k/n [GeV] (对齐两图的近似分箱中心)
    std::vector<double> energy = {0.5, 1.0, 1.5, 2.2, 2.8, 3.5, 4.2, 5.2, 6.2, 7.5, 9.0, 10.6, 12.5, 14.7, 17.2, 20.1};
    
    // 图1 (Beryllium) 的 ISS/MC 及其误差
    std::vector<double> y1 = {0.82, 0.77, 0.62, 0.69, 0.56, 0.62, 0.62, 0.49, 0.48, 0.76, 0.55, 0.80, 0.70, 0.53, 0.85, 0.74};
    std::vector<double> e1 = {0.03, 0.04, 0.08, 0.10, 0.11, 0.05, 0.05, 0.06, 0.08, 0.12, 0.15, 0.22, 0.28, 0.45, 0.50, 0.60};

    // 图2 (Boron -> 10Be) 的 ISS/MC 及其误差
    std::vector<double> y2 = {0.85, 0.72, 0.64, 0.78, 0.58, 0.61, 0.54, 0.38, 0.38, 0.74, 0.61, 0.65, 0.71, 0.38, 1.15, 1.88};
    std::vector<double> e2 = {0.04, 0.10, 0.12, 0.12, 0.15, 0.08, 0.08, 0.08, 0.10, 0.12, 0.22, 0.30, 0.40, 1.00, 1.20, 1.50};

    int n = energy.size();
    std::vector<double> ratio(n), ratio_err(n), zero_err(n, 0.0);

    // --- 计算比例 (Fig2 / Fig1) ---
    for(int i=0; i<n; ++i) {
        ratio[i] = y2[i] / y1[i];
        // 误差传递公式: R = y2/y1 -> sigma_R = R * sqrt((e1/y1)^2 + (e2/y2)^2)
        ratio_err[i] = ratio[i] * std::sqrt(std::pow(e1[i]/y1[i], 2) + std::pow(e2[i]/y2[i], 2));
    }

    // --- 绘图 ---
    TCanvas *c1 = new TCanvas("c1", "Ratio Plot", 800, 600);
    c1->SetGrid();

    TGraphErrors *gr = new TGraphErrors(n, &energy[0], &ratio[0], &zero_err[0], &ratio_err[0]);
    gr->SetTitle("Ratio of ISS/MC (Fig 2 / Fig 1) vs Energy;E_{k}/n [GeV];Ratio (10Be / Be)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    
    gr->Draw("AP");

    // 设置坐标轴范围
    gr->GetYaxis()->SetRangeUser(0.0, 4.0);
    gr->GetXaxis()->SetLimits(0.0, 22.0);

    // 添加基准线 (ratio = 1)
    TLine *line = new TLine(0, 1, 22, 1);
    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw();

    c1->SaveAs("ISS_MC_Ratio.pdf");
    c1->SaveAs("ISS_MC_Ratio.png");
}