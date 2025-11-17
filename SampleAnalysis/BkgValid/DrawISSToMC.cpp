#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <map>
#include <utility> // for std::pair

// --- 全局样式和常量 ---
const TString path_noBkgCut = "/eos/user/z/zixuan/Isotope/BkgValid/";
const TString path_withBkgCut = "/eos/user/z/zixuan/Isotope/BkgValid_B/";
const TString plot_output_path = "/eos/user/z/zixuan/Isotope/BkgValid/Plots/ratio/";

// 颜色顺序 (使用您更新的)
const std::vector<int> colors = {kGreen + 2, kRed, kBlue, kOrange + 1, kMagenta, kCyan + 2, kGray + 2, kBlack};
// Marker 样式
const int marker_style = 20; // 实心圆 (kFullCircle)
// 轴范围 (使用您更新的)
const double x_min = 0.1;
const double x_max = 21.5;
const double y_ranges[2][2] = {{0.2, 2}, {0., 4.5}};

/**
 * @brief 辅助函数：从文件中获取 ISS/MC Ratio 直方图
 * (此函数保持不变)
 */
TH1D* getRatio(TString filePath, TString histName_iss, TString histName_mc) {
    TFile* f = TFile::Open(filePath);
    if (!f || f->IsZombie()) {
        return nullptr;
    }
    TH1* h_iss = (TH1*)f->Get(histName_iss);
    TH1* h_mc = (TH1*)f->Get(histName_mc);
    if (!h_iss) {
        f->Close(); delete f;
        return nullptr;
    }
    if (!h_mc) {
        f->Close(); delete f;
        return nullptr;
    }
    TH1D* h_ratio = (TH1D*)h_iss->Clone(TString(histName_iss) + "_ratio");
    h_ratio->SetDirectory(0); 
    h_ratio->Sumw2(kTRUE);
    h_mc->Sumw2(kTRUE);
    h_ratio->Divide(h_mc);
    for (int i = 0; i <= h_ratio->GetNbinsX() + 1; ++i) {
        h_ratio->SetBinError(i, 0);
    }
    f->Close();
    delete f;
    return h_ratio;
}

/**
 * @brief 辅助函数：设置直方图样式
 * (此函数保持不变)
 */
void styleHist(TH1* hist, int color, int markerStyle) {
    if (!hist) return;
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(1.2);
    hist->GetXaxis()->SetRangeUser(x_min, x_max);
    hist->GetXaxis()->SetTitle("Measured Ek/n[GeV]"); 
    hist->GetYaxis()->SetTitle("ISS / MC Ratio");
    hist->SetStats(0); // 关闭统计框
}

/**
 * @brief 辅助函数：创建并设置图例
 * (此函数保持不变)
 */
TLegend* createLegend() {
    TLegend* leg = new TLegend(0.7, 0.72, 0.88, 0.87); // 调整了位置，避免遮挡
    leg->SetFillStyle(0);   // 透明背景
    leg->SetBorderSize(1);  // 有边框 (根据您的代码)
    leg->SetTextSize(0.03);
    return leg;
}

/**
 * @brief 新增辅助函数：查找直方图在 X 轴范围内的 Y 轴极值
 */
std::pair<double, double> findMinMaxInRange(TH1* hist, double xmin, double xmax) {
    if (!hist) return {0, 0};

    int first_bin = hist->GetXaxis()->FindBin(xmin);
    int last_bin = hist->GetXaxis()->FindBin(xmax);

    double local_min = 1e30;
    double local_max = -1e30;
    bool data_found = false;

    for (int i = first_bin; i <= last_bin; ++i) {
        double val = hist->GetBinContent(i);
        // 仅在 bin 有内容时才更新极值
        if (val != 0 || hist->GetBinError(i) != 0) { 
            if (val < local_min) local_min = val;
            if (val > local_max) local_max = val;
            data_found = true;
        }
    }

    if (!data_found) return {0, 0}; // 如果范围内没有数据
    return {0.8*local_min, 1.2*local_max};
}


// --- 主函数 (函数名保持不变) ---
void DrawISSToMC() {
    gROOT->SetBatch(kTRUE); // 后台运行，不显示画布
    gStyle->SetOptStat(0);   // 全局关闭统计框

    gSystem->Exec("mkdir -p " + plot_output_path);

    TCanvas* c = new TCanvas("c", "Bkg Ratios", 800, 600);
    c->SetGrid(); // 设置网格 (保持)
    TLegend* leg;

    // --- 组 1: Boron -> Beryllium (isotopes) ---
    {
        std::vector<TString> isotopes_be = {"Be7", "Be9", "Be10", "total"};
        std::vector<TString> labels_be = {"^{7}Be", "^{9}Be", "^{10}Be", "Total (B#rightarrowBe)"};
        TString filename = "Boron_to_Beryllium_L1Inner_Validation.root";

        for (int i_range = 0; i_range < 2; ++i_range) {
            
            // --- 1.1: No BkgCut (重构) ---
            {
                c->Clear();
                c->SetLogx(0); c->SetLogy(0);
                leg = createLegend();
                std::vector<TH1D*> hists_to_draw;
                double global_min = 1e30, global_max = -1e30;

                // 循环 1: 获取数据并计算极值
                for (size_t j = 0; j < isotopes_be.size(); ++j) {
                    TH1D* h = getRatio(path_noBkgCut + filename, "h_iss_ratio_" + isotopes_be[j], "h_mc_ratio_" + isotopes_be[j]);
                    if (!h) continue;
                    hists_to_draw.push_back(h);
                    std::pair<double, double> minmax = findMinMaxInRange(h, x_min, x_max);
                    if (minmax.first < global_min) global_min = minmax.first;
                    if (minmax.second > global_max) global_max = minmax.second;
                }

                // 计算最终 Y 轴范围
                double preset_min = y_ranges[i_range][0], preset_max = y_ranges[i_range][1];
                double final_min = preset_min, final_max = preset_max;
                if (global_min > preset_min && global_min < preset_max) final_min = global_min;
                if (global_max < preset_max && global_max > preset_min) final_max = global_max;
                if (final_min >= final_max) { final_min = preset_min; final_max = preset_max; } // 防止范围无效

                // 循环 2: 绘制
                bool firstHist = true;
                for (size_t j = 0; j < hists_to_draw.size(); ++j) {
                    TH1D* h = hists_to_draw[j];
                    styleHist(h, colors[j], marker_style);
                    h->GetYaxis()->SetRangeUser(final_min, final_max);
                    h->SetTitle("Bkg: B #rightarrow Be (No BkgCut)");
                    if (firstHist) { h->Draw("P"); firstHist = false; }
                    else { h->Draw("P SAME"); }
                    leg->AddEntry(h, labels_be[j], "p");
                }
                
                leg->Draw();
                c->SaveAs(plot_output_path + "B_to_Be_isotopes_range" + (i_range + 1) + "_NoBkgCut.png");
                delete leg;
                for (auto h : hists_to_draw) delete h; // 清理内存
            }

            // --- 1.2: With BkgCut (重构) ---
            {
                c->Clear();
                c->SetLogx(0); c->SetLogy(0);
                leg = createLegend();
                std::vector<TH1D*> hists_to_draw;
                double global_min = 1e30, global_max = -1e30;

                for (size_t j = 0; j < isotopes_be.size(); ++j) {
                    TH1D* h = getRatio(path_withBkgCut + filename, "h_iss_ratio_" + isotopes_be[j], "h_mc_ratio_" + isotopes_be[j]);
                    if (!h) continue;
                    hists_to_draw.push_back(h);
                    std::pair<double, double> minmax = findMinMaxInRange(h, x_min, x_max);
                    if (minmax.first < global_min) global_min = minmax.first;
                    if (minmax.second > global_max) global_max = minmax.second;
                }

                double preset_min = y_ranges[i_range][0], preset_max = y_ranges[i_range][1];
                double final_min = preset_min, final_max = preset_max;
                if (global_min > preset_min && global_min < preset_max) final_min = global_min;
                if (global_max < preset_max && global_max > preset_min) final_max = global_max;
                if (final_min >= final_max) { final_min = preset_min; final_max = preset_max; }

                bool firstHist = true;
                for (size_t j = 0; j < hists_to_draw.size(); ++j) {
                    TH1D* h = hists_to_draw[j];
                    styleHist(h, colors[j], marker_style);
                    h->GetYaxis()->SetRangeUser(final_min, final_max);
                    h->SetTitle("Bkg: B #rightarrow Be (With BkgCut)");
                    if (firstHist) { h->Draw("P"); firstHist = false; }
                    else { h->Draw("P SAME"); }
                    leg->AddEntry(h, labels_be[j], "p");
                }
                
                leg->Draw();
                c->SaveAs(plot_output_path + "B_to_Be_isotopes_range" + (i_range + 1) + "_WithBkgCut.png");
                delete leg;
                for (auto h : hists_to_draw) delete h;
            }
        }
    }

    // --- 组 2: C,N,O -> Beryllium (total) ---
    {
        std::vector<TString> sources_cno = {"Carbon", "Nitrogen", "Oxygen"};
        std::vector<TString> labels_cno = {"C #rightarrow Be (Total)", "N #rightarrow Be (Total)", "O #rightarrow Be (Total)"};

        for (int i_range = 0; i_range < 2; ++i_range) {
            // --- 2.1: No BkgCut (重构) ---
            {
                c->Clear();
                c->SetLogx(0); c->SetLogy(0);
                leg = createLegend();
                std::vector<TH1D*> hists_to_draw;
                double global_min = 1e30, global_max = -1e30;

                for (size_t j = 0; j < sources_cno.size(); ++j) {
                    TString filename = sources_cno[j] + "_to_Beryllium_L1Inner_Validation.root";
                    TH1D* h = getRatio(path_noBkgCut + filename, "h_iss_ratio_total", "h_mc_ratio_total");
                    if (!h) continue;
                    hists_to_draw.push_back(h);
                    std::pair<double, double> minmax = findMinMaxInRange(h, x_min, x_max);
                    if (minmax.first < global_min) global_min = minmax.first;
                    if (minmax.second > global_max) global_max = minmax.second;
                }

                double preset_min = y_ranges[i_range][0], preset_max = y_ranges[i_range][1];
                double final_min = preset_min, final_max = preset_max;
                if (global_min > preset_min && global_min < preset_max) final_min = global_min;
                if (global_max < preset_max && global_max > preset_min) final_max = global_max;
                if (final_min >= final_max) { final_min = preset_min; final_max = preset_max; }

                bool firstHist = true;
                for (size_t j = 0; j < hists_to_draw.size(); ++j) {
                    TH1D* h = hists_to_draw[j];
                    styleHist(h, colors[j], marker_style);
                    h->GetYaxis()->SetRangeUser(final_min, final_max);
                    h->SetTitle("Bkg: C,N,O #rightarrow Be (No BkgCut)");
                    if (firstHist) { h->Draw("P"); firstHist = false; }
                    else { h->Draw("P SAME"); }
                    leg->AddEntry(h, labels_cno[j], "p");
                }
                
                leg->Draw();
                c->SaveAs(plot_output_path + "CNO_to_Be_total_range" + (i_range + 1) + "_NoBkgCut.png");
                delete leg;
                for (auto h : hists_to_draw) delete h;
            }

            // --- 2.2: With BkgCut (重构) ---
            {
                c->Clear();
                c->SetLogx(0); c->SetLogy(0);
                leg = createLegend();
                std::vector<TH1D*> hists_to_draw;
                double global_min = 1e30, global_max = -1e30;

                for (size_t j = 0; j < sources_cno.size(); ++j) {
                    TString filename = sources_cno[j] + "_to_Beryllium_L1Inner_Validation.root";
                    TH1D* h = getRatio(path_withBkgCut + filename, "h_iss_ratio_total", "h_mc_ratio_total");
                    if (!h) continue;
                    hists_to_draw.push_back(h);
                    std::pair<double, double> minmax = findMinMaxInRange(h, x_min, x_max);
                    if (minmax.first < global_min) global_min = minmax.first;
                    if (minmax.second > global_max) global_max = minmax.second;
                }

                double preset_min = y_ranges[i_range][0], preset_max = y_ranges[i_range][1];
                double final_min = preset_min, final_max = preset_max;
                if (global_min > preset_min && global_min < preset_max) final_min = global_min;
                if (global_max < preset_max && global_max > preset_min) final_max = global_max;
                if (final_min >= final_max) { final_min = preset_min; final_max = preset_max; }

                bool firstHist = true;
                for (size_t j = 0; j < hists_to_draw.size(); ++j) {
                    TH1D* h = hists_to_draw[j];
                    styleHist(h, colors[j], marker_style);
                    h->GetYaxis()->SetRangeUser(final_min, final_max);
                    h->SetTitle("Bkg: C,N,O #rightarrow Be (With BkgCut)");
                    if (firstHist) { h->Draw("P"); firstHist = false; }
                    else { h->Draw("P SAME"); }
                    leg->AddEntry(h, labels_cno[j], "p");
                }
                
                leg->Draw();
                c->SaveAs(plot_output_path + "CNO_to_Be_total_range" + (i_range + 1) + "_WithBkgCut.png");
                delete leg;
                for (auto h : hists_to_draw) delete h;
            }
        }
    }

    // --- 组 3: Carbon -> Boron (isotopes) ---
    {
        std::vector<TString> isotopes_b = {"B10", "B11", "total"};
        std::vector<TString> labels_b = {"^{10}B", "^{11}B", "Total (C#rightarrowB)"};
        TString filename = "Carbon_to_Boron_L1Inner_Validation.root";

        for (int i_range = 0; i_range < 2; ++i_range) {
            // --- 3.1: No BkgCut (重构) ---
            {
                c->Clear();
                c->SetLogx(0); c->SetLogy(0);
                leg = createLegend();
                std::vector<TH1D*> hists_to_draw;
                double global_min = 1e30, global_max = -1e30;

                for (size_t j = 0; j < isotopes_b.size(); ++j) {
                    TH1D* h = getRatio(path_noBkgCut + filename, "h_iss_ratio_" + isotopes_b[j], "h_mc_ratio_" + isotopes_b[j]);
                    if (!h) continue;
                    hists_to_draw.push_back(h);
                    std::pair<double, double> minmax = findMinMaxInRange(h, x_min, x_max);
                    if (minmax.first < global_min) global_min = minmax.first;
                    if (minmax.second > global_max) global_max = minmax.second;
                }

                double preset_min = y_ranges[i_range][0], preset_max = y_ranges[i_range][1];
                double final_min = preset_min, final_max = preset_max;
                if (global_min > preset_min && global_min < preset_max) final_min = global_min;
                if (global_max < preset_max && global_max > preset_min) final_max = global_max;
                if (final_min >= final_max) { final_min = preset_min; final_max = preset_max; }

                bool firstHist = true;
                for (size_t j = 0; j < hists_to_draw.size(); ++j) {
                    TH1D* h = hists_to_draw[j];
                    styleHist(h, colors[j], marker_style);
                    h->GetYaxis()->SetRangeUser(final_min, final_max);
                    h->SetTitle("Bkg: C #rightarrow B (No BkgCut)");
                    if (firstHist) { h->Draw("P"); firstHist = false; }
                    else { h->Draw("P SAME"); }
                    leg->AddEntry(h, labels_b[j], "p");
                }
                
                leg->Draw();
                c->SaveAs(plot_output_path + "C_to_B_isotopes_range" + (i_range + 1) + "_NoBkgCut.png");
                delete leg;
                for (auto h : hists_to_draw) delete h;
            }

            // --- 3.2: With BkgCut (重构) ---
            {
                c->Clear();
                c->SetLogx(0); c->SetLogy(0);
                leg = createLegend();
                std::vector<TH1D*> hists_to_draw;
                double global_min = 1e30, global_max = -1e30;

                for (size_t j = 0; j < isotopes_b.size(); ++j) {
                    TH1D* h = getRatio(path_withBkgCut + filename, "h_iss_ratio_" + isotopes_b[j], "h_mc_ratio_" + isotopes_b[j]);
                    if (!h) continue;
                    hists_to_draw.push_back(h);
                    std::pair<double, double> minmax = findMinMaxInRange(h, x_min, x_max);
                    if (minmax.first < global_min) global_min = minmax.first;
                    if (minmax.second > global_max) global_max = minmax.second;
                }

                double preset_min = y_ranges[i_range][0], preset_max = y_ranges[i_range][1];
                double final_min = preset_min, final_max = preset_max;
                if (global_min > preset_min && global_min < preset_max) final_min = global_min;
                if (global_max < preset_max && global_max > preset_min) final_max = global_max;
                if (final_min >= final_max) { final_min = preset_min; final_max = preset_max; }

                bool firstHist = true;
                for (size_t j = 0; j < hists_to_draw.size(); ++j) {
                    TH1D* h = hists_to_draw[j];
                    styleHist(h, colors[j], marker_style);
                    h->GetYaxis()->SetRangeUser(final_min, final_max);
                    h->SetTitle("Bkg: C #rightarrow B (With BkgCut)");
                    if (firstHist) { h->Draw("P"); firstHist = false; }
                    else { h->Draw("P SAME"); }
                    leg->AddEntry(h, labels_b[j], "p");
                }
                
                leg->Draw();
                c->SaveAs(plot_output_path + "C_to_B_isotopes_range" + (i_range + 1) + "_WithBkgCut.png");
                delete leg;
                for (auto h : hists_to_draw) delete h;
            }
        }
    }

    // --- 组 4: N,O -> Boron (total) ---
    {
        std::vector<TString> sources_no = {"Nitrogen", "Oxygen"};
        std::vector<TString> labels_no = {"N #rightarrow B (Total)", "O #rightarrow B (Total)"};

        for (int i_range = 0; i_range < 2; ++i_range) {
            // --- 4.1: No BkgCut (重构) ---
            {
                c->Clear();
                c->SetLogx(0); c->SetLogy(0);
                leg = createLegend();
                std::vector<TH1D*> hists_to_draw;
                double global_min = 1e30, global_max = -1e30;

                for (size_t j = 0; j < sources_no.size(); ++j) {
                    TString filename = sources_no[j] + "_to_Boron_L1Inner_Validation.root";
                    TH1D* h = getRatio(path_noBkgCut + filename, "h_iss_ratio_total", "h_mc_ratio_total");
                    if (!h) continue;
                    hists_to_draw.push_back(h);
                    std::pair<double, double> minmax = findMinMaxInRange(h, x_min, x_max);
                    if (minmax.first < global_min) global_min = minmax.first;
                    if (minmax.second > global_max) global_max = minmax.second;
                }

                double preset_min = y_ranges[i_range][0], preset_max = y_ranges[i_range][1];
                double final_min = preset_min, final_max = preset_max;
                if (global_min > preset_min && global_min < preset_max) final_min = global_min;
                if (global_max < preset_max && global_max > preset_min) final_max = global_max;
                if (final_min >= final_max) { final_min = preset_min; final_max = preset_max; }

                bool firstHist = true;
                for (size_t j = 0; j < hists_to_draw.size(); ++j) {
                    TH1D* h = hists_to_draw[j];
                    styleHist(h, colors[j], marker_style);
                    h->GetYaxis()->SetRangeUser(final_min, final_max);
                    h->SetTitle("Bkg: N,O #rightarrow B (No BkgCut)");
                    if (firstHist) { h->Draw("P"); firstHist = false; }
                    else { h->Draw("P SAME"); }
                    leg->AddEntry(h, labels_no[j], "p");
                }
                
                leg->Draw();
                c->SaveAs(plot_output_path + "NO_to_B_total_range" + (i_range + 1) + "_NoBkgCut.png");
                delete leg;
                for (auto h : hists_to_draw) delete h;
            }

            // --- 4.2: With BkgCut (重构) ---
            {
                c->Clear();
                c->SetLogx(0); c->SetLogy(0);
                leg = createLegend();
                std::vector<TH1D*> hists_to_draw;
                double global_min = 1e30, global_max = -1e30;

                for (size_t j = 0; j < sources_no.size(); ++j) {
                    TString filename = sources_no[j] + "_to_Boron_L1Inner_Validation.root";
                    TH1D* h = getRatio(path_withBkgCut + filename, "h_iss_ratio_total", "h_mc_ratio_total");
                    if (!h) continue;
                    hists_to_draw.push_back(h);
                    std::pair<double, double> minmax = findMinMaxInRange(h, x_min, x_max);
                    if (minmax.first < global_min) global_min = minmax.first;
                    if (minmax.second > global_max) global_max = minmax.second;
                }

                double preset_min = y_ranges[i_range][0], preset_max = y_ranges[i_range][1];
                double final_min = preset_min, final_max = preset_max;
                if (global_min > preset_min && global_min < preset_max) final_min = global_min;
                if (global_max < preset_max && global_max > preset_min) final_max = global_max;
                if (final_min >= final_max) { final_min = preset_min; final_max = preset_max; }

                bool firstHist = true;
                for (size_t j = 0; j < hists_to_draw.size(); ++j) {
                    TH1D* h = hists_to_draw[j];
                    styleHist(h, colors[j], marker_style);
                    h->GetYaxis()->SetRangeUser(final_min, final_max);
                    h->SetTitle("Bkg: N,O #rightarrow B (With BkgCut)");
                    if (firstHist) { h->Draw("P"); firstHist = false; }
                    else { h->Draw("P SAME"); }
                    leg->AddEntry(h, labels_no[j], "p");
                }
                
                leg->Draw();
                c->SaveAs(plot_output_path + "NO_to_B_total_range" + (i_range + 1) + "_WithBkgCut.png");
                delete leg;
                for (auto h : hists_to_draw) delete h;
            }
        }
    }

    delete c;
    gROOT->SetBatch(kFALSE);
    std::cout << "--- 所有 16 张图已生成完毕并保存至: " << plot_output_path << " (Y轴已自动缩放) ---" << std::endl;
}