#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TF1.h>
#include <TBox.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>

// 全局常量定义
const TString path_noBkgCut = "/eos/user/z/zixuan/Isotope/BkgValid/";
const TString path_withBkgCut = "/eos/user/z/zixuan/Isotope/BkgValid_B/";
const TString plot_output_path = "/eos/user/z/zixuan/Isotope/BkgValid/Plots/ratio_final_plots_fully_corrected/";

const int marker_style = 20;
const int data_color = kBlue;
const int fit_color_linear = kMagenta;
const int fit_color_const = kOrange + 1;

const double x_min = 0.1;
const double x_max = 21.5;
const double fit_linear_min = 0.3;
const double fit_linear_max = 5.0;
const double fit_const_min = 5.0;
const double fit_const_max = 21.5;

/**
 * @brief 计算68%置信区间的不确定度（针对一次函数）
 */
double calculate68PercentUncertainty_Linear(TH1D* hist, TF1* fit_func, double fit_min, double fit_max) {
    if (!hist || !fit_func) return 0.0;
    
    // 收集拟合范围内的所有有效数据点
    int total_points = 0;
    
    for (int bin = hist->GetXaxis()->FindBin(fit_min); bin <= hist->GetXaxis()->FindBin(fit_max); ++bin) {
        double x = hist->GetBinCenter(bin);
        double y = hist->GetBinContent(bin);
        double y_err = hist->GetBinError(bin);
        
        if (y != 0 && !std::isnan(y) && !std::isinf(y) && y_err > 0) {
            total_points++;
        }
    }
    
    if (total_points == 0) {
        std::cout << "Warning: No valid points for uncertainty calculation" << std::endl;
        return 0.0;
    }
    
    std::cout << "Total valid points in range: " << total_points << std::endl;
    
    // 从0开始递增常数，直到68%的点落在 f(x)±delta 范围内
    double delta = 0.0;
    double step = 0.01;
    int target_points = static_cast<int>(std::ceil(0.68 * total_points));
    
    std::cout << "Target: " << target_points << " points out of " << total_points << " (68%)" << std::endl;
    
    while (delta < 10.0) {
        int points_within = 0;
        
        for (int bin = hist->GetXaxis()->FindBin(fit_min); bin <= hist->GetXaxis()->FindBin(fit_max); ++bin) {
            double x = hist->GetBinCenter(bin);
            double y = hist->GetBinContent(bin);
            double y_err = hist->GetBinError(bin);
            
            if (y != 0 && !std::isnan(y) && !std::isinf(y) && y_err > 0) {
                double fit_val = fit_func->Eval(x);
                double upper = fit_val + delta;
                double lower = fit_val - delta;
                
                if (y >= lower && y <= upper) {
                    points_within++;
                }
            }
        }
        
        if (points_within >= target_points) {
            std::cout << "Found delta = " << delta << " with " << points_within << " points within band" << std::endl;
            return delta;
        }
        
        delta += step;
    }
    
    std::cout << "Warning: Could not find suitable uncertainty band (reached limit)" << std::endl;
    return delta;
}

/**
 * @brief 计算68%置信区间的不确定度（针对常数函数）
 */
double calculate68PercentUncertainty_Const(TH1D* hist, TF1* fit_func, double fit_min, double fit_max) {
    if (!hist || !fit_func) return 0.0;
    
    int total_points = 0;
    
    for (int bin = hist->GetXaxis()->FindBin(fit_min); bin <= hist->GetXaxis()->FindBin(fit_max); ++bin) {
        double x = hist->GetBinCenter(bin);
        double y = hist->GetBinContent(bin);
        double y_err = hist->GetBinError(bin);
        
        if (y != 0 && !std::isnan(y) && !std::isinf(y) && y_err > 0) {
            total_points++;
        }
    }
    
    if (total_points == 0) {
        std::cout << "Warning: No valid points for uncertainty calculation" << std::endl;
        return 0.0;
    }
    
    std::cout << "Total valid points in range: " << total_points << std::endl;
    
    double delta = 0.0;
    double step = 0.01;
    int target_points = static_cast<int>(std::ceil(0.68 * total_points));
    
    std::cout << "Target: " << target_points << " points out of " << total_points << " (68%)" << std::endl;
    
    while (delta < 10.0) {
        int points_within = 0;
        
        for (int bin = hist->GetXaxis()->FindBin(fit_min); bin <= hist->GetXaxis()->FindBin(fit_max); ++bin) {
            double x = hist->GetBinCenter(bin);
            double y = hist->GetBinContent(bin);
            double y_err = hist->GetBinError(bin);
            
            if (y != 0 && !std::isnan(y) && !std::isinf(y) && y_err > 0) {
                double fit_val = fit_func->Eval(x);
                double upper = fit_val + delta;
                double lower = fit_val - delta;
                
                if (y >= lower && y <= upper) {
                    points_within++;
                }
            }
        }
        
        if (points_within >= target_points) {
            std::cout << "Found delta = " << delta << " with " << points_within << " points within band" << std::endl;
            return delta;
        }
        
        delta += step;
    }
    
    std::cout << "Warning: Could not find suitable uncertainty band (reached limit)" << std::endl;
    return delta;
}

/**
 * @brief 从文件中获取两个本身就是比例的直方图，计算它们的比值
 */
TH1D* getRatio(TString fp, TString his_iss, TString his_mc) {
    TFile* f = TFile::Open(fp);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open file: " << fp << std::endl;
        return nullptr;
    }
    TH1* h_iss = (TH1*)f->Get(his_iss);
    TH1* h_mc = (TH1*)f->Get(his_mc);
    if (!h_iss || !h_mc) {
        std::cerr << "Error: Cannot find histograms '" << his_iss << "' or '" << his_mc << "' in file: " << fp << std::endl;
        f->Close();
        delete f;
        return nullptr;
    }

    h_iss->Sumw2(kTRUE);
    h_mc->Sumw2(kTRUE);

    TH1D* h_ratio = (TH1D*)h_iss->Clone(TString(his_iss) + "_ratio");
    h_ratio->SetDirectory(0);
    h_ratio->Sumw2(kTRUE);

    h_ratio->Divide(h_iss, h_mc, 1.0, 1.0, "");

    f->Close();
    delete f;
    return h_ratio;
}

/**
 * @brief 在指定x范围内查找直方图的y值最小和最大值
 */
std::pair<double, double> findMinMaxInRange(TH1* hist, double xmin, double xmax) {
    if (!hist) return {0, 2};
    int first = hist->GetXaxis()->FindBin(xmin);
    int last = hist->GetXaxis()->FindBin(xmax);
    double local_min = 1e30, local_max = -1e30;
    bool found = false;
    for (int i = first; i <= last; ++i) {
        double val = hist->GetBinContent(i);
        if (std::isnan(val) || std::isinf(val) || val == 0) continue;
        
        local_min = std::min(local_min, val);
        local_max = std::max(local_max, val);
        found = true;
    }
    if (!found) return {0, 2};

    double margin = (local_max - local_min) * 0.2;
    if (margin < 0.1) margin = 0.1;
    double y_min_final = local_min - margin;
    double y_max_final = local_max + margin;

    if (std::abs(y_max_final - y_min_final) < 1e-6) {
        y_max_final = y_min_final + 0.5;
    }

    return {y_min_final, y_max_final};
}

/**
 * @brief 创建一个标准化的图例
 */
TLegend* createLegend() {
    TLegend* leg = new TLegend(0.15, 0.72, 0.7, 0.95);
    leg->SetFillStyle(0);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.028);
    return leg;
}

/**
 * @brief 设置直方图的通用样式
 */
void styleHist(TH1D* h) {
    if (!h) return;
    h->SetLineColor(data_color);
    h->SetMarkerColor(data_color);
    h->SetMarkerStyle(marker_style);
    h->SetMarkerSize(1.0);
    h->GetXaxis()->SetRangeUser(x_min, x_max);
    h->GetXaxis()->SetTitle("Measured Ek/n [GeV]");
    h->GetYaxis()->SetTitle("ISS/MC");
    h->SetStats(0);
}

/**
 * @brief 绘制单个比值图，进行双段拟合并计算68%置信区间
 */
void drawSinglePlot(TCanvas* c, TH1D* h, const TString& title, const TString& output_fn)
{
    if (!h) {
        std::cerr << "Warning: Histogram for " << output_fn << " is null. Skipping." << std::endl;
        return;
    }
    c->Clear();
    
    std::pair<double, double> y_range = findMinMaxInRange(h, x_min, x_max);
    double ymin = y_range.first;
    double ymax = y_range.second;

    styleHist(h);
    h->GetYaxis()->SetRangeUser(ymin, ymax);
    h->SetTitle("");

    std::cout << "\n--- Processing: " << output_fn << " ---" << std::endl;
    
    h->Draw("E1");

    TLegend* leg = createLegend();
    leg->AddEntry(h, title, "pe");

    // ====== 一次函数拟合 (0.3 - 5 GeV) ======
    bool hasDataLinear = false;
    for (int bin = h->GetXaxis()->FindBin(fit_linear_min); bin <= h->GetXaxis()->FindBin(fit_linear_max); ++bin) {
        if (h->GetBinContent(bin) != 0) {
            hasDataLinear = true;
            break;
        }
    }

    TF1* fit_linear = nullptr;
    double delta_linear = 0.0;
    
    if (hasDataLinear) {
        std::cout << "Linear fit range [" << fit_linear_min << ", " << fit_linear_max << "]" << std::endl;
        
        TString fit_name_linear = TString(h->GetName()) + "_fit_linear";
        fit_linear = new TF1(fit_name_linear, "pol1", fit_linear_min, fit_linear_max);
        fit_linear->SetLineColor(fit_color_linear);
        fit_linear->SetLineWidth(2);
        
        TFitResultPtr r_linear = h->Fit(fit_linear, "SRQ");
        
        if (r_linear->IsValid() && r_linear->Status() == 0) {
            double p0 = fit_linear->GetParameter(0);
            double p1 = fit_linear->GetParameter(1);
            double chi2 = r_linear->Chi2();
            int ndf = r_linear->Ndf();
            
            std::cout << "Linear fit: p0 = " << p0 << ", p1 = " << p1 << std::endl;
            std::cout << "Chi2/NDF = " << chi2 << " / " << ndf << std::endl;
            
            // 计算68%置信区间
            delta_linear = calculate68PercentUncertainty_Linear(h, fit_linear, fit_linear_min, fit_linear_max);
            
            fit_linear->Draw("SAME");
            
            // 绘制一次函数的不确定度带（使用TGraph）
            const int n_points = 100;
            double x_arr[n_points * 2];
            double y_arr[n_points * 2];
            
            for (int i = 0; i < n_points; ++i) {
                double x = fit_linear_min + i * (fit_linear_max - fit_linear_min) / (n_points - 1);
                double y_center = fit_linear->Eval(x);
                
                // 上边界
                x_arr[i] = x;
                y_arr[i] = y_center + delta_linear;
                
                // 下边界（反向填充）
                x_arr[2 * n_points - 1 - i] = x;
                y_arr[2 * n_points - 1 - i] = y_center - delta_linear;
            }
            
            TGraph* uncertainty_band_linear = new TGraph(2 * n_points, x_arr, y_arr);
            uncertainty_band_linear->SetFillColorAlpha(fit_color_linear, 0.25);
            uncertainty_band_linear->SetFillStyle(1001);
            uncertainty_band_linear->SetLineColor(fit_color_linear - 7);
            uncertainty_band_linear->SetLineWidth(1);
            uncertainty_band_linear->Draw("F SAME");
            
            fit_linear->Draw("SAME");  // 重新画拟合线，确保在前面
            
            TString fit_label = TString::Format("Linear fit (0.3-5 GeV): %.2f + %.2f#timesE #pm %.2f", p0, p1, delta_linear);
            leg->AddEntry(uncertainty_band_linear, fit_label, "lf");
            
            if (ndf > 0) {
                TString chi2_label = TString::Format("#chi^{2}/NDF = %.1f / %d = %.2f", chi2, ndf, chi2 / ndf);
                leg->AddEntry((TObject*)0, chi2_label, "");
            }
        } else {
            std::cout << "Warning: Linear fit failed. Status: " << r_linear->Status() << std::endl;
        }
    }

    // ====== 常数拟合 (5 - 21.5 GeV) ======
    bool hasDataConst = false;
    for (int bin = h->GetXaxis()->FindBin(fit_const_min); bin <= h->GetXaxis()->FindBin(fit_const_max); ++bin) {
        if (h->GetBinContent(bin) != 0) {
            hasDataConst = true;
            break;
        }
    }

    TF1* fit_const = nullptr;
    double delta_const = 0.0;
    
    if (hasDataConst) {
        std::cout << "Constant fit range [" << fit_const_min << ", " << fit_const_max << "]" << std::endl;
        
        TString fit_name_const = TString(h->GetName()) + "_fit_const";
        fit_const = new TF1(fit_name_const, "pol0", fit_const_min, fit_const_max);
        fit_const->SetLineColor(fit_color_const);
        fit_const->SetLineWidth(2);
        
        TFitResultPtr r_const = h->Fit(fit_const, "SRQ+");
        
        if (r_const->IsValid() && r_const->Status() == 0) {
            double c = fit_const->GetParameter(0);
            double chi2 = r_const->Chi2();
            int ndf = r_const->Ndf();
            
            std::cout << "Constant fit: c = " << c << std::endl;
            std::cout << "Chi2/NDF = " << chi2 << " / " << ndf << std::endl;
            
            // 计算68%置信区间
            delta_const = calculate68PercentUncertainty_Const(h, fit_const, fit_const_min, fit_const_max);
            
            fit_const->Draw("SAME");
            
            // 绘制常数的不确定度带（矩形）
            TBox* error_box_const = new TBox(fit_const_min, c - delta_const, 
                                            fit_const_max, c + delta_const);
            error_box_const->SetFillColorAlpha(fit_color_const, 0.25);
            error_box_const->SetFillStyle(1001);
            error_box_const->SetLineColor(fit_color_const - 7);
            error_box_const->SetLineWidth(1);
            error_box_const->Draw("SAME");
            
            fit_const->Draw("SAME");  // 重新画拟合线
            
            TString fit_label = TString::Format("Constant fit (5-21.5 GeV): %.2f #pm %.2f", c, delta_const);
            leg->AddEntry(error_box_const, fit_label, "lf");
            
            if (ndf > 0) {
                TString chi2_label = TString::Format("#chi^{2}/NDF = %.1f / %d = %.2f", chi2, ndf, chi2 / ndf);
                leg->AddEntry((TObject*)0, chi2_label, "");
            }
        } else {
            std::cout << "Warning: Constant fit failed. Status: " << r_const->Status() << std::endl;
        }
    }

    leg->Draw();
    c->SaveAs(plot_output_path + output_fn + ".png");
}

/**
 * @brief 处理来自单个文件的多个碎片的绘图任务
 */
void processGroup(TCanvas* c, const TString& file_path_base, const TString& file_name_base,
                  const std::vector<TString>& frag_names, const std::vector<TString>& labels,
                  const TString& title_pfx, const TString& output_pfx)
{
    TString cut_sfx_label = (file_path_base.Contains("BkgValid_B")) ? " (With BkgCut)" : " (No BkgCut)";
    TString cut_sfx_fn = (file_path_base.Contains("BkgValid_B")) ? "WithBkgCut" : "NoBkgCut";

    for (size_t j = 0; j < frag_names.size(); ++j) {
        TH1D* h = getRatio(file_path_base + file_name_base, "h_iss_ratio_" + frag_names[j], "h_mc_ratio_" + frag_names[j]);
        
        TString current_title = title_pfx + " (" + labels[j] + ")" + cut_sfx_label;
        TString output_fn = output_pfx + "_" + frag_names[j] + "_" + cut_sfx_fn;
        
        drawSinglePlot(c, h, current_title, output_fn);
        
        if (h) delete h;
    }
}

/**
 * @brief 处理单个源核到单个碎片目标的绘图任务
 */
void processSingleSourceToFragPlot(TCanvas* c, const TString& file_path_base, 
                                   const TString& src_name, const TString& frag_name, 
                                   const TString& label, const TString& output_pfx)
{
    TString cut_sfx_label = (file_path_base.Contains("BkgValid_B")) ? " (With BkgCut)" : " (No BkgCut)";
    TString cut_sfx_fn = (file_path_base.Contains("BkgValid_B")) ? "WithBkgCut" : "NoBkgCut";

    TString current_fn = src_name + "_to_" + frag_name + "_L1Inner_Validation.root";
    TH1D* h = getRatio(file_path_base + current_fn, "h_iss_ratio_total", "h_mc_ratio_total");

    TString current_title = "" + label + cut_sfx_label;
    TString output_fn = output_pfx + "_" + cut_sfx_fn;
    
    drawSinglePlot(c, h, current_title, output_fn);

    if (h) delete h;
}

/**
 * @brief 主函数
 */
void plotR() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gSystem->Exec("mkdir -p " + plot_output_path);
    TCanvas* c = new TCanvas("c", "Bkg Ratios", 800, 600);
    c->SetTopMargin(0.05);
    c->SetGrid();

    // --- Group 1: B -> Be ---
    std::vector<TString> is_be = {"Be7", "Be9", "Be10", "Be"};
    std::vector<TString> lb_be = {"^{7}Be", "^{9}Be", "^{10}Be", "Be"};
    TString fn_be = "Boron_to_Beryllium_L1Inner_Validation.root";
    processGroup(c, path_noBkgCut, fn_be, is_be, lb_be, "B #rightarrow Be", "B_to_Be");
    processGroup(c, path_withBkgCut, fn_be, is_be, lb_be, "B #rightarrow Be", "B_to_Be");

    // --- Group 2: C -> B ---
    std::vector<TString> is_b = {"B10", "B11", "B"};
    std::vector<TString> lb_b = {"^{10}B", "^{11}B", "B"};
    TString fn_b = "Carbon_to_Boron_L1Inner_Validation.root";
    processGroup(c, path_noBkgCut, fn_b, is_b, lb_b, "C #rightarrow B", "C_to_B");
    processGroup(c, path_withBkgCut, fn_b, is_b, lb_b, "C #rightarrow B", "C_to_B");

    // --- Group 3: Heavier nuclei fragmentation ---
    processSingleSourceToFragPlot(c, path_noBkgCut, "Carbon", "Beryllium", "C #rightarrow Be", "C_to_Be_total");
    processSingleSourceToFragPlot(c, path_withBkgCut, "Carbon", "Beryllium", "C #rightarrow Be", "C_to_Be_total");
    
    processSingleSourceToFragPlot(c, path_noBkgCut, "Nitrogen", "Beryllium", "N #rightarrow Be", "N_to_Be_total");
    processSingleSourceToFragPlot(c, path_withBkgCut, "Nitrogen", "Beryllium", "N #rightarrow Be", "N_to_Be_total");

    processSingleSourceToFragPlot(c, path_noBkgCut, "Oxygen", "Beryllium", "O #rightarrow Be", "O_to_Be_total");
    processSingleSourceToFragPlot(c, path_withBkgCut, "Oxygen", "Beryllium", "O #rightarrow Be", "O_to_Be_total");

    processSingleSourceToFragPlot(c, path_noBkgCut, "Nitrogen", "Boron", "N #rightarrow B", "N_to_B_total");
    processSingleSourceToFragPlot(c, path_withBkgCut, "Nitrogen", "Boron", "N #rightarrow B", "N_to_B_total");

    processSingleSourceToFragPlot(c, path_noBkgCut, "Oxygen", "Boron", "O #rightarrow B", "O_to_B_total");
    processSingleSourceToFragPlot(c, path_withBkgCut, "Oxygen", "Boron", "O #rightarrow B", "O_to_B_total");

    delete c;
    gROOT->SetBatch(kFALSE);
    
    std::cout << "\n=== 完成！已生成所有图 (双段拟合 + 68%置信区间) ===" << std::endl;
}