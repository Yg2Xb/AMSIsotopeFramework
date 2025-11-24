#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iomanip>
#include <sstream>
#include <TStyle.h>
#include <TLine.h>
#include "../Tool.h" // 确保 DoGausPlusAsymGausFit 函数在此头文件中

using namespace std;
using namespace AMS_Iso;

// --- 结构体定义 (更新：为 LR 和 RR 增加误差字段) ---
struct HistInfo {
    string y_axis_label;
    string output_suffix;
    string title_description;
};

struct FitResult {
    // 核心高斯的参数 (用于最终存储和绘图)
    double mean;
    double mean_err;
    double sigma;
    double sigma_err;
    // 附加的拟合信息 (用于PDF绘图/调试)
    double chi2;
    double ndf;
    double LR; 
    double LR_err; // 新增 LR 误差
    double RR; 
    double RR_err; // 新增 RR 误差
};

using H4ResultsMap = map<string, FitResult>;

// --- 常量定义 ---
const vector<double> H5_RIG_BINS_EDGES = {30.0, 50.0, 80.0, 120.0, 160.0, 240.0};
const vector<string> H5_RIG_LABELS = {"30-50 GV", "50-80 GV", "80-120 GV", "120-160 GV", "160-240 GV"};
const int N_RIG_BINS = H5_RIG_BINS_EDGES.size() - 1;

const vector<string> H5_CHARGE_LABELS = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
const vector<double> H5_CHARGE_BINS_EDGES = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
const int N_CHARGE_BINS = H5_CHARGE_BINS_EDGES.size() - 1;

// --- 函数实现：查找拟合范围 ---
void findFitRange(TH1* hist, double coverage, double center_x, double& x_min, double& x_max) {
    if (!hist || hist->GetEntries() == 0) return;
    double total_integral = hist->Integral();
    if (total_integral <= 0) return;
    
    // 优先使用用户提供的中心点，否则使用最大 bin
    int center_bin = hist->GetXaxis()->FindFixBin(center_x);
    if (center_bin < 1 || center_bin > hist->GetNbinsX()) {
        center_bin = hist->GetMaximumBin(); 
        center_x = hist->GetXaxis()->GetBinCenter(center_bin);
    }
    
    double required_integral = total_integral * coverage;
    double current_integral = hist->GetBinContent(center_bin);
    int low_bin = center_bin - 1;
    int high_bin = center_bin + 1;
    x_min = hist->GetXaxis()->GetBinLowEdge(center_bin);
    x_max = hist->GetXaxis()->GetBinUpEdge(center_bin);

    // 从中心 bin 向两侧扩展，直到覆盖所需积分比例
    while (current_integral < required_integral) {
        bool extended = false;
        double content_low = (low_bin >= 1) ? hist->GetBinContent(low_bin) : 0;
        double content_high = (high_bin <= hist->GetNbinsX()) ? hist->GetBinContent(high_bin) : 0;

        if (low_bin >= 1 && high_bin <= hist->GetNbinsX()) {
            current_integral += content_low + content_high;
            x_min = hist->GetXaxis()->GetBinLowEdge(low_bin);
            x_max = hist->GetXaxis()->GetBinUpEdge(high_bin);
            low_bin--;
            high_bin++;
            extended = true;
        } else if (low_bin >= 1) {
            current_integral += content_low;
            x_min = hist->GetXaxis()->GetBinLowEdge(low_bin);
            low_bin--;
            extended = true;
        } else if (high_bin <= hist->GetNbinsX()) {
            current_integral += content_high;
            x_max = hist->GetXaxis()->GetBinUpEdge(high_bin);
            high_bin++;
            extended = true;
        }
        if (!extended) break;
    }
    
    // 确保拟合范围关于最终确定的中心点对称
    double final_center = (x_min + x_max) / 2.0; // 使用当前找到的范围的中心作为新的中心点
    double half_range = max(abs(x_max - final_center), abs(x_min - final_center));
    x_min = final_center - half_range;
    x_max = final_center + half_range;
}

// --- 函数实现：获取全局 Y 轴范围 ---
void getGlobalYRange(const vector<TGraphErrors*>& graphs, double x_min, double x_max, double& y_min, double& y_max) {
    y_min = 1e10;
    y_max = -1e10;
    bool found_point = false;
    for (const auto& g : graphs) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            if (x >= x_min && x <= x_max) {
                double y_low = y - g->GetErrorY(i);
                double y_high = y + g->GetErrorY(i);
                y_min = min(y_min, y_low);
                y_max = max(y_max, y_high);
                found_point = true;
            }
        }
    }
    if (found_point) {
        double range = y_max - y_min;
        if (range == 0.0) {
            range = abs(y_min * 0.1);
            if (range == 0) range = 0.01;
        }
        double buffer = range * 0.10;
        y_min -= buffer;
        y_max += buffer;
    } else {
        y_min = -0.01;
        y_max = 0.01;
    }
}

// --- 函数实现：获取直方图信息 ---
HistInfo getHistInfo(const string& suffix) {
    HistInfo info;
    if (suffix == "ID_H5a") { info = {"Rigidity [GV]", "Rigidity", "NaF-Tracker #Delta(1/#beta)"}; }
    else if (suffix == "ID_H5b") { info = {"Rigidity [GV]", "Rigidity", "AGL-Tracker #Delta(1/#beta)"}; }
    else if (suffix == "ID_H4a") { info = {"Rigidity [GV]", "Rigidity", "NaF 1/#beta (Rig > 80GV)"}; }
    else if (suffix == "ID_H4b") { info = {"Rigidity [GV]", "Rigidity", "AGL 1/#beta (Rig > 150GV)"}; }
    else { info = {"Y Variable", "YVariable", "Unknown Delta Beta"}; }
    return info;
}

// --- 函数实现：绘制和保存图表 (修改1: 移除 logx 选项) ---
void DrawAndSaveGraphs(vector<TGraphErrors*>& graphs, const vector<string>& legend_labels,
                         const string& x_axis_title, const string& y_axis_title,
                         const string& canvas_title_prefix, const string& output_dir,
                         const string& output_name_base, bool logx, TFile* save_to_root,
                         const string& graph_name_suffix) {
    
    if (graphs.empty()) return;
    string canvas_name = "c_" + output_name_base + graph_name_suffix;
    TCanvas* c = new TCanvas(canvas_name.c_str(), (canvas_title_prefix + graph_name_suffix).c_str(), 700, 500);
    c->SetGrid();
    // 原始代码: if (logx) c->SetLogx(1); // 移除或忽略logx设置

    // 假设所有图表的 X 轴范围相同，使用第一个图表的范围
    double x_min = graphs[0]->GetXaxis()->GetXmin();
    double x_max = graphs[0]->GetXaxis()->GetXmax();
    double y_min, y_max;
    getGlobalYRange(graphs, x_min, x_max, y_min, y_max);

    TLegend* leg = new TLegend(0.8, 0.8, 0.99, 0.99);
    leg->SetFillStyle(0);
    leg->SetBorderSize(1);

    for (size_t i = 0; i < graphs.size(); ++i) {
        string draw_opt = (i == 0) ? "APZ" : "PZ same";
        graphs[i]->Draw(draw_opt.c_str());
        
        if (i == 0) {
            graphs[i]->SetTitle((canvas_title_prefix + " vs " + x_axis_title).c_str());
            graphs[i]->GetXaxis()->SetTitle(x_axis_title.c_str());
            graphs[i]->GetYaxis()->SetTitle(y_axis_title.c_str());
            graphs[i]->GetXaxis()->SetRangeUser(x_min, x_max);
            graphs[i]->GetYaxis()->SetRangeUser(y_min, y_max);
        }

        if (i < legend_labels.size()) {
            leg->AddEntry(graphs[i], legend_labels[i].c_str(), "p");
        }

        if (save_to_root) {
            save_to_root->cd();
            graphs[i]->SetName(("g_" + output_name_base + graph_name_suffix + "_" + to_string(i)).c_str());
            graphs[i]->Write();
        }
    }

    leg->Draw();
    c->SaveAs((output_dir + output_name_base + graph_name_suffix + ".png").c_str());
    delete leg;
    delete c;
}

// --- 函数实现：绘制 H4 拟合结果 (图例标签改为 "ISS") ---
void DrawH4Results(const H4ResultsMap& results, const string& output_dir, const HistInfo& info,
                    const string& output_name_base, const string& x_axis_label,
                    const vector<string>& particles, TFile* save_to_root) {
    
    TGraphErrors* g_mean_z = new TGraphErrors();
    TGraphErrors* g_sigma_z = new TGraphErrors();
    g_mean_z->SetMarkerStyle(20); g_sigma_z->SetMarkerStyle(20);
    g_mean_z->SetMarkerSize(1.0); g_sigma_z->SetMarkerSize(1.0);
    g_mean_z->SetLineColor(kBlack); g_sigma_z->SetLineColor(kBlack);
    g_mean_z->SetMarkerColor(kBlack); g_sigma_z->SetMarkerColor(kBlack);
    
    for (size_t p = 0; p < particles.size(); ++p) {
        const string& particle = particles[p];
        if (p >= H5_CHARGE_BINS_EDGES.size() - 1) continue;
        double z_center = (H5_CHARGE_BINS_EDGES[p] + H5_CHARGE_BINS_EDGES[p+1]) / 2.0; 

        if (results.count(particle)) {
            const auto& res = results.at(particle);
            // 存储时只使用核心高斯的 mean 和 sigma
            if (res.mean_err > 0 || res.sigma_err > 0) {
                int n_mean = g_mean_z->GetN();
                g_mean_z->SetPoint(n_mean, z_center, res.mean);
                g_mean_z->SetPointError(n_mean, 0.0, res.mean_err);
                int n_sigma = g_sigma_z->GetN();
                g_sigma_z->SetPoint(n_sigma, z_center, res.sigma);
                g_sigma_z->SetPointError(n_sigma, 0.0, res.sigma_err);
            }
        }
    }
    
    double min_z = H5_CHARGE_BINS_EDGES.front();
    double max_z = H5_CHARGE_BINS_EDGES.back();
    g_mean_z->GetXaxis()->SetRangeUser(min_z, max_z);
    g_sigma_z->GetXaxis()->SetRangeUser(min_z, max_z);

    // H4图例标签改为 "ISS"
    vector<string> single_legend = {"ISS"};
    
    vector<TGraphErrors*> graphs_mean_z = {g_mean_z};
    DrawAndSaveGraphs(graphs_mean_z, single_legend, "Charge (Z)", ("#mu_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_name_base, false, 
                      save_to_root, "_Mean_vs_Charge");

    vector<TGraphErrors*> graphs_sigma_z = {g_sigma_z};
    DrawAndSaveGraphs(graphs_sigma_z, single_legend, "Charge (Z)", ("#sigma_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_name_base, false, 
                      save_to_root, "_Sigma_vs_Charge");

    delete g_mean_z; 
    delete g_sigma_z;
}

// --- 函数实现：绘制 H5 拟合结果 ---
void DrawH5ResultsFromTH2(TH2F* h2_mean, TH2F* h2_sigma, const string& output_dir, const HistInfo& info,
                          const string& output_name_base, const string& x_axis_label, TFile* save_to_root) {
    
    vector<int> colors = {kYellow+2, kMagenta, kBlack, kRed, kBlue, kGreen + 2, kOrange + 1};
    vector<int> rig_colors = {kBlack, kRed, kBlue, kGreen+2, kMagenta}; 

    vector<double> rig_centers;
    for (int i = 0; i < N_RIG_BINS; ++i) {
        rig_centers.push_back((H5_RIG_BINS_EDGES[i] + H5_RIG_BINS_EDGES[i+1]) / 2.0);
    }
    double rig_min = H5_RIG_BINS_EDGES.front();
    double rig_max = H5_RIG_BINS_EDGES.back() * 1.5; // 确保最右侧 bin 有显示空间

    vector<double> charge_centers;
    for (int i = 0; i < N_CHARGE_BINS; ++i) {
        charge_centers.push_back((H5_CHARGE_BINS_EDGES[i] + H5_CHARGE_BINS_EDGES[i+1]) / 2.0);
    }
    double charge_min = H5_CHARGE_BINS_EDGES.front();
    double charge_max = H5_CHARGE_BINS_EDGES.back();

    // 绘制 Mean/Sigma vs Rigidity (按 Z 分类)
    vector<TGraphErrors*> graphs_mean_rig, graphs_sigma_rig;
    for (int i = 1; i <= N_CHARGE_BINS; ++i) {
        TGraphErrors* g_mean = new TGraphErrors();
        TGraphErrors* g_sigma = new TGraphErrors();
        int color = colors[i-1];
        g_mean->SetMarkerStyle(20); g_sigma->SetMarkerStyle(20);
        g_mean->SetMarkerSize(1.0); g_sigma->SetMarkerSize(1.0);
        g_mean->SetMarkerColor(color); g_sigma->SetMarkerColor(color);
        g_mean->SetLineColor(color); g_sigma->SetLineColor(color);

        for (int j = 1; j <= N_RIG_BINS; ++j) {
            double mean = h2_mean->GetBinContent(j, i);
            double mean_err = h2_mean->GetBinError(j, i);
            double sigma = h2_sigma->GetBinContent(j, i);
            double sigma_err = h2_sigma->GetBinError(j, i);
            
            // 存储时只使用核心高斯的 mean 和 sigma
            if (mean_err > 0 || sigma_err > 0) {
                double rig_center = rig_centers[j-1];
                g_mean->SetPoint(g_mean->GetN(), rig_center, mean);
                g_mean->SetPointError(g_mean->GetN()-1, 0.0, mean_err); 
                g_sigma->SetPoint(g_sigma->GetN(), rig_center, sigma);
                g_sigma->SetPointError(g_sigma->GetN()-1, 0.0, sigma_err); 
            }
        }
        g_mean->GetXaxis()->SetRangeUser(rig_min, rig_max);
        g_sigma->GetXaxis()->SetRangeUser(rig_min, rig_max);
        graphs_mean_rig.push_back(g_mean);
        graphs_sigma_rig.push_back(g_sigma);
    }

    // mean vs Rigidity, logx=false
    DrawAndSaveGraphs(graphs_mean_rig, H5_CHARGE_LABELS, "Rigidity [GV]", ("#mu_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_name_base, false, 
                      save_to_root, "_Mean_vs_Rigidity");

    // sigma vs Rigidity, logx=false
    DrawAndSaveGraphs(graphs_sigma_rig, H5_CHARGE_LABELS, "Rigidity [GV]", ("#sigma_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_name_base, false, 
                      save_to_root, "_Sigma_vs_Rigidity");

    for (auto g : graphs_mean_rig) delete g;
    for (auto g : graphs_sigma_rig) delete g;

    // 绘制 Mean/Sigma vs Charge (按 Rigidity 分类)
    vector<TGraphErrors*> graphs_mean_z, graphs_sigma_z;
    for (int j = 1; j <= N_RIG_BINS; ++j) {
        TGraphErrors* g_mean_z = new TGraphErrors();
        TGraphErrors* g_sigma_z = new TGraphErrors();
        int color = rig_colors[j-1];
        g_mean_z->SetMarkerStyle(20); g_sigma_z->SetMarkerStyle(20);
        g_mean_z->SetMarkerSize(1.0); g_sigma_z->SetMarkerSize(1.0);
        g_mean_z->SetMarkerColor(color); g_sigma_z->SetMarkerColor(color);
        g_mean_z->SetLineColor(color); g_sigma_z->SetLineColor(color);

        for (int i = 1; i <= N_CHARGE_BINS; ++i) {
            double mean = h2_mean->GetBinContent(j, i);
            double mean_err = h2_mean->GetBinError(j, i);
            double sigma = h2_sigma->GetBinContent(j, i);
            double sigma_err = h2_sigma->GetBinError(j, i);

            // 存储时只使用核心高斯的 mean 和 sigma
            if (mean_err > 0 || sigma_err > 0) {
                double z_center = charge_centers[i-1];
                g_mean_z->SetPoint(g_mean_z->GetN(), z_center, mean);
                g_mean_z->SetPointError(g_mean_z->GetN()-1, 0.0, mean_err);
                g_sigma_z->SetPoint(g_sigma_z->GetN(), z_center, sigma);
                g_sigma_z->SetPointError(g_sigma_z->GetN()-1, 0.0, sigma_err);
            }
        }
        g_mean_z->GetXaxis()->SetRangeUser(charge_min, charge_max);
        g_sigma_z->GetXaxis()->SetRangeUser(charge_min, charge_max);
        graphs_mean_z.push_back(g_mean_z);
        graphs_sigma_z.push_back(g_sigma_z);
    }

    // mean vs Charge, logx=false
    DrawAndSaveGraphs(graphs_mean_z, H5_RIG_LABELS, "Charge (Z)", ("#mu_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Mean", output_dir, output_name_base, false, 
                      save_to_root, "_Mean_vs_Charge");

    // sigma vs Charge, logx=false
    DrawAndSaveGraphs(graphs_sigma_z, H5_RIG_LABELS, "Charge (Z)", ("#sigma_{" + x_axis_label + "}").c_str(), 
                      info.title_description + " Sigma", output_dir, output_name_base, false, 
                      save_to_root, "_Sigma_vs_Charge");

    for (auto g : graphs_mean_z) delete g;
    for (auto g : graphs_sigma_z) delete g;
}

// --- 函数实现：两轮拟合 (更新：存储和显示 LR/RR 误差) ---
FitResult twoStepGaussianFit(TH1* hist, const string& x_axis_label, TCanvas* c_fit, const string& title_prefix, double default_center_x, const string& rig_label) {
    // 初始化结果结构体
    FitResult result = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (!hist || hist->GetEntries() < 10) return result;

    // --- 第一轮拟合：使用 80% 积分范围确定初始参数 (普通高斯) ---
    double x_min_fit1 = 0.0, x_max_fit1 = 0.0;
    findFitRange(hist, 0.80, default_center_x, x_min_fit1, x_max_fit1); // 80% 覆盖率
    if (x_max_fit1 <= x_min_fit1) return result;

    TF1* f_gaus1 = new TF1("f_gaus1", "gaus", x_min_fit1, x_max_fit1);
    f_gaus1->SetParameters(hist->GetMaximum(), hist->GetMean(), hist->GetRMS());
    int fit_status1 = hist->Fit(f_gaus1, "QRS");
    
    if (fit_status1 != 0) {
        delete f_gaus1;
        return result; 
    }

    double mean0 = f_gaus1->GetParameter(1);
    double sigma0 = f_gaus1->GetParameter(2);
    delete f_gaus1; // 释放内存

    // --- 第二轮拟合：使用 DoGausPlusAsymGausFit (核心高斯 +/- 3*sigma 范围) ---
    // 使用第一轮的高斯结果确定拟合范围
    double x_min_fit2 = mean0 - 4.0 * abs(sigma0);
    double x_max_fit2 = mean0 + 4.0 * abs(sigma0);
    //if Beryllium and 30-50GV, using narrower range:2.0
    if(rig_label == "30-50 GV" && title_prefix.find("Beryllium") != string::npos){
        x_min_fit2 = mean0 - 2.5 * abs(sigma0);
        x_max_fit2 = mean0 + 2.5 * abs(sigma0);
    }

    if (x_max_fit2 <= x_min_fit2) return result;
    
    // 调用封装的拟合函数
    // Row 0: Values [Mean, Sigma_core, LR, RR, Chi2]
    // Row 1: Errors & NDF [Mean_err, Sigma_core_err, LR_err, RR_err, NDF]
    hist->SetTitle(title_prefix.c_str());
    hist->GetYaxis()->SetTitle("Events");
    //hist->GetXaxis()->SetTitle(x_axis_label.c_str());
    double buffer = 0.6 * (x_max_fit2 - x_min_fit2);
    hist->GetXaxis()->SetRangeUser(x_min_fit2 - buffer, x_max_fit2 + buffer);
    vector<vector<double>> fit_data = DoGausPlusAsymGausFit(hist, x_min_fit2, x_max_fit2, c_fit, true);
    
    if (fit_data[1][4] > 0) { // NDF > 0 表示拟合成功
        // 成功拟合，将结果填充到 FitResult 结构体
        result.mean = fit_data[0][0];
        result.mean_err = fit_data[1][0];
        result.sigma = fit_data[0][1];
        result.sigma_err = fit_data[1][1];
        result.LR = fit_data[0][2];
        result.LR_err = fit_data[1][2]; // 存储 LR 误差
        result.RR = fit_data[0][3];
        result.RR_err = fit_data[1][3]; // 存储 RR 误差
        result.chi2 = fit_data[0][4];
        result.ndf = fit_data[1][4];
        c_fit->cd();
        
        // 标记拟合范围
        TLine* l_min = new TLine(x_min_fit2, 0, x_min_fit2, hist->GetMaximum() * 1.);
        TLine* l_max = new TLine(x_max_fit2, 0, x_max_fit2, hist->GetMaximum() * 1.);
        l_min->SetLineStyle(2); l_max->SetLineStyle(2);
        l_min->SetLineColor(kRed); l_max->SetLineColor(kRed);
        l_min->Draw("same"); l_max->Draw("same");

        c_fit->Update();
    } else {
        // 拟合失败
        c_fit->cd();
        hist->Draw("hist");
        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
        latex.DrawLatex(0.6, 0.85, "Fit Failed");
        if (!rig_label.empty()) {
            stringstream ss; ss << "Rig: " << rig_label;
            latex.DrawLatex(0.15, 0.85, ss.str().c_str());
        }
        c_fit->Update();
    }
    
    return result;
}

// --- 函数实现：分析主逻辑 ---
void Analyze(TFile* file, const string& suffix, const string& output_dir,
             const vector<string>& particles, const vector<int>& colors, TFile* output_root_file, bool is_h5) {

    cout << "\n--- Analyzing " << (is_h5 ? "H5" : "H4") << " Series: " << suffix << " ---" << endl;
    gStyle->SetOptFit(0);
    
    const string prefix = "UnbiasedL1Inner_";
    const string x_axis_label = is_h5 ? "#Delta(1/#beta)" : "1/#beta";
    HistInfo info = getHistInfo(suffix);
    string output_name_base = info.output_suffix + suffix.substr(3);
    string fit_pdf_path = output_dir + "FitResults_" + output_name_base + ".pdf";
    TCanvas* c_fit = new TCanvas("c_fit", "Gaussian Fit Results", 800, 600);
    c_fit->SetLogy(0);
    c_fit->Print((fit_pdf_path + "[").c_str(), "pdf");

    H4ResultsMap h4_results;
    TH2F* h2_mean = nullptr;
    TH2F* h2_sigma = nullptr;
    if (is_h5) {
        h2_mean = new TH2F(("h2_mean_" + suffix).c_str(), (info.title_description + " Mean;Rigidity [GV];Charge (Z)").c_str(), N_RIG_BINS, &H5_RIG_BINS_EDGES[0], N_CHARGE_BINS, &H5_CHARGE_BINS_EDGES[0]);
        h2_sigma = new TH2F(("h2_sigma_" + suffix).c_str(), (info.title_description + " Sigma;Rigidity [GV];Charge (Z)").c_str(), N_RIG_BINS, &H5_RIG_BINS_EDGES[0], N_CHARGE_BINS, &H5_CHARGE_BINS_EDGES[0]);
    }

    for (size_t p = 0; p < particles.size(); ++p) {
        const string& particle = particles[p];
        string full_name = prefix + particle + "_" + suffix;
        
        TObject* obj = file->Get(full_name.c_str());
        if (!obj) continue;

        if (is_h5) {
            TH2F* h2 = dynamic_cast<TH2F*>(obj);
            if (!h2) continue;

            for (size_t i = 0; i < N_RIG_BINS; ++i) {
                double rig_min = H5_RIG_BINS_EDGES[i];
                double rig_max = H5_RIG_BINS_EDGES[i+1];
                int rig_bin_index = i + 1;
                string rig_label = H5_RIG_LABELS[i];

                // 2D 上的 Y 轴是 Rigidity，X 轴是要拟合的量 (Delta 1/beta)
                int bin_y_min = h2->GetYaxis()->FindFixBin(rig_min);
                // 确保最后一个 bin 包含最后一个边界
                int bin_y_max = (i == N_RIG_BINS - 1) ? h2->GetNbinsY() : h2->GetYaxis()->FindFixBin(rig_max - 1e-6);

                TH1D* h1_proj = h2->ProjectionX(("proj_" + full_name + "_" + to_string(i)).c_str(), bin_y_min, bin_y_max);
                if (h1_proj->GetEntries() == 0) { delete h1_proj; continue; }
                
                // 尝试 Rebin 直到最大值大于60
                for(int r=1 ; r<=10; ++r) { if (h1_proj->GetMaximum() >= 80 && h1_proj->GetBinWidth(1)>0.000001) break; h1_proj->Rebin(2); }
                
                string title = particle + " - " + info.output_suffix + " [" + rig_label + "]";

                FitResult fit_res = twoStepGaussianFit(h1_proj, x_axis_label, c_fit, title, 0.0, rig_label);

                // 存储的时候只用存核心高斯的 mean 和 sigma
                if (fit_res.mean_err > 0 || fit_res.sigma_err > 0) {
                    h2_mean->SetBinContent(rig_bin_index, p + 1, fit_res.mean);
                    h2_mean->SetBinError(rig_bin_index, p + 1, fit_res.mean_err);
                    h2_sigma->SetBinContent(rig_bin_index, p + 1, fit_res.sigma);
                    h2_sigma->SetBinError(rig_bin_index, p + 1, fit_res.sigma_err);
                    c_fit->Print(fit_pdf_path.c_str(), "pdf");
                }
                
                delete h1_proj;
            }
        } else {
            TH1F* h1 = dynamic_cast<TH1F*>(obj);
            if (!h1) continue;
            TH1F* h1_clone = (TH1F*)h1->Clone(("clone_" + full_name).c_str());
            h1_clone->SetMarkerStyle(20);
            h1_clone->SetMarkerSize(1.2);

            // 尝试 Rebin 直到最大值大于60
            for(int r=1 ; r<=10; ++r) { if (h1_clone->GetMaximum() >= 80 && h1_clone->GetBinWidth(1)>0.000001) break; h1_clone->Rebin(2); }
            
            string title = particle + " - " + info.title_description;
            // 默认中心点 1.0 用于 H4 拟合
            FitResult fit_res = twoStepGaussianFit(h1_clone, x_axis_label, c_fit, title, 1.0, ""); 

            if (fit_res.mean_err > 0 || fit_res.sigma_err > 0) {
                h4_results[particle] = fit_res;
                c_fit->Print(fit_pdf_path.c_str(), "pdf");
            }

            delete h1_clone;
        }
    }
    
    c_fit->Print((fit_pdf_path + "]").c_str(), "pdf");
    delete c_fit;

    if (is_h5) {
        output_root_file->cd();
        h2_mean->Write();
        h2_sigma->Write();
        DrawH5ResultsFromTH2(h2_mean, h2_sigma, output_dir, info, output_name_base, x_axis_label, output_root_file);
        delete h2_mean;
        delete h2_sigma;
    } else {
        // H4 Results with "ISS" legend
        DrawH4Results(h4_results, output_dir, info, output_name_base, x_axis_label, particles, output_root_file);
    }
}


// --- ROOT 宏主函数 ---
void BetaFit() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetErrorX(0); 

    const string input_file_path = "/eos/user/z/zixuan/Isotope/Add/Be_frag4.root";
    const string output_dir = "/eos/user/z/zixuan/Isotope/Beta/ISS/";
    const string output_root_path = output_dir + "ISS_RICHBetaStudy.root";
    
    vector<string> particles = {"Helium","Lithium","Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
    vector<int> colors = {kYellow+2, kMagenta, kBlack, kRed, kBlue, kGreen + 2, kOrange + 1}; 

    if (gSystem->AccessPathName(output_dir.c_str())) {
        cout << "Creating output directory: " << output_dir << endl;
        gSystem->mkdir(output_dir.c_str(), kTRUE);
    }

    TFile* file = TFile::Open(input_file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "ERROR: Cannot open input file: " << input_file_path << endl;
        return;
    }

    TFile* output_root_file = TFile::Open(output_root_path.c_str(), "RECREATE");
    if (!output_root_file || output_root_file->IsZombie()) {
        cerr << "ERROR: Cannot create output ROOT file: " << output_root_path << endl;
        file->Close();
        delete file;
        return;
    }
    
    // H5 系列分析 (is_h5 = true)
    vector<string> h5_suffixes = {"ID_H5a", "ID_H5b"};
    for (const string& suffix : h5_suffixes) {
        Analyze(file, suffix, output_dir, particles, colors, output_root_file, true);
    }

    // H4 系列分析 (is_h5 = false)
    vector<string> h4_suffixes = {"ID_H4a", "ID_H4b"};
    for (const string& suffix : h4_suffixes) {
        Analyze(file, suffix, output_dir, particles, colors, output_root_file, false);
    }

    output_root_file->Close();
    file->Close();
    delete output_root_file;
    delete file;

    gROOT->SetBatch(kFALSE);
    cout << "\nAnalysis finished. Results saved to: " << output_dir << endl;
}