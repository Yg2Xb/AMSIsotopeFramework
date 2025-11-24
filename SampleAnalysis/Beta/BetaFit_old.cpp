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
#include <TRegexp.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iomanip>
#include <sstream>
#include <TStyle.h>
#include <TLine.h>

using namespace std;

// --- 辅助结构和函数定义 ---

struct HistInfo {
    string y_axis_label;
    string output_suffix;
    string title_description;
};

struct FitResult {
    double mean;
    double mean_err;
    double sigma;
    double sigma_err;
};

// 用 map 存储 H4 拟合结果： key=particle name
using H4ResultsMap = map<string, FitResult>; 

// 仅用于 H5 系列的固定 Rigidity Bin (X轴)
const vector<double> H5_RIG_BINS_EDGES = {30.0, 50.0, 80.0, 120.0, 160.0, 240.0};
const vector<string> H5_RIG_LABELS = {"30-50 GV", "50-80 GV", "80-120 GV", "120-160 GV", ">160 GV"};
const int N_RIG_BINS = H5_RIG_BINS_EDGES.size() - 1;

// 用于 H5 系列的固定 Charge/Particle Bin (Y轴)
const vector<string> H5_CHARGE_LABELS = {"Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"};
const vector<double> H5_CHARGE_BINS_EDGES = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5}; // Z=2 到 Z=8
const int N_CHARGE_BINS = H5_CHARGE_BINS_EDGES.size() - 1;


// --- 修改后的 findFitRange 函数：基于指定中心点对称扩展 ---
void findFitRange(TH1* hist, double coverage, double center_x, double& x_min, double& x_max) {
    if (!hist || hist->GetEntries() == 0) return;

    double total_integral = hist->Integral();
    if (total_integral <= 0) return;
    
    // 1. 找到初始中心点所在的 bin
    int center_bin = hist->GetXaxis()->FindFixBin(center_x);
    if (center_bin < 1 || center_bin > hist->GetNbinsX()) {
        // 如果指定中心点不在范围内，回退到最大值 bin (作为起点)
        center_bin = hist->GetMaximumBin(); 
    }

    double required_integral = total_integral * coverage;
    double current_integral = hist->GetBinContent(center_bin);
    
    int low_bin = center_bin - 1;
    int high_bin = center_bin + 1;
    
    x_min = hist->GetXaxis()->GetBinLowEdge(center_bin);
    x_max = hist->GetXaxis()->GetBinUpEdge(center_bin);

    // 2. 对称地向外扩展，直到达到所需覆盖率
    while (current_integral < required_integral) {
        bool extended = false;
        double content_low = (low_bin >= 1) ? hist->GetBinContent(low_bin) : 0;
        double content_high = (high_bin <= hist->GetNbinsX()) ? hist->GetBinContent(high_bin) : 0;
        
        // 尝试向两侧扩展
        if (low_bin >= 1 && high_bin <= hist->GetNbinsX()) {
            current_integral += content_low + content_high;
            x_min = hist->GetXaxis()->GetBinLowEdge(low_bin);
            x_max = hist->GetXaxis()->GetBinUpEdge(high_bin);
            low_bin--;
            high_bin++;
            extended = true;
        } 
        // 如果只剩下一边可以扩展
        else if (low_bin >= 1) {
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
    
    // 3. 确保拟合范围相对于指定的中心点 (center_x) 对称
    double final_center = center_x;
    double half_range = max(abs(x_max - final_center), abs(x_min - final_center));
    
    x_min = final_center - half_range;
    x_max = final_center + half_range;
}

// 计算 TGraph 集合在特定 X 范围内的全局 Y 范围
void getGlobalYRange(const vector<TGraphErrors*>& graphs, double x_min, double x_max, double& y_min, double& y_max) {
    y_min = 1e10;
    y_max = -1e10;
    bool found_point = false;

    for (const auto& g : graphs) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);

            if (x >= x_min && x <= x_max) {
                double y_low = y;
                double y_high = y;
                
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
        
        double buffer = range * 0.10; // 10% 扩充
        y_min -= buffer;
        y_max += buffer;
    } else {
        y_min = -0.01;
        y_max = 0.01;
    }
}

HistInfo getHistInfo(const string& suffix) {
    HistInfo info;
    
    if (suffix == "ID_H5a") { info = {"Rigidity [GV]", "Rigidity", "NaF-Tracker #Delta(1/#beta)"}; }
    else if (suffix == "ID_H5b") { info = {"Rigidity [GV]", "Rigidity", "AGL-Tracker #Delta(1/#beta)"}; }
    else if (suffix == "ID_H4a") { info = {"Rigidity [GV]", "Rigidity", "NaF 1/#beta (Tracker Rig)"}; }
    else if (suffix == "ID_H4b") { info = {"Rigidity [GV]", "Rigidity", "AGL 1/#beta (Tracker Rig)"}; }
    else { info = {"Y Variable", "YVariable", "Unknown Delta Beta"}; }
    
    return info;
}


// --- 新增 H4 系列绘图函数：Mean/Sigma vs Charge (Z) ---
void DrawH4Results(const H4ResultsMap& results, const string& output_dir, const HistInfo& info, 
                   const string& output_name_base, const string& x_axis_label, 
                   const vector<string>& particles, const vector<int>& colors) {
    
    // 创建 TGraphErrors  
    TGraphErrors* g_mean_z = new TGraphErrors();
    TGraphErrors* g_sigma_z = new TGraphErrors();

    g_mean_z->SetMarkerStyle(20); g_sigma_z->SetMarkerStyle(20);
    g_mean_z->SetMarkerSize(1.0); g_sigma_z->SetMarkerSize(1.0);
    g_mean_z->SetLineColor(kBlack); g_sigma_z->SetLineColor(kBlack);
    g_mean_z->SetMarkerColor(kBlack); g_sigma_z->SetMarkerColor(kBlack);
    
    vector<TGraphErrors*> graphs_mean_z_vec = {g_mean_z}; // 用于 getGlobalYRange
    vector<TGraphErrors*> graphs_sigma_z_vec = {g_sigma_z};

    double min_z = H5_CHARGE_BINS_EDGES.front();
    double max_z = H5_CHARGE_BINS_EDGES.back();

    for (size_t p = 0; p < particles.size(); ++p) {
        const string& particle = particles[p];
        // 假设 H5_CHARGE_BINS_EDGES 的定义是从 Z=2 开始的
        double z_center = (H5_CHARGE_BINS_EDGES[p] + H5_CHARGE_BINS_EDGES[p+1]) / 2.0; 

        if (results.count(particle)) {
            const auto& res = results.at(particle);
            if (res.mean_err > 0 || res.sigma_err > 0) { // 仅绘制有效点
                g_mean_z->SetPoint(g_mean_z->GetN(), z_center, res.mean);
                g_mean_z->SetPointError(g_mean_z->GetN()-1, 0.0, res.mean_err); // X Error = 0
                g_sigma_z->SetPoint(g_sigma_z->GetN(), z_center, res.sigma);
                g_sigma_z->SetPointError(g_sigma_z->GetN()-1, 0.0, res.sigma_err); // X Error = 0
            }
        }
    }

    // --- 绘制 Mean vs Charge ---
    TCanvas* c_mean_charge = new TCanvas(("c_H4_mean_charge_" + output_name_base).c_str(), "H4 Mean vs Charge", 700, 500);
    c_mean_charge->SetGrid();
    double mean_z_min, mean_z_max;
    getGlobalYRange(graphs_mean_z_vec, min_z, max_z, mean_z_min, mean_z_max);

    g_mean_z->Draw("APZ");
    g_mean_z->SetTitle((info.title_description + " Mean vs Charge").c_str());
    g_mean_z->GetXaxis()->SetTitle("Charge (Z)");
    g_mean_z->GetYaxis()->SetTitle(("#mu_{" + x_axis_label + "}").c_str());
    g_mean_z->GetXaxis()->SetRangeUser(min_z, max_z);
    g_mean_z->GetYaxis()->SetRangeUser(mean_z_min, mean_z_max);
    
    c_mean_charge->SaveAs((output_dir + output_name_base + "_Mean_vs_Charge.png").c_str());
    
    // --- 绘制 Sigma vs Charge ---
    TCanvas* c_sigma_charge = new TCanvas(("c_H4_sigma_charge_" + output_name_base).c_str(), "H4 Sigma vs Charge", 700, 500);
    c_sigma_charge->SetGrid();
    double sigma_z_min, sigma_z_max;
    getGlobalYRange(graphs_sigma_z_vec, min_z, max_z, sigma_z_min, sigma_z_max);

    g_sigma_z->Draw("APZ");
    g_sigma_z->SetTitle((info.title_description + " Sigma vs Charge").c_str());
    g_sigma_z->GetXaxis()->SetTitle("Charge (Z)");
    g_sigma_z->GetYaxis()->SetTitle(("#sigma_{" + x_axis_label + "}").c_str());
    g_sigma_z->GetXaxis()->SetRangeUser(min_z, max_z);
    g_sigma_z->GetYaxis()->SetRangeUser(sigma_z_min, sigma_z_max);

    c_sigma_charge->SaveAs((output_dir + output_name_base + "_Sigma_vs_Charge.png").c_str());

    // 清理
    delete g_mean_z; delete g_sigma_z;
    delete c_mean_charge; delete c_sigma_charge;
}

// --- H5 系列绘图函数 (保持不变) ---
void DrawH5ResultsFromTH2(TH2F* h2_mean, TH2F* h2_sigma, const string& output_dir, const HistInfo& info, const string& suffix, const string& output_name_base, const string& x_axis_label) {

    vector<int> colors = {kYellow+2, kMagenta, kBlack, kRed, kBlue, kGreen + 2, kOrange + 1};
    vector<int> rig_colors = {kBlack, kRed, kBlue, kGreen+2, kMagenta}; // Rigidity 边框颜色

    // 1. Mean/Sigma vs Rigidity (X轴) - 7 条曲线
    vector<TGraphErrors*> graphs_mean;
    vector<TGraphErrors*> graphs_sigma;

    // 计算 Rigidity Bin 中心值
    vector<double> rig_centers;
    for (int i = 0; i < N_RIG_BINS; ++i) {
        rig_centers.push_back((H5_RIG_BINS_EDGES[i] + H5_RIG_BINS_EDGES[i+1]) / 2.0);
    }

    // 循环粒子 (Y轴 Bin) 准备 TGraph
    for (int i = 1; i <= N_CHARGE_BINS; ++i) {
        TGraphErrors* g_mean = new TGraphErrors();
        TGraphErrors* g_sigma = new TGraphErrors();
        int color = colors[i-1];
        
        // 统一 Marker Style 和 Size
        g_mean->SetMarkerStyle(20); g_sigma->SetMarkerStyle(20);
        g_mean->SetMarkerSize(1.0); g_sigma->SetMarkerSize(1.0);
        g_mean->SetMarkerColor(color); g_sigma->SetMarkerColor(color);
        g_mean->SetLineColor(color); g_sigma->SetLineColor(color);

        for (int j = 1; j <= N_RIG_BINS; ++j) { // 循环 Rigidity Bin (X轴 Bin)
            double mean = h2_mean->GetBinContent(j, i);
            double mean_err = h2_mean->GetBinError(j, i);
            double sigma = h2_sigma->GetBinContent(j, i);
            double sigma_err = h2_sigma->GetBinError(j, i);
            
            if (mean_err > 0 || sigma_err > 0) { // 仅绘制有效点
                double rig_center = rig_centers[j-1];
                // X轴误差为0，但为了保持 TGraphErrors 的数据结构，我们只将 XError 设为 0
                
                g_mean->SetPoint(g_mean->GetN(), rig_center, mean);
                g_mean->SetPointError(g_mean->GetN()-1, 0.0, mean_err); // X Error = 0
                g_sigma->SetPoint(g_sigma->GetN(), rig_center, sigma);
                g_sigma->SetPointError(g_sigma->GetN()-1, 0.0, sigma_err); // X Error = 0
            }
        }
        graphs_mean.push_back(g_mean);
        graphs_sigma.push_back(g_sigma);
    }

    // --- 1.1 绘制 Mean vs Rigidity ---
    TCanvas* c_mean_rig = new TCanvas(("c_mean_rig_" + output_name_base).c_str(), "Mean vs Rigidity", 700, 500);
    c_mean_rig->SetLogx(0); c_mean_rig->SetGrid();
    double mean_y_min, mean_y_max;
    getGlobalYRange(graphs_mean, H5_RIG_BINS_EDGES.front(), H5_RIG_BINS_EDGES.back() + 100.0, mean_y_min, mean_y_max);
    
    TLine* line_h5_ref = new TLine(H5_RIG_BINS_EDGES.front(), 0.0, H5_RIG_BINS_EDGES.back() + 100.0, 0.0);
    line_h5_ref->SetLineColor(kBlue); line_h5_ref->SetLineStyle(2);

    TLegend* leg_mean = new TLegend(0.8, 0.8, 0.99, 0.99); leg_mean->SetFillStyle(0); leg_mean->SetBorderSize(1);
    for (size_t i = 0; i < graphs_mean.size(); ++i) {
        string draw_opt = (i == 0) ? "APZ" : "PZ same";
        graphs_mean[i]->Draw(draw_opt.c_str());
        leg_mean->AddEntry(graphs_mean[i], H5_CHARGE_LABELS[i].c_str(), "p");

        if (i == 0) {
            graphs_mean[i]->SetTitle((info.title_description + " Mean vs Rigidity").c_str());
            graphs_mean[i]->GetXaxis()->SetTitle("Rigidity [GV]");
            graphs_mean[i]->GetYaxis()->SetTitle(("#mu_{" + x_axis_label + "}").c_str());
            graphs_mean[i]->GetXaxis()->SetRangeUser(H5_RIG_BINS_EDGES.front(), H5_RIG_BINS_EDGES.back() + 100.0);
            graphs_mean[i]->GetYaxis()->SetRangeUser(mean_y_min, mean_y_max);
        }
    }
    leg_mean->Draw();
    c_mean_rig->SaveAs((output_dir + output_name_base + "_Mean_vs_Rigidity.png").c_str());
    delete line_h5_ref;

    // --- 1.2 绘制 Sigma vs Rigidity ---
    TCanvas* c_sigma_rig = new TCanvas(("c_sigma_rig_" + output_name_base).c_str(), "Sigma vs Rigidity", 700, 500);
    c_sigma_rig->SetLogx(1); c_sigma_rig->SetGrid();
    double sigma_y_min, sigma_y_max;
    getGlobalYRange(graphs_sigma, H5_RIG_BINS_EDGES.front(), H5_RIG_BINS_EDGES.back() + 100.0, sigma_y_min, sigma_y_max);
    
    TLegend* leg_sigma = new TLegend(0.8, 0.8, 0.99, 0.99); leg_sigma->SetFillStyle(0); leg_sigma->SetBorderSize(1);
    for (size_t i = 0; i < graphs_sigma.size(); ++i) {
        string draw_opt = (i == 0) ? "APZ" : "PZ same";
        graphs_sigma[i]->Draw(draw_opt.c_str());
        leg_sigma->AddEntry(graphs_sigma[i], H5_CHARGE_LABELS[i].c_str(), "p");

        if (i == 0) {
            graphs_sigma[i]->SetTitle((info.title_description + " Sigma vs Rigidity").c_str());
            graphs_sigma[i]->GetXaxis()->SetTitle("Rigidity [GV]");
            graphs_sigma[i]->GetYaxis()->SetTitle(("#sigma_{" + x_axis_label + "}").c_str());
            graphs_sigma[i]->GetXaxis()->SetRangeUser(H5_RIG_BINS_EDGES.front(), H5_RIG_BINS_EDGES.back() + 100.0);
            graphs_sigma[i]->GetYaxis()->SetRangeUser(sigma_y_min, sigma_y_max);
        }
    }
    leg_sigma->Draw();
    c_sigma_rig->SaveAs((output_dir + output_name_base + "_Sigma_vs_Rigidity.png").c_str());

    // 清理 Rigidity plots
    for (auto g : graphs_mean) delete g;
    for (auto g : graphs_sigma) delete g;
    delete leg_mean; delete leg_sigma;
    delete c_mean_rig; delete c_sigma_rig;

    // --- 2. Mean/Sigma vs Charge (Y轴) - 5 (或 N_RIG_BINS) 条曲线 (保持不变) ---
    vector<TGraphErrors*> graphs_mean_z;
    vector<TGraphErrors*> graphs_sigma_z;

    // 循环 Rigidity Bin (X轴 Bin) 准备 TGraph
    for (int j = 1; j <= N_RIG_BINS; ++j) {
        TGraphErrors* g_mean_z = new TGraphErrors();
        TGraphErrors* g_sigma_z = new TGraphErrors();
        int color = rig_colors[j-1];
        
        // 统一 Marker Style 和 Size
        g_mean_z->SetMarkerStyle(20); g_sigma_z->SetMarkerStyle(20);
        g_mean_z->SetMarkerSize(1.0); g_sigma_z->SetMarkerSize(1.0);
        g_mean_z->SetMarkerColor(color); g_sigma_z->SetMarkerColor(color);
        g_mean_z->SetLineColor(color); g_sigma_z->SetLineColor(color);

        for (int i = 1; i <= N_CHARGE_BINS; ++i) { // 循环 Charge Bin (Y轴 Bin)
            double mean = h2_mean->GetBinContent(j, i);
            double mean_err = h2_mean->GetBinError(j, i);
            double sigma = h2_sigma->GetBinContent(j, i);
            double sigma_err = h2_sigma->GetBinError(j, i);

            if (mean_err > 0 || sigma_err > 0) { // 仅绘制有效点
                double z_center = (H5_CHARGE_BINS_EDGES[i-1] + H5_CHARGE_BINS_EDGES[i])/2.0;
                
                g_mean_z->SetPoint(g_mean_z->GetN(), z_center, mean);
                g_mean_z->SetPointError(g_mean_z->GetN()-1, 0.0, mean_err); // X Error = 0
                g_sigma_z->SetPoint(g_sigma_z->GetN(), z_center, sigma);
                g_sigma_z->SetPointError(g_sigma_z->GetN()-1, 0.0, sigma_err); // X Error = 0
            }
        }
        graphs_mean_z.push_back(g_mean_z);
        graphs_sigma_z.push_back(g_sigma_z);
    }

    // --- 2.1 绘制 Mean vs Charge ---
    TCanvas* c_mean_charge = new TCanvas(("c_mean_charge_" + output_name_base).c_str(), "Mean vs Charge", 700, 500);
    c_mean_charge->SetGrid();
    double mean_z_min, mean_z_max;
    getGlobalYRange(graphs_mean_z, H5_CHARGE_BINS_EDGES.front(), H5_CHARGE_BINS_EDGES.back(), mean_z_min, mean_z_max);

    TLegend* leg_mean_z = new TLegend(0.8, 0.8, 0.99, 0.99); leg_mean_z->SetFillStyle(0); leg_mean_z->SetBorderSize(1);
    for (size_t i = 0; i < graphs_mean_z.size(); ++i) {
        string draw_opt = (i == 0) ? "APZ" : "PZ same";
        graphs_mean_z[i]->Draw(draw_opt.c_str());
        leg_mean_z->AddEntry(graphs_mean_z[i], H5_RIG_LABELS[i].c_str(), "p");

        if (i == 0) {
            graphs_mean_z[i]->SetTitle((info.title_description + " Mean vs Charge").c_str());
            graphs_mean_z[i]->GetXaxis()->SetTitle("Charge (Z)");
            graphs_mean_z[i]->GetYaxis()->SetTitle(("#mu_{" + x_axis_label + "}").c_str());
            graphs_mean_z[i]->GetXaxis()->SetRangeUser(H5_CHARGE_BINS_EDGES.front(), H5_CHARGE_BINS_EDGES.back());
            graphs_mean_z[i]->GetYaxis()->SetRangeUser(mean_z_min, mean_z_max);
        }
    }
    leg_mean_z->Draw();
    c_mean_charge->SaveAs((output_dir + output_name_base + "_Mean_vs_Charge.png").c_str());

    // --- 2.2 绘制 Sigma vs Charge ---
    TCanvas* c_sigma_charge = new TCanvas(("c_sigma_charge_" + output_name_base).c_str(), "Sigma vs Charge", 700, 500);
    c_sigma_charge->SetGrid();
    double sigma_z_min, sigma_z_max;
    getGlobalYRange(graphs_sigma_z, H5_CHARGE_BINS_EDGES.front(), H5_CHARGE_BINS_EDGES.back(), sigma_z_min, sigma_z_max);

    TLegend* leg_sigma_z = new TLegend(0.8, 0.8, 0.99, 0.99); leg_sigma_z->SetFillStyle(0); leg_sigma_z->SetBorderSize(1);
    for (size_t i = 0; i < graphs_sigma_z.size(); ++i) {
        string draw_opt = (i == 0) ? "APZ" : "PZ same";
        graphs_sigma_z[i]->Draw(draw_opt.c_str());
        leg_sigma_z->AddEntry(graphs_sigma_z[i], H5_RIG_LABELS[i].c_str(), "p");

        if (i == 0) {
            graphs_sigma_z[i]->SetTitle((info.title_description + " Sigma vs Charge").c_str());
            graphs_sigma_z[i]->GetXaxis()->SetTitle("Charge (Z)");
            graphs_sigma_z[i]->GetYaxis()->SetTitle(("#sigma_{" + x_axis_label + "}").c_str());
            graphs_sigma_z[i]->GetXaxis()->SetRangeUser(H5_CHARGE_BINS_EDGES.front(), H5_CHARGE_BINS_EDGES.back());
            graphs_sigma_z[i]->GetYaxis()->SetRangeUser(sigma_z_min, sigma_z_max);
        }
    }
    leg_sigma_z->Draw();
    c_sigma_charge->SaveAs((output_dir + output_name_base + "_Sigma_vs_Charge.png").c_str());

    // 清理 Charge plots
    for (auto g : graphs_mean_z) delete g;
    for (auto g : graphs_sigma_z) delete g;
    delete leg_mean_z; delete leg_sigma_z;
    delete c_mean_charge; delete c_sigma_charge;
}

// --- 核心分析函数：H5 系列 (TH2F 投影) ---
void AnalyzeH5(TFile* file, const string& suffix, const string& output_dir,
    const vector<string>& particles) {

    cout << "\n--- Analyzing H5 Series: " << suffix << " ---" << endl;
    gStyle->SetOptFit(0);
    
    const string prefix = "L1Inner_";
    const string x_axis_label = "#Delta(1/#beta)";
    HistInfo info = getHistInfo(suffix);
    string output_name_base = info.output_suffix + suffix.substr(3);

    // 拟合结果 PDF 文件
    string fit_pdf_path = output_dir + "FitResults_" + output_name_base + ".pdf";
    TCanvas* c_fit = new TCanvas("c_fit_H5", "Gaussian Fit Results (H5)", 800, 600);
    c_fit->SetLogy(0);
    c_fit->Print((fit_pdf_path + "[").c_str(), "pdf");

    // 创建最终的 TH2F 结果存储
    // X轴: Rigidity (5 bins), Y轴: Charge (7 bins)
    TH2F* h2_mean = new TH2F(("h2_mean_" + suffix).c_str(), 
        (info.title_description + " Mean;Rigidity [GV];Charge (Z)").c_str(), 
        N_RIG_BINS, &H5_RIG_BINS_EDGES[0], N_CHARGE_BINS, &H5_CHARGE_BINS_EDGES[0]);
    TH2F* h2_sigma = new TH2F(("h2_sigma_" + suffix).c_str(), 
        (info.title_description + " Sigma;Rigidity [GV];Charge (Z)").c_str(), 
        N_RIG_BINS, &H5_RIG_BINS_EDGES[0], N_CHARGE_BINS, &H5_CHARGE_BINS_EDGES[0]);

    // --- 循环粒子 ---
    for (size_t p = 0; p < particles.size(); ++p) {
        const string& particle = particles[p];
        int charge_bin_index = p + 1; // Z=2 -> 1, Z=3 -> 2, ...
        string full_name = prefix + particle + "_" + suffix;
        
        TH2F* h2 = dynamic_cast<TH2F*>(file->Get(full_name.c_str()));
        if (!h2) continue;

        // --- 循环 Rigidity Bin (自定义 5 个) ---
        for (size_t i = 0; i < N_RIG_BINS; ++i) {
            double rig_min = H5_RIG_BINS_EDGES[i];
            double rig_max = H5_RIG_BINS_EDGES[i+1];
            int rig_bin_index = i + 1;
            string rig_label = H5_RIG_LABELS[i];

            // 1. 投影到 1D 直方图
            int bin_y_min = h2->GetYaxis()->FindFixBin(rig_min);
            int bin_y_max = h2->GetYaxis()->FindFixBin(rig_max - 1e-6);
            if (i == N_RIG_BINS - 1) { bin_y_max = h2->GetNbinsY(); }

            TH1D* h1_proj = h2->ProjectionX(("proj_" + full_name + "_" + to_string(i)).c_str(), bin_y_min, bin_y_max);

            if (h1_proj->GetEntries() == 0) { delete h1_proj; continue; }

            // 2. Rebin 直到 Max > 60
            for(int r=1 ; r<=10; ++r) {
                if (h1_proj->GetMaximum() >= 60) break;
                h1_proj->Rebin(2);
            }

            // 3. 拟合
            double x_min_fit = 0.0, x_max_fit = 0.0;
            // H5 系列: 中心点设为 0.0
            findFitRange(h1_proj, 0.95, 0.0, x_min_fit, x_max_fit); 

            if (x_max_fit <= x_min_fit) { delete h1_proj; continue; }

            TF1* f_gaus = new TF1("f_gaus", "gaus", x_min_fit, x_max_fit);
            f_gaus->SetParameters(h1_proj->GetMaximum(), h1_proj->GetMean(), h1_proj->GetRMS());
            f_gaus->SetLineColor(kRed);

            int fit_status = h1_proj->Fit(f_gaus, "QRS");

            if (fit_status == 0) {
                double mean = f_gaus->GetParameter(1);
                double mean_err = f_gaus->GetParError(1);
                double sigma = f_gaus->GetParameter(2);
                double sigma_err = f_gaus->GetParError(2);
                double chi2 = f_gaus->GetChisquare();
                int ndf = f_gaus->GetNDF();
                double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                
                // 4. 记录数据到 TH2F
                h2_mean->SetBinContent(rig_bin_index, charge_bin_index, mean);
                h2_mean->SetBinError(rig_bin_index, charge_bin_index, mean_err);
                h2_sigma->SetBinContent(rig_bin_index, charge_bin_index, sigma);
                h2_sigma->SetBinError(rig_bin_index, charge_bin_index, sigma_err);

                // 5. 绘制拟合结果到 PDF
                c_fit->cd();
                string title = particle + " - " + info.output_suffix + " [" + rig_label + "]";
                h1_proj->SetTitle(title.c_str());
                h1_proj->GetXaxis()->SetTitle(x_axis_label.c_str());
                h1_proj->GetXaxis()->SetRangeUser(mean - 10*sigma, mean + 10*sigma);
                h1_proj->Draw("hist");
                f_gaus->Draw("SAME");
                
                TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
                stringstream ss; ss << fixed << setprecision(5);
                ss.str(""); ss << "#mu = " << mean << " #pm " << mean_err;
                latex.DrawLatex(0.6, 0.85, ss.str().c_str());
                ss.str(""); ss << "#sigma = " << sigma << " #pm " << sigma_err;
                latex.DrawLatex(0.6, 0.80, ss.str().c_str());
                ss.str(""); ss << fixed << setprecision(2) << "#chi^{2}/NDF = " << chi2_ndf;
                latex.DrawLatex(0.6, 0.75, ss.str().c_str());
                ss.str(""); ss << "Rig: " << rig_label;
                latex.DrawLatex(0.15, 0.85, ss.str().c_str());

                c_fit->Update();
                c_fit->Print(fit_pdf_path.c_str(), "pdf");
            }
            
            delete f_gaus;
            delete h1_proj;
        }
    }
    
    // PDF 文件结束写入
    c_fit->Print((fit_pdf_path + "]").c_str(), "pdf");
    delete c_fit;

    // 6. 最终绘图
    DrawH5ResultsFromTH2(h2_mean, h2_sigma, output_dir, info, suffix, output_name_base, x_axis_label);

    // 清理
    delete h2_mean;
    delete h2_sigma;
}

// --- 核心分析函数：H4 系列 (TH1F 直接拟合) ---
void AnalyzeH4(TFile* file, const string& suffix, const string& output_dir,
    const vector<string>& particles, const vector<int>& colors) {

    cout << "\n--- Analyzing H4 Series: " << suffix << " ---" << endl;
    gStyle->SetOptFit(0);
    
    const string prefix = "L1Inner_";
    const string x_axis_label = "1/#beta";
    HistInfo info = getHistInfo(suffix);
    string output_name_base = info.output_suffix + suffix.substr(3);
    
    // 引入 H4 结果存储 map
    H4ResultsMap h4_results;

    // 拟合结果 PDF 文件
    string fit_pdf_path = output_dir + "FitResults_" + output_name_base + ".pdf";
    TCanvas* c_fit = new TCanvas("c_fit_H4", "Gaussian Fit Results (H4)", 800, 600);
    c_fit->SetLogy(0);
    c_fit->Print((fit_pdf_path + "[").c_str(), "pdf");

    // --- 循环粒子 ---
    for (size_t p = 0; p < particles.size(); ++p) {
        const string& particle = particles[p];
        int color = colors[p];
        string full_name = prefix + particle + "_" + suffix;
        
        TH1F* h1 = dynamic_cast<TH1F*>(file->Get(full_name.c_str()));
        if (!h1) continue;

        TH1F* h1_clone = (TH1F*)h1->Clone(("clone_" + full_name).c_str());
        h1_clone->SetMarkerStyle(20);
        h1_clone->SetMarkerSize(1.2);
        //h1_clone->SetMarkerColor(color);

        // 1. Rebin 直到 Max > 60
            for(int r=1 ; r<=10; ++r) {
                if (h1_clone->GetMaximum() >= 60) break;
                h1_clone->Rebin(2);
            }

        // 2. 拟合
        double x_min_fit = 0.0, x_max_fit = 0.0;
        // H4 系列: 中心点设为 1.0
        findFitRange(h1_clone, 0.95, 1.0, x_min_fit, x_max_fit); 

        if (x_max_fit <= x_min_fit) { delete h1_clone; continue; }

        TF1* f_gaus = new TF1("f_gaus", "gaus", x_min_fit, x_max_fit);
        f_gaus->SetParameters(h1_clone->GetMaximum(), h1_clone->GetMean(), h1_clone->GetRMS());
        f_gaus->SetLineColor(kRed);

        int fit_status = h1_clone->Fit(f_gaus, "QRS");

        // 3. 绘制拟合结果到 PDF
        c_fit->cd();
        string title = particle + " - " + info.title_description;
        h1_clone->SetTitle(title.c_str());
        h1_clone->GetXaxis()->SetTitle(x_axis_label.c_str());
        h1_clone->GetXaxis()->SetRangeUser(f_gaus->GetParameter(1) - 10*f_gaus->GetParameter(2),
            f_gaus->GetParameter(1) + 10*f_gaus->GetParameter(2));
        h1_clone->Draw("hist");
        f_gaus->Draw("SAME");
        
        // 显示参数
        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035);
        stringstream ss; ss << fixed << setprecision(5);
        
        if (fit_status == 0) {
            double mean = f_gaus->GetParameter(1);
            double mean_err = f_gaus->GetParError(1);
            double sigma = f_gaus->GetParameter(2);
            double sigma_err = f_gaus->GetParError(2);
            double chi2 = f_gaus->GetChisquare();
            int ndf = f_gaus->GetNDF();
            double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0.0;
            
            // 记录数据到 H4ResultsMap
            h4_results[particle] = {mean, mean_err, sigma, sigma_err};

            ss.str(""); ss << "#mu = " << mean << " #pm " << mean_err;
            latex.DrawLatex(0.6, 0.85, ss.str().c_str());
            ss.str(""); ss << "#sigma = " << sigma << " #pm " << sigma_err;
            latex.DrawLatex(0.6, 0.80, ss.str().c_str());
            ss.str("");
            ss << fixed << setprecision(2) << "#chi^{2}/NDF = " << chi2_ndf;
            latex.DrawLatex(0.6, 0.75, ss.str().c_str());
        } else {
            latex.DrawLatex(0.6, 0.85, "Fit Failed");
        }

        c_fit->Update();
        c_fit->Print(fit_pdf_path.c_str(), "pdf");
        
        delete f_gaus;
        delete h1_clone;
    }
    
    // PDF 文件结束写入
    c_fit->Print((fit_pdf_path + "]").c_str(), "pdf");
    delete c_fit;

    // 4. 最终绘图 (H4 Mean/Sigma vs Charge)
    DrawH4Results(h4_results, output_dir, info, output_name_base, x_axis_label, particles, colors);
}


// --- 主函数 ---
void BetaFit_old() {
    gROOT->SetBatch(kTRUE);
    
    // 强制设置 TGraph/TGraphErrors 的 X 轴误差为 0
    gStyle->SetErrorX(0); 

    const string input_file_path = "/eos/user/z/zixuan/Isotope/Add/Be_frag4.root";
    const string output_dir = "/eos/user/z/zixuan/Isotope/Beta/";
    
    // 粒子和颜色定义
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

    // --- 1. H5 系列 (TH2F) 分析 ---
    vector<string> h5_suffixes = {"ID_H5a", "ID_H5b"};
    for (const string& suffix : h5_suffixes) {
        AnalyzeH5(file, suffix, output_dir, particles);
    }

    // --- 2. H4 系列 (TH1F) 分析 ---
    vector<string> h4_suffixes = {"ID_H4a", "ID_H4b"};
    for (const string& suffix : h4_suffixes) {
        AnalyzeH4(file, suffix, output_dir, particles, colors);
    }

    file->Close();
    delete file;
    gROOT->SetBatch(kFALSE);
}