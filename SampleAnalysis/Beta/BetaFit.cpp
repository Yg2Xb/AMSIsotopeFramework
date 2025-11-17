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
#include <TLatex.h>   // 新增
#include <iomanip>   // 新增
#include <sstream>   // 新增

using namespace std;

// --- 辅助结构和函数定义 (保持不变) ---

void findFitRange(TH1* hist, double coverage, double& x_min, double& x_max) {
    if (!hist || hist->GetEntries() == 0) return;
    
    int max_bin = hist->GetMaximumBin();
    double total_integral = hist->Integral();
    if (total_integral <= 0) return;

    double required_integral = total_integral * coverage;
    double current_integral = hist->GetBinContent(max_bin);
    
    x_min = hist->GetXaxis()->GetBinLowEdge(max_bin);
    x_max = hist->GetXaxis()->GetBinUpEdge(max_bin);

    int low_bin = max_bin - 1;
    int high_bin = max_bin + 1;
    
    while (current_integral < required_integral) {
        bool extended = false;
        double content_low = (low_bin >= 1) ? hist->GetBinContent(low_bin) : 0;
        double content_high = (high_bin <= hist->GetNbinsX()) ? hist->GetBinContent(high_bin) : 0;

        if (content_low >= content_high && low_bin >= 1 && current_integral + content_low <= total_integral) {
            current_integral += content_low;
            x_min = hist->GetXaxis()->GetBinLowEdge(low_bin);
            low_bin--;
            extended = true;
        } else if (content_high > content_low && high_bin <= hist->GetNbinsX() && current_integral + content_high <= total_integral) {
            current_integral += content_high;
            x_max = hist->GetXaxis()->GetBinUpEdge(high_bin);
            high_bin++;
            extended = true;
        } else if (low_bin >= 1 && current_integral + content_low <= total_integral) {
            current_integral += content_low;
            x_min = hist->GetXaxis()->GetBinLowEdge(low_bin);
            low_bin--;
            extended = true;
        } else if (high_bin <= hist->GetNbinsX() && current_integral + content_high <= total_integral) {
            current_integral += content_high;
            x_max = hist->GetXaxis()->GetBinUpEdge(high_bin);
            high_bin++;
            extended = true;
        }

        if (!extended) break;
    }
}

struct HistInfo {
    string y_axis_label;
    string output_suffix;
    string title_description;
};

struct PlotConfig {
    HistInfo info;
    int skip_points;
    double x_min_plot;
    double x_max_plot;
};

HistInfo getHistInfo(const string& suffix) {
    HistInfo info;
    
    if (suffix == "ID_H5a") { info = {"Rigidity [GV]", "Rigidity", "NaF-Tracker #Delta(1/#beta)"}; } 
    else if (suffix == "ID_H5b") { info = {"Rigidity [GV]", "Rigidity", "AGL-Tracker #Delta(1/#beta)"}; } 
    else if (suffix == "ID_H5a2") { info = {"gene Rigidity [GV]", "geneRigidity", "NaF-Tracker #Delta(1/#beta)"}; } 
    else if (suffix == "ID_H5b2") { info = {"gene Rigidity [GV]", "geneRigidity", "AGL-Tracker #Delta(1/#beta)"}; } 
    else if (suffix == "ID_H5a3") { info = {"gene Rigidity [GV]", "geneRigidity", "NaF-true #Delta(1/#beta)"}; } 
    else if (suffix == "ID_H5b3") { info = {"gene Rigidity [GV]", "geneRigidity", "AGL-True #Delta(1/#beta)"}; } 
    else if (suffix == "ID_H6a") { info = {"NaF E_{k}/n [GeV/n]", "EkPerN", "TOF-NaF #Delta(1/#beta)"}; } 
    else if (suffix == "ID_H6b") { info = {"AGL E_{k}/n [GeV/n]", "EkPerN", "TOF-AGL #Delta(1/#beta)"}; } 
    else if (suffix == "ID_H7a") { info = {"NaF #beta Rig[GV]", "BetaRig", "TOF-NaF #Delta(1/#beta)"}; } 
    else if (suffix == "ID_H7b") { info = {"AGL #beta Rig[GV]", "BetaRig", "TOF-AGL #Delta(1/#beta)"}; } 
    else { info = {"Y Variable", "YVariable", "Unknown Delta Beta"}; }
    
    if (info.y_axis_label.rfind("NaF ", 0) == 0) { info.y_axis_label = info.y_axis_label.substr(4); } 
    else if (info.y_axis_label.rfind("AGL ", 0) == 0) { info.y_axis_label = info.y_axis_label.substr(4); }
    
    return info;
}

PlotConfig getPlotConfig(const string& suffix) {
    PlotConfig config;
    config.info = getHistInfo(suffix);
    config.x_min_plot = 0.1; 
    config.x_max_plot = 1000.0; 
    config.skip_points = 0;   

    // 1. 设置跳过点数 (保持注释，使用用户提供的版本)
    //if (suffix.find('a') != string::npos) { config.skip_points = 2; } 
    //else if (suffix.find('b') != string::npos) { config.skip_points = 6; } 

    // 2. 设置 X 轴范围 (Rigidity/Ek/n) (保持用户提供的版本)
    if (suffix == "ID_H7a") {
        config.x_min_plot = 2.0; config.x_max_plot = 20.0;
    } else if (suffix == "ID_H7b") {
        config.x_min_plot = 4.0; config.x_max_plot = 20.0;
    } else if (suffix == "ID_H6a") {
        config.x_min_plot = 0.8; config.x_max_plot = 20.0;
    } else if (suffix == "ID_H6b") {
        config.x_min_plot = 2.0; config.x_max_plot = 30.0;
    } else if (suffix.find("ID_H5a") != string::npos) {
        config.x_min_plot = 20.; config.x_max_plot = 150.0;
    } else if (suffix.find("ID_H5b") != string::npos) {
        config.x_min_plot = 20; config.x_max_plot = 200.0;
    }
    
    return config;
}

// 计算 TGraph 集合在特定 X 范围内的全局 Y 范围 (保持不变)
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
        
        double buffer = range * 0.10; // 使用用户提供的 10% 扩充
        y_min -= buffer;
        y_max += buffer;
    } else {
        y_min = -0.01;
        y_max = 0.01;
    }
}

// 图例位置函数 (保持不变)
void setLegendPosition(const string& suffix, bool is_mean, double& x1, double& y1, double& x2, double& y2) {
    // 默认图例大小 (相对坐标)
    const double w = 0.25; // 宽度
    const double h = 0.20; // 高度

    // 默认图例边距
    const double margin = 0.03;

    if (suffix.find("ID_H5") != string::npos) {
        // H5 系列 (5a, 5b, 5a2, 5b2, 5a3, 5b3)
        if (is_mean) { // Mean: 右下
            x1 = 1.0 - 2*margin - w; y1 = 5*margin;
            x2 = 1.0 - 2*margin;     y2 = 5*margin + h;
        } else { // Sigma: 右上
            x1 = 1.0 - 2*margin - w; y1 = 1.0 - 5*margin - h;
            x2 = 1.0 - 2*margin;     y2 = 1.0 - 5*margin;
        }
    } else if (suffix.find("ID_H6") != string::npos) {
        // H6 系列 (6a, 6b)
        if (is_mean) { // Mean: 右下
            x1 = 1.0 - 2*margin - w; y1 = 5*margin;
            x2 = 1.0 - 2*margin;     y2 = 5*margin + h;
        } else { // Sigma: 左上
            x1 = 5*margin;           y1 = 1.0 - 5*margin - h;
            x2 = 5*margin + w;       y2 = 1.0 - 5*margin;
        }
    } else if (suffix == "ID_H7a") {
        // H7a: Mean 右下, Sigma 右上
        if (is_mean) { // Mean: 右下
            x1 = 1.0 - 2*margin - w; y1 = 5*margin;
            x2 = 1.0 - 2*margin;     y2 = 5*margin + h;
        } else { // Sigma: 右上
            x1 = 1.0 - 2*margin - w; y1 = 1.0 - 5*margin - h;
            x2 = 1.0 - 2*margin;     y2 = 1.0 - 5*margin;
        }
    } else if (suffix == "ID_H7b") {
        // H7b: Mean 右下, Sigma 左上
        if (is_mean) { // Mean: 右下
            x1 = 1.0 - 2*margin - w; y1 = 5*margin;
            x2 = 1.0 - 2*margin;     y2 = 5*margin + h;
        } else { // Sigma: 左上
            x1 = 5*margin;           y1 = 1.0 - 5*margin - h;
            x2 = 5*margin + w;       y2 = 1.0 - 5*margin;
        }
    } else {
        // 默认: 右上 
        x1 = 0.7; y1 = 0.7;
        x2 = 0.9; y2 = 0.9;
    }
}


// --- 主分析函数 ---
void BetaFit() {
    gROOT->SetBatch(kTRUE); 
    
    const string input_file_path = "/eos/user/z/zixuan/Isotope/Add/Be_frag4.root";
    const string output_dir = "/eos/user/z/zixuan/Isotope/Beta/";
    const string prefix = "L1Inner_";
    const string x_axis_label = "#Delta(1/#beta)"; 
    
    vector<string> particles = {"Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen"}; 
    vector<int> colors = {kBlack, kRed, kBlue, kGreen + 2, kOrange + 1}; 
    vector<string> id_suffixes = {
        "ID_H5a", "ID_H5b", "ID_H5a2", "ID_H5b2", "ID_H5a3", "ID_H5b3",
        "ID_H6a", "ID_H6b",
        "ID_H7a", "ID_H7b"
    };

    if (gSystem->AccessPathName(output_dir.c_str())) {
        cout << "Creating output directory: " << output_dir << endl;
        gSystem->mkdir(output_dir.c_str(), kTRUE);
    }

    TFile* file = TFile::Open(input_file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
        cerr << "ERROR: Cannot open input file: " << input_file_path << endl;
        return;
    }

    // **开始循环处理**
    for (const string& suffix : id_suffixes) {
        
        PlotConfig config = getPlotConfig(suffix);
        const HistInfo& info = config.info;
        string output_name_base = info.output_suffix + suffix.substr(3); 
        
        // **新增: 拟合结果输出 PDF 文件名**
        string fit_pdf_path = output_dir + "FitResults_" + output_name_base + ".pdf";
        TCanvas* c_fit = new TCanvas("c_fit", "Gaussian Fit Results", 800, 600);
        
        // **PDF 文件开启写入 (使用 [.c_str() 结束方括号)]**
        c_fit->Print((fit_pdf_path + "[").c_str(), "pdf"); 

        TCanvas* c_mean = new TCanvas(("c_mean_" + output_name_base).c_str(), "Mean Plot", 800, 600);
        TCanvas* c_sigma = new TCanvas(("c_sigma_" + output_name_base).c_str(), "Sigma Plot", 800, 600);
        
        // 动态创建和定位图例 (保持不变)
        double mean_x1, mean_y1, mean_x2, mean_y2;
        setLegendPosition(suffix, true, mean_x1, mean_y1, mean_x2, mean_y2);
        TLegend* leg_mean = new TLegend(mean_x1, mean_y1, mean_x2, mean_y2);
        leg_mean->SetFillStyle(0); leg_mean->SetBorderSize(0); 

        double sigma_x1, sigma_y1, sigma_x2, sigma_y2;
        setLegendPosition(suffix, false, sigma_x1, sigma_y1, sigma_x2, sigma_y2);
        TLegend* leg_sigma = new TLegend(sigma_x1, sigma_y1, sigma_x2, sigma_y2);
        leg_sigma->SetFillStyle(0); leg_sigma->SetBorderSize(0);

        c_mean->SetGrid(); c_sigma->SetGrid();
        c_mean->SetLogx(1); c_sigma->SetLogx(1);

        vector<TGraphErrors*> graphs_mean;
        vector<TGraphErrors*> graphs_sigma;

        // --- 提取数据并拟合 ---
        for (size_t p = 0; p < particles.size(); ++p) {
            const string& particle = particles[p];
            int color = colors[p];
            string full_name = prefix + particle + "_" + suffix;
            
            TH2D* h2 = dynamic_cast<TH2D*>(file->Get(full_name.c_str()));
            
            if (!h2) continue;

            TH2D* h2_rebin = (TH2D*)h2->Clone((full_name + "_rebin").c_str());
            h2_rebin->RebinY(2);
            int new_bins_y = h2_rebin->GetNbinsY(); 
            
            if (new_bins_y <= config.skip_points) {
                 delete h2_rebin;
                 continue;
            }

            TGraphErrors* g_mean = new TGraphErrors(new_bins_y);
            TGraphErrors* g_sigma = new TGraphErrors(new_bins_y);

            g_mean->SetMarkerStyle(20); g_sigma->SetMarkerStyle(20);
            g_mean->SetMarkerSize(1.0); g_sigma->SetMarkerSize(1.0);
            g_mean->SetMarkerColor(color); g_sigma->SetMarkerColor(color);
            g_mean->SetLineColor(color); g_sigma->SetLineColor(color);
            
            for (int i = 1; i <= new_bins_y; ++i) {
                
                if (i <= config.skip_points) { continue; }

                double y_center = h2_rebin->GetYaxis()->GetBinCenter(i);
                TH1D* h1_proj = h2_rebin->ProjectionX(("proj_" + full_name + "_" + to_string(i)).c_str(), i, i);

                // 检查 Y 轴范围是否在用户设定的 X 轴绘图范围内（只拟合要绘制的点）
                if (y_center < config.x_min_plot || y_center > config.x_max_plot) {
                    delete h1_proj; 
                    continue;
                }

                if (h1_proj->GetEntries() < 50 || h1_proj->GetMaximum() <= 0) { 
                    delete h1_proj; continue;
                }

                double x_min_fit = 0.0, x_max_fit = 0.0;
                findFitRange(h1_proj, 0.95, x_min_fit, x_max_fit);

                if (x_max_fit <= x_min_fit) { delete h1_proj; continue; }

                TF1* f_gaus = new TF1("f_gaus", "gaus", x_min_fit, x_max_fit);
                f_gaus->SetParameters(h1_proj->GetMaximum(), h1_proj->GetMean(), h1_proj->GetRMS());
                f_gaus->SetLineColor(kRed);

                int fit_status = h1_proj->Fit(f_gaus, "QRS"); // "S" 选项获取拟合结果

                if (fit_status == 0) { 
                    double mean = f_gaus->GetParameter(1);
                    double mean_err = f_gaus->GetParError(1);
                    double sigma = f_gaus->GetParameter(2);
                    double sigma_err = f_gaus->GetParError(2);
                    double chi2 = f_gaus->GetChisquare();
                    int ndf = f_gaus->GetNDF();
                    double chi2_ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                    
                    // 1. 记录数据点 (保持不变)
                    int point_index = g_mean->GetN();
                    g_mean->SetPoint(point_index, y_center, mean);
                    g_mean->SetPointError(point_index, 0.0, mean_err); 
                    g_sigma->SetPoint(point_index, y_center, sigma);
                    g_sigma->SetPointError(point_index, 0.0, sigma_err); 

                    // 2. **新增: 绘制拟合结果到 PDF**
                    c_fit->cd();
                    
                    // **修正：删除了 SetTitle 外部多余的括号**
                    // 设置直方图标题，包括 Rigidity/Ek/n 范围
                    h1_proj->SetTitle((particle + " - " + info.output_suffix + " [" + to_string(static_cast<int>(y_center)) + "]").c_str());
                    h1_proj->GetXaxis()->SetRangeUser(mean - 10*sigma, mean + 10*sigma);
                    h1_proj->SetMarkerSize(1.0);
                    h1_proj->SetMarkerStyle(20);
                    h1_proj->Draw("PZ"); 
                    f_gaus->Draw("SAME");
                    
                    // **显示参数**
                    TLatex latex;
                    latex.SetNDC();
                    latex.SetTextSize(0.035);

                    stringstream ss;
                    ss << fixed << setprecision(5);
                    
                    // Mean
                    ss.str(""); ss << "#mu = " << mean << " #pm " << mean_err;
                    latex.DrawLatex(0.6, 0.85, ss.str().c_str());
                    
                    // Sigma
                    ss.str(""); ss << "#sigma = " << sigma << " #pm " << sigma_err;
                    latex.DrawLatex(0.6, 0.80, ss.str().c_str());

                    // Chi2/NDF
                    ss.str(""); 
                    ss << fixed << setprecision(2) << "#chi^{2}/NDF = " << chi2_ndf;
                    latex.DrawLatex(0.6, 0.75, ss.str().c_str());

                    // Y Center (Rigidity/Ek/n)
                    ss.str(""); 
                    ss << info.y_axis_label << ": " << y_center;
                    latex.DrawLatex(0.15, 0.85, ss.str().c_str());

                    c_fit->Update();
                    c_fit->Print(fit_pdf_path.c_str(), "pdf"); // 输出当前页

                }
                
                delete f_gaus;
                delete h1_proj;
            } 

            if (g_mean->GetN() > 0) {
                 graphs_mean.push_back(g_mean);
                 graphs_sigma.push_back(g_sigma);
                 leg_mean->AddEntry(g_mean, particle.c_str(), "p"); 
                 leg_sigma->AddEntry(g_sigma, particle.c_str(), "p");
            } else {
                 delete g_mean; delete g_sigma;
            }

            delete h2_rebin; 
        } 
        
        // **PDF 文件结束写入**
        c_fit->Print((fit_pdf_path + "]").c_str(), "pdf"); 
        delete c_fit;

        // 5. **绘图输出 (Mean/Sigma) (保持不变)**
        double plot_x_min = config.x_min_plot;
        double plot_x_max = config.x_max_plot; 
        double mean_y_min, mean_y_max;
        double sigma_y_min, sigma_y_max;
        
        getGlobalYRange(graphs_mean, plot_x_min, plot_x_max, mean_y_min, mean_y_max);
        getGlobalYRange(graphs_sigma, plot_x_min, plot_x_max, sigma_y_min, sigma_y_max);


        // Mean 绘图
        c_mean->cd();
        for (size_t i = 0; i < graphs_mean.size(); ++i) {
            string draw_opt = (i == 0) ? "A P" : "P SAME"; 
            graphs_mean[i]->Draw(draw_opt.c_str());
            
            if (i == 0) {
                string full_title = info.title_description + " Mean vs " + info.y_axis_label;
                graphs_mean[i]->SetTitle(full_title.c_str());
                graphs_mean[i]->GetXaxis()->SetTitle(info.y_axis_label.c_str());
                graphs_mean[i]->GetYaxis()->SetTitle(("#mu_{" + x_axis_label + "}").c_str());
                
                graphs_mean[i]->GetXaxis()->SetRangeUser(plot_x_min, plot_x_max); 
                graphs_mean[i]->GetYaxis()->SetRangeUser(mean_y_min, mean_y_max);
                c_mean->Update(); 
            }
        }
        if (graphs_mean.size() > 0) {
            leg_mean->Draw();
            string mean_output_path = output_dir + output_name_base + "_Mean.png";
            c_mean->SaveAs(mean_output_path.c_str());
        }

        // Sigma 绘图
        c_sigma->cd();
        for (size_t i = 0; i < graphs_sigma.size(); ++i) {
            string draw_opt = (i == 0) ? "A P" : "P SAME";
            graphs_sigma[i]->Draw(draw_opt.c_str());
            
            if (i == 0) {
                string full_title = info.title_description + " Sigma vs " + info.y_axis_label;
                graphs_sigma[i]->SetTitle(full_title.c_str());
                graphs_sigma[i]->GetXaxis()->SetTitle(info.y_axis_label.c_str());
                graphs_sigma[i]->GetYaxis()->SetTitle(("#sigma_{" + x_axis_label + "}").c_str());
                
                graphs_sigma[i]->GetXaxis()->SetRangeUser(plot_x_min, plot_x_max); 
                graphs_sigma[i]->GetYaxis()->SetRangeUser(sigma_y_min, sigma_y_max);
                c_sigma->Update(); 
            }
        }
        if (graphs_sigma.size() > 0) {
            leg_sigma->Draw();
            string sigma_output_path = output_dir + output_name_base + "_Sigma.png";
            c_sigma->SaveAs(sigma_output_path.c_str());
        }
        
        // 清理 (Mean/Sigma)
        for (auto g : graphs_mean) delete g;
        for (auto g : graphs_sigma) delete g;
        delete c_mean; delete c_sigma;
        delete leg_mean; delete leg_sigma;

    } // End of suffix loop

    file->Close();
    delete file;
    gROOT->SetBatch(kFALSE);
}