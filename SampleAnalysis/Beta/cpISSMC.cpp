#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

using namespace std;

struct ElementInfo {
    int Z;
    string element_label;
    int iss_index; // Z - 2
    vector<string> mc_nuclides;
};

const map<string, int> MC_NUCLIDE_INDICES = {
    {"B10", 0}, {"B11", 1}, {"Be10", 2}, {"Be7", 3}, 
    {"Be9", 4}, {"C12", 5}, {"N15", 6}, {"O16", 7}
};

const vector<ElementInfo> ALL_ELEMENTS = {
    {4, "Be", 2, {"Be7", "Be9", "Be10"}},
    {5, "B", 3, {"B10", "B11"}},
    {6, "C", 4, {"C12"}},
    {7, "N", 5, {"N15"}},
    {8, "O", 6, {"O16"}}
};

void DrawAndSaveComparison(TFile* iss_file, TFile* mc_file, 
                           const ElementInfo& element_info, 
                           const string& h5_suffix, const string& type, 
                           const string& output_dir) {
    
    string canvas_name = "c_" + element_info.element_label + "_" + h5_suffix + "_" + type;
    gROOT->cd(); 
    TCanvas* c = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 800, 600);
    c->SetGrid();

    string rich_label = (h5_suffix == "ID_H5a") ? "NaF-Tracker #Delta(1/#beta)" : "AGL-Tracker #Delta(1/#beta)";
    string y_title = (type == "Mean") ? "#mu_{#Delta(1/#beta)}" : "#sigma_{#Delta(1/#beta)}";
    string graph_name_base = "g_Rigidity" + h5_suffix.substr(3) + "_" + type + "_vs_Rigidity";

    TLegend* leg = new TLegend(0.70, 0.70, 0.95, 0.95);
    leg->SetFillStyle(0);
    leg->SetBorderSize(1);

    vector<TGraphErrors*> graphs;
    vector<string> labels;
    vector<int> colors = {kRed, kBlue, kGreen + 2, kMagenta, kOrange + 1, kViolet, kBlack};
    int color_idx = 0;

    // 1. Get ISS Graph (Data)
    string iss_graph_name = graph_name_base + "_" + to_string(element_info.iss_index);
    TGraphErrors* g_iss = dynamic_cast<TGraphErrors*>(iss_file->Get(iss_graph_name.c_str()));
    
    if (g_iss) {
        g_iss->SetMarkerStyle(20);
        g_iss->SetMarkerColor(kBlack); 
        g_iss->SetLineColor(kBlack);
        graphs.push_back(g_iss);
        labels.push_back("ISS Data " + element_info.element_label);
    }

    // 2. Get MC Graphs (Simulation)
    for (const string& nuclide : element_info.mc_nuclides) {
        if (!MC_NUCLIDE_INDICES.count(nuclide)) continue;
        int mc_index = MC_NUCLIDE_INDICES.at(nuclide);
        string mc_graph_name = "g_Rigidity" + h5_suffix.substr(3) + "_MC_Combined_" + type + "_vs_Rigidity_" + to_string(mc_index);
        
        TGraphErrors* g_mc = dynamic_cast<TGraphErrors*>(mc_file->Get(mc_graph_name.c_str()));
        
        if (g_mc) {
            int color = colors[color_idx % colors.size()];
            g_mc->SetMarkerStyle(20);
            g_mc->SetMarkerColor(color);
            g_mc->SetLineColor(color);
            graphs.push_back(g_mc);
            labels.push_back("MC " + nuclide);
            color_idx++;
        }
    }

    if (graphs.empty()) {
        cout << "Warning: No graphs found for " << element_info.element_label << " " << h5_suffix << " " << type << endl;
        return;
    }

    // Draw
    double x_min = 30.0;
    double x_max = 240.0;

    double y_min, y_max;
    // 使用 DrawAndSaveGraphs 里的逻辑获取全局 Y 范围 (为了代码简洁，直接在 Draw 逻辑中实现)
    y_min = 1e10; y_max = -1e10;
    bool first_draw = true;
    for (TGraphErrors* g : graphs) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            if (x >= x_min && x <= x_max) {
                double y_low = y - g->GetErrorY(i);
                double y_high = y + g->GetErrorY(i);
                y_min = min(y_min, y_low);
                y_max = max(y_max, y_high);
            }
        }
    }
    double range = y_max - y_min;
    double buffer = (range == 0.0) ? 0.01 : range * 0.10; 
    y_min -= buffer;
    y_max += buffer;
    if (y_min >= y_max) { y_min = -0.01; y_max = 0.01; }

    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->Draw(first_draw ? "AP" : "P same");
        if (first_draw) {
            graphs[i]->SetTitle((element_info.element_label + " " + type + " vs Rigidity;Rigidity [GV];" + y_title).c_str());
            graphs[i]->GetXaxis()->SetTitle("Rigidity [GV]");
            graphs[i]->GetYaxis()->SetTitle(y_title.c_str());
            graphs[i]->GetXaxis()->SetRangeUser(x_min, x_max);
            graphs[i]->GetYaxis()->SetRangeUser(y_min, y_max);
            first_draw = false;
        }
        leg->AddEntry(graphs[i], labels[i].c_str(), "lp");
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.15, 0.85, rich_label.c_str());

    leg->Draw();
    c->SaveAs((output_dir + element_info.element_label + "_" + h5_suffix + "_" + type + "_Comp.png").c_str());
    
    // 不手动 delete TGraphErrors, TCanvas, TLegend
}

void cpISSMC() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetErrorX(0); 

    const string iss_path = "/eos/user/z/zixuan/Isotope/Beta/ISS/ISS_RICHBetaStudy.root";
    const string mc_path = "/eos/user/z/zixuan/Isotope/Beta/MC/MC_RICHBetaStudy.root";
    const string output_dir = "/eos/user/z/zixuan/Isotope/Beta/";

    if (gSystem->AccessPathName(output_dir.c_str())) {
        gSystem->mkdir(output_dir.c_str(), kTRUE);
    }

    TFile* iss_file = TFile::Open(iss_path.c_str(), "READ");
    TFile* mc_file = TFile::Open(mc_path.c_str(), "READ");

    if (!iss_file || iss_file->IsZombie()) {
        cerr << "ERROR: Cannot open ISS file: " << iss_path << endl;
        return;
    }
    if (!mc_file || mc_file->IsZombie()) {
        cerr << "ERROR: Cannot open MC file: " << mc_path << endl;
        iss_file->Close();
        return;
    }

    vector<string> suffixes = {"ID_H5a", "ID_H5b"};
    vector<string> types = {"Mean", "Sigma"};

    for (const auto& element : ALL_ELEMENTS) {
        for (const string& suffix : suffixes) {
            for (const string& type : types) {
                DrawAndSaveComparison(iss_file, mc_file, element, suffix, type, output_dir);
            }
        }
    }

    iss_file->Close();
    mc_file->Close();
    
    // 不手动 delete TFile*
    gROOT->SetBatch(kFALSE);
    cout << "Comparison plots finished. Saved to: " << output_dir << endl;
}