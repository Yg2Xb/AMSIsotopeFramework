#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TAxis.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm> // For std::sort and std::unique

// --- 辅助函数 ---

// 通用画布设置函数
void setupCanvas(TCanvas* canvas) {
    canvas->SetLeftMargin(0.12); // 稍微加宽以容纳Y轴标签
    canvas->SetRightMargin(0.03);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.15);
}
    auto setupHistogram = [](TH1D* hist) {
        // 去掉X和Y轴标题
        hist->GetXaxis()->SetTitle("");
        hist->GetYaxis()->SetTitle("");
        
        // 调大标签大小
        hist->GetXaxis()->SetLabelSize(0.1);
        hist->GetYaxis()->SetLabelSize(0.1);
        
        // 设置Y轴刻度分隔为505
        hist->GetYaxis()->SetNdivisions(505);
    };


// 数据拷贝辅助函数：将源直方图在指定刚度范围内的数据拷贝到目标直方图
void copyData(TH1D* target, TH1D* source, double minR, double maxR) {
    if (!source) return;
    for (int i = 1; i <= source->GetNbinsX(); ++i) {
        double rigidity = source->GetBinCenter(i);
        if (rigidity >= minR && rigidity < maxR) {
            int targetBin = target->FindBin(rigidity);
            target->SetBinContent(targetBin, source->GetBinContent(i));
            target->SetBinError(targetBin, source->GetBinError(i));
        }
    }
}


// --- 主绘图函数 ---

void DrawFrac() {
    // --- 1. 参数设置 ---
    const std::string sourceElement = "Carbon";
    const std::string fragElement   = "Boron";
    const std::string trackerLayer  = "L1Inner";
    const std::string outputPath = "/eos/user/z/zixuan/Isotope/ChargeTemp/";

    // 设置全局ROOT样式
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // --- 2. 元素和样式定义 ---
    // 根据源元素确定左右邻居
    std::string elementCenter = sourceElement;
    std::string elementLeft, elementRight;
    std::map<std::string, int> chargeMap;

    if (sourceElement == "Boron") {
        elementLeft = "Beryllium";
        elementRight = "Carbon";
        chargeMap[elementLeft] = 4;
        chargeMap[elementCenter] = 5;
        chargeMap[elementRight] = 6;
    } else if (sourceElement == "Beryllium") {
        elementLeft = "Lithium";
        elementRight = "Boron";
        chargeMap[elementLeft] = 3;
        chargeMap[elementCenter] = 4;
        chargeMap[elementRight] = 5;
    } else if (sourceElement == "Carbon") {
        elementLeft = "Boron";
        elementRight = "Nitrogen";
        chargeMap[elementLeft] = 5;
        chargeMap[elementCenter] = 6;
        chargeMap[elementRight] = 7;
    } else {
        std::cerr << "Error: Unknown source element '" << sourceElement << "'. Please add it to the logic." << std::endl;
        return;
    }
    
    const std::vector<std::string> elementsToPlot = {elementLeft, elementCenter, elementRight};

    // --- 3. 打开输入文件 ---
    TString inputFilePath = TString::Format(
        "root://eoshome-z.cern.ch//eos/user/z/zixuan/Isotope/ChargeTemp/QFit_%s_to_%s_%s.root",
        sourceElement.c_str(),
        fragElement.c_str(),
        trackerLayer.c_str()
    );
    
    std::cout << "Processing file: " << inputFilePath << std::endl;
    
    TFile* file = TFile::Open(inputFilePath, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << inputFilePath << std::endl;
        return;
    }

    // --- 4. 循环处理每个元素，生成拼接图 ---
    for (const auto& element : elementsToPlot) {
        
        std::cout << "\nProcessing element: " << element << std::endl;

        // --- 4.1 获取三个探测器的原始直方图 ---
        TH1D *h_tof = nullptr, *h_naf = nullptr, *h_agl = nullptr;
        file->GetObject(TString::Format("h_narrowfrac_%s_TOF", element.c_str()), h_tof);
        file->GetObject(TString::Format("h_narrowfrac_%s_NaF", element.c_str()), h_naf);
        file->GetObject(TString::Format("h_narrowfrac_%s_AGL", element.c_str()), h_agl);

        if (!h_tof || !h_naf || !h_agl) {
            std::cerr << "Warning: Missing one or more detector histograms for element '" << element << "'. Skipping." << std::endl;
            continue;
        }

        // --- 4.2 创建一个新的、包含所有bins的组合直方图 ---
        // 收集所有原始直方图的bin边界
        std::vector<double> bin_edges;
        TAxis* axes[] = {h_tof->GetXaxis(), h_naf->GetXaxis(), h_agl->GetXaxis()};
        for (TAxis* axis : axes) {
            for (int i = 1; i <= axis->GetNbins(); ++i) {
                bin_edges.push_back(axis->GetBinLowEdge(i));
            }
            bin_edges.push_back(axis->GetBinUpEdge(axis->GetNbins()));
        }
        // 排序并移除重复的bin边界
        std::sort(bin_edges.begin(), bin_edges.end());
        bin_edges.erase(std::unique(bin_edges.begin(), bin_edges.end()), bin_edges.end());

        // 创建组合直方图
        TString combinedHistName = TString::Format("h_combined_narrowfrac_%s", element.c_str());
        TH1D* h_combined = new TH1D(combinedHistName, "", bin_edges.size() - 1, &bin_edges[0]);

        // --- 4.3 根据刚度范围，从原始直方图填充组合直方图 ---
        // TOF: R < 1.28 (包含1.11之前的所有部分)
        copyData(h_combined, h_tof, 0.0, 1.28);
        // NaF: 1.28 <= R < 3.06
        copyData(h_combined, h_naf, 1.28, 3.06);
        // AGL: R >= 3.06
        copyData(h_combined, h_agl, 3.06, 1e9); // 使用一个很大的数作为上限

        // --- 4.4 设置样式和Y轴范围 ---
        // 颜色规则: 电荷从小到大 -> 绿, 蓝, 品红
        int color = kBlack; // 默认颜色
        int charge = chargeMap.count(element) ? chargeMap[element] : 0;
        if (charge == chargeMap[elementLeft]) color = kGreen + 2;
        else if (charge == chargeMap[elementCenter]) color = kBlue;
        else if (charge == chargeMap[elementRight]) color = kMagenta;
        
        h_combined->SetLineColor(color);
        h_combined->SetMarkerColor(color);
        h_combined->SetMarkerStyle(20);
        h_combined->SetMarkerSize(1.2);

        h_combined->GetXaxis()->SetRangeUser(0.41, 21);
        // Y轴范围规则
        if (element == sourceElement) {
            h_combined->GetYaxis()->SetRangeUser(0.985, 1.005);
        } else if (element == "Beryllium") {
            h_combined->GetYaxis()->SetRangeUser(0, 0.006);
        } else if (element == "Carbon" || element == "Nitrogen") {
            h_combined->GetYaxis()->SetRangeUser(0, 0.0005);
        } else {
            // 为其他可能的元素设置一个默认范围
            h_combined->GetYaxis()->SetRangeUser(0, 0.01);
        }

        // --- 4.5 创建画布并绘制 ---
        TString canvasName = TString::Format("canvas_combined_%s", element.c_str());
        TCanvas* canvas = new TCanvas(canvasName, "", 900, 300); 
        setupCanvas(canvas);
        
        // 应用通用样式，并传入元素名以设置坐标轴标题
        setupHistogram(h_combined);
        
        h_combined->Draw("p0");

        // --- 4.6 保存图像 ---
        TString outputFileName = TString::Format(
            "%s/Combined_Frac_%s_to_%s_%s.png",
            outputPath.c_str(),
            sourceElement.c_str(),
            fragElement.c_str(),
            element.c_str()
        );
        
        canvas->SaveAs(outputFileName);
        std::cout << "  -> Saved combined plot: " << outputFileName << std::endl;

        // 清理内存
        delete canvas;
        delete h_combined;
    }

    // --- 5. 清理 ---
    file->Close();
    delete file;

    std::cout << "\nAll combined fraction plots saved successfully." << std::endl;
}