#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TLegend.h> // 需要 TLedgend
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>

// --- 辅助函数 ---

// 通用画布设置函数 (LogY)
void setupCanvas(TCanvas* canvas) {
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.1);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.15);
    canvas->SetLogy(); // *** 启用 Log Y 轴 ***
}

// 通用直方图样式和轴设置函数
void setupHistogramForYield(TH1D* hist, const std::string& element, bool isFirst) {
    // 调大标签大小
    hist->GetXaxis()->SetLabelSize(0.045);
    hist->GetYaxis()->SetLabelSize(0.045);
    
    // 设置 X, Y 轴标题
    if (isFirst) {
        hist->GetXaxis()->SetTitle("Measured Ek/n(GeV/n)");
        hist->GetYaxis()->SetTitle("Events in L1Q 4.8-5.4");
    } else {
        // 非第一个绘制的直方图，不显示标题，避免重复
        hist->GetXaxis()->SetTitle("");
        hist->GetYaxis()->SetTitle("");
    }
    
    // 设置 Y 轴范围 (基于 Yield 的对数刻度)
    // 假设最大 Yield 接近 1e6，最小 Yield 接近 10
    hist->GetYaxis()->SetRangeUser(1, 2e6); 
    
    // 设置刻度分隔
    hist->GetYaxis()->SetNdivisions(505);
}


// 数据拷贝辅助函数：将源直方图在指定刚度范围内的数据拷贝到目标直方图
void copyData(TH1D* target, TH1D* source, double minR, double maxR) {
    if (!source) return;
    for (int i = 1; i <= source->GetNbinsX(); ++i) {
        double rigidity = source->GetBinCenter(i);
        if (rigidity >= minR && rigidity < maxR) {
            int targetBin = target->FindBin(rigidity);
            if (targetBin > 0 && targetBin <= target->GetNbinsX()) {
                target->SetBinContent(targetBin, source->GetBinContent(i));
                target->SetBinError(targetBin, source->GetBinError(i));
            }
        }
    }
}


// --- 主绘图函数 (合并 Yield) ---

void DrawFrac() {
    // --- 1. 参数设置 ---
    const std::string windowElement = "Boron";
    const std::string trackerLayer  = "UnbiasedL1Inner";
    const std::string outputPath    = "/eos/user/z/zixuan/Isotope/ChargeTemp/";
    const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
    const std::vector<std::string> elementsToPlot = {"Boron", "Beryllium", "Carbon"}; // 绘制顺序: 主成分优先
    
    // 设置全局ROOT样式
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    // --- 2. 元素和样式定义 ---
    std::map<std::string, int> colorMap = {
        {"Beryllium", kBlue},
        {"Boron",     kOrange - 3},
        {"Carbon",    kGreen+2}
    };
    std::map<std::string, std::string> labelMap = {
        {"Beryllium", "Be"},
        {"Boron",     "B"},
        {"Carbon",    "C"}
    };
    
    // 用于存储拼接后的直方图
    std::vector<std::unique_ptr<TH1D>> combinedHists;

    // --- 3. 打开输入文件 ---
    TString inputFilePath = TString::Format(
        "root://eoshome-z.cern.ch//eos/user/z/zixuan/Isotope/ChargeTemp/PureL1TempFit_AllL1Q%s.root",
        trackerLayer.c_str()
    );
    std::cout << "Processing file: " << inputFilePath << std::endl;
    std::unique_ptr<TFile> file(TFile::Open(inputFilePath, "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << inputFilePath << std::endl;
        return;
    }

    // --- 4. 循环处理每个元素，生成拼接的 Yield 直方图 ---
    std::vector<double> all_bin_edges;
    for (const auto& element : elementsToPlot) {
        
        std::cout << "\nProcessing element (Yield): " << element << std::endl;

        // 4.1 获取并收集 Bin 边界 (只收集一次所有存在的边界)
        std::map<std::string, TH1D*> hists;
        bool found_any_hist = false;
        
        for (const auto& det : detectors) {
            // *** 修改为 h_yield_in_... ***
            TString histName = TString::Format("h_yield_in_%s_from_%s_%s", windowElement.c_str(), element.c_str(), det.c_str());
            TH1D* h = nullptr;
            file->GetObject(histName, h); 
            hists[det] = h;
            
            if (h) {
                found_any_hist = true;
                TAxis* axis = h->GetXaxis();
                for (int i = 1; i <= axis->GetNbins(); ++i) {
                    all_bin_edges.push_back(axis->GetBinLowEdge(i));
                }
                all_bin_edges.push_back(axis->GetBinUpEdge(axis->GetNbins()));
            } else {
                 std::cerr << "Warning: Missing yield hist '" << histName << "'. Skipping this detector range." << std::endl;
            }
        }

        if (!found_any_hist) {
            std::cerr << "Error: No yield histograms found for element " << element << ". Skipping." << std::endl;
            continue;
        }

        // --- 4.2 创建组合直方图（只需为每个元素创建一次） ---
        // 这一步在 4.3 之前执行，用于创建直方图的容器，但 binning 在循环结束后统一确定。
        // 为了简化，我们先跳过这个创建，将所有数据统一到 Boron 的 binning 上，或者使用一个统一的 binning
    }
    
    // 确保 bin 边界唯一且排序
    std::sort(all_bin_edges.begin(), all_bin_edges.end());
    all_bin_edges.erase(std::unique(all_bin_edges.begin(), all_bin_edges.end()), all_bin_edges.end());
    
    if (all_bin_edges.size() < 2) {
        std::cerr << "Error: Not enough bin boundaries found for combined plot." << std::endl;
        return;
    }


    // --- 5. 再次循环，创建组合图并填充数据 ---
    for (const auto& element : elementsToPlot) {
        
        // 5.1 重新获取数据 (这里可以优化，但为清晰保持结构)
        std::map<std::string, TH1D*> hists;
        bool found_any_hist = false;
        for (const auto& det : detectors) {
            TString histName = TString::Format("h_yield_in_%s_from_%s_%s", windowElement.c_str(), element.c_str(), det.c_str());
            TH1D* h = nullptr;
            file->GetObject(histName, h); 
            hists[det] = h;
            if (h) found_any_hist = true;
        }
        if (!found_any_hist) continue;

        // 5.2 创建新的组合直方图 (使用统一的 binning)
        TString combinedHistName = TString::Format("h_combined_yield_in_%s_from_%s", windowElement.c_str(), element.c_str());
        auto h_combined = std::make_unique<TH1D>(combinedHistName, "", 
                                                 all_bin_edges.size() - 1, &all_bin_edges[0]);

        // 5.3 填充数据
        copyData(h_combined.get(), hists["TOF"], 0.0, 1.1);
        copyData(h_combined.get(), hists["NaF"], 1.1, 3.06);
        copyData(h_combined.get(), hists["AGL"], 3.06, 1e9); 

        // 5.4 设置样式
        int color = colorMap.count(element) ? colorMap[element] : kBlack;
        h_combined->SetLineColor(color);
        h_combined->SetMarkerColor(color);
        h_combined->SetMarkerStyle(20);
        h_combined->SetMarkerSize(1.2);
        
        // 5.5 存储以供绘制
        combinedHists.push_back(std::move(h_combined));
    }


    // --- 6. 创建画布并绘制所有组合图 ---
    TString canvasName = TString::Format("canvas_combined_yield_in_%s", windowElement.c_str());
    std::unique_ptr<TCanvas> canvas(new TCanvas(canvasName, "", 800, 600)); // 调高画布
    setupCanvas(canvas.get());
    
    // 创建图例
    std::unique_ptr<TLegend> legend(new TLegend(0.8, 0.7, 0.95, 0.88));
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    
    bool isFirst = true;
    for (const auto& h_ptr : combinedHists) {
        std::string element = h_ptr->GetName(); // 从名称中提取元素名
        // 提取元素名 (从 "h_combined_yield_in_Boron_from_ELEMENT" 中)
        element = element.substr(element.rfind('_') + 1); 

        // 6.1 设置轴和样式
        setupHistogramForYield(h_ptr.get(), element, isFirst);
        h_ptr->GetXaxis()->SetRangeUser(0.41, 21);
        h_ptr->GetYaxis()->SetRangeUser(10, 1000000);
        
        // 6.2 绘制
        if (isFirst) {
            h_ptr->Draw("p0");
            isFirst = false;
        } else {
            h_ptr->Draw("p0 same"); // 叠加绘制
        }
        
        // 6.3 添加图例
        legend->AddEntry(h_ptr.get(), labelMap.at(element).c_str(), "p");
    }
    
    if (!combinedHists.empty()) {
        legend->Draw();
    }
    
    // --- 7. 保存图像 ---
    TString outputFileName = TString::Format(
        "%s/Combined_Yield_in_%s_%s.png",
        outputPath.c_str(),
        windowElement.c_str(),
        trackerLayer.c_str()
    );
    
    canvas->SaveAs(outputFileName);
    std::cout << "\n-> Saved combined yield plot: " << outputFileName << std::endl;
    std::cout << "\nAll combined yield plots saved successfully." << std::endl;
}