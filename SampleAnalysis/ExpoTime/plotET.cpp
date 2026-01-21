#include "../Tool.h"

using namespace AMS_Iso;

struct IsotopeConfig {
    std::string name;
    int charge;
    std::vector<int> masses;
    std::vector<std::string> legends;
    std::vector<int> colors;

    IsotopeConfig(const std::string& elementName) {
        // 从IsotopeData查找对应元素
        for (const auto& iso : IsotopeData) {
            if (ConvertElementName(iso.name_, false) == elementName) {
                name = elementName;
                charge = iso.charge_;
                // 只添加Mass7
                for (int i = 0; i < iso.isotope_count_; ++i) {
                    if (iso.mass_[i] == 7) {  // 只处理Mass7
                        masses.push_back(iso.mass_[i]);
                        legends.push_back(elementName + std::to_string(iso.mass_[i]));
                        colors.push_back(iso.color_[i]);
                        break;  // 找到Mass7就退出
                    }
                }
                break;
            }
        }
    }
};

void plotIsoET(const std::string& element, const std::string& inputFile) {
    // 获取同位素配置
    IsotopeConfig config(element);
    if (config.masses.empty()) {
        std::cerr << "Error: Cannot find Mass7 for element " << element << std::endl;
        return;
    }

    const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
    const int nMasses = 1;  // 只有Mass7
    const int nDetectors = detectors.size();

    // 打开输入文件
    TFile* inputFilePtr = TFile::Open(inputFile.c_str());
    if (!inputFilePtr || inputFilePtr->IsZombie()) {
        std::cerr << "Error: Unable to open file " << inputFile << std::endl;
        return;
    }

    std::cout << "Successfully opened file: " << inputFile << std::endl;
    
    // 创建输出文件
    std::string outputRootFile = "/eos/user/z/zixuan/Isotope/ExpoTime/" + 
                                 element + "_Mass7_ExposureTime.root";
    TFile* outputFile = new TFile(outputRootFile.c_str(), "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Unable to create output file " << outputRootFile << std::endl;
        inputFilePtr->Close();
        return;
    }

    std::cout << "Created output file: " << outputRootFile << std::endl;

    // 创建存储直方图的二维数组
    std::vector<std::vector<TH1F*>> histograms(nMasses);
    for (int i = 0; i < nMasses; ++i) {
        histograms[i].resize(nDetectors);
        for (int j = 0; j < nDetectors; ++j) {
            std::string histName = Form("ISS_FLUX_H3_%s_Mass%dBin", 
                                       detectors[j].c_str(), config.masses[i]);
            
            std::cout << "Looking for histogram: " << histName << std::endl;
            
            TH1F* tempHist = (TH1F*)inputFilePtr->Get(histName.c_str());
            
            if (!tempHist) {
                std::cerr << "Error: Cannot find histogram " << histName 
                          << " in file " << inputFile << std::endl;
                continue;
            }
            
            // 克隆直方图并保存到输出文件
            std::string cloneName = Form("%s_Mass7_%s", element.c_str(), detectors[j].c_str());
            histograms[i][j] = (TH1F*)tempHist->Clone(cloneName.c_str());
            
            if (histograms[i][j]) {
                histograms[i][j]->SetTitle(Form("%s Mass7 %s Exposure Time;E_{k}/n [GeV/n];Exposure Time [s]", 
                                               element.c_str(), detectors[j].c_str()));
                histograms[i][j]->SetLineColor(config.colors[i]);
                histograms[i][j]->SetLineWidth(2);
                histograms[i][j]->SetStats(0);
                
                // 写入输出文件
                outputFile->cd();
                histograms[i][j]->Write();
                
                std::cout << "Saved histogram: " << cloneName << std::endl;
            }
        }
    }

    // 创建画布
    TCanvas *c1 = new TCanvas("c1", "Exposure Time for Mass7", 1200, 400);
    c1->SetLeftMargin(0.10);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.12);
    c1->SetBottomMargin(0.15);
    c1->Divide(3, 1);

    // 从BetaTypes获取能量范围
    std::vector<std::pair<double, double>> xRanges;
    for (const auto& beta : Detector::BetaTypes) {
        xRanges.push_back({std::max(beta.Ekn_range_[0],0.5), beta.Ekn_range_[1]});
    }

    // 绘制直方图
    for (int det = 0; det < nDetectors; ++det) {
        c1->cd(det + 1);
        
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.03);
        gPad->SetLogx();
        
        // 获取当前探测器的直方图
        TH1F* hist = nullptr;
        if (nMasses > 0 && histograms[0].size() > det) {
            hist = histograms[0][det];
        }
        
        if (!hist) {
            std::cerr << "Warning: No histogram for detector " << detectors[det] << std::endl;
            continue;
        }

        // 设置坐标轴
        hist->GetXaxis()->SetRangeUser(det == 0 ? 0.45 : xRanges[det].first, xRanges[det].second);
        hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.2);
        
        hist->GetXaxis()->SetTitleOffset(1.0);
        hist->GetYaxis()->SetTitleOffset(1.4);
        hist->GetXaxis()->SetLabelSize(0.07);
        hist->GetYaxis()->SetLabelSize(0.07);
        hist->GetXaxis()->SetTitleSize(0.08);
        hist->GetYaxis()->SetTitleSize(0.08);
        
        hist->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
        hist->GetYaxis()->SetTitle("Exposure Time [s]");
        
        // 绘制
        hist->Draw("hist");
        
        // 添加探测器标签
        TLatex latex;
        latex.SetTextSize(0.09);
        latex.SetTextAlign(22);
        latex.DrawLatexNDC(0.5, 0.92, detectors[det].c_str());
        
        // 添加同位素标签
        latex.SetTextSize(0.07);
        latex.SetTextAlign(12);
        latex.DrawLatexNDC(0.18, 0.82, Form("%s7", element.c_str()));
    }

    // 添加总标题
    c1->cd(0);
    TLatex title;
    title.SetTextSize(0.03);
    title.SetTextAlign(22);
    title.DrawLatexNDC(0.5, 0.98, Form("%s Mass7 Exposure Time Comparison", element.c_str()));

    // 保存图像
    std::string outputPdf = "/eos/user/z/zixuan/Isotope/ExpoTime/" + 
                           element + "_Mass7_ExposureTime.pdf";
    c1->Print(outputPdf.c_str());
    std::cout << "Saved PDF: " << outputPdf << std::endl;

    // 保存画布到输出文件
    outputFile->cd();
    c1->Write("Canvas_Mass7");
    
    // 关闭并保存输出文件
    outputFile->Close();
    delete outputFile;
    
    // 关闭输入文件
    inputFilePtr->Close();
    delete inputFilePtr;

    // 清理内存
    delete c1;
    for (auto& histVec : histograms) {
        for (auto* hist : histVec) {
            if (hist) delete hist;
        }
    }

    std::cout << "All done! Results saved to:" << std::endl;
    std::cout << "  PDF: " << outputPdf << std::endl;
    std::cout << "  ROOT: " << outputRootFile << std::endl;
}

void plotET()
{
    // 处理Be的Mass7
    std::string beFile = "/eos/user/z/zixuan/Isotope/Add/Be_frag4_withBkg_Tune_full.root";
    
    // 检查文件是否存在
    if (gSystem->AccessPathName(beFile.c_str())) {
        std::cerr << "Error: File not found: " << beFile << std::endl;
        return;
    }
    
    std::cout << "Processing Be Mass7..." << std::endl;
    plotIsoET("Be", beFile);
    
    // 如果需要处理B的Mass10（示例）
    // std::string bFile = "/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Bor_BeToC.root";
    // plotIsoET("B", bFile);
}