#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLine.h>
#include <TMath.h>

// Configuration
const std::string baseDir = "/eos/user/z/zixuan/Isotope/Add/";
const std::string outputDir = "/eos/user/z/zixuan/Isotope/";

const std::vector<std::string> elements = {"Carbon", "Nitrogen", "Oxygen"};
const std::vector<std::string> detectors = {"TOF", "NaF", "AGL"};
const std::vector<std::string> chargeTypes = {"UnbiasedL1Inner", "L1Inner"};
const std::vector<std::string> mcFiles = {"C12_all.root", "N14_all.root", "O16_all.root"};

// Detector boundaries for stitching
const double boundary1 = 1.1130274; // GetBinLowEdge(7)
const double boundary2 = 2.4349398; // GetBinLowEdge(10)

class RatioAnalyzer {
public:
    RatioAnalyzer();
    ~RatioAnalyzer() = default;
    
    void runAnalysis();
    
private:
    std::map<std::string, TFile*> mcFiles_;
    TFile* issFile_;
    
    void openFiles();
    void closeFiles();
    
    // MC analysis
    TH1F* processMCElement(const std::string& element, const std::string& chargeType);
    TH1F* getMCSource(const std::string& element, const std::string& chargeType);
    TH1F* getMCFragmentation(const std::string& element, const std::string& chargeType);
    
    // ISS analysis  
    TH1F* processISSElement(const std::string& element, const std::string& chargeType);
    TH1F* getISSSource(const std::string& element, const std::string& chargeType);
    TH1F* getISSFragmentation(const std::string& element, const std::string& chargeType);
    
    // Utility functions
    TH1F* stitchDetectors(TH1F* tof, TH1F* naf, TH1F* agl, const std::string& name);
    TH1F* calculateRatio(TH1F* numerator, TH1F* denominator, const std::string& name);
    void createComparisonPlot(const std::string& element, const std::string& chargeType, 
                            TH1F* mcRatio, TH1F* issRatio);
    void setupHistogramStyle(TH1F* hist, int color, int style = 1);
    std::pair<double, double> calculateYRange(TH1F* h1, TH1F* h2, double xmin, double xmax);
};

RatioAnalyzer::RatioAnalyzer() : issFile_(nullptr) {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
}

void RatioAnalyzer::openFiles() {
    // Open MC files
    for (size_t i = 0; i < elements.size(); ++i) {
        std::string filepath = baseDir + mcFiles[i];
        mcFiles_[elements[i]] = TFile::Open(filepath.c_str());
        if (!mcFiles_[elements[i]] || mcFiles_[elements[i]]->IsZombie()) {
            std::cerr << "Error: Cannot open MC file " << filepath << std::endl;
            return;
        }
    }
    
    // Open ISS file
    std::string issPath = baseDir + "Be_all.root";
    issFile_ = TFile::Open(issPath.c_str());
    if (!issFile_ || issFile_->IsZombie()) {
        std::cerr << "Error: Cannot open ISS file " << issPath << std::endl;
        return;
    }
}

void RatioAnalyzer::closeFiles() {
    for (auto& [key, file] : mcFiles_) {
        if (file) file->Close();
    }
    if (issFile_) issFile_->Close();
}

TH1F* RatioAnalyzer::getMCSource(const std::string& element, const std::string& chargeType) {
    TFile* file = mcFiles_[element];
    
    std::vector<TH1F*> detHists;
    for (const auto& det : detectors) {
        std::string histName = chargeType + "_MC_BKG_H1_" + det;
        TH1F* hist = (TH1F*)file->Get(histName.c_str());
        if (!hist) {
            std::cerr << "Cannot find histogram: " << histName << std::endl;
            return nullptr;
        }
        hist->Rebin(2);
        detHists.push_back(hist);
    }
    
    return stitchDetectors(detHists[0], detHists[1], detHists[2], 
                          element + "_" + chargeType + "_MC_Source");
}

TH1F* RatioAnalyzer::getMCFragmentation(const std::string& element, const std::string& chargeType) {
    TFile* file = mcFiles_[element];
    
    std::vector<TH1F*> totalFragHists;
    
    for (const auto& det : detectors) {
        TH1F* detTotal = nullptr;
        
        // Sum all three Be isotopes (Mass7, Mass9, Mass10)
        for (int mass : {7, 9, 10}) {
            std::string histName = chargeType + "_MC_BKG_H2_" + det + "_Z4_Mass" + std::to_string(mass);
            TH1F* hist = (TH1F*)file->Get(histName.c_str());
            if (!hist) {
                std::cerr << "Cannot find histogram: " << histName << std::endl;
                continue;
            }
            hist->Rebin(2);
            
            if (!detTotal) {
                detTotal = (TH1F*)hist->Clone((det + "_total").c_str());
            } else {
                detTotal->Add(hist);
            }
        }
        
        if (detTotal) {
            totalFragHists.push_back(detTotal);
        }
    }
    
    if (totalFragHists.size() != 3) {
        std::cerr << "Error: Missing fragmentation histograms for " << element << std::endl;
        return nullptr;
    }
    
    return stitchDetectors(totalFragHists[0], totalFragHists[1], totalFragHists[2],
                          element + "_" + chargeType + "_MC_Frag");
}

TH1F* RatioAnalyzer::processMCElement(const std::string& element, const std::string& chargeType) {
    TH1F* source = getMCSource(element, chargeType);
    TH1F* fragmentation = getMCFragmentation(element, chargeType);
    
    if (!source || !fragmentation) {
        std::cerr << "Error: Failed to get MC histograms for " << element << std::endl;
        return nullptr;
    }
    
    return calculateRatio(fragmentation, source, element + "_" + chargeType + "_MC_Ratio");
}

TH1F* RatioAnalyzer::getISSSource(const std::string& element, const std::string& chargeType) {
    std::vector<TH1F*> detHists;
    
    for (const auto& det : detectors) {
        std::string histName = chargeType + "_ISS_BKG_H1_" + element + "_" + det;
        TH1F* hist = (TH1F*)issFile_->Get(histName.c_str());
        if (!hist) {
            std::cerr << "Cannot find histogram: " << histName << std::endl;
            return nullptr;
        }
        hist->Rebin(2);
        detHists.push_back(hist);
    }
    
    return stitchDetectors(detHists[0], detHists[1], detHists[2],
                          element + "_" + chargeType + "_ISS_Source");
}

TH1F* RatioAnalyzer::getISSFragmentation(const std::string& element, const std::string& chargeType) {
    std::vector<TH1F*> detHists;
    
    for (const auto& det : detectors) {
        std::string histName = chargeType + "_ISS_BKG_H3_" + element + "_" + det;
        TH1F* hist = (TH1F*)issFile_->Get(histName.c_str());
        if (!hist) {
            std::cerr << "Cannot find histogram: " << histName << std::endl;
            return nullptr;
        }
        hist->Rebin(2);
        detHists.push_back(hist);
    }
    
    return stitchDetectors(detHists[0], detHists[1], detHists[2],
                          element + "_" + chargeType + "_ISS_Frag");
}

TH1F* RatioAnalyzer::processISSElement(const std::string& element, const std::string& chargeType) {
    TH1F* source = getISSSource(element, chargeType);
    TH1F* fragmentation = getISSFragmentation(element, chargeType);
    
    if (!source || !fragmentation) {
        std::cerr << "Error: Failed to get ISS histograms for " << element << std::endl;
        return nullptr;
    }
    
    std::cout << "  " << element << " " << chargeType << " - Source entries: " << source->GetEntries() 
              << ", Fragmentation entries: " << fragmentation->GetEntries() << std::endl;
    
    return calculateRatio(fragmentation, source, element + "_" + chargeType + "_ISS_Ratio");
}

TH1F* RatioAnalyzer::stitchDetectors(TH1F* tof, TH1F* naf, TH1F* agl, const std::string& name) {
    if (!tof || !naf || !agl) return nullptr;
    
    // Create combined histogram with same binning as input
    TH1F* combined = (TH1F*)tof->Clone(name.c_str());
    combined->Reset();
    
    int nBins = tof->GetNbinsX();
    
    for (int i = 1; i <= nBins; ++i) {
        double binCenter = tof->GetBinCenter(i);
        double content = 0, error2 = 0;
        
        if (binCenter < boundary1) {
            // Use TOF
            content = tof->GetBinContent(i);
            error2 = pow(tof->GetBinError(i), 2);
        } else if (binCenter < boundary2) {
            // Use NaF
            content = naf->GetBinContent(i);
            error2 = pow(naf->GetBinError(i), 2);
        } else {
            // Use AGL
            content = agl->GetBinContent(i);
            error2 = pow(agl->GetBinError(i), 2);
        }
        
        combined->SetBinContent(i, content);
        combined->SetBinError(i, sqrt(error2));
    }
    
    return combined;
}

TH1F* RatioAnalyzer::calculateRatio(TH1F* numerator, TH1F* denominator, const std::string& name) {
    if (!numerator || !denominator) return nullptr;
    
    TH1F* ratio = (TH1F*)numerator->Clone(name.c_str());
    
    // Simple bin-by-bin division without any corrections
    for (int i = 1; i <= ratio->GetNbinsX(); ++i) {
        double num = numerator->GetBinContent(i);
        double den = denominator->GetBinContent(i);
        double num_err = numerator->GetBinError(i);
        double den_err = denominator->GetBinError(i);
        
        if (den > 0) {
            double r = num / den;
            // Error propagation for ratio: sqrt((σ_num/den)² + (num*σ_den/den²)²)
            double err = 0;
            if (num > 0) {
                err = r * sqrt(pow(num_err/num, 2) + pow(den_err/den, 2));
            } else {
                err = num_err / den;
            }
            
            ratio->SetBinContent(i, r);
            ratio->SetBinError(i, err);
        } else {
            ratio->SetBinContent(i, 0);
            ratio->SetBinError(i, 0);
        }
    }
    
    return ratio;
}

void RatioAnalyzer::setupHistogramStyle(TH1F* hist, int color, int style) {
    if (!hist) return;
    
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.0);
    hist->SetLineWidth(2);
    hist->SetLineStyle(style);
    
    // Keep X title, remove Y title and other labels
    hist->SetTitle("");
    hist->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
    hist->GetYaxis()->SetTitle("");
    hist->GetXaxis()->SetTitleSize(0.12);
    hist->GetXaxis()->SetLabelSize(0.12);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
}

std::pair<double, double> RatioAnalyzer::calculateYRange(TH1F* h1, TH1F* h2, double xmin, double xmax) {
    double ymin = 1e10, ymax = -1e10;
    
    // Find min and max in the specified x range
    for (int i = 1; i <= h1->GetNbinsX(); ++i) {
        double x = h1->GetBinCenter(i);
        if (x < xmin || x > xmax) continue;
        
        double y1 = h1->GetBinContent(i);
        double y2 = h2->GetBinContent(i);
        
        if (y1 > 0) {
            ymin = TMath::Min(ymin, y1);
            ymax = TMath::Max(ymax, y1);
        }
        if (y2 > 0) {
            ymin = TMath::Min(ymin, y2);
            ymax = TMath::Max(ymax, y2);
        }
    }
    
    // Add some margin
    if (ymin < 1e9 && ymax > -1e9) {
        ymin *= 0.5;  // Lower by factor of 2
        ymax *= 2.0;  // Upper by factor of 2
    } else {
        ymin = 1e-4;
        ymax = 1e-1;
    }
    
    return {0, 0.6*ymax};
}

void RatioAnalyzer::createComparisonPlot(const std::string& element, const std::string& chargeType,
                                       TH1F* mcRatio, TH1F* issRatio) {
    if (!mcRatio || !issRatio) {
        std::cerr << "Error: Missing histograms for " << element << " " << chargeType << std::endl;
        return;
    }
    
    // Setup styles
    setupHistogramStyle(mcRatio, kRed, 1);
    setupHistogramStyle(issRatio, kBlack, 1);
    
    // Create canvas
    TCanvas* canvas = new TCanvas(("c_" + element + "_" + chargeType).c_str(),
                                 "", 800, 800);
    canvas->Divide(1, 2);
    
    // Calculate Y range for ratio plot based on data in x range [0.5, 16.3]
    auto [ratioYMin, ratioYMax] = calculateYRange(mcRatio, issRatio, 0.5, 16.3);
    
    // Upper pad - ratio comparison
    TPad* pad1 = (TPad*)canvas->cd(1);
    pad1->SetLogx();
    pad1->SetPad(0, 0.3, 1, 1);
    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.02);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.02);
    pad1->SetGridy();
    pad1->SetGridx();
    
    mcRatio->GetYaxis()->SetRangeUser(ratioYMin, ratioYMax);
    mcRatio->GetXaxis()->SetRangeUser(0.5, 16.3);
    mcRatio->GetXaxis()->SetLabelSize(0);
    
    mcRatio->Draw("EP");
    issRatio->Draw("EP SAME");
    
    std::cout << "  Ratio Y range: [" << ratioYMin << ", " << ratioYMax << "]" << std::endl;
    
    // Create ISS/MC ratio histogram
    TH1F* dataToMC = (TH1F*)issRatio->Clone("dataToMC");
    dataToMC->Divide(mcRatio);
    
    // Calculate Y range for ISS/MC ratio in x range [0.5, 16.3]
    double issmc_ymin = 1e10, issmc_ymax = -1e10;
    for (int i = 1; i <= dataToMC->GetNbinsX(); ++i) {
        double x = dataToMC->GetBinCenter(i);
        if (x < 0.5 || x > 16.3) continue;
        
        double y = dataToMC->GetBinContent(i);
        if (y > 0 && y < 1e9) {  // Exclude unreasonable values
            issmc_ymin = TMath::Min(issmc_ymin, y);
            issmc_ymax = TMath::Max(issmc_ymax, y);
        }
    }
    
    // Add margin to ISS/MC range
    if (issmc_ymin < 1e9 && issmc_ymax > -1e9) {
        double center = (issmc_ymin + issmc_ymax) / 2;
        double range = issmc_ymax - issmc_ymin;
        issmc_ymin = center - range * 0.6;  // 20% extra margin
        issmc_ymax = center + range * 0.6;
        
        // Ensure reasonable bounds
        issmc_ymin = TMath::Max(issmc_ymin, 0.1);
        issmc_ymax = TMath::Min(issmc_ymax, 10.0);
    } else {
        issmc_ymin = 0.5;
        issmc_ymax = 1.5;
    }
    
    // Lower pad - ISS/MC ratio
    TPad* pad2 = (TPad*)canvas->cd(2);
    pad2->SetLogx();
    pad2->SetPad(0, 0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.02);
    pad2->SetGridy();
    
    dataToMC->SetTitle("");
    dataToMC->GetYaxis()->SetTitle("");
    dataToMC->GetXaxis()->SetTitle("E_{k}/n [GeV/n]");
    dataToMC->GetYaxis()->SetRangeUser(0.2,1.5);
    dataToMC->GetXaxis()->SetRangeUser(0.5, 16.3);
    dataToMC->GetYaxis()->SetLabelSize(0.08);
    dataToMC->GetXaxis()->SetLabelSize(0.08);
    dataToMC->GetYaxis()->SetTitleSize(0);
    dataToMC->GetXaxis()->SetTitleSize(0.1);
    dataToMC->GetXaxis()->SetTitleOffset(1.2);
    dataToMC->GetYaxis()->SetNdivisions(505);
    
    setupHistogramStyle(dataToMC, kBlack, 1);
    dataToMC->SetMarkerSize(0.6);
    dataToMC->Draw("EP");
    
    std::cout << "  ISS/MC Y range: [" << issmc_ymin << ", " << issmc_ymax << "]" << std::endl;
    
    // Add unity line
    TLine* unityLine = new TLine(0.5, 1, 16.3, 1);
    unityLine->SetLineStyle(2);
    unityLine->SetLineColor(kGray + 1);
    unityLine->SetLineWidth(2);
    unityLine->Draw("SAME");
    
    // Save plot
    std::string outputName = outputDir + "FragmentationRatio_Clean_" + element + "_" + chargeType + ".png";
    canvas->SaveAs(outputName.c_str());
    
    outputName = outputDir + "FragmentationRatio_Clean_" + element + "_" + chargeType + ".pdf";
    canvas->SaveAs(outputName.c_str());
    
    std::cout << "Saved clean plot: " << outputName << std::endl;
    
    delete canvas;
    delete unityLine;
}

void RatioAnalyzer::runAnalysis() {
    std::cout << "Starting clean fragmentation ratio analysis..." << std::endl;
    
    openFiles();
    
    for (const auto& chargeType : chargeTypes) {
        std::cout << "\nProcessing " << chargeType << " acceptance..." << std::endl;
        
        for (const auto& element : elements) {
            std::cout << "  Processing " << element << "..." << std::endl;
            
            // Process MC
            TH1F* mcRatio = processMCElement(element, chargeType);
            if (!mcRatio) {
                std::cerr << "    Failed to process MC for " << element << std::endl;
                continue;
            }
            
            // Process ISS
            TH1F* issRatio = processISSElement(element, chargeType);
            if (!issRatio) {
                std::cerr << "    Failed to process ISS for " << element << std::endl;
                delete mcRatio;
                continue;
            }
            
            // Create comparison plot
            createComparisonPlot(element, chargeType, mcRatio, issRatio);
            
            // Cleanup
            delete mcRatio;
            delete issRatio;
            
            std::cout << "    Completed " << element << std::endl;
        }
    }
    
    closeFiles();
    std::cout << "\nClean fragmentation ratio analysis completed!" << std::endl;
}

// Entry point
void BkgRatioCp() {
    RatioAnalyzer analyzer;
    analyzer.runAnalysis();
}