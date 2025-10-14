#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <map>
#include <string>
#include <numeric>

// --- Helper: ValueWithError ---
struct ValueWithError {
    double val = 0.0;
    double err = 0.0;

    ValueWithError() = default;
    ValueWithError(double v, double e) : val(v), err(e) {}

    ValueWithError operator+(const ValueWithError& o) const {
        return {val + o.val, std::hypot(err, o.err)};
    }
    ValueWithError operator-(const ValueWithError& o) const {
        return {val - o.val, std::hypot(err, o.err)};
    }
    ValueWithError operator*(double scalar) const {
        return {val * scalar, std::abs(err * scalar)};
    }
    // Standard error propagation for Z = X*Y, assuming X and Y are uncorrelated:
    // err_Z^2 = (Y*err_X)^2 + (X*err_Y)^2
    ValueWithError operator*(const ValueWithError& o) const {
        double p_val = val * o.val;
        double p_err = std::hypot(err * o.val, val * o.err);
        return {p_val, p_err};
    }
    // Standard error propagation for Z = X/Y, assuming X and Y are uncorrelated:
    // err_Z^2 = (err_X/Y)^2 + (X*err_Y/Y^2)^2
    ValueWithError operator/(const ValueWithError& o) const {
        if (o.val == 0.0) {
            throw std::runtime_error("Division by zero in ValueWithError (denominator is zero)");
        }
        double q_val = val / o.val;
        double term1_err_div_y = err / o.val;
        double term2_x_erry_div_y2 = (val * o.err) / (o.val * o.val);
        double q_err = std::hypot(term1_err_div_y, term2_x_erry_div_y2);
        return {q_val, q_err};
    }
};
ValueWithError operator*(double scalar, const ValueWithError& ve) {
    return ve * scalar;
}

// --- Helper: File Operations ---
std::unique_ptr<TFile> openFile(const char* filename, const char* option = "READ") {
    auto file = std::unique_ptr<TFile>(TFile::Open(filename, option));
    if (!file || file->IsZombie()) {
        throw std::runtime_error("Error opening file: " + std::string(filename) + " with option " + std::string(option));
    }
    return file;
}

TH1F* getHistFromFile(TFile* file, const TString& histName) {
    if (!file) {
         throw std::runtime_error("Input TFile pointer is null when trying to get histogram: " + std::string(histName.Data()));
    }
    TH1F* hist = nullptr;
    file->GetObject(histName, hist);
    if (!hist) {
        throw std::runtime_error("Could not find histogram: " + std::string(histName.Data()) + " in file " + file->GetName());
    }
    return hist;
}

TH1D* getHistFromFileD(TFile* file, const TString& histName) {
    if (!file) {
        throw std::runtime_error("Input TFile pointer is null when trying to get histogram: " + std::string(histName.Data()));
    }
    TH1D* hist = nullptr;  // T 被替换为 TH1D
    file->GetObject(histName, hist);
    if (!hist) {
        throw std::runtime_error("Could not find histogram: " + std::string(histName.Data()) + " in file " + file->GetName());
    }
    return hist;
}
// --- Helper: Histogram Creation & Styling ---
TH1F* createTH1F(const char* name, const char* title, int nBins, const double* bins, std::vector<TObject*>& objList) {
    auto hist = new TH1F(name, title, nBins, bins);
    objList.push_back(hist);
    return hist;
}

void styleHist(TH1F* hist, EColor color, EMarkerStyle markerStyle = (EMarkerStyle)20) {
    if (!hist) return;
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(0.8);
    hist->SetStats(0);
}

void fillHistBin(TH1F* hist, int bin_idx, const ValueWithError& vw) {
    if(hist) {
        hist->SetBinContent(bin_idx, vw.val);
        hist->SetBinError(bin_idx, vw.err);
    }
}

// --- Detector Logic ---
enum class Detector { TOF, NaF, AGL };

Detector getDetectorForEnergy(double energy) {
    if (energy < 1.17) return Detector::TOF;
    if (energy < 3.23) return Detector::NaF;
    return Detector::AGL;
}

std::string getDetectorName(Detector detector) {
    switch (detector) {
        case Detector::TOF: return "TOF";
        case Detector::NaF: return "NaF";
        case Detector::AGL: return "AGL";
        default: throw std::runtime_error("Unknown detector type");
    }
}

// --- Main Calculation Function ---
void calculateFragmentedBeryllium(bool isISS) {
    const std::array<double, 15> WideBins {{
        0.61, 0.86, 1.17, 1.55, 2.01, 2.57, 3.23, 4.00, 4.91, 5.99, 7.18, 8.60, 10.25, 12.13, 16.38
    }};
    const std::array<double, 31> NarrowBins {{
        0.42, 0.51, 0.61, 0.73, 0.86, 1.00, 1.17, 1.35, 1.55, 1.77, 2.01, 2.28, 2.57, 2.88, 3.23, 3.60, 4.00, 4.44, 4.91, 5.42, 5.99, 6.56, 7.18, 7.86, 8.60, 9.40, 10.25, 11.16, 12.13, 14.35, 16.38
    }};
    const int nBins = WideBins.size() - 1;
    const std::array<std::string, 3> isotopes = {"Be7", "Be9", "Be10"};
    const std::array<std::string, 2> primary_isotopes = {"Be7", "Be9"};

    std::vector<TObject*> allObjectsToSave;
    
    //bkgBBeFinal
    const char *FileSuffix = isISS ? "BtoBe_SDIATBeBcut" : "MCB11toBeIso";
    cout<<FileSuffix<<endl;
    const char *FilePrefix = isISS ? "ISS" : "MCB11";
    const char *CountsFileName = isISS ? "Bor_temp_narrow_bkg_unbL1_SDIATBeBcut" : "B11_temp_wide_bkg_frag";

    auto fragmentBeFitFile = openFile(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/MassTF_Be_wide_%s_Use7.root", FileSuffix));
    auto normalBeFitFile   = openFile("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/TempFit/MassTF_Be_wide_normalcut_Use7.root");
    auto CountsFile        = openFile(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/%s.root", CountsFileName));
    auto Be_ratioFile      = openFile("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGE_detfuncT_unb_wide.root");
    auto B_ratioFile       = openFile("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/chargeTempFit/chargeFit_B_EGE_detfuncT_unb_wide.root");

    auto B11IsoFracFile = openFile("/eos/ams/user/z/zetong/Boron_isotope/pre/B11f1.root"); 
    
    // 新增：打开B10 MC文件
    auto B10MCFile = openFile("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/B10_temp_narrow_frag_unbL1_SDIATBeBcut.root");
    
    auto outputFile = openFile(Form("/eos/ams/user/z/zuhao/yanzx/Isotope/IsoResults/IssValidateL12/%s_FragmentBer_MassFitAnalysis_unb_perfectQfit.root",FilePrefix), "RECREATE");

    TH1D* h_be_ratio_src    = getHistFromFileD(Be_ratioFile.get(), "h_be_fraction");
    h_be_ratio_src->Rebin(2); 
    TH1D* h_boron_ratio_src = getHistFromFileD(B_ratioFile.get(), "h_b_fraction");
    h_boron_ratio_src->Rebin(2); 
    
    TH1D* h_B11_Frac = (TH1D*)B11IsoFracFile->Get("B11");
    
    struct OutputHistograms {
        TH1F* fragmented_counts = nullptr;
        TH1F* fragmented_fraction = nullptr;
        TH1F* normal_fraction = nullptr;
        TH1F* realfrac_totalfrac_ratio = nullptr;
        TH1F* raw_isotope_counts = nullptr;
        TH1F* BeIso_B11Source_Ratio = nullptr;
        TH1F* BeIso_B10Source_Ratio = nullptr;
        TH1F* BeIso_BorSource_Ratio = nullptr;
    };
    std::map<std::string, OutputHistograms> outHists;

    for (const auto& iso : isotopes) {
        std::string iso_num_str = (iso.length() > 2) ? iso.substr(2) : "?";
        outHists[iso].fragmented_counts = createTH1F(Form("h_%s_fragmented", iso.c_str()), Form("^{%s}Be from Boron fragmentation;E_{k}/n [GeV/n];Counts", iso_num_str.c_str()), nBins, WideBins.data(), allObjectsToSave);
        outHists[iso].fragmented_fraction = createTH1F(Form("h_fragment_%s_frac", iso.c_str()), Form("Fragmented ^{%s}Be Fraction;E_{k}/n [GeV/n];Fraction", iso_num_str.c_str()), nBins, WideBins.data(), allObjectsToSave);
        outHists[iso].normal_fraction = createTH1F(Form("h_normal_%s_frac", iso.c_str()), Form("Normal ^{%s}Be Fraction;E_{k}/n [GeV/n];Fraction", iso_num_str.c_str()), nBins, WideBins.data(), allObjectsToSave);
        outHists[iso].realfrac_totalfrac_ratio = createTH1F(Form("realfrac_sel_%s_ratio", iso.c_str()), Form("^{%s}Be exclude contamination / selectde ^{%s}Be ;E_{k}/n [GeV/n];Ratio", iso_num_str.c_str(), iso_num_str.c_str()), nBins, WideBins.data(), allObjectsToSave);
        outHists[iso].raw_isotope_counts = createTH1F(Form("h_raw_%s_counts", iso.c_str()), Form("Raw ^{%s}Be Counts;E_{k}/n [GeV/n];Counts", iso_num_str.c_str()), nBins, WideBins.data(), allObjectsToSave);
        //final result: fragmented Be Isotope(exclude contami.) / Boron Isotope in L1
        outHists[iso].BeIso_B11Source_Ratio = createTH1F(Form("BeIso_B11Source_Ratio_%s", iso.c_str()), Form("Fragmented ^{%s}Be(exclude contami.) / Boron11 in L1;E_{k}/n [GeV/n];Ratio", iso_num_str.c_str()), nBins, WideBins.data(), allObjectsToSave);
        outHists[iso].BeIso_B10Source_Ratio = createTH1F(Form("BeIso_B10Source_Ratio_%s", iso.c_str()), Form("Fragmented ^{%s}Be(exclude contami.) / Boron10 in L1;E_{k}/n [GeV/n];Ratio", iso_num_str.c_str()), nBins, WideBins.data(), allObjectsToSave);
        outHists[iso].BeIso_BorSource_Ratio = createTH1F(Form("BeIso_BorSource_Ratio_%s", iso.c_str()), Form("Fragmented ^{%s}Be(exclude contami.) / Boron in L1;E_{k}/n [GeV/n];Ratio", iso_num_str.c_str()), nBins, WideBins.data(), allObjectsToSave);
    }
    
    TH1F* h_L1SourceCounts     = createTH1F("h_L1SourceCounts", "L1 Source Counts;E_{k}/n [GeV/n];Counts", nBins, WideBins.data(), allObjectsToSave);
    TH1F* h_contamination_Be = createTH1F("h_contamination_Be", "Contamination Be Counts;E_{k}/n [GeV/n];Counts", nBins, WideBins.data(), allObjectsToSave);
    TH1F* h_boron_counts_out = createTH1F("h_boron_counts_out", "L1 Source Boron(Calculate);E_{k}/n [GeV/n];Counts", nBins, WideBins.data(), allObjectsToSave);
    TH1F* h_fragment_Be_raw  = createTH1F("h_fragment_Be_raw", "Raw Fragmented Be (Total);E_{k}/n [GeV/n];Counts", nBins, WideBins.data(), allObjectsToSave);
    
    // 新增：用于存储B11碎裂Be10与B11碎裂Be总数的比值
    TH1F* h_Be10_from_B11_ratio = createTH1F("h_Be10_from_B11_ratio", "Be10 from B11 / Total Be from B11;E_{k}/n [GeV/n];Ratio", nBins, WideBins.data(), allObjectsToSave);

    std::map<std::string, TH1F*> inputFragFrac_HistMap; 
    std::map<std::string, TH1F*> inputNormFrac_HistMap;
    std::map<std::string, TH1F*> inputL1SourceCounts_HistMap;
    std::map<std::string, TH1F*> inputFragNucCounts_HistMap; 
    
    // 新增：B10 MC数据的直方图
    std::map<std::string, TH1F*> B10MC_Be7_HistMap;
    std::map<std::string, TH1F*> B10MC_Be9_HistMap;
    std::map<std::string, TH1F*> B10MC_Source_HistMap;

    const char *L1SourcePrefix = isISS ? "Cut" : "Truth"; 
    for (const auto* det_cstr : {"TOF", "NaF", "AGL"}) {
        std::string det_str(det_cstr);
        for (const auto& iso : primary_isotopes) { 
            inputFragFrac_HistMap[iso + "_" + det_str] = getHistFromFile(fragmentBeFitFile.get(), Form("h_best_%s_frac_%s", iso.c_str(), det_cstr));
            inputNormFrac_HistMap[iso + "_" + det_str] = getHistFromFile(normalBeFitFile.get(), Form("h_best_%s_frac_%s", iso.c_str(), det_cstr));
        }
        inputL1SourceCounts_HistMap[det_str] = getHistFromFile(CountsFile.get(), Form("%sL1SourceCounts_%s", L1SourcePrefix, det_cstr));
        inputL1SourceCounts_HistMap[det_str]->Rebin(2);
        inputFragNucCounts_HistMap[det_str] = getHistFromFile(CountsFile.get(), Form("FragNucCounts_%s", det_cstr));
        inputFragNucCounts_HistMap[det_str]->Rebin(2);
        
        // 新增：读取B10 MC数据
        B10MC_Be7_HistMap[det_str] = getHistFromFile(B10MCFile.get(), Form("SourceToIsoCounts_%s_Be7", det_cstr));
        B10MC_Be9_HistMap[det_str] = getHistFromFile(B10MCFile.get(), Form("SourceToIsoCounts_%s_Be9", det_cstr));
        B10MC_Source_HistMap[det_str] = getHistFromFile(B10MCFile.get(), Form("TruthL1SourceCounts_%s", det_cstr));
        // Rebin以匹配主数据
        B10MC_Be7_HistMap[det_str]->Rebin(2);
        B10MC_Be9_HistMap[det_str]->Rebin(2);
        B10MC_Source_HistMap[det_str]->Rebin(2);
    }

    for (int bin_idx = 1; bin_idx <= nBins; ++bin_idx) {
        double binCenter = outHists["Be7"].fragmented_counts->GetBinCenter(bin_idx);
        std::string detName = getDetectorName(getDetectorForEnergy(binCenter));
        cout<<detName<<endl;

        std::map<std::string, ValueWithError> frag_fracs_vw, norm_fracs_vw, raw_iso_counts_vw, pure_iso_counts_vw;
        std::map<std::string, ValueWithError> pure_iso_ratios_vw_map;

        for (const auto& iso : primary_isotopes) {
            frag_fracs_vw[iso] = {inputFragFrac_HistMap[iso + "_" + detName]->GetBinContent(bin_idx), inputFragFrac_HistMap[iso + "_" + detName]->GetBinError(bin_idx)};
            norm_fracs_vw[iso] = {inputNormFrac_HistMap[iso + "_" + detName]->GetBinContent(bin_idx), inputNormFrac_HistMap[iso + "_" + detName]->GetBinError(bin_idx)};
        }
        frag_fracs_vw["Be10"] = {1.0 - frag_fracs_vw["Be7"].val - frag_fracs_vw["Be9"].val, std::hypot(frag_fracs_vw["Be7"].err, frag_fracs_vw["Be9"].err)};
        norm_fracs_vw["Be10"] = {1.0 - norm_fracs_vw["Be7"].val - norm_fracs_vw["Be9"].val, std::hypot(norm_fracs_vw["Be7"].err, norm_fracs_vw["Be9"].err)};

        double raw_FracNuc_val = inputFragNucCounts_HistMap[detName]->GetBinContent(bin_idx+1);
        cout<<raw_FracNuc_val<<endl;
        ValueWithError raw_FragNucCount_vw(raw_FracNuc_val, (raw_FracNuc_val >=0 ? std::sqrt(raw_FracNuc_val) : 0.0) );

        double raw_L1Souce_val = inputL1SourceCounts_HistMap[detName]->GetBinContent(bin_idx+1);
        cout<<inputL1SourceCounts_HistMap[detName]->GetBinLowEdge(bin_idx+1)<<" "<<inputL1SourceCounts_HistMap[detName]->GetBinContent(bin_idx+1)<<endl;
        ValueWithError raw_L1Souce_vw (raw_L1Souce_val, (raw_L1Souce_val >=0 ? std::sqrt(raw_L1Souce_val) : 0.0) );
        
        ValueWithError be_ratio_vw(h_be_ratio_src->GetBinContent(bin_idx+1), h_be_ratio_src->GetBinError(bin_idx+1));
        ValueWithError boron_ratio_vw(h_boron_ratio_src->GetBinContent(bin_idx+1), h_boron_ratio_src->GetBinError(bin_idx+1));
        cout<<h_be_ratio_src->GetBinLowEdge(bin_idx+1)<<" "<<h_be_ratio_src->GetBinContent(bin_idx+1)<<endl;
        
        ValueWithError contam_Be_vw = isISS ? raw_L1Souce_vw * be_ratio_vw : ValueWithError(0.0, 0.0);
        ValueWithError boron_count_final_vw = isISS ? raw_L1Souce_vw * boron_ratio_vw : raw_L1Souce_vw;

        // 新增：计算B10贡献的Be7和Be9
        // 从MC获取B10碎裂产生的Be7和Be9
        double mc_B10_Be7_val = B10MC_Be7_HistMap[detName]->GetBinContent(bin_idx+1);
        double mc_B10_Be9_val = B10MC_Be9_HistMap[detName]->GetBinContent(bin_idx+1);
        double mc_B10_source_val = B10MC_Source_HistMap[detName]->GetBinContent(bin_idx+1);
        
        ValueWithError mc_B10_Be7_vw(mc_B10_Be7_val, 0.0);
        ValueWithError mc_B10_Be9_vw(mc_B10_Be9_val, 0.0);
        ValueWithError mc_B10_source_vw(mc_B10_source_val, 0.0);
        
        // 假设B10占总Boron的30%
        const double B10_fraction = 0.3;
        ValueWithError data_B10_count_vw = boron_count_final_vw * B10_fraction;
        
        // 按比例计算实际数据中B10贡献的Be7和Be9
        ValueWithError B10_contrib_Be7_vw(0.0, 0.0);
        ValueWithError B10_contrib_Be9_vw(0.0, 0.0);
        
        if (mc_B10_source_vw.val > 0) {
            ValueWithError scale_factor = data_B10_count_vw / mc_B10_source_vw;
            B10_contrib_Be7_vw = mc_B10_Be7_vw * scale_factor;
            B10_contrib_Be9_vw = mc_B10_Be9_vw * scale_factor;
        }
        
        // 计算B11贡献的Be总数（扣除B10贡献）
        ValueWithError total_Be_from_B10 = B10_contrib_Be7_vw + B10_contrib_Be9_vw;
        
        for (const auto& iso : isotopes) {
            raw_iso_counts_vw[iso] = raw_FragNucCount_vw * frag_fracs_vw[iso];
            
            ValueWithError contamination_term_for_iso = contam_Be_vw * norm_fracs_vw[iso];
            pure_iso_counts_vw[iso] = raw_iso_counts_vw[iso] - contamination_term_for_iso;
            
            ValueWithError current_pure_iso_ratio;
            if (raw_iso_counts_vw[iso].val != 0.0) {
                ValueWithError B_over_A = contamination_term_for_iso / raw_iso_counts_vw[iso];
                current_pure_iso_ratio.val = 1.0 - B_over_A.val;
                current_pure_iso_ratio.err = B_over_A.err;
            } else {
                current_pure_iso_ratio = {0.0, 0.0}; 
            }
            pure_iso_ratios_vw_map[iso] = current_pure_iso_ratio;

            double B11Frac = h_B11_Frac->GetBinContent(bin_idx);
            double B11FracError = h_B11_Frac->GetBinError(bin_idx);
            if(B11Frac <= 0) {B11Frac = 0.7; B11FracError = 0;}
            ValueWithError b11_frac_vw(B11Frac, 0);
            ValueWithError b10_frac_vw(1-B11Frac, 0);

            ValueWithError BeIso_B11Source_vw = isISS ? pure_iso_counts_vw[iso] / (b11_frac_vw*boron_count_final_vw) : pure_iso_counts_vw[iso] / (boron_count_final_vw); 
            ValueWithError BeIso_B10Source_vw = isISS ? pure_iso_counts_vw[iso] / (b10_frac_vw*boron_count_final_vw) : pure_iso_counts_vw[iso] / (boron_count_final_vw); 
            ValueWithError BeIso_BorSource_vw = isISS ? pure_iso_counts_vw[iso] / (boron_count_final_vw) : pure_iso_counts_vw[iso] / (boron_count_final_vw); 

            fillHistBin(outHists[iso].fragmented_counts, bin_idx, pure_iso_counts_vw[iso]);
            fillHistBin(outHists[iso].fragmented_fraction, bin_idx, frag_fracs_vw[iso]);
            fillHistBin(outHists[iso].normal_fraction, bin_idx, norm_fracs_vw[iso]);
            fillHistBin(outHists[iso].realfrac_totalfrac_ratio, bin_idx, pure_iso_ratios_vw_map[iso]);
            fillHistBin(outHists[iso].raw_isotope_counts, bin_idx, raw_iso_counts_vw[iso]);
            fillHistBin(outHists[iso].BeIso_B11Source_Ratio, bin_idx, BeIso_B11Source_vw);
            fillHistBin(outHists[iso].BeIso_B10Source_Ratio, bin_idx, BeIso_B10Source_vw);
            fillHistBin(outHists[iso].BeIso_BorSource_Ratio, bin_idx, BeIso_BorSource_vw);
        }
        
        // 新增：计算Be10 from B11 / Total Be from B11的比值
        ValueWithError Be10_from_B11 = pure_iso_counts_vw["Be10"];  // Be10全部来自B11
        ValueWithError Be_total_from_B11 = pure_iso_counts_vw["Be7"] + pure_iso_counts_vw["Be9"] + pure_iso_counts_vw["Be10"] - total_Be_from_B10;
        
        ValueWithError Be10_B11_ratio(0.0, 0.0);
        if (Be_total_from_B11.val > 0) {
            Be10_B11_ratio = Be10_from_B11 / Be_total_from_B11;
        }
        
        fillHistBin(h_Be10_from_B11_ratio, bin_idx, Be10_B11_ratio);
        
        fillHistBin(h_L1SourceCounts, bin_idx, raw_L1Souce_vw);
        fillHistBin(h_contamination_Be, bin_idx, contam_Be_vw);
        fillHistBin(h_boron_counts_out, bin_idx, boron_count_final_vw);
        fillHistBin(h_fragment_Be_raw, bin_idx, raw_FragNucCount_vw);

        std::cout << "Energy bin: " << outHists["Be7"].fragmented_counts->GetXaxis()->GetBinLowEdge(bin_idx) << " - " 
                  << outHists["Be7"].fragmented_counts->GetXaxis()->GetBinUpEdge(bin_idx) << " GeV/n" << std::endl;
        std::cout << "Detector: " << detName << std::endl;
        std::cout << "Total Be and B: " << raw_L1Souce_vw.val << " +/- " << raw_L1Souce_vw.err << std::endl;
        std::cout << "Raw Be Count: " << raw_FragNucCount_vw.val << std::endl;
        std::cout << "Contamination Be: " << contam_Be_vw.val << " +/- " << contam_Be_vw.err << std::endl;
        std::cout << "Boron count: " << boron_count_final_vw.val << " +/- " << boron_count_final_vw.err << std::endl;
        std::cout << "B10 contribution to Be7: " << B10_contrib_Be7_vw.val << " +/- " << B10_contrib_Be7_vw.err << std::endl;
        std::cout << "B10 contribution to Be9: " << B10_contrib_Be9_vw.val << " +/- " << B10_contrib_Be9_vw.err << std::endl;
        std::cout << "Be10 from B11 / Total Be from B11: " << Be10_B11_ratio.val << " +/- " << Be10_B11_ratio.err << std::endl;
        for (const auto& iso : isotopes) {
            std::cout << iso << " fragmented: " << pure_iso_counts_vw[iso].val << " +/- " << pure_iso_counts_vw[iso].err
                      << " (" << (pure_iso_ratios_vw_map[iso].val * 100.0) << "%)" << std::endl;
        }
        std::cout << "-----------------------------------------------------" << std::endl;
    }

    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    auto createCanvas = [&](const char* name, const char* title, int divX, int divY) {
        auto canvas = new TCanvas(name, title, 350 * divX, 300 * divY);
        canvas->Divide(divX, divY);
        allObjectsToSave.push_back(canvas);
        return canvas;
    };
    
    std::map<std::string, EColor> iso_colors = {
        {"Be7", (EColor)(kGreen + 2)},
        {"Be9", kBlue},
        {"Be10", kMagenta}
    };

    TCanvas* c_frag_counts = createCanvas("canvas_frag_counts", "Fragmented Beryllium Isotopes Counts", 1, 3);
    for (size_t i = 0; i < isotopes.size(); ++i) {
        c_frag_counts->cd(i + 1);
        styleHist(outHists[isotopes[i]].fragmented_counts, iso_colors[isotopes[i]]);
        outHists[isotopes[i]].fragmented_counts->Draw("E");
    }

    TCanvas* c_pure_ratios = createCanvas("canvas_pure_ratios", "Pure Fragmented Beryllium Ratio", 1, 3);
    for (size_t i = 0; i < isotopes.size(); ++i) {
        c_pure_ratios->cd(i + 1);
        styleHist(outHists[isotopes[i]].realfrac_totalfrac_ratio, iso_colors[isotopes[i]]);
        outHists[isotopes[i]].realfrac_totalfrac_ratio->Draw("E");
    }
    
    TCanvas* c_counts_intermediate = createCanvas("canvas_counts_intermediate", "Intermediate Counts", 2, 2);
    c_counts_intermediate->cd(1); styleHist(h_L1SourceCounts, kBlack); h_L1SourceCounts->Draw("E");
    c_counts_intermediate->cd(2); styleHist(h_contamination_Be, kRed); h_contamination_Be->Draw("E");
    c_counts_intermediate->cd(3); styleHist(h_boron_counts_out, (EColor)(kOrange+7)); h_boron_counts_out->Draw("E");
    c_counts_intermediate->cd(4); styleHist(h_fragment_Be_raw, kGreen); h_fragment_Be_raw->Draw("E");

    styleHist(outHists["Be10"].BeIso_B11Source_Ratio, isISS ? (EColor)kBlack : (EColor)(kOrange+1));

    TCanvas* c_fractions = createCanvas("canvas_fractions", "Isotope Fractions", 2, 3);
    for (size_t i = 0; i < isotopes.size(); ++i) {
        const auto& iso = isotopes[i];
        c_fractions->cd(i * 2 + 1); 
        styleHist(outHists[iso].fragmented_fraction, iso_colors[iso], (EMarkerStyle)20);
        outHists[iso].fragmented_fraction->Draw("E");
        
        c_fractions->cd(i * 2 + 2); 
        styleHist(outHists[iso].normal_fraction, iso_colors[iso], (EMarkerStyle)24);
        outHists[iso].normal_fraction->Draw("E");
    }
    
    TCanvas* c_raw_iso_counts = createCanvas("canvas_raw_iso_counts", "Raw Isotope Counts", 1, 3);
    for (size_t i = 0; i < isotopes.size(); ++i) {
        c_raw_iso_counts->cd(i + 1);
        styleHist(outHists[isotopes[i]].raw_isotope_counts, iso_colors[isotopes[i]]);
        outHists[isotopes[i]].raw_isotope_counts->Draw("E");
    }
    
    // 新增：绘制Be10 from B11 / Total Be from B11比值
    TCanvas* c_Be10_B11_ratio = createCanvas("canvas_Be10_B11_ratio", "Be10 from B11 Ratio", 1, 1);
    c_Be10_B11_ratio->SetLogx();
    c_Be10_B11_ratio->cd(1);
    styleHist(h_Be10_from_B11_ratio, kRed);
    h_Be10_from_B11_ratio->GetYaxis()->SetRangeUser(0,1);
    //h_Be10_from_B11_ratio->GetXaxis()->SetRangeUser(0.01,16.5);
    h_Be10_from_B11_ratio->Draw("E");

    outputFile->cd();
    for (TObject* obj : allObjectsToSave) {
        if (obj) obj->Write(obj->GetName(), TObject::kOverwrite);
    }
    
    std::cout << "Results saved to " << outputFile->GetName() << std::endl;
}

void CalFragBe() {
    try {
        calculateFragmentedBeryllium(true);
        //calculateFragmentedBeryllium(false);
    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
    }
}