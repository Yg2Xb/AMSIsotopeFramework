#ifndef HIST_TEMP_FIT_H
#define HIST_TEMP_FIT_H

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooHistPdf.h>
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include "RooAbsReal.h"
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <string>
#include <vector>
#include <map>
#include <memory>

// Configuration for each element to be fitted
struct IsotopeConfig {
    std::string name;                   // Element name, e.g., "Li", "Be"
    int charge;                         // Charge of the element
    std::vector<int> masses;            // Isotope masses, e.g., {7, 9, 10} for Be
    std::vector<double> initFractions;  // Initial fractions for the fit
    std::vector<double> fitRangeLow;    // Fit range lower bound for [TOF, NaF, AGL]
    std::vector<double> fitRangeUp;     // Fit range upper bound for [TOF, NaF, AGL]
    int nFitParams;                     // Number of free fraction parameters (usually masses.size() - 1)
};

// Constants and configurations for the fitting process
class IsoFitConstants {
public:
    static const std::map<std::string, IsotopeConfig> configs;
};

// Main function declaration
void HistTempFit(const std::string& isotype = "Be", 
                 int UseMass = 7, 
                 bool useUnbiasedChain = true, 
                 bool usePureTemplates = true,
                 int rebinX = 1, 
                 int ProNbin = 1, 
                 bool fitFragMass = false,
                 const string& sourceName = ""
                );

// Entry point function
void HistTempFit_hist();

// ==================== Global Configurations ====================

// Defines the properties for each element available for fitting
const char* DetName[] = {"TOF", "NaF", "AGL"};
const double DetRanges[3][2] = {{0.25, 1.5}, {0.61, 6.10}, {2.50, 23.0}};

const std::map<std::string, IsotopeConfig> IsoFitConstants::configs = {
    {"Li", {"Li", 3, {6, 7}, {0.5}, 
           {0.095, 0.095, 0.095}, {0.24, 0.242, 0.242}, 1}},

    {"Be", {"Be", 4, {7, 9, 10}, {0.6, 0.33}, 
           {0.05, 0.05, 0.05}, {0.23, 0.23, 0.23}, 2}},

    {"B", {"B", 5, {10, 11}, {0.3}, 
          {0.05, 0.05, 0.05}, {0.15, 0.15, 0.15}, 1}}
};

// ==================== Plotting Helper Functions ====================

inline void setLegend(TLegend *legend) {
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.035);
}

inline void setPaveText(TPaveText *pt) {
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.03);
}

#endif // HIST_TEMP_FIT_H