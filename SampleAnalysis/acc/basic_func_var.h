#ifndef BASIC_FUNC_VAR_H
#define BASIC_FUNC_VAR_H

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <array>
#include "RootFunc.h"

// Style settings
void acc_goodlogon() {
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(1.);
    gStyle->SetLineWidth(1);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
}

// Physics conversion functions
double convertBetaToRig(double beta, int Z, double a) {
    if (fabs(beta) >= 1) return -100000.;
    double mass = a * 0.931;
    return beta * mass / (sqrt(1 - beta * beta) * Z);
}

double convertEkToBeta(double ek) {
    double gamma = ek / 0.931 + 1;
    return sqrt(1 - 1 / (gamma * gamma));
}

double convertBetaToEk(double beta) {
    if (beta < -1 || beta > 1) {
        std::cout << "Wrong beta input for convertBetaToEk()!!!" << std::endl;
        return 0;
    }
    double gamma = 1 / sqrt(1 - beta * beta);
    return (gamma - 1) * 0.931;
}

// Constants and configurations
const double _xlow = 1.;
const double _xup = 2000.;

// Detector configurations
struct DetectorConfig {
    const char* name;
    int color;
    double betaCut;
    double ekMin;
    double ekMax;
    int index;
};

const std::array<DetectorConfig, 3> DETECTORS = {{
    {"tof", kBlue,      0.5,   0.5,  1.8, 0},
    {"NaF", kOrange+1,  0.75,  1.0,  5.0, 1},
    {"AGL", kGreen+2,   0.953, 2.24, 14.0, 2}
}};

// Isotope configurations
struct IsotopeConfig {
    std::string name;       // Li, B, Be etc.
    std::string histPrefix; // Lithium, Boron, Berlium
    int Z;                  // atomic number
    std::vector<int> masses;// mass numbers
    const double* bins;     // energy bins
    int nBins;             // number of bins
};

const double Be_bins[] = {0.08, 0.13, 0.17, 0.21, 0.27, 0.33, 0.41, 0.49, 0.59, 0.70, 0.82, 
    0.96, 1.11, 1.28, 1.47, 1.68, 1.91, 2.16, 2.44, 2.73, 3.06, 3.41, 3.79, 4.20, 4.65, 5.14, 
    5.64, 6.18, 6.78, 7.42, 8.12, 8.86, 9.66, 10.51, 11.45, 12.45, 13.50, 14.65, 15.84, 17.14, 
    18.54, 20.04, 21.64, 23.34, 25.19, 27.13, 29.23, 31.48, 33.93, 36.53, 39.33, 42.33, 45.58, 
    49.08, 53.08, 57.08, 61.58, 66.58, 72.57, 79.07, 86.57, 95.07, 104.57, 115.57, 128.57, 
    144.57, 164.07, 188.57, 219.57, 261.57, 329.07, 439.07, 649.07, 1649.07};
const double LiB_bins[] = {0.08, 0.13, 0.17, 0.21, 0.27, 0.33, 0.41, 0.49, 0.59, 0.70, 0.82, 
    0.96, 1.11, 1.28, 1.47, 1.68, 1.91, 2.16, 2.44, 2.73, 3.06, 3.41, 3.79, 4.20, 4.65, 5.14, 
    5.64, 6.18, 6.78, 7.42, 8.12, 8.86, 9.66, 10.51, 11.45, 12.45, 13.50, 14.65, 15.84, 17.14, 
    18.54, 20.04, 21.64, 23.34, 25.19, 27.13, 29.23, 31.48, 33.93, 36.53, 39.33, 42.33, 45.58, 
    49.08, 53.08, 57.08, 61.58, 66.58, 72.57, 79.07, 86.57, 95.07, 104.57, 115.57, 128.57, 
    144.57, 164.07, 188.57, 219.57, 261.57, 329.07, 439.07, 649.07, 1649.07};

double TofEdgeEk = convertBetaToEk(0.5);
double NaFEdgeEk = convertBetaToEk(0.75);
double AGLEdgeEk = convertBetaToEk(0.953);
std::array<std::vector<double>, 3> xpointsEk = {{{0.25, 0.4, 0.63, 0.89, 1.41, 2.00, 3.16, 5.62, 10.00, 30, 60, 100, 300, 650},
                                                 {1.4 * NaFEdgeEk, 0.75, 0.92, 1.25, 1.41, 2.00, 3.16, 5.62, 10.00, 30, 60, 100, 300, 650},
                                                 {1.1 * AGLEdgeEk, 2.3, 3.2, 3.8, 4.5, 5.62, 10.00, 30, 60, 100, 300, 650}}};

const std::map<std::string, IsotopeConfig> isotope_configs = {
    {"Li", {"Li", "Lithium", 3, {6, 7}, LiB_bins, 74}},
    {"Be", {"Be", "Beryllium", 4, {7, 9, 10}, Be_bins, 74}},
    {"B",  {"B",  "Boron",   5, {10, 11}, LiB_bins, 74}}
};

const int colors[] = {600, 880, 632, 800-3, 900+10, 416-1, 1, 416-3, 840+3, 632, 416 + 2};
const int nColors = sizeof(colors) / sizeof(colors[0]);

#endif // BASIC_FUNC_VAR_H