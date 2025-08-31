/***********************************************************
 * File: Tools.h
 *
 * Modern C++ header file for AMS analysis tools.
 * Consolidates functions from legacy Tool1.h, Tool2.h,
 * and RootFunc.h.
 *
 * History:
 * 20250831 - Ported and refactored by Gemini
 ***********************************************************/

#pragma once

#include <string>
#include <vector>
#include <map>
#include <optional>
#include <array>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <sys/resource.h>

// ROOT includes
#include <TF1.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <IsoToolbox/PhysicsConstants.h>


namespace IsoToolbox {
namespace Tools {

// DATA STRUCTURES
struct Point2D {
    double x;
    double y;
};

struct MassResult {
    double beta;
    double gamma;
    double ek;
    double invMass;
};

// PHYSICS & KINEMATICS
double betaToKineticEnergy(double beta);
double kineticEnergyToBeta(double kineticEnergy);
double rigidityToBeta(double rigidity, int charge, int mass, bool isElectron = false);
double betaToRigidity(double beta, int charge, int mass, bool isElectron = false);
double kineticEnergyToRigidity(double ek_per_nucleon, int z, int a);
double rigidityToKineticEnergy(double rig_gv, int z, int a);
double dR_dEk(double ek_per_nucleon, int z, int a);
double dEk_dR(double rig_gv, int z, int a);
MassResult calculateMass(double beta, double alpha, double innerRig, int charge);
double calculateWeight(double mmom, int charge, bool isISS);

// EVENT SELECTION & BINNING
bool isValidBeta(double beta);
int findBin(const std::vector<double>& bins, double value);
int findRigidityLowerBound(double rigidity, double *xAxis, int range);

// GEOMETRY & UTILITIES
std::string convertElementName(const std::string& fullName, bool isISS);
double calculateAverage(const double* values, int count, double ignoreValue = -1000000.001);
void modifyPositionByZ(double targetZ, std::array<double, 3>& position, double theta, double phi);
std::optional<Point2D> calculateXYAtZ(const Float_t positions[9][3], const Float_t directions[9][3], double zpl, int trackIndex = 0);
std::optional<Point2D> calculateXYAtNaFZ(const std::array<double, 3>& pos1, const std::array<double, 3>& pos2, double nafZ);
double getCurrentRSS_MB();
std::unique_ptr<TH1D> getEnergySlice(int plusIbin, int Nrebin, TH2D* h2d, int binIndex, const std::string& name);

// FITTING FUNCTIONS
Double_t landauGaussianFit(Double_t *x, Double_t *par);
double extendedExponentialGaussianFit(double *x, double *p);

// PLOTTING & ANALYSIS
void calculatePull(RooPlot* frame, TGraphErrors* pullGraph, double fit_min, double fit_max);
void setupPullPlot(TGraphErrors* pullGraph, double fit_min, double fit_max);
std::unique_ptr<TLegend> createLegend();
std::unique_ptr<TH1D> extendHistogram(const TH1D* hist, double new_min, double new_max);

} // namespace Tools
} // namespace IsoToolbox
