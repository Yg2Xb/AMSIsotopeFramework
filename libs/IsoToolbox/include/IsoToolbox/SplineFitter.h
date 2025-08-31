/***********************************************************
 * File: SplineFitter.h
 *
 * Header file for the SplineFitter class, a tool for
 * performing spline fits on ROOT histograms and graphs.
 * Ported from legacy QSplineFit.C.
 *
 * History:
 * 20250831 - Ported and refactored by Gemini
 ***********************************************************/

#pragma once

#include "TH1.h"
#include "TF1.h"
#include "TSpline.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TDirectory.h"
#include "TVirtualFitter.h"
#include <iostream>

namespace IsoToolbox {
namespace Tools {

///--Spline Fit Tool
class SplineFitter {
public:
    enum FitMode {
        LinearXY = 0x1,      // y=f(X)
        LogX = 0x2,          // y=f(Log(X))
        LogY = 0x4,          // Log(Y)=f(X)
        LogX2 = 0x8,         // y=f(Log(X)) using geometric mean for bin center
        ExtrapolateLB = 0x10, // Linear Extrapolate Begin
        ExtrapolateLE = 0x20, // Linear Extrapolate End
        ExtrapolateFlux = 0x40 // Y=CX^a
    };

public:
    SplineFitter();
    SplineFitter(const char *name, TH1 *h, int Nnodes, double *nodes, int sfitmode = SplineFitter::LinearXY, const char *sopt = "b1e1", double svalbeg = 0, double svalend = 0, double serror = 1e-12);
    SplineFitter(const char *name, TGraphErrors *gr, int Nnodes, double *nodes, int sfitmode = SplineFitter::LinearXY, const char *sopt = "b1e1", double svalbeg = 0, double svalend = 0, double serror = 1e-12);
    SplineFitter(const char *name, int Nnodes, double *nodes, double *ynodes, int sfitmode = SplineFitter::LinearXY, const char *sopt = "b1e1", double svalbeg = 0, double svalend = 0, double serror = 1e-12); //Fix Nodes
    virtual ~SplineFitter();

    void SetName(const char *name);
    void SetHisto(TH1 *h);
    void SetGraph(TGraphErrors *gr);
    void SetNode(int Nnodes, double *nodes, double *ynodes = 0);
    int GetNodeX(double *xnodes);
    void SetSplinePar(int sfitmode = SplineFitter::LinearXY, const char *sopt = "b1e1", double svalbeg = 0, double svalend = 0, double serror = 1e-12);
    void SetVerbose(int verbose = 0);
    TF1 *GetFun(double xmin = 0, double xmax = 0);
    virtual TFitResultPtr DoFit(Option_t *option = "", Option_t *goption = "", Double_t xmin = 0, Double_t xmax = 0);
    int FillFlux(TH1 *hflux, double scaleindex = 0, double scalenorm = 1);
    double operator()(double *x, double *par);

private:
    void Init();
    void Clear();
    TGraphErrors* HistoToGraph(int cb = 0);

    const char *qname;
    int qnnodes;
    double *qxnodes;
    double *qynodes;

    const char *_sopt;
    double _serror;
    double _svalbeg;
    double _svalend;
    int _sfitmode;
    int _verbose;

    TH1 *qh;
    TF1 *qfun;
    TSpline3 *qspline;
    TGraphErrors *qgr;
};

// Helper function for easier spline fitting, ported from legacy RootFunc.h
TF1* splineFit(TH1* hist, double* xpoints, int npoint, int qOption, TString splineOption = "b1e1", TString name = "SplineFunc", double fitRangeLow = -9999, double fitRangeHigh = -9999);
TF1* splineFit(TGraph* graph, double* xpoints, int npoint, int qOption, TString splineOption = "b1e1", TString name = "SplineFunc", double fitRangeLow = -9999, double fitRangeHigh = -9999);

}
} 
