#pragma once
#include <fstream>
#include "TSystem.h"
#include <vector>
#include <map>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TCollection.h>
#include <TChainElement.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TString.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSpline.h>
#include "QSplineFit.C"
//#include "../FitFunc.h"
using namespace std;
TString TxtFunc_Dir = "./tmpstuff/";
TString TxtFunc_Name = "RootFunc.txt";
vector<int> emptyintvector;
vector<double> emptydoublevector;
vector<TH1F*> emptyTH1Fvector;
const double NANNumber = -9999;
TSpline3* globalSpline = NULL;

TF1* SplineFit(TH1* haim, double* xpoints, int npoint, int QOption, TString SplineOption = "b1e1", TString Name = "SplineFunc", double FitRangeLow = NANNumber, double FitRangeHigh = NANNumber)
{
    QSplineFit* spfit = new QSplineFit("spfit",(TH1*)haim,npoint,xpoints,QOption,SplineOption);
    TF1* spfit_func;

    if(FitRangeLow == NANNumber && FitRangeHigh == NANNumber)
    {
        spfit->DoFit();
        spfit_func = spfit->GetFun();
    }
    else
    {
        spfit->DoFit("","",FitRangeLow,FitRangeHigh);
        spfit_func = spfit->GetFun(FitRangeLow,FitRangeHigh);
    }

    spfit_func->SetName(Name);

    return spfit_func;
}

TF1* SplineFit(TGraph* gaim, double* xpoints, int npoint, int QOption, TString SplineOption = "b1e1", TString Name = "SplineFunc", double FitRangeLow = NANNumber, double FitRangeHigh = NANNumber)
{
    QSplineFit* spfit = new QSplineFit("spfit",(TGraphErrors*)gaim,npoint,xpoints,QOption,SplineOption);
    TF1* spfit_func;

    if(FitRangeLow == NANNumber && FitRangeHigh == NANNumber)
    {
        spfit->DoFit();
        spfit_func = spfit->GetFun();
    }
    else
    {
        spfit->DoFit("","",FitRangeLow,FitRangeHigh);
        spfit_func = spfit->GetFun(FitRangeLow,FitRangeHigh);
    }

    spfit_func->SetName(Name);

    return spfit_func;
}

int findRigidityLowerBound(double rigidity, double *xAxis, int range)
{
    int pos(-1);
    for (int i = 0; i != range; ++i)
    {
        if (rigidity >= xAxis[i] and rigidity < xAxis[i + 1])
        {
            pos = i;
            break;
        }
    }
    return pos;
}
