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
#include <Math/GSLIntegrator.h>
#include <TInterpreter.h>
#include "./QSplineFit.C"

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

TH1* HistRatio(TH1* hnum, TH1* hden, TH1* hratio = NULL, TString Option = "", TString Title = "", bool DrawFlag = true) //numerator denominator fraction , Option: '' or 'noerr'
{
    int nhb_num = hnum->GetNbinsX(); int nhb_den = hden->GetNbinsX();
    if(!hratio)
    {
        hratio = (TH1F*)hden->Clone("hratio");
    }

    for(int ihb_den=1; ihb_den <= nhb_den; ihb_den++)
    {
        int ihb_num = hnum->FindBin(hden->GetBinCenter(ihb_den));

        bool Flowunder = false; bool Flowover = false; bool Flowin = false;
        if(ihb_num > 1 && ihb_den == 1) Flowunder = true;
        if(ihb_num < nhb_num && ihb_den == nhb_den) Flowover = true;
        if(ihb_num <= 0 || ihb_num >= nhb_num+1) Flowin = true; 

        double xlow_num, xlow_den, xup_num, xup_den;
        double y_num, y_den, y_ratio;
        double yerr_num, yerr_den, yerr_ratio;

        xlow_num = hnum->GetBinLowEdge(ihb_num);
        xup_num = hnum->GetBinWidth(ihb_num)+xlow_num;
        xlow_den = hden->GetBinLowEdge(ihb_den);
        xup_den = hden->GetBinWidth(ihb_den)+xlow_den;

        y_num = hnum->GetBinContent(ihb_num); yerr_num = hnum->GetBinError(ihb_num);
        y_den = hden->GetBinContent(ihb_den); yerr_den = hden->GetBinError(ihb_den);
        y_ratio = y_num/y_den; yerr_ratio = y_ratio*TMath::Sqrt((yerr_num*yerr_num)/(y_num*y_num) + (yerr_den*yerr_den)/(y_den*y_den));
        if(y_num == 0 || y_den == 0 || Flowin) { y_ratio = 0; yerr_ratio = 0; } 

        hratio->SetBinContent(ihb_den,y_ratio); hratio->SetBinError(ihb_den,yerr_ratio);

        if(Option.Contains("noerr")) hratio->SetBinError(ihb_den,0);


        if(Flowunder) {
            for(int itmp = 1; itmp < ihb_num; itmp++)
            {
                double tmplow_num; double tmpup_num;
                double tmpy_num; double tmpyerr_num;
                tmplow_num = hnum->GetBinLowEdge(itmp); tmpup_num = hnum->GetBinWidth(itmp)+tmplow_num;
                tmpy_num = hnum->GetBinContent(itmp); tmpyerr_num = hnum->GetBinError(itmp);
                if(DrawFlag) cout << Form("Numerator (%d) [%g,%g] (%g #pm %g) v.s Denominator (Non) [Non,Non] Non",itmp,tmplow_num,tmpup_num,tmpy_num,tmpyerr_num) << endl;
            }
        }

        if(DrawFlag)
        {
            if(Flowin) cout << Form("Numerator (Non) [Non,Non] Non v.s Denominator (%d) [%g,%g] (%g #pm %g)",ihb_den,xlow_den,xup_den,y_den,yerr_den) << endl;
            else cout << Form("Numerator (%d) [%g,%g] (%g #pm %g) v.s Denominator (%d) [%g,%g] (%g #pm %g)",ihb_num,xlow_num,xup_num,y_num,yerr_num,ihb_den,xlow_den,xup_den,y_den,yerr_den) << endl;
        }

        if(Flowover) {
            for(int itmp = ihb_num+1; itmp <= nhb_num; itmp++)
            {
                double tmplow_num; double tmpup_num;
                double tmpy_num; double tmpyerr_num;
                tmplow_num = hnum->GetBinLowEdge(itmp); tmpup_num = hnum->GetBinWidth(itmp)+tmplow_num;
                tmpy_num = hnum->GetBinContent(itmp); tmpyerr_num = hnum->GetBinError(itmp);
                if(DrawFlag) cout << Form("Numerator (%d) [%g,%g] (%g#pm%g) v.s Denominator (Non) [Non,Non]",itmp,tmplow_num,tmpup_num,tmpy_num,tmpyerr_num) << endl;
            }
        }
    }

    hratio->SetLineWidth(3); hratio->SetMarkerStyle(8);
    hratio->SetLineColor(4); hratio->SetMarkerColor(4);
    
    hratio->GetYaxis()->SetTitle("ratio");
    if(Title != "") {
        hratio->SetTitle(Title);
    }
    else {
        TString t1 = hnum->GetTitle();
        TString t2 = hden->GetTitle();
        hratio->SetTitle(Form("%s <Divide> %s",t1.Data(),t2.Data()));
    }

    if(DrawFlag == true)
    {
        new TCanvas(); 
        hratio->Draw();
    }

    return hratio;
}
TH1* HistOperate(TH1* hnum, TH1* hden, TString FuncString, TString FuncErrString = "", TH1* hratio = NULL, TString Option = "", TString Title = "", bool DrawFlag = true) // FuncString: TF1 formula string, hnum relate to [0], hden relate to [1]; Option: '' or 'noerr' ; FuncErrString: hnum value relate to [0], hden value relate to [1], hnum error relate to [2], hden error relate to [3];
//Example: positron fraction: FuncString:"[0]/([0]+[1])", FuncErrString:"TMath::Sqrt( ([1]*[1]*[0]*[0]) + ([0]*[0]*[1]*[1]) ) / ( ([0]+[1])*([0]+[1]) )"
{
    TString OutputFuncFormat;
    TF1* ffunc = new TF1("ffunc",FuncString,0,10);
    TF1* ferrfunc;

    OutputFuncFormat = "Func: "+FuncString;
    if(FuncErrString == "")
    {
        ferrfunc = NULL;
        OutputFuncFormat += " ; FuncErr: Set0";
    }
    else
    {
        ferrfunc = new TF1("ferrfunc",FuncErrString,0,10);
        OutputFuncFormat += " ; FuncErr: "+FuncErrString;
    }

    OutputFuncFormat.ReplaceAll("[0]","val0");
    OutputFuncFormat.ReplaceAll("[1]","val1");
    OutputFuncFormat.ReplaceAll("[2]","err0");
    OutputFuncFormat.ReplaceAll("[3]","err1");

    int nhb_num = hnum->GetNbinsX(); int nhb_den = hden->GetNbinsX();
    if(!hratio)
    {
        hratio = new TH1F((TH1F)*((TH1F*)hden));
        hratio->SetName("hratio");
    }

    for(int ihb_den=1; ihb_den <= nhb_den; ihb_den++)
    {
        int ihb_num = hnum->FindBin(hden->GetBinCenter(ihb_den));

        bool Flowunder = false; bool Flowover = false; bool Flowin = false;
        if(ihb_num > 1 && ihb_den == 1) Flowunder = true;
        if(ihb_num < nhb_num && ihb_den == nhb_den) Flowover = true;
        if(ihb_num <= 0 || ihb_num >= nhb_num+1) Flowin = true; 

        double xlow_num, xlow_den, xup_num, xup_den;
        double y_num, y_den, y_ratio;
        double yerr_num, yerr_den, yerr_ratio;

        xlow_num = hnum->GetBinLowEdge(ihb_num);
        xup_num = hnum->GetBinWidth(ihb_num)+xlow_num;
        xlow_den = hden->GetBinLowEdge(ihb_den);
        xup_den = hden->GetBinWidth(ihb_den)+xlow_den;

        y_num = hnum->GetBinContent(ihb_num); yerr_num = hnum->GetBinError(ihb_num);
        y_den = hden->GetBinContent(ihb_den); yerr_den = hden->GetBinError(ihb_den);
////        y_ratio = y_num/y_den; yerr_ratio = y_ratio*TMath::Sqrt((yerr_num*yerr_num)/(y_num*y_num) + (yerr_den*yerr_den)/(y_den*y_den));
        ffunc->SetParameters(y_num,y_den); y_ratio = ffunc->Eval(1); // ffunc is a const value function, don't change with x;
        if(ferrfunc)
        {
            ferrfunc->SetParameters(y_num,y_den,yerr_num,yerr_den);
            yerr_ratio = ferrfunc->Eval(1);
        }
        else
        {
            yerr_ratio = 0;
        }
        if(y_num == 0 || y_den == 0 || Flowin) { y_ratio = 0; yerr_ratio = 0; } 

        hratio->SetBinContent(ihb_den,y_ratio); hratio->SetBinError(ihb_den,yerr_ratio);

        if(Option.Contains("noerr")) hratio->SetBinError(ihb_den,0);


        if(Flowunder) {
            for(int itmp = 1; itmp < ihb_num; itmp++)
            {
                double tmplow_num; double tmpup_num;
                double tmpy_num; double tmpyerr_num;
                tmplow_num = hnum->GetBinLowEdge(itmp); tmpup_num = hnum->GetBinWidth(itmp)+tmplow_num;
                tmpy_num = hnum->GetBinContent(itmp); tmpyerr_num = hnum->GetBinError(itmp);
                if(DrawFlag) cout << Form("h0 (%d) [%g,%g] (%g #pm %g) v.s h1 (Non) [Non,Non] Non ; %s",itmp,tmplow_num,tmpup_num,tmpy_num,tmpyerr_num,OutputFuncFormat.Data()) << endl;
            }
        }

        if(DrawFlag)
        {
            if(Flowin) cout << Form("h0 (Non) [Non,Non] Non v.s h1 (%d) [%g,%g] (%g #pm %g) ; %s",ihb_den,xlow_den,xup_den,y_den,yerr_den,OutputFuncFormat.Data()) << endl;
            else cout << Form("h0 (%d) [%g,%g] (%g #pm %g) v.s h1 (%d) [%g,%g] (%g #pm %g) ; %s",ihb_num,xlow_num,xup_num,y_num,yerr_num,ihb_den,xlow_den,xup_den,y_den,yerr_den,OutputFuncFormat.Data()) << endl;
        }

        if(Flowover) {
            for(int itmp = ihb_num+1; itmp <= nhb_num; itmp++)
            {
                double tmplow_num; double tmpup_num;
                double tmpy_num; double tmpyerr_num;
                tmplow_num = hnum->GetBinLowEdge(itmp); tmpup_num = hnum->GetBinWidth(itmp)+tmplow_num;
                tmpy_num = hnum->GetBinContent(itmp); tmpyerr_num = hnum->GetBinError(itmp);
                if(DrawFlag) cout << Form("h0 (%d) [%g,%g] (%g#pm%g) v.s h1 (Non) [Non,Non] ; %s",itmp,tmplow_num,tmpup_num,tmpy_num,tmpyerr_num,OutputFuncFormat.Data()) << endl;
            }
        }
    }

    hratio->SetLineWidth(3); hratio->SetMarkerStyle(8);
    hratio->SetLineColor(4); hratio->SetMarkerColor(4);
    
    hratio->GetYaxis()->SetTitle("ratio");
    if(Title != "") {
        hratio->SetTitle(Title);
    }
    else {
        TString t1 = hnum->GetTitle();
        TString t2 = hden->GetTitle();
        hratio->SetTitle(Form("%s <Divide> %s",t1.Data(),t2.Data()));
    }

    if(DrawFlag == true)
    {
        new TCanvas(); 
        hratio->Draw();
    }

    return hratio;
}

TGraphErrors* GraphRatio(TGraphErrors* gnum, TGraphErrors* gden, TGraphErrors* gratio = NULL, TString Option = "", TString Title = "", bool DrawFlag = true, bool Estimate1percentage = true) //numerator denominator fraction, Option: '' or 'noerr'
{
    if(!gratio) {gratio = (TGraphErrors*)gden->Clone("gratio");}
    int N_num = gnum->GetN();
    int N_den = gden->GetN();
    int N_ratio = gratio->GetN();

    for(int ip=0; ip < N_ratio; ip++){
        gratio->RemovePoint(0);
    }

    int imatch_num = -1;
    double x_den; double y_den; double xerr_den = 0; double yerr_den = 0;
    double x_num; double y_num; double xerr_num = 0; double yerr_num = 0;
    double y_ratio, yerr_ratio;
    int ip_ratio = -1;

    for(int ip_den=0; ip_den < N_den; ip_den++)
    {
        gden->GetPoint(ip_den,x_den,y_den);
        xerr_den = gden->GetErrorX(ip_den);
        yerr_den = gden->GetErrorY(ip_den);

        if(x_den != x_den || y_den != y_den || xerr_den != xerr_den || yerr_den != yerr_den) continue;

        imatch_num = -1;
        for(int ip_num=0; ip_num < N_num; ip_num++)
        {
            gnum->GetPoint(ip_num,x_num,y_num);
            xerr_num = gnum->GetErrorX(ip_num);
            yerr_num = gnum->GetErrorY(ip_num);

            if(x_num != x_num || y_num != y_num || xerr_num != xerr_num || yerr_num != yerr_num) continue;

            if(x_num == x_den || ( Estimate1percentage && (x_num/x_den) > 0.95 && (x_num/x_den) < 1.05 ) ) 
            {
                imatch_num = ip_num;
                break;
            }
        }

        if(DrawFlag == true)
        {
            if (imatch_num > 0 && imatch_num <= N_num-1 && ip_den == 0) { 
                for(int itmp_num=0; itmp_num < imatch_num; itmp_num++)
                {
                    double xtmp_num, ytmp_num, xtmperr_num, ytmperr_num;

                    gnum->GetPoint(itmp_num,xtmp_num,ytmp_num);
                    xtmperr_num = gnum->GetErrorX(itmp_num);
                    ytmperr_num = gnum->GetErrorY(itmp_num);
                    cout << Form("Numerator (%d) (%g +- %g , %g +- %g) v.s Denominator (Non) Non,Non ",itmp_num,xtmp_num,xtmperr_num,ytmp_num,ytmperr_num) << endl;
                }
            }
        }
        if(imatch_num >= 0)
        {
            ip_ratio++;
            y_ratio = y_num/y_den; yerr_ratio = y_ratio*TMath::Sqrt((yerr_num*yerr_num)/(y_num*y_num) + (yerr_den*yerr_den)/(y_den*y_den));
            if(y_num == 0 || y_den == 0) {y_ratio = 0; yerr_ratio = 0;}

            gratio->SetPoint(ip_ratio,x_den,y_ratio);
            gratio->SetPointError(ip_ratio,0,yerr_ratio);
            if(Option.Contains("noerr")) gratio->SetPointError(ip_ratio,0,0);

            if(DrawFlag) cout << Form("Numerator (%d) (%g +- %g , %g +- %g) v.s Denominator (%d) (%g +- %g , %g +- %g)",imatch_num,x_num,xerr_num,y_num,yerr_num,ip_den,x_den,xerr_den,y_den,yerr_den) << endl;
        }

        if(imatch_num == -1 && DrawFlag == true) {
            cout << Form("Numerator (Non) Non,Non v.s Denominator (%d) (%g +- %g , %g +- %g) ",ip_den,x_den,xerr_den,y_den,yerr_den) << endl;
        }

        if (imatch_num >= 0 && imatch_num < N_num-1 && ip_den == N_den-1 && DrawFlag == true) {
            for(int itmp_num=imatch_num+1; itmp_num < N_num; itmp_num++)
            {
                double xtmp_num, ytmp_num, xtmperr_num, ytmperr_num;

                gnum->GetPoint(itmp_num,xtmp_num,ytmp_num);
                xtmperr_num = gnum->GetErrorX(itmp_num);
                ytmperr_num = gnum->GetErrorY(itmp_num);
                cout << Form("Numerator (%d) (%g +- %g , %g +- %g) v.s Denominator (Non) Non,Non ",itmp_num,xtmp_num,xtmperr_num,ytmp_num,ytmperr_num) << endl;
            }
        }
    }

    gratio->GetYaxis()->SetTitle("ratio");
    if(Title != "") {
        gratio->SetTitle(Title);
    }
    else {
        TString t1 = gnum->GetTitle();
        TString t2 = gden->GetTitle();
//        gratio->SetTitle(Form("%s <Divide> %s",t1.Data(),t2.Data()));
    }

    if(DrawFlag == true){
        new TCanvas();
        gratio->SetMarkerStyle(8);
        gratio->SetLineWidth(3);
        gratio->Draw("AP");
    }

    return gratio;
}
TGraphErrors* GraphOperate(TGraphErrors* gnum, TGraphErrors* gden, TString FuncString, TString FuncErrString = "", TGraphErrors* gratio = NULL, TString Option = "", TString Title = "", bool DrawFlag = true, bool Estimate1percentage = true) // FuncString: TF1 formula string, hnum relate to [0], hden relate to [1]; Option: '' or 'noerr' ; FuncErrString: hnum value relate to [0], hden value relate to [1], hnum error relate to [2], hden error relate to [3];
//Example: positron fraction: FuncString:"[0]/([0]+[1])", FuncErrString:"TMath::Sqrt( ([1]*[1]*[0]*[0]) + ([0]*[0]*[1]*[1]) ) / ( ([0]+[1])*([0]+[1]) )"
{
    TString OutputFuncFormat;
    TF1* ffunc = new TF1("ffunc",FuncString,0,10);
    TF1* ferrfunc;

    OutputFuncFormat = "Func: "+FuncString;
    if(FuncErrString == "")
    {
        ferrfunc = NULL;
        OutputFuncFormat += " ; FuncErr: Set0";
    }
    else
    {
        ferrfunc = new TF1("ferrfunc",FuncErrString,0,10);
        OutputFuncFormat += " ; FuncErr: "+FuncErrString;
    }

    OutputFuncFormat.ReplaceAll("[0]","val0");
    OutputFuncFormat.ReplaceAll("[1]","val1");
    OutputFuncFormat.ReplaceAll("[2]","err0");
    OutputFuncFormat.ReplaceAll("[3]","err1");

    if(!gratio) {gratio = (TGraphErrors*)gden->Clone("gratio");}
    int N_num = gnum->GetN();
    int N_den = gden->GetN();
    int N_ratio = gratio->GetN();

    for(int ip=0; ip < N_ratio; ip++){
        gratio->RemovePoint(0);
    }

    int imatch_num = -1;
    double x_den = 0; double y_den = 0; double xerr_den = 0; double yerr_den = 0;
    double x_num = 0; double y_num = 0; double xerr_num = 0; double yerr_num = 0;
    double y_ratio, yerr_ratio;
    int ip_ratio = -1;

    for(int ip_den=0; ip_den < N_den; ip_den++)
    {
        gden->GetPoint(ip_den,x_den,y_den);
        xerr_den = gden->GetErrorX(ip_den);
        yerr_den = gden->GetErrorY(ip_den);

        if(x_den != x_den || y_den != y_den || xerr_den != xerr_den || yerr_den != yerr_den) continue;

        imatch_num = -1;
        for(int ip_num=0; ip_num < N_num; ip_num++)
        {
            gnum->GetPoint(ip_num,x_num,y_num);
            xerr_num = gnum->GetErrorX(ip_num);
            yerr_num = gnum->GetErrorY(ip_num);

            if(x_num != x_num || y_num != y_num || xerr_num != xerr_num || yerr_num != yerr_num) continue;

            if(x_num == x_den || ( Estimate1percentage && (x_num/x_den) > 0.99 && (x_num/x_den) < 1.01 ) ) 
            {
                imatch_num = ip_num;
                break;
            }
        }

        if(DrawFlag == true)
        {
            if (imatch_num > 0 && imatch_num <= N_num-1 && ip_den == 0) { 
                for(int itmp_num=0; itmp_num < imatch_num; itmp_num++)
                {
                    double xtmp_num, ytmp_num, xtmperr_num, ytmperr_num;

                    gnum->GetPoint(itmp_num,xtmp_num,ytmp_num);
                    xtmperr_num = gnum->GetErrorX(itmp_num);
                    ytmperr_num = gnum->GetErrorY(itmp_num);
                    cout << Form("Numerator (%d) (%g +- %g , %g +- %g) v.s Denominator (Non) Non,Non ; %s",itmp_num,xtmp_num,xtmperr_num,ytmp_num,ytmperr_num,OutputFuncFormat.Data()) << endl;
                }
            }
        }
        if(imatch_num >= 0)
        {
            ip_ratio++;
            ////y_ratio = y_num/y_den; yerr_ratio = y_ratio*TMath::Sqrt((yerr_num*yerr_num)/(y_num*y_num) + (yerr_den*yerr_den)/(y_den*y_den));

            ffunc->SetParameters(y_num,y_den); y_ratio = ffunc->Eval(1); // ffunc is a const value function, don't change with x;
            if(ferrfunc)
            {
                ferrfunc->SetParameters(y_num,y_den,yerr_num,yerr_den);
                yerr_ratio = ferrfunc->Eval(1);
            }
            else
            {
                yerr_ratio = 0;
            }

            if(y_num == 0 || y_den == 0) {y_ratio = 0; yerr_ratio = 0;}

            gratio->SetPoint(ip_ratio,x_den,y_ratio);
            gratio->SetPointError(ip_ratio,0,yerr_ratio);
            if(Option.Contains("noerr")) gratio->SetPointError(ip_ratio,0,0);

            if(DrawFlag) cout << Form("Numerator (%d) (%g +- %g , %g +- %g) v.s Denominator (%d) (%g +- %g , %g +- %g) ; %s",imatch_num,x_num,xerr_num,y_num,yerr_num,ip_den,x_den,xerr_den,y_den,yerr_den,OutputFuncFormat.Data()) << endl;
        }

        if(imatch_num == -1 && DrawFlag == true) {
            cout << Form("Numerator (Non) Non,Non v.s Denominator (%d) (%g +- %g , %g +- %g) ; %s ",ip_den,x_den,xerr_den,y_den,yerr_den, OutputFuncFormat.Data()) << endl;
        }

        if (imatch_num >= 0 && imatch_num < N_num-1 && ip_den == N_den-1 && DrawFlag == true) {
            for(int itmp_num=imatch_num+1; itmp_num < N_num; itmp_num++)
            {
                double xtmp_num, ytmp_num, xtmperr_num, ytmperr_num;

                gnum->GetPoint(itmp_num,xtmp_num,ytmp_num);
                xtmperr_num = gnum->GetErrorX(itmp_num);
                ytmperr_num = gnum->GetErrorY(itmp_num);
                cout << Form("Numerator (%d) (%g +- %g , %g +- %g) v.s Denominator (Non) Non,Non ; %s",itmp_num,xtmp_num,xtmperr_num,ytmp_num,ytmperr_num, OutputFuncFormat.Data()) << endl;
            }
        }
    }

    gratio->GetYaxis()->SetTitle("ratio");
    if(Title != "") {
        gratio->SetTitle(Title);
    }
    else {
        TString t1 = gnum->GetTitle();
        TString t2 = gden->GetTitle();
        gratio->SetTitle(Form("%s <Divide> %s",t1.Data(),t2.Data()));
    }

    if(DrawFlag == true){
        new TCanvas();
        gratio->SetMarkerStyle(8);
        gratio->SetLineWidth(3);
        gratio->Draw("AP");
    }

    return gratio;
}
/*TGraph* GraphSetLog(TGraph* graw, TString Option="|Logx|Logy|NewCreate|")//Option: |Logx| |Logy| |NewCreate| |Raw|
{
    bool FlagErr = false;
    TGraph* gaim = NULL;
    double* xerr = NULL; double* yerr = NULL;
    int N = graw->GetN();
    double* x = graw->GetX();
    double* y = graw->GetY();

    xerr = graw->GetEX();
    yerr = graw->GetEY();

    if(xerr != NULL && yerr != NULL) FlagErr = true;

    if(Option.Contains("|Raw|"))
    {
        gaim = graw;
    }
//    if(Option.Contains("|NewCreate|"))
    else
    {
        gaim = new TGraphErrors();
        TString Name = graw->GetName();
        gaim->SetName(Name+"_Log");
    }

    for(int ip=0; ip < N; ip++)
    {
        double xaim = NANNumber; double yaim = NANNumber; double xerraim = NANNumber; double yerraim = NANNumber;

        if(Option.Contains("|Logx|"))
        {
            xaim = TMath::Log(x[ip]);

            if(FlagErr) {xerraim = xerr[ip]/x[ip]*xaim;}

            if(xaim != xaim) {xaim = NANNumber; if(FlagErr){xerraim = NANNumber;} }
        }
        else
        {
            xaim = x[ip]; if(FlagErr){xerraim = xerr[ip];}
        }

        if(Option.Contains("|Logy|"))
        {
            yaim = TMath::Log(y[ip]);

            if(FlagErr){yerraim = yerr[ip]/y[ip]*yaim;}

            if(yaim != yaim) {yaim = NANNumber; if(FlagErr){yerraim = NANNumber;} }
        }
        else
        {
            yaim = y[ip]; if(FlagErr){yerraim = yerr[ip];}
        }

        gaim->SetPoint(ip,xaim,yaim);
        if(FlagErr) ((TGraphErrors*)gaim)->SetPointError(ip,xerraim,yerraim);
    }

    return gaim;
}*/

TGraph* GraphReCount(TGraph* graw, TString XCountFunc = "", TString YCountFunc = "", TString Option="|NewCreate|", TString SPSetName = "")//Option: |NewCreate| |Raw|
{
    bool FlagErr = false;
    TGraph* gaim = NULL;
    double* xerr = NULL; double* yerr = NULL;
    TF1* yfunc = NULL; TF1* xfunc = NULL;

    if(YCountFunc != "") yfunc = new TF1("yfunc",YCountFunc,0,1);
    if(XCountFunc != "") xfunc = new TF1("xfunc",XCountFunc,0,1); 

    int N = graw->GetN();
    double* x = graw->GetX();
    double* y = graw->GetY();

    xerr = graw->GetEX();
    yerr = graw->GetEY();

    if(xerr != NULL && yerr != NULL) FlagErr = true;

    if(Option.Contains("|Raw|"))
    {
        gaim = graw;
    }
//    if(Option.Contains("|NewCreate|"))
    else
    {
        gaim = new TGraphErrors();
        if(SPSetName == "")
        {
            TString Name = graw->GetName();
            gaim->SetName(Name+"_ReCount");
        }
        else{
            gaim->SetName(SPSetName);
        }
    }

    for(int ip=0; ip < N; ip++)
    {
        double xaim = NANNumber; double yaim = NANNumber; double xerraim = NANNumber; double yerraim = NANNumber;

        if(XCountFunc != "")
        {
            if(XCountFunc.Contains("[0]")) {xfunc->SetParameter(0,y[ip]);}
            xaim = xfunc->Eval(x[ip]);

            if(FlagErr) {xerraim = xerr[ip]/x[ip]*xaim;}

            if(xaim != xaim) {xaim = NANNumber; if(FlagErr){xerraim = NANNumber;} }
        }
        else
        {
            xaim = x[ip]; if(FlagErr){xerraim = xerr[ip];}
        }

        if(YCountFunc != "")
        {
            if(YCountFunc.Contains("[0]")) {yfunc->SetParameter(0,x[ip]);}
            yaim = yfunc->Eval(y[ip]);

            if(FlagErr){yerraim = yerr[ip]/y[ip]*yaim;}

            if(yaim != yaim) {yaim = NANNumber; if(FlagErr){yerraim = NANNumber;} }
        }
        else
        {
            yaim = y[ip]; if(FlagErr){yerraim = yerr[ip];}
        }

        gaim->SetPoint(ip,xaim,yaim);
        if(FlagErr) ((TGraphErrors*)gaim)->SetPointError(ip,xerraim,yerraim);
    }

    return gaim;
}

void CanvasPrint(TCanvas* cc, TString OutRawDir, TString DrawTxt = "", TString OutKind = "all")// ep: OutDir: basic/July/test , OutKind: all/ps/root/gif
{
    ////TString basicdir = "/publicfs/cms/data/Publics/zhangcheng/WeekReportPlots";
    TString basicdir = "/afs/cern.ch/user/c/czhang/work/private/Jobfit/WeekReportPlots";
    //TString KindName[]={"root","ps","gif"};
    TString OutName="";
    TString OutDir="";

    TObjArray* InfoSeperate = OutRawDir.Tokenize("/");
    int NInfo = InfoSeperate->GetEntries();
    
    for(int i=0; i < NInfo-1; i++)
    {
        OutDir += ((TString)((TObjString*)InfoSeperate->At(i))->String())+"/";
    }

    OutName = ((TString)((TObjString*)InfoSeperate->At(NInfo-1))->String());

    if(OutDir.Contains("basic"))
    {
        OutDir.ReplaceAll("basic",basicdir);
    }

    system(Form("source /afs/cern.ch/user/c/czhang/useful/BashFunction.sh ; mkdir_check %s",OutDir.Data()));

    if(DrawTxt != "")
    {
        ofstream fout; fout.open(Form("%s/%s.txt",OutDir.Data(),OutName.Data()));
        fout << "The RootCommand Line:  " << endl;

        if(DrawTxt == "<opentxt>")
        {
            fout.close();
            system(Form("gvim %s/%s.txt",OutDir.Data(),OutName.Data()));
        }
        else if(DrawTxt == "<functxt>")
        {
            ifstream fin; fin.open(TxtFunc_Dir+TxtFunc_Name);
            fout << fin.rdbuf();
            fout.close();
        }
        else
        {
            fout << DrawTxt << endl;
            fout.close();
        }
    }

    if(OutKind == "all")
    {
        cc->Print(Form("%s/%s.root",OutDir.Data(),OutName.Data()));
        cc->Print(Form("%s/%s.ps",OutDir.Data(),OutName.Data()));
        cc->Print(Form("%s/%s.gif",OutDir.Data(),OutName.Data()));
    }
    else
    {
        cc->Print(Form("%s/%s.%s",OutDir.Data(),OutName.Data(),OutKind.Data()));
    }
}

TH1* GetErrorHist(TH1* hraw, TH1* herror = NULL,TString Title = "", bool DrawFlag = true)
{
    int nhb = hraw->GetNbinsX();
    if(!herror)
    {
        herror = (TH1*)hraw->Clone("herror");
    }
    if(Title != "") {
        herror->SetTitle(Title);
    }

    for(int ihb=1; ihb <= nhb; ihb++)
    {
        //double tmpx;
        double tmpy, tmpyerr;
        //tmpx = hraw->GetBinCenter(ihb);
        tmpy = hraw->GetBinContent(ihb);
        tmpyerr = hraw->GetBinError(ihb);

        if(tmpy != 0) tmpyerr /= tmpy;
        else tmpyerr = -9999;

        herror->SetBinContent(ihb,tmpyerr);
        herror->SetBinError(ihb,0);

    }

    if(DrawFlag)
    {
        for(int ihb=1; ihb <= nhb; ihb++)
        {
            double tmpxlow_raw, tmpxhigh_raw, tmpy_raw, tmpyerr_raw;
            double tmpxlow_error, tmpxhigh_error, tmpy_error, tmpyerr_error;

            tmpxlow_raw = hraw->GetBinLowEdge(ihb); tmpxhigh_raw = hraw->GetBinWidth(ihb)+tmpxlow_raw;
            tmpy_raw = hraw->GetBinContent(ihb); tmpyerr_raw = hraw->GetBinError(ihb);
            tmpxlow_error = herror->GetBinLowEdge(ihb); tmpxhigh_error = herror->GetBinWidth(ihb)+tmpxlow_error;
            tmpy_error = herror->GetBinContent(ihb); tmpyerr_error = herror->GetBinError(ihb);
            cout << Form("RawHist (%d) [%g,%g] (%g #pm %g) v.s ErrorHist (%d) [%g,%g] (%g #pm %g)",ihb,tmpxlow_raw,tmpxhigh_raw,tmpy_raw,tmpyerr_raw,ihb,tmpxlow_error,tmpxhigh_error,tmpy_error,tmpyerr_error) << endl;
        }

        new TCanvas();
        herror->SetMarkerStyle(8);
        herror->SetLineWidth(3);
        herror->Draw();
    }

    return herror;
}
int SetHistErr(TH1* haim, TH1* href)
{
    int Naim = haim->GetNbinsX();
    int Nref = href->GetNbinsX();

    if(Nref < Naim) { cout << "HRef are less than HAim!!! Can't Deal" << endl; return -1;}

    for(int ih=1; ih <= Naim; ih++)
    {
        double xaim = haim->GetBinCenter(ih);
        double yaim = haim->GetBinContent(ih);
        int setih = href->FindBin(xaim);

        double xref = href->GetBinCenter(setih);
        if(xaim != xref) { cout << "HRef and HAim have different binning !!! Can't Deal" << endl; return -2;}

        double yvalue = href->GetBinContent(setih); double yerr = href->GetBinError(setih);
        double relateerr = yerr/yvalue;

        double yaimerr = yaim*relateerr;

        haim->SetBinError(ih,yaimerr);
    }

    return 0;
}
TGraphErrors* GetErrorGraph(TGraphErrors* graw, TGraphErrors* gerror = NULL,TString Title = "", bool DrawFlag = true)
{
    int NPoint = graw->GetN();
    if(!gerror) {gerror = (TGraphErrors*)graw->Clone("gerror");}
    if(Title != "") {
        gerror->SetTitle(Title);
    }

    for(int ip=0; ip < NPoint; ip++){
        gerror->RemovePoint(0);
    }

    int icount = 0;
    for(int ip=0; ip < NPoint; ip++)
    {
        double tmpx, tmpy, tmpyerr;
        graw->GetPoint(ip,tmpx,tmpy);
        tmpyerr = graw->GetErrorY(ip);

        if(tmpx != tmpx || tmpy != tmpy || tmpyerr != tmpyerr) continue;

        if(tmpy != 0) tmpyerr /= tmpy;
        else tmpyerr = -9999;

        gerror->SetPoint(icount,tmpx,tmpyerr);
        gerror->SetPointError(icount,0,0);

        icount++;
    }

    if(DrawFlag)
    {
        for(int ip=0; ip < NPoint; ip++)
        {
            double tmpx_raw, tmpy_raw, tmpyerr_raw;
            double tmpx_error, tmpy_error, tmpyerr_error;

            graw->GetPoint(ip,tmpx_raw,tmpy_raw);
            tmpyerr_raw = graw->GetErrorY(ip);
            gerror->GetPoint(ip,tmpx_error,tmpy_error);
            tmpyerr_error = gerror->GetErrorY(ip);

            cout << Form("RawGraph (%d) x: %g (%g #pm %g) v.s ErrorGraph (%d) x: %g (%g #pm %g)",ip, tmpx_raw, tmpy_raw,tmpyerr_raw, ip, tmpx_error, tmpy_error,tmpyerr_error) << endl;
        }
        new TCanvas();
        gerror->SetMarkerStyle(8);
        gerror->SetLineWidth(3);
        gerror->Draw("AP");
    }
    return gerror;
}
TH1F* GetErrorGraphAsHist(TGraphErrors* graw, TH1F* htmp, TGraphErrors* gerror = NULL,TString Title = "", bool DrawFlag = true)
{
    TGraphErrors* gerr = GetErrorGraph(graw,gerror,Title,DrawFlag);
    TH1F* houtput = (TH1F*)htmp->Clone("herr");
    
    int Nhist = houtput->GetNbinsX();
    for(int ihist=1; ihist <= Nhist; ihist++)
    {
        houtput->SetBinContent(ihist,-999);
        houtput->SetBinError(ihist,0);
    }

    int NPoint = gerr->GetN();
    for(int ip=0; ip < NPoint; ip++)
    {
        double x; double y;
        gerr->GetPoint(ip,x,y);
        
        int iset = houtput->FindBin(x);
        houtput->SetBinContent(iset,y);
    }

    return houtput;
}

TGraphErrors* HistToGraph(TH1* hraw, TString Option = "CenterAsX", TString Name = "") // Option: "CenterAsX","CenterLogAsX"
{
    TGraphErrors *gaim = new TGraphErrors();

    if(Name == "")
    {
        Name = hraw->GetName();
        gaim->SetName("g_"+Name);

    }
    else 
    {
        gaim->SetName(Name);
    }

    int Nhist = hraw->GetNbinsX();
    for(int ib=1;ib<=Nhist;ib++){
        double xv = 0;
        if(Option.Contains("CenterLogAsX")){
            xv=hraw->GetXaxis()->GetBinCenterLog(ib);//sqrt(lb*hb)
        }
        else if(Option.Contains("CenterAsX")) {
            xv=hraw->GetXaxis()->GetBinCenter(ib);//(lb+hb)/2
        }
        double yv=hraw->GetBinContent(ib);
        double ey=hraw->GetBinError(ib);
        int ip=gaim->GetN();
        gaim->SetPoint(ip,xv,yv);
        gaim->SetPointError(ip,0,ey);
    }
    return gaim;
}

vector<TString> GetSeperateTString(TString SRaw, TString SMark = ",")
{
    vector<TString> SOut;
    TObjArray* InfoAll = SRaw.Tokenize(SMark);
    int NString = InfoAll->GetEntries();
    TString tmps;
    for(int is=0; is < NString; is++)
    {
        tmps = (TString)((TObjString*)InfoAll->At(is))->String();
        SOut.push_back(tmps);
    }

    return SOut;
}
int GetSeperateTString_PointTrans(TString SRaw, vector<TString>* SOut, TString SMark = ",")
{
////    vector<TString> SOut;

    cout << "test" << endl;
    cout << "In String vector, SOut address: " << SOut << endl;
    cout << "SRaw: " << SRaw << endl;

    TObjArray* InfoAll = SRaw.Tokenize(SMark);
    int NString = InfoAll->GetEntries();
    cout << "InfoAll Num: " << NString << endl;
    TString tmps;
    for(int is=0; is < NString; is++)
    {
        cout << "In Loop: " << is << endl;
        tmps = (TString)((TObjString*)InfoAll->At(is))->String();
        cout << "tmps: " << tmps << endl;
        SOut->push_back(tmps);
        cout << "test again" << endl;
        cout << "SOut value: " << (*SOut)[is] << endl;
    }

    ////return SOut;
    return 1;
}
double* OptionStringSPSet_RangeGet(TString Option)
{
    double range[2];
    vector<TString> string_stp1 = GetSeperateTString(Option, "()");
    vector<TString> string_stp2 = GetSeperateTString(string_stp1[1], ",");

    range[0] = string_stp2[0].Atof(); range[1] = string_stp2[1].Atof();
    
    return &(range[0]);
}

double* GetHistArray_LogBins(int Narray,double genrange[2], double setlow = -999, double sethigh = -999)
{
    double* array = new double[Narray+1];
    double xlow = genrange[0]; double xhigh = genrange[1];

    if(setlow != -999 && setlow > genrange[0]) { xlow = setlow;}
    if(sethigh != -999 && sethigh < genrange[1]) { xhigh = sethigh;}

    for(int il=0; il <= Narray; il++)
    {
        array[il] = TMath::Power(xhigh/xlow,il*1.0/(Narray))*xlow;
    }

    return array;
}
int VectorMinMax(vector<double> raw, TString Option, double& minvalue, double& maxvalue ) // Option: Abs , NoZero, ...
{
    double tmpmin = 0; double tmpmax = 0; int id_tmpmin = 0; int id_tmpmax =0;
    bool flagAbs = false; bool flagNoZero = false; bool flagPositive = false;

    if(Option.Contains("Abs")) flagAbs = true;
    if(Option.Contains("NoZero")) flagNoZero = true;
    if(Option.Contains("JustPositive")) flagPositive = true;

    int count = -1;
    for(unsigned int i=0; i < raw.size(); i++)
    {
        double tmpvalue;
        if(flagAbs) tmpvalue = TMath::Abs(raw[i]);
        else tmpvalue = raw[i];

        if(flagNoZero && tmpvalue ==0) continue;
        if(flagPositive && tmpvalue < 0) continue;

        count++;
        if(count == 0) {tmpmin = raw[i]; tmpmax = raw[i];}

        if(tmpvalue > tmpmax) {tmpmax = tmpvalue; id_tmpmax = i;}
        if(tmpvalue < tmpmin) {tmpmin = tmpvalue; id_tmpmin = i;}
        cout << Form("loop: %d, min: %g, max: %g",i,tmpmin,tmpmax) << endl; 
    }

    minvalue = tmpmin; maxvalue = tmpmax;
    cout << "IDMin: " << id_tmpmin << " and IDMax: " << id_tmpmax << endl;
    return 0;
}

vector<double> ArraytoVector(double* array, int NLength)
{
    vector<double> out;
    for(int i=0; i < NLength; i++)
    {
        out.push_back(array[i]);
    }
    return out;
}
void VectortoArray(vector<double> vet, double out[100])
{
    if(vet.size() > 100) {cout << "Error vector size out of limit(100)" << endl;}
    for(unsigned int i=0; i < 100; i++)
    {
        if(i < vet.size()) out[i] = vet[i];
        else out[i] = 0;
    }
}
vector<bool> ArraytoVector(bool* array, int NLength)
{
    vector<bool> out;
    for(int i=0; i < NLength; i++)
    {
        out.push_back(array[i]);
    }
    return out;
}
int GraphMinMax(TGraph* graw, TString Option, double& minvalue, double& maxvalue ) // Option: Abs , NoZero, ..., Range(a,b)'==> a,b is the hist binvalue range for get the maxvalue
{
//    double tmpmin; double tmpmax; int id_tmpmin; int id_tmpmax;
    bool flagAbs = false; bool flagNoZero = false; bool flagPositive = false;
    double rangelow = -999; double rangehigh = -999;
    vector<double> raw;

    if(Option.Contains("Abs")) flagAbs = true;
    if(Option.Contains("NoZero")) flagNoZero = true;
    if(Option.Contains("JustPositive")) flagPositive = true;

    cout << "Flags: " << flagAbs << " and " << flagNoZero << " and " << flagPositive << endl;

    if(Option.Contains("Range"))
    {
        vector<TString> string_stp1 = GetSeperateTString(Option, "()");
        vector<TString> string_stp2 = GetSeperateTString(string_stp1[1], ",");

        rangelow = string_stp2[0].Atof(); rangehigh = string_stp2[1].Atof();
    }

    for(int ig=0; ig < graw->GetN(); ig++)
    {
        double gbinx; double gbiny;
        graw->GetPoint(ig,gbinx,gbiny);

        if(rangelow != -999 && rangehigh != -999)
        {
            if(gbinx >= rangelow && gbinx < rangehigh)
            {
                raw.push_back(gbiny);
            }
        }
    }
/*
       Old Design, no Range option set
       int count = -1;
       double* raw_array = graw->GetY();
       int N = graw->GetN();
       vector<double> raw = ArraytoVector(raw_array,N);
*/
    int I_return = VectorMinMax(raw,Option,minvalue,maxvalue);
    return I_return;
}
int HistMinMax(TH1* hraw, TString Option, double& minvalue, double& maxvalue ) ////Option: Range(a,b), Abs, NoZero, ...
{
//    double tmpmin; double tmpmax; int id_tmpmin; int id_tmpmax;
    bool flagAbs = false; bool flagNoZero = false; bool flagPositive = false;
    double rangelow = -999; double rangehigh = -999;
    vector<double> raw;
    if(Option.Contains("Abs")) flagAbs = true;
    if(Option.Contains("NoZero")) flagNoZero = true;
    if(Option.Contains("JustPositive")) flagPositive = true;

    cout << "Flags: " << flagAbs << " and " << flagNoZero << " and " << flagPositive << endl;

    if(Option.Contains("Range"))
    {
        vector<TString> string_stp1 = GetSeperateTString(Option, "()");
        vector<TString> string_stp2 = GetSeperateTString(string_stp1[1], ",");

        rangelow = string_stp2[0].Atof(); rangehigh = string_stp2[1].Atof();
    }

    for(int ih=1; ih <= hraw->GetNbinsX(); ih++)
    {
        double hbinx; double hbiny;
        hbinx = hraw->GetBinCenter(ih); hbiny = hraw->GetBinContent(ih);

        if(rangelow != -999 && rangehigh != -999)
        {
            if(hbinx >= rangelow && hbinx < rangehigh)
            {
                raw.push_back(hbiny);
            }
        }
    }

    int I_return = VectorMinMax(raw,Option,minvalue,maxvalue);
    return I_return;
}
TH1* HistNorm(TH1* HRaw, TString normkind = "EventsNorm", TH1* HNorm = NULL, bool ChangeRawHist = true, double normvalue = 1) // 'EventsNorm','MaxValueNorm','MaxValueNormRange(a,b)'==> a,b is the hist binvalue range for get the maxvalue
{
    double valuefornorm = -1;

    if(!ChangeRawHist)
    {
        if (!HNorm)
        {
            TString Name = HRaw->GetName();
            HNorm = (TH1*)HRaw->Clone(Name);
        }
    }
    else
    {
        HNorm = HRaw;
    }
    if(normkind == "EventsNorm")
    {
        valuefornorm = HNorm->Integral();
        HNorm->Scale(normvalue/valuefornorm);
    }
    else if(normkind.Contains("MaxValueNorm"))
    {
        if(normkind == "MaxValueNorm")
        {
            valuefornorm = HNorm->GetMaximum();
            HNorm->Scale(normvalue/valuefornorm);
        }
        else if(normkind.Contains("MaxValueNormRange"))
        {
            double tmpmin; double tmpmax;

            HistMinMax(HNorm, normkind, tmpmin, tmpmax );
            valuefornorm = tmpmax;
            HNorm->Scale(normvalue/valuefornorm);
        }
    }

    return HNorm;
}

void SetHistStyle(TH1* h1draw, int colorset=1, int markerstyle = 8, int linestyle = 1, bool FlagDraw = false, bool AddPadSet = true)
{
    ////if(FlagDraw) { new TCanvas(); }
    if(AddPadSet && FlagDraw) { new TCanvas(); }

    if(AddPadSet && FlagDraw) {
        gPad->SetLeftMargin(.15);
        gPad->SetBottomMargin(.15);
    }

    h1draw->GetXaxis()->SetLabelSize(.06);
    h1draw->GetYaxis()->SetLabelSize(.06);
    h1draw->GetXaxis()->SetTitleSize(.07);
    h1draw->GetYaxis()->SetTitleSize(.07);
    h1draw->GetXaxis()->CenterTitle();
    h1draw->GetYaxis()->CenterTitle();
    h1draw->GetXaxis()->SetNdivisions(505);
    h1draw->GetYaxis()->SetNdivisions(505);
    h1draw->GetXaxis()->SetTitleOffset(.9);
    h1draw->GetYaxis()->SetTitleOffset(.9);

    h1draw->SetLineWidth(3);
    h1draw->SetLineColor(colorset);
    h1draw->SetLineStyle(linestyle);
    h1draw->SetMarkerColor(colorset);
    h1draw->SetMarkerStyle(markerstyle);

    if(FlagDraw) { h1draw->Draw(); }
}
void SetGraphStyle(TGraph* h1draw, int colorset = 1, int markerstyle = 8, int linestyle = 1, bool FlagDraw = false, bool AddPadSet = true)
{
    ////if(FlagDraw) { new TCanvas(); }
    if(AddPadSet && FlagDraw) { new TCanvas(); }
    if(AddPadSet && FlagDraw) {
    gPad->SetLeftMargin(.15);
    gPad->SetBottomMargin(.15);
    }

   h1draw->GetXaxis()->SetLabelSize(.06);
   h1draw->GetYaxis()->SetLabelSize(.06);
   h1draw->GetXaxis()->SetTitleSize(.07);
   h1draw->GetYaxis()->SetTitleSize(.07);
   h1draw->GetXaxis()->CenterTitle();
   h1draw->GetYaxis()->CenterTitle();
   h1draw->GetXaxis()->SetNdivisions(505);
   h1draw->GetYaxis()->SetNdivisions(505);
   h1draw->GetXaxis()->SetTitleOffset(.9);
   h1draw->GetYaxis()->SetTitleOffset(.9);

   h1draw->SetLineWidth(3);
   h1draw->SetLineColor(colorset);
   h1draw->SetLineStyle(linestyle);
   h1draw->SetMarkerColor(colorset);
   h1draw->SetMarkerStyle(markerstyle);

   if(FlagDraw) { h1draw->Draw("AP"); }
}
TH1* HistGetFromND(TH1* hHist,TString OutName = "",TString ProjectDim = "x",double SetDim2Value = -1) //ProjectDim: 'x' or 'y'
{
    TH1D* hdim1d = NULL;

    int NDim = hHist->GetDimension();
    if(OutName == "") OutName = hHist->GetName();

    if(NDim == 1)
    {
        hdim1d = (TH1D*)hHist->Clone(OutName);
    }
    else if(NDim == 2)
    {
        int ProjectID = -1; TH1D* htmp = NULL; 
        if(ProjectDim == "x") ProjectID = 0;
        else if(ProjectDim == "y") ProjectID = 1;

        if(ProjectID == 0) htmp = ((TH2F*)hHist)->ProjectionX("htmp");
        else if(ProjectID == 1) htmp = ((TH2F*)hHist)->ProjectionY("htmp");

        hdim1d = (TH1D*)htmp->Clone(OutName);

        if(SetDim2Value != -1)
        {
            int nhb = hdim1d->GetNbinsX();
            for(int ihb=1; ihb <= nhb; ihb++)
            {
                double binvalue; double bincenter; double bingid = 0;
                bincenter = hdim1d->GetBinCenter(ihb);

                if(ProjectID == 0) bingid = hHist->FindBin(bincenter,SetDim2Value);
                if(ProjectID == 1) bingid = hHist->FindBin(SetDim2Value,bincenter);

                binvalue = hHist->GetBinContent(bingid);
                hdim1d->SetBinContent(ihb,binvalue);
            }
        }
    }
    return hdim1d;
}
TH1* HistGetFromND(TString InName,TString OutName = "",TString ProjectDim = "x",double SetDim2Value = -1) //ProjectDim: 'x' or 'y'
{
    TH1D* hdim1d;

    if(OutName == "") OutName = InName;
    TH1* hHist = (TH1*)(gROOT->FindObject(InName));
    
//    cout << "GetFunciton Name: " << InName << endl;
//    if(hHist) cout << "hHist has get" << endl;
//    else cout << "hHist NOT get" << endl;

    if(!hHist) return hHist;

    hdim1d = (TH1D*)HistGetFromND(hHist,OutName,ProjectDim,SetDim2Value);
    return hdim1d;
}
double HistNDChi2Count(TH1* Hist_data,TH1* Hist_fit)
{
    double totalvalue = 0;
    
    int Nhbx_data = Hist_data->GetNbinsX();
    int Nhby_data = Hist_data->GetNbinsY();
    int Nhbz_data = Hist_data->GetNbinsZ();
    int Nhbx_fit = Hist_fit->GetNbinsX();
    int Nhby_fit = Hist_fit->GetNbinsY();
    int Nhbz_fit = Hist_fit->GetNbinsZ();

    if(Nhbx_data != Nhbx_fit || Nhby_data != Nhby_fit || Nhbz_data != Nhbz_fit) 
    {
        cout << "===!!! In Func 'Hist1DChi2Count': Nhb_data != Nhb_fit !!!===" << endl;
        return -999;
    }

    for(int ihbx=1; ihbx <= Nhbx_data; ihbx++)
    {
        for(int ihby=1; ihby <= Nhby_data; ihby++)
        {
            for(int ihbz=1; ihbz <= Nhbz_data; ihbz++)
            {
                double countvalue;
                double value_data; double value_fit; double err_data;
                value_data = Hist_data->GetBinContent(ihbx,ihby,ihbz);
                err_data = Hist_data->GetBinError(ihbx,ihby,ihbz);
                value_fit = Hist_fit->GetBinContent(ihbx,ihby,ihbz);

                if(value_data > 0)
                {
                    cout << Form("For ih(%d,%d,%d), the value_data is: %g, the value_fit is: %g, the err_data is: %g",ihbx,ihby,ihbz,value_data,value_fit,err_data) << endl;

                    countvalue = (value_data-value_fit)*(value_data-value_fit)/(err_data*err_data);
                    totalvalue += countvalue;

                    if(err_data <= 0) {
                        countvalue = -999; totalvalue = -99999;
                    }
                }
            }
        }
    }

    totalvalue /= (Nhbx_data*Nhby_data*Nhbz_data-1);
    if(totalvalue <= 0) totalvalue = -999;

    return totalvalue;
}

void readDataFrom(TChain *fChain,TString fname)
{
    int StartRunNo = 0;
    int StartEvtNo = 0;

    TString prefix = "root://castorpublic.cern.ch//";
    TString prefix1 = "root://eosams.cern.ch//";
    TString server = "?svcClass=amsuser";
    if( fname.Contains( "root" ) ){
        TString TmpBase = gSystem->BaseName(fname);
        StartRunNo = TString(TmpBase(0,10)).Atoi();
        StartEvtNo = TString(TmpBase(11,8)).Atoi();
        //cout<<StartRunNo<<"."<<StartEvtNo<<".root"<<endl;
        if( fname.Contains("castor") ){
            fname = prefix + fname + server;
        }
        else if (fname.Contains("eos")){
            fname = prefix1 + fname;
        }
        fChain->Add( fname );
    }
    else{
        ifstream fin;
        fin.open( fname );
        TString fnm;
        while(1){
            fin >> fnm;
            if( fin.eof() ) break;
            if( fnm.Contains("#") ) continue;
            if( fnm.Contains("root") ){
                TString TmpBase = gSystem->BaseName(fname);
                StartRunNo = TString(TmpBase(0,10)).Atoi();
                StartEvtNo = TString(TmpBase(11,8)).Atoi();
                if( fnm.Contains("castor") ){
                    fnm = prefix + fnm + server;
                }
                else if (fnm.Contains("eos")){
                    fnm = prefix1 + fnm;
                }
            }
            fChain->Add( fnm );

        }
    }
    cout << "totally " << fChain->GetEntries() << " data events " << endl;
    cout << "StartRunNo: " << StartRunNo << " and StartEvtNo: " << StartEvtNo << endl;
}
void SetOfficalAMSFileRead(TString Opt = "offical") // offical / local
{
    system("source ~/useful/vdev_icc.sh");
    if(Opt == "offical")
    {
    gSystem->Load("/cvmfs/ams.cern.ch/Offline/vdev/lib/linuxx8664icc5.34/libntuple_slc6_PG_dynamic.so");
    }
    else if(Opt == "local")
    {
    gSystem->Load("/afs/cern.ch/user/c/czhang/useful/WangXueQiang_ntuple_slc6_PG.so");
    }
}

vector<TFile*> GetFileFromChain(TChain* t)
{
    vector<TFile*> files;
    TObjArray* fileElements = t->GetListOfFiles();
    TIter next(fileElements);
    TChainElement* chE1 = 0;
    while (( chE1 = (TChainElement*)next() )) {
        files.push_back( new TFile(chE1->GetTitle()) );
    }
    return files;
}

vector<TH1F*> RootDrawHist(vector<TChain*>* t, vector<TString>* DrawValueSet, TString BasicCut, vector<TString>* AddCut, vector<TString>* Name_t = NULL, TString NormStyle = "EventsNorm") //'EventsNorm','MaxValueNorm',""
{
    int Color[] = {2,4,1,6,8,9,7,5,3};
    const int NMax = sizeof(Color)/sizeof(Color[0]);
    vector<TH1F*> h; vector<TString> String_Draw;
    TString LegendName[NMax];
    ofstream fout;
    ////fout.open("./RootFunc_RootDrawHist_record.txt");
    fout.open(TxtFunc_Dir+TxtFunc_Name);

    int Nf; int NCut; int NDrawValue; int Ntotal;
 //   int i_Nf = 0; int i_NCut = 0;
    vector<TString> DrawValue; vector<TString> HistSet;

    NDrawValue = (*DrawValueSet).size();

    for(int idraw=0; idraw < NDrawValue; idraw++)
    {
        String_Draw = GetSeperateTString((*DrawValueSet)[idraw],";");
        DrawValue[idraw] = String_Draw[0];
        HistSet[idraw] = String_Draw[1];
    }

/*
    while (t[i_Nf] != NULL) { i_Nf++; }
    while (AddCut[i_NCut] != "NULL") { i_NCut++; }
    Nf = i_Nf; NCut = i_NCut;
*/
    Nf = (*t).size(); NCut = (*AddCut).size();

    Ntotal = Nf*NCut*NDrawValue;
    cout << "NChain: " << Nf << " and NDrawValue: " << NDrawValue << " and NCut: " << NCut << endl;
    if(Ntotal > NMax) cout << "Error for RootDrawHist, InputN: " << Ntotal << " is larger than Max: " << NMax << endl;
    cout << "BasicCut is: " << BasicCut << endl;

    cout << "DrawValues are: ";
    for(int i=0; i < NDrawValue; i++){ cout << DrawValue[i] << " , ";};
    cout << endl;

    cout << "HistSets are: ";
    for(int i=0; i < NDrawValue; i++){ cout << HistSet[i] << " , ";};
    cout << endl;
    cout << "Add Cuts are: ";
    for(int i=0; i < NCut; i++){ cout << (*AddCut)[i] << " , "; }; 
    cout << endl;

    //======================== For record output ==================
    for(int i=0; i < Nf; i++)
    {
        fout << Form("TChain* t%d = new TChain(\"%s\"); t%d->Add(\" \");",i,(*t)[i]->GetName(),i) << endl; 
    }
    fout << Form("vector<TChain*> tall = {");
    for(int i=0; i < Nf; i++)
    {
        if(i != 0) fout << ",";
        fout << Form(" t%d ",i); 
    }
    fout << "};" << endl;
    fout << "TString DrawValue; TString BasicCut; vector<TH1F*> hall;" << endl;
    fout << Form("BasicCut = \"%s\"",BasicCut.Data()) << endl;
    fout << Form("vector<TString> DrawValue = {");
    for(int i=0; i < NCut; i++)
    {
        if(i != 0) fout << ",";
        fout << Form(" \"%s\" ",(*DrawValueSet)[i].Data()) << endl;
    }
    fout << "};" << endl;

    fout << Form("vector<TString> AddCut = {");
    for(int i=0; i < NCut; i++)
    {
        if(i != 0) fout << ",";
        fout << Form(" \"%s\" ",(*AddCut)[i].Data()); 
    }
    fout << "};" << endl;
    if(Name_t == NULL)
    {
        fout << Form("hall = RootDrawHist(&tall,&DrawValue,BasicCut,&AddCut,NULL,\"%s\")",NormStyle.Data()) << endl;
    }
    else{
        fout << Form("vector<TString> Name_t = {");
        for(int i=0; i < Nf; i++)
        {
            if(i != 0) fout << ",";
            fout << Form(" \"%s\" ",(*Name_t)[i].Data()); 
        }
        fout << "};" << endl;
        fout << Form("hall = RootDrawHist(&tall,DrawValue,BasicCut,&AddCut,&Name_t,\"%s\")",NormStyle.Data()) << endl;
    }
    fout << "The root lists: " << endl;
    for(int i=0; i < Nf; i++)
    {
        vector<TFile*> fadd;
        fout << Form(" t%d: ",i) << endl; 
        fadd = GetFileFromChain((*t)[i]);
        int Nfile = fadd.size();
        TString tmps;
        for(int is=0; is < Nfile; is++)
        {
            tmps = (TString)(fadd[is]->GetName());
            fout << tmps << "\t";
            delete fadd[is];
        }
        fout << endl;
        fadd.clear();
    }
    fout.close();

    //=============================================================

    TCanvas* call = new TCanvas("call","call",1600,1200);
    int NCanvas_Sepe = (int)((TMath::Sqrt(Ntotal)*2+1.9999)/2);
    call->Divide(NCanvas_Sepe,NCanvas_Sepe);

    int i_Canvas = 0;
    for(int idf=0; idf < Nf; idf++)
    {
        for(int idv=0; idv < NDrawValue; idv++)
        {
            for(int idc=0; idc < NCut; idc++)
            {
                call->cd(i_Canvas+1);
                cout << "Pass TChain: " << idf << " and DrawValue: " << idv << " and Cut: " << idc << endl;
                (*t)[idf]->Draw(Form("%s >> h%d%s",DrawValue[idv].Data(),i_Canvas,HistSet[idv].Data()),Form("%s && %s",BasicCut.Data(),(*AddCut)[idc].Data())); 
                h.push_back((TH1F*)gROOT->FindObject(Form("h%d",i_Canvas)));
                if(Nf > 1) { 
                    if(Name_t != NULL)  LegendName[i_Canvas] = (*Name_t)[idf]+" "+DrawValue[idv]+" "+(*AddCut)[idc]; 
                    else if(Name_t == NULL) LegendName[i_Canvas] = Form("F%d ",idf)+DrawValue[idv]+(*AddCut)[idc];
                }
                else { LegendName[i_Canvas] = DrawValue[idv]+(*AddCut)[idc]; }

                i_Canvas++; 
            }
        }
    }
    call->cd(0); call->Draw();

    TLegend* gl = new TLegend(0.75,0.75,0.9,0.9);

    TCanvas* cc = new TCanvas("cc","cc");
    for(int id=0; id < Ntotal; id++)
    {
        SetHistStyle(h[id]);
        if(NormStyle != "") { HistNorm(h[id],NormStyle); }

        h[id]->SetLineColor(Color[id]);
        h[id]->SetLineWidth(3);
        h[id]->SetMarkerColor(Color[id]);
        h[id]->SetMarkerStyle(8);
        h[id]->GetXaxis()->SetTitle(DrawValue[0]);

        gl->AddEntry(h[id],LegendName[id],"l");

        if(id ==0) h[id]->Draw();
        else h[id]->Draw("same");
    }
    gl->SetFillColor(0);
    gl->Draw("same");
    cc->Draw();

    return h;
}

int FileStuffReplace(TString FileAim,TString FileSource, TString StuffName_Aim, TString StuffName_Source = "", TString StuffKind = "Hist1D" )
{
    TH1F* haim = NULL; TH1F* hsource = NULL;
    if(StuffName_Source == "") { StuffName_Source = StuffName_Aim; }
    TFile* faim = TFile::Open(FileAim,"update");
    TFile* fsource = TFile::Open(FileSource);

    fsource->cd();
    if(StuffKind == "Hist1D")
    {
        hsource = (TH1F*)gROOT->FindObject(StuffName_Source);
    }
    faim->Delete(StuffName_Aim);
    faim->cd();
    hsource->SetName(StuffName_Aim);
    hsource->Write();
    faim->Close();
    fsource->Close();

    faim = TFile::Open(FileAim);
    fsource = TFile::Open(FileSource);
    if(StuffKind == "Hist1D")
    {
        fsource->cd();
        hsource = (TH1F*)gROOT->FindObject(StuffName_Source);
        faim->cd();
        haim = (TH1F*)gROOT->FindObject(StuffName_Aim);
    }
    TH1F* hratio = (TH1F*)HistRatio(haim,hsource);
    faim->Close();
    fsource->Close();
    return 0;
}
void HistNameAdd(TH1F** haim, TString NameAdd)
{
    if((*haim)) 
    {
        TString NameRaw = (*haim)->GetName();
        (*haim)->SetName(NameRaw+NameAdd);
    }
}
TH1F* VariableRejectionPower(TH2F* hsig, TH2F* hbkg) // hsig,hbkg: EvsAimVariable
{
    TH1F* hresult;

    int Nh = hsig->GetNbinsX();
    hresult = (TH1F*)hsig->ProjectionX("hrej"); // Just to Get the Ebins Setting.

    for(int ih=1; ih <= Nh; ih++)
    {
        double xcenter = hresult->GetBinCenter(ih);
        TH1F* h1Dsig = (TH1F*)HistGetFromND(hsig,"tmpsig","y",xcenter);
        TH1F* h1Dbkg = (TH1F*)HistGetFromND(hbkg,"tmpsig","y",xcenter);

        int Nh_y = h1Dsig->GetNbinsX();
        double totalN = h1Dsig->Integral();
        double totalN_proton = h1Dbkg->Integral();
        //int setib = -1;
        double ratioreject;
        for(int ihy = 1; ihy <= Nh_y; ihy++)
        {
            double countN0 = h1Dsig->Integral(ihy,Nh_y);
            double countN1 = h1Dsig->Integral(ihy+1,Nh_y);
            if(countN0*1.0/totalN >= 0.9 && countN1*1.0/totalN <= 0.9) 
            {
                double protoncount = h1Dbkg->Integral(ihy,Nh_y);
                cout << Form("Set CutValue: %g, NSigreal: %g, NSigAll: %g, NProtonreal: %g, NProtonAll: %g",h1Dsig->GetBinCenter(ihy),countN0,totalN,protoncount,totalN_proton) << endl;
                if(protoncount != 0)
                {
                    ratioreject = (totalN_proton*1.0)/protoncount;
                }
                else
                {
                    ratioreject = 0;
                }
                hresult->SetBinContent(ih,ratioreject);

                break;
            }
        }
    }

    return hresult;
}
TH1F* CombineHist(TH1F* hlow, TH1F* hhigh, double SetXValue, TString SetName = "") // Need Add Option based on HighPart, or based on LowPart; Currently the code is based on HighPart.
{
    if(SetName == "") SetName = hhigh->GetName();
    TH1F* hcombine = (TH1F*)hhigh->Clone(SetName);
    TString HistTitle = hhigh->GetTitle();
    hcombine->SetTitle(HistTitle);

    double xhigh, xlow, yhigh, ylow, ylowerr, yhigherr;
    int Nhist = hhigh->GetNbinsX();
    for(int ih=1; ih <= Nhist; ih++)
    {
        xhigh = hhigh->GetBinCenter(ih);
        yhigh = hhigh->GetBinContent(ih);

        if(xhigh < SetXValue)
        {
            int setih = hlow->FindBin(xhigh);
            xlow = hlow->GetBinCenter(setih);
            ylow = hlow->GetBinContent(setih);
            ylowerr = hlow->GetBinError(setih);

            hcombine->SetBinContent(ih,ylow);
            hcombine->SetBinError(ih,ylowerr);
        }
    }
 
    return hcombine;
}
void HistClear(TH1F* haim, TString Option = "All") // All, JustError, ...
{
    if(Option == "All")
    {
        haim->Reset("M");
    }
    else if(Option == "JustError")
    {
        int N = haim->GetNbinsX();
        for(int i=1; i <= N; i++)
        {
            haim->SetBinError(i,0);
        }
    }
}

double* AxisBinArrayGet(TAxis* haim, int size)
{
    TArrayD* harray = (TArrayD*)(haim->GetXbins());
    double* arraypoint;
    arraypoint = harray->GetArray();
    return arraypoint;
}
double* HistBinArrayGet(TH1F* haim,int size)
{
    TArrayD* harray = (TArrayD*)(haim->GetXaxis()->GetXbins());
    double* arraypoint;
    arraypoint = harray->GetArray();
    return arraypoint;
}

vector<double> HistBinArrayGet(TH1F* haim) // Need change from 1D Hist x bins to ND Hist x/y/z all axis bins
{
    vector<double> array;
    int nhx = haim->GetNbinsX();

    double * arraypoint = HistBinArrayGet(haim,nhx);
    for(int i=0; i <= nhx; i++)
    {
        array.push_back(arraypoint[i]);
    }

    return array;
}
vector<double> AxisBinArrayGet(TAxis* haim)
{
    vector<double> array;

    int size = haim->GetNbins();
    double * arraypoint = AxisBinArrayGet(haim,size);
    for(int i=0; i <= size; i++)
    {
        array.push_back(arraypoint[i]);
    }

    return array;
}

TH1* HistPartGet(TH1* haim_raw, double LowBinValueSet, double HighBinValueSet, TString Option = "noreset", double* NDLowBinValueSet = NULL, double* NDHighBinValueSet = NULL) // "reset" or "noreset"
{
    TH1* hget = NULL; TH1* haim = NULL;
    double TotalBinValueSet[3][2]; // Max 3D, low and high

    int ndim = haim_raw->GetDimension();

    if(ndim == 1) { haim = (TH1F*)(haim_raw); }
    if(ndim == 2) { haim = (TH2F*)(haim_raw); }
    if(ndim == 3) { haim = (TH3F*)(haim_raw); }

    if(NDLowBinValueSet != NULL && NDHighBinValueSet != NULL)
    {
        if(ndim >= 1) { TotalBinValueSet[0][0] = NDLowBinValueSet[0]; TotalBinValueSet[0][1] = NDHighBinValueSet[0]; }
        if(ndim >= 2) { TotalBinValueSet[1][0] = NDLowBinValueSet[1]; TotalBinValueSet[1][1] = NDHighBinValueSet[1]; }
        if(ndim >= 3) { TotalBinValueSet[2][0] = NDLowBinValueSet[2]; TotalBinValueSet[2][1] = NDHighBinValueSet[2]; }
    }
    else{
        if(ndim >= 1) { TotalBinValueSet[0][0] = LowBinValueSet; TotalBinValueSet[0][1] = HighBinValueSet; }
        if(ndim >= 2) { TotalBinValueSet[1][0] = LowBinValueSet; TotalBinValueSet[1][1] = HighBinValueSet; }
        if(ndim >= 3) { TotalBinValueSet[2][0] = LowBinValueSet; TotalBinValueSet[2][1] = HighBinValueSet; }
    }

    vector<double> xbins_raw[3];
    vector<double> xbins_part[3]; double array_part[3][100];
    vector<int> xbins_ihbs[3];

    for(int idm = 0; idm < ndim; idm++)
    {
        TAxis* haxis = NULL;
        if(idm+1 == 1) haxis = (TAxis*)haim->GetXaxis();
        else if(idm+1 == 2) haxis = (TAxis*)haim->GetYaxis();
        else if(idm+1 == 3) haxis = (TAxis*)haim->GetZaxis();

        unsigned int ibhlow; unsigned int ibhhigh;
        ibhlow = haxis->FindBin(TotalBinValueSet[idm][0]);
        ibhhigh = haxis->FindBin(TotalBinValueSet[idm][1]);
        xbins_raw[idm] = AxisBinArrayGet(haxis);


        for(unsigned int ihb=0; ihb < xbins_raw[idm].size(); ihb++)
        {
////            cout << Form("Dim: %d, For ihb: %d, binvalue: %g",idm,ihb+1,xbins_raw[idm][ihb]) << endl;
            if(ihb+1 >= ibhlow && ihb+1 <= ibhhigh+1)
            {
                xbins_part[idm].push_back(xbins_raw[idm][ihb]);
                xbins_ihbs[idm].push_back(ihb+1);
////                cout << Form("Set In part vector: Dim: %d, For ihb: %d, binvalue: %g",idm,ihb+1,xbins_raw[idm][ihb]) << endl;
////                cout << Form("ihbLow: %d, LowBinCenter: %g, ihbHigh: %d, HighBinCenter: %g",ibhlow,haxis->GetBinCenter(ibhlow),ibhhigh,haxis->GetBinCenter(ibhhigh)) << endl;
            }
        }

        VectortoArray(xbins_part[idm],&(array_part[idm][0]));
    }
        TString HistTitle = haim->GetTitle();
        TString HistName = haim->GetName();
        TString XTitle = haim->GetXaxis()->GetName();
        TString YTitle = haim->GetYaxis()->GetName();
        TString ZTitle = haim->GetZaxis()->GetName();

        if(ndim == 1) hget = new TH1F(Form("%s_part",HistName.Data()),HistTitle,(int)(xbins_part[0].size())-1,&(array_part[0][0]));
        if(ndim == 2) hget = new TH2F(Form("%s_part",HistName.Data()),HistTitle,(int)(xbins_part[0].size())-1,&(array_part[0][0]),(int)(xbins_part[1].size())-1,&(array_part[1][0]));
        if(ndim == 3) hget = new TH3F(Form("%s_part",HistName.Data()),HistTitle,(int)(xbins_part[0].size())-1,&(array_part[0][0]),(int)(xbins_part[1].size())-1,&(array_part[1][0]),(int)(xbins_part[2].size())-1,&(array_part[2][0]));

        if(Option == "noreset")
        {
            double value; double err;
            for(unsigned int ihb_x=0; ihb_x < xbins_ihbs[0].size(); ihb_x++)
            {
                if(ndim < 1) continue;
                if(ndim == 1)
                {
                    value = haim->GetBinContent(xbins_ihbs[0][ihb_x]);
                    err = haim->GetBinError(xbins_ihbs[0][ihb_x]);
                    hget->SetBinContent(ihb_x+1,value);
                    hget->SetBinError(ihb_x+1,err);
                }

                for(unsigned int ihb_y=0; ihb_y < xbins_ihbs[1].size(); ihb_y++)
                {
                    if(ndim < 2) continue;
                    if(ndim == 2)
                    {
                        value = haim->GetBinContent(xbins_ihbs[0][ihb_x],xbins_ihbs[1][ihb_y]);
                        err = haim->GetBinError(xbins_ihbs[0][ihb_x],xbins_ihbs[1][ihb_y]);
                        hget->SetBinContent(ihb_x+1,ihb_y+1,value);
                        hget->SetBinError(ihb_x+1,ihb_y+1,err);
                    }

                    for(unsigned int ihb_z=0; ihb_z < xbins_ihbs[2].size(); ihb_z++)
                    {
                        if(ndim < 3) continue;
                        if(ndim == 3)
                        {
                            value = haim->GetBinContent(xbins_ihbs[0][ihb_x],xbins_ihbs[1][ihb_y],xbins_ihbs[2][ihb_z]);
                            err = haim->GetBinError(xbins_ihbs[0][ihb_x],xbins_ihbs[1][ihb_y],xbins_ihbs[2][ihb_z]);
                            hget->SetBinContent(ihb_x+1,ihb_y+1,ihb_z+1,value);
                            hget->SetBinError(ihb_x+1,ihb_y+1,ihb_z+1,err);
                        }
                    }
                }
            }
        }

//            for(int ihb=0; ihb < xbins_ihbs[0].size(); ihb++)
//            {
//                cout << Form("%d, rawhistbin_id: %d, rawhistbin_value: %g",ihb,xbins_ihbs[0][ihb],xbins_part[0][ihb]) << endl;
//            }
        return hget;
}

TGraphErrors* FluxCompare_Graph(TGraphErrors* gnum, TGraphErrors* gden, TH1F* hline = NULL, TH1F* hlineall = NULL, TCanvas* c = NULL, TString Option = "", double SetYRangeLow = -999, double SetYRangeHigh = -999) // Option: line1 or line2; comp or avege ; '' or nodraw
{
    TGraphErrors* gratio;
    if(c==NULL)
    {
        c = new TCanvas("c","c",1200,900);
    }

    gratio = GraphRatio(gnum,gden,NULL,"noerr");

    TH1F* h_high; TH1F* h_low;

    h_high = (TH1F*)hline->Clone("h_high");
    h_low = (TH1F*)hline->Clone("h_low");

    int N = hline->GetNbinsX();

    for(int ihb=1; ihb <= N; ihb++)
    {
        double value = hline->GetBinContent(ihb);

        h_high->SetBinContent(ihb,1+value);
        h_low->SetBinContent(ihb,1-value);
    }

    TH1F* hall_high = NULL; TH1F* hall_low = NULL;
    if(hlineall != NULL)
    {
        hall_high = (TH1F*)hlineall->Clone("hall_high");
        hall_low = (TH1F*)hlineall->Clone("hall_low");

        int Nall = hlineall->GetNbinsX();

        for(int ihb=1; ihb <= Nall; ihb++)
        {
            double value = hlineall->GetBinContent(ihb);

            hall_high->SetBinContent(ihb,1+value);
            hall_low->SetBinContent(ihb,1-value);
        }
    }

    SetGraphStyle(gratio);

    if(!Option.Contains("avege"))
    {
        c->cd();

        gratio->SetLineColor(4); gratio->SetMarkerColor(4); gratio->SetMarkerStyle(8);
        gratio->SetLineWidth(3); gratio->SetMarkerSize(1.3);
        h_high->SetLineWidth(3); h_high->SetLineColor(8);
        h_low->SetLineWidth(3); h_low->SetLineColor(8);

        if(hlineall != NULL)
        {
            hall_high->SetLineWidth(3); hall_high->SetLineColor(6);
            hall_low->SetLineWidth(3); hall_low->SetLineColor(6);
        }

        gPad->SetLogx(1);
        gratio->Draw("AP"); 
        gratio->GetYaxis()->SetRangeUser(0.8,1.2);

        if(SetYRangeLow != -999) {gratio->GetYaxis()->SetLimits(SetYRangeLow,SetYRangeHigh);}
        h_high->Draw("same");
        h_low->Draw("same");
    }

    if(Option.Contains("avege"))
    {
        TGraphErrors* gnum_ave = (TGraphErrors*)gratio->Clone("gnum_ave");

        TGraphErrors* gden_ave = (TGraphErrors*)gratio->Clone("gden_ave");

        int Nratio = gratio->GetN();

        for(int i=0; i<Nratio; i++)
        {
            double ynum_ave; double yden_ave; double yratio; double xratio;

            gratio->GetPoint(i,xratio,yratio);

            ynum_ave = 2.0/(1+1./yratio);
            yden_ave = 2.0/(1+yratio);

            if(yratio == 0) {ynum_ave = 0; yden_ave = 0;}

            gnum_ave->SetPoint(i,xratio,ynum_ave);
            gden_ave->SetPoint(i,xratio,yden_ave);
        }

        c->cd();
        gnum_ave->SetLineColor(2); gnum_ave->SetMarkerColor(2); gnum_ave->SetMarkerStyle(8);
        gnum_ave->SetLineWidth(3); gnum_ave->SetMarkerSize(1.3);
        gden_ave->SetLineColor(4); gden_ave->SetMarkerColor(4); gden_ave->SetMarkerStyle(8);
        gden_ave->SetLineWidth(3); gden_ave->SetMarkerSize(1.3);
        h_high->SetLineWidth(3); h_high->SetLineColor(8);
        h_low->SetLineWidth(3); h_low->SetLineColor(8);
 
        gPad->SetLogx(1);
        gnum_ave->GetYaxis()->SetRangeUser(0.8,1.2);
        if(SetYRangeLow != -999) {gnum_ave->GetXaxis()->SetLimits(SetYRangeLow,SetYRangeHigh);}
        gnum_ave->Draw("AP"); 
        gden_ave->Draw("P");
        h_high->Draw("same");
        h_low->Draw("same");

        if(hlineall != NULL)
        {
            hall_high->SetLineWidth(3); hall_high->SetLineColor(6);
            hall_low->SetLineWidth(3); hall_low->SetLineColor(6);
            hall_high->Draw("same");
            hall_low->Draw("same");
        }
        TLegend* gl = new TLegend(0.25,0.25,0.5,0.5);
        TString Title_Num = gnum->GetTitle();
        TString Title_Den = gden->GetTitle();
        gl->AddEntry(gnum_ave,Form("%s/Average",Title_Num.Data()),"p");
        gl->AddEntry(gden_ave,Form("%s/Average",Title_Den.Data()),"p");
        gl->SetFillColor(0);
        gl->Draw("same");
    }
    return gratio;
}

TGraphErrors* GraphPointShift(TGraphErrors* graw, TH1F* htmp = NULL, TString Option = "", double SetShiftRatio = 0.05) // Option: '' or changeraw
{
    TGraphErrors* g;

    if(Option.Contains("changeraw"))
    {
        g = graw;
    }
    else{
        TString name = graw->GetName();
        g = (TGraphErrors*)graw->Clone(Form("%s_ShiftX",name.Data()));
    }

    int N = g->GetN();
    double x; double y; double xshift;
    double xerr; double yerr; double xerrshift;
    for(int ip=0; ip < N; ip++)
    {
        g->GetPoint(ip,x,y);
        yerr = g->GetErrorY(ip); xerr = g->GetErrorX(ip);
        if(htmp == NULL)
        {
            xshift = x*(1+SetShiftRatio);
        }
        else {
            int setihb = htmp->FindBin(x);
            double binwidth = htmp->GetBinWidth(setihb);
            xshift = x + binwidth*SetShiftRatio*2;
        }

        xerrshift = xerr/x*xshift;
        g->SetPoint(ip,xshift,y);
        g->SetPointError(ip,xerrshift,yerr);
    }

    return g;
}

double SigvsBKGRatio(TH1F* hsig, TH1F* hbkg, double SetSigEff = 0.9)
{
    double NumAllsig = hsig->Integral();
    double NumAllbkg = hbkg->Integral();
    double Numsig=0; double Numbkg=0; double ratio_sigeff=0; 
    int setih_left=-1; int setih_right = -1;

    int Nhsig = hsig->GetNbinsX();
    int Nhbkg = hbkg->GetNbinsX();

    double SetSigEff_ForCut = (1-SetSigEff)/2.0;
    int ih_SetSigEff_left = -1; double xcenter_SetSigEff_left = 99999;
    int ih_SetSigEff_right = -1; double xcenter_SetSigEff_right = 99999;
    double NumSig_SetSigEff = 0; double NumBkg_SetSigEff = 0;
    double SigvsBKGRatio=0;

    for(int ih=1; ih <= Nhsig; ih++)
    {
        Numsig = hsig->Integral(1,ih);
        ratio_sigeff = Numsig/NumAllsig;

        if(ratio_sigeff >= SetSigEff_ForCut) 
        {
            ih_SetSigEff_left = ih;
            xcenter_SetSigEff_left = hsig->GetBinCenter(ih_SetSigEff_left);
            break;
        }
    }
    for(int ih=Nhsig; ih >= 1; ih--)
    {
        Numsig = hsig->Integral(ih,Nhsig);
        ratio_sigeff = Numsig/NumAllsig;

        if(ratio_sigeff >= SetSigEff_ForCut) 
        {
            ih_SetSigEff_right = ih;
            xcenter_SetSigEff_right = hsig->GetBinCenter(ih_SetSigEff_right);
            break;
        }
    }
    
    NumSig_SetSigEff = hsig->Integral(ih_SetSigEff_left,ih_SetSigEff_right);
    double RatioSigForCheck = NumSig_SetSigEff/NumAllsig;
    if(TMath::Abs(RatioSigForCheck-SetSigEff) > 0.1*SetSigEff)
    {
        NumSig_SetSigEff = -999;
        cout << "Error For Get " << SetSigEff << " percent signal events, the real ratio is: " << RatioSigForCheck << " and the leftxcenter is: " << xcenter_SetSigEff_left << " and the rightxcenter is: " << xcenter_SetSigEff_right << endl;
    }

    setih_left = hbkg->FindBin(xcenter_SetSigEff_left);
    setih_right = hbkg->FindBin(xcenter_SetSigEff_right);
    Numbkg = hbkg->Integral(setih_left,setih_right);

    NumBkg_SetSigEff = Numbkg;

    SigvsBKGRatio = NumSig_SetSigEff/NumBkg_SetSigEff;

    if(NumBkg_SetSigEff == 0) SigvsBKGRatio = 0;

    cout << Form("ih_SetSigEff: [%d,%d], xcenter_SetSigEff: [%g,%g], NumSig_SetSigEff: %g, NumBkg_SetSigEff: %g",ih_SetSigEff_left,ih_SetSigEff_right,xcenter_SetSigEff_left,xcenter_SetSigEff_right,NumSig_SetSigEff,NumBkg_SetSigEff) << endl;
    return SigvsBKGRatio;
}
//Need Update to be usedful at multiple cuts and histograms
//Use Guide: HCutStatis =  LoopCutSaveHist(Cuts_Used,NCuts_Used,rigidity,DB::NBin,DB::Rbins);  At line:271, /afs/cern.ch/user/c/czhang/work/private/Jobfit/Nuclei/FluxCount/FluxCount.C 
//
/*TH2F* LoopCutSaveHist(bool* LCSH_cuts, int LCSH_nlengthcut, double LCSH_xvalue = 0, int LCSH_nx = 1, double* LCSH_binsx = NULL ) // 2D Hist, Y: list of cuts, X: normally rigidity bins
{
    double LCSH_nox[2] = {0,1};
    if(LCSH_binsx == NULL)
    {
        LCSH_binsx = LCSH_nox;
    }

    double LCSH_histxn = -1; double* LCSH_histxbins = NULL; double LCSH_histyn = -1; double LCSH_histylow = -1; double LCSH_histyhigh = -1;

    LCSH_histxn = LCSH_nx;  LCSH_histxbins = LCSH_binsx; 
    LCSH_histyn = LCSH_nlengthcut; LCSH_histylow = 0.5; LCSH_histyhigh = LCSH_nlengthcut+0.5;

    //////static TH2F* LCSH_h = new TH2F("hstatis"," statis of cuts", LCSH_nx, LCSH_binsx, LCSH_nlengthcut,0.5,LCSH_nlengthcut+0.5);
    static TH2F* LCSH_h = new TH2F("hstatis"," statis of cuts", LCSH_histxn, LCSH_histxbins, LCSH_histyn,LCSH_histylow,LCSH_histyhigh);

    for(int iy=0; iy < LCSH_nlengthcut; iy++)
    {
        if(LCSH_cuts[iy]) LCSH_h->Fill(LCSH_xvalue,iy+1);
    }

    return LCSH_h;
}
*/
TH2F* LoopCutSaveHist(TString Name, bool* LCSH_cuts = NULL, int LCSH_nlengthcut = 0, double LCSH_xvalue = 0, int LCSH_nx = 1, double* LCSH_binsx = NULL ) // 2D Hist, Y: list of cuts, X: normally rigidity bins
{
    double LCSH_nox[2] = {0,1};
    if(LCSH_binsx == NULL)
    {
        LCSH_binsx = LCSH_nox;
    }

    double LCSH_histxn = -1; double* LCSH_histxbins = NULL; double LCSH_histyn = -1; double LCSH_histylow = -1; double LCSH_histyhigh = -1;

    LCSH_histxn = LCSH_nx;  LCSH_histxbins = LCSH_binsx; 
    LCSH_histyn = LCSH_nlengthcut; LCSH_histylow = 0.5; LCSH_histyhigh = LCSH_nlengthcut+0.5;

    static map<TString,TH2F*> hmap;
    map<TString,TH2F*>::iterator hit;

    hit = hmap.find(Name);
    if(hit == hmap.end())//No Find, create new
    {
        hmap[Name] = new TH2F(Name," statis of cuts", LCSH_histxn, LCSH_histxbins, LCSH_histyn,LCSH_histylow,LCSH_histyhigh);
    }
    if(LCSH_cuts != NULL) //input exist
    {
        for(int iy=0; iy < LCSH_nlengthcut; iy++)
        {
            if(LCSH_cuts[iy]) hmap[Name]->Fill(LCSH_xvalue,iy+1);
        }
    }

    return hmap[Name];
}
//=============================================N-Dimension Array Save And Read =======================================
////vector<int>* ArraySizeSet(int size_a, int size_b =0, int size_c = 0, int size_d = 0, int size_e = 0, int size_f =0) //6-dimension Max
vector<int> ArraySizeSet(int size_a, int size_b =0, int size_c = 0, int size_d = 0, int size_e = 0, int size_f =0) //6-dimension Max
{
    vector<int> arraysize;
    arraysize.push_back(size_a);
    if(size_b > 0) arraysize.push_back(size_b);
    if(size_c > 0) arraysize.push_back(size_c);
    if(size_d > 0) arraysize.push_back(size_d);
    if(size_e > 0) arraysize.push_back(size_e);
    if(size_f > 0) arraysize.push_back(size_f);

    ////return &(arraysize);
    return arraysize;
}
TString ArrayIDCount(int itotal, vector<int> parraysize)
{
    int N = parraysize.size();
    vector<int> unitvalue;
    vector<int> idvalue;
    TString arrayid = "";
    int itotal_copy = itotal;

    for(int in=0; in < N; in++)
    {
        int tmpvalue = 1;
        for(int isubn=1+in; isubn < N; isubn++)
        {
            tmpvalue *= parraysize[isubn];
        }
        unitvalue.push_back(tmpvalue);
    }
    for(int in=0; in < N; in++)
    {
        int tmpid = 0;
        tmpid = (int)(itotal_copy/unitvalue[in]);
        idvalue.push_back(tmpid);
        arrayid += Form("%d",idvalue[in]);

        itotal_copy -= tmpid*unitvalue[in];
    }

    return arrayid;
}

int ArrrayTotalICount(TString ArrayID, vector<int> parraysize)
{
    int N = parraysize.size();
    vector<int> unitarrayvalue;
    vector<int> unit10value;
    int arraytotali = 0;
    for(int in=0; in < N; in++)
    {
        int tmparrayvalue = 1;
        int tmp10value = 1;
        for(int isubn=1+in; isubn < N; isubn++)
        {
            tmparrayvalue *= parraysize[isubn];
            tmp10value *= 10;
        }
        unitarrayvalue.push_back(tmparrayvalue);
        unit10value.push_back(tmp10value);
    }

    int tmpvalue = ArrayID.Atoi();

    for(int in=0; in < N; in++)
    {
        int tmpid = (int)(tmpvalue/unit10value[in]);
        arraytotali += tmpid*unitarrayvalue[in];

        tmpvalue -= tmpid*unit10value[in];
    }

    return arraytotali;
}

double ArraySaveRead(TString Name, double* array=NULL, vector<int> array_NLength = emptyintvector) //SaveModel: Name: arrayname, double* array and vector<int>* array_NLength need fill ; ReadModel: Name arrayname+id in each dimension, double* array and vector<int>* array_NLength just keep NULL
{
    static map<TString,double> hmap;
    map<TString,double>::iterator hit;

    double outvalue = 0;

    ////if(array != NULL && array_NLength != NULL)
    if(array != NULL)
    {
        int N = array_NLength.size();
        int NTotalArray = 1;
        for(int in=0; in < N; in++)
        {
            NTotalArray *= array_NLength[in];
        }

        for(int iarray=0; iarray < NTotalArray; iarray++)
        {
            TString Name_ArrayID = ArrayIDCount(iarray,array_NLength);
            TString Name_Array = Name+Name_ArrayID;

            ////hit = hmap.find(Name_Array);
            hmap[Name_Array] = array[iarray];
        }
    }
    else
    {
        hit = hmap.find(Name);

        if(hit == hmap.end())//No Find, create new
        {
            outvalue = -9999;
        }
        else
        {
            outvalue = hmap[Name];
        }
    }

    return outvalue;
}
//========================================================================================================
//=======================================About MC Acceptance Count =======================================
double GetReweightIndex_Oxygen(double rigidity_gen)
{
    TString Func_SpectralIndex_Oxygen="[0]*TMath::Log(x)*TMath::Log(x)+[1]*TMath::Log(x)+[2]-[3]*TMath::Power(x,-3)";
    ////double FuncPar_SpectralIndex_Oxygen[4]={0.0496814,-0.5085,-1.426,-185.394};
    ////double FuncXLimit_SpectralIndex_Oxygen[2]={12,1500};
    double FuncPar_SpectralIndex_Oxygen[4]={0.06484,-0.6644,-1.039,-15.04};
    double FuncXLimit_SpectralIndex_Oxygen[2]={2.4,1200};

    static TF1* f = new TF1("SpectralIndexFit",Func_SpectralIndex_Oxygen,FuncXLimit_SpectralIndex_Oxygen[0],FuncXLimit_SpectralIndex_Oxygen[1]);
    static bool Flag_SetParameter = true;
    if(Flag_SetParameter) f->SetParameters(FuncPar_SpectralIndex_Oxygen[0],FuncPar_SpectralIndex_Oxygen[1],FuncPar_SpectralIndex_Oxygen[2],FuncPar_SpectralIndex_Oxygen[3]);

    Flag_SetParameter=false;
    double Index;
    if(rigidity_gen < FuncXLimit_SpectralIndex_Oxygen[0]) {Index = f->Eval(FuncXLimit_SpectralIndex_Oxygen[0]);}
    else if(rigidity_gen > FuncXLimit_SpectralIndex_Oxygen[1]) {Index = f->Eval(FuncXLimit_SpectralIndex_Oxygen[1]);}
    else Index = f->Eval(rigidity_gen); 

    /////return TMath::Power(rigidity_gen,Index);
    return Index;
}
TF1* TF1SaveRead(TString Name, TString FuncString = "", vector<double> FuncPar = emptydoublevector, double funcrangelow = 0, double funcrangehigh = 0)
{
    static map<TString,TF1*> hmap;
    map<TString,TF1*>::iterator hit;

    if(FuncString != "")
    {
        hmap[Name] = new TF1(Name,FuncString,funcrangelow,funcrangehigh);
        int N = FuncPar.size();
        for(int ipar=0; ipar < N; ipar++)
        {
            hmap[Name]->SetParameter(ipar,FuncPar[ipar]);
        }
        return hmap[Name];
    }
    else{
        hit = hmap.find(Name);

        if(hit == hmap.end())//No Find, create new
        {
            return NULL;
        }
        else
        {
            return hmap[Name];
        }
    }

}

double GetReweight(double rigidity_gen, TString ReWeightModel = "InputTF1", TF1* Fit_Set = NULL, TSpline3* Spline_Set = NULL, double RawGenIndex = -1) // May Need MCGenRange[2] For Future to add Normalization function
// ReWeightModel: "Index_-1", "Index_-2.7", "Index_func_oxygen", "InputTF1", "InputTSpline3"
{
    double AimDerivative = NANNumber;
    double RawGenDerivative = NANNumber;

    RawGenDerivative = TMath::Power(rigidity_gen,RawGenIndex); // when generate flux is E^index function

    if(ReWeightModel=="Index_-1") {AimDerivative = TMath::Power(rigidity_gen,-1);}
    else if(ReWeightModel == "Index_-2.7") {AimDerivative = TMath::Power(rigidity_gen,-2.7);}
    else if(ReWeightModel == "Index_func_oxygen") {AimDerivative = TMath::Power(rigidity_gen,GetReweightIndex_Oxygen(rigidity_gen));}
    else if(ReWeightModel == "InputTF1" && Fit_Set != NULL) {AimDerivative = Fit_Set->Eval(rigidity_gen);}
    else if(ReWeightModel == "InputTSpline3" && Spline_Set != NULL) {
        double xinput = rigidity_gen;
        double xmin = Spline_Set->GetXmin(); double xmax = Spline_Set->GetXmax();
        if(xinput < xmin) xinput = xmin;
        if(xinput > xmax) xinput = xmax;
        AimDerivative = Spline_Set->Eval(xinput);
    }
    double reweightvalue = 0;

    reweightvalue = AimDerivative/RawGenDerivative;

    return reweightvalue;
}


TGraph* GetRunMaxEventGraph(TString name, int run = -1, int event = -1) // 2D Hist, Y: list of cuts, X: normally rigidity bins
{
    //    static vector<TGraph> gsave;
    static map<TString,TGraph*> gmap;
    map<TString,TGraph*>::iterator git;

    git = gmap.find(name);
    if(git != gmap.end())
    {
        if(run != -1 && event != -1)
        {
            int N = gmap[name]->GetN();
            int setig = -1;
            for(int ig=0; ig < N; ig++)
            {
                double run_g; double event_g;
                gmap[name]->GetPoint(ig,run_g,event_g);
                if((double)run == run_g)
                {
                    setig = ig;
                    if((double)event > event_g)
                    {
                        gmap[name]->SetPoint(ig,(double)run_g,(double)event);
                    }
                }
            }
            if(setig == -1)
            {
                setig = N;
                gmap[name]->SetPoint(setig,(double)run,(double)event);
            }
        }
    }
    else
    {
        gmap[name] = new TGraph();
        gmap[name]->SetName(name);
        if(run != -1 && event != -1)
        {
            gmap[name]->SetPoint(0,(double)run,(double)event);
        }
    }

    return gmap[name];
}
long int GetTotalEventAllRunGraph(TGraph* gtmp) //TGraph has the point set as the GetRunMaxEventGraph output
{
    double xtmp; double ytmp;
    int run_tmp; int event_tmp;

    int run_dump = -99; int event_dump = 0; // use to correct the dumplicate points: same run number, different event value.
    long int EventTotal = 0; // Output EventsValue
    gtmp->Sort(); // sort points in order, just nearby points could have the same run number.

    int Npoint = gtmp->GetN();

    for(int ip=0; ip < Npoint; ip++)
    {
        gtmp->GetPoint(ip,xtmp,ytmp);
        run_tmp = (int)xtmp; event_tmp = (int)ytmp;

        if(run_tmp == run_dump)
        {
            if(event_tmp > event_dump)
            {
                EventTotal += event_tmp - event_dump;
            }
            else
            {
                continue;
            }
        }
        else if(run_tmp != run_dump)
        {
            EventTotal += event_tmp;
            run_dump = run_tmp; event_dump = event_tmp;
        }
    }
    return EventTotal;
}
void GetTotalEventCountHist(TH1F* histcreate, long int totalevent, double genrange[2], TString ReWeightModel = "InputTF1", TF1* Fit_Set = NULL, double RawGenIndex = -1) // ReWeightModel: "Index_-1", "Index_-2.7", "Index_func_oxygen", "InputTF1"
{
    double scale = 1;
    if(RawGenIndex == -1)
    {
        scale = totalevent/TMath::Log(genrange[1]/genrange[0]);
    }
    else{
        cout << "Error RawGenIndex != -1 still not available" << endl;
    }

    int N = histcreate->GetNbinsX();
    for(int i=1; i <= N; i++)
    {
        double xlow = 0; double xhigh = 0; double y = 0; double yerr = 0; double xcenter = 0;
        xlow = histcreate->GetXaxis()->GetBinLowEdge(i);
        xhigh = histcreate->GetXaxis()->GetBinUpEdge(i);
////        xcenter = histcreate->GetBinCenter(i);

        if(xlow < genrange[0]) xlow = genrange[0];
        if(xhigh > genrange[1]) xhigh = genrange[1];

        xcenter = TMath::Sqrt(xlow*xhigh); //CenterLog

        double index = 0; double value_integralbin = 0;

        if(ReWeightModel=="Index_-1") {index = -1; value_integralbin = TMath::Log(xhigh/xlow);}
        else if(ReWeightModel == "Index_-2.7") {index = -2.7; value_integralbin = (TMath::Power(xhigh,(index+1))-TMath::Power(xlow,(index+1)))/(index+1);}
        else if(ReWeightModel == "Index_func_oxygen") {index = GetReweightIndex_Oxygen(xcenter); value_integralbin = (TMath::Power(xhigh,(index+1))-TMath::Power(xlow,(index+1)))/(index+1);}
        else if(ReWeightModel == "InputTF1") {
            if(Fit_Set != NULL)
            {
                value_integralbin = Fit_Set->Integral(xlow,xhigh); 
            }
            else{
                cout << "Fit_Set Should not be NULL, when ReweightModel is InputTF1" << endl;
            }
        }

        y = scale*value_integralbin;

        //if(y >= 0) yerr = TMath::Sqrt(y);
        //else yerr = 0;
        yerr = 0;

        histcreate->SetBinContent(i,y); histcreate->SetBinError(i,yerr);
    }
}
TH1F* AcceptanceCount(TH1F* hnumerator, long int TotalGenNum, double genrange[2], TString ReWeightModel = "InputTF1", TF1* Fit_Set = NULL, TString HistName = "HAccp", double RawGenIndex = -1, vector<TH1F*>* HAccepRelate = NULL)
// ReWeightModel: "Index_-1", "Index_-2.7", "Index_func_oxygen", "InputTF1"
// HAccepRelate outputs: HAcceptance, Hnumerator, Hdenominator, when input is not NULL, save in. 
{
    const double MCGenAccp = TMath::Pi()*3.9*3.9;
    TH1F* HNum = new TH1F((TH1F)*hnumerator);
    HNum->SetName(HistName+"_Numerator");
    TH1F* HDen = new TH1F((TH1F)*hnumerator);
    HDen->SetName(HistName+"_Deonominator");
    HistClear(HDen);

    GetTotalEventCountHist(HDen,TotalGenNum,genrange,ReWeightModel,Fit_Set,RawGenIndex);

    TH1F* HAccep = (TH1F*)HistOperate(HNum,HDen,Form("%g*[0]/[1]",MCGenAccp),"0");
    HAccep->SetName(HistName);

    cout << "In AcceptanceCount, Before Judge HAccepRelate" << endl;
    if(HAccepRelate != NULL)
    {
        cout << "In AcceptanceCount, The HAccepRelate is not empty" << endl;
        HAccepRelate->clear();
        HAccepRelate->push_back(HAccep);
        HAccepRelate->push_back(HNum);
        HAccepRelate->push_back(HDen);
    }

    return HAccep;
}
//=======================================End MC Acceptance Count =======================================
double CountBinXWithPowerConst(double xlow, double xhigh, double index) // the index of reference flux. Exp: positron xpoint count from E^-3, so index = -3.
{
    double power = -1*(index+1);
    double mainpart = (TMath::Power(xlow,-1*power)-TMath::Power(xhigh,-1*power))/(xhigh-xlow)/power;
    if(mainpart < 0) cout << "Error in Function CountBinXWithPowerConst, mainpart < 0 " << endl;

    double xlw = TMath::Exp(-1./(power+1)*TMath::Log(mainpart));

    return xlw;
}
//====================================== Count Flux ====================================================
TH1F* Hist_CountBinXWithPowerConst(TH1F* hbin, double index, TString Name = "Xflux")
{
    TString Namehbin = hbin->GetName();
    TH1F* haim = (TH1F*)hbin->Clone(Namehbin+Name);
    HistClear(haim);

    int Nhist = hbin->GetNbinsX();
    for(int ihist = 1; ihist <= Nhist; ihist++)
    {
        double xlow; double xhigh; double xfluxcount;
        xlow = hbin->GetXaxis()->GetBinLowEdge(ihist);
        xhigh = hbin->GetXaxis()->GetBinUpEdge(ihist);

        xfluxcount = CountBinXWithPowerConst(xlow,xhigh,index);

        haim->SetBinContent(ihist,xfluxcount);
    }

    return haim;
}
TH1F* TFitTransHist(TF1* fit, TH1F* histxcount, TString HistName) //Attention, histxcount should create from function "Hist_CountBinXWithPowerConst()"
{
    //////TH1F* haim = (TH1F*)histxcount->Clone(Form("%s",HistName.Data()));
    TH1F* haim = new TH1F((TH1F)*((TH1F*)histxcount));
    haim->SetName(HistName);
    HistClear(haim);

    int Nhist = histxcount->GetNbinsX();
    for(int ihist=1; ihist <= Nhist; ihist++)
    {
        double xset; double yset;
        xset = histxcount->GetBinContent(ihist);

        yset = fit->Eval(xset);

        haim->SetBinContent(ihist,yset);
    }

    return haim;
}
TH1F* TMultiFitTransHist(vector<TF1*> fit, TH1F* histxcount, TString HistName)
{
    int Nfunc = fit.size();
    int Nhist = histxcount->GetNbinsX();

    TH1F* haim = new TH1F((TH1F)*((TH1F*)histxcount));
    haim->SetName(HistName);
    HistClear(haim);
    
    for(int ihist=1; ihist <= Nhist; ihist++)
    {
        double xset = 0; double yset = 1;
        xset = histxcount->GetBinContent(ihist);

        yset = 1;
        for(int ifit=0; ifit < Nfunc; ifit++)
        {
            yset *= fit[ifit]->Eval(xset);
        }

        haim->SetBinContent(ihist,yset);
    }

    return haim;
}

TGraphErrors* GetFluxResult(TH1F* hnum, TH1F* hExp, TH1F* hTrigger, TH1F* hAccp, TH1F* hEffCorr, double Xfluxcountindex = 9999, TString GraphName = "") //If use -2.7 count Xflux, input -2.7 to Xfluxcountindex 
{
    bool Flag_XfluxCountIndex = false;
    TH1F* hxflux = NULL;
    if(Xfluxcountindex != 9999)
    {
        Flag_XfluxCountIndex = true;
        hxflux = Hist_CountBinXWithPowerConst(hnum,Xfluxcountindex);
    }

    if(GraphName == "") { GraphName = hnum->GetName(); GraphName = "g"+GraphName; }
    TGraphErrors* gflux = new TGraphErrors();
    gflux->SetName(GraphName);

    int Nhist = hnum->GetNbinsX();
    for(int ihist=1; ihist <= Nhist; ihist++)
    {
        double xset = 0; int i_exp = -1; int i_trigger = -1; int i_accp = -1; int i_effcorr = -1;
        double y_num = 0; double y_err; double y_exp = 0; double y_trigger = 0; double y_accp = 0; double y_effcorr = 0; double y_binwidth = 0; double y_flux = 0; double y_fluxerr = 0;

        if(Flag_XfluxCountIndex)
        {
            xset = hxflux->GetBinContent(ihist); // x count with powerconst flux
        }
        else{
            xset = hnum->GetBinCenter(ihist); // x count with (xlow+xhigh)/2
        }

        y_num = hnum->GetBinContent(ihist);
        y_err = hnum->GetBinError(ihist);
        y_binwidth = hnum->GetBinWidth(ihist);

        i_exp = hExp->FindBin(xset);
        y_exp = hExp->GetBinContent(i_exp);

        i_trigger = hTrigger->FindBin(xset);
        y_trigger = hTrigger->GetBinContent(i_trigger);

        i_accp = hAccp->FindBin(xset);
        y_accp = hAccp->GetBinContent(i_accp);

        i_effcorr = hEffCorr->FindBin(xset);
        y_effcorr = hEffCorr->GetBinContent(i_effcorr);

        y_flux = y_num/(y_exp*y_trigger*y_accp*y_effcorr*y_binwidth);

        y_fluxerr = (y_err/y_num)*y_flux;

        gflux->SetPoint(ihist-1,xset,y_flux);
        gflux->SetPointError(ihist-1,0,y_fluxerr);
    }

    return gflux;
}
//====================================== End Count Flux ================================================
//====================================== About Unfolding ===============================================
int GetGraphPointID(TGraph* g, double xpoint)
{
    int seti = -1;
    int N = g->GetN();
    double xleft = NANNumber; double yleft = NANNumber; int ileft = NANNumber;
    double xright = NANNumber; double yright = NANNumber; int iright = NANNumber;

    for(int ip=0; ip < N; ip++)
    {
        double xtmp; double ytmp;
        g->GetPoint(ip,xtmp,ytmp);
        if(xtmp < xpoint)
        {
            if(xleft == NANNumber || xleft < xtmp)
            {
                xleft = xtmp;
                yleft = ytmp;
                ileft = ip;
            }
        }
        else if(xtmp >= xpoint)
        {
            if(xright == NANNumber || xright > xtmp)
            {
                xright = xtmp;
                yright = ytmp;
                iright = ip;
            }
        }
    }
    if(TMath::Abs(xpoint-xleft) > TMath::Abs(xpoint-xright))
    {
        seti = iright;
    }
    else
    {
        seti = ileft;
    }
    return seti;
}

TSpline3* SplineReCount(TSpline3* fraw, TString XCountFunc = "", TString YCountFunc = "", int Npx = 100000)
{
    TF1* yfunc = NULL; TF1* xfunc = NULL;

    TGraph* g = new TGraph();

    if(YCountFunc != "") yfunc = new TF1("yfunc",YCountFunc,0,1);
    if(XCountFunc != "") xfunc = new TF1("xfunc",XCountFunc,0,1); 

    double xmin = fraw->GetXmin(); double xmax = fraw->GetXmax();
    for(int iloop=0; iloop < Npx; iloop++)
    {
        double xset = NANNumber; double yset = NANNumber; double xrecount = NANNumber; double yrecount = NANNumber;
        xset = xmin + (xmax-xmin)*iloop/(1.0*(Npx-1));
        yset = fraw->Eval(xset);

        if(iloop == 0 || iloop == Npx-1) cout << "xset: " << xset << endl;

        if(YCountFunc != "") { 
            if(YCountFunc.Contains("[0]")) {yfunc->SetParameter(0,xset);}
            yrecount = yfunc->Eval(yset);
        }
        else yrecount = yset;

        if(XCountFunc != "") { 
            if(XCountFunc.Contains("[0]")) {xfunc->SetParameter(0,yset);}
            xrecount = xfunc->Eval(xset);
        }
        else xrecount = xset;

        g->SetPoint(iloop,xrecount,yrecount);
    }

    TString Name = fraw->GetName();
    TSpline3* frecount = new TSpline3(Name+"_ReCount",g,"e2b2");

    return frecount;
}

TGraph* TF1ReCountToGraph(TF1* fraw, TString XCountFunc = "", TString YCountFunc = "", int Npx = 100000)
{
    TF1* yfunc = NULL; TF1* xfunc = NULL;

    TGraph* g = new TGraph();

    if(YCountFunc != "") yfunc = new TF1("yfunc",YCountFunc,0,1);
    if(XCountFunc != "") xfunc = new TF1("xfunc",XCountFunc,0,1); 

    double xmin = fraw->GetXmin(); double xmax = fraw->GetXmax();
    for(int iloop=0; iloop < Npx; iloop++)
    {
        double xset = NANNumber; double yset = NANNumber; double xrecount = NANNumber; double yrecount = NANNumber;
        xset = xmin + (xmax-xmin)*iloop/(1.0*(Npx-1));
        yset = fraw->Eval(xset);

        if(iloop == 0 || iloop == Npx-1) cout << "xset: " << xset << endl;

        if(YCountFunc != "") { 
            if(YCountFunc.Contains("[0]")) {yfunc->SetParameter(0,xset);}
            yrecount = yfunc->Eval(yset);
        }
        else yrecount = yset;

        if(XCountFunc != "") { 
            if(XCountFunc.Contains("[0]")) {xfunc->SetParameter(0,yset);}
            xrecount = xfunc->Eval(xset);
        }
        else xrecount = xset;

        g->SetPoint(iloop,xrecount,yrecount);
    }

    TString Name = fraw->GetName();
    g->SetName("g_"+Name+"_ReCount");

    return g;
}

/*TSpline3* GetFluxFuncForReweight(TGraph* gflux, TString Option = "|AllPoints|") // |AllPoints|, |SetPointsNuclei|, |SetPointsNuclei|Logx| 
{
    double Esetpoint[]={2.2,4.2,10.2,20.2,102.2,300,600,1000,1500,2500};
    const int Nsetpoint = sizeof(Esetpoint)/sizeof(Esetpoint[0]);
    double xsetpoint[Nsetpoint]; double ysetpoint[Nsetpoint];
    TSpline3* funcflux = NULL;

    if(Option.Contains("|AllPoints|"))
    {
        funcflux = new TSpline3("funcflux",gflux,"e2b2");
    }
    else if(Option.Contains("|SetPointsNuclei|"))
    {
        for(int isetp=0; isetp < Nsetpoint; isetp++)
        {
            double xinput;
            if(Option.Contains("|Logx|"))
            {
                xinput = TMath::Log(Esetpoint[isetp]);
            }
            else{
                xinput = Esetpoint[isetp];
            }
            int tmpi = GetGraphPointID(gflux,xinput);
            gflux->GetPoint(tmpi,xsetpoint[isetp],ysetpoint[isetp]);
        }

        funcflux = new TSpline3("funcflux",xsetpoint,ysetpoint,Nsetpoint,"e2b2");
    }

    return funcflux;
}
*/
double funcfromspline(double* x, double* t)
{
    return globalSpline->Eval(x[0]);
}
TF1* TransSplineToTF1(TSpline3* spline)
{
    TString Name = spline->GetName();
    double xmin = spline->GetXmin(); double xmax = spline->GetXmax();
    globalSpline = spline;
    TF1* fline = new TF1(Name,funcfromspline,xmin,xmax,0);

    ////globalSpline = NULL;

    return fline;
}
TSpline3* GetFluxFuncForReweight(TGraph* gflux, TString Option = "|AllPoints|", TString Name = "funcflux", bool Debug = false) // |AllPoints|, |SetPointsNuclei|, |SetPointsNuclei|Spline-LogXLogY| 
{
    double Esetpoint[]={2.2,4.2,10.2,20.2,102.2,300,600,1000,1500,2500};
    const int Nsetpoint = sizeof(Esetpoint)/sizeof(Esetpoint[0]);
    double xsetpoint[Nsetpoint]; double ysetpoint[Nsetpoint];
    int idsetpoint[Nsetpoint];
    TString StringCountFunc[2][2] = { //[GraphTrans,SplineTransBack][x,y]
        "TMath::Log(x)","TMath::Log(x)",
        "TMath::Exp(x)","TMath::Exp(x)"
    };

    TSpline3* funcflux = NULL;
    TSpline3* funcfluxOut = NULL;
    TGraph* gaim = NULL;

    for(int ip=0; ip < Nsetpoint; ip++)
    {
        idsetpoint[ip] = GetGraphPointID(gflux,Esetpoint[ip]);
    }

    if(Option.Contains("|Spline-LogXLogY|"))
    {
        gaim = GraphReCount(gflux,StringCountFunc[0][0],StringCountFunc[0][1]);
    }
    else 
    {
        gaim = gflux;
    }

    if(Option.Contains("|AllPoints|"))
    {
        funcflux = new TSpline3(Name,gaim,"e2b2");
    }
    else if(Option.Contains("|SetPointsNuclei|"))
    {
        for(int isetp=0; isetp < Nsetpoint; isetp++)
        {
            gaim->GetPoint(idsetpoint[isetp],xsetpoint[isetp],ysetpoint[isetp]);
        }

        funcflux = new TSpline3(Name,xsetpoint,ysetpoint,Nsetpoint,"e2b2");
    }

    if(Option.Contains("|Spline-LogXLogY|"))
    {
        funcfluxOut = SplineReCount(funcflux,StringCountFunc[1][0],StringCountFunc[1][1]);
    }
    else
    {
        funcfluxOut = funcflux;
    }

    if(Debug)
    {
        new TCanvas();
        SetGraphStyle(gaim);
        gaim->Draw("AP");
        funcflux->Draw("Lsame");
    }

    return funcfluxOut;
}
double FluxIterationCheck(TGraphErrors* gafter, TGraphErrors* gpre, double rangelow = NANNumber, double rangehigh = NANNumber)
{
    int N = gafter->GetN(); int Npre = gpre->GetN();
    if(N!=Npre)
    {
        cout << "The Npoints of gafter and gpre are not the same, can't deal" << endl;
        return -1;
    }

    double ydiffertotal = 0;
    int ncount = 0;
    for(int ip=0; ip < N; ip++)
    {
        double xpre; double ypre; double xafter; double yafter;

        gafter->GetPoint(ip,xafter,yafter);
        gpre->GetPoint(ip,xpre,ypre);

        if(xpre != xafter)
        {
            cout << "gafter and gpre not in the same x points, can't deal" << endl;
            return -1;
        }

        if(rangelow != NANNumber)
        {
            if(xpre < rangelow) {
                cout << "Set RangeLow: " << rangelow << " and Graph PointX: " << xpre << ", Not count in, pass" << endl;
                continue;
            }
        }

        if(rangehigh != NANNumber)
        {
            if(xpre > rangehigh) {
                cout << "Set RangeHigh: " << rangehigh << " and Graph PointX: " << xpre << ", Not count in, pass" << endl;
                continue;
            }
        }

        ydiffertotal += TMath::Abs(yafter-ypre)/ypre;
        ncount++;
    }

    ydiffertotal /= ncount*1.0;

    return ydiffertotal;
}

TGraph* GetIndexFromFlux(TGraphErrors* GFlux, double FluxDrawScale = 0, TString Option = "")
{
    TGraph* gindex = new TGraph();

    int N = GFlux->GetN();
    TString Name = GFlux->GetName();
    gindex->SetName(Form("%s_IndexFlux",Name.Data()));

    for(int ip=0; ip < N-1; ip++)
    {
        double y0; double y1; double x0; double x1;
        double y_index; double x_index;

        GFlux->GetPoint(ip,x0,y0);
        GFlux->GetPoint(ip+1,x1,y1);

        x_index = (x1+x0)/2;
        y_index = (TMath::Log(y1)-TMath::Log(y0))/(TMath::Log(x1)-TMath::Log(x0));
        y_index = y_index-FluxDrawScale;

        gindex->SetPoint(ip,x_index,y_index);
    }

    return gindex;
}

//====================================== End Unfolding =================================================
/*TGraph* GetRunMaxEventGraph_useVector(TString name, int run = -1, int event = -1) // 2D Hist, Y: list of cuts, X: normally rigidity bins
  {
  static vector<TGraph*> gsave;
  static vector<TString> stringsave;

  int setil = -1;
  for(int il=0; il < (int)stringsave.size(); il++)
  {
  if(stringsave[il] == name) setil = il;
  }

  if(setil == -1)
  {
  setil = stringsave.size();
  gsave.push_back(new TGraph());
  stringsave.push_back(name);

  gsave[setil]->SetName(name);

  if(run != -1 && event != -1) 
  {
  gsave[setil]->SetPoint(0,run,event);
  }
  }
  else if(run != -1 && event != -1)
  {
  int N=gsave[setil]->GetN();
  int setig = -1;
  for(int ig=0; ig < N; ig++)
  {
  double run_g; double event_g;
  gsave[setil]->GetPoint(ig,run_g,event_g);
  if(run_g == run)
  {
  setig = ig;
  if(event > event_g)
  {
  gsave[setil]->SetPoint(ig,run_g,event);
  }
  }
  }
  if(setig == -1)
  {
  setig = N;
  gsave[setil]->SetPoint(setig,run,event);
  }
  }

  return gsave[setil];
  }
  */
vector<TString> ReadInput(int argc, char* argv[], TString KeyWord, TString* ReturnString = NULL) // List: -OutDir:,-OutName:,-InFile:,-Option:,All
{//When RootFunc.h has been writen, public TString[] could be {"outdir","outname","infile","option"} and could be reset through function()

    vector<TString> allreturn;
    TString List[]={"-OutDir:","-OutName:","-InFile:","-Option:"};
    int NList = sizeof(List)/sizeof(List[0]);

    if(KeyWord != "All")
    {
        for(int c=1; c < argc; c++)
        {
            TString tmpstring = argv[c];

            if(tmpstring.Contains(KeyWord))
            {
                (*ReturnString) = tmpstring.ReplaceAll(KeyWord,"");
                allreturn.push_back((*ReturnString));
            }
        }
    }
    else
    {
        for(int il=0; il < NList; il++)
        {
            bool Flag_Get = false;
            TString tmpkeyword = List[il];
            for(int c=1; c < argc; c++)
            {
                TString tmpstring = argv[c];

                if(tmpstring.Contains(tmpkeyword))
                {
                    TString tmpreturnstring = tmpstring.ReplaceAll(tmpkeyword,"");
                    allreturn.push_back(tmpreturnstring);
                    Flag_Get = true;
                    break;
                }
            }
            if(!Flag_Get)
            {
                TString nonfind = "NULL";
                allreturn.push_back(nonfind);
            }
        }
    }

    return allreturn;
}

double GetCorrectedGeomRigidity(double rigidity, int igeom = 1, int opt = 0, int version = 0)
{
     double geomcorr = 0;
     if(igeom == 0||igeom == 1)geomcorr += 1.0*0.0001;
     double nrig = rigidity;
     if(opt==0)nrig = 1./(1./rigidity + geomcorr);
     else nrig = 1./(1./rigidity - geomcorr);
 
     return nrig;
}

/*double SigvsBKGRatio(TH1F* hsig, TH1F* hbkg, double SetSigEff = 0.9)
  {
  double NumAllsig = hsig->Integral();
  double NumAllbkg = hbkg->Integral();
  double Numsig=0; double Numbkg=0; double ratio_sigeff=0; int setih=-1;

  int Nhsig = hsig->GetNbinsX();
  int Nhbkg = hbkg->GetNbinsX();

  int ih_SetSigEff = -1; double xcenter_SetSigEff = 99999;
  double NumSig_SetSigEff = 0; double NumBkg_SetSigEff = 0;
  double SigvsBKGRatio=0;

  for(int ih=1; ih <= Nhsig; ih++)
  {
  Numsig = hsig->Integral(1,ih);
  ratio_sigeff = Numsig/NumAllsig;

  if(ratio_sigeff >= SetSigEff) 
  {
  ih_SetSigEff = ih;
  xcenter_SetSigEff = hsig->GetBinCenter(ih_SetSigEff);
  xcenter_SetSigEff = hsig->GetBinCenter(ih_SetSigEff);
  NumSig_SetSigEff = Numsig;
  break;
  }
  }

  setih = hbkg->FindBin(xcenter_SetSigEff);
  Numbkg = hbkg->Integral(1,setih);

  NumBkg_SetSigEff = Numbkg;

  SigvsBKGRatio = NumSig_SetSigEff/NumBkg_SetSigEff;

  if(NumBkg_SetSigEff == 0) SigvsBKGRatio = 0;

  cout << Form("ih_SetSigEff: %d, xcenter_SetSigEff: %g, NumSig_SetSigEff: %g, NumBkg_SetSigEff: %g",ih_SetSigEff,xcenter_SetSigEff,NumSig_SetSigEff,NumBkg_SetSigEff) << endl;
  return SigvsBKGRatio;
  }
  */


/*TH1F* HistPartGet(TH1F* haim, double LowBinValueSet, double HighBinValueSet, TString Option = "noreset") // "reset" or "noreset"
  {
  TH1F* hget;
  int ibhlow; int ibhhigh;
  ibhlow = haim->FindBin(LowBinValueSet);
  ibhhigh = haim->FindBin(HighBinValueSet);

  vector<double> xbins_raw = HistBinArrayGet(haim);
  vector<double> xbins_part; double array_part[100];
  vector<int> xbins_ihbs;

  for(int ihb=0; ihb < xbins_raw.size(); ihb++)
  {
  cout << Form("For ihb: %d, binvalue: %g",ihb+1,xbins_raw[ihb]) << endl;
  if(ihb+1 >= ibhlow && ihb+1 <= ibhhigh)
  {
  xbins_part.push_back(xbins_raw[ihb]);
  xbins_ihbs.push_back(ihb+1);
  cout << Form("Set In part vector: For ihb: %d, binvalue: %g",ihb+1,xbins_raw[ihb]) << endl;
  }
  }

  VectortoArray(xbins_part,array_part);

  TString HistTitle = haim->GetTitle();
  TString HistName = haim->GetName();
  TString XTitle = haim->GetXaxis()->GetName();
  TString YTitle = haim->GetYaxis()->GetName();

  hget = new TH1F(Form("%s_part",HistName.Data()),HistTitle,(int)(xbins_part.size())-1,&(array_part[0]));

  if(Option == "noreset")
  {
  double value; double err;
  for(int ihb=0; ihb < xbins_ihbs.size(); ihb++)
  {
  value = haim->GetBinContent(xbins_ihbs[ihb]);
  err = haim->GetBinError(xbins_ihbs[ihb]);

  hget->SetBinContent(ihb+1,value);
  hget->SetBinError(ihb+1,err);
  }
  }

  return hget;
  }
  */
/*double HistNDChi2Count(TH1* Hist_data,TH1* Hist_fit,TString ProjectDim = "xy",double SetDim2Value = -1,double* ProjectChi2Value = NULL) // The Return is the combine chi2 value // ProjectChi2Value is used to get the seperately x and y project chi2 value //
  {
  double totalvalue = 0;
  int NDim_data = Hist_data->GetDimension();
  int NDim_fit = Hist_fit->GetDimension();

  if(NDim_data != NDim_fit) 
  {
  cout << "!!!=== Error in Func 'HistNDChi2Count': NDim_data != NDim_fit ===!!!" << endl;
  return -999;
  }

  if(NDim_data == 1)
  {
  totalvalue = Hist1DChi2Count(Hist_data,Hist_fit);
  return totalvalue;
  }
  else if(NDim_data == 2)
  {
  TH1* htmp_data; TH1* htmp_fit;
  TString ProjectName[2] = {"x","y"};
  for(int id=0; id < 2; id++)
  {
  if(ProjectDim.Contains(ProjectName[id]))
  {
  htmp_data = HistGetFromND(Hist_data,"",ProjectName[id],SetDim2Value)

  }
  }

  }

  }
  */
