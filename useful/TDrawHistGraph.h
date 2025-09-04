#include <iostream>
#include "THStack.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TText.h"
#include "TFile.h"
class TDrawHistGraph{
    private:
    const static int BasicDrawNum = 6;
    const static int DrawNumMax = 9;
    TString GlobalPlotDir;
    TLegend * gl;
    THStack * hs;
    TMultiGraph *gs;
    TFile* fsave;
//    int ColorList[];
//    int StyleList[];

    public:
    TCanvas * call;
    TCanvas * cc;
    int CanvasDivideNum = -1;
    bool Flag_CanvasAutoPrint = true;

    TH1F * htemp[BasicDrawNum];
    TGraphErrors * gtemp[BasicDrawNum];

    TH1F * htemp_add[BasicDrawNum];
    TGraphErrors * gtemp_add[BasicDrawNum];

    TString legendtitle[BasicDrawNum];
    TString legendtitle_add[BasicDrawNum];
    TString labeltitle[3];
    double *legendcoord; // [xstart,xend, ystart, yend]
    bool FlagHistSet, FlagGraphSet, FlagSetLogy, FlagSetLogx;
    TString PlotDir, GifShowName;
    double HistScaleValue; //double MarkerStyle;
    double SetYMin, SetYMax, SetXMin, SetXMax;
    int LineWidth, MarkerStyle;
    int OptStat, OptFit;

    int ColorList[DrawNumMax];
    int StyleList[DrawNumMax];
    int ConstColorList[DrawNumMax];
    int ConstStyleList[DrawNumMax];

    TString ShadowName;
    double ShadowRange[BasicDrawNum][2];
    TH1F* hshadow[BasicDrawNum];
    bool ShadowFlag[BasicDrawNum];
    double ShadowRatio[BasicDrawNum];
    double ShadowRatioCood[BasicDrawNum][2];

    bool FileSaveFlag; 
    bool HistErrFlag;


    TDrawHistGraph(TString AllPlotDir = "");
    ~TDrawHistGraph();
    void SetHist(TH1 * hh);
    void SetHStack(THStack * hh);
    void SetGraph(TGraphErrors * hh);
    void SetMultiGraph(TMultiGraph * hh);

    TH1F* SetShadowRangeHist(TH1F* htemp, TString Title = "copy", int color = 1, double Rlow = -9999, double Rup = 9999);
    void InitalHist();
    void SaveDataHistInital();
    void InitalGraph();

    int ClassHistDraw();
    int BasicHistDraw(TH1F * Htempdraw[BasicDrawNum], TString LegendTitle[BasicDrawNum], TString LabelTitle[3], double* LegendCoord, bool FlagHistSet = true, TString PlotDir = "", double HistScaleValue = -1, bool FlagSetLogy = false, bool FlagSetLogx = false, double SetYMin = -1, double SetYMax = -1, double SetXMin = -1, double SetXMax = -1, TString GifShowName = "", int LineWidth = -1, int OptStat = 0, int OptFit = 0);

    int ClassGraphDraw();
    int BasicGraphDraw(TGraphErrors * Gtempdraw[BasicDrawNum], TString LegendTitle[BasicDrawNum], TString LabelTitle[3], double* LegendCoord, bool FlagGraphSet = true, TString PlotDir = "", int MarkerStyle = -1, bool FlagSetLogy = false, bool FlagSetLogx = false, double SetYMin = -1, double SetYMax = -1, double SetXMin = -1, double SetXMax = -1, TString GifShowName = "", int LineWidth = -1, int OptStat = 0, int OptFit = 0);
};

TH1F* TDrawHistGraph::SetShadowRangeHist(TH1F* htemp, TString Title, int color, double Rlow, double Rup)
{
    TH1F* hcopy = (TH1F*)(htemp->Clone(Title));

    int nlength = htemp->GetNbinsX();
    for(int il=0; il < nlength; il++)
    {
        double xbinwidth = htemp->GetBinWidth(il+1);
        double xbinlow = htemp->GetBinLowEdge(il+1);
        double xbinup = xbinlow+xbinwidth;

        if(xbinup <= Rlow && Rlow != -9999) hcopy->SetBinContent(il+1,0);
        if(xbinlow >= Rup && Rup != 9999) hcopy->SetBinContent(il+1,0);
    }

    hcopy->SetFillColor(color);
    hcopy->SetFillStyle(3004);

    return &(hcopy[0]);
}

void TDrawHistGraph::SaveDataHistInital()
{
    for(int il=0; il < BasicDrawNum; il++)
    {
        ShadowRatio[il] = 0;
    }
}
void TDrawHistGraph::InitalHist()
{
    for(int il=0; il < BasicDrawNum; il++)
    {
        htemp[il] = NULL; //Input Histogram
        htemp_add[il] = NULL; //Input Histogram
        legendtitle[il] = ""; //Legend of each histogram
        legendtitle_add[il] = ""; //Legend of each histogram

        ColorList[il] = ConstColorList[il];
        StyleList[il] = ConstStyleList[il];

        ShadowFlag[il] = false;
        hshadow[il] = NULL;
        ShadowRange[il][0] = -9999; ShadowRange[il][1] = 9999;
        ShadowRatioCood[il][0] = 0.16; ShadowRatioCood[il][1] = 0.86;
    }
    labeltitle[0] = ""; labeltitle[1] = ""; labeltitle[2] = ""; // X title, Y title, PlotTitle
    //legendcoord[0][0] = 0;legendcoord[0][1] = 0;legendcoord[1][0] = 0;legendcoord[1][1] = 0; // Legend coordinate [x,y][start,end]
    legendcoord = NULL; // Legend coordinate [xstart,xend,ystart,yend]
    FlagHistSet=true; // true: Make the basic settting of histogram; false: if the histogram has been used and basic setted once. 
    PlotDir=GlobalPlotDir; GifShowName="";// PlotDir: the ${PlotDir}.ps file you want to write in; GifShowName: the ${GifShowName}.gif file you want to create
    HistScaleValue=-1; // the value of the integral of histogram after scale. -1: No Scale
    FlagSetLogy=false; FlagSetLogx=false; // flase: Normal axis; true: Log axis
    SetYMin=-1; SetYMax=-1; SetXMin=-1; SetXMax=-1; // The Setting Range of Pave
    LineWidth=-1; // The LineWidth you want to set 
    OptStat=0; OptFit=0; // Set Pave Stat and Fit Info
    ShadowName = "shadow";

    FileSaveFlag = false;
    HistErrFlag = false;

    cc = NULL;
    CanvasDivideNum = -1;
    Flag_CanvasAutoPrint = true;
}

void TDrawHistGraph::InitalGraph()
{
    for(int il=0; il < BasicDrawNum; il++)
    {
        gtemp[il] = NULL; //Input Histogram
        gtemp_add[il] = NULL; //Input Histogram
        legendtitle[il] = ""; //Legend of each histogram
        legendtitle_add[il] = ""; //Legend of each histogram

        ColorList[il] = ConstColorList[il];
        StyleList[il] = ConstStyleList[il];
    }
    labeltitle[0] = ""; labeltitle[1] = ""; labeltitle[2] = ""; // X title, Y title, PlotTitle
////    legendcoord[0][0] = 0;legendcoord[0][1] = 0;legendcoord[1][0] = 0;legendcoord[1][1] = 0; // Legend coordinate [x,y][start,end]
    legendcoord = NULL; // Legend coordinate [xstart,xend,ystart,yend]
    FlagGraphSet=true; // true: Make the basic settting of histogram; false: if the histogram has been used and basic setted once.
    PlotDir=GlobalPlotDir; GifShowName=""; // PlotDir: the ${PlotDir}.ps file you want to write in; GifShowName: the ${GifShowName}.gif file you want to create
    MarkerStyle=-1; LineWidth=-1; // The LineWidth and MarkerStyle you want to set (with this set, all Graph members has the same Marker Style. -1: different Graph has different Marker Style 
    FlagSetLogy=false; FlagSetLogx=false; // flase: Normal axis; true: Log axis
    SetYMin=-1; SetYMax=-1; SetXMin=-1; SetXMax=-1; // The Setting Range of Pave
    OptStat=0; OptFit=0; // Set Pave Stat and Fit Info

    cc = NULL;
    CanvasDivideNum = -1;
    Flag_CanvasAutoPrint = true;
}

void TDrawHistGraph::SetHist(TH1 * hh)
{
    hh->GetXaxis()->SetLabelSize(.06);
    hh->GetYaxis()->SetLabelSize(.06);
    hh->GetXaxis()->SetTitleSize(.07);
    hh->GetYaxis()->SetTitleSize(.07);
    hh->GetXaxis()->CenterTitle();
    hh->GetYaxis()->CenterTitle();
    hh->GetXaxis()->SetNdivisions(505);
    hh->GetYaxis()->SetNdivisions(505);
    hh->GetXaxis()->SetTitleOffset(.9);
    hh->GetYaxis()->SetTitleOffset(.9);
    hh->SetLineWidth(3);
}
void TDrawHistGraph::SetHStack(THStack * hh)
{
    hh->GetXaxis()->SetLabelSize(.06);
    hh->GetYaxis()->SetLabelSize(.06);
    hh->GetXaxis()->SetTitleSize(.07);
    hh->GetYaxis()->SetTitleSize(.07);
    hh->GetXaxis()->CenterTitle();
    hh->GetYaxis()->CenterTitle();
    hh->GetXaxis()->SetNdivisions(505);
    hh->GetYaxis()->SetNdivisions(505);
    hh->GetXaxis()->SetTitleOffset(.9);
    hh->GetYaxis()->SetTitleOffset(.9);
    //    hh->SetLineWidth(3);
}
void TDrawHistGraph::SetGraph(TGraphErrors * hh)
{
    hh->GetXaxis()->SetLabelSize(.06);
    hh->GetYaxis()->SetLabelSize(.06);
    hh->GetXaxis()->SetTitleSize(.07);
    hh->GetYaxis()->SetTitleSize(.07);
    hh->GetXaxis()->CenterTitle();
    hh->GetYaxis()->CenterTitle();
    hh->GetXaxis()->SetNdivisions(505);
    hh->GetYaxis()->SetNdivisions(505);
    hh->GetXaxis()->SetTitleOffset(.9);
    hh->GetYaxis()->SetTitleOffset(.9);
    hh->SetLineWidth(3);
}
void TDrawHistGraph::SetMultiGraph(TMultiGraph * hh)
{
    hh->GetXaxis()->SetLabelSize(.06);
    hh->GetYaxis()->SetLabelSize(.06);
    hh->GetXaxis()->SetTitleSize(.07);
    hh->GetYaxis()->SetTitleSize(.07);
    hh->GetXaxis()->CenterTitle();
    hh->GetYaxis()->CenterTitle();
    hh->GetXaxis()->SetNdivisions(505);
    hh->GetYaxis()->SetNdivisions(505);
    hh->GetXaxis()->SetTitleOffset(.9);
    hh->GetYaxis()->SetTitleOffset(.9);
    //    hh->SetLineWidth(3);
}
//int TDrawHistGraph::BasicHistDraw(TH1F * Htempdraw[BasicDrawNum], TString LegendTitle[BasicDrawNum], TString LabelTitle[3], double LegendCoord[2][2], bool FlagHistSet = true, TString PlotDir = "", double HistScaleValue = -1, bool FlagSetLogy = false, bool FlagSetLogx = false, double SetYMin = -1, double SetYMax = -1, double SetXMin = -1, double SetXMax = -1, TString GifShowName = "", int LineWidth = -1, int OptStat = 0, int OptFit = 0)
int TDrawHistGraph::BasicHistDraw(TH1F * Htempdraw[BasicDrawNum], TString LegendTitle[BasicDrawNum], TString LabelTitle[3], double* LegendCoord, bool FlagHistSet, TString PlotDir, double HistScaleValue, bool FlagSetLogy, bool FlagSetLogx, double SetYMin, double SetYMax, double SetXMin, double SetXMax, TString GifShowName, int LineWidth, int OptStat, int OptFit)
{
    SaveDataHistInital();
    ////    int ColorList[] = {4,2,1,3,6,9,5,7,8};
    ////    int ColorList[] = {4,2,1,8,6,7,5,9,3};
    ////    const int NColorNum = sizeof(ColorList)/sizeof(ColorList[0]);
    ////    if(NColorNum < BasicDrawNum) {std::cout << Form("Error for BasicHistDraw Number: %d, which is less than ColorNumber: %d",BasicDrawNum, NColorNum); return -1;}

 //   cout << "~~~~ HistScaleValue: " << HistScaleValue << " ~~~~" << endl;

    double HistEntries[BasicDrawNum];
    bool StatFlag = false;
    if(OptStat != 0 || Htempdraw[1] == NULL) StatFlag = true;

    memset(HistEntries,0,sizeof(HistEntries));

    if(cc == NULL)
    {
        cc = new TCanvas("cc",Form("Canvas for Histogram Drawing"),1200,800);
    }
    else
    {
        Flag_CanvasAutoPrint = false;
        if(CanvasDivideNum == -1) {cc->cd();}
        else {cc->cd(CanvasDivideNum);}
    }
    gPad->SetLeftMargin(.15);
    gPad->SetBottomMargin(.15);
    gPad->SetLogy(FlagSetLogy);
    gPad->SetLogx(FlagSetLogx);
    gStyle->SetOptStat(OptStat);
    ////    if( !Htempdraw[1] ) {gStyle->SetOptStat(1111111); std::cout << "=== Settting OptStat === " << std::endl;}
    if( StatFlag ) {gStyle->SetOptStat(1111111); std::cout << "=== Settting OptStat " << OptStat << std::endl;}
    else gStyle->SetOptStat(OptStat);

    gStyle->SetOptFit(OptFit);

////    if(LegendCoord[0][0] == 0 && LegendCoord[0][1] == 0) gl = new TLegend(0.75,0.75,0.99,0.9);
////    else gl = new TLegend(LegendCoord[0][0],LegendCoord[1][0],LegendCoord[0][1],LegendCoord[1][1]);
    if(LegendCoord == NULL) gl = new TLegend(0.75,0.75,0.99,0.9);
    else gl = new TLegend(LegendCoord[2*0+0],LegendCoord[2*1+0],LegendCoord[2*0+1],LegendCoord[2*1+1]);
    hs = new THStack("hs","Comparison of Histogram");

    for(int il=0; il < BasicDrawNum; il++)
    {
        if(!Htempdraw[il]) continue;

        std::cout << "Pass Hist " << il << std::endl;
        if(FlagHistSet)
        {
            SetHist(Htempdraw[il]);
            Htempdraw[il]->GetYaxis()->SetTitle(LabelTitle[1]);
            Htempdraw[il]->GetXaxis()->SetTitle(LabelTitle[0]);
            Htempdraw[il]->SetTitle(LabelTitle[2]);
            Htempdraw[il]->SetLineWidth(3);
        }


        if(HistScaleValue != -1 && HistScaleValue != 0) {
            HistEntries[il] = Htempdraw[il]->Integral();
            Htempdraw[il]->Scale(HistScaleValue/HistEntries[il]);

//            cout << Form("~~~ HistScaleValue: %g, Listnum: %d, totalEvents: %g, ScaleRatio: %g",HistScaleValue,il,HistEntries[il],HistScaleValue/HistEntries[il]) << endl;
        }

        if(LineWidth != -1)
        {
            Htempdraw[il]->SetLineWidth(LineWidth);
        }

        Htempdraw[il]->SetLineColor(ColorList[il]);
        if( StatFlag ) {Htempdraw[il]->SetStats(1111111); std::cout << "=== Settting HistStat === " << std::endl;}
        if( !StatFlag ) {Htempdraw[il]->SetStats(0); std::cout << "===No Settting HistStat === " << std::endl;}

        if(HistErrFlag == true)
        {
            gStyle->SetErrorX(0);
            Htempdraw[il]->SetMarkerStyle(8);
            Htempdraw[il]->SetMarkerSize(1.3);
            Htempdraw[il]->SetMarkerColor(ColorList[il]);
        }

        if(SetYMin != -1 || SetYMax != -1) {Htempdraw[il]->GetYaxis()->SetRangeUser(SetYMin,SetYMax);}
        if(SetXMin != -1 || SetXMax != -1) {Htempdraw[il]->GetXaxis()->SetRangeUser(SetXMin,SetXMax);}
        //if(SetYMin != -1 || SetYMax != -1) {Htempdraw[il]->SetMinimum(SetYMin); Htempdraw[il]->SetMaximum(SetYMax);}
        //if(SetXMin != -1 || SetXMax != -1) {Htempdraw[il]->GetXaxis()->SetRangeUser(SetXMin,SetXMax);}

        if(LegendTitle[il] != "") gl->AddEntry(Htempdraw[il],LegendTitle[il],"l");

        hs->Add(Htempdraw[il]);
    }

    if(HistErrFlag == true) hs->Draw("enostack");
    else hs->Draw("nostack");
    
    SetHStack(hs);
    hs->GetXaxis()->SetTitle(LabelTitle[0]);
    hs->GetYaxis()->SetTitle(LabelTitle[1]);
    hs->SetTitle(LabelTitle[2]);

    //if(SetYMin != -1 || SetYMax != -1) {hs->GetYaxis()->SetRangeUser(SetYMin,SetYMax);}
    if(SetYMin != -1 || SetYMax != -1) {hs->SetMinimum(SetYMin); hs->SetMaximum(SetYMax);}
    if(SetXMin != -1 || SetXMax != -1) {hs->GetXaxis()->SetRangeUser(SetXMin,SetXMax);}

    std::cout << "SetYmin: " << SetYMin << std::endl;
    std::cout << "SetYmax: " << SetYMax << std::endl;

    for(int il=0; il < BasicDrawNum; il++)
    {
        if(!Htempdraw[il]) continue;
        if(ShadowFlag[il] == true) 
        {
            hshadow[il] = SetShadowRangeHist(Htempdraw[il],Form("%s_%d",ShadowName.Data(),il),ColorList[il],ShadowRange[il][0],ShadowRange[il][1]);
            ShadowRatio[il] = hshadow[il]->Integral()/Htempdraw[il]->Integral();
            hshadow[il]->Draw("same");

            TText* tmptext  = new TText();
//            tmptext->DrawTextNDC(0.16,0.86,Form("ResEleRatio: %.2g%%",ShadowRatio[il]*100));
            tmptext->DrawTextNDC(ShadowRatioCood[il][0],ShadowRatioCood[il][1],Form("ResEleRatio: %.2g%%",ShadowRatio[il]*100));
////            tmptext->DrawTextNDC(0.16,0.86,Form("ResEleRatio: %g",ShadowRatio[il]));
            std::cout << "-----ID: " << il << " has proceding!! ----" << std::endl;
        }
    }

    for(int il=0; il < BasicDrawNum; il++)
    {
        if(!gtemp_add[il]) continue;
        SetGraph(gtemp_add[il]);
        gtemp_add[il]->SetMarkerStyle(8);
        if(legendtitle_add[il] && legendtitle_add[il] != "") gl->AddEntry(gtemp_add[il],legendtitle_add[il],"p");
        gtemp_add[il]->Draw("P");
    }

    gl->SetFillColor(0);
    if(Htempdraw[1]) gl->Draw();

    if(Flag_CanvasAutoPrint)
    {
        cc->Print(Form("%s.ps",PlotDir.Data()));
        if(GifShowName != "") cc->Print(Form("%s.gif",GifShowName.Data()));
        if(GifShowName != "") cc->Print(Form("%s.ps",GifShowName.Data()));
        if(GifShowName != "" && FileSaveFlag == false) cc->Print(Form("%s.root",GifShowName.Data()));
    }

    if(FileSaveFlag)
    {
        fsave = new TFile(Form("%s.root",GifShowName.Data()),"recreate");
        for(int il=0; il < BasicDrawNum; il++)
        {
            if(!Htempdraw[il]) continue;
            Htempdraw[il]->Write();
            if(ShadowFlag[il]) hshadow[il]->Write();
        }
        hs->Write();
        cc->Write();

        fsave->Close();
    }
    //    if(GifShowName != "") cc->Print(Form("%s_ROOTSave.root",GifShowName.Data()));

    if(HistScaleValue != -1 && HistScaleValue != 0)
    {
        for(int il=0; il < BasicDrawNum; il++)
        {
            if(!Htempdraw[il]) continue;
            Htempdraw[il]->Scale(HistEntries[il]/HistScaleValue);
        }
    }

    if(Flag_CanvasAutoPrint)
    {
        delete hs;
    }
    InitalHist();

    return 0;
}
//int TDrawHistGraph::BasicGraphDraw(TGraphErrors * Gtempdraw[BasicDrawNum], TString LegendTitle[BasicDrawNum], TString LabelTitle[3], double* LegendCoord, bool FlagGraphSet = true, TString PlotDir = "", int MarkerStyle = -1, bool FlagSetLogy = false, bool FlagSetLogx = false, double SetYMin = -1, double SetYMax = -1, double SetXMin = -1, double SetXMax = -1, TString GifShowName = "", int LineWidth = -1, int OptStat = 0, int OptFit = 0)
int TDrawHistGraph::BasicGraphDraw(TGraphErrors * Gtempdraw[BasicDrawNum], TString LegendTitle[BasicDrawNum], TString LabelTitle[3], double* LegendCoord, bool FlagGraphSet, TString PlotDir, int MarkerStyle, bool FlagSetLogy, bool FlagSetLogx, double SetYMin, double SetYMax, double SetXMin, double SetXMax, TString GifShowName, int LineWidth, int OptStat, int OptFit)
{
    ////    int ColorList[] = {4,2,1,8,6,7,5,9,3};
    ////    int StyleList[] = {20,21,22,23,24,25,26,32,27};
    const int NColorNum = sizeof(ColorList)/sizeof(ColorList[0]);

    if(NColorNum < BasicDrawNum) {std::cout << Form("Error for BasicDraw Number: %d, which is less than ColorNumber: %d",BasicDrawNum, NColorNum); return -1;}

    ////cc = new TCanvas("cc",Form("Canvas for Histogram Drawing"),1200,800);
    if(cc == NULL)
    {
        cc = new TCanvas("cc",Form("Canvas for Histogram Drawing"),1200,800);
    }
    else
    {
        Flag_CanvasAutoPrint = false;
        if(CanvasDivideNum == -1) {cc->cd();}
        else {cc->cd(CanvasDivideNum);}
    }

    gPad->SetLeftMargin(.15);
    gPad->SetBottomMargin(.15);
    gPad->SetLogy(FlagSetLogy);
    gPad->SetLogx(FlagSetLogx);
    gStyle->SetOptStat(OptStat);
    gStyle->SetOptFit(OptFit);

    ////    if(LegendCoord[0][0] == 0 && LegendCoord[0][1] == 0) gl = new TLegend(0.75,0.75,0.99,0.9);
    ////    else gl = new TLegend(LegendCoord[0][0],LegendCoord[1][0],LegendCoord[0][1],LegendCoord[1][1]);
    if(LegendCoord == NULL) gl = new TLegend(0.75,0.75,0.99,0.9);
    else gl = new TLegend(LegendCoord[2*0+0],LegendCoord[2*1+0],LegendCoord[2*0+1],LegendCoord[2*1+1]);

    gs = new TMultiGraph("gs",LabelTitle[2]);

    for(int il=0; il < BasicDrawNum; il++)
    {
        if(!Gtempdraw[il]) continue;

        if(FlagGraphSet)
        {
            SetGraph(Gtempdraw[il]);
            Gtempdraw[il]->GetYaxis()->SetTitle(LabelTitle[1]);
            Gtempdraw[il]->GetXaxis()->SetTitle(LabelTitle[0]);
            Gtempdraw[il]->SetTitle(LabelTitle[2]);
            Gtempdraw[il]->SetLineWidth(3);
            Gtempdraw[il]->SetMarkerSize(1);
        }

        Gtempdraw[il]->SetLineColor(ColorList[il]);
        Gtempdraw[il]->SetMarkerColor(ColorList[il]);
        if(MarkerStyle == -1) Gtempdraw[il]->SetMarkerStyle(StyleList[il]);
        else Gtempdraw[il]->SetMarkerStyle(MarkerStyle);

        if(LineWidth != -1) Gtempdraw[il]->SetLineWidth(LineWidth);

        gl->AddEntry(Gtempdraw[il],LegendTitle[il],"p");
        gs->Add(Gtempdraw[il]);
    }

    gs->Draw("AP");
    gs->GetXaxis()->SetTitle(LabelTitle[0]);
    gs->GetYaxis()->SetTitle(LabelTitle[1]);
    gs->SetTitle(LabelTitle[2]);

    if(SetYMin != -1 || SetYMax != -1) {gs->GetYaxis()->SetRangeUser(SetYMin,SetYMax);}
    if(SetXMin != -1 || SetXMax != -1) {gs->GetXaxis()->SetRangeUser(SetXMin,SetXMax);}
    SetMultiGraph(gs);

    for(int il=0; il < BasicDrawNum; il++)
    {
        if(!htemp_add[il]) continue;
        SetHist(htemp_add[il]);
        htemp_add[il]->SetLineWidth(3);
        if(legendtitle_add[il] && legendtitle_add[il] != "") gl->AddEntry(htemp_add[il],legendtitle_add[il],"l");
        htemp_add[il]->Draw("same");
    }

    gl->SetFillColor(0);
    if(Gtempdraw[1]) gl->Draw();

    if(Flag_CanvasAutoPrint)
    {
        cc->Print(Form("%s.ps",PlotDir.Data()));
        if(GifShowName != "") cc->Print(Form("%s.gif",GifShowName.Data()));
        if(GifShowName != "") cc->Print(Form("%s.ps",GifShowName.Data()));
    }
    //    if(GifShowName != "") cc->Print(Form("%s_ROOTSave.root",GifShowName.Data()));

    InitalGraph();

    return 0;
}
//TDrawHistGraph::TDrawHistGraph(int &StatusFlag)
////TDrawHistGraph::TDrawHistGraph(TString AllPlotDir = "")
TDrawHistGraph::TDrawHistGraph(TString AllPlotDir)
{

    int StatusFlag = 0;
    int InitColorList[] = {4,2,1,8,6,7,5,9,3};
    int InitStyleList[] = {20,21,22,23,24,25,26,32,27};

    const int NColorNum = sizeof(InitColorList)/sizeof(InitColorList[0]);
    if(NColorNum < BasicDrawNum) {std::cout << Form("Error for BasicDraw Number: %d, which is larger than ColorNumber: %d",BasicDrawNum, NColorNum); StatusFlag = -1;}
    if(NColorNum > DrawNumMax) {std::cout << Form("Error for BasicDrawMax Number: %d, which is less than ColorNumber: %d",DrawNumMax, NColorNum); StatusFlag = -2;}
    for(int i=0; i < NColorNum; i++)
    {
        ConstColorList[i] = InitColorList[i];
        ConstStyleList[i] = InitStyleList[i];
    }

    GlobalPlotDir = "";

    if(AllPlotDir != "")
    {
        GlobalPlotDir=AllPlotDir;
        call = new TCanvas("call","Create ps file",1200,800);
        call->Print(Form("%s.ps[",GlobalPlotDir.Data()));
    }

    InitalHist();
    InitalGraph();
}
TDrawHistGraph::~TDrawHistGraph()
{
    if(GlobalPlotDir != "")
    {
        call->Print(Form("%s.ps]",GlobalPlotDir.Data()));
    }
}
int TDrawHistGraph::ClassHistDraw()
{
    int tmpstat =  BasicHistDraw(htemp, legendtitle, labeltitle, legendcoord, FlagHistSet, PlotDir, HistScaleValue, FlagSetLogy, FlagSetLogx, SetYMin, SetYMax, SetXMin, SetXMax, GifShowName, LineWidth, OptStat, OptFit);

    return tmpstat;
}
int TDrawHistGraph::ClassGraphDraw()
{
    //    int tmpstat =  TDrawHistGraph::BasicHistDraw(htemp, legendtitle, labeltitle, legendcoord, FlagHistSet, PlotDir, HistScaleValue, FlagSetLogy, FlagSetLogx, SetYMin, SetYMax, SetXMin, SetXMax, GifShowName, LineWidth, OptStat, OptFit);
    int tmpstat = BasicGraphDraw(gtemp, legendtitle, labeltitle, legendcoord, FlagGraphSet, PlotDir, MarkerStyle, FlagSetLogy, FlagSetLogx, SetYMin, SetYMax, SetXMin, SetXMax, GifShowName, LineWidth, OptStat, OptFit);

    return tmpstat;
}
