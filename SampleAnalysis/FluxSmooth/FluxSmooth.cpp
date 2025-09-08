#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TROOT.h"

#include "../Tool.h"  // 提供 SplineFit()

using namespace std;

// --------- 分段点生成函数 ---------
static vector<double> generateXPoints(const TGraphAsymmErrors* g, int nTail=3) {
    vector<double> pts;
    int n = g->GetN();
    if (n <= 0) return pts;
    double* xs = g->GetX();

    int i = 0;
    while (i < n - nTail) {
        pts.push_back(xs[i]);
        cout << "[DEBUG] section point i=" << i << "  x=" << xs[i] << endl;

        if (i < 10) i += 5;        // 低能密集
        else if (i < 50) i += 6;  // 中段适中
        else if (i < 200) i += 4;  // 中段适中
        else i += 3;              // 高能稀疏
    }

    // 保证最后一段用最后一个点收尾
    pts.push_back(xs[n-3]);
    cout << "[DEBUG] section point i=" << n-2 << "  x=" << xs[n-2] << endl;
    pts.push_back(xs[n-1]);
    //cout << "[DEBUG] final section point i=" << n-1 << "  x=" << xs[n-1] << endl;

    return pts;
}

// --------- spline拟合 wrapper ---------
static TF1* performSplineFit(TGraphAsymmErrors* g, const string& baseName) {
    auto xpts = generateXPoints(g);
    cout << "[DEBUG] Performing spline fit for " << baseName 
         << " with " << xpts.size() << " segments" << endl;

    try {
        auto fit = SplineFit(g, xpts.data(), xpts.size(), 0x38, "b2e1",
                             baseName.c_str(), xpts.front(), xpts.back());
        fit->SetLineWidth(2);
        fit->SetLineColor(kBlue+1);

        cout << "[DEBUG] Fit ready for " << baseName 
             << "  Chisq=" << fit->GetChisquare()
             << "  NDF=" << fit->GetNDF()
             << "  Chi2/NDF=" << (fit->GetNDF()>0 ? fit->GetChisquare()/fit->GetNDF() : -1)
             << endl;

        return fit;
    } catch (...) {
        cerr << "[ERROR] SplineFit failed for " << baseName << endl;
        return nullptr;
    }
}

// --------- Info box ---------
static void drawInfoBox(TF1* f, double x1=0.55,double y1=0.75,double x2=0.89,double y2=0.89) {
    auto box = new TPaveText(x1,y1,x2,y2,"NDC");
    box->SetFillStyle(0); box->SetBorderSize(0); box->SetTextSize(0.035);
    int ndf = f->GetNDF();
    double chi2 = f->GetChisquare();
    box->AddText(Form("chi2/ndf = %.2f/%d = %.2f", chi2, ndf, (ndf>0?chi2/ndf:0)));
    box->Draw();
}

// ========== 主执行函数 ==========
void FluxSmooth() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    string indir = "/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg";
    vector<pair<string,string>> files = {
        {"Be", indir + "/AMS2011to2018PhysReport_BerylliumFlux.root"},
        {"B",  indir + "/AMS2011to2018PhysReport_BoronFlux.root"},
        {"C",  indir + "/AMS2011to2018PhysReport_CarbonFlux.root"},
        {"N",  indir + "/AMS2011to2018PhysReport_NitrogenFlux.root"},
        {"O",  indir + "/AMS2011to2018PhysReport_OxygenFlux.root"}
    };

    TF1 *fBe=nullptr,*fB=nullptr,*fC=nullptr,*fN=nullptr,*fO=nullptr;
    TF1 *fBe7=nullptr,*fBe9=nullptr,*fBe10=nullptr,*fB10=nullptr,*fB11=nullptr;

    string pdfOut = "/eos/user/z/zixuan/Isotope/FluxSmooth/FluxSmooth.pdf";
    auto c0 = new TCanvas("c0","Flux",800,600);
    c0->Print((pdfOut+"[").c_str());

    for (auto& [name,fpath]:files) {
        cout << "\n[INFO] Processing " << name << " from " << fpath << endl;
        TFile fin(fpath.c_str(),"READ");
        auto g = dynamic_cast<TGraphAsymmErrors*>(fin.Get("graph1"));
        if (!g) { cerr << "[ERROR] Missing graph1 in " << fpath << endl; continue; }

        auto f = performSplineFit(g, name+"_spline");
        if (!f) continue;
        f->SetNpx(3000);
        f->SetRange(1.0,3300.0);

        auto cv = new TCanvas(("c_"+name).c_str(),name.c_str(),800,600);
	cv->SetLogx();
	cv->SetLogy();
        g->SetMarkerStyle(20);
        g->Draw("AP");
	g->GetYaxis()->SetRangeUser(1e-10,10);
        f->Draw("same");
        drawInfoBox(f);

        if (name=="Be") {
            fBe=f;
            fBe7=new TF1("Be7",[f](double* x,double*){return 0.7*f->Eval(x[0]);},1,3300,0);
            fBe9=new TF1("Be9",[f](double* x,double*){return 0.2*f->Eval(x[0]);},1,3300,0);
            fBe10=new TF1("Be10",[f](double* x,double*){return 0.1*f->Eval(x[0]);},1,3300,0);
            fBe7->SetLineColor(kRed); fBe9->SetLineColor(kGreen+2); fBe10->SetLineColor(kMagenta);
            fBe7->Draw("same"); fBe9->Draw("same"); fBe10->Draw("same");
        }
        else if (name=="B") {
            fB=f;
            fB10=new TF1("B10",[f](double* x,double*){return 0.3*f->Eval(x[0]);},1,3300,0);
            fB11=new TF1("B11",[f](double* x,double*){return 0.7*f->Eval(x[0]);},1,3300,0);
            fB10->SetLineColor(kRed); fB11->SetLineColor(kGreen+2);
            fB10->Draw("same"); fB11->Draw("same");
        }
        else if (name=="C") { fC=f; }
        else if (name=="N") { fN=f; }
        else if (name=="O") { fO=f; }

        cv->Print(pdfOut.c_str());
        delete cv;
    }

    c0->Print((pdfOut+"]").c_str());
    delete c0;

    TFile fout("/eos/user/z/zixuan/Isotope/FluxSmooth/FluxSmooth.root","RECREATE");
    if (fBe) fBe->Write("Be");
    if (fBe7) fBe7->Write("Be7");
    if (fBe9) fBe9->Write("Be9");
    if (fBe10) fBe10->Write("Be10");
    if (fB) fB->Write("B");
    if (fB10) fB10->Write("B10");
    if (fB11) fB11->Write("B11");
    if (fC) fC->Write("C");
    if (fN) fN->Write("N");
    if (fO) fO->Write("O");
    fout.Close();

    cout << "\n[INFO] Done. ROOT saved to FluxSmooth.root, PDF to FluxSmooth.pdf" << endl;
}
