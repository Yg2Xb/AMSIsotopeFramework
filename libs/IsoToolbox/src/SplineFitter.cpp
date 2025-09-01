/***********************************************************
 * File: SplineFitter.cpp
 *
 * Implementation file for the SplineFitter class.
 *
 * History:
 * 20250831 - Ported and refactored by Gemini
 ***********************************************************/

#include <IsoToolbox/SplineFitter.h>
#include <string>

namespace IsoToolbox {
namespace Tools {

void SplineFitter::Init() {
    qname = 0;
    qnnodes = 0;
    qxnodes = qynodes = 0;
    qh = 0;
    qgr = 0;
    qfun = 0;
    qspline = 0;
    _verbose = 0;
}

void SplineFitter::Clear() {
    if (qfun) {
        delete qfun;
        qfun = 0;
    }
    if (qspline) {
        delete qspline;
        qspline = 0;
    }
}

SplineFitter::SplineFitter() {
    Init();
}

SplineFitter::SplineFitter(const char *name, TH1 *h, int Nnodes, double *nodes, int sfitmode, const char *sopt, double svalbeg, double svalend, double serror) {
    Init();
    SetName(name);
    SetHisto(h);
    SetSplinePar(sfitmode, sopt, svalbeg, svalend, serror);
    SetNode(Nnodes, nodes);
}

SplineFitter::SplineFitter(const char *name, TGraphErrors *gr, int Nnodes, double *nodes, int sfitmode, const char *sopt, double svalbeg, double svalend, double serror) {
    Init();
    SetName(name);
    SetGraph(gr);
    SetSplinePar(sfitmode, sopt, svalbeg, svalend, serror);
    SetNode(Nnodes, nodes);
}

SplineFitter::SplineFitter(const char *name, int Nnodes, double *nodes, double *ynodes, int sfitmode, const char *sopt, double svalbeg, double svalend, double serror) {
    Init();
    SetName(name);
    SetSplinePar(sfitmode, sopt, svalbeg, svalend, serror);
    SetNode(Nnodes, nodes, ynodes);
}

SplineFitter::~SplineFitter() {
    Clear();
    if(qxnodes) delete [] qxnodes;
    if(qynodes) delete [] qynodes;
}


void SplineFitter::SetName(const char *name) {
    qname = name;
    int i = 0;
    while (gDirectory->FindObject(Form("%s_fun", qname)) != 0) {
        qname = Form("%s_f%d", name, i);
        i++;
    }
    Clear();
}

void SplineFitter::SetHisto(TH1 *h) {
    qh = h;
    qgr = 0;
    Clear();
}

void SplineFitter::SetGraph(TGraphErrors *gr) {
    qgr = gr;
    qh = 0;
    Clear();
}

void SplineFitter::SetNode(int Nnodes, double *nodes, double *ynodes) {
    if (qxnodes) delete[] qxnodes;
    if (qynodes) delete[] qynodes;
    qnnodes = Nnodes;
    qxnodes = new double[Nnodes];
    qynodes = new double[Nnodes];
    for (int i = 0; i < qnnodes; i++) {
        qxnodes[i] = nodes[i];
        if (ynodes != 0) qynodes[i] = ynodes[i];
    }
    if (ynodes != 0) {
        qh = 0;
        qgr = 0;
    }
    Clear();
}

int SplineFitter::GetNodeX(double *xnodes) {
    for (int i = 0; i < qnnodes; i++) {
        xnodes[i] = qxnodes[i];
    }
    return qnnodes;
}

void SplineFitter::SetSplinePar(int sfitmode, const char *sopt, double svalbeg, double svalend, double serror) {
    _sopt = sopt;
    _svalbeg = svalbeg;
    _svalend = svalend;
    _serror = serror;
    _sfitmode = sfitmode;
    Clear();
}

void SplineFitter::SetVerbose(int verbose) {
    _verbose = verbose;
}

double SplineFitter::operator()(double *x, double *par) {
    double xv = x[0];
    if (_sfitmode & (LogX | LogX2)) {
        xv = log(fabs(x[0]));
    }

    int update = 0;
    if (qspline == 0) update = 1;
    else {
        for (int ino = 0; ino < qnnodes; ino++) {
            double parn = par[ino];
            if (_sfitmode & LogY) {
                parn = log(fabs(parn));
            }
            if (fabs(qynodes[ino] - parn) > _serror * fabs(qynodes[ino])) {
                update = 1;
                break;
            }
        }
    }
    if (update) {
        if (qspline) delete qspline;
        double *qxnodesn = new double[qnnodes];
        for (int ino = 0; ino < qnnodes; ino++) {
            qxnodesn[ino] = qxnodes[ino];
            qynodes[ino] = par[ino];
            if (_sfitmode & (LogX | LogX2)) qxnodesn[ino] = log(fabs(qxnodes[ino]));
            if (_sfitmode & LogY) qynodes[ino] = log(fabs(par[ino]));
        }
        qspline = new TSpline3(Form("%s_spline", qname), qxnodesn, qynodes, qnnodes, _sopt, _svalbeg, _svalend);
        delete[] qxnodesn;
    }

    if (_sfitmode & ExtrapolateFlux) {
        int isx0 = 0;
        double x0 = xv;
        if (xv > log(fabs(qxnodes[qnnodes - 1]))) {
            x0 = log(fabs(qxnodes[qnnodes - 1]));
            isx0 = 2;
        } //End
        else if (xv < log(fabs(qxnodes[0]))) {
            x0 = log(fabs(qxnodes[0]));
            isx0 = 1;
        } //Begin
        double kx = qspline->Derivative(x0);
        if (_verbose >= 1) {
            std::cout << "xv=" << x[0] << " log(xv)=" << xv << " index=" << kx << std::endl;
        }
        double y0 = qspline->Eval(x0);
        double b = y0 - kx * x0;
        if (isx0) return exp(kx * xv + b);
    }

    for (int ilh = 0; ilh < 2; ilh++) {
        bool islinear = false;
        if (ilh == 0) islinear = (_sfitmode & ExtrapolateLB);
        else islinear = (_sfitmode & ExtrapolateLE);
        if (!islinear) continue;
        double x0 = 0, y0 = 0;
        if (ilh == 0) qspline->GetKnot(0, x0, y0);
        else qspline->GetKnot(qspline->GetNp() - 1, x0, y0);
        double kx = qspline->Derivative(x0);
        double b = y0 - kx * x0;
        double yv = 0;
        if (ilh == 0 && xv < x0) {
            yv = kx * xv + b;
        } else if (ilh == 1 && xv > x0) {
            yv = kx * xv + b;
        } else continue;
        if (_sfitmode & LogY) yv = exp(yv);
        return yv;
    }

    double yv = qspline->Eval(xv);
    if (_sfitmode & LogY) yv = exp(qspline->Eval(xv));
    return yv;
}

TF1* SplineFitter::GetFun(double xmin, double xmax) {
    if (qfun) {
        if (xmin != 0 || xmax != 0) qfun->SetRange(xmin, xmax);
        return qfun;
    }
    qfun = new TF1(Form("%s_fun", qname), this, xmin, xmax, qnnodes);
    for (int ino = 0; ino < qnnodes; ino++) {
        if (qgr) qfun->SetParameter(ino, qgr->Eval(qxnodes[ino]));
        else if (qh) {
            double y0 = qh->Interpolate(qxnodes[ino]);
            if (_sfitmode & LogY) {
                if (y0 < 1e-12) y0 = 1e-12;
            }
            qfun->SetParameter(ino, y0);
        } else qfun->SetParameter(ino, qynodes[ino]);
    }
    qfun->SetNpx(100000);
    return qfun;
}

TGraphErrors* SplineFitter::HistoToGraph(int cb) {
    TGraphErrors *gr = new TGraphErrors();
    for (int ib = 1; ib <= qh->GetNbinsX(); ib++) {
        double xv = qh->GetXaxis()->GetBinCenterLog(ib);
        if (cb == 1) xv = qh->GetXaxis()->GetBinCenter(ib);
        double yv = qh->GetBinContent(ib);
        double ey = qh->GetBinError(ib);
        int ip = gr->GetN();
        gr->SetPoint(ip, xv, yv);
        gr->SetPointError(ip, 0, ey);
    }
    return gr;
}

TFitResultPtr SplineFitter::DoFit(Option_t *option, Option_t *goption, Double_t xxmin, Double_t xxmax) {
    if (qgr) {
        return qgr->Fit(GetFun(xxmin, xxmax), option, goption, xxmin, xxmax);
    } else if (_sfitmode & LogX2) {
        qgr = HistoToGraph(0);
        TFitResultPtr fitresult = qgr->Fit(GetFun(xxmin, xxmax), option, goption, xxmin, xxmax);
        std::string opt(option);
        if (opt.find("+") == std::string::npos) qh->GetListOfFunctions()->Clear();
        
        qh->GetListOfFunctions()->Add(GetFun(xxmin, xxmax));
        (TVirtualFitter::GetFitter())->SetObjectFit(qh);
        delete qgr;
        qgr = 0;
        return fitresult;
    } else {
        return qh->Fit(GetFun(xxmin, xxmax), option, goption, xxmin, xxmax);
    }
}

int SplineFitter::FillFlux(TH1 *hflux, double scaleindex, double scalenorm) {
    for (int ib = 1; ib <= hflux->GetNbinsX(); ib++) {
        double xv = hflux->GetXaxis()->GetBinCenterLog(ib);
        double yv = GetFun()->Eval(xv) * pow(xv, scaleindex) * scalenorm;
        hflux->SetBinContent(ib, yv);
        hflux->SetBinError(ib, 0);
    }
    return 0;
}

TF1* splineFit(TH1* hist, double* xpoints, int npoint, int qOption, TString splineOption, TString name, double fitRangeLow, double fitRangeHigh) {
    SplineFitter* spfit = new SplineFitter("spfit", hist, npoint, xpoints, qOption, splineOption.Data());
    TF1* spfit_func;

    if (fitRangeLow == -9999 && fitRangeHigh == -9999) {
        spfit->DoFit();
        spfit_func = spfit->GetFun();
    } else {
        spfit->DoFit("","", fitRangeLow, fitRangeHigh);
        spfit_func = spfit->GetFun(fitRangeLow, fitRangeHigh);
    }

    spfit_func->SetName(name);
    // The SplineFitter object is managed by ROOT's list of functions,
    // but it's safer to delete it if we created it.
    // However, the function object `spfit_func` uses `spfit` internally, so we cannot delete it here.
    // This will lead to a memory leak. Proper management is needed.
    // For now, let's assume ROOT handles it or accept the leak as in legacy code.
    return spfit_func;
}

TF1* splineFit(TGraph* graph, double* xpoints, int npoint, int qOption, TString splineOption, TString name, double fitRangeLow, double fitRangeHigh) {
    SplineFitter* spfit = new SplineFitter("spfit", (TGraphErrors*)graph, npoint, xpoints, qOption, splineOption.Data());
    TF1* spfit_func;

    if (fitRangeLow == -9999 && fitRangeHigh == -9999) {
        spfit->DoFit();
        spfit_func = spfit->GetFun();
    } else {
        spfit->DoFit("","", fitRangeLow, fitRangeHigh);
        spfit_func = spfit->GetFun(fitRangeLow, fitRangeHigh);
    }

    spfit_func->SetName(name);
    // Same memory leak issue as above.
    return spfit_func;
}

} // namespace Tools
} // namespace IsoToolbox
