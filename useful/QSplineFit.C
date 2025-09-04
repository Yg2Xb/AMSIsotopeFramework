///--Author Q.Yan(Spline Tool)
#ifndef _QSpline_
#define _QSpline_
#include "TH1.h"
#include "TF1.h"
#include "TSpline.h"
#include "TGraphErrors.h"
#include "TFitResultPtr.h"
#include "TDirectory.h"
#include "TVirtualFitter.h"
#include <iostream>
///--Spline Fit Tool
class QSplineFit{
 public:
   enum FitMode{
      LinearXY=0x1, //y=f(X)
      LogX=0x2,// y=f(Log(X))
      LogY=0x4,// Log(Y)=f(X)
      LogX2=0x8,// y=f(Log(X))//sqrt(bl*bh)
      ExtrapolateLB=0x10,//Linear Extrapolate Begin
      ExtrapolateLE=0x20,//Linear Extrapolate End
      ExtrapolateFlux=0x40,//Y=CX^a
    }; 
 public:
  QSplineFit();
  QSplineFit(const char *name,TH1 *h,int Nnodes,double *nodes,int sfitmode=QSplineFit::LinearXY,const char *sopt="b1e1",double svalbeg=0, double svalend=0,double serror=1e-12);
  QSplineFit(const char *name,TGraphErrors *gr,int Nnodes,double *nodes,int sfitmode=QSplineFit::LinearXY,const char *sopt="b1e1",double svalbeg=0, double svalend=0,double serror=1e-12);
  QSplineFit(const char *name,int Nnodes,double *nodes,double *ynodes,int sfitmode=QSplineFit::LinearXY,const char *sopt="b1e1",double svalbeg=0, double svalend=0,double serror=1e-12);//Fix Nodes
  virtual ~QSplineFit(){};
  void SetName(const char *name);
  void SetHisto(TH1 *h);//SetHisto
  void SetGraph(TGraphErrors *gr);
  void SetNode(int Nnodes,double *nodes,double *ynodes=0);//SetNode
  int  GetNodeX(double *xnodes);
  void SetSplinePar(int sfitmode=QSplineFit::LinearXY,const char *sopt="b1e1",double svalbeg=0, double svalend=0,double serror=1e-12);//(a+bx)'=0 //Set Spline Par
  void SetVerbose(int verbose=0);
  TF1* GetFun(double xmin=0,double xmax=0);
  virtual TFitResultPtr DoFit(Option_t* option = "", Option_t* goption = "", Double_t xmin = 0, Double_t xmax = 0);
  int  FillFlux(TH1 *hflux,double scaleindex=0,double scalenorm=1);
  double operator()(double *x, double *par);
 private:
  int Init();///Init  for Begin
  int Clear();//Clear qspline+qfun
  TGraphErrors* HistoToGraph(int cb=0);///centerbin or logbin
  const char *qname;
  int     qnnodes;
  double *qxnodes;
  double *qynodes;
//---
  const char *_sopt;
  double  _serror;
  double _svalbeg;
  double _svalend;
  int    _sfitmode;
  int    _verbose;
//---
  TH1* qh;
  TF1 *qfun;
  TSpline3 *qspline;
  TGraphErrors *qgr;
};

///---
int  QSplineFit::Init(){
  qname=0;
  qnnodes=0;
  qxnodes=qynodes=0;
  qh=0;
  qgr=0;
  qfun=0;
  qspline=0;
 _verbose=0;
 return 0;
}

///---
int  QSplineFit::Clear(){
  if(qfun)delete qfun;
  if(qspline)delete qspline;
  qfun=0;
  qspline=0;
  return 0;
}

///---
QSplineFit::QSplineFit(){
  Init();
}

///---
QSplineFit::QSplineFit(const char *name,TH1 *h,int Nnodes,double *nodes,int sfitmode,const char *sopt,double svalbeg, double svalend,double serror){//Constant
  Init();
  SetName(name);
  SetHisto(h);
  SetSplinePar(sfitmode,sopt,svalbeg,svalend,serror);//SplinePar
  SetNode(Nnodes,nodes);
}

///---
QSplineFit::QSplineFit(const char *name,TGraphErrors *gr,int Nnodes,double *nodes,int sfitmode,const char *sopt,double svalbeg, double svalend,double serror){//Constant
  Init();
  SetName(name);
  SetGraph(gr);
  SetSplinePar(sfitmode,sopt,svalbeg,svalend,serror);//SplinePar
  SetNode(Nnodes,nodes);
}

///---
QSplineFit::QSplineFit(const char *name,int Nnodes,double *nodes,double *ynodes,int sfitmode,const char *sopt,double svalbeg, double svalend,double serror){//Constant
  Init();
  SetName(name);
  SetSplinePar(sfitmode,sopt,svalbeg,svalend,serror);//SplinePar
  SetNode(Nnodes,nodes,ynodes);
}

///---
void QSplineFit::SetName(const char *name){
  qname=name;
  int i=0;
  while (gDirectory->FindObject(Form("%s_fun",qname))!=0){
    qname=Form("%s_f%d",name,i);
    i++;
  }
  Clear();//Clear spline+qfun 
}

///---
void QSplineFit::SetHisto(TH1 *h){
  qh=h;
  qgr=0;
  Clear();//Clear spline+qfun 
}

///---
void QSplineFit::SetGraph(TGraphErrors *gr){
  qgr=gr;
  qh=0;
  Clear();
}

///-----
void QSplineFit::SetNode(int Nnodes,double *nodes,double *ynodes){
  if(qxnodes) delete [] qxnodes;
  if(qynodes) delete [] qynodes;
  qnnodes=Nnodes;
  qxnodes=new double[Nnodes];
  qynodes=new double[Nnodes]; 
  for(int i=0;i<qnnodes;i++){
    qxnodes[i]=nodes[i];
    if(ynodes!=0)qynodes[i]=ynodes[i];
  } 
  if(ynodes!=0){qh=0;qgr=0;}//Fix-Nodes
  Clear();//Clear spline+qfun
};

///-----
int  QSplineFit::GetNodeX(double *xnodes){
  for(int i=0;i<qnnodes;i++){
     xnodes[i]=qxnodes[i];
  }
  return qnnodes; 
}

///-----
void QSplineFit::SetSplinePar(int sfitmode,const char *sopt,double svalbeg,double svalend,double serror){
  _sopt=sopt;
  _svalbeg=svalbeg;
  _svalend=svalend;
  _serror=serror;
  _sfitmode=sfitmode;
   Clear();
}

///---
void QSplineFit::SetVerbose(int verbose){
  _verbose=verbose;
}

///--
double QSplineFit::operator()(double *x, double *par){
//---LogX
   double xv=x[0];
   if(_sfitmode&(LogX|LogX2)){xv=log(fabs(x[0]));}
//---
   int update=0;
   if(qspline==0)update=1;
   else {
      for(int ino=0;ino<qnnodes;ino++){//Nodes
         double parn=par[ino];
         if(_sfitmode&LogY){parn=log(fabs(parn));}
         if(fabs(qynodes[ino]-parn)>_serror*fabs(qynodes[ino])){update=1;break;}
      }
   }
   if(update){
      if(qspline)delete qspline;
      double *qxnodesn=new double[qnnodes];
      for(int ino=0;ino<qnnodes;ino++){
        qxnodesn[ino]=qxnodes[ino];
        qynodes[ino]=par[ino]; 
        if(_sfitmode&(LogX|LogX2))qxnodesn[ino]=log(fabs(qxnodes[ino]));
        if(_sfitmode&LogY)qynodes[ino]=log(fabs(par[ino]));
      }
      qspline = new TSpline3(Form("%s_spline",qname),qxnodesn,qynodes,qnnodes,_sopt,_svalbeg,_svalend);
      delete [] qxnodesn;
  }
//---Add Flux Extrapolate 
  if(_sfitmode&ExtrapolateFlux){
    int    isx0=0;
    double x0=xv;
    if     (xv>log(fabs(qxnodes[qnnodes-1]))){x0=log(fabs(qxnodes[qnnodes-1]));isx0=2;}//End
    else if(xv<log(fabs(qxnodes[0])))        {x0=log(fabs(qxnodes[0]));isx0=1;}//Begin
    double kx=qspline->Derivative(x0);
    if(_verbose>=1){
       std::cout<<"xv="<<x[0]<<" log(xv)="<<xv<<" index="<<kx<<std::endl;
    } 
    double y0=qspline->Eval(x0); 
    double b=y0-kx*x0;
    if(isx0)return exp(kx*xv+b);
  }
 
//---Low+High Flux Linear Extrapolate
 for(int ilh=0;ilh<2;ilh++){
    bool islinear=0;
    if(ilh==0)islinear=(_sfitmode&ExtrapolateLB);
    else      islinear=(_sfitmode&ExtrapolateLE);
    if(!islinear)continue;
    double x0=0,y0=0;
    if(ilh==0)qspline->GetKnot(0,x0,y0);
    else      qspline->GetKnot(qspline->GetNp()-1,x0,y0);
    double kx=qspline->Derivative(x0);
    double b=y0-kx*x0;
    double yv=0;
    if     (ilh==0&&xv<x0){yv=kx*xv+b;}
    else if(ilh==1&&xv>x0){yv=kx*xv+b;}
    else  continue;
    if(_sfitmode&LogY)yv=exp(yv);
    return yv;
 }
//---

//---End
  double yv=qspline->Eval(xv);
  if(_sfitmode&LogY)yv=exp(qspline->Eval(xv));
  return yv;
}
 
///---
TF1 *QSplineFit::GetFun(double xmin,double xmax){
  if(qfun){
    if(xmin!=0||xmax!=0)qfun->SetRange(xmin,xmax);
    return qfun;
  }
  qfun=new TF1(Form("%s_fun",qname),this,xmin,xmax,qnnodes);
  for(int ino=0;ino<qnnodes;ino++){
    if     (qgr)qfun->SetParameter(ino,qgr->Eval(qxnodes[ino]));//Graph First Order
    else if(qh) {
       double y0=qh->Interpolate(qxnodes[ino]);
       if(_sfitmode&LogY){if(y0<1e-12)y0=1e-12;}
       qfun->SetParameter(ino,y0);//Histo Second Order
     }
    else        qfun->SetParameter(ino,qynodes[ino]);//ynodes last
  }
  qfun->SetNpx(100000);
  return qfun;
}

///--
TGraphErrors* QSplineFit::HistoToGraph(int cb){
   TGraphErrors *gr=new TGraphErrors();
   for(int ib=1;ib<=qh->GetNbinsX();ib++){
     double xv=qh->GetXaxis()->GetBinCenterLog(ib);//sqrt(lb*hb)
     if(cb==1)xv=qh->GetXaxis()->GetBinCenter(ib);//(lb+hb)/2
     double yv=qh->GetBinContent(ib);
     double ey=qh->GetBinError(ib);
     int ip=gr->GetN();
     gr->SetPoint(ip,xv,yv);
     gr->SetPointError(ip,0,ey);
  }
  return gr;
}
  
///--
TFitResultPtr QSplineFit::DoFit(Option_t *option ,Option_t *goption, Double_t xxmin, Double_t xxmax){
  if(qgr){
    return qgr->Fit(GetFun(xxmin,xxmax),option,goption,xxmin,xxmax);
  }
  else if(_sfitmode&LogX2){
    qgr=HistoToGraph(0);//Convert to Graph sqrt(lr*hr)
    TFitResultPtr fitresult=qgr->Fit(GetFun(xxmin,xxmax),option,goption,xxmin,xxmax);
    std::string opt(option),gopt(goption);
    if(opt.find("+")==std::string::npos)qh->GetListOfFunctions()->Clear();
//h->GetListOfFunctions()->Detele()
//h.GetFunction("myFunction")->Delete();
//h.GetFunction("myFunction")->SetBit(TF1::kNotDraw);
//--
    qh->GetListOfFunctions()->Add(GetFun(xxmin,xxmax));//Add To result Histo
    (TVirtualFitter::GetFitter())->SetObjectFit(qh);
    delete qgr;qgr=0;
    return fitresult; 
  }
  else { 
    return qh->Fit(GetFun(xxmin,xxmax),option,goption,xxmin,xxmax);
  }
}


int QSplineFit::FillFlux(TH1 *hflux,double scaleindex,double scalenorm){
//---
 for(int ib=1;ib<=hflux->GetNbinsX();ib++){
   double xv=hflux->GetXaxis()->GetBinCenterLog(ib);
   double yv=GetFun()->Eval(xv)*pow(xv,scaleindex)*scalenorm;
   hflux->SetBinContent(ib,yv);
   hflux->SetBinError(ib,0);
 }
 return 0;
}

#endif
