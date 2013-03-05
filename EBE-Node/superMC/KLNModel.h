#ifndef KLNModel_h
#define KLNModel_h

#include "ParamDefs.h"
#include "Bases.h"
#include "UnintegPartonDist.h"
#include "OverLap.h"
#include <cmath>

//... kt factorization model for small-x gluon production

class KLNModel : public Bases
{
private:
  double siginNN, siginNN200;
  UnintegPartonDist* wavefunc;
  double ecm;
  int    model;
  int    dNdeta;       // =1: calculate dN/deta
  int    dEtdy;        // =1: calculate dE/dy
  double partonPt;    // produce parton with fixed pt
  double g2hfac;      // g-->hadr conversion factor
  double alphaS;      // alpha strong for fixed coupling case
  double lambda;
  double lambdaQCD;
  double lambdaQCD2;
  double mHadron;
  double Ptmin, Ptmax;
  double Npart1, Npart2;
  int    optFixAlpha;

  double rapidity;
  double transEtaY;


 public:
  static double Nc;
  static double Nf;
  static double CF;
  static double FacSQ;
  static double Norm;
  static double Beta0;

 public:
  KLNModel(double srt,int mode,UnintegPartonDist* f);
  virtual ~KLNModel();

  double getdNdy(double y,double npart1,double npart2, double pt=-1.0);
  inline double getAlphaStrong(const double q2);

  double SaturationScale(double x,double npart);
  UnintegPartonDist* getwavefunc() {return wavefunc;}
  double getJacobian(double h, double pt2);
  void FFconv(int iFF, int bins, double ptmin, double dpt, double *dNdyd2pt);

  void setNormalization(double c)      {g2hfac *= c;}
  void setg2hfac(double c)             {g2hfac *= c;}
  int  getModel()                      {return model;}
  void setdNdeta(int i)                {dNdeta=i;}
  int  probeEtaFlag()                  {return dNdeta;}
  void setdEtdy(int i)                 {dEtdy=i;}
  int  probeEtFlag()                   {return dEtdy;}
  void setLambdaQCD(double l)          {lambdaQCD=l; lambdaQCD2=l*l;}
  double getLambdaQCD()                {return lambdaQCD;}
  void setLambdaSaturation(double l)   {lambda = l;}
  void getPtRange(double *a, double *b){*a=Ptmin; *b=Ptmax;}
  void setPtRange(double a, double b)  {Ptmin=a; Ptmax=b;}
  double getMHadron()                  {return mHadron;}
  double getEcm()                      {return ecm;}
  void setOptFixAlpha(int i)           {optFixAlpha=i;}
  double Dgl(int iFF, double z, double Q);

 protected:
  double ktF_MCintegral();
  double waveFunction(double qssq,double x,double kt2);
  double func(double* x);

 private:
  double pow__dd(double *x, double *y) {return pow(*x,*y);}
  void kkpFF(int *ih, int *iset, double *x, double *qs, double *dh);
};


inline double KLNModel::getAlphaStrong(const double q2)
{
  if(optFixAlpha) return alphaS;
  if(q2 <= lambdaQCD2) return alphaS;
  return std::min(alphaS,  1.0/( Beta0 * log( q2/lambdaQCD2 ) ) );
}

#endif
