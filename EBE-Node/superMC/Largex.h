#ifndef ForwardJet_h
#define ForwardJet_h

#include "KLNModel.h"
#include <cmath>
#include "UnintegPartonDist.h"
#include "rcBKfunc.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


//... production of large-x partons with x>x0

class Large_x
{
private:
  KLNModel* kln;
  UnintegPartonDist* rcBKugd;
  double Norm;  // normalization factor for large-xF parton production
  int dNdeta;   // flag for dNdeta vs dNdy
  int dEtdy;    // flag for Et vs N
  double rapidity;
  double mHadron;
  double ecm;
  double Ptmin, Ptmax;
  double transEtaY;
  gsl_spline *ptspline;  // for hadron pt distribution
  gsl_interp_accel *ptacc;
  double PartPtDistr[4][40];
  double HadrPtDistr[40], ptPoints[40];
  int DonePtDistr;

  double DHJfunc(double x1, double x2, double mt, double Q0sq, int flav=-1);
  double qDHJfunc(double x1, double x2, double pt2, double Q2, double Q0sq, int flav=0);
  double gDHJfunc(double x1, double x2, double pt2, double Q2, double Q0sq);
  UnintegPartonDist *wavefunc;
  double *xg, *wg;
  double grv94l_q(double x, double q2, int flav=0), grv94l_g(double, double);
  inline double grvv (double x, double n, double ak, double bk,
			     double a, double b, double c, double d);
  inline double grvs (double x, double s, double sth, double al,
			       double be, double ak, double ag, double b, double d, double e, double es);
  inline double grvw (double x, double s, double al, double be, double ak,
		      double bk, double a, double b, double c, double d, double e, double es);


public:
  Large_x(KLNModel* klnX, UnintegPartonDist* rcBKugdX);
  virtual ~Large_x();

  void lgXClearPtTable() {
    DonePtDistr=0;  
    if (ptspline) gsl_spline_free(ptspline);
    if (ptacc)	gsl_interp_accel_free(ptacc);
    ptacc = 0; ptspline = 0;
  }
  void setNormalization(double c)      {Norm=c;}
  double getdNdy(double, double, double);
  double getdNdyd2pt(double, double, double, double);
};


#endif


