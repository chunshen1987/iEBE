
#include "Largex.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "KLNModel.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


using namespace std;

Large_x::Large_x(KLNModel* klnX, UnintegPartonDist* rcBKugdX)
{
  kln = klnX;
  rcBKugd = rcBKugdX;

  wavefunc = kln->getwavefunc(); // get pointer to ugd
  dNdeta = kln->probeEtaFlag(); // read out flags that have already been set
  dEtdy  = kln->probeEtFlag();  // for small-x part
  // we treat large-x (x>xCut) partons here, not in kt factorization part
  // wavefunc->setlgXflag(0);

  kln->getPtRange(&Ptmin,&Ptmax);
  if (Ptmax > 0.66*rcBKugd->getMax_kt())
    cout << "#  WARNING in Large_x(): max hadron pt=" << Ptmax <<
      ",  max gluon kt=" << rcBKugd->getMax_kt() << endl;

  ecm = kln->getEcm();
  mHadron = kln->getMHadron();  // mass scale for eta <--> y transformation
  Norm = 2.5;   // initialize normalization factor for large-x contribution
  DonePtDistr = 0;  // flag for pt distribution
  ptspline = 0;  ptacc=0;  // spline pointers

  // for Gaussian integration
  xg = new double [38];
  wg = new double [38];
  OverLap::Gauss38(0.0,1.0,xg,wg);

}


Large_x::~Large_x()
{
  delete xg;
  delete wg;
  if (ptspline) gsl_spline_free(ptspline);
  if (ptacc)	gsl_interp_accel_free(ptacc);
}



// large-x (x>x0) production
//    returns dN/dyd2r [1/fm^2] summed over parton flavors
double Large_x::getdNdy(double y, double Npart1, double Npart2)
{
  if (kln->getModel() < rcBK) {
    cout << "WARNING: large-x parton production only in rcBK mode !\n" ;
    return 0.;
  }
  rapidity=y;

  double trEtaY=1.0;
  if(dNdeta<0) {  // global (pt-independent) y->eta Jacobian
    double avg_pt = 0.13 + .32*pow(200./1000.,0.115);
      // this call modifies variable "rapidity" from eta to y !
    trEtaY=kln->getJacobian(rapidity,avg_pt*avg_pt);
  }

  // now integrate over pt to obtain dN/d2r dy in forw/backw regions
  int N2=38;  // # of integration points
  double dN=0.;
  for(int ipt=0; ipt<N2; ipt++) { // sum over pt
    double pt=Ptmin+xg[ipt]*(Ptmax-Ptmin);
    double dpt=wg[ipt]*(Ptmax-Ptmin);

    transEtaY=1.0;
    double mt=pt;
    double rap_save = rapidity;
    if(dNdeta>0) {
      double pt2 = pt*pt;
      // this call modifies variable "rapidity" from eta to y !
      transEtaY=kln->getJacobian(rapidity,pt2);
      mt=sqrt(pt2+mHadron*mHadron);
    }
    double x1,x2,Q0sq;
    if (rapidity<0.0) {
      x1=mt/ecm/exp(rapidity);
      x2=mt/ecm*exp(rapidity);
      // get Qs(x0)^2 of target (right-mover)
      Q0sq = kln->SaturationScale(x2,Npart1);
    } else {
      x1=mt/ecm*exp(rapidity);
      x2=mt/ecm/exp(rapidity);
      // get Qs(x0)^2 of target (left-mover)
      Q0sq = kln->SaturationScale(x2,Npart2);
    }
    rapidity = rap_save;
    // cut on x1>x0, x2<x0
    if (x1 <= UnintegPartonDist::xCut || x2 >= UnintegPartonDist::xCut ||
	x1 > 1.0 || x2 > 1.0) continue;

    if (dEtdy)
      dN += DHJfunc(x1,x2,mt,Q0sq) *pt*dpt*transEtaY*mt;
    else
      dN += DHJfunc(x1,x2,mt,Q0sq) *pt*dpt*transEtaY;
  } // sum over pt

  // # of nucl / sigma0 in right / left mover
  dN *= (rapidity>0) ? Npart1 : Npart2;
  dN *= Norm * 2.*M_PI;  // 2pi from integration over azimuth of pt
  return (dN*trEtaY);
}


// maps (0=g, 1=u, 2=d, 3=s) flavor coding to KKP (see kkp.cxx)
#define KKPflavors(flav)  (flav ? (2*flav-1) : flav)

//    returns dN/dyd2ptd2r [1/(GeV fm)^2] for charged hadrons
// (ad. & YN)
//   ATTN: you need to call lgXClearPtTable() before calling with new
//         values for y, Npart1, Npart2 !!
double Large_x::getdNdyd2pt(double y, double Npart1, double Npart2, double pt)
{
  if (pt > Ptmax || pt < Ptmin) return 0.;
  if (DonePtDistr)   // hadron pt distribution is available
    return gsl_spline_eval(ptspline, pt, ptacc);

  // ... if not, need to calculate it
  if (kln->getModel() < rcBK) {
    cout << "WARNING: large-x parton production only in rcBK mode !\n" ;
    return 0.;
  }
  rapidity=y;

  double trEtaY=1.0;
  if(dNdeta<0) {  // global (pt-independent) y->eta Jacobian
    double avg_pt = 0.13 + .32*pow(200./1000.,0.115);
      // this call modifies variable "rapidity" from eta to y !
    trEtaY=kln->getJacobian(rapidity,avg_pt*avg_pt);
  }

  for (int iqt=0; iqt<40; iqt++) {
    for (int flav=0; flav<4; flav++)
      PartPtDistr[flav][iqt] = 0.;
    // init parton pt points
    ptPoints[iqt] = Ptmin + iqt*(rcBKugd->getMax_kt()-Ptmin)/40.;
  }

  // generate parton qt-distributions, each flavor individually
  for (int iqt=0; iqt<40; iqt++) {  // loop over parton qt
    double x1,x2,Q0sq,pt_p;
    pt_p = ptPoints[iqt];
    if (rapidity<0.0) {
      x1=pt_p/ecm/exp(rapidity);
      x2=pt_p/ecm*exp(rapidity);
      // get Qs(x0)^2 of target (right-mover)
      Q0sq = kln->SaturationScale(x2,Npart1);
    } else {
      x1=pt_p/ecm*exp(rapidity);
      x2=pt_p/ecm/exp(rapidity);
      // get Qs(x0)^2 of target (left-mover)
      Q0sq = kln->SaturationScale(x2,Npart2);
    }
    // cut on large x1, small x2
    if (x1 <= UnintegPartonDist::xCut) continue;
    if (x1 >= 1.0 || x2 >= UnintegPartonDist::xCut) break;

    for (int flav=0; flav<4; flav++)    // loop over parton flavor
      PartPtDistr[flav][iqt] = DHJfunc(x1,x2,pt_p,Q0sq,flav);
  }  // loop over parton qt
  // --- parton qt distributions/tables done ---

  // --- now generate hadron pt tables & GSL splines ---
  for (int ipt_h=0; ipt_h<40; ipt_h++)
    HadrPtDistr[ipt_h]=0.;

  for (int flav=0; flav<4; flav++) {   // loop over parton flavor
    ptacc = gsl_interp_accel_alloc();
    ptspline = gsl_spline_alloc(gsl_interp_cspline, 40);
    // init GSL spline for parton dN/dyd2pt
    gsl_spline_init (ptspline, ptPoints, &PartPtDistr[flav][0], 40);

    for (int ipt_h=0; ipt_h<40; ipt_h++) {  // loop over hadron pt
      double pt_h = Ptmin + ipt_h*(Ptmax-Ptmin)/40.;
      // integrate over z to obtain dN_h/d2pt dy
      double zmax=1.0, zmin;
      zmin=pt_h/rcBKugd->getMax_kt();
      if (zmin<0.05) zmin = 0.05; // very small z may violate mom. sum rule
      if (zmin >= zmax) break;
      int N2=38;  // # of integration points
      for(int iz=0; iz<N2; iz++) { // sum over z
	double z=zmin+xg[iz]*(zmax-zmin);
	double qt = pt_h/z;
	if (qt>rcBKugd->getMax_kt()) continue;
	double dz=wg[iz]*(zmax-zmin);
	HadrPtDistr[ipt_h] += dz/z/z
	  * kln->Dgl(KKPflavors(flav), z, qt)
	  * gsl_spline_eval(ptspline, qt, ptacc);
      }
    }
    gsl_spline_free(ptspline);
    gsl_interp_accel_free(ptacc);
    ptspline = 0;  ptacc=0;  // clear spline pointers
  }

  for (int ipt_h=0; ipt_h<40; ipt_h++) {
    ptPoints[ipt_h] = Ptmin + ipt_h*(Ptmax-Ptmin)/40.;

    // # of nucl / sigma0 in right / left mover
    HadrPtDistr[ipt_h] *= (rapidity>0) ? Npart1 : Npart2;
    HadrPtDistr[ipt_h] *= Norm;  // apply "K factor" for large-x contribution
    HadrPtDistr[ipt_h] *= trEtaY;  // eta <--> y Jacobian
  }

  ptacc = gsl_interp_accel_alloc();
  ptspline = gsl_spline_alloc(gsl_interp_cspline, 40);
  // init GSL spline for hadron dN/dyd2pt
  gsl_spline_init (ptspline, ptPoints, HadrPtDistr, 40);
  DonePtDistr = 1;
  return gsl_spline_eval(ptspline, pt, ptacc);
}



// DHJ formula for parton production; from eq.(22) of hep-ph/0506308
//  returns dsigma/dyd2ptd2r [1/GeV^2] summed over parton flavors
double Large_x::DHJfunc(double x1, double x2, double mt, double Q0sq, int flav)
{
  double Q2 = wavefunc->getQs(Q0sq,x2);  Q2 = Q2*Q2;
  double pt2 = mt*mt;
  if (pt2 > Q2) Q2 = pt2;     // Q^2 = max(pt^2 , Qs^2)

  switch (flav) {
  case -1 : // sum all partons
    return gDHJfunc(x1,x2,pt2,Q2,Q0sq)+qDHJfunc(x1,x2,pt2,Q2,Q0sq);
  case 0 : return gDHJfunc(x1,x2,pt2,Q2,Q0sq);  // g
  case 1 :
  case 2 :
  case 3 : return qDHJfunc(x1,x2,pt2,Q2,Q0sq,flav);  // q+qbar, flavor "flav"
  default :   // ?
    cout << "# WARNING in Large_x::DHJfunc(), unknown parton flavor "
	 << flav << ", ignoring...\n";
    return 0.0;
  }
}


// dsigma/dyd2ptd2r [1/GeV^2] for quarks
double Large_x::qDHJfunc(double x1, double x2, double pt2, double Q2, double Q0sq, int flav)
{
  double dNdyd2pt = 
    grv94l_q(x1,Q2,flav) * wavefunc->getFuncNF(Q0sq,x2,pt2);
  return dNdyd2pt/pow(2.*M_PI,2.)*pow(1.-x2,4.);
}


// dsigma/dyd2ptd2r [1/GeV^2] for gluons
double Large_x::gDHJfunc(double x1, double x2, double pt2, double Q2, double Q0sq)
{
  double dNdyd2pt = 
    // note: N_A(pt) = Phi(pt) * alpha_s * 6pi^3/pt^2  for rcBK()
    grv94l_g(x1,Q2) * wavefunc->getFunc(Q0sq,x2,pt2,pt2/(6.*M_PI*M_PI*M_PI));
  return dNdyd2pt/pow(2.*M_PI,2.)*pow(1.-x2,4.);
}




// ******** DGLAP PDFs ********

//...GRV 94 L (leading order) quark (val+sea) distribution function set
//...in parametrized form; xq(x,Q2)
//...Authors: M. Glueck, E. Reya and A. Vogt.
double Large_x::grv94l_q(double x, double q2, int flav)
{
 
//...Common expressions.
    double mu2  = 0.23;
    double lam2 = 0.2322 * 0.2322;
    double s  = log (log(q2/lam2) / log(mu2/lam2));
    double  ds = sqrt (s);
    double  s2 = s * s;
    double  s3 = s2 * s;
 
//    UV  = U(VAL) = U - U(BAR)
//    DV  = D(VAL) = D - D(BAR)
//    DEL = D(BAR) - U(BAR)
//    UDB = U(BAR) + D(BAR) 
//    SB  = S = S(BAR)
//    GL  = GLUON
//...uv :
      double nu  =  2.284 + 0.802 * s + 0.055 * s2;
      double aku =  0.590 - 0.024 * s;
      double bku =  0.131 + 0.063 * s;
      double au  = -0.449 - 0.138 * s - 0.076 * s2;
      double bu  =  0.213 + 2.669 * s - 0.728 * s2;
      double cu  =  8.854 - 9.135 * s + 1.979 * s2;
      double du  =  2.997 + 0.753 * s - 0.076 * s2;
      double uv  = grvv (x, nu, aku, bku, au, bu, cu, du);
 
//...dv :
      double nd  =  0.371 + 0.083 * s + 0.039* s2;
      double akd =  0.376;
      double bkd =  0.486 + 0.062 * s;
      double ad  = -0.509 + 3.310 * s - 1.248 * s2;
      double bd  =  12.41 - 10.52 * s + 2.267 * s2;
      double cd  =  6.373 - 6.208 * s + 1.418 * s2;
      double dd  =  3.691 + 0.799 * s - 0.071 * s2;
      double dv  = grvv (x, nd, akd, bkd, ad, bd, cd, dd);
 
//...del :
      double ne  =  0.082 + 0.014 * s + 0.008 * s2;
      double ake =  0.409 - 0.005 * s;
      double bke =  0.799 + 0.071 * s;
      double ae  = -38.07 + 36.13 * s - 0.656 * s2;
      double be  =  90.31 - 74.15 * s + 7.645 * s2;
      double ce  =  0.0;
      double de  =  7.486 + 1.217 * s - 0.159 * s2;
      double del = grvv (x, ne, ake, bke, ae, be, ce, de);
 
//...udb :
      double alx =  1.451;
      double bex =  0.271;
      double akx =  0.410 - 0.232 * s;
      double bkx =  0.534 - 0.457 * s;
      double agx =  0.890 - 0.140 * s;
      double bgx = -0.981;
      double cx  =  0.320 + 0.683 * s;
      double dx  =  4.752 + 1.164 * s + 0.286 * s2;
      double ex  =  4.119 + 1.713 * s;
      double esx =  0.682 + 2.978 * s;
      double udb = grvw (x, s, alx, bex, akx, bkx, agx, bgx, cx, dx, ex, esx);
 
//...sb :
      double sts =  0.0;
      double als =  0.914;
      double bes =  0.577;
      double aks =  1.798 - 0.596 * s;
      double as  = -5.548 + 3.669 * ds - 0.616 * s;
      double bs  =  18.92 - 16.73 * ds + 5.168 * s;
      double dst =  6.379 - 0.350 * s + 0.142 * s2;
      double est =  3.981 + 1.638 * s;
      double ess =  6.402;
      double sb  = grvs (x, s, sts, als, bes, aks, as, bs, dst, est, ess);
 
      switch (flav) {
      case 0 : return uv+dv+2.*udb+2.*sb;  // u+ubar+d+dbar+s+sbar
      case 1 : return uv+udb-del;  // u+ubar
      case 2 : return dv+udb+del;  // d+dbar
      case 3 : return 2.*sb;  // s+sbar
      default : return 0.;
      }
}


//...GRV 94 L (leading order) gluon distribution function set
//...in parametrized form; xg(x,Q2) 
//...Authors: M. Glueck, E. Reya and A. Vogt.
double Large_x::grv94l_g(double x, double q2)
{
 
//...Common expressions.
    double mu2  = 0.23;
    double lam2 = 0.2322 * 0.2322;
    double s  = log (log(q2/lam2) / log(mu2/lam2));
    double  ds = sqrt (s);
    double  s2 = s * s;
    double  s3 = s2 * s;
 
//...gl :
      double alg =  0.524;
      double beg =  1.088;
      double akg =  1.742 - 0.930 * s;
      double bkg =                         - 0.399 * s2;
      double ag  =  7.486 - 2.185 * s;
      double bg  =  16.69 - 22.74 * s  + 5.779 * s2;
      double cg  = -25.59 + 29.71 * s  - 7.296 * s2;
      double dg  =  2.792 + 2.215 * s  + 0.422 * s2 - 0.104 * s3;
      double eg  =  0.807 + 2.005 * s;
      double esg =  3.841 + 0.316 * s;
      double gl  = grvw (x, s, alg, beg, akg, bkg, ag, bg, cg, dg, eg, esg);

      return gl;
}


//*********************************************************************
//...Auxiliary for the GRV 94 parton distribution functions
//...for u and d valence and d-u sea.
//...Authors: M. Glueck, E. Reya and A. Vogt.
inline double Large_x::grvv (double x, double n, double ak, double bk,
       	double a, double b, double c, double d)
{
    double dx = sqrt (x);
    return n * pow(x,ak) * (1.0+ a*pow(x,bk) + x * (b + c*dx)) * pow(1.0- x,d);
}
 
//...Auxiliary for the GRV 94 parton distribution functions
//...for d+u sea and gluon.
//...Authors: M. Glueck, E. Reya and A. Vogt.
inline double Large_x::grvw (double x, double s, double al, double be, double ak,
       	double bk, double a, double b, double c, double d, double e, double es)
{
    double lx = log (1.0/x);
    return (pow(x,ak) * (a + x * (b + x*c)) * pow(lx,bk) + pow(s,al)
           * exp (-e + sqrt (es * pow(s,be) * lx))) * pow(1.0- x,d);
}
 
//...Auxiliary for the GRV 94 parton distribution functions
//...for s, c and b sea.
//...Authors: M. Glueck, E. Reya and A. Vogt.
inline double Large_x::grvs (double x, double s, double sth, double al,
       double be, double ak, double ag, double b, double d, double e, double es)
{
 
    if(s <= sth) return 0.0;

    double dx = sqrt (x);
    double lx = log (1.0/x);
    return pow(s - sth,al) / pow(lx,ak) * (1.0+ ag*dx + b*x) *
           pow(1.0- x,d) * exp (-e + sqrt (es * pow(s,be) * lx));
 
}
