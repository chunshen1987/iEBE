#ifndef rcBKfunc_h
#define rcBKfunc_h

#include <fstream>
#include <iostream>
#include <cstdlib>
#include "ParamDefs.h"
#include "UnintegPartonDist.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>



class rcBKfunc: public UnintegPartonDist
{
private:
  std::ifstream ifs;
  double *table;    // ugd Phi(x,kt)
  double *tableNF;  // N_F(x,kt)   (fund. dipole)
  double *ktPoint;  // index of kt
  double *QsTab;    // table of Qs(x;Q0)
  int maxQ0; // number of Qs_0 bins, initialized in rcBKfunc.cxx
  int maxY;  // # of rapidity bins
  int maxKt;  // # of kt bins
  double UGD_maxKt;  // max kt[GeV] for UGD table
  double dQ0; // difference of Qs_0^2 from one table to the next
  int uGDflag;
  gsl_interp_accel ***acc;
  gsl_spline ***spline;
  gsl_interp_accel ***accNF;
  gsl_spline ***splineNF;


 public:
  rcBKfunc(int uGDmodel);
  
  ~rcBKfunc() {
    delete table;
    delete ktPoint;
    delete QsTab;
    delete tableNF;
    for(int iq=0;iq<maxQ0;iq++) {
      for(int iy=0;iy<maxY;iy++) {
	gsl_spline_free(spline[iq][iy]);
	gsl_spline_free(splineNF[iq][iy]);
	gsl_interp_accel_free(acc[iq][iy]);
	gsl_interp_accel_free(accNF[iq][iy]);
      }
      delete acc[iq];  delete accNF[iq];
      delete spline[iq]; delete splineNF[iq];
    }
    delete acc; delete spline;
    delete accNF; delete splineNF;
  }

  double getMax_kt() {return UGD_maxKt;}

  double getTable(int iq, int iy, int ik) {return table[ik + maxKt*(iy + maxY*iq)];}
  double getKt(int iq, int iy, int ik)  {return ktPoint[ik + maxKt*(iy + maxY*iq)];}
  double getQsTab(int iq, int iy)  {return QsTab[iy + maxY*iq];}


  // --- returns rcBK unintegrated gluon distribution ---
  double getFunc(double qs0_2, double x, double kt2, double alp) {

    if (x<0. || x>1. || qs0_2<0) {
      std::cout << "WARNING: x or Qs_0^2 out of bounds!\n";
      return 0.;
    }

    // get the rcBK-uGD from the table.
    
    double dY=0.1;     // rapidity evolution step
    double Y = log(x0/x);
    if (Y<0.0) {
      // no large-x (x>xCut) extrapolation
      if (probelgXflag() == 0 && x>xCut) return 0.;
      // extrap. to Y<0 assuming Qs^2(Y) = Q0^2 exp(lambda Y)
      qs0_2 *= exp(lgXlambda*Y);
      Y=0.;  // ini. cond. is MV model
    }

    int iy= (int) (Y/dY+.5);  // round to nearest table entry
    if (iy>=maxY) {
      std::cout << "# WARNING: x too small, Y overflow !" << std::endl;      
      iy=maxY-1;
    }

    // convert to Q0 for quark, tables labelled by Q0 for fund rep
    double Q02 = 4./8.*qs0_2;
    // find table entry just below actual Q0^2
    int iq= (int) (Q02/dQ0); 
    iq -= 1;
    int iqoffset =1;
    if (uGDflag==rcBKalbacete) {
      iq -= 1;  // first table corr. to Q0^2 = 0.2GeV^2  <--> iq=0
      iqoffset++;
    }
    if (iq == maxQ0-1) iq--;  // right on last table
    else if (iq > maxQ0-1) {
      std::cout << "# WARNING: Q0=" << sqrt(Q02) 
		<< "  overflow, performing extrapolation !" << std::endl;
      iq=maxQ0-2;
    }

    // lin. interpolation in Q0^2
    double val1 = 0.;
    double val2 = 0.;
    // this defines Phi(x,kt) = N_A(x,pt) * C_F * kt^2 / (2pi)^3 alpha(kt)
    double fac = kt2/(6.*M_PI*M_PI*M_PI)/alp;
    if (iq>=0) {
      val1 = gsl_spline_eval(spline[iq][iy], sqrt(kt2), acc[iq][iy]);
      val2 = gsl_spline_eval(spline[iq+1][iy], sqrt(kt2), acc[iq+1][iy]);
      return  (val1+(val2-val1)*(Q02-(iq+iqoffset)*dQ0)/dQ0) *fac;
    }
    else { // interpolate between Q0^2 = 0 and first table
      val2 = gsl_spline_eval(spline[0][iy], sqrt(kt2), acc[0][iy]);
      return  (val2*Q02/(iqoffset*dQ0)) *fac;
    }      
  }



  // --- returns N_F(kt) ---
  double getFuncNF(double qs0_2, double x, double kt2) {

    if (x<0. || x>1. || qs0_2<0) {
      std::cout << "WARNING: x or Qs_0^2 out of bounds!\n";
      return 0.;
    }

    // get N_F from the table.
    
    double dY=0.1;     // rapidity evolution step
    double Y = log(x0/x);
    if (Y<0.0) {
      // extrap. to Y<0 assuming Qs^2(Y) = Q0^2 exp(lambda Y)
      qs0_2 *= exp(lgXlambda*Y);
      Y=0.;  // ini. cond. is MV model
    }

    int iy= (int) (Y/dY+.5);  // round to nearest table entry
    if (iy>=maxY) {
      std::cout << "# WARNING: x too small, Y overflow !" << std::endl;      
      iy=maxY-1;
    }

    // convert to Q0 for quark, tables labelled by Q0 for fund rep
    double Q02 = 4./8.*qs0_2;
    // find table entry just below actual Q0^2
    int iq= (int) (Q02/dQ0); 
    iq -= 1;
    int iqoffset =1;
    if (uGDflag==rcBKalbacete) {
      iq -= 1;  // first table corr. to Q0^2 = 0.2GeV^2  <--> iq=0
      iqoffset++;
    }
    if (iq == maxQ0-1) iq--;  // right on last table
    else if (iq > maxQ0-1) {
      std::cout << "# WARNING: Q0=" << sqrt(Q02) 
		<< "  overflow, performing extrapolation !" << std::endl;
      iq=maxQ0-2;
    }

    // lin. interpolation in Q0^2
    double val1 = 0.;
    double val2 = 0.;
    if (iq>=0) {
      val1 = gsl_spline_eval(splineNF[iq][iy], sqrt(kt2), accNF[iq][iy]);
      val2 = gsl_spline_eval(splineNF[iq+1][iy], sqrt(kt2), accNF[iq+1][iy]);
      return  (val1+(val2-val1)*(Q02-(iq+iqoffset)*dQ0)/dQ0);
    }
    else { // interpolate between Q0^2 = 0 and first table
      val2 = gsl_spline_eval(splineNF[0][iy], sqrt(kt2), accNF[0][iy]);
      return  (val2*Q02/(iqoffset*dQ0));
    }      
  }




  // returns Qs[GeV] at x for given Qs(x0) as obtained from max of N_F*kt2
  //    ATTN: in most applications, this would correspond to the TARGET-x
  //    for production of a large-x_p _projectile_
  //          parton with rapidity y and tr. mom. pt:  
  //          x = x_p / exp(2y) = pt^2/(x_p s)

  double getQs(double qs0_2, double x) {

    if (x<0. || x>1. || qs0_2<0) {
      std::cout << "WARNING: Qs_0^2 out of bounds in getQs() !\n";
      return 0.;
    }

    double dY=0.1;     // rapidity evolution step
    double Y = log(x0/x);
    if (Y<0.0) {
      // extrap. to Y<0 assuming Qs^2(Y) = Q0^2 exp(lambda Y)
      qs0_2 *= exp(lgXlambda*Y);
      Y=0.;  // ini. cond. is MV model
    }

    int iy= (int) (Y/dY+.5);  // round to nearest table entry
    if (iy>=maxY) {
      std::cout << "# WARNING: x too small, Y overflow !" << std::endl;      
      iy=maxY-1;
    }

    // convert to Q0 for quark, tables labelled by Q0 for fund rep
    double Q02 = 4./8.*qs0_2;
    // find table entry just below actual Q0^2
    int iq= (int) (Q02/dQ0); 
    iq -= 1;
    int iqoffset =1;
    if (uGDflag==rcBKalbacete) {
      iq -= 1;  // first table corr. to Q0^2 = 0.2GeV^2  <--> iq=0
      iqoffset++;
    }
    if (iq == maxQ0-1) iq--;  // right on last table
    else if (iq > maxQ0-1) {
      std::cout << "# WARNING: Q0=" << sqrt(Q02) 
		<< "  overflow, performing extrapolation !" << std::endl;
      iq=maxQ0-2;
    }

    // lin. interpolation in Q0^2
    double val1 = 0.;
    double val2 = 0.;
    if (iq>=0) {   // Q0^2 equal to first table or larger
      val1 = getQsTab(iq,iy);
      val2 = getQsTab(iq+1,iy);
      return  (val1+(val2-val1)*(Q02-(iq+iqoffset)*dQ0)/dQ0);
    }
    else { // interpolate between Q0^2 = 0 and first table
      val2 = getQsTab(0,iy);
      return  (val2*Q02/(iqoffset*dQ0));
    }      
  }  // end getQs()


};
#endif

