#ifndef KLNfunc_h
#define KLNfunc_h
#include "KLNModel.h"
#include "UnintegPartonDist.h"

class KLNfunc: public UnintegPartonDist
{
private:
    double fac;
public:
    KLNfunc() {
      //fac = 1.0/(2*M_PI*M_PI/KLNModel::CF)/M_PI;
      fac = KLNModel::CF *2./(3.*M_PI*M_PI);
    }
    double getFunc(double qs2, double x, double kt2, double alp) {
      if(kt2 <= qs2) return fac/alp;
      return fac*qs2/kt2/alp;
    }
    double getFuncNF(double qs2, double x, double kt2) {return 0.;}
    double getQs(double qs0_2, double x) {return 0.;}
    double getMax_kt() {return 1e38;}
};
#endif
