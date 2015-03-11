#include <vector> 
#include "RandomVariable.h"

#ifndef HulthenFunc_h
#define HulthenFunc_h

class HulthenFunc: public RandomVariable
{
    private:
       // allocate Hulthen Wave Function's cummulative distribution
        double (HulthenFunc::*pointCDF)(double);
    public: 
      HulthenFunc();
     ~HulthenFunc();
      static double CDF(double);
      double invCDF(double);
      double rand();  // returns random number x governed by Hulthen Function
};

#endif         


