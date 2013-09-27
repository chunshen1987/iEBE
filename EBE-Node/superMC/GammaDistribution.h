#include <vector>
#include "RandomVariable.h"

#ifndef GammaDistribution_h
#define GammaDistribution_h

class GammaDistribution: public RandomVariable
{
    private:
       double alpha, beta;  // parameter for Gamma distribution: alpha > 0 shape, beta > 0 rate

    public:
       GammaDistribution(double alpha_in = 1, double beta_in = 1);
       ~GammaDistribution();

       double mean, std;      // mean and standard deviation
       double mode;              // maximum position of pdf
       double pdf(double x);
       void resetDistribution(double alpha_in, double beta_in);
       
       double rand();           // return random integer according to Poisson distribution
       double rand(double alpha_in, double beta_in);    // return random integer according to Poisson distribution with given parameter lambda

};

#endif
