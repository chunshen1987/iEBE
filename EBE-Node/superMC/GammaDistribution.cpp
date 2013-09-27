#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "arsenal.h"
#include "GammaDistribution.h"

#define ZERO 1e-15

using namespace std;

GammaDistribution::GammaDistribution(double alpha_in, double beta_in)
{  
   resetDistribution(alpha_in, beta_in);
}

GammaDistribution::~GammaDistribution()
{  
}

double GammaDistribution::pdf(double x)
{
  if (x<=0) return 0;
  double result = exp(alpha*log(beta) - lgamma(alpha) + (alpha - 1)*log(x)
                      - beta*x);
  return(result);
}

void GammaDistribution::resetDistribution(double alpha_in, double beta_in)
// reset parameter for pdf with alpha_in and beta_in
// calculate mean, std, and mode for the given pdf
// it then construct the envelop functions for better sampling efficiencies.
{
   if(alpha_in <= 0 || beta_in <=0)
   {
      cout << "Error: parameter of Poisson distribution is negative!" << endl;
      exit(1);
   }
   alpha = alpha_in;
   beta = beta_in;
   mean = alpha/beta;
   std = sqrt(alpha/(beta*beta));
   
   if(alpha > 1)
      mode = (alpha - 1.)/beta;
   else
      mode = 0.0;

   // fill in envelop pdf and inverse CDF
   // first test left boundary
   int nstep = 20;
   int step_left, step_right;
   double step_width = std/2.;
   step_right = nstep;  // no constrain boundary on the right side
   //find the boundary at 0 on the left side
   for (step_left=nstep; step_left>0; step_left--) if (mode-step_width*step_left>=0) break;
   // then fill envelop functions
   constructEnvelopTab(mode, step_width, step_left, step_right);
   return;
}

double GammaDistribution::rand()
{
  return sampleUsingPDFAndEnvelopFunc();
}

double GammaDistribution::rand(double alpha_in, double beta_in)
{
   resetDistribution(alpha_in, beta_in);
   return rand();
}
