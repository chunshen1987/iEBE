
#include <vector> 
#include "HulthenFunc.h"
#include "arsenal.h"
#include <cmath>
#include <math.h>
#include <time.h>

#define ZERO 1e-15

using namespace std; 

HulthenFunc::HulthenFunc()
{
  double (HulthenFunc::*pointCDF)(double);
}

HulthenFunc::~HulthenFunc()
{
}

double HulthenFunc::rand()
{
	return sampleUsingInvCDF(0,1);
}

double HulthenFunc::CDF(double r)
{
  double alpha = .228;
  double beta = 1.18;
  if (r<=0) return 0;
	/* define CDF function */
  
  const double c = (alpha*beta*(alpha+beta))/((alpha-beta)*(alpha-beta));
  return 2*c*(2*(exp(-r*(alpha+beta))/(alpha+beta))-.5*exp(-2*alpha*r)/alpha-.5*exp(-2*beta*r)/beta+.5/alpha+.5/beta-2/(alpha+beta));
}

double HulthenFunc::invCDF(double x)
{
  return invertFunc(&HulthenFunc::CDF,x,0,100.0,0.001,1.0,0.001);
}
