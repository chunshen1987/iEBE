#include "GaussianDistribution.h"
#include <cmath>
#include "arsenal.h"
#include "iostream"

GaussianDistribution::GaussianDistribution(double m, double s)
{
	mu = m;
	sigma = s;
      CDF_array_size = 1000;
      CDF = new vector<double>;
      CDF_x = new vector<double>;
      initCDF();
}

GaussianDistribution::~GaussianDistribution()
{
      delete CDF;
      delete CDF_x;
}


void GaussianDistribution::initCDF()
{
	double dx = 12.*sigma/(CDF_array_size - 1);
	double sum = 0.0;
	for(int i = 0; i < CDF_array_size; i++)
      {
            double x_local = mu - 6*sigma + i*dx;
		sum += pdf(x_local)*dx;
		CDF->push_back(sum);
		CDF_x->push_back(x_local);
	}
	for(int i = 0; i < CDF_array_size; i++)
            (*CDF)[i] = (*CDF)[i]/sum;
}

double GaussianDistribution::pdf(double x)
{
      return exp(-(x-mu)*(x-mu)/(2*sigma*sigma))/sqrt(2*M_PI*sigma*sigma);
}

double GaussianDistribution::rand()
{
	return sampleUsingInvCDF(0.0, 1.0);
}

double GaussianDistribution::invCDF(double y)
{
   long idx;
   double x;
   idx = binarySearch(CDF, y, true);
   x = (*CDF_x)[idx] + ((*CDF_x)[idx+1] - (*CDF_x)[idx])/((*CDF)[idx+1] - (*CDF)[idx])*(y - (*CDF)[idx]);
   return(x);
}

