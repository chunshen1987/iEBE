#ifndef GaussianDistribution_h
#define GaussianDistribution_h

#include "RandomVariable.h"
#include <vector>
using namespace std;


class GaussianDistribution: public RandomVariable
{

private:
	double mu, sigma;
      int CDF_array_size;
	vector<double>* CDF;
	vector<double>* CDF_x;

public:
	GaussianDistribution(double mu,double sigma);
	~GaussianDistribution();

	void initCDF();

	double pdf(double);
	double rand();
      double invCDF(double y);
};

#endif
