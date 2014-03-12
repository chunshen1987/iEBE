#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <cstdarg>

#include "mistools.h"
using namespace std;

double Linear1dInterp(double x, const double dx, double y0, double y1)
{

	/*Linear interpolation on 1d function f(x)
	f(x)=y0+(x-x0)*(y1-y0)/(x1-x0), here x1-x0=1.
	*/
	if(dx!=1.) 
    {
    	cout << "Linear1dInterp: spacing is wrong!" << endl;
    	exit(-1);
    }
	return y0+(x-(int) x)*(y1-y0);


}

double Bilinear2dInterp(double xi, double yi, const double dx, const double dy, double y00, double y01, double y11, double y10)
{

    /* It is used for interpolate 2d function f(x,y) using form
    f(x,y) = (1-t)*(1-u)*y00+t*(1-u)*y01+t*u*y11+(1-t)*u*y10
    Where:
    f(0,0)=v00, f(0,1)=y01, f(1,1)=y11, f(1,0)=y10
    Then find the value in the grid forming by these four points.
  */
	if(dx!=1. || dy!=1.) 
    {
    	cout << "Linear1dInterp: spacing is wrong!" << endl;
    	exit(-1);
    }
	long x1a, x2a;
	double t, u;

	x1a=long (xi);
	x2a=long (yi);

	t=(xi-x1a);
	u=(yi-x2a);

	return (1-t)*(1-u)*y00+t*(1-u)*y01+t*u*y11+(1-t)*u*y10;
}


double stepfunc(double x)
{
	return (x>=0)*1.0+(x<0)*0.0;
};
