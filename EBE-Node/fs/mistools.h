#ifndef mistools_h
#define mistools_h

#include <string>
#include <vector>


double Bilinear2dInterp(double x, double y, const double dx, const double dy,
    double y00, double y01, double y11, double y10);

double Linear1dInterp(double ipts, const double dpt, double y00, double y01);
double stepfunc(double x);


#endif
