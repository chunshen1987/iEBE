#include <iostream>
#include <cmath>
#include "MathBasics.h"

using namespace std;

//----------------------------------------------------------------------

Point2D::Point2D(double x0, double y0)
{
  x = x0; y = y0;
}

void Point2D::rotate(double theta)
// Rotate (x,y) by theta counter-clockwisely
{
  double x0=x, y0=y;
  x = cos(theta)*x0 - sin(theta)*y0;
  y = sin(theta)*x0 + cos(theta)*y0;
}

void Point2D::shift(double shift_x, double shift_y)
// Shift (x,y) by (shift_x,shift_y)
{
  x += shift_x; y += shift_y;
}

void Point2D::printMe()
// Print out (x,y)
{
  cout << "(" << x << "," << y << ")";
}

//----------------------------------------------------------------------

Point3D::Point3D(double x0, double y0, double z0)
{
  x = x0; y = y0; z = z0;
}

void Point3D::rotate(double costheta, double phi)
// Rotate (x,y,z) by theta and phi, the first variable is cos(theta)
{
  double x0=x, y0=y, z0=z;
  double cth = costheta, cphi = cos(phi);
  double sth = sqrt(1.-cth*cth), sphi = sqrt(1.-cphi*cphi);
  x = cth*cphi*x0 - sphi*y0 + sth*cphi*z0;
  y = cth*sphi*x0 + cphi*y0 + sth*sphi*z0;
  z = -sth    *x0 + 0.  *y0 + cth     *z0;
}

void Point3D::shift(double shift_x, double shift_y, double shift_z)
// Shift (x,y,z) by (shift_x, shift_y, shift_z)
{
  x += shift_x; y += shift_y; z += shift_z;
}

void Point3D::printMe()
// Print out (x,y,z)
{
  cout << "(" << x << "," << y << "," << z << ")";
}
