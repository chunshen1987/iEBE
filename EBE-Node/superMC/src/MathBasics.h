// Define basic mathematics objects and implement basic operations
// Version 0.1
// 09-12-2011, Zhi Qiu

#ifndef MathBasicsHeader
#define MathBasicsHeader

#include <vector>
#include <string>

using namespace std;

class Point2D
{
  public:
    double x, y;
    Point2D(double x0, double y0);
    void rotate(double theta); // rotate the point counter-clockwisely by theta angle (in radian)
    void shift(double shift_x, double shift_y); // shift the point by shift_x and shift_y
    void printMe(); // print out (x,y)
};

class Point3D
{
  public:
    double x, y, z;
    Point3D(double x0, double y0, double z0);
    void rotate(double costheta, double phi); // rotate the point by angle theta and phi
    void shift(double shift_x, double shift_y, double shift_z); // shift the point by shift_x, shift_y, shift_z
    void printMe(); // print out (x,y,z)
};

#endif



/*----------------------------------------------------------------------
 Change logs:



-----------------------------------------------------------------------*/
