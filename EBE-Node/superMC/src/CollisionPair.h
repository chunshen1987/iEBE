#ifndef CollisionPair_h
#define CollisionPair_h

class CollisionPair
{
 protected:
  double x,y;
  double fluctfactor;
 public:
  double additional_weight; // store additional weight used in Uli-Glb model
  CollisionPair(double x0,double y0) {
    x = x0; y = y0; additional_weight = 0.;
  }
  ~CollisionPair() {};
  double getX() {return x;}
  double getY() {return y;}
  void setX(double a) {x=a;}
  void setY(double a) {y=a;}

  void setfluctfactor(double fluct) {fluctfactor = fluct;}
  double getfluctfactor() {return fluctfactor;}
};
#endif
