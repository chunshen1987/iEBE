#ifndef PARTICLE_h
#define PARTICLE_h

#include <iostream>
#include <cstdlib>
#include <cmath>

class Particle
{
 protected:
  double x,y,z;
  int numberOfCollision;

 public:
  Particle(double x0,double y0, double z0) {
    x = x0; y = y0; z = z0;
    numberOfCollision=0;
  }
  ~Particle() {};

  double getX() {return x;}
  double getY() {return y;}
  double getZ() {return z;}
  void   setX(double a) {x=a;}
  void   setY(double a) {y=a;}
  void   setZ(double a) {z=a;}

  int    getNumberOfCollision() {return numberOfCollision;}
  void   setNumberOfCollision() {numberOfCollision++;}
  void   setNumberOfCollision(int i) {numberOfCollision=i;}
};

#endif
