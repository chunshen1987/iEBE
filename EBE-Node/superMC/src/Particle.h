#ifndef PARTICLE_h
#define PARTICLE_h

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "GaussianDistribution.h"

class Particle
{
 protected:
  double x,y,z;
  int numberOfCollision;
  
  // added by Kevin Welsh
  bool QuarkGen;
  double ValenceQuarks[3][3];
  double flucfactors[3];

 public:
  Particle(double x0,double y0, double z0);
  ~Particle();

  double getX() {return x;}
  double getY() {return y;}
  double getZ() {return z;}
  void   setX(double a) {x=a;}
  void   setY(double a) {y=a;}
  void   setZ(double a) {z=a;}

  int    getNumberOfCollision() {return numberOfCollision;}
  void   setNumberOfCollision() {numberOfCollision++;}
  void   setNumberOfCollision(int i) {numberOfCollision=i;}

  // functions for nucleon substructure added by Kevin Welsh
  void setfluctfactorQuarks(double f1, double f2, double f3);
  void getQuarkPos(double Q[][3], GaussianDistribution *gaussDist);
  double getInternalStructDensity(double xg, double yg, double quarkWidth, GaussianDistribution *gaussDist);

};

#endif
