#include <cmath>
#include <fstream>
#include <istream>
#include <iomanip>

#include "Particle.h"
#include "GaussianDistribution.h"
using namespace std;

Particle::Particle(double x0, double y0, double z0)
{
   x = x0; 
   y = y0; 
   z = z0;
   numberOfCollision=0;
   
   QuarkGen = false;
   flucfactors[0] = 1;
   flucfactors[1] = 1;
   flucfactors[2] = 1;
}

Particle::~Particle()
{
}

double Particle::getInternalStructDensity(double xg, double yg, double quarkWidth, GaussianDistribution *gaussDist)
{
   double dens = 0;

   if(!QuarkGen)
     getQuarkPos(ValenceQuarks,gaussDist);
   
   double d;
   for(int i(0);i<3;i++)
   {
   	d = pow(ValenceQuarks[i][0]+x-xg,2)+pow(ValenceQuarks[i][1]+y-yg,2); //The squared distance between the Quark and the grid
   	dens += flucfactors[i]*(1./(2*M_PI*quarkWidth*quarkWidth))*exp(-d/(2*quarkWidth*quarkWidth))/3; //Divide a total of 1 density between 3 quarks
   }
   return dens;
}

void Particle::getQuarkPos(double Q[][3], GaussianDistribution *gaussDist)
{
   double x, y;

   int i=0;
   while(i<3)
   {
     x = gaussDist->rand();
     y = gaussDist->rand();

     Q[i][0] = x;
     Q[i][1] = y;
     
     i++;
   }

   QuarkGen = true;
}

void Particle::setfluctfactorQuarks(double f1, double f2, double f3)
{
   flucfactors[0] = f1;
   flucfactors[1] = f2;
   flucfactors[2] = f3;
}
