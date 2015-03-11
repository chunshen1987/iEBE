#ifndef Spectator_h
#define Spectator_h

#include <vector>
#include "Particle.h"

class Spectator
{
   private:
      double x, y, rapidity_Y;

   public:
      Spectator(double xi, double yi, double rapi)
      {
         x = xi;
         y = yi;
         rapidity_Y = rapi;
      };
      ~Spectator() {};
      
      double getX() {return(x);};
      double getY() {return(y);};
      double getRapidity_Y() {return(rapidity_Y);};
      
};

#endif
