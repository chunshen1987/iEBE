#ifndef Participant_h
#define Participant_h

#include <vector>
#include "Particle.h"

class Participant
{
protected:
    int nucl;
    Particle* part;
    double xsave,ysave;
    double fluctfactor;
public:
    std::vector<int> who_hit_me; // those nucleons collided with this nucleon, their indices in the binaryCollision array is stored here. This is not against the capsulate rule since the interpretation of these vector replied on external programs, and this vector is just a storage space for additional info.
    Participant(Particle* part0,int i) {
      part=part0;
      nucl=i;
      xsave = part->getX();
      ysave = part->getY();
    }
    ~Participant() {};

    Particle* getParticle() {return part;}
    int getNumberOfCollision() {return part->getNumberOfCollision();}
    double getX() {return part->getX();}
    double getY() {return part->getY();}
    double getXorg() {return xsave;}
    double getYorg() {return ysave;}
    void setX(double a) {part->setX(a);}
    void setY(double a) {part->setY(a);}
    int    isNucl() {return nucl;}
    void   setNucl(int i) {nucl=i;}
    void resetCoordinate() {part->setX(xsave);part->setY(ysave);}

    void setfluctfactor(double fluct) {fluctfactor = fluct;}
    double getfluctfactor() {return fluctfactor;}
};

#endif
