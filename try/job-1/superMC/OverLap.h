#ifndef OVERLAP_h
#define OVERLAP_h

#include <iostream>
#include <cstdlib>
#include <cmath>

class OverLap
{
protected:
    int deformed;
    double rad,rmaxCut,rwMax;
    double dr;
    double density0;
    double A;  // mass number of nucleus as double
    int atomic;// same as int
    double  sig;  // inelastic NN cross section [fm^2].

    // working area.
    double zini;
    double zfin;
    double z[38],zw[38];

    double beta2,beta4; //deformation parameters 05032010 by TH
    double ctr, phir; // ctr = cos(theta)

public:
    // a=atomic number, b=impact parameter [fm], sigin: [mb]
    OverLap(int a, double signn, int deformed=0);
    virtual ~OverLap();
    int    getAtomic() {return atomic;}

    // main function: returns a random coordinate of nucleon
    void getRandomWS(double& x, double& y, double& z);

    static void Gauss38(double xini,double xfin,double* xn,double* wn);

    //Deformation
    void getDeformRandomWS(double& x, double& y, double& z);
    void setRotation(double costheta, double phi) {ctr=costheta; phir=phi;}
    double SphericalHarmonics(int l, double theta);

};

#endif
