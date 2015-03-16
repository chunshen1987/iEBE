#ifndef OVERLAP_h
#define OVERLAP_h

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include "ParameterReader.h"
#include "HulthenFunc.h"

class OverLap
{
private:
    ParameterReader* paraRdr;

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

    HulthenFunc sample_deuteron;
    vector< vector<double> > triton_pos;

    int flag_NN_correlation;
    int n_configuration;
    double*** nucleon_pos_array;

public:
    // a=atomic number, b=impact parameter [fm], sigin: [mb]
    OverLap(ParameterReader* paraRdr_in, int a, double signn, int deformed=0);
    virtual ~OverLap();
    int    getAtomic() {return atomic;}

    // main function: returns a random coordinate of nucleon
    void getRandomWS(double& x, double& y, double& z);

    static void Gauss38(double xini,double xfin,double* xn,double* wn);

    //Deformation
    void getDeformRandomWS(double& x, double& y, double& z);
    void setRotation(double costheta, double phi) {ctr=costheta; phir=phi;}
    double SphericalHarmonics(int l, double theta);
    void readin_nucleon_positions();
    void get_nucleon_position_with_NN_correlation(double **nucleon_ptr);

    // nucleon positions for light nuclei
    void GetDeuteronPosition(double& x1,double& y1,double& z1,double& x2,double& y2,double& z2);
    void readin_triton_position();
    void GetTritonPosition(double& x1,double& y1,double& z1,double &x2,double& y2,double& z2, double &x3, double &y3, double &z3);

};

#endif
