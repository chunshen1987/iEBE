#ifndef GaussianNucleonsCal_h

#define GaussianNucleonsCal_h

#include "ParameterReader.h"

class GaussianNucleonsCal
{
  public:
    ParameterReader *paraRdr;
    double lambda; // lambda=sigma_gg/(4*pi*width^2)
    double sigma_gg, width; // sigma_gg and gaussian width

    double entropy_gaussian_width; // (fm), used to determine the shape of the gaussian for entropy dumped into the medium from MC-Glb model; the height is a free paramter controlled by final multiplicity

    GaussianNucleonsCal(ParameterReader*);

    bool testCollision(double b);

    static double get2DHeightFromWidth(double w);

    double getSigEff(double siginNN, double width);

};
#endif

/*----------------------------------------------------------------------
 Required parameters used for Gaussian-shaped nucleons are calcualted in
 the constructor; parameters are quoted directly later (public).
----------------------------------------------------------------------*/
