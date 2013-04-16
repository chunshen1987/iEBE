#ifndef GlueDensity_h
#define GlueDensity_h

#include <vector>
#include "Participant.h"
#include "CollisionPair.h"

class GlueDensity
{
protected:
    double ***density;
    double Xmax,Ymax,Xmin,Ymin,dx,dy;
    int    nRap;
    double rapMin, rapMax;
    double *dNdy;
    double *Xcm, *Ycm, *AngleG;
    double *Xcm2, *Ycm2, *XYcm;
    double sigX, sigY, sigXY;
    int    Maxx, Maxy;
public:
    GlueDensity(double xmax, double ymax,double dx0,double dy0,
		    int nrap, double rmin, double rmax);
    ~GlueDensity();

    void setDensity(int i, int ix, int iy,double a) {density[i][ix][iy]=a;}
    double getDensity(int iy, int i, int j) {return density[iy][i][j];}
    double getXcm(int i) {return Xcm[i];}
    double getYcm(int i) {return Ycm[i];}
    double getXcm2(int i) {return Xcm2[i];}
    double getYcm2(int i) {return Ycm2[i];}
    double getXYcm(int i) {return XYcm[i];}
    double getSigX() {return sigX;}
    double getSigY() {return sigY;}
    double getSigXY() {return sigXY;}

    //void setXcm(double a) {Xcm=a;}
    //void setYcm(double a) {Ycm=a;}
    //void setAngle(double a) {AngleG=a;}

    void getCMAngle(const int iy, int n=2);
    //void rotateGrid(const int ix, const int iy);
    void rotateParticle(std::vector<Participant*> participant,
			std::vector<CollisionPair*> binaryCollision,
			const int iy);

};

#endif
