#include <cmath>
#include "GlueDensity.h"
#include "MCnucl.h"

using namespace std;

GlueDensity::GlueDensity(double xmax, double ymax, double ptmin, double ptmax, double dx0,double dy0,double dpt0, 
    int ny, double rapmin, double rapmax)
{
    Xmax=xmax;
    Ymax=ymax;
    Xmin=-xmax;
    Ymin=-ymax;
    dx=dx0;
    dy=dy0;
    Maxx=(int)((Xmax-Xmin)/dx+0.1)+1;
    Maxy=(int)((Ymax-Ymin)/dy+0.1)+1;
    nRap=ny;
    rapMin=rapmin;
    rapMax=rapmax;

    PTmax=ptmax;
    PTmin=ptmin;
    dpt=dpt0;
    MaxPT=(int)((PTmax-PTmin)/dpt+0.1)+1;

    dNdy = new double [nRap];
    Xcm = new double [nRap];
    Ycm = new double [nRap];
    AngleG = new double [nRap];
    Xcm2 = new double [nRap];
    Ycm2 = new double [nRap];
    XYcm = new double [nRap];

    density  = new double** [nRap];
    for(int iy=0;iy<nRap;iy++) {
	density[iy] =  new double* [Maxx];
	for(int i=0;i<Maxx;i++) density[iy][i] = new double[Maxy];
	dNdy[iy]=0.0; Xcm[iy]=0.0; Ycm[iy]=0.0; AngleG[iy]=0.0;
	Xcm2[iy]=0.0; Ycm2[iy]=0.0; XYcm[iy]=0.0;
    }

    //density with pt dependence
    densitypt  = new double*** [nRap];
    for(int iy=0;iy<nRap;iy++) {
    densitypt[iy] =  new double** [Maxx];
    for(int i=0;i<Maxx;i++) {
        densitypt[iy][i] = new double* [Maxy];
        for(int j=0;j<Maxy;j++) {
            densitypt[iy][i][j]=new double[MaxPT];
                      for(int ipt=0;ipt<MaxPT; ipt++)
                            densitypt[iy][i][j][ipt]=0;}
        }
    }    
}


GlueDensity::~GlueDensity()
{
    for(int iy=0;iy<nRap;iy++) {
	for(int j=0;j<Maxy;j++) delete [] density[iy][j];
	delete [] density[iy];
    }
    delete [] density;

//clean densitypt
    for(int iy=0;iy<nRap;iy++) {
    for(int i=0;i<Maxx;i++) {
        for(int j=0;j<Maxy;j++) delete [] densitypt[iy][i][j];
        delete [] densitypt[iy][i];
        }
    delete [] densitypt[iy];
    }
    delete [] densitypt;     
}


void GlueDensity::getCMAngle(const int iy, int n)
{
    Xcm[iy]=0.0, Ycm[iy]=0.0, Xcm2[iy]=0.0, Ycm2[iy]=0.0,XYcm[iy]=0.0;
    double weight=0.0;
    // for un-recenterd profile
    double x, y, wei, th;
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++) {
        x = Xmin + i*dx;
        y = Ymin + j*dy;
        wei = density[iy][i][j]*dx*dy;
        weight += wei;
        Xcm[iy]  += x*wei;
        Ycm[iy]  += y*wei;
        Xcm2[iy] += x*x*wei;
        Ycm2[iy] += y*y*wei;
        XYcm[iy] += x*y*wei;
    }

    Xcm[iy] /= weight;
    Ycm[iy] /= weight;
    Xcm2[iy] /= weight;
    Ycm2[iy] /= weight;
    XYcm[iy] /= weight;
    dNdy[iy] = weight;
    sigX = Xcm2[iy] - Xcm[iy]*Xcm[iy];
    sigY = Ycm2[iy] - Ycm[iy]*Ycm[iy];
    sigXY = XYcm[iy] - Xcm[iy]*Ycm[iy];
    AngleG[iy]=atan2(2*sigXY,sigY-sigX)/2;
    // cout << "old: " << AngleG[iy] << endl;

    // for re-centerd profile
    double Num_real=0.0, Num_imag=0.0; // numberator of eccentricity
    for(int i=0;i<Maxx;i++)
    for(int j=0;j<Maxy;j++) {
        x = Xmin + i*dx - Xcm[iy];
        y = Ymin + j*dy - Ycm[iy];
        th = atan2(y,x);
        wei = density[iy][i][j];
        double rwei = pow(sqrt(x*x+y*y), n);
        Num_real += rwei*cos(n*th)*wei; //r^n weighted eccentricity
        Num_imag += rwei*sin(n*th)*wei;
    }
    if(weight <1e-15) {
        cout << "GlueDensity::getCMAngle weight =0 ? " << iy
            << " weight= " << weight
            << endl;
            exit(0);
    }
    int rand_orientation = rand() % n; //random integer from 0 to n-1
    AngleG[iy] = -atan2(-Num_imag, -Num_real)/n + 2*M_PI*rand_orientation/n; //AngleG takes the range from -pi to pi
    // cout << "imag=" << Num_imag << "," << "real=" << Num_real << endl;
    // cout << "new: " << AngleG[iy] << "," << "order=" << n << endl;
}


void GlueDensity::rotateParticle(vector<Participant*> participant,
    vector<CollisionPair*> binaryCollision, const int iy)
{
    double xcm=Xcm[iy];
    double ycm=Ycm[iy];
    double ang0 = AngleG[iy];

    int npart = participant.size();
    for(int i=0;i<npart;i++) {
      double x = participant[i]->getX()-xcm;
      double y = participant[i]->getY()-ycm;
      double ang = MCnucl::Angle(x,y);
      double r=sqrt(x*x+y*y);
      double x0 = r*cos(ang+ang0);
      double y0 = r*sin(ang+ang0);
      participant[i]->setX(x0);
      participant[i]->setY(y0);
    }

    int ncoll=binaryCollision.size();
    for(int icoll=0;icoll<ncoll;icoll++) {
      double x = binaryCollision[icoll]->getX()-xcm;
      double y = binaryCollision[icoll]->getY()-ycm;
      double ang = MCnucl::Angle(x,y);
      double r=sqrt(x*x+y*y);
      double x0 = r*cos(ang+ang0);
      double y0 = r*sin(ang+ang0);
      binaryCollision[icoll]->setX(x0);
      binaryCollision[icoll]->setY(y0);
    }

}
