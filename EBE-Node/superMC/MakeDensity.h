#ifndef MAKEDENSITY_h
#define MAKEDENSITY_h
#include <fstream>
#include "OverLap.h"
#include "KLNModel.h"
#include "OverLap.h"
#include "MCnucl.h"
#include "UnintegPartonDist.h"
#include "EOS.h"

#include "ParameterReader.h"

class MakeDensity
{
protected:
  OverLap* proj;
  OverLap* targ;
  UnintegPartonDist* wf;
  KLNModel* kln;
  MCnucl* mc;
  GlueDensity *rho;
  Large_x *val;

  double ecm,Anucl1,Anucl2;
  double bmin, bmax, Npmin, Npmax;
  double siginNN;
  int Maxx, Maxy;
  double Xmin, Ymin, Xmax, Ymax;
  double dx,dy;
  int    binRapidity;
  double rapMin, rapMax;
  double Npart;
  double Nbin;
  double Alpha;
  double finalFactor;
  EOS eos;

  ParameterReader *paraRdr;

 public:
  MakeDensity(ParameterReader*);
  ~MakeDensity();

  void averageDensity(const int iy, const double wei);
  void dumpDensity(std::ofstream& o, const int iy, const double wei);
  void dumpDensityBlock(char filename[], double *** data, const int iy);
  void dumpDensity4Col(char filename[], double *** data, const int iy);
  void generate(int nevent);
  void generate_profile_ebe(int nevent);
  void generate_profile_average(int nevent);
  void generateEccTable(int nevent);
  void dumpEccentricities(char*, double***, const int, int, int, double, double, double);
  void setSd(double*** data, int iy);
  void setEd(double*** data, int iy);
};

#endif
