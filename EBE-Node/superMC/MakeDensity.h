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
  int cutdSdy;
  double cutdSdy_lowerBound, cutdSdy_upperBound;
  double siginNN;
  int Maxx, Maxy, MaxPT;
  double Xmin, Ymin, Xmax, Ymax;
  double dx,dy;
  int    binRapidity;
  double rapMin, rapMax;
  //double cx, phi; //andy
  double Npart;
  double Nbin;
  double Alpha;
  double finalFactor;
  double PTinte, PTmax, PTmin, dpt, MixedMode;
  EOS eos;
  int Operation;

  ParameterReader *paraRdr;

 public:
  MakeDensity(ParameterReader*);
  ~MakeDensity();

  void averageDensity(const int iy, const double wei);
  void dumpDensity(std::ofstream& o, const int iy, const double wei);
  void dumpDensityBlock(char filename[], double *** data, const int iy);
  void dumpDensity4Col(char filename[], double *** data, const int iy);
  void dumpDensity5Col(char filename[], double **** data, const int iy);
  void dumpDesityptCol(char filename[], double **** data, const int iy);
  void generate(int nevent);
  void generate_profile_ebe_Jet(int nevent);
  void generate_profile_ebe(int nevent);
  void generate_profile_average(int nevent);
  void generateEccTable(int nevent);
  void dumpEccentricities(char* base_filename, double*** density, const int iy, int from_order, int to_order, double Npart_current, double Nbin_current, double b);
  void setSd(double*** data, int iy);
  void setEd(double*** data, int iy);
  void setSd(double**** data, int iy);
  void setEd(double**** data, int iy);
  double gettotaldSdy(const int iy);
};

#endif
