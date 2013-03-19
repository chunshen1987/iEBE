#include "rcBKfunc.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
using namespace std;



// functions to read tables for rcBK-uGD, and to set up interpolations

rcBKfunc::rcBKfunc(int uGDmodel) {

  // table: y kt phi1 phi2 phi3

  uGDflag = uGDmodel;  // save uGD flag

  string filename[100];
  string tmpfilename[]={ "ft_rcbk_mv_qs02_02_ad.dat",
			 "ft_rcbk_mv_qs02_03_ad.dat",
			 "ft_rcbk_mv_qs02_04_ad.dat",
			 "ft_rcbk_mv_qs02_05_ad.dat",
			 "ft_rcbk_mv_qs02_06_ad.dat",
			 "ft_rcbk_mv_qs02_07_ad.dat",
			 "ft_rcbk_mv_qs02_08_ad.dat",
			 "ft_rcbk_mv_qs02_09_ad.dat",
			 "ft_rcbk_mv_qs02_1_ad.dat",
			 "ft_rcbk_mv_qs02_11_ad.dat",
			 "ft_rcbk_mv_qs02_12_ad.dat",
			 "ft_rcbk_mv_qs02_13_ad.dat",
			 "ft_rcbk_mv_qs02_14_ad.dat",
			 "ft_rcbk_mv_qs02_15_ad.dat",
			 "ft_rcbk_mv_qs02_16_ad.dat",
			 "ft_rcbk_mv_qs02_17_ad.dat",
			 "ft_rcbk_mv_qs02_18_ad.dat",
			 "ft_rcbk_mv_qs02_19_ad.dat",
			 "ft_rcbk_mv_qs02_2_ad.dat",
			 "ft_rcbk_mv_qs02_21_ad.dat",
			 "ft_rcbk_mv_qs02_22_ad.dat",
			 "ft_rcbk_mv_qs02_23_ad.dat",
			 "ft_rcbk_mv_qs02_24_ad.dat",
			 "ft_rcbk_mv_qs02_25_ad.dat",
			 "ft_rcbk_mv_qs02_26_ad.dat",
			 "ft_rcbk_mv_qs02_27_ad.dat",
			 "ft_rcbk_mv_qs02_28_ad.dat",
			 "ft_rcbk_mv_qs02_29_ad.dat",
			 "ft_rcbk_mv_qs02_3_ad.dat",
			 "ft_rcbk_mv_qs02_31_ad.dat",
			 "ft_rcbk_mv_qs02_32_ad.dat",
			 "ft_rcbk_mv_qs02_33_ad.dat",
			 "ft_rcbk_mv_qs02_34_ad.dat",
			 "ft_rcbk_mv_qs02_35_ad.dat",
			 "ft_rcbk_mv_qs02_36_ad.dat",
			 "ft_rcbk_mv_qs02_37_ad.dat",
			 "ft_rcbk_mv_qs02_38_ad.dat",
			 "ft_rcbk_mv_qs02_39_ad.dat",
			 "ft_rcbk_mv_qs02_4_ad.dat",
			 "ft_rcbk_mv_qs02_41_ad.dat",
			 "ft_rcbk_mv_qs02_42_ad.dat",
			 "ft_rcbk_mv_qs02_43_ad.dat",
			 "ft_rcbk_mv_qs02_44_ad.dat",
			 "ft_rcbk_mv_qs02_45_ad.dat",
			 "ft_rcbk_mv_qs02_46_ad.dat",
			 "ft_rcbk_mv_qs02_47_ad.dat",
			 "ft_rcbk_mv_qs02_48_ad.dat",
			 "ft_rcbk_mv_qs02_49_ad.dat",
			 "ft_rcbk_mv_qs02_5_ad.dat",
			 "ft_rcbk_mv_qs02_51_ad.dat",
			 "ft_rcbk_mv_qs02_52_ad.dat",
			 "ft_rcbk_mv_qs02_53_ad.dat",
			 "ft_rcbk_mv_qs02_54_ad.dat",
			 "ft_rcbk_mv_qs02_55_ad.dat",
			 "ft_rcbk_mv_qs02_56_ad.dat",
			 "ft_rcbk_mv_qs02_57_ad.dat",
			 "ft_rcbk_mv_qs02_58_ad.dat",
			 "ft_rcbk_mv_qs02_59_ad.dat",
			 "ft_rcbk_mv_qs02_6_ad.dat"
  };
  string tmpfile2[]={ "ft_rcbk_mv_qs02_0168_g1_119_1.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_2.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_3.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_4.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_5.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_6.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_7.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_8.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_9.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_10.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_11.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_12.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_13.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_14.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_15.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_16.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_17.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_18.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_19.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_20.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_21.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_22.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_23.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_24.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_25.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_26.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_27.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_28.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_29.dat",
		      "ft_rcbk_mv_qs02_0168_g1_119_30.dat",
  };


  switch (uGDmodel) {
  case rcBKalbacete : 
    maxQ0=59;  // number of Qs_0 bins = number of files
    dQ0 = 0.1; // difference of Qs_0^2 from one table to the next
    for (int i=0; i<maxQ0; i++) filename[i] = tmpfilename[i];
    UGD_maxKt = 20.0;  // table up to kt=20GeV
    break;
  case rcBKalbaceteSet2 : 
    maxQ0=30;    // number of Qs_0 bins = number of files
    dQ0 = 0.168; // difference of Qs_0^2 from one table to the next
    for (int i=0; i<maxQ0; i++) filename[i] = tmpfile2[i];
    UGD_maxKt = 20.0;  // table up to kt=20GeV
    break;
  default :
    cout << "ERROR in rcBKfunc::rcBKfunc() : unknown uGD model !\n";
    exit(0);
  }
  maxY=121;  // number of Y bins (from 0..12, dY=0.1)
  maxKt=101; // number of phi(kt) points in each Y bin


  table   = new double [maxQ0*maxY*maxKt];
  tableNF = new double [maxQ0*maxY*maxKt];
  ktPoint = new double [maxQ0*maxY*maxKt];
  QsTab   = new double [maxQ0*maxY];
  double *tdumy=table;      // save pointers to first entry of tables
  double *tdumyNF=tableNF;
  double *ktdumy=ktPoint;
  double *Qsdumy = QsTab;

  for(int iq=0; iq<maxQ0; iq++) {

    string fname = "javier/"+filename[iq];
    ifs.open(fname.c_str());
    if(!ifs) { 
      cerr << "Error unable to open file " << fname << endl;
      exit(1);
    }
    double dumy1;    // dummy variable for reading tables

    for(int iy=0; iy<maxY; iy++) {
      *QsTab = 0.0;
      for(int ikt=0; ikt<maxKt; ikt++) {
	//      Y         kt       N_F(kt)  N_A(kt)
	ifs >> dumy1 >> *ktPoint >> *tableNF >> *table;
	if (ifs.eof()) {
	  cerr << "ERROR reading phi(x,kt) tables, too few entries !\n";
	  exit(1);
	}
	// Qs(y;Q0) is defined by max of N_F(kt)*kt^2
	if (ikt)
	  if (pow(*ktPoint,2.) * (*tableNF) > 
	      pow(*(ktPoint-1),2.) * (*(tableNF-1)))   *QsTab = *ktPoint;
	ktPoint++;  table++;  tableNF++;
      }
      QsTab++;
    }
    ifs >> dumy1;
    if (!ifs.eof()) { // read beyond last entry to raise eof
      cerr << "ERROR reading phi(x,kt) tables, too many entries !\n";
      exit(1);
    }
    ifs.close();
  }
  table=tdumy;   // reset pointers to beginning of tables
  ktPoint=ktdumy;
  QsTab = Qsdumy;
  tableNF = tdumyNF;

  /* allocate accelerator and spline object for every value of Y and Q0 */
  acc = new gsl_interp_accel**[maxQ0];
  spline = new gsl_spline**[maxQ0];
  accNF = new gsl_interp_accel**[maxQ0];
  splineNF = new gsl_spline**[maxQ0];
  for(int iq=0;iq<maxQ0;iq++) {
    acc[iq] = new gsl_interp_accel*[maxY];
    spline[iq] = new gsl_spline *[maxY];
    accNF[iq] = new gsl_interp_accel*[maxY];
    splineNF[iq] = new gsl_spline *[maxY];
  }

  for(int iq=0; iq<maxQ0; iq++)
    for(int iy=0; iy<maxY; iy++) {
      acc[iq][iy] = gsl_interp_accel_alloc();
      accNF[iq][iy] = gsl_interp_accel_alloc();
      spline[iq][iy] = gsl_spline_alloc(gsl_interp_cspline, maxKt);
      splineNF[iq][iy] = gsl_spline_alloc(gsl_interp_cspline, maxKt);
      // init GSL spline for N_A(kt) and N_F(kt)
      gsl_spline_init (spline[iq][iy], &ktPoint[0 + maxKt*(iy + maxY*iq)],
		       &table[0 + maxKt*(iy + maxY*iq)], maxKt);
      gsl_spline_init (splineNF[iq][iy], &ktPoint[0 + maxKt*(iy + maxY*iq)],
		       &tableNF[0 + maxKt*(iy + maxY*iq)], maxKt);
    }
}


//#define MAIN
#ifdef MAIN
//g++  -lgsl -lgslcblas
int main()
{
    rcBKfunc *rcbk = new rcBKfunc();

    double qs2=6.;
    double x=0.0001;
    double alp=0.2;

    for (double kt2=1.0; kt2<10.; kt2+=1.)
      cout << sqrt(kt2) << " " << rcbk->getFunc(qs2, x, kt2,  alp) << endl;
    return 0;
}
#endif
