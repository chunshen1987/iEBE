#ifndef Freestreaming_h
#define Freestreaming_h

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include "stdlib.h"

/*
revise history:
Apr.02, 2013
Add the function to calculate spatial eccentricity

Mar.22, 2013
1. To maintain higher precision, use cubic interpolation again, but set negative 
   number to zero by hand;
2. Add dEdy table, function generateEdTable() generateEpx() to calculate eccentricity in free-streaming
stage

Mar.18, 2013, 
In function switch to linear interpolation to avoid negative value returned by cubic method
the output file of both methods have been compaired:
relative error 2.76%, difference of total entropy densities -2.7345e-04;
*/

using namespace std;


class FreeStrm
{
protected:
	double ***shiftedTable;
	double ***unshiftedTable;
	double Xmax,Ymax,Xmin,Ymin,dx,dy;
	int    nRap;
	double rapMin, rapMax;
	int    Maxx, Maxy;
	double Taui, Tauf;
	double PTmin, dpt, PTmax, MaxPT;
	double Xcm, Ycm;    //center of the energy density profile  

public:
	FreeStrm(double xmax, double ymax, double dx0,double dy0,
			    int nrap, double rmin, double rmax, double taui, double tauf);
	~FreeStrm();

    void CopyTable(double ***source);  //copy data table from outside
	double GetDensity(int iy, int i, int j, double phip);  //debugging
	void ShiftDensity(const int iy, double phip);   //free stream density to another coordinate with a specific angle

	void OutputTable(const char* filename, const int iy);  //output free steamed and pt integrated gluon density
	void dumpBlockTable(const char *filename, double ***data, const int iy);

	void CreateDataTable(const char *filename, const int iy=0); //testing
	double GaussProfile(int irap, int i, int j, int ipt);  		//generating test profile
	double BoxProfile(const int iRap, int i, int j, int ipt);
	
	double getEpx(int n, const int iRap=0);    //calculate eccentricity in the free-streaming stage
};


#endif
