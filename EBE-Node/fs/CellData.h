#ifndef CellData_h
#define CellData_h

#include <iostream>
#include <string>
#include <cstring>
#include "stdlib.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

using namespace std;

typedef struct cell0
{
	gsl_matrix *Tmn;
	gsl_matrix *Pi_mn;
	gsl_vector *Um;
	double Bulk_Pi;
	double Pres;
	double Ed;
	double Vis_ratio;
	double PiSq;
	double upVal;
	double Sd;
}cell0, *cell0ptr;   //data stored in one fluid cell


class CellData
{
protected:	
	cell0ptr ***ptsData;   //pointers matrix to structures
	double Xmax,Ymax,Xmin,Ymin,dx,dy;
	int    nRap;
	double rapMin, rapMax;
	int    Maxx, Maxy;

public:
	CellData(double xmax, double ymax, double dx0,double dy0,
		    int ny, double rapmin, double rapmax);
	~CellData();

	void   SetTmn(int iy, int idx_x, int idx_y, int row, int col, double data) 
	      { gsl_matrix_set(ptsData[iy][idx_x][idx_y]->Tmn, row, col, data); }
    double GetTmn(int iy, int idx_x, int idx_y, int row, int col)
          {return gsl_matrix_get(ptsData[iy][idx_x][idx_y]->Tmn, row, col);}
    void Tmn_matrix_get(gsl_matrix *dest, int iy, int i, int j);

	void   SetPi_mn(int iy,int idx_x, int idx_y, int row, int col, double data) { gsl_matrix_set(ptsData[iy][idx_x][idx_y]->Pi_mn, row, col, data);} 
	//set pi_/mu/nu at a point
    double GetPi_mn(int iy,int idx_x, int idx_y, int row, int col) {return gsl_matrix_get(ptsData[iy][idx_x][idx_y]->Pi_mn, row, col);}
    void Pimn_matrix_get(gsl_matrix *dest, int iy, int i, int j);

    void   SetUm(int iy, int idx_x, int idx_y, int col, double data)   {gsl_vector_set(ptsData[iy][idx_x][idx_y]->Um, col, data);}  //set u^/mu at a point
    double GetUm(int iy, int idx_x, int idx_y, int col) {return gsl_vector_get(ptsData[iy][idx_x][idx_y]->Um, col);}
    void Um_vector_get(gsl_vector *dest, int iy, int i, int j);

    void   SetEd(int iy, int idx_x, int idx_y, double data)  { ptsData[iy][idx_x][idx_y]->Ed=data;}
    double GetEd(int iy, int idx_x, int idx_y)   {return ptsData[iy][idx_x][idx_y]->Ed;}

    void   SetSd(int iy, int idx_x, int idx_y, double data)  { ptsData[iy][idx_x][idx_y]->Sd=data;}
    double GetSd(int iy, int idx_x, int idx_y)   {return ptsData[iy][idx_x][idx_y]->Sd;}

    void   SetBulk_Pi(int iy, int idx_x, int idx_y, double data) { ptsData[iy][idx_x][idx_y]->Bulk_Pi=data;}
    double GetBulk_Pi(int iy, int idx_x, int idx_y)   {return ptsData[iy][idx_x][idx_y]->Bulk_Pi;}

    void   SetPres(int iy, int idx_x, int idx_y, double data) { ptsData[iy][idx_x][idx_y]->Pres=data;}
    double GetPres(int iy, int idx_x, int idx_y)   {return ptsData[iy][idx_x][idx_y]->Pres;}

    void   SetVisRatio(int iy, int idx_x, int idx_y, double data) { ptsData[iy][idx_x][idx_y]->Vis_ratio=data;}
    double GetVisRatio(int iy, int idx_x, int idx_y)   {return ptsData[iy][idx_x][idx_y]->Vis_ratio;}
    
// Store and fetch Pi_mn*Pi^mn
    void   SetPiSq(int iy, int idx_x, int idx_y, double data) { ptsData[iy][idx_x][idx_y]->PiSq=data;}
    double GetPiSq(int iy, int idx_x, int idx_y)   {return ptsData[iy][idx_x][idx_y]->PiSq;}
// Store and fetch Integrate[(u*p)^2 * f]
    void   SetupVal(int iy, int idx_x, int idx_y, double data) { ptsData[iy][idx_x][idx_y]->upVal=data;}
    double GetupVal(int iy, int idx_x, int idx_y)   {return ptsData[iy][idx_x][idx_y]->upVal;}


//for debugging
    void demonstration(void);
};

#endif
