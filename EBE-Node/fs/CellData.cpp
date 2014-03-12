#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include "CellData.h"

using namespace std;

CellData::CellData(double xmax, double ymax, double dx0,double dy0,	    
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

    ptsData  = new cell0ptr** [nRap];

    for(int iy=0;iy<nRap;iy++) 
    {
    	ptsData[iy]= new cell0ptr* [Maxx];

    	for(int i=0;i<Maxx;i++)
    	{
    		ptsData[iy][i] = new cell0ptr [Maxy];

    		for(int j=0;j<Maxy;j++)
			{
				ptsData[iy][i][j]= (cell0*) malloc(sizeof(cell0));
				ptsData[iy][i][j]-> Tmn  = gsl_matrix_alloc (4, 4);
				ptsData[iy][i][j]-> Pi_mn= gsl_matrix_alloc (4, 4);
				ptsData[iy][i][j]-> Um   = gsl_vector_alloc (4);

				gsl_matrix_set_zero (ptsData[iy][i][j]-> Tmn);
				gsl_matrix_set_zero (ptsData[iy][i][j]-> Pi_mn);
				gsl_vector_set_zero (ptsData[iy][i][j]-> Um);
				ptsData[iy][i][j]-> Bulk_Pi =0.0;
				ptsData[iy][i][j]-> Ed      =0.0;
                ptsData[iy][i][j]-> Sd      =0.0;
				ptsData[iy][i][j]-> Pres    =0.0;
                ptsData[iy][i][j]-> Vis_ratio=0.0;
                ptsData[iy][i][j]-> PiSq    =0.0;
                ptsData[iy][i][j]-> upVal   =0.0;
				// cout<<gsl_matrix_get(ptsData[iy][i][j]-> Tmn, 1,1)<<endl;
			}  
    	}
    }
}


CellData::~CellData()
{
	for(int iy=0;iy<nRap;iy++)  
    {
        for(int i=0;i<Maxx;i++)
        { 
        	for(int j=0;j<Maxy;j++)
        	{
        		gsl_matrix_free(ptsData[iy][i][j]-> Tmn);
        		gsl_matrix_free(ptsData[iy][i][j]-> Pi_mn);
        		gsl_vector_free(ptsData[iy][i][j]-> Um);
                free((void *)(ptsData[iy][i][j]));  //call free(*void) to destroy if the 
                                                    //memory is allocated by malloc()
        	}
        delete [] ptsData[iy][i];  //if use new to create, use delete to destroy
        }
        delete [] ptsData[iy];
    }
    delete [] ptsData;
}

void CellData::Tmn_matrix_get(gsl_matrix *dest, int iy, int i, int j)
{
    gsl_matrix_memcpy(dest, ptsData[iy][i][j]-> Tmn); 
    // cout<<"Copy from the following matrix"<<endl;
    // gsl_matrix_fprintf(stdout, ptsData[iy][i][j]-> Tmn, "%.18f");
    // cout<<"The result of copy: "<<endl;
    // gsl_matrix_fprintf(stdout, dest, "%.18f");   
}

void CellData::Um_vector_get(gsl_vector *dest, int iy, int i, int j)
{
    gsl_vector_memcpy(dest, ptsData[iy][i][j]-> Um); 
    // cout<<"Copy from the following vector"<<endl;
    // gsl_vector_fprintf(stdout, ptsData[iy][i][j]-> Um, "%.18f");
    // cout<<"The result of copy: "<<endl;
    // gsl_vector_fprintf(stdout, dest, "%.18f"); 
}

void CellData::Pimn_matrix_get(gsl_matrix *dest, int iy, int i, int j)
{
    gsl_matrix_memcpy(dest, ptsData[iy][i][j]-> Pi_mn); 
    // cout<<"Copy from the following matrix"<<endl;
    // gsl_matrix_fprintf(stdout, ptsData[iy][i][j]-> Pi_mn, "%.18f");
    // cout<<"The result of copy: "<<endl;
    // gsl_matrix_fprintf(stdout, dest, "%.18f");  
}

void  CellData::demonstration(void)
{
    cout<<"Grid info:"<<endl
        <<"X points: "<<Maxx<<"  "<<"y points: "<<Maxy<<endl;


    for(int iy=0;iy<nRap;iy++)
        for(int i=0;i<Maxx;i++) 
        {   cout<<"Row: "<<i+1<<"  ";

            for(int j=0;j<Maxy;j++)
                cout<<setw(5)
                    <<GetTmn(iy, i, j, 0, 0);
            cout<<endl;
        }
    
    for(int iy=0;iy<nRap;iy++)
        for(int i=0;i<Maxx;i++)
            for(int j=0;j<Maxy;j++)
                SetTmn(iy, i, j, 0, 0, i+j);

    for(int iy=0;iy<nRap;iy++)
        for(int i=0;i<Maxx;i++)
            {for(int j=0;j<Maxy;j++)
                cout<<setw(5)<<GetTmn(iy, i, j, 0, 0);
            cout<<endl;}



}
