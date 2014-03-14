/* For given gluon density distribution, dN/dy_d^2r_d^2Pt, calculate energy momentum tensor
T^{mu nu}. Then use Landau matching conditions to get energy density, velocity profile, 
pressure at z=0 plane. To accomplish all, an input of EOS is necessary. 
*/

#ifndef LdMatching_h
#define LdMatching_h

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdlib.h>
#include "Freestreaming.h"
#include "CellData.h"
#include "EOS.h"
#include "ParameterReader.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

using namespace std;

class LdMatching
{
private:	
    FreeStrm *Streaming;
    CellData *DataTable;
    ParameterReader *lm_params;
    double*** dNd2rdyTable;
	double Xmax,Ymax,Xmin,Ymin,dx,dy;
    double Xcm, Ycm;  //center of the energy density 
    double edMax;  //maximum energy density in the lab frame
	int    nRap;
	double rapMin, rapMax;
	int    Maxx, Maxy;
	double Tau0, Taui, Tauf, Dtau, delta_tau;
	int EOS_type;
	EOS eos;
    bool outputData;
    string Result_Dir, Dst_Folder;
    double event_phi2, event_phi3; //store the spatial angle phi2 and phi3 

    void findCM_ed(const int iRap = 0);  //find the center of the profile
    //A general function to shift energy density table to a new center (x0, y0) and/or rotate phi angle clockwisely
    double getShiftedProfile(double ***data, int i, int j, double x0, double y0, 
        double phi = 0, bool limit=false, const int iRap=0); 
    void regulateDiluteRegion(const int iRap=0); //regulate u^\mu, set it to 1 in dilute region


public:
	LdMatching(ParameterReader *params_in, string result_dir);
	~LdMatching();
	void echo();
	void Diagnostic(int iRap, int i, int j);
    void MultiMatching(string filename);

	void ReadTable(string filename);
    void CalTmunu(const int iRap); 
    void Matching_eig(const int nrap=1);
    //useful routines for matrix manipulations
    void Lower_matrix_single(gsl_matrix *dest, gsl_matrix *src);
    void Lower_matrix_double(gsl_matrix *dest, gsl_matrix *src);  //give g_ma * S^ab * g_bn
    double Contract_matrix(gsl_matrix *upper); //give the contraction T^mn * T_mn

    double GetIdealPressure(double edens);

    void CalBulkVis(const int nrap=1);
    void CalShearVis(const int nrap=1);
    void CalPresTable(const int nrap=1);
    void CalVis2Ideal_Ratio(const int iRap=0);
    void GenerateSdTable(const int nrap = 1);
    double getEpx(int n, const int iRap = 0);         //calculate eccentricity in the free-streaming stage

    // void OutputTable(const char *filename, const int iRap);
    void OutputTable_ux(const char *filename, const int iRap=0);
    void OutputTable_uy(const char *filename, const int iRap=0);
    void OutputTable_uz(const char *filename, const int iRap=1);
    void Output4colTable_ed(const char *filename, const int iRap=0);

    //for debugging
    void OutputTmnTable(const char *filename,const int iRap, const int mu, const int nu);
    void Output4colTable_visratio(const char *filename, const int iRap);

    void  OutputTable_ed(const char *filename, const int iRap);
    void  OutputTable_Sd(const char *filename, const int iRap=0);
    void  OutputTable_pressure(const char *filename, const int iRap=0);
    void  OutputTable_BulkPi(const char *filename, const int iRap=0);

    // void  CalPiSquare(const int iRap=0);
    // void  Test_piSq_part(const int iRap=0);
    void  Output_picontract_comp(const char *filename, const int iRap=0);
    void  OutputTables_Pimn(const int iRap=0);

    double getEpx_test(const int irap=1);

};

#endif
