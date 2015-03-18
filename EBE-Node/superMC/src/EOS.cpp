// Ver. 1.2
// Starting from version 1.2, this EOS is only valid for s95p-PCE at large energy
// because when the values are out of table a fitted formula is used instead.

// The EOS data file is a 4-column double data file (see loadEOSFromFile).
// The energy density is assumed to be equal-spaced; if not, change all
// the interpCubicDirect fucntions to interpCubic functions and all the
// invertTableDirect functions to invertTable functions in the following.
// Note that using a equal-spaced energy density increases the efficiency
// due to the fact otherwise the interpolation method involves an inversion
// operation.

#include <fstream>
#include <iostream>
#include <string>
#include "stdlib.h"
#include <cmath>
#include "arsenal.h"
#include "EOS.h"

#define INIT_GUESS 40
#define ACCURACY 1e-4

EOS *zq_global_eos;

using namespace std;

//----------------------------------------------------------------------
EOS::EOS() {};


//----------------------------------------------------------------------
EOS::EOS(char* filename) { loadEOSFromFile(filename); };

//----------------------------------------------------------------------
EOS::EOS(char* data_filename, char* coeff_filename) {
    loadEOSFromFile(data_filename, coeff_filename);
};

//----------------------------------------------------------------------
void EOS::loadEOSFromFile(char* data_filename)
// The EOS data file (data_filename) is assumed to be a 4 column file:
// 1st column: the energy density
// 2nd column: the pressure
// 3rd column: the entropy density
// 4th column: the temprature
// The units are not relavent here; the external program is responsible
// for converting to the correct units in according to the EOS data file.
{
    fstream fs(data_filename);
    if (fs.is_open()==false)
    {
        cout << "EOS::loadEOSFromFile error: the EOS data file cannot be opened." << endl;
        exit(-1);
    }
    vector< vector<double>* >* data = readBlockData(fs);
    if ((*data).size()!=4)
    {
        cout << "EOS::loadEOSFromFile error: the EOS data file is not a valid 4-column data file." << endl;
        exit(-1);
    }
    ed_table = (*data)[0];
    p_table = (*data)[1];
    sd_table = (*data)[2];
    T_table = (*data)[3];
    table_length = ed_table->size();
    delta_ed = (*ed_table)[1]-(*ed_table)[0];
    max_ed = (*ed_table)[table_length-1];
    zq_global_eos = this;
    delete data;
};


//----------------------------------------------------------------------
void EOS::loadEOSFromFile(char* data_filename, char* coeff_filename)
// The version takes an additional file coeff_filename, containing 6 numbers
// p1,p2,s1,s2,T1,T2, and they are assumed to give the relations:
// p = p1*ed^p2; s = s1*ed^s2; T = T1*ed^T2;
// And they are used for extrapolation when EOS table is not long enough.
{
    loadEOSFromFile(data_filename); // setup the EOS table
    fstream fs(coeff_filename);
    if (fs.is_open()==false)
    {
        cout << "EOS::loadEOSFromFile error: the coeff file cannot be opened." << endl;
        exit(-1);
    }
    fs >> p1 >> p2 >> s1 >> s2 >> T1 >> T2;
    //cout << p1 <<" "<< p2 <<" "<< s1 <<" "<< s2 <<" "<< T1 <<" "<< T2 << endl; // for debug
}

//----------------------------------------------------------------------
double EOS::p(double ed0)
// Return the pressure from the energy density ed0.
{
    if (ed0<=(*ed_table)[0]) return ed0/(*ed_table)[0]*(*p_table)[0]; // linear interpolation to 0
    if (ed0>=max_ed) return p1*pow(ed0,p2); // log-fitted
    return interpCubicDirect(ed_table, p_table, ed0);
};


//----------------------------------------------------------------------
double EOS::sd(double ed0)
// Return the entropy density from the energy density ed0.
{
    if (ed0<=(*ed_table)[0]) return ed0/(*ed_table)[0]*(*sd_table)[0]; // linear interpolation to 0
    if (ed0>=max_ed) return s1*pow(ed0,s2); // log-fitted
    return interpCubicDirect(ed_table, sd_table, ed0);
};


//----------------------------------------------------------------------
double EOS::T(double ed0)
// Return the temperature from the energy density ed0.
{
    if (ed0<=(*ed_table)[0]) return ed0/(*ed_table)[0]*(*T_table)[0]; // linear interpolation to 0
    if (ed0>=max_ed) return T1*pow(ed0,T2); // log-fitted
    return interpCubicDirect(ed_table, T_table, ed0);
};


//----------------------------------------------------------------------
double zq_global_edFromP_hook(double ee) { return zq_global_eos->p(ee); }
double EOS::edFromP(double p0)
// Return the energe density from given pressure p0.
{
    return invertFunc(&zq_global_edFromP_hook, p0, (*ed_table)[0], (*ed_table)[table_length-1], delta_ed, (*ed_table)[INIT_GUESS], ACCURACY);
}


//----------------------------------------------------------------------
double zq_global_edFromSd_hook(double ee) { return zq_global_eos->sd(ee); }
double EOS::edFromSd(double sd0)
// Return the energe density from given entropy density sd0.
{
    return invertFunc(&zq_global_edFromSd_hook, sd0, (*ed_table)[0], (*ed_table)[table_length-1], delta_ed, (*ed_table)[INIT_GUESS], ACCURACY);
};


//----------------------------------------------------------------------
double zq_global_edFromT_hook(double ee) { return zq_global_eos->T(ee); }
double EOS::edFromT(double T0)
// Return the energe density from given temperature T0.
{
    return invertFunc(&zq_global_edFromT_hook, T0, (*ed_table)[0], (*ed_table)[table_length-1], delta_ed, (*ed_table)[INIT_GUESS], ACCURACY);
};

// Ver 1.1:
// -- In the xxFromEd functions, the invertFunc functions are used instead
//    of invertTable functions to increase stability of the inverting algorithm.
// -- In loadEOSFromFile, the code "delete data;" added to avoid memory leakage.
// Ver 1.2:
// -- The variable max_ed initialized in loadEOSFromFile function.
// -- For p, sd, and T functions, EOS is extended by a formula from log-fit,
//    of which the coefficients are read from a coeff file. The coeff file should
//    contain 6 number, see loadEOSFromFile(char*, char*) for details.
