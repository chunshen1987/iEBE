/////////////////////////////////////////////////////////////////////////
//                      hydrodynamics analysis
//                          photon emission 
//
//              author: Chun Shen <shen@mps.ohio-state.edu>
//
//  This program calculates the photon emission from the relativistic
//  heavy ion collision. It reads in viscous hydrodynamics results in 
//  OSCAR format.
//  
//
//  To do in the future:
/////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>

#include "PhotonEmission.h"
#include "Hydroinfo_h5.h"
#include "Stopwatch.h"
#include "Arsenal.h"
#include "ParameterReader.h"
#include "gauss_quadrature.h"

using namespace std;

int main(int argc, char** argv)
{
  Stopwatch sw;

  sw.tic();
  ParameterReader* paraRdr = new ParameterReader();
  paraRdr->readFromFile("parameters.dat");
  paraRdr->readFromArguments(argc, argv);

  int bufferSize = paraRdr->getVal("HydroinfoBuffersize");
  int hydroInfoVisflag = paraRdr->getVal("HydroinfoVisflag");
  HydroinfoH5* hydroinfo_ptr = new HydroinfoH5();
  hydroinfo_ptr->readHydroinfoH5("results/JetData.h5", bufferSize, hydroInfoVisflag); //hydro data file pointer
  int neta = paraRdr->getVal("neta");
  double eta_i = paraRdr->getVal("eta_i");
  double eta_f = paraRdr->getVal("eta_f");
  double* eta_ptr = new double [neta];
  double* etaweight_ptr = new double [neta];
  gauss_quadrature(neta, 1, 0.0, 0.0, eta_i, eta_f, eta_ptr, etaweight_ptr);

  PhotonEmission thermalPhotons(paraRdr);

  thermalPhotons.calPhotonemission(hydroinfo_ptr, eta_ptr, etaweight_ptr);

  thermalPhotons.calPhoton_SpvnpT_individualchannel();
  thermalPhotons.calPhoton_total_SpMatrix();
  thermalPhotons.calPhoton_total_Spvn();

  thermalPhotons.outputPhotonSpvn();

  sw.toc();
  cout << "totally takes : " << sw.takeTime() << " seconds." << endl;

  return(0);
}

