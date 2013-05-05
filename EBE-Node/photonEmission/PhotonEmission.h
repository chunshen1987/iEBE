#ifndef PHOTONEMISSION_H
#define PHOTONEMISSION_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "Hydroinfo_h5.h"
#include "ThermalPhoton.h"
#include "ParameterReader.h"

using namespace std;

class PhotonEmission
{
   private:
      ParameterReader *paraRdr;
      string output_path;

      int neta;
      int np, nphi, nrapidity;
      int norder;

      double gridDx, gridDy, gridDtau;
      double gridX0, gridY0, gridTau0;
      int gridNx, gridNy;

      double T_dec, T_sw_high, T_sw_low;

      int calHGIdFlag;

      double** lambda; // Lorentz boost transverse only
      double* Eq_localrest_Tb;
      double* pi_photon_Tb;

      double ***dNd2pTdphidy_eq, *dNd2pT_eq;
      double **vnpT_cos_eq, **vnpT_sin_eq;
      double ***dNd2pTdphidy;
      double *dNd2pT;
      double **vnpT_cos, **vnpT_sin;

      //photon production processes
      ThermalPhoton* photon_QGP;
      ThermalPhoton* photon_HG;

      ThermalPhoton* photon_pirho;
      ThermalPhoton* photon_KstarK;
      ThermalPhoton* photon_piK;
      ThermalPhoton* photon_piKstar;
      ThermalPhoton* photon_pipi;
      ThermalPhoton* photon_rhoK;
      ThermalPhoton* photon_rho;
      ThermalPhoton* photon_pirho_omegat;

   public:
      PhotonEmission(ParameterReader* paraRdr_in);
      ~PhotonEmission();
      
      void set_hydroGridinfo();
      void print_hydroGridinfo();
      void InitializePhotonEmissionRateTables();
      void calPhotonemission(HydroinfoH5* hydroinfo_ptr, double* eta_ptr, double* etaweight_ptr);
      void calPhoton_total_SpMatrix();
      void calPhoton_SpvnpT_individualchannel();
      void calPhoton_total_Spvn();
      void outputPhoton_total_SpvnpT(string );
      void outputPhotonSpvn();
};

#endif
