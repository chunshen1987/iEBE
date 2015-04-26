#ifndef THERMALPHOTON_H
#define THERMALPHOTON_H

#include <string>
#include<fstream>

#include "Arsenal.h"
#include "Table2D.h"
#include "ParameterReader.h"

using namespace std;


class ThermalPhoton
{
   private:
      ParameterReader* paraRdr;

      int np, nphi, nrapidity;
      int norder;
      int neta;
      string rate_path;

      //photon emission rate
      Table2D* Photonemission_eqrateTable_ptr;
      Table2D* Photonemission_viscous_rateTable_ptr;
      Table2D* Photonemission_bulkvis_rateTable_ptr;
      
      double** Emission_eqrateTb_ptr;
      double** Emission_viscous_rateTb_ptr;
      double** Emission_bulkvis_rateTb_ptr;
      double* EmissionrateTb_Yidxptr;
      double EmissionrateTb_Xmin;
      double EmissionrateTb_Ymin;
      int EmissionrateTb_sizeX;
      int EmissionrateTb_sizeY;
      double EmissionrateTb_dX;
      double EmissionrateTb_dY;


      //photon spectra parameters
      string emissionProcess_name;
      double *p, *p_weight;
      double *phi, *phi_weight;
      double *y;
      double *theta;

      double ***dNd2pTdphidy_eq, ***dNd2pTdphidy_vis, ***dNd2pTdphidy_tot;
      double ***dNd2pTdphidy_bulkvis;
      double *dNd2pT_eq, **vnpT_cos_eq, **vnpT_sin_eq;
      double *dNd2pT_vis, **vnpT_cos_vis, **vnpT_sin_vis;
      double *dNd2pT_bulkvis, **vnpT_cos_bulkvis, **vnpT_sin_bulkvis;
      double *dNd2pT_tot, **vnpT_cos_tot, **vnpT_sin_tot;

      double dNdy_eq, dNdy_vis, dNdy_tot, dNdy_bulkvis;
      double *vn_cos_eq, *vn_sin_eq;
      double *vn_cos_vis, *vn_sin_vis;
      double *vn_cos_bulkvis, *vn_sin_bulkvis;
      double *vn_cos_tot, *vn_sin_tot;
      
      //matrix for cuts on temperature and proper time
      int nTcut, nTaucut;
      double Tcut_high, Tcut_low;
      double Taucut_high, Taucut_low;
      double *****dNd2pTdphidydTdtau_eq, *****dNd2pTdphidydTdtau_tot;
      double *****dNd2pTdphidydTdtau_vis, *****dNd2pTdphidydTdtau_bulkvis;
      double **dNdydTdtau_eq, **dNdydTdtau_tot;
      double **dNdydTdtau_vis, **dNdydTdtau_bulkvis;
      double ***vndTdtau_cos_eq, ***vndTdtau_sin_eq;
      double ***vndTdtau_cos_vis, ***vndTdtau_sin_vis;
      double ***vndTdtau_cos_bulkvis, ***vndTdtau_sin_bulkvis;
      double ***vndTdtau_cos_tot, ***vndTdtau_sin_tot;

   public:
      ThermalPhoton(ParameterReader* paraRdr_in);
      ~ThermalPhoton();

      void setupEmissionrate(string emissionProcess, double Xmin, double dX,  double Ymin, double dY);
      void readEmissionrate(string);

      Table2D* get_eqRatetableptr() {return(Photonemission_eqrateTable_ptr);};
      Table2D* get_visRatetableptr() {return(Photonemission_viscous_rateTable_ptr);};
      double getPhotonp(int i) {return(p[i]);};
      double getPhoton_pweight(int i) {return(p_weight[i]);};
      double getPhotonphi(int i) {return(phi[i]);};
      double getPhoton_phiweight(int i) {return(phi_weight[i]);};
      double getPhotontheta(int i) {return(theta[i]);};
      double getPhotonrapidity(int i) {return(y[i]);};
      double getPhotonSpMatrix_eq(int i, int j, int k) {return(dNd2pTdphidy_eq[i][j][k]);};
      double getPhotonSpMatrix_tot(int i, int j, int k) {return(dNd2pTdphidy_tot[i][j][k]);};

      void getPhotonemissionRate(double* Eq, double* pi_zz, double* bulkPi, int Eq_length, double T, double* eqrate_ptr, double* visrate_ptr, double* bulkvis_ptr);

      void calThermalPhotonemission(double* Eq, double* pi_zz, double* bulkPi, int Tb_length, double T, double* volume, double fraction);
      void calThermalPhotonemissiondTdtau(double* Eq, double* pi_zz, double* bulkPi, int Tb_length, double T, double tau, double* volume, double fraction);
      void calPhoton_SpvnpT();
      void calPhoton_SpvnpT_dTdtau();
      void outputPhoton_SpvnpT(string path);
      void outputPhoton_SpvnpTdTdtau(string path);
      void interpolation2D_bilinear(double varX, double* varY, int Y_length, double** Table2D_ptr, double* results);

      void update_rates_with_polyakov_suppression();
      double get_polyakov_suppression_factor(double T_in_GeV);

};
#endif
