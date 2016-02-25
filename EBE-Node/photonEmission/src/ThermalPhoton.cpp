/////////////////////////////////////////////////////////////////////////
//  To do in the future:
//      change the integration routines into gaussian integration
/////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "Table2D.h"
#include "ThermalPhoton.h"
#include "ParameterReader.h"
#include "gauss_quadrature.h"

using namespace std;

ThermalPhoton::ThermalPhoton(ParameterReader* paraRdr_in)
{
    paraRdr = paraRdr_in;

    neta = paraRdr->getVal("neta");
    np = paraRdr->getVal("np");
    nphi = paraRdr->getVal("nphi");
    nrapidity = paraRdr->getVal("nrapidity");
    norder = paraRdr->getVal("norder");
    rate_path = "ph_rates/";

    //initial variables for photon spectra 
    double p_i = paraRdr->getVal("photon_q_i"); 
    double p_f = paraRdr->getVal("photon_q_f");
    double phi_i = paraRdr->getVal("photon_phi_q_i");
    double phi_f = paraRdr->getVal("photon_phi_q_f");
    double y_i = paraRdr->getVal("photon_y_i");
    double y_f = paraRdr->getVal("photon_y_f");
    double dy = (y_f - y_i)/(nrapidity - 1 + 1e-100);
    
    p = new double [np];
    p_weight = new double [np];
    phi = new double [nphi];
    phi_weight = new double [nphi];

    gauss_quadrature(np, 1, 0.0, 0.0, p_i, p_f, p, p_weight);
    gauss_quadrature(nphi, 1, 0.0, 0.0, phi_i, phi_f, phi, phi_weight);

    y = new double [nrapidity];
    theta = new double [nrapidity];
    for(int i=0;i<nrapidity;i++) 
    {
        y[i] = y_i + i*dy;
        theta[i] = acos(tanh(y[i]));  //rapidity's corresponding polar angle
    }
    
    dNd2pT_eq = new double [np];
    dNd2pT_vis = new double [np];
    dNd2pT_vis_deltaf_restricted = new double [np];
    dNd2pT_bulkvis = new double [np];
    dNd2pT_bulkvis_deltaf_restricted = new double [np];
    dNd2pT_tot = new double [np];
    dNd2pTdphidy_eq = new double** [np];
    dNd2pTdphidy_vis = new double** [np];
    dNd2pTdphidy_vis_deltaf_restricted = new double** [np];
    dNd2pTdphidy_bulkvis = new double** [np];
    dNd2pTdphidy_bulkvis_deltaf_restricted = new double** [np];
    dNd2pTdphidy_tot = new double** [np];
    for(int i=0;i<np;i++)
    {
      dNd2pT_eq[i] = 0.0;
      dNd2pT_vis[i] = 0.0;
      dNd2pT_vis_deltaf_restricted[i] = 0.0;
      dNd2pT_bulkvis[i] = 0.0;
      dNd2pT_bulkvis_deltaf_restricted[i] = 0.0;
      dNd2pT_tot[i] = 0.0;
      dNd2pTdphidy_eq[i] = new double* [nphi];
      dNd2pTdphidy_vis[i] = new double* [nphi];
      dNd2pTdphidy_vis_deltaf_restricted[i] = new double* [nphi];
      dNd2pTdphidy_bulkvis[i] = new double* [nphi];
      dNd2pTdphidy_bulkvis_deltaf_restricted[i] = new double* [nphi];
      dNd2pTdphidy_tot[i] = new double* [nphi];
      for(int j=0;j<nphi;j++)
      {
         dNd2pTdphidy_eq[i][j] = new double [nrapidity];
         dNd2pTdphidy_vis[i][j] = new double [nrapidity];
         dNd2pTdphidy_vis_deltaf_restricted[i][j] = new double [nrapidity];
         dNd2pTdphidy_bulkvis[i][j] = new double [nrapidity];
         dNd2pTdphidy_bulkvis_deltaf_restricted[i][j] = new double [nrapidity];
         dNd2pTdphidy_tot[i][j] = new double [nrapidity];
         for(int k=0;k<nrapidity;k++)
         {
            dNd2pTdphidy_eq[i][j][k] = 0.0;
            dNd2pTdphidy_vis[i][j][k] = 0.0;
            dNd2pTdphidy_vis_deltaf_restricted[i][j][k] = 0.0;
            dNd2pTdphidy_bulkvis[i][j][k] = 0.0;
            dNd2pTdphidy_bulkvis_deltaf_restricted[i][j][k] = 0.0;
            dNd2pTdphidy_tot[i][j][k] = 0.0;
         }
      }
    }
    vnpT_cos_eq = new double* [norder];
    vnpT_sin_eq = new double* [norder];
    vnpT_cos_vis = new double* [norder];
    vnpT_sin_vis = new double* [norder];
    vnpT_cos_vis_deltaf_restricted = new double* [norder];
    vnpT_sin_vis_deltaf_restricted = new double* [norder];
    vnpT_cos_bulkvis = new double* [norder];
    vnpT_sin_bulkvis = new double* [norder];
    vnpT_cos_bulkvis_deltaf_restricted = new double* [norder];
    vnpT_sin_bulkvis_deltaf_restricted = new double* [norder];
    vnpT_cos_tot = new double* [norder];
    vnpT_sin_tot = new double* [norder];
    for(int order=0; order<norder; order++)
    {
      vnpT_cos_eq[order] = new double [np];
      vnpT_sin_eq[order] = new double [np];
      vnpT_cos_vis[order] = new double [np];
      vnpT_sin_vis[order] = new double [np];
      vnpT_cos_vis_deltaf_restricted[order] = new double [np];
      vnpT_sin_vis_deltaf_restricted[order] = new double [np];
      vnpT_cos_bulkvis[order] = new double [np];
      vnpT_sin_bulkvis[order] = new double [np];
      vnpT_cos_bulkvis_deltaf_restricted[order] = new double [np];
      vnpT_sin_bulkvis_deltaf_restricted[order] = new double [np];
      vnpT_cos_tot[order] = new double [np];
      vnpT_sin_tot[order] = new double [np];
      for(int i =0; i < np; i++)
      {
         vnpT_cos_eq[order][i] = 0.0;
         vnpT_cos_vis[order][i] = 0.0;
         vnpT_cos_vis_deltaf_restricted[order][i] = 0.0;
         vnpT_cos_bulkvis[order][i] = 0.0;
         vnpT_cos_bulkvis_deltaf_restricted[order][i] = 0.0;
         vnpT_cos_tot[order][i] = 0.0;
         vnpT_sin_eq[order][i] = 0.0;
         vnpT_sin_vis[order][i] = 0.0;
         vnpT_sin_vis_deltaf_restricted[order][i] = 0.0;
         vnpT_sin_bulkvis[order][i] = 0.0;
         vnpT_sin_bulkvis_deltaf_restricted[order][i] = 0.0;
         vnpT_sin_tot[order][i] = 0.0;
      }
    }

    dNdy_eq = 0.0;
    dNdy_vis = 0.0;
    dNdy_vis_deltaf_restricted = 0.0;
    dNdy_bulkvis = 0.0;
    dNdy_bulkvis_deltaf_restricted = 0.0;
    dNdy_tot = 0.0;
    vn_cos_eq = new double [norder];
    vn_sin_eq = new double [norder];
    vn_cos_vis = new double [norder];
    vn_sin_vis = new double [norder];
    vn_cos_vis_deltaf_restricted = new double [norder];
    vn_sin_vis_deltaf_restricted = new double [norder];
    vn_cos_bulkvis = new double [norder];
    vn_sin_bulkvis = new double [norder];
    vn_cos_bulkvis_deltaf_restricted = new double [norder];
    vn_sin_bulkvis_deltaf_restricted = new double [norder];
    vn_cos_tot = new double [norder];
    vn_sin_tot = new double [norder];
    for(int i = 0; i < norder; i++)
    {
       vn_cos_eq[i] = 0.0;
       vn_sin_eq[i] = 0.0;
       vn_cos_vis[i] = 0.0;
       vn_sin_vis[i] = 0.0;
       vn_cos_vis_deltaf_restricted[i] = 0.0;
       vn_sin_vis_deltaf_restricted[i] = 0.0;
       vn_cos_bulkvis[i] = 0.0;
       vn_sin_bulkvis[i] = 0.0;
       vn_cos_bulkvis_deltaf_restricted[i] = 0.0;
       vn_sin_bulkvis_deltaf_restricted[i] = 0.0;
       vn_cos_tot[i] = 0.0;
       vn_sin_tot[i] = 0.0;
    }


    int diff_flag = paraRdr->getVal("differential_flag");

    if(diff_flag == 1 or diff_flag > 10)
    {
       nTcut = paraRdr->getVal("nTcut");
       nTaucut = paraRdr->getVal("nTaucut");

       Tcut_high = paraRdr->getVal("T_cuthigh");
       Tcut_low = paraRdr->getVal("T_cutlow");
       Taucut_high = paraRdr->getVal("tau_end");
       Taucut_low = paraRdr->getVal("tau_start");

       dNd2pTdphidydTdtau_eq = new double**** [nTcut];
       dNd2pTdphidydTdtau_vis = new double**** [nTcut];
       dNd2pTdphidydTdtau_bulkvis = new double**** [nTcut];
       dNd2pTdphidydTdtau_tot = new double**** [nTcut];

       dNdydTdtau_eq = new double* [nTcut];
       dNdydTdtau_vis = new double* [nTcut];
       dNdydTdtau_bulkvis = new double* [nTcut];
       dNdydTdtau_tot = new double* [nTcut];

       vndTdtau_cos_eq = new double** [nTcut];
       vndTdtau_sin_eq = new double** [nTcut];
       vndTdtau_cos_vis = new double** [nTcut];
       vndTdtau_sin_vis = new double** [nTcut];
       vndTdtau_cos_bulkvis = new double** [nTcut];
       vndTdtau_sin_bulkvis = new double** [nTcut];
       vndTdtau_cos_tot = new double** [nTcut];
       vndTdtau_sin_tot = new double** [nTcut];

       for(int i = 0; i < nTcut; i++)
       {
          dNd2pTdphidydTdtau_eq[i] = new double*** [nTaucut];
          dNd2pTdphidydTdtau_vis[i] = new double*** [nTaucut];
          dNd2pTdphidydTdtau_bulkvis[i] = new double*** [nTaucut];
          dNd2pTdphidydTdtau_tot[i] = new double*** [nTaucut];

          dNdydTdtau_eq[i] = new double [nTaucut];
          dNdydTdtau_vis[i] = new double [nTaucut];
          dNdydTdtau_bulkvis[i] = new double [nTaucut];
          dNdydTdtau_tot[i] = new double [nTaucut];

          vndTdtau_cos_eq[i] = new double* [nTaucut];
          vndTdtau_sin_eq[i] = new double* [nTaucut];
          vndTdtau_cos_vis[i] = new double* [nTaucut];
          vndTdtau_sin_vis[i] = new double* [nTaucut];
          vndTdtau_cos_bulkvis[i] = new double* [nTaucut];
          vndTdtau_sin_bulkvis[i] = new double* [nTaucut];
          vndTdtau_cos_tot[i] = new double* [nTaucut];
          vndTdtau_sin_tot[i] = new double* [nTaucut];

          for(int j = 0; j < nTaucut; j++)
          {
             dNd2pTdphidydTdtau_eq[i][j] = new double** [np];
             dNd2pTdphidydTdtau_vis[i][j] = new double** [np];
             dNd2pTdphidydTdtau_bulkvis[i][j] = new double** [np];
             dNd2pTdphidydTdtau_tot[i][j] = new double** [np];

             vndTdtau_cos_eq[i][j] = new double [norder];
             vndTdtau_sin_eq[i][j] = new double [norder];
             vndTdtau_cos_vis[i][j] = new double [norder];
             vndTdtau_sin_vis[i][j] = new double [norder];
             vndTdtau_cos_bulkvis[i][j] = new double [norder];
             vndTdtau_sin_bulkvis[i][j] = new double [norder];
             vndTdtau_cos_tot[i][j] = new double [norder];
             vndTdtau_sin_tot[i][j] = new double [norder];

             dNdydTdtau_eq[i][j] = 0.0;
             dNdydTdtau_vis[i][j] = 0.0;
             dNdydTdtau_bulkvis[i][j] = 0.0;
             dNdydTdtau_tot[i][j] = 0.0;

             for(int jj = 0; jj < norder; jj++)
             {
                vndTdtau_cos_eq[i][j][jj] = 0.0;
                vndTdtau_sin_eq[i][j][jj] = 0.0;
                vndTdtau_cos_vis[i][j][jj] = 0.0;
                vndTdtau_sin_vis[i][j][jj] = 0.0;
                vndTdtau_cos_bulkvis[i][j][jj] = 0.0;
                vndTdtau_sin_bulkvis[i][j][jj] = 0.0;
                vndTdtau_cos_tot[i][j][jj] = 0.0;
                vndTdtau_sin_tot[i][j][jj] = 0.0;
             }
             for(int k = 0; k < np; k++)
             {
                dNd2pTdphidydTdtau_eq[i][j][k] = new double* [nphi];
                dNd2pTdphidydTdtau_vis[i][j][k] = new double* [nphi];
                dNd2pTdphidydTdtau_bulkvis[i][j][k] = new double* [nphi];
                dNd2pTdphidydTdtau_tot[i][j][k] = new double* [nphi];
                for(int l = 0; l < nphi; l++)
                {
                   dNd2pTdphidydTdtau_eq[i][j][k][l] = new double [nrapidity];
                   dNd2pTdphidydTdtau_vis[i][j][k][l] = new double [nrapidity];
                   dNd2pTdphidydTdtau_bulkvis[i][j][k][l] = new double [nrapidity];
                   dNd2pTdphidydTdtau_tot[i][j][k][l] = new double [nrapidity];
                   for(int m = 0; m < nrapidity; m++)
                   {
                      dNd2pTdphidydTdtau_eq[i][j][k][l][m] = 0.0;
                      dNd2pTdphidydTdtau_vis[i][j][k][l][m] = 0.0;
                      dNd2pTdphidydTdtau_bulkvis[i][j][k][l][m] = 0.0;
                      dNd2pTdphidydTdtau_tot[i][j][k][l][m] = 0.0;
                   }
                }
             }
          }
       }

    }
    
    if(diff_flag == 2 or diff_flag > 10)
    {
       n_xperp_cut = paraRdr->getVal("n_xperp_cut");
       n_tau_cut_xtau = paraRdr->getVal("nTaucut");

       xperp_high = paraRdr->getVal("xperp_cuthigh");
       xperp_low = paraRdr->getVal("xperp_cutlow");
       tau_cut_high = paraRdr->getVal("tau_end");
       tau_cut_low = paraRdr->getVal("tau_start");

       dNd2pTdphidydxperpdtau_eq = new double**** [n_xperp_cut];
       dNd2pTdphidydxperpdtau_vis = new double**** [n_xperp_cut];
       dNd2pTdphidydxperpdtau_bulkvis = new double**** [n_xperp_cut];
       dNd2pTdphidydxperpdtau_tot = new double**** [n_xperp_cut];

       dNdydxperpdtau_eq = new double* [n_xperp_cut];
       dNdydxperpdtau_vis = new double* [n_xperp_cut];
       dNdydxperpdtau_bulkvis = new double* [n_xperp_cut];
       dNdydxperpdtau_tot = new double* [n_xperp_cut];

       vndxperpdtau_cos_eq = new double** [n_xperp_cut];
       vndxperpdtau_sin_eq = new double** [n_xperp_cut];
       vndxperpdtau_cos_vis = new double** [n_xperp_cut];
       vndxperpdtau_sin_vis = new double** [n_xperp_cut];
       vndxperpdtau_cos_bulkvis = new double** [n_xperp_cut];
       vndxperpdtau_sin_bulkvis = new double** [n_xperp_cut];
       vndxperpdtau_cos_tot = new double** [n_xperp_cut];
       vndxperpdtau_sin_tot = new double** [n_xperp_cut];

       for(int i = 0; i < n_xperp_cut; i++)
       {
          dNd2pTdphidydxperpdtau_eq[i] = new double*** [n_tau_cut_xtau];
          dNd2pTdphidydxperpdtau_vis[i] = new double*** [n_tau_cut_xtau];
          dNd2pTdphidydxperpdtau_bulkvis[i] = new double*** [n_tau_cut_xtau];
          dNd2pTdphidydxperpdtau_tot[i] = new double*** [n_tau_cut_xtau];

          dNdydxperpdtau_eq[i] = new double [n_tau_cut_xtau];
          dNdydxperpdtau_vis[i] = new double [n_tau_cut_xtau];
          dNdydxperpdtau_bulkvis[i] = new double [n_tau_cut_xtau];
          dNdydxperpdtau_tot[i] = new double [n_tau_cut_xtau];

          vndxperpdtau_cos_eq[i] = new double* [n_tau_cut_xtau];
          vndxperpdtau_sin_eq[i] = new double* [n_tau_cut_xtau];
          vndxperpdtau_cos_vis[i] = new double* [n_tau_cut_xtau];
          vndxperpdtau_sin_vis[i] = new double* [n_tau_cut_xtau];
          vndxperpdtau_cos_bulkvis[i] = new double* [n_tau_cut_xtau];
          vndxperpdtau_sin_bulkvis[i] = new double* [n_tau_cut_xtau];
          vndxperpdtau_cos_tot[i] = new double* [n_tau_cut_xtau];
          vndxperpdtau_sin_tot[i] = new double* [n_tau_cut_xtau];

          for(int j = 0; j < n_tau_cut_xtau; j++)
          {
             dNd2pTdphidydxperpdtau_eq[i][j] = new double** [np];
             dNd2pTdphidydxperpdtau_vis[i][j] = new double** [np];
             dNd2pTdphidydxperpdtau_bulkvis[i][j] = new double** [np];
             dNd2pTdphidydxperpdtau_tot[i][j] = new double** [np];

             vndxperpdtau_cos_eq[i][j] = new double [norder];
             vndxperpdtau_sin_eq[i][j] = new double [norder];
             vndxperpdtau_cos_vis[i][j] = new double [norder];
             vndxperpdtau_sin_vis[i][j] = new double [norder];
             vndxperpdtau_cos_bulkvis[i][j] = new double [norder];
             vndxperpdtau_sin_bulkvis[i][j] = new double [norder];
             vndxperpdtau_cos_tot[i][j] = new double [norder];
             vndxperpdtau_sin_tot[i][j] = new double [norder];

             dNdydxperpdtau_eq[i][j] = 0.0;
             dNdydxperpdtau_vis[i][j] = 0.0;
             dNdydxperpdtau_bulkvis[i][j] = 0.0;
             dNdydxperpdtau_tot[i][j] = 0.0;

             for(int jj = 0; jj < norder; jj++)
             {
                vndxperpdtau_cos_eq[i][j][jj] = 0.0;
                vndxperpdtau_sin_eq[i][j][jj] = 0.0;
                vndxperpdtau_cos_vis[i][j][jj] = 0.0;
                vndxperpdtau_sin_vis[i][j][jj] = 0.0;
                vndxperpdtau_cos_bulkvis[i][j][jj] = 0.0;
                vndxperpdtau_sin_bulkvis[i][j][jj] = 0.0;
                vndxperpdtau_cos_tot[i][j][jj] = 0.0;
                vndxperpdtau_sin_tot[i][j][jj] = 0.0;
             }
             for(int k = 0; k < np; k++)
             {
                dNd2pTdphidydxperpdtau_eq[i][j][k] = new double* [nphi];
                dNd2pTdphidydxperpdtau_vis[i][j][k] = new double* [nphi];
                dNd2pTdphidydxperpdtau_bulkvis[i][j][k] = new double* [nphi];
                dNd2pTdphidydxperpdtau_tot[i][j][k] = new double* [nphi];
                for(int l = 0; l < nphi; l++)
                {
                   dNd2pTdphidydxperpdtau_eq[i][j][k][l] = new double [nrapidity];
                   dNd2pTdphidydxperpdtau_vis[i][j][k][l] = new double [nrapidity];
                   dNd2pTdphidydxperpdtau_bulkvis[i][j][k][l] = new double [nrapidity];
                   dNd2pTdphidydxperpdtau_tot[i][j][k][l] = new double [nrapidity];
                   for(int m = 0; m < nrapidity; m++)
                   {
                      dNd2pTdphidydxperpdtau_eq[i][j][k][l][m] = 0.0;
                      dNd2pTdphidydxperpdtau_vis[i][j][k][l][m] = 0.0;
                      dNd2pTdphidydxperpdtau_bulkvis[i][j][k][l][m] = 0.0;
                      dNd2pTdphidydxperpdtau_tot[i][j][k][l][m] = 0.0;
                   }
                }
             }
          }
       }
    }
    return;
}

ThermalPhoton::~ThermalPhoton()
{
    int TbsizeX = Photonemission_eqrateTable_ptr->getTbsizeX();
    for(int i=0; i<TbsizeX; i++)
    {
       delete [] Emission_eqrateTb_ptr[i];
       delete [] Emission_viscous_rateTb_ptr[i];
       delete [] Emission_bulkvis_rateTb_ptr[i];
    }
    delete [] Emission_eqrateTb_ptr;
    delete [] Emission_viscous_rateTb_ptr;
    delete [] Emission_bulkvis_rateTb_ptr;
    delete [] EmissionrateTb_Yidxptr;
    delete Photonemission_eqrateTable_ptr;
    delete Photonemission_viscous_rateTable_ptr;

    delete [] p;
    delete [] p_weight;
    delete [] phi;
    delete [] phi_weight;
    delete [] y;
    delete [] theta;

    delete [] dNd2pT_eq;
    delete [] dNd2pT_vis;
    delete [] dNd2pT_vis_deltaf_restricted;
    delete [] dNd2pT_bulkvis;
    delete [] dNd2pT_bulkvis_deltaf_restricted;
    delete [] dNd2pT_tot;

    for(int i = 0; i < np; i++)
    {
       for(int j = 0; j < nphi; j++)
       {
          delete [] dNd2pTdphidy_eq[i][j];
          delete [] dNd2pTdphidy_vis[i][j];
          delete [] dNd2pTdphidy_vis_deltaf_restricted[i][j];
          delete [] dNd2pTdphidy_bulkvis[i][j];
          delete [] dNd2pTdphidy_bulkvis_deltaf_restricted[i][j];
          delete [] dNd2pTdphidy_tot[i][j];
       }
       delete [] dNd2pTdphidy_eq[i];
       delete [] dNd2pTdphidy_vis[i];
       delete [] dNd2pTdphidy_vis_deltaf_restricted[i];
       delete [] dNd2pTdphidy_bulkvis[i];
       delete [] dNd2pTdphidy_bulkvis_deltaf_restricted[i];
       delete [] dNd2pTdphidy_tot[i];
    }
    delete [] dNd2pTdphidy_eq;
    delete [] dNd2pTdphidy_vis;
    delete [] dNd2pTdphidy_vis_deltaf_restricted;
    delete [] dNd2pTdphidy_bulkvis;
    delete [] dNd2pTdphidy_bulkvis_deltaf_restricted;
    delete [] dNd2pTdphidy_tot;

    for(int i = 0; i < norder; i++)
    {
       delete [] vnpT_cos_eq[i];
       delete [] vnpT_sin_eq[i];
       delete [] vnpT_cos_vis[i];
       delete [] vnpT_sin_vis[i];
       delete [] vnpT_cos_vis_deltaf_restricted[i];
       delete [] vnpT_sin_vis_deltaf_restricted[i];
       delete [] vnpT_cos_bulkvis[i];
       delete [] vnpT_sin_bulkvis[i];
       delete [] vnpT_cos_bulkvis_deltaf_restricted[i];
       delete [] vnpT_sin_bulkvis_deltaf_restricted[i];
       delete [] vnpT_cos_tot[i];
       delete [] vnpT_sin_tot[i];
    }
    delete [] vnpT_cos_eq;
    delete [] vnpT_sin_eq;
    delete [] vnpT_cos_vis;
    delete [] vnpT_sin_vis;
    delete [] vnpT_cos_vis_deltaf_restricted;
    delete [] vnpT_sin_vis_deltaf_restricted;
    delete [] vnpT_cos_bulkvis;
    delete [] vnpT_sin_bulkvis;
    delete [] vnpT_cos_bulkvis_deltaf_restricted;
    delete [] vnpT_sin_bulkvis_deltaf_restricted;
    delete [] vnpT_cos_tot;
    delete [] vnpT_sin_tot;

    delete [] vn_cos_eq;
    delete [] vn_sin_eq;
    delete [] vn_cos_vis;
    delete [] vn_sin_vis;
    delete [] vn_cos_vis_deltaf_restricted;
    delete [] vn_sin_vis_deltaf_restricted;
    delete [] vn_cos_bulkvis;
    delete [] vn_sin_bulkvis;
    delete [] vn_cos_bulkvis_deltaf_restricted;
    delete [] vn_sin_bulkvis_deltaf_restricted;
    delete [] vn_cos_tot;
    delete [] vn_sin_tot;

    int diff_flag = paraRdr->getVal("differential_flag");
    if(diff_flag == 1 or diff_flag > 10)
    {
       for(int i = 0; i < nTcut; i++)
       {
          for(int j = 0; j < nTaucut; j++)
          {
             for(int k = 0; k < np; k++)
             {
                for(int l = 0; l < nphi; l++)
                {
                   delete[] dNd2pTdphidydTdtau_eq[i][j][k][l];
                   delete[] dNd2pTdphidydTdtau_vis[i][j][k][l];
                   delete[] dNd2pTdphidydTdtau_bulkvis[i][j][k][l];
                   delete[] dNd2pTdphidydTdtau_tot[i][j][k][l];
                }
                delete[] dNd2pTdphidydTdtau_eq[i][j][k];
                delete[] dNd2pTdphidydTdtau_vis[i][j][k];
                delete[] dNd2pTdphidydTdtau_bulkvis[i][j][k];
                delete[] dNd2pTdphidydTdtau_tot[i][j][k];
             }
             delete[] dNd2pTdphidydTdtau_eq[i][j];
             delete[] dNd2pTdphidydTdtau_vis[i][j];
             delete[] dNd2pTdphidydTdtau_bulkvis[i][j];
             delete[] dNd2pTdphidydTdtau_tot[i][j];
             delete[] vndTdtau_cos_eq[i][j];
             delete[] vndTdtau_sin_eq[i][j];
             delete[] vndTdtau_cos_vis[i][j];
             delete[] vndTdtau_sin_vis[i][j];
             delete[] vndTdtau_cos_bulkvis[i][j];
             delete[] vndTdtau_sin_bulkvis[i][j];
             delete[] vndTdtau_cos_tot[i][j];
             delete[] vndTdtau_sin_tot[i][j];
          }
          delete[] dNd2pTdphidydTdtau_eq[i];
          delete[] dNd2pTdphidydTdtau_vis[i];
          delete[] dNd2pTdphidydTdtau_bulkvis[i];
          delete[] dNd2pTdphidydTdtau_tot[i];

          delete[] dNdydTdtau_eq[i];
          delete[] dNdydTdtau_vis[i];
          delete[] dNdydTdtau_bulkvis[i];
          delete[] dNdydTdtau_tot[i];

          delete[] vndTdtau_cos_eq[i];
          delete[] vndTdtau_sin_eq[i];
          delete[] vndTdtau_cos_vis[i];
          delete[] vndTdtau_sin_bulkvis[i];
          delete[] vndTdtau_cos_tot[i];
          delete[] vndTdtau_sin_tot[i];
       }
       delete[] dNd2pTdphidydTdtau_eq;
       delete[] dNd2pTdphidydTdtau_vis;
       delete[] dNd2pTdphidydTdtau_bulkvis;
       delete[] dNd2pTdphidydTdtau_tot;

       delete[] dNdydTdtau_eq;
       delete[] dNdydTdtau_vis;
       delete[] dNdydTdtau_bulkvis;
       delete[] dNdydTdtau_tot;

       delete[] vndTdtau_cos_eq;
       delete[] vndTdtau_sin_eq;
       delete[] vndTdtau_cos_vis;
       delete[] vndTdtau_sin_vis;
       delete[] vndTdtau_cos_bulkvis;
       delete[] vndTdtau_sin_bulkvis;
       delete[] vndTdtau_cos_tot;
       delete[] vndTdtau_sin_tot;
    }

    if(diff_flag == 2 or diff_flag > 10)
    {
       for(int i = 0; i < n_xperp_cut; i++)
       {
          for(int j = 0; j < n_tau_cut_xtau; j++)
          {
             for(int k = 0; k < np; k++)
             {
                for(int l = 0; l < nphi; l++)
                {
                   delete[] dNd2pTdphidydxperpdtau_eq[i][j][k][l];
                   delete[] dNd2pTdphidydxperpdtau_vis[i][j][k][l];
                   delete[] dNd2pTdphidydxperpdtau_bulkvis[i][j][k][l];
                   delete[] dNd2pTdphidydxperpdtau_tot[i][j][k][l];
                }
                delete[] dNd2pTdphidydxperpdtau_eq[i][j][k];
                delete[] dNd2pTdphidydxperpdtau_vis[i][j][k];
                delete[] dNd2pTdphidydxperpdtau_bulkvis[i][j][k];
                delete[] dNd2pTdphidydxperpdtau_tot[i][j][k];
             }
             delete[] dNd2pTdphidydxperpdtau_eq[i][j];
             delete[] dNd2pTdphidydxperpdtau_vis[i][j];
             delete[] dNd2pTdphidydxperpdtau_bulkvis[i][j];
             delete[] dNd2pTdphidydxperpdtau_tot[i][j];
             delete[] vndxperpdtau_cos_eq[i][j];
             delete[] vndxperpdtau_sin_eq[i][j];
             delete[] vndxperpdtau_cos_vis[i][j];
             delete[] vndxperpdtau_sin_vis[i][j];
             delete[] vndxperpdtau_cos_bulkvis[i][j];
             delete[] vndxperpdtau_sin_bulkvis[i][j];
             delete[] vndxperpdtau_cos_tot[i][j];
             delete[] vndxperpdtau_sin_tot[i][j];
          }
          delete[] dNd2pTdphidydxperpdtau_eq[i];
          delete[] dNd2pTdphidydxperpdtau_vis[i];
          delete[] dNd2pTdphidydxperpdtau_bulkvis[i];
          delete[] dNd2pTdphidydxperpdtau_tot[i];

          delete[] dNdydxperpdtau_eq[i];
          delete[] dNdydxperpdtau_vis[i];
          delete[] dNdydxperpdtau_bulkvis[i];
          delete[] dNdydxperpdtau_tot[i];

          delete[] vndxperpdtau_cos_eq[i];
          delete[] vndxperpdtau_sin_eq[i];
          delete[] vndxperpdtau_cos_vis[i];
          delete[] vndxperpdtau_sin_bulkvis[i];
          delete[] vndxperpdtau_cos_tot[i];
          delete[] vndxperpdtau_sin_tot[i];
       }
       delete[] dNd2pTdphidydxperpdtau_eq;
       delete[] dNd2pTdphidydxperpdtau_vis;
       delete[] dNd2pTdphidydxperpdtau_bulkvis;
       delete[] dNd2pTdphidydxperpdtau_tot;

       delete[] dNdydxperpdtau_eq;
       delete[] dNdydxperpdtau_vis;
       delete[] dNdydxperpdtau_bulkvis;
       delete[] dNdydxperpdtau_tot;

       delete[] vndxperpdtau_cos_eq;
       delete[] vndxperpdtau_sin_eq;
       delete[] vndxperpdtau_cos_vis;
       delete[] vndxperpdtau_sin_vis;
       delete[] vndxperpdtau_cos_bulkvis;
       delete[] vndxperpdtau_sin_bulkvis;
       delete[] vndxperpdtau_cos_tot;
       delete[] vndxperpdtau_sin_tot;
    }
}

void ThermalPhoton::setupEmissionrate(string emissionProcess, double Xmin, double dX,  double Ymin, double dY)
{
     EmissionrateTb_Xmin = Xmin;
     EmissionrateTb_Ymin = Ymin;
     EmissionrateTb_dX = dX;
     EmissionrateTb_dY = dY;
     readEmissionrate(emissionProcess);
}

void ThermalPhoton::readEmissionrate(string emissionProcess)
{
    emissionProcess_name = emissionProcess;
    ostringstream eqrate_filename_stream;
    ostringstream visrate_filename_stream;
    ostringstream bulkvisrate_filename_stream;
    eqrate_filename_stream << rate_path << "rate_" << emissionProcess << "_eqrate.dat";
    visrate_filename_stream << rate_path << "rate_" << emissionProcess << "_viscous.dat";
    bulkvisrate_filename_stream << rate_path << "rate_" << emissionProcess << "_bulkvis.dat";
    Photonemission_eqrateTable_ptr = new Table2D(eqrate_filename_stream.str().c_str());
    Photonemission_viscous_rateTable_ptr = new Table2D(visrate_filename_stream.str().c_str());
    Photonemission_bulkvis_rateTable_ptr = new Table2D(bulkvisrate_filename_stream.str().c_str());

    EmissionrateTb_sizeX = Photonemission_eqrateTable_ptr->getTbsizeX();
    EmissionrateTb_sizeY = Photonemission_eqrateTable_ptr->getTbsizeY();
    Emission_eqrateTb_ptr = new double* [EmissionrateTb_sizeX];
    Emission_viscous_rateTb_ptr = new double* [EmissionrateTb_sizeX];
    Emission_bulkvis_rateTb_ptr = new double* [EmissionrateTb_sizeX];
    for(int i=0; i<EmissionrateTb_sizeX; i++)
    {
       Emission_eqrateTb_ptr[i] = new double [EmissionrateTb_sizeY];
       Emission_viscous_rateTb_ptr[i] = new double [EmissionrateTb_sizeY];
       Emission_bulkvis_rateTb_ptr[i] = new double [EmissionrateTb_sizeY];
    }

    EmissionrateTb_Yidxptr = new double [EmissionrateTb_sizeY];
    for(int i=0; i<EmissionrateTb_sizeY; i++)
       EmissionrateTb_Yidxptr[i] = EmissionrateTb_Ymin + i*EmissionrateTb_dY;

    // take log for the equilibrium emission rates for better interpolation precisions; take ratio of viscous rates w.r.t eq. rates for faster performance
    for(int i=0; i<EmissionrateTb_sizeX; i++)
       for(int j=0; j<EmissionrateTb_sizeY; j++)
       {
           Emission_eqrateTb_ptr[i][j] = log(Photonemission_eqrateTable_ptr->getTbdata(i,j) + 1e-30);  
           Emission_viscous_rateTb_ptr[i][j] = Photonemission_viscous_rateTable_ptr->getTbdata(i,j)/(Photonemission_eqrateTable_ptr->getTbdata(i,j) + 1e-30);
           Emission_bulkvis_rateTb_ptr[i][j] = Photonemission_bulkvis_rateTable_ptr->getTbdata(i,j)/(Photonemission_eqrateTable_ptr->getTbdata(i,j) + 1e-30);
       }

    return;
}

void ThermalPhoton::getPhotonemissionRate(double* Eq, double* pi_zz, double* bulkPi, int Eq_length, double T, double* eqrate_ptr, double* visrate_ptr, double* bulkvis_ptr)
{
    //interpolate equilibrium rate
    interpolation2D_bilinear(T, Eq, Eq_length, Emission_eqrateTb_ptr, eqrate_ptr);
    //interpolate viscous rate
    interpolation2D_bilinear(T, Eq, Eq_length, Emission_viscous_rateTb_ptr, visrate_ptr);
    //interpolate bulk viscous rate
    interpolation2D_bilinear(T, Eq, Eq_length, Emission_bulkvis_rateTb_ptr, bulkvis_ptr);
    for(int i=0; i<Eq_length; i++)
    {
       eqrate_ptr[i] = exp(eqrate_ptr[i]);
       visrate_ptr[i] = pi_zz[i]*visrate_ptr[i]*eqrate_ptr[i];
       bulkvis_ptr[i] = bulkPi[i]*bulkvis_ptr[i]*eqrate_ptr[i];
    }

    return;
}

void ThermalPhoton::calThermalPhotonemission(double* Eq, double* pi_zz, double* bulkPi, int Tb_length, double T, double* volume, double fraction)
{
    double* em_eqrate = new double [Tb_length];   //photon emission equilibrium rate at local rest cell
    double* em_visrate = new double [Tb_length];   //photon emission viscous correction at local rest cell
    double* em_bulkvis = new double [Tb_length];   //photon emission bulk viscous correction at local rest cell
    getPhotonemissionRate(Eq, pi_zz, bulkPi, Tb_length, T, em_eqrate, em_visrate, em_bulkvis);

    int n_pt_point = nrapidity*np*nphi;
    
    double temp_eq_sum, temp_vis_sum, temp_bulkvis_sum;
    double temp_vis_deltaf_restricted_sum, temp_bulkvis_deltaf_restricted_sum;
    int idx=0;
    for(int k=0; k<nrapidity; k++)
    {
      for(int m=0; m<nphi; m++)
      {
        for(int l=0; l<np; l++)
        {
          temp_eq_sum = 0.0;
          temp_vis_sum = 0.0;
          temp_bulkvis_sum = 0.0;
          temp_vis_deltaf_restricted_sum = 0.0;
          temp_bulkvis_deltaf_restricted_sum = 0.0;
          for(int i=0; i < neta; i++)
          {
             double local_eq = em_eqrate[idx + i*n_pt_point];
             double local_vis = em_visrate[idx + i*n_pt_point];
             double local_bulk = em_bulkvis[idx + i*n_pt_point];

             temp_eq_sum += local_eq*volume[i]*fraction;
             temp_vis_sum += local_vis*volume[i]*fraction;
             temp_bulkvis_sum += local_bulk*volume[i]*fraction;
             
             double ratio = fabs(local_vis + local_bulk)/(local_eq + 1e-20);
             if(ratio > 1.0)
             {
                 local_vis = local_vis/ratio;
                 local_bulk = local_bulk/ratio;
             }
             
             temp_vis_deltaf_restricted_sum += local_vis*volume[i]*fraction;
             temp_bulkvis_deltaf_restricted_sum += local_bulk*volume[i]*fraction;
          }
          dNd2pTdphidy_eq[l][m][k] += temp_eq_sum;
          dNd2pTdphidy_vis[l][m][k] += temp_eq_sum + temp_vis_sum;
          dNd2pTdphidy_vis_deltaf_restricted[l][m][k] += temp_eq_sum + temp_vis_deltaf_restricted_sum;
          dNd2pTdphidy_bulkvis[l][m][k] += temp_eq_sum + temp_bulkvis_sum;
          dNd2pTdphidy_bulkvis_deltaf_restricted[l][m][k] += temp_eq_sum + temp_bulkvis_deltaf_restricted_sum;
          dNd2pTdphidy_tot[l][m][k] += temp_eq_sum + temp_vis_sum + temp_bulkvis_sum;
          idx++;
        }
      }
    }
    delete[] em_eqrate;
    delete[] em_visrate;
    delete[] em_bulkvis;
    /*if(isnan(em_rate))
    {
        cout<<" em_rate is nan" << endl;
        exit(1);
    }*/
    return;
}

void ThermalPhoton::calThermalPhotonemissiondTdtau(double* Eq, double* pi_zz, double* bulkPi, int Tb_length, double T, double tau, double* volume, double fraction)
{
    double* em_eqrate = new double [Tb_length];   //photon emission equilibrium rate at local rest cell
    double* em_visrate = new double [Tb_length];   //photon emission viscous correction at local rest cell
    double* em_bulkvis = new double [Tb_length];   //photon emission bulk viscous correction at local rest cell
    getPhotonemissionRate(Eq, pi_zz, bulkPi, Tb_length, T, em_eqrate, em_visrate, em_bulkvis);

    int n_pt_point = nrapidity*np*nphi;

    double eps = 1e-15;
    double dT = (Tcut_high - Tcut_low)/(nTcut - 1);
    double dtau = (Taucut_high - Taucut_low)/(nTaucut - 1);
    int idx_T = (int)((T - Tcut_low)/dT + eps);
    int idx_tau = (int)((tau - Taucut_low)/dtau + eps);

    double temp_eq_sum, temp_vis_sum, temp_bulkvis_sum;
    int idx=0;
    for(int k=0; k<nrapidity; k++)
    {
      for(int m=0; m<nphi; m++)
      {
        for(int l=0; l<np; l++)
        {
          temp_eq_sum = 0.0;
          temp_vis_sum = 0.0;
          temp_bulkvis_sum = 0.0;
          for(int i=0; i < neta; i++)
          {
             temp_eq_sum += em_eqrate[idx + i*n_pt_point]*volume[i]*fraction;
             temp_vis_sum += em_visrate[idx + i*n_pt_point]*volume[i]*fraction;
             temp_bulkvis_sum += em_bulkvis[idx + i*n_pt_point]*volume[i]*fraction;
          }
          dNd2pTdphidydTdtau_eq[idx_T][idx_tau][l][m][k] += temp_eq_sum;
          dNd2pTdphidydTdtau_vis[idx_T][idx_tau][l][m][k] += temp_eq_sum + temp_vis_sum;
          dNd2pTdphidydTdtau_bulkvis[idx_T][idx_tau][l][m][k] += temp_eq_sum + temp_bulkvis_sum;
          dNd2pTdphidydTdtau_tot[idx_T][idx_tau][l][m][k] += temp_eq_sum + temp_vis_sum + temp_bulkvis_sum;
          idx++;
        }
      }
    }

    delete[] em_eqrate;
    delete[] em_visrate;
    delete[] em_bulkvis;
    return;
}

void ThermalPhoton::calThermalPhotonemissiondxperpdtau(double* Eq, double* pi_zz, double* bulkPi, int Tb_length, double T, double x_local, double tau, double* volume, double fraction)
{
    double* em_eqrate = new double [Tb_length];   //photon emission equilibrium rate at local rest cell
    double* em_visrate = new double [Tb_length];   //photon emission viscous correction at local rest cell
    double* em_bulkvis = new double [Tb_length];   //photon emission bulk viscous correction at local rest cell
    getPhotonemissionRate(Eq, pi_zz, bulkPi, Tb_length, T, em_eqrate, em_visrate, em_bulkvis);

    int n_pt_point = nrapidity*np*nphi;

    double eps = 1e-15;
    double dxperp = (xperp_high - xperp_low)/(n_xperp_cut - 1);
    double dtau = (tau_cut_high - tau_cut_low)/(n_tau_cut_xtau - 1);
    int idx_xperp = (int)((x_local - xperp_low)/dxperp + eps);
    int idx_tau = (int)((tau - tau_cut_low)/dtau + eps);

    double temp_eq_sum, temp_vis_sum, temp_bulkvis_sum;
    int idx=0;
    for(int k=0; k<nrapidity; k++)
    {
      for(int m=0; m<nphi; m++)
      {
        for(int l=0; l<np; l++)
        {
          temp_eq_sum = 0.0;
          temp_vis_sum = 0.0;
          temp_bulkvis_sum = 0.0;
          for(int i=0; i < neta; i++)
          {
             temp_eq_sum += em_eqrate[idx + i*n_pt_point]*volume[i]*fraction;
             temp_vis_sum += em_visrate[idx + i*n_pt_point]*volume[i]*fraction;
             temp_bulkvis_sum += em_bulkvis[idx + i*n_pt_point]*volume[i]*fraction;
          }
          dNd2pTdphidydxperpdtau_eq[idx_xperp][idx_tau][l][m][k] += temp_eq_sum;
          dNd2pTdphidydxperpdtau_vis[idx_xperp][idx_tau][l][m][k] += temp_eq_sum + temp_vis_sum;
          dNd2pTdphidydxperpdtau_bulkvis[idx_xperp][idx_tau][l][m][k] += temp_eq_sum + temp_bulkvis_sum;
          dNd2pTdphidydxperpdtau_tot[idx_xperp][idx_tau][l][m][k] += temp_eq_sum + temp_vis_sum + temp_bulkvis_sum;
          idx++;
        }
      }
    }

    delete[] em_eqrate;
    delete[] em_visrate;
    delete[] em_bulkvis;
    return;
}


void ThermalPhoton::calPhoton_SpvnpT()
//calculate the photon spectra and differential vn at mid-rapidity
{
   int k = 0;  //calculate at y = 0
   double eps = 1e-15;
   for(int i=0;i<np;i++)
   {
       for(int j=0;j<nphi;j++)
       {
         dNd2pT_eq[i] += dNd2pTdphidy_eq[i][j][k]*phi_weight[j];  //dN/(dy pT dpT)
         dNd2pT_vis[i] += dNd2pTdphidy_vis[i][j][k]*phi_weight[j];
         dNd2pT_vis_deltaf_restricted[i] += dNd2pTdphidy_vis_deltaf_restricted[i][j][k]*phi_weight[j];
         dNd2pT_bulkvis[i] += dNd2pTdphidy_bulkvis[i][j][k]*phi_weight[j];
         dNd2pT_bulkvis_deltaf_restricted[i] += dNd2pTdphidy_bulkvis_deltaf_restricted[i][j][k]*phi_weight[j];
         dNd2pT_tot[i] += dNd2pTdphidy_tot[i][j][k]*phi_weight[j];
         for(int order=0; order<norder; order++)
         {
            vnpT_cos_eq[order][i] += dNd2pTdphidy_eq[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_cos_vis[order][i] += dNd2pTdphidy_vis[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_cos_vis_deltaf_restricted[order][i] += dNd2pTdphidy_vis_deltaf_restricted[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_cos_bulkvis[order][i] += dNd2pTdphidy_bulkvis[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_cos_bulkvis_deltaf_restricted[order][i] += dNd2pTdphidy_bulkvis_deltaf_restricted[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_cos_tot[order][i] += dNd2pTdphidy_tot[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_sin_eq[order][i] += dNd2pTdphidy_eq[i][j][k]*sin(order*phi[j])*phi_weight[j];
            vnpT_sin_vis[order][i] += dNd2pTdphidy_vis[i][j][k]*sin(order*phi[j])*phi_weight[j];
            vnpT_sin_vis_deltaf_restricted[order][i] += dNd2pTdphidy_vis_deltaf_restricted[i][j][k]*sin(order*phi[j])*phi_weight[j];
            vnpT_sin_bulkvis[order][i] += dNd2pTdphidy_bulkvis[i][j][k]*sin(order*phi[j])*phi_weight[j];
            vnpT_sin_bulkvis_deltaf_restricted[order][i] += dNd2pTdphidy_bulkvis_deltaf_restricted[i][j][k]*sin(order*phi[j])*phi_weight[j];
            vnpT_sin_tot[order][i] += dNd2pTdphidy_tot[i][j][k]*sin(order*phi[j])*phi_weight[j];
         }
       }
       dNdy_eq += dNd2pT_eq[i]*p[i]*p_weight[i];   //dN/dy
       dNdy_vis += dNd2pT_vis[i]*p[i]*p_weight[i];
       dNdy_vis_deltaf_restricted += dNd2pT_vis_deltaf_restricted[i]*p[i]*p_weight[i];
       dNdy_bulkvis += dNd2pT_bulkvis[i]*p[i]*p_weight[i];
       dNdy_bulkvis_deltaf_restricted += dNd2pT_bulkvis_deltaf_restricted[i]*p[i]*p_weight[i];
       dNdy_tot += dNd2pT_tot[i]*p[i]*p_weight[i];
       for(int order = 0; order < norder ; order++)
       {
          vn_cos_eq[order] += vnpT_cos_eq[order][i]*p[i]*p_weight[i];
          vn_sin_eq[order] += vnpT_sin_eq[order][i]*p[i]*p_weight[i];
          vn_cos_vis[order] += vnpT_cos_vis[order][i]*p[i]*p_weight[i];
          vn_sin_vis[order] += vnpT_sin_vis[order][i]*p[i]*p_weight[i];
          vn_cos_vis_deltaf_restricted[order] += vnpT_cos_vis_deltaf_restricted[order][i]*p[i]*p_weight[i];
          vn_sin_vis_deltaf_restricted[order] += vnpT_sin_vis_deltaf_restricted[order][i]*p[i]*p_weight[i];
          vn_cos_bulkvis[order] += vnpT_cos_bulkvis[order][i]*p[i]*p_weight[i];
          vn_sin_bulkvis[order] += vnpT_sin_bulkvis[order][i]*p[i]*p_weight[i];
          vn_cos_bulkvis_deltaf_restricted[order] += vnpT_cos_bulkvis_deltaf_restricted[order][i]*p[i]*p_weight[i];
          vn_sin_bulkvis_deltaf_restricted[order] += vnpT_sin_bulkvis_deltaf_restricted[order][i]*p[i]*p_weight[i];
          vn_cos_tot[order] += vnpT_cos_tot[order][i]*p[i]*p_weight[i];
          vn_sin_tot[order] += vnpT_sin_tot[order][i]*p[i]*p_weight[i];
          
          //vn(pT)
          vnpT_cos_eq[order][i] = vnpT_cos_eq[order][i]/dNd2pT_eq[i];
          vnpT_cos_vis[order][i] = vnpT_cos_vis[order][i]/dNd2pT_vis[i];
          vnpT_cos_vis_deltaf_restricted[order][i] = vnpT_cos_vis_deltaf_restricted[order][i]/dNd2pT_vis_deltaf_restricted[i];
          vnpT_cos_bulkvis[order][i] = vnpT_cos_bulkvis[order][i]/dNd2pT_bulkvis[i];
          vnpT_cos_bulkvis_deltaf_restricted[order][i] = vnpT_cos_bulkvis_deltaf_restricted[order][i]/dNd2pT_bulkvis_deltaf_restricted[i];
          vnpT_cos_tot[order][i] = vnpT_cos_tot[order][i]/dNd2pT_tot[i];
          vnpT_sin_eq[order][i] = vnpT_sin_eq[order][i]/dNd2pT_eq[i];
          vnpT_sin_vis[order][i] = vnpT_sin_vis[order][i]/dNd2pT_vis[i];
          vnpT_sin_vis_deltaf_restricted[order][i] = vnpT_sin_vis_deltaf_restricted[order][i]/dNd2pT_vis_deltaf_restricted[i];
          vnpT_sin_bulkvis[order][i] = vnpT_sin_bulkvis[order][i]/dNd2pT_bulkvis[i];
          vnpT_sin_bulkvis_deltaf_restricted[order][i] = vnpT_sin_bulkvis_deltaf_restricted[order][i]/dNd2pT_bulkvis_deltaf_restricted[i];
          vnpT_sin_tot[order][i] = vnpT_sin_tot[order][i]/dNd2pT_tot[i];
       }
       dNd2pT_eq[i] = dNd2pT_eq[i]/(2*M_PI);  //dN/(2pi dy pT dpT)
       dNd2pT_vis[i] = dNd2pT_vis[i]/(2*M_PI);
       dNd2pT_vis_deltaf_restricted[i] = dNd2pT_vis_deltaf_restricted[i]/(2*M_PI);
       dNd2pT_bulkvis[i] = dNd2pT_bulkvis[i]/(2*M_PI);
       dNd2pT_bulkvis_deltaf_restricted[i] = dNd2pT_bulkvis_deltaf_restricted[i]/(2*M_PI);
       dNd2pT_tot[i] = dNd2pT_tot[i]/(2*M_PI);
   }
   for(int order = 1; order < norder ; order++)
   {
       //vn
       vn_cos_eq[order] = vn_cos_eq[order]/(dNdy_eq + eps);
       vn_sin_eq[order] = vn_sin_eq[order]/(dNdy_eq + eps);
       vn_cos_vis[order] = vn_cos_vis[order]/(dNdy_vis + eps);
       vn_sin_vis[order] = vn_sin_vis[order]/(dNdy_vis + eps);
       vn_cos_vis_deltaf_restricted[order] = vn_cos_vis_deltaf_restricted[order]/(dNdy_vis_deltaf_restricted + eps);
       vn_sin_vis_deltaf_restricted[order] = vn_sin_vis_deltaf_restricted[order]/(dNdy_vis_deltaf_restricted + eps);
       vn_cos_bulkvis[order] = vn_cos_bulkvis[order]/(dNdy_vis + eps);
       vn_sin_bulkvis[order] = vn_sin_bulkvis[order]/(dNdy_vis + eps);
       vn_cos_bulkvis_deltaf_restricted[order] = vn_cos_bulkvis_deltaf_restricted[order]/(dNdy_vis_deltaf_restricted + eps);
       vn_sin_bulkvis_deltaf_restricted[order] = vn_sin_bulkvis_deltaf_restricted[order]/(dNdy_vis_deltaf_restricted + eps);
       vn_cos_tot[order] = vn_cos_tot[order]/(dNdy_tot + eps);
       vn_sin_tot[order] = vn_sin_tot[order]/(dNdy_tot + eps);
   }
   return;
}

void ThermalPhoton::output_photon_spectra_dTdtau(string path)
//calculate the inverse slope of the photon spectra at T-tau interval
{
   ostringstream filename_SpdTdtau_eq;
   ostringstream filename_SpdTdtau_vis;
   ostringstream filename_SpdTdtau_bulkvis;
   ostringstream filename_SpdTdtau_tot;
   filename_SpdTdtau_eq << path << emissionProcess_name << "_SpdTdtau_eq.dat";
   filename_SpdTdtau_vis << path << emissionProcess_name << "_SpdTdtau_vis.dat";
   filename_SpdTdtau_bulkvis << path << emissionProcess_name << "_SpdTdtau_bulkvis.dat";
   filename_SpdTdtau_tot << path << emissionProcess_name << "_SpdTdtau_tot.dat";

   ofstream ofeq(filename_SpdTdtau_eq.str().c_str());
   ofstream ofvis(filename_SpdTdtau_vis.str().c_str());
   ofstream ofbulkvis(filename_SpdTdtau_bulkvis.str().c_str());
   ofstream oftot(filename_SpdTdtau_tot.str().c_str());

   int irap = 0;  //calculate at y = 0
   double dT = (Tcut_high - Tcut_low)/(nTcut - 1);
   double dtau = (Taucut_high - Taucut_low)/(nTaucut - 1);
   for(int i = 0; i < nTcut; i++)
   {
       double T_local = Tcut_low + i*dT;
       for(int j = 0; j < nTaucut; j++)
       {
           double tau_local = Taucut_low + j*dtau;

           ofeq << scientific << setw(18) << setprecision(8) 
                << T_local << "   "  << tau_local << "   ";
           ofvis << scientific << setw(18) << setprecision(8) 
                 << T_local << "   "  << tau_local << "   ";
           ofbulkvis << scientific << setw(18) << setprecision(8) 
                     << T_local << "   "  << tau_local << "   ";
           oftot << scientific << setw(18) << setprecision(8) 
                 << T_local << "   "  << tau_local << "   ";
           for(int k = 0; k < np; k++)
           {
               double temp_dNdypTdpT_eq = 0.0;
               double temp_dNdypTdpT_vis = 0.0;
               double temp_dNdypTdpT_bulkvis = 0.0;
               double temp_dNdypTdpT_tot = 0.0;
               for(int l = 0; l < nphi; l++)
               {
                   temp_dNdypTdpT_eq += dNd2pTdphidydTdtau_eq[i][j][k][l][irap]*phi_weight[l];
                   temp_dNdypTdpT_vis += dNd2pTdphidydTdtau_vis[i][j][k][l][irap]*phi_weight[l];
                   temp_dNdypTdpT_bulkvis += dNd2pTdphidydTdtau_bulkvis[i][j][k][l][irap]*phi_weight[l];
                   temp_dNdypTdpT_tot += dNd2pTdphidydTdtau_tot[i][j][k][l][irap]*phi_weight[l];
               }
               ofeq << temp_dNdypTdpT_eq << "   ";
               ofvis << temp_dNdypTdpT_vis << "   ";
               ofbulkvis << temp_dNdypTdpT_bulkvis << "   ";
               oftot << temp_dNdypTdpT_tot << "   ";
           }
           ofeq << endl;
           ofvis << endl;
           ofbulkvis << endl;
           oftot << endl;
       }
   }
   ofeq.close();
   ofvis.close();
   ofbulkvis.close();
   oftot.close();
   return;
}

void ThermalPhoton::calPhoton_SpvnpT_dTdtau()
//calculate the photon spectra and differential vn at mid-rapidity
{
   int irap = 0;  //calculate at y = 0
   double eps = 1e-15;
   for(int i = 0; i < nTcut; i++)
   {
       for(int j = 0; j < nTaucut; j++)
       {
          for(int k = 0; k < np; k++)
          {
             for(int l = 0; l < nphi; l++)
             {
                dNdydTdtau_eq[i][j] += dNd2pTdphidydTdtau_eq[i][j][k][l][irap]*p[k]*p_weight[k]*phi_weight[l];
                dNdydTdtau_vis[i][j] += dNd2pTdphidydTdtau_vis[i][j][k][l][irap]*p[k]*p_weight[k]*phi_weight[l];
                dNdydTdtau_bulkvis[i][j] += dNd2pTdphidydTdtau_bulkvis[i][j][k][l][irap]*p[k]*p_weight[k]*phi_weight[l];
                dNdydTdtau_tot[i][j] += dNd2pTdphidydTdtau_tot[i][j][k][l][irap]*p[k]*p_weight[k]*phi_weight[l];
                for(int order=0; order<norder; order++)
                {
                   vndTdtau_cos_eq[i][j][order] += dNd2pTdphidydTdtau_eq[i][j][k][l][irap]*p[k]*p_weight[k]*cos(order*phi[l])*phi_weight[l];
                   vndTdtau_sin_eq[i][j][order] += dNd2pTdphidydTdtau_eq[i][j][k][l][irap]*p[k]*p_weight[k]*sin(order*phi[l])*phi_weight[l];
                   vndTdtau_cos_vis[i][j][order] += dNd2pTdphidydTdtau_vis[i][j][k][l][irap]*p[k]*p_weight[k]*cos(order*phi[l])*phi_weight[l];
                   vndTdtau_sin_vis[i][j][order] += dNd2pTdphidydTdtau_vis[i][j][k][l][irap]*p[k]*p_weight[k]*sin(order*phi[l])*phi_weight[l];
                   vndTdtau_cos_bulkvis[i][j][order] += dNd2pTdphidydTdtau_bulkvis[i][j][k][l][irap]*p[k]*p_weight[k]*cos(order*phi[l])*phi_weight[l];
                   vndTdtau_sin_bulkvis[i][j][order] += dNd2pTdphidydTdtau_bulkvis[i][j][k][l][irap]*p[k]*p_weight[k]*sin(order*phi[l])*phi_weight[l];
                   vndTdtau_cos_tot[i][j][order] += dNd2pTdphidydTdtau_tot[i][j][k][l][irap]*p[k]*p_weight[k]*cos(order*phi[l])*phi_weight[l];
                   vndTdtau_sin_tot[i][j][order] += dNd2pTdphidydTdtau_tot[i][j][k][l][irap]*p[k]*p_weight[k]*sin(order*phi[l])*phi_weight[l];
                }
             }
          }
          for(int order = 1; order < norder ; order++)
          {
              vndTdtau_cos_eq[i][j][order] = vndTdtau_cos_eq[i][j][order]/(dNdydTdtau_eq[i][j] + eps);
              vndTdtau_sin_eq[i][j][order] = vndTdtau_sin_eq[i][j][order]/(dNdydTdtau_eq[i][j] + eps);
              vndTdtau_cos_vis[i][j][order] = vndTdtau_cos_vis[i][j][order]/(dNdydTdtau_vis[i][j] + eps);
              vndTdtau_sin_vis[i][j][order] = vndTdtau_sin_vis[i][j][order]/(dNdydTdtau_vis[i][j] + eps);
              vndTdtau_cos_bulkvis[i][j][order] = vndTdtau_cos_bulkvis[i][j][order]/(dNdydTdtau_bulkvis[i][j] + eps);
              vndTdtau_sin_bulkvis[i][j][order] = vndTdtau_sin_bulkvis[i][j][order]/(dNdydTdtau_bulkvis[i][j] + eps);
              vndTdtau_cos_tot[i][j][order] = vndTdtau_cos_tot[i][j][order]/(dNdydTdtau_tot[i][j] + eps);
              vndTdtau_sin_tot[i][j][order] = vndTdtau_sin_tot[i][j][order]/(dNdydTdtau_tot[i][j] + eps);
          }
       }
   }
   return;
}

void ThermalPhoton::calPhoton_SpvnpT_dxperpdtau()
//calculate the photon spectra and differential vn at mid-rapidity
{
   int irap = 0;  //calculate at y = 0
   double eps = 1e-15;
   for(int i = 0; i < n_xperp_cut; i++)
   {
       for(int j = 0; j < n_tau_cut_xtau; j++)
       {
          for(int k = 0; k < np; k++)
          {
             for(int l = 0; l < nphi; l++)
             {
                dNdydxperpdtau_eq[i][j] += dNd2pTdphidydxperpdtau_eq[i][j][k][l][irap]*p[k]*p_weight[k]*phi_weight[l];
                dNdydxperpdtau_vis[i][j] += dNd2pTdphidydxperpdtau_vis[i][j][k][l][irap]*p[k]*p_weight[k]*phi_weight[l];
                dNdydxperpdtau_bulkvis[i][j] += dNd2pTdphidydxperpdtau_bulkvis[i][j][k][l][irap]*p[k]*p_weight[k]*phi_weight[l];
                dNdydxperpdtau_tot[i][j] += dNd2pTdphidydxperpdtau_tot[i][j][k][l][irap]*p[k]*p_weight[k]*phi_weight[l];
                for(int order=0; order<norder; order++)
                {
                   vndxperpdtau_cos_eq[i][j][order] += dNd2pTdphidydxperpdtau_eq[i][j][k][l][irap]*p[k]*p_weight[k]*cos(order*phi[l])*phi_weight[l];
                   vndxperpdtau_sin_eq[i][j][order] += dNd2pTdphidydxperpdtau_eq[i][j][k][l][irap]*p[k]*p_weight[k]*sin(order*phi[l])*phi_weight[l];
                   vndxperpdtau_cos_vis[i][j][order] += dNd2pTdphidydxperpdtau_vis[i][j][k][l][irap]*p[k]*p_weight[k]*cos(order*phi[l])*phi_weight[l];
                   vndxperpdtau_sin_vis[i][j][order] += dNd2pTdphidydxperpdtau_vis[i][j][k][l][irap]*p[k]*p_weight[k]*sin(order*phi[l])*phi_weight[l];
                   vndxperpdtau_cos_bulkvis[i][j][order] += dNd2pTdphidydxperpdtau_bulkvis[i][j][k][l][irap]*p[k]*p_weight[k]*cos(order*phi[l])*phi_weight[l];
                   vndxperpdtau_sin_bulkvis[i][j][order] += dNd2pTdphidydxperpdtau_bulkvis[i][j][k][l][irap]*p[k]*p_weight[k]*sin(order*phi[l])*phi_weight[l];
                   vndxperpdtau_cos_tot[i][j][order] += dNd2pTdphidydxperpdtau_tot[i][j][k][l][irap]*p[k]*p_weight[k]*cos(order*phi[l])*phi_weight[l];
                   vndxperpdtau_sin_tot[i][j][order] += dNd2pTdphidydxperpdtau_tot[i][j][k][l][irap]*p[k]*p_weight[k]*sin(order*phi[l])*phi_weight[l];
                }
             }
          }
          for(int order = 1; order < norder ; order++)
          {
              vndxperpdtau_cos_eq[i][j][order] = vndxperpdtau_cos_eq[i][j][order]/(dNdydxperpdtau_eq[i][j] + eps);
              vndxperpdtau_sin_eq[i][j][order] = vndxperpdtau_sin_eq[i][j][order]/(dNdydxperpdtau_eq[i][j] + eps);
              vndxperpdtau_cos_vis[i][j][order] = vndxperpdtau_cos_vis[i][j][order]/(dNdydxperpdtau_vis[i][j] + eps);
              vndxperpdtau_sin_vis[i][j][order] = vndxperpdtau_sin_vis[i][j][order]/(dNdydxperpdtau_vis[i][j] + eps);
              vndxperpdtau_cos_bulkvis[i][j][order] = vndxperpdtau_cos_bulkvis[i][j][order]/(dNdydxperpdtau_bulkvis[i][j] + eps);
              vndxperpdtau_sin_bulkvis[i][j][order] = vndxperpdtau_sin_bulkvis[i][j][order]/(dNdydxperpdtau_bulkvis[i][j] + eps);
              vndxperpdtau_cos_tot[i][j][order] = vndxperpdtau_cos_tot[i][j][order]/(dNdydxperpdtau_tot[i][j] + eps);
              vndxperpdtau_sin_tot[i][j][order] = vndxperpdtau_sin_tot[i][j][order]/(dNdydxperpdtau_tot[i][j] + eps);
          }
       }
   }
   return;
}

void ThermalPhoton::outputPhoton_SpvnpT(string path)
{
    ostringstream filename_stream_SpMatrix_eq;
    ostringstream filename_stream_SpMatrix_vis;
    ostringstream filename_stream_SpMatrix_bulkvis;
    ostringstream filename_stream_SpMatrix_vis_deltaf_restricted;
    ostringstream filename_stream_SpMatrix_bulkvis_deltaf_restricted;
    ostringstream filename_stream_SpMatrix_tot;
    ostringstream filename_stream_Spvn_eq;
    ostringstream filename_stream_Spvn_vis;
    ostringstream filename_stream_Spvn_bulkvis;
    ostringstream filename_stream_Spvn_vis_deltaf_restricted;
    ostringstream filename_stream_Spvn_bulkvis_deltaf_restricted;
    ostringstream filename_stream_Spvn_tot;
    ostringstream filename_stream_inte_Spvn_eq;
    ostringstream filename_stream_inte_Spvn_vis;
    ostringstream filename_stream_inte_Spvn_bulkvis;
    ostringstream filename_stream_inte_Spvn_vis_deltaf_restricted;
    ostringstream filename_stream_inte_Spvn_bulkvis_deltaf_restricted;
    ostringstream filename_stream_inte_Spvn_tot;

    filename_stream_SpMatrix_eq << path << emissionProcess_name << "_SpMatrix_eq.dat";
    filename_stream_SpMatrix_vis << path << emissionProcess_name << "_SpMatrix_vis.dat";
    filename_stream_SpMatrix_vis_deltaf_restricted << path << emissionProcess_name << "_SpMatrix_vis_deltaf_restricted.dat";
    filename_stream_SpMatrix_bulkvis << path << emissionProcess_name << "_SpMatrix_bulkvis.dat";
    filename_stream_SpMatrix_bulkvis_deltaf_restricted << path << emissionProcess_name << "_SpMatrix_bulkvis_deltaf_restricted.dat";
    filename_stream_SpMatrix_tot << path << emissionProcess_name << "_SpMatrix_tot.dat";
    filename_stream_Spvn_eq << path << emissionProcess_name << "_Spvn_eq.dat";
    filename_stream_Spvn_vis << path << emissionProcess_name << "_Spvn_vis.dat";
    filename_stream_Spvn_vis_deltaf_restricted << path << emissionProcess_name << "_Spvn_vis_deltaf_restricted.dat";
    filename_stream_Spvn_bulkvis << path << emissionProcess_name << "_Spvn_bulkvis.dat";
    filename_stream_Spvn_bulkvis_deltaf_restricted << path << emissionProcess_name << "_Spvn_bulkvis_deltaf_restricted.dat";
    filename_stream_Spvn_tot << path << emissionProcess_name << "_Spvn_tot.dat";
    filename_stream_inte_Spvn_eq << path << emissionProcess_name << "_Spvn_eq_inte.dat";
    filename_stream_inte_Spvn_vis << path << emissionProcess_name << "_Spvn_vis_inte.dat";
    filename_stream_inte_Spvn_vis_deltaf_restricted << path << emissionProcess_name << "_Spvn_vis_deltaf_restricted_inte.dat";
    filename_stream_inte_Spvn_bulkvis << path << emissionProcess_name << "_Spvn_bulkvis_inte.dat";
    filename_stream_inte_Spvn_bulkvis_deltaf_restricted << path << emissionProcess_name << "_Spvn_bulkvis_deltaf_restricted_inte.dat";
    filename_stream_inte_Spvn_tot << path << emissionProcess_name << "_Spvn_tot_inte.dat";

    ofstream fphotonSpMatrix_eq(filename_stream_SpMatrix_eq.str().c_str());
    ofstream fphotonSpMatrix_vis(filename_stream_SpMatrix_vis.str().c_str());
    ofstream fphotonSpMatrix_vis_deltaf_restricted(filename_stream_SpMatrix_vis_deltaf_restricted.str().c_str());
    ofstream fphotonSpMatrix_bulkvis(filename_stream_SpMatrix_bulkvis.str().c_str());
    ofstream fphotonSpMatrix_bulkvis_deltaf_restricted(filename_stream_SpMatrix_bulkvis_deltaf_restricted.str().c_str());
    ofstream fphotonSpMatrix_tot(filename_stream_SpMatrix_tot.str().c_str());
    ofstream fphotonSpvn_eq(filename_stream_Spvn_eq.str().c_str());
    ofstream fphotonSpvn_vis(filename_stream_Spvn_vis.str().c_str());
    ofstream fphotonSpvn_vis_deltaf_restricted(filename_stream_Spvn_vis_deltaf_restricted.str().c_str());
    ofstream fphotonSpvn_bulkvis(filename_stream_Spvn_bulkvis.str().c_str());
    ofstream fphotonSpvn_bulkvis_deltaf_restricted(filename_stream_Spvn_bulkvis_deltaf_restricted.str().c_str());
    ofstream fphotonSpvn_tot(filename_stream_Spvn_tot.str().c_str());
    ofstream fphotoninteSpvn_eq(filename_stream_inte_Spvn_eq.str().c_str());
    ofstream fphotoninteSpvn_vis(filename_stream_inte_Spvn_vis.str().c_str());
    ofstream fphotoninteSpvn_vis_deltaf_restricted(filename_stream_inte_Spvn_vis_deltaf_restricted.str().c_str());
    ofstream fphotoninteSpvn_bulkvis(filename_stream_inte_Spvn_bulkvis.str().c_str());
    ofstream fphotoninteSpvn_bulkvis_deltaf_restricted(filename_stream_inte_Spvn_bulkvis_deltaf_restricted.str().c_str());
    ofstream fphotoninteSpvn_tot(filename_stream_inte_Spvn_tot.str().c_str());

    for(int i=0;i<nphi;i++)
    {
      fphotonSpMatrix_eq << phi[i] << "  ";
      fphotonSpMatrix_vis << phi[i] << "  ";
      fphotonSpMatrix_bulkvis << phi[i] << "  ";
      fphotonSpMatrix_vis_deltaf_restricted << phi[i] << "  ";
      fphotonSpMatrix_bulkvis_deltaf_restricted << phi[i] << "  ";
      fphotonSpMatrix_tot << phi[i] << "  ";
      for(int j=0;j<np;j++)
        for(int k=0;k<nrapidity;k++)
        {
           fphotonSpMatrix_eq << scientific << setprecision(6) << setw(16) 
                              << dNd2pTdphidy_eq[j][i][k] << "  ";
           fphotonSpMatrix_vis << scientific << setprecision(6) << setw(16) 
                               << dNd2pTdphidy_vis[j][i][k] << "  ";
           fphotonSpMatrix_vis_deltaf_restricted << scientific << setprecision(6) << setw(16) 
                               << dNd2pTdphidy_vis_deltaf_restricted[j][i][k] << "  ";
           fphotonSpMatrix_bulkvis << scientific << setprecision(6) << setw(16) 
                                   << dNd2pTdphidy_bulkvis[j][i][k] << "  ";
           fphotonSpMatrix_bulkvis_deltaf_restricted << scientific << setprecision(6) << setw(16) 
                                   << dNd2pTdphidy_bulkvis_deltaf_restricted[j][i][k] << "  ";
           fphotonSpMatrix_tot << scientific << setprecision(6) << setw(16) 
                               << dNd2pTdphidy_tot[j][i][k] << "  ";
        }
      fphotonSpMatrix_eq << endl;
      fphotonSpMatrix_vis << endl;
      fphotonSpMatrix_vis_deltaf_restricted << endl;
      fphotonSpMatrix_bulkvis << endl;
      fphotonSpMatrix_bulkvis_deltaf_restricted << endl;
      fphotonSpMatrix_tot << endl;
    }
    for(int i=0;i<np;i++)
    {
      fphotonSpvn_eq << scientific << setprecision(6) << setw(16) 
                     << p[i] << "  " << dNd2pT_eq[i] << "  " ;
      fphotonSpvn_vis << scientific << setprecision(6) << setw(16) 
                      << p[i] << "  " << dNd2pT_vis[i] << "  " ;
      fphotonSpvn_vis_deltaf_restricted << scientific << setprecision(6) << setw(16) 
                      << p[i] << "  " << dNd2pT_vis_deltaf_restricted[i] << "  " ;
      fphotonSpvn_bulkvis << scientific << setprecision(6) << setw(16) 
                          << p[i] << "  " << dNd2pT_bulkvis[i] << "  " ;
      fphotonSpvn_bulkvis_deltaf_restricted << scientific << setprecision(6) << setw(16) 
                          << p[i] << "  " << dNd2pT_bulkvis_deltaf_restricted[i] << "  " ;
      fphotonSpvn_tot << scientific << setprecision(6) << setw(16) 
                      << p[i] << "  " << dNd2pT_tot[i] << "  " ;
      for(int order=1; order<norder; order++)
      {
         fphotonSpvn_eq << scientific << setprecision(6) << setw(16) 
                        << vnpT_cos_eq[order][i] << "  " 
                        << vnpT_sin_eq[order][i] << "  "
                        << sqrt(pow(vnpT_cos_eq[order][i], 2) + pow(vnpT_sin_eq[order][i], 2)) << "  ";
         fphotonSpvn_vis << scientific << setprecision(6) << setw(16) 
                         << vnpT_cos_vis[order][i] << "  "
                         << vnpT_sin_vis[order][i] << "  "
                         << sqrt(pow(vnpT_cos_vis[order][i], 2) + pow(vnpT_sin_vis[order][i], 2)) << "  ";
         fphotonSpvn_vis_deltaf_restricted << scientific << setprecision(6) << setw(16) 
                         << vnpT_cos_vis_deltaf_restricted[order][i] << "  "
                         << vnpT_sin_vis_deltaf_restricted[order][i] << "  "
                         << sqrt(pow(vnpT_cos_vis_deltaf_restricted[order][i], 2) + pow(vnpT_sin_vis_deltaf_restricted[order][i], 2)) << "  ";
         fphotonSpvn_bulkvis << scientific << setprecision(6) << setw(16) 
                             << vnpT_cos_bulkvis[order][i] << "  "
                             << vnpT_sin_bulkvis[order][i] << "  "
                             << sqrt(pow(vnpT_cos_bulkvis[order][i], 2) + pow(vnpT_sin_bulkvis[order][i], 2)) << "  ";
         fphotonSpvn_bulkvis_deltaf_restricted << scientific << setprecision(6) << setw(16) 
                             << vnpT_cos_bulkvis_deltaf_restricted[order][i] << "  "
                             << vnpT_sin_bulkvis_deltaf_restricted[order][i] << "  "
                             << sqrt(pow(vnpT_cos_bulkvis_deltaf_restricted[order][i], 2) + pow(vnpT_sin_bulkvis_deltaf_restricted[order][i], 2)) << "  ";
         fphotonSpvn_tot << scientific << setprecision(6) << setw(16) 
                         << vnpT_cos_tot[order][i] << "  "
                         << vnpT_sin_tot[order][i] << "  "
                         << sqrt(pow(vnpT_cos_tot[order][i], 2) + pow(vnpT_sin_tot[order][i], 2)) << "  ";
      }
      fphotonSpvn_eq << endl;
      fphotonSpvn_vis << endl;
      fphotonSpvn_vis_deltaf_restricted << endl;
      fphotonSpvn_bulkvis << endl;
      fphotonSpvn_bulkvis_deltaf_restricted << endl;
      fphotonSpvn_tot << endl;
    }

    for(int order = 0; order < norder; order++)
    {
        fphotoninteSpvn_eq << scientific << setprecision(6) << setw(16) 
                           << order << "   " << vn_cos_eq[order] << "   " 
                           << vn_sin_eq[order] << "   " 
                           << sqrt(pow(vn_cos_eq[order], 2) + pow(vn_sin_eq[order], 2)) 
                           << endl;
        fphotoninteSpvn_vis << scientific << setprecision(6) << setw(16) 
                            << order << "   " << vn_cos_vis[order] << "   " 
                            << vn_sin_vis[order] << "   " 
                            << sqrt(pow(vn_cos_vis[order], 2) + pow(vn_sin_vis[order], 2)) 
                            << endl;
        fphotoninteSpvn_vis_deltaf_restricted << scientific << setprecision(6) << setw(16) 
                            << order << "   " << vn_cos_vis_deltaf_restricted[order] << "   " 
                            << vn_sin_vis_deltaf_restricted[order] << "   " 
                            << sqrt(pow(vn_cos_vis_deltaf_restricted[order], 2) + pow(vn_sin_vis_deltaf_restricted[order], 2)) 
                            << endl;
        fphotoninteSpvn_bulkvis << scientific << setprecision(6) << setw(16) 
                                << order << "   " << vn_cos_bulkvis[order] << "   " 
                                << vn_sin_bulkvis[order] << "   " 
                                << sqrt(pow(vn_cos_bulkvis[order], 2) + pow(vn_sin_bulkvis[order], 2)) 
                                << endl;
        fphotoninteSpvn_bulkvis_deltaf_restricted << scientific << setprecision(6) << setw(16) 
                                << order << "   " << vn_cos_bulkvis_deltaf_restricted[order] << "   " 
                                << vn_sin_bulkvis_deltaf_restricted[order] << "   " 
                                << sqrt(pow(vn_cos_bulkvis_deltaf_restricted[order], 2) + pow(vn_sin_bulkvis_deltaf_restricted[order], 2)) 
                                << endl;
        fphotoninteSpvn_tot << scientific << setprecision(6) << setw(16) 
                            << order << "   " << vn_cos_tot[order] << "   " 
                            << vn_sin_tot[order] << "   " 
                            << sqrt(pow(vn_cos_tot[order], 2) + pow(vn_sin_tot[order], 2)) 
                            << endl;
    }

    fphotonSpMatrix_eq.close();
    fphotonSpMatrix_vis.close();
    fphotonSpMatrix_vis_deltaf_restricted.close();
    fphotonSpMatrix_bulkvis.close();
    fphotonSpMatrix_bulkvis_deltaf_restricted.close();
    fphotonSpMatrix_tot.close();
    fphotonSpvn_eq.close();
    fphotonSpvn_vis.close();
    fphotonSpvn_vis_deltaf_restricted.close();
    fphotonSpvn_bulkvis.close();
    fphotonSpvn_bulkvis_deltaf_restricted.close();
    fphotonSpvn_tot.close();
    fphotoninteSpvn_eq.close();
    fphotoninteSpvn_vis.close();
    fphotoninteSpvn_vis_deltaf_restricted.close();
    fphotoninteSpvn_bulkvis.close();
    fphotoninteSpvn_bulkvis_deltaf_restricted.close();
    fphotoninteSpvn_tot.close();

    return;
}

void ThermalPhoton::outputPhoton_SpvnpTdTdtau(string path)
{
    double dT = (Tcut_high - Tcut_low)/(nTcut - 1);
    double dtau = (Taucut_high - Taucut_low)/(nTaucut - 1);
    ostringstream filename_stream_dNdydTdtau_eq;
    ostringstream filename_stream_dNdydTdtau_vis;
    ostringstream filename_stream_dNdydTdtau_bulkvis;
    ostringstream filename_stream_dNdydTdtau_tot;

    filename_stream_dNdydTdtau_eq << path << emissionProcess_name << "_dNdydTdtau_eq.dat";
    filename_stream_dNdydTdtau_vis << path << emissionProcess_name << "_dNdydTdtau_vis.dat";
    filename_stream_dNdydTdtau_bulkvis << path << emissionProcess_name << "_dNdydTdtau_bulkvis.dat";
    filename_stream_dNdydTdtau_tot << path << emissionProcess_name << "_dNdydTdtau_tot.dat";

    ofstream fphotondNdy_eq(filename_stream_dNdydTdtau_eq.str().c_str());
    ofstream fphotondNdy_vis(filename_stream_dNdydTdtau_vis.str().c_str());
    ofstream fphotondNdy_bulkvis(filename_stream_dNdydTdtau_bulkvis.str().c_str());
    ofstream fphotondNdy_tot(filename_stream_dNdydTdtau_tot.str().c_str());

    for(int i = 0; i < nTcut; i++)
    {
       for(int j = 0; j < nTaucut; j++)
       {
          fphotondNdy_eq << dNdydTdtau_eq[i][j]/dT/dtau << "    ";
          fphotondNdy_vis << dNdydTdtau_vis[i][j]/dT/dtau << "    ";
          fphotondNdy_bulkvis << dNdydTdtau_bulkvis[i][j]/dT/dtau << "    ";
          fphotondNdy_tot << dNdydTdtau_tot[i][j]/dT/dtau << "    ";
       }
       fphotondNdy_eq << endl;
       fphotondNdy_vis << endl;
       fphotondNdy_bulkvis << endl;
       fphotondNdy_tot << endl;
    }
    fphotondNdy_eq.close();
    fphotondNdy_vis.close();
    fphotondNdy_bulkvis.close();
    fphotondNdy_tot.close();

    for(int order = 1; order < norder; order++)
    {
       ostringstream filename_stream_vncosdTdtau_eq;
       ostringstream filename_stream_vncosdTdtau_vis;
       ostringstream filename_stream_vncosdTdtau_bulkvis;
       ostringstream filename_stream_vncosdTdtau_tot;
       ostringstream filename_stream_vnsindTdtau_eq;
       ostringstream filename_stream_vnsindTdtau_vis;
       ostringstream filename_stream_vnsindTdtau_bulkvis;
       ostringstream filename_stream_vnsindTdtau_tot;
       filename_stream_vncosdTdtau_eq << path << emissionProcess_name << "_v_" << order << "_cos_dTdtau_eq.dat";
       filename_stream_vncosdTdtau_vis << path << emissionProcess_name << "_v_" << order << "_cos_dTdtau_vis.dat";
       filename_stream_vncosdTdtau_bulkvis << path << emissionProcess_name << "_v_" << order << "_cos_dTdtau_bulkvis.dat";
       filename_stream_vncosdTdtau_tot << path << emissionProcess_name << "_v_" << order << "_cos_dTdtau_tot.dat";
       filename_stream_vnsindTdtau_eq << path << emissionProcess_name << "_v_" << order << "_sin_dTdtau_eq.dat";
       filename_stream_vnsindTdtau_vis << path << emissionProcess_name << "_v_" << order << "_sin_dTdtau_vis.dat";
       filename_stream_vnsindTdtau_bulkvis << path << emissionProcess_name << "_v_" << order << "_sin_dTdtau_bulkvis.dat";
       filename_stream_vnsindTdtau_tot << path << emissionProcess_name << "_v_" << order << "_sin_dTdtau_tot.dat";

       ofstream fphotonvncos_eq(filename_stream_vncosdTdtau_eq.str().c_str());
       ofstream fphotonvncos_vis(filename_stream_vncosdTdtau_vis.str().c_str());
       ofstream fphotonvncos_bulkvis(filename_stream_vncosdTdtau_bulkvis.str().c_str());
       ofstream fphotonvncos_tot(filename_stream_vncosdTdtau_tot.str().c_str());
       ofstream fphotonvnsin_eq(filename_stream_vnsindTdtau_eq.str().c_str());
       ofstream fphotonvnsin_vis(filename_stream_vnsindTdtau_vis.str().c_str());
       ofstream fphotonvnsin_bulkvis(filename_stream_vnsindTdtau_bulkvis.str().c_str());
       ofstream fphotonvnsin_tot(filename_stream_vnsindTdtau_tot.str().c_str());
       for(int i = 0; i < nTcut; i++)
       {
          for(int j = 0; j < nTaucut; j++)
          {
             fphotonvncos_eq << vndTdtau_cos_eq[i][j][order]/dT/dtau << "    ";
             fphotonvncos_vis << vndTdtau_cos_vis[i][j][order]/dT/dtau << "    ";
             fphotonvncos_bulkvis << vndTdtau_cos_bulkvis[i][j][order]/dT/dtau << "    ";
             fphotonvncos_tot << vndTdtau_cos_tot[i][j][order]/dT/dtau << "    ";
             fphotonvnsin_eq << vndTdtau_sin_eq[i][j][order]/dT/dtau << "    ";
             fphotonvnsin_vis << vndTdtau_sin_vis[i][j][order]/dT/dtau << "    ";
             fphotonvnsin_bulkvis << vndTdtau_sin_bulkvis[i][j][order]/dT/dtau << "    ";
             fphotonvnsin_tot << vndTdtau_sin_tot[i][j][order]/dT/dtau << "    ";
          }
          fphotonvncos_eq << endl;
          fphotonvncos_vis << endl;
          fphotonvncos_bulkvis << endl;
          fphotonvncos_tot << endl;
          fphotonvnsin_eq << endl;
          fphotonvnsin_vis << endl;
          fphotonvnsin_bulkvis << endl;
          fphotonvnsin_tot << endl;
       }
       fphotonvncos_eq.close();
       fphotonvnsin_eq.close();
       fphotonvncos_vis.close();
       fphotonvnsin_vis.close();
       fphotonvncos_bulkvis.close();
       fphotonvnsin_bulkvis.close();
       fphotonvncos_tot.close();
       fphotonvnsin_tot.close();
    }

    return;
}

void ThermalPhoton::outputPhoton_SpvnpTdxperpdtau(string path)
{
    double dxperp = (xperp_high - xperp_low)/(n_xperp_cut - 1);
    double dtau = (tau_cut_high - tau_cut_low)/(n_tau_cut_xtau - 1);
    ostringstream filename_stream_dNdydxperpdtau_eq;
    ostringstream filename_stream_dNdydxperpdtau_vis;
    ostringstream filename_stream_dNdydxperpdtau_bulkvis;
    ostringstream filename_stream_dNdydxperpdtau_tot;

    filename_stream_dNdydxperpdtau_eq << path << emissionProcess_name << "_dNdydxperpdtau_eq.dat";
    filename_stream_dNdydxperpdtau_vis << path << emissionProcess_name << "_dNdydxperpdtau_vis.dat";
    filename_stream_dNdydxperpdtau_bulkvis << path << emissionProcess_name << "_dNdydxperpdtau_bulkvis.dat";
    filename_stream_dNdydxperpdtau_tot << path << emissionProcess_name << "_dNdydxperpdtau_tot.dat";

    ofstream fphotondNdy_eq(filename_stream_dNdydxperpdtau_eq.str().c_str());
    ofstream fphotondNdy_vis(filename_stream_dNdydxperpdtau_vis.str().c_str());
    ofstream fphotondNdy_bulkvis(filename_stream_dNdydxperpdtau_bulkvis.str().c_str());
    ofstream fphotondNdy_tot(filename_stream_dNdydxperpdtau_tot.str().c_str());

    for(int i = 0; i < n_xperp_cut; i++)
    {
       for(int j = 0; j < n_tau_cut_xtau; j++)
       {
          fphotondNdy_eq << dNdydxperpdtau_eq[i][j]/dxperp/dtau << "    ";
          fphotondNdy_vis << dNdydxperpdtau_vis[i][j]/dxperp/dtau << "    ";
          fphotondNdy_bulkvis << dNdydxperpdtau_bulkvis[i][j]/dxperp/dtau << "    ";
          fphotondNdy_tot << dNdydxperpdtau_tot[i][j]/dxperp/dtau << "    ";
       }
       fphotondNdy_eq << endl;
       fphotondNdy_vis << endl;
       fphotondNdy_bulkvis << endl;
       fphotondNdy_tot << endl;
    }
    fphotondNdy_eq.close();
    fphotondNdy_vis.close();
    fphotondNdy_bulkvis.close();
    fphotondNdy_tot.close();

    for(int order = 1; order < norder; order++)
    {
       ostringstream filename_stream_vncosdxperpdtau_eq;
       ostringstream filename_stream_vncosdxperpdtau_vis;
       ostringstream filename_stream_vncosdxperpdtau_bulkvis;
       ostringstream filename_stream_vncosdxperpdtau_tot;
       ostringstream filename_stream_vnsindxperpdtau_eq;
       ostringstream filename_stream_vnsindxperpdtau_vis;
       ostringstream filename_stream_vnsindxperpdtau_bulkvis;
       ostringstream filename_stream_vnsindxperpdtau_tot;
       filename_stream_vncosdxperpdtau_eq << path << emissionProcess_name << "_v_" << order << "_cos_dxperpdtau_eq.dat";
       filename_stream_vncosdxperpdtau_vis << path << emissionProcess_name << "_v_" << order << "_cos_dxperpdtau_vis.dat";
       filename_stream_vncosdxperpdtau_bulkvis << path << emissionProcess_name << "_v_" << order << "_cos_dxperpdtau_bulkvis.dat";
       filename_stream_vncosdxperpdtau_tot << path << emissionProcess_name << "_v_" << order << "_cos_dxperpdtau_tot.dat";
       filename_stream_vnsindxperpdtau_eq << path << emissionProcess_name << "_v_" << order << "_sin_dxperpdtau_eq.dat";
       filename_stream_vnsindxperpdtau_vis << path << emissionProcess_name << "_v_" << order << "_sin_dxperpdtau_vis.dat";
       filename_stream_vnsindxperpdtau_bulkvis << path << emissionProcess_name << "_v_" << order << "_sin_dxperpdtau_bulkvis.dat";
       filename_stream_vnsindxperpdtau_tot << path << emissionProcess_name << "_v_" << order << "_sin_dxperpdtau_tot.dat";

       ofstream fphotonvncos_eq(filename_stream_vncosdxperpdtau_eq.str().c_str());
       ofstream fphotonvncos_vis(filename_stream_vncosdxperpdtau_vis.str().c_str());
       ofstream fphotonvncos_bulkvis(filename_stream_vncosdxperpdtau_bulkvis.str().c_str());
       ofstream fphotonvncos_tot(filename_stream_vncosdxperpdtau_tot.str().c_str());
       ofstream fphotonvnsin_eq(filename_stream_vnsindxperpdtau_eq.str().c_str());
       ofstream fphotonvnsin_vis(filename_stream_vnsindxperpdtau_vis.str().c_str());
       ofstream fphotonvnsin_bulkvis(filename_stream_vnsindxperpdtau_bulkvis.str().c_str());
       ofstream fphotonvnsin_tot(filename_stream_vnsindxperpdtau_tot.str().c_str());
       for(int i = 0; i < n_xperp_cut; i++)
       {
          for(int j = 0; j < n_tau_cut_xtau; j++)
          {
             fphotonvncos_eq << vndxperpdtau_cos_eq[i][j][order]/dxperp/dtau << "    ";
             fphotonvncos_vis << vndxperpdtau_cos_vis[i][j][order]/dxperp/dtau << "    ";
             fphotonvncos_bulkvis << vndxperpdtau_cos_bulkvis[i][j][order]/dxperp/dtau << "    ";
             fphotonvncos_tot << vndxperpdtau_cos_tot[i][j][order]/dxperp/dtau << "    ";
             fphotonvnsin_eq << vndxperpdtau_sin_eq[i][j][order]/dxperp/dtau << "    ";
             fphotonvnsin_vis << vndxperpdtau_sin_vis[i][j][order]/dxperp/dtau << "    ";
             fphotonvnsin_bulkvis << vndxperpdtau_sin_bulkvis[i][j][order]/dxperp/dtau << "    ";
             fphotonvnsin_tot << vndxperpdtau_sin_tot[i][j][order]/dxperp/dtau << "    ";
          }
          fphotonvncos_eq << endl;
          fphotonvncos_vis << endl;
          fphotonvncos_bulkvis << endl;
          fphotonvncos_tot << endl;
          fphotonvnsin_eq << endl;
          fphotonvnsin_vis << endl;
          fphotonvnsin_bulkvis << endl;
          fphotonvnsin_tot << endl;
       }
       fphotonvncos_eq.close();
       fphotonvnsin_eq.close();
       fphotonvncos_vis.close();
       fphotonvnsin_vis.close();
       fphotonvncos_bulkvis.close();
       fphotonvnsin_bulkvis.close();
       fphotonvncos_tot.close();
       fphotonvnsin_tot.close();
    }

    return;
}

void ThermalPhoton::interpolation2D_bilinear(double varX, double* varY, int Y_length, double** Table2D_ptr, double* results)
//this function is used the most frequent one, it needs to be as fast as possible
{
     double varX_min = EmissionrateTb_Xmin;
     double varY_min = EmissionrateTb_Ymin;
     double dX = EmissionrateTb_dX;
     double dY = EmissionrateTb_dY;
     //int tb_sizeX = EmissionrateTb_sizeX;
     int tb_sizeY = EmissionrateTb_sizeY;

     double dXdYinverse = 1/(dX*dY);
     double dYinverse = 1/dY;

     int idx_Y;
     double f00, f10, f01, f11;
     double ddy;
     double dYhigh;
     
     int idx_X;
     double ddx;
     double dXhigh;
     idx_X = (int)((varX-varX_min)/dX);
     ddx = (varX - (varX_min + idx_X*dX));
     dXhigh = (dX - ddx)*dXdYinverse;
     ddx = ddx*dXdYinverse;


     for(int ii=0; ii<Y_length; ii++)
     {
        idx_Y = (int)((varY[ii]-varY_min)*dYinverse);

     //if it is out the table, doing linear extrapolation
    /* if (idx_X<0 || idx_X>=tb_sizeX-1)
     {
        if(Iwarning == 1)
        {
           cout<<"interpolation2D_bilinear: varX out of bounds!"<<endl;
           cout<<"varX= " << varX<<"idx_X= "<<idx_X<<endl;
        }
        if(idx_X < 0) idx_X = 0;
        if(idx_X >= tb_sizeX-1) idx_X = tb_sizeX - 1;
     }*/
     
        /*if(Iwarning == 1)
        {
           cout<<"interpolation2D_bilinear: varY out of bounds!"<<endl;
           cout<<"varY= " << varY<<"idx_Y= "<<idx_Y<<endl;
        }*/
        if(idx_Y < 0) idx_Y = 0;
        if(idx_Y >= tb_sizeY-1) idx_Y = tb_sizeY - 2;
       
        ddy = varY[ii] - EmissionrateTb_Yidxptr[idx_Y];
        dYhigh = dY - ddy;

        f00 = Table2D_ptr[idx_X][idx_Y];
        f10 = Table2D_ptr[idx_X+1][idx_Y];
        f01 = Table2D_ptr[idx_X][idx_Y+1];
        f11 = Table2D_ptr[idx_X+1][idx_Y+1];
       
        results[ii] = ((f00*dYhigh + f01*ddy)*dXhigh + (f10*dYhigh + f11*ddy)*ddx);
     }

     return;
}

void ThermalPhoton::update_rates_with_polyakov_suppression()
{
     for(int i=0; i<EmissionrateTb_sizeX; i++)
     {
         double T_local = EmissionrateTb_Xmin + i*EmissionrateTb_dX;
         double suppression_factor = get_polyakov_suppression_factor(T_local);
         for(int j=0; j<EmissionrateTb_sizeY; j++)
         {
             Emission_eqrateTb_ptr[i][j] += log(suppression_factor);
         }
    }
}

double ThermalPhoton::get_polyakov_suppression_factor(double T_in_GeV)
{
    double T_in_MeV=T_in_GeV*1e3;
    const double a = 1.49201e-9;
    const double b = -7.48088e-7;
    const double c = -0.000480142;
    const double d = 0.420208;
    double Qratio=a*pow(T_in_MeV,3) + b*T_in_MeV*T_in_MeV + c*T_in_MeV + d;
    double f_suppression=10./3.*Qratio*Qratio-4*Qratio+1;
    return(f_suppression);
}
