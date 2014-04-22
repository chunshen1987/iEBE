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
    dNd2pT_bulkvis = new double [np];
    dNd2pT_tot = new double [np];
    dNd2pTdphidy_eq = new double** [np];
    dNd2pTdphidy_vis = new double** [np];
    dNd2pTdphidy_bulkvis = new double** [np];
    dNd2pTdphidy_tot = new double** [np];
    for(int i=0;i<np;i++)
    {
      dNd2pT_eq[i] = 0.0;
      dNd2pT_vis[i] = 0.0;
      dNd2pT_bulkvis[i] = 0.0;
      dNd2pT_tot[i] = 0.0;
      dNd2pTdphidy_eq[i] = new double* [nphi];
      dNd2pTdphidy_vis[i] = new double* [nphi];
      dNd2pTdphidy_bulkvis[i] = new double* [nphi];
      dNd2pTdphidy_tot[i] = new double* [nphi];
      for(int j=0;j<nphi;j++)
      {
         dNd2pTdphidy_eq[i][j] = new double [nrapidity];
         dNd2pTdphidy_vis[i][j] = new double [nrapidity];
         dNd2pTdphidy_bulkvis[i][j] = new double [nrapidity];
         dNd2pTdphidy_tot[i][j] = new double [nrapidity];
         for(int k=0;k<nrapidity;k++)
         {
            dNd2pTdphidy_eq[i][j][k] = 0.0;
            dNd2pTdphidy_vis[i][j][k] = 0.0;
            dNd2pTdphidy_bulkvis[i][j][k] = 0.0;
            dNd2pTdphidy_tot[i][j][k] = 0.0;
         }
      }
    }
    vnpT_cos_eq = new double* [norder];
    vnpT_sin_eq = new double* [norder];
    vnpT_cos_vis = new double* [norder];
    vnpT_sin_vis = new double* [norder];
    vnpT_cos_bulkvis = new double* [norder];
    vnpT_sin_bulkvis = new double* [norder];
    vnpT_cos_tot = new double* [norder];
    vnpT_sin_tot = new double* [norder];
    for(int order=0; order<norder; order++)
    {
      vnpT_cos_eq[order] = new double [np];
      vnpT_sin_eq[order] = new double [np];
      vnpT_cos_vis[order] = new double [np];
      vnpT_sin_vis[order] = new double [np];
      vnpT_cos_bulkvis[order] = new double [np];
      vnpT_sin_bulkvis[order] = new double [np];
      vnpT_cos_tot[order] = new double [np];
      vnpT_sin_tot[order] = new double [np];
      for(int i =0; i < np; i++)
      {
         vnpT_cos_eq[order][i] = 0.0;
         vnpT_cos_vis[order][i] = 0.0;
         vnpT_cos_bulkvis[order][i] = 0.0;
         vnpT_cos_tot[order][i] = 0.0;
         vnpT_sin_eq[order][i] = 0.0;
         vnpT_sin_vis[order][i] = 0.0;
         vnpT_sin_bulkvis[order][i] = 0.0;
         vnpT_sin_tot[order][i] = 0.0;
      }
    }

    dNdy_eq = 0.0;
    dNdy_vis = 0.0;
    dNdy_bulkvis = 0.0;
    dNdy_tot = 0.0;
    vn_cos_eq = new double [norder];
    vn_sin_eq = new double [norder];
    vn_cos_vis = new double [norder];
    vn_sin_vis = new double [norder];
    vn_cos_bulkvis = new double [norder];
    vn_sin_bulkvis = new double [norder];
    vn_cos_tot = new double [norder];
    vn_sin_tot = new double [norder];
    for(int i = 0; i < norder; i++)
    {
       vn_cos_eq[i] = 0.0;
       vn_sin_eq[i] = 0.0;
       vn_cos_vis[i] = 0.0;
       vn_sin_vis[i] = 0.0;
       vn_cos_bulkvis[i] = 0.0;
       vn_sin_bulkvis[i] = 0.0;
       vn_cos_tot[i] = 0.0;
       vn_sin_tot[i] = 0.0;
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
    delete [] dNd2pT_tot;

    for(int i = 0; i < np; i++)
    {
       for(int j = 0; j < nphi; j++)
       {
          delete [] dNd2pTdphidy_eq[i][j];
          delete [] dNd2pTdphidy_vis[i][j];
          delete [] dNd2pTdphidy_bulkvis[i][j];
          delete [] dNd2pTdphidy_tot[i][j];
       }
       delete [] dNd2pTdphidy_eq[i];
       delete [] dNd2pTdphidy_vis[i];
       delete [] dNd2pTdphidy_bulkvis[i];
       delete [] dNd2pTdphidy_tot[i];
    }
    delete [] dNd2pTdphidy_eq;
    delete [] dNd2pTdphidy_vis;
    delete [] dNd2pTdphidy_bulkvis;
    delete [] dNd2pTdphidy_tot;

    for(int i = 0; i < norder; i++)
    {
       delete [] vnpT_cos_eq[i];
       delete [] vnpT_sin_eq[i];
       delete [] vnpT_cos_vis[i];
       delete [] vnpT_sin_vis[i];
       delete [] vnpT_cos_bulkvis[i];
       delete [] vnpT_sin_bulkvis[i];
       delete [] vnpT_cos_tot[i];
       delete [] vnpT_sin_tot[i];
    }
    delete [] vnpT_cos_eq;
    delete [] vnpT_sin_eq;
    delete [] vnpT_cos_vis;
    delete [] vnpT_sin_vis;
    delete [] vnpT_cos_bulkvis;
    delete [] vnpT_sin_bulkvis;
    delete [] vnpT_cos_tot;
    delete [] vnpT_sin_tot;

    delete [] vn_cos_eq;
    delete [] vn_sin_eq;
    delete [] vn_cos_vis;
    delete [] vn_sin_vis;
    delete [] vn_cos_bulkvis;
    delete [] vn_sin_bulkvis;
    delete [] vn_cos_tot;
    delete [] vn_sin_tot;
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
           Emission_eqrateTb_ptr[i][j] = log(Photonemission_eqrateTable_ptr->getTbdata(i,j));  
           Emission_viscous_rateTb_ptr[i][j] = Photonemission_viscous_rateTable_ptr->getTbdata(i,j)/Photonemission_eqrateTable_ptr->getTbdata(i,j);
           Emission_bulkvis_rateTb_ptr[i][j] = Photonemission_bulkvis_rateTable_ptr->getTbdata(i,j)/Photonemission_eqrateTable_ptr->getTbdata(i,j);
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
          dNd2pTdphidy_eq[l][m][k] += temp_eq_sum;
          dNd2pTdphidy_vis[l][m][k] += temp_eq_sum + temp_vis_sum;
          dNd2pTdphidy_bulkvis[l][m][k] += temp_eq_sum + temp_bulkvis_sum;
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


void ThermalPhoton::calPhoton_SpvnpT()
//calculate the photon spectra and differential vn at mid-rapidity
{
   int k = 0;  //calculate at y = 0
   for(int i=0;i<np;i++)
   {
       for(int j=0;j<nphi;j++)
       {
         dNd2pT_eq[i] += dNd2pTdphidy_eq[i][j][k]*phi_weight[j];  //dN/(dy pT dpT)
         dNd2pT_vis[i] += dNd2pTdphidy_vis[i][j][k]*phi_weight[j];
         dNd2pT_bulkvis[i] += dNd2pTdphidy_bulkvis[i][j][k]*phi_weight[j];
         dNd2pT_tot[i] += dNd2pTdphidy_tot[i][j][k]*phi_weight[j];
         for(int order=0; order<norder; order++)
         {
            vnpT_cos_eq[order][i] += dNd2pTdphidy_eq[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_cos_vis[order][i] += dNd2pTdphidy_vis[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_cos_bulkvis[order][i] += dNd2pTdphidy_bulkvis[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_cos_tot[order][i] += dNd2pTdphidy_tot[i][j][k]*cos(order*phi[j])*phi_weight[j];
            vnpT_sin_eq[order][i] += dNd2pTdphidy_eq[i][j][k]*sin(order*phi[j])*phi_weight[j];
            vnpT_sin_vis[order][i] += dNd2pTdphidy_vis[i][j][k]*sin(order*phi[j])*phi_weight[j];
            vnpT_sin_bulkvis[order][i] += dNd2pTdphidy_bulkvis[i][j][k]*sin(order*phi[j])*phi_weight[j];
            vnpT_sin_tot[order][i] += dNd2pTdphidy_tot[i][j][k]*sin(order*phi[j])*phi_weight[j];
         }
       }
       dNdy_eq += dNd2pT_eq[i]*p[i]*p_weight[i];   //dN/dy
       dNdy_vis += dNd2pT_vis[i]*p[i]*p_weight[i];
       dNdy_bulkvis += dNd2pT_bulkvis[i]*p[i]*p_weight[i];
       dNdy_tot += dNd2pT_tot[i]*p[i]*p_weight[i];
       for(int order = 0; order < norder ; order++)
       {
          vn_cos_eq[order] += vnpT_cos_eq[order][i]*p[i]*p_weight[i];
          vn_sin_eq[order] += vnpT_sin_eq[order][i]*p[i]*p_weight[i];
          vn_cos_vis[order] += vnpT_cos_vis[order][i]*p[i]*p_weight[i];
          vn_sin_vis[order] += vnpT_sin_vis[order][i]*p[i]*p_weight[i];
          vn_cos_bulkvis[order] += vnpT_cos_bulkvis[order][i]*p[i]*p_weight[i];
          vn_sin_bulkvis[order] += vnpT_sin_bulkvis[order][i]*p[i]*p_weight[i];
          vn_cos_tot[order] += vnpT_cos_tot[order][i]*p[i]*p_weight[i];
          vn_sin_tot[order] += vnpT_sin_tot[order][i]*p[i]*p_weight[i];
          
          //vn(pT)
          vnpT_cos_eq[order][i] = vnpT_cos_eq[order][i]/dNd2pT_eq[i];
          vnpT_cos_vis[order][i] = vnpT_cos_vis[order][i]/dNd2pT_vis[i];
          vnpT_cos_bulkvis[order][i] = vnpT_cos_bulkvis[order][i]/dNd2pT_bulkvis[i];
          vnpT_cos_tot[order][i] = vnpT_cos_tot[order][i]/dNd2pT_tot[i];
          vnpT_sin_eq[order][i] = vnpT_sin_eq[order][i]/dNd2pT_eq[i];
          vnpT_sin_vis[order][i] = vnpT_sin_vis[order][i]/dNd2pT_vis[i];
          vnpT_sin_bulkvis[order][i] = vnpT_sin_bulkvis[order][i]/dNd2pT_bulkvis[i];
          vnpT_sin_tot[order][i] = vnpT_sin_tot[order][i]/dNd2pT_tot[i];
       }
       dNd2pT_eq[i] = dNd2pT_eq[i]/(2*M_PI);  //dN/(2pi dy pT dpT)
       dNd2pT_vis[i] = dNd2pT_vis[i]/(2*M_PI);
       dNd2pT_bulkvis[i] = dNd2pT_bulkvis[i]/(2*M_PI);
       dNd2pT_tot[i] = dNd2pT_tot[i]/(2*M_PI);
   }
   for(int order = 1; order < norder ; order++)
   {
       //vn
       vn_cos_eq[order] = vn_cos_eq[order]/dNdy_eq;
       vn_sin_eq[order] = vn_sin_eq[order]/dNdy_eq;
       vn_cos_vis[order] = vn_cos_vis[order]/dNdy_vis;
       vn_sin_vis[order] = vn_sin_vis[order]/dNdy_vis;
       vn_cos_bulkvis[order] = vn_cos_bulkvis[order]/dNdy_vis;
       vn_sin_bulkvis[order] = vn_sin_bulkvis[order]/dNdy_vis;
       vn_cos_tot[order] = vn_cos_tot[order]/dNdy_tot;
       vn_sin_tot[order] = vn_sin_tot[order]/dNdy_tot;
   }
   return;
}

void ThermalPhoton::outputPhoton_SpvnpT(string path)
{
    ostringstream filename_stream_SpMatrix_eq;
    ostringstream filename_stream_SpMatrix_vis;
    ostringstream filename_stream_SpMatrix_bulkvis;
    ostringstream filename_stream_SpMatrix_tot;
    ostringstream filename_stream_Spvn_eq;
    ostringstream filename_stream_Spvn_vis;
    ostringstream filename_stream_Spvn_bulkvis;
    ostringstream filename_stream_Spvn_tot;
    ostringstream filename_stream_inte_Spvn_eq;
    ostringstream filename_stream_inte_Spvn_vis;
    ostringstream filename_stream_inte_Spvn_bulkvis;
    ostringstream filename_stream_inte_Spvn_tot;

    filename_stream_SpMatrix_eq << path << emissionProcess_name << "_SpMatrix_eq.dat";
    filename_stream_SpMatrix_vis << path << emissionProcess_name << "_SpMatrix_vis.dat";
    filename_stream_SpMatrix_bulkvis << path << emissionProcess_name << "_SpMatrix_bulkvis.dat";
    filename_stream_SpMatrix_tot << path << emissionProcess_name << "_SpMatrix_tot.dat";
    filename_stream_Spvn_eq << path << emissionProcess_name << "_Spvn_eq.dat";
    filename_stream_Spvn_vis << path << emissionProcess_name << "_Spvn_vis.dat";
    filename_stream_Spvn_bulkvis << path << emissionProcess_name << "_Spvn_bulkvis.dat";
    filename_stream_Spvn_tot << path << emissionProcess_name << "_Spvn_tot.dat";
    filename_stream_inte_Spvn_eq << path << emissionProcess_name << "_Spvn_eq_inte.dat";
    filename_stream_inte_Spvn_vis << path << emissionProcess_name << "_Spvn_vis_inte.dat";
    filename_stream_inte_Spvn_bulkvis << path << emissionProcess_name << "_Spvn_bulkvis_inte.dat";
    filename_stream_inte_Spvn_tot << path << emissionProcess_name << "_Spvn_tot_inte.dat";

    ofstream fphotonSpMatrix_eq(filename_stream_SpMatrix_eq.str().c_str());
    ofstream fphotonSpMatrix_vis(filename_stream_SpMatrix_vis.str().c_str());
    ofstream fphotonSpMatrix_bulkvis(filename_stream_SpMatrix_bulkvis.str().c_str());
    ofstream fphotonSpMatrix_tot(filename_stream_SpMatrix_tot.str().c_str());
    ofstream fphotonSpvn_eq(filename_stream_Spvn_eq.str().c_str());
    ofstream fphotonSpvn_vis(filename_stream_Spvn_vis.str().c_str());
    ofstream fphotonSpvn_bulkvis(filename_stream_Spvn_bulkvis.str().c_str());
    ofstream fphotonSpvn_tot(filename_stream_Spvn_tot.str().c_str());
    ofstream fphotoninteSpvn_eq(filename_stream_inte_Spvn_eq.str().c_str());
    ofstream fphotoninteSpvn_vis(filename_stream_inte_Spvn_vis.str().c_str());
    ofstream fphotoninteSpvn_bulkvis(filename_stream_inte_Spvn_bulkvis.str().c_str());
    ofstream fphotoninteSpvn_tot(filename_stream_inte_Spvn_tot.str().c_str());

    for(int i=0;i<nphi;i++)
    {
      fphotonSpMatrix_eq << phi[i] << "  ";
      fphotonSpMatrix_vis << phi[i] << "  ";
      fphotonSpMatrix_bulkvis << phi[i] << "  ";
      fphotonSpMatrix_tot << phi[i] << "  ";
      for(int j=0;j<np;j++)
        for(int k=0;k<nrapidity;k++)
        {
           fphotonSpMatrix_eq << scientific << setprecision(6) << setw(16) 
                              << dNd2pTdphidy_eq[j][i][k] << "  ";
           fphotonSpMatrix_vis << scientific << setprecision(6) << setw(16) 
                               << dNd2pTdphidy_vis[j][i][k] << "  ";
           fphotonSpMatrix_bulkvis << scientific << setprecision(6) << setw(16) 
                                   << dNd2pTdphidy_bulkvis[j][i][k] << "  ";
           fphotonSpMatrix_tot << scientific << setprecision(6) << setw(16) 
                               << dNd2pTdphidy_tot[j][i][k] << "  ";
        }
      fphotonSpMatrix_eq << endl;
      fphotonSpMatrix_vis << endl;
      fphotonSpMatrix_bulkvis << endl;
      fphotonSpMatrix_tot << endl;
    }
    for(int i=0;i<np;i++)
    {
      fphotonSpvn_eq << scientific << setprecision(6) << setw(16) 
                     << p[i] << "  " << dNd2pT_eq[i] << "  " ;
      fphotonSpvn_vis << scientific << setprecision(6) << setw(16) 
                      << p[i] << "  " << dNd2pT_vis[i] << "  " ;
      fphotonSpvn_bulkvis << scientific << setprecision(6) << setw(16) 
                          << p[i] << "  " << dNd2pT_bulkvis[i] << "  " ;
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
         fphotonSpvn_bulkvis << scientific << setprecision(6) << setw(16) 
                             << vnpT_cos_bulkvis[order][i] << "  "
                             << vnpT_sin_bulkvis[order][i] << "  "
                             << sqrt(pow(vnpT_cos_bulkvis[order][i], 2) + pow(vnpT_sin_bulkvis[order][i], 2)) << "  ";
         fphotonSpvn_tot << scientific << setprecision(6) << setw(16) 
                         << vnpT_cos_tot[order][i] << "  "
                         << vnpT_sin_tot[order][i] << "  "
                         << sqrt(pow(vnpT_cos_tot[order][i], 2) + pow(vnpT_sin_tot[order][i], 2)) << "  ";
      }
      fphotonSpvn_eq << endl;
      fphotonSpvn_vis << endl;
      fphotonSpvn_bulkvis << endl;
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
        fphotoninteSpvn_bulkvis << scientific << setprecision(6) << setw(16) 
                                << order << "   " << vn_cos_bulkvis[order] << "   " 
                                << vn_sin_bulkvis[order] << "   " 
                                << sqrt(pow(vn_cos_bulkvis[order], 2) + pow(vn_sin_bulkvis[order], 2)) 
                                << endl;
        fphotoninteSpvn_tot << scientific << setprecision(6) << setw(16) 
                            << order << "   " << vn_cos_tot[order] << "   " 
                            << vn_sin_tot[order] << "   " 
                            << sqrt(pow(vn_cos_tot[order], 2) + pow(vn_sin_tot[order], 2)) 
                            << endl;
    }

    fphotonSpMatrix_eq.close();
    fphotonSpMatrix_vis.close();
    fphotonSpMatrix_bulkvis.close();
    fphotonSpMatrix_tot.close();
    fphotonSpvn_eq.close();
    fphotonSpvn_vis.close();
    fphotonSpvn_bulkvis.close();
    fphotonSpvn_tot.close();
    fphotoninteSpvn_eq.close();
    fphotoninteSpvn_vis.close();
    fphotoninteSpvn_bulkvis.close();
    fphotoninteSpvn_tot.close();

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

