#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include "main.h"
#include "readindata.h"
#include "emissionfunction.h"
#include "Stopwatch.h"
#include "arsenal.h"
#include "ParameterReader.h"

#define AMOUNT_OF_OUTPUT 0 // smaller value means less outputs

using namespace std;


// Class EmissionFunctionArray ------------------------------------------
EmissionFunctionArray::EmissionFunctionArray(ParameterReader* paraRdr_in, double particle_y_in, Table* chosen_particles_in, Table* pT_tab_in, Table* phi_tab_in, Table* eta_tab_in, particle_info* particles_in, int Nparticles_in, FO_surf* FOsurf_ptr_in, long FO_length_in)
{
  paraRdr = paraRdr_in;
  particle_y = particle_y_in;
  pT_tab = pT_tab_in; pT_tab_length = pT_tab->getNumberOfRows();
  phi_tab = phi_tab_in; phi_tab_length = phi_tab->getNumberOfRows();
  eta_tab = eta_tab_in; eta_tab_length = eta_tab->getNumberOfRows();

  dN_ptdptdphidy = new Table(pT_tab_length, phi_tab_length);
  dN_ptdptdphidy_filename = "results/dN_ptdptdphidy.dat";

  // get control parameters
  CALCULATEDED3P = paraRdr->getVal("calculate_dEd3p");
  INCLUDE_BULKDELTAF = paraRdr->getVal("turn_on_bulk");
  INCLUDE_MUB = paraRdr->getVal("turn_on_muB");
  INCLUDE_DELTAF = paraRdr->getVal("turn_on_shear");
  GROUPING_PARTICLES = paraRdr->getVal("grouping_particles");
  PARTICLE_DIFF_TOLERANCE = paraRdr->getVal("particle_diff_tolerance");
  USE_HISTORIC_FORMAT = paraRdr->getVal("use_historic_format");
  F0_IS_NOT_SMALL = paraRdr->getVal("f0_is_not_small");
  bulk_deltaf_kind = paraRdr->getVal("bulk_deltaf_kind");

  if(CALCULATEDED3P == 1)
  {
     dE_ptdptdphidy = new Table(pT_tab_length, phi_tab_length);
     dE_ptdptdphidy_filename = "results/dE_ptdptdphidy.dat";
  }

  particles = particles_in;
  Nparticles = Nparticles_in;

  FOsurf_ptr = FOsurf_ptr_in;
  FO_length = FO_length_in;

  number_of_chosen_particles = chosen_particles_in->getNumberOfRows();

  chosen_particles_01_table = new int[Nparticles];
  for (int n=0; n<Nparticles; n++) chosen_particles_01_table[n]=0;
  for (int m=0; m<number_of_chosen_particles; m++)
  {
    int monval = chosen_particles_in->get(1,m+1);
    for (int n=0; n<Nparticles; n++)
    {
      if (particles[n].monval==monval)
      {
        chosen_particles_01_table[n]=1;
        break;
      }
    }
  }
  // next, for sampling processes
  chosen_particles_sampling_table = new int[number_of_chosen_particles];
  // first copy the chosen_particles table, but now using indices instead of monval
  int current_idx = 0;
  for (int m=0; m<number_of_chosen_particles; m++)
  {
    int monval = chosen_particles_in->get(1,m+1);
    for (int n=0; n<Nparticles; n++)
    {
      if (particles[n].monval==monval)
      {
        chosen_particles_sampling_table[current_idx] = n;
        current_idx ++;
        break;
      }
    }
  }
  // next re-order them so that particles with similar mass are adjacent
  if (GROUPING_PARTICLES == 1) // sort particles according to their mass; bubble-sorting
  {
    for (int m=0; m<number_of_chosen_particles; m++)
      for (int n=0; n<number_of_chosen_particles-m-1; n++)
        if (particles[chosen_particles_sampling_table[n]].mass > particles[chosen_particles_sampling_table[n+1]].mass)
        {
          // swap them
          int particle_idx = chosen_particles_sampling_table[n+1];
          chosen_particles_sampling_table[n+1] = chosen_particles_sampling_table[n];
          chosen_particles_sampling_table[n] = particle_idx;
        }
  }

  flow_differential_filename_old =  "results/v2data.dat";
  flow_integrated_filename_old = "results/v2data-inte.dat";
  flow_differential_filename = "results/thermal_%d_vndata.dat";
  flow_integrated_filename = "results/thermal_%d_integrated_vndata.dat";
  if(CALCULATEDED3P == 1)
  {
     energyflow_differential_filename_old =  "results/ET_v2data.dat";
     energyflow_integrated_filename_old = "results/ET_v2data-inte.dat";
     energyflow_differential_filename = "results/thermal_%d_ET_vndata.dat";
     energyflow_integrated_filename = "results/thermal_%d_ET_integrated_vndata.dat";
  }
  last_particle_idx = -1;

  //arrays for bulk delta f coefficients
  bulkdf_coeff = new Table ("EOS/BulkDf_Coefficients_Hadrons_s95p-v0-PCE.dat");
}


EmissionFunctionArray::~EmissionFunctionArray()
{
  delete dN_ptdptdphidy;
  if(CALCULATEDED3P == 1) delete dE_ptdptdphidy;
  delete[] chosen_particles_01_table;
  delete[] chosen_particles_sampling_table;
  delete bulkdf_coeff;
}


void EmissionFunctionArray::calculate_dN_ptdptdphidy(int particle_idx)
// Calculate dN_xxx array.
{
  last_particle_idx = particle_idx;
  Stopwatch sw;
  sw.tic();

  double y = particle_y;
  particle_info* particle;
  particle = &particles[particle_idx];

  double mass = particle->mass;
  double sign = particle->sign;
  double degen = particle->gspin;
  int baryon = particle->baryon;

  double prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;

  FO_surf* surf = &FOsurf_ptr[0];

  double *bulkvisCoefficients;
  if(bulk_deltaf_kind == 0)
      bulkvisCoefficients = new double [3];
  else
      bulkvisCoefficients = new double [2];


  // for intermediate results
  double dN_ptdptdphidy_tab[pT_tab_length][phi_tab_length];
  double dE_ptdptdphidy_tab[pT_tab_length][phi_tab_length];
  for (int i=0; i<pT_tab_length; i++)
  for (int j=0; j<phi_tab_length; j++)
  {
    dN_ptdptdphidy_tab[i][j] = 0.0;
    dE_ptdptdphidy_tab[i][j] = 0.0;
  }

  // pre-calculated variables //!!!!!!
  double trig_phi_table[phi_tab_length][2]; // 2: 0,1-> cos,sin
  for (int j=0; j<phi_tab_length; j++)
  {
    double phi = phi_tab->get(1,j+1);
    trig_phi_table[j][0] = cos(phi);
    trig_phi_table[j][1] = sin(phi);
  }

  double hypertrig_etas_table[eta_tab_length][2]; // 2: 0,1-> cosh,sinh
  double delta_eta_tab[eta_tab_length]; // cache it
  for (int k=0; k<eta_tab_length; k++)
  {
    double eta_s = eta_tab->get(1,k+1);
    hypertrig_etas_table[k][0] = cosh(y-eta_s);
    hypertrig_etas_table[k][1] = sinh(y-eta_s);
    delta_eta_tab[k] = eta_tab->get(2,k+1);
  }

  //---------------------------
  // THE main summation loop
  //---------------------------
  double progress_total = pT_tab_length*phi_tab_length;
  if (AMOUNT_OF_OUTPUT>0) print_progressbar(-1);

  for (int i=0; i<pT_tab_length; i++)
  {
      double pT = pT_tab->get(1,i+1); 
      double mT = sqrt(mass*mass + pT*pT);

      for (int j=0; j<phi_tab_length; j++)
      {
          double px = pT*trig_phi_table[j][0];
          double py = pT*trig_phi_table[j][1];

          double dN_ptdptdphidy_tmp = 0.0;
          double dE_ptdptdphidy_tmp = 0.0;

          for (long l=0; l<FO_length; l++)
          {
              surf = &FOsurf_ptr[l];

              double Tdec = surf->Tdec;
              double Pdec = surf->Pdec;
              double Edec = surf->Edec;
              double mu = surf->particle_mu[last_particle_idx];
              double tau = surf->tau;
              double gammaT = surf->u0;
              double ux = surf->u1;
              double uy = surf->u2;
              double da0 = surf->da0;
              double da1 = surf->da1;
              double da2 = surf->da2;
              double pi00 = surf->pi00;
              double pi01 = surf->pi01;
              double pi02 = surf->pi02;
              double pi11 = surf->pi11;
              double pi12 = surf->pi12;
              double pi22 = surf->pi22;
              double pi33 = surf->pi33;
              double muB = surf->muB;
              double bulkPi = 0.0;
              double deltaf_prefactor = 0.0;
              if(INCLUDE_DELTAF)
                  deltaf_prefactor = 1.0/(2.0*Tdec*Tdec*(Edec+Pdec));
              if(INCLUDE_BULKDELTAF == 1)
              {
                  if(bulk_deltaf_kind == 0)
                      bulkPi = surf->bulkPi;
                  else
                      bulkPi = surf->bulkPi/hbarC;   // need unit in fm^-4 for the parameterization
                  getbulkvisCoefficients(Tdec, bulkvisCoefficients);
              }

              for (int k=0; k<eta_tab_length; k++)
              {
                  double delta_eta = delta_eta_tab[k];

                  double pt = mT*hypertrig_etas_table[k][0];
                  double pz = mT*hypertrig_etas_table[k][1];

                  double pdotu = pt*gammaT - px*ux - py*uy;
                  double expon = (pdotu - mu - baryon*muB) / Tdec;
                  double f0 = 1./(exp(expon)+sign);       //thermal equilibrium distributions

                  // Must adjust this to be correct for the p*del \tau term.  The plus sign is
                  // due to the fact that the DA# variables are for the covariant surface integration
                  double pdsigma = pt*da0 + px*da1 + py*da2;

                  //viscous corrections
                  double delta_f_shear = 0.0;
                  if(INCLUDE_DELTAF)
                  {
                      double Wfactor = pt*pt*pi00 - 2.0*pt*px*pi01 - 2.0*pt*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
                      delta_f_shear = (1 - F0_IS_NOT_SMALL*sign*f0)*Wfactor*deltaf_prefactor;
                  }
                  double delta_f_bulk = 0.0;
                  if (INCLUDE_BULKDELTAF == 1)
                  {
                      if(bulk_deltaf_kind == 0)
                          delta_f_bulk = -(1. - F0_IS_NOT_SMALL*sign*f0)*bulkPi*(bulkvisCoefficients[0]*mass*mass + bulkvisCoefficients[1]*pdotu + bulkvisCoefficients[2]*pdotu*pdotu);
                      else if (bulk_deltaf_kind == 1)
                      {

                          double E_over_T = pdotu/Tdec;
                          double mass_over_T = mass/Tdec;
                          delta_f_bulk = -1.0*(1.-sign*f0)/E_over_T*bulkvisCoefficients[0]*(mass_over_T*mass_over_T/3. - bulkvisCoefficients[1]*E_over_T*E_over_T)*bulkPi;
                      }
                      else if (bulk_deltaf_kind == 2)
                      {
                          double E_over_T = pdotu/Tdec;
                          delta_f_bulk = -1.*(1.-sign*f0)*(-bulkvisCoefficients[0] + bulkvisCoefficients[1]*E_over_T)*bulkPi;
                      }
                      else if (bulk_deltaf_kind == 3)
                      {
                          double E_over_T = pdotu/Tdec;
                          delta_f_bulk = -1.0*(1.-sign*f0)/sqrt(E_over_T)*(-bulkvisCoefficients[0] + bulkvisCoefficients[1]*E_over_T)*bulkPi;
                      }
                      else if (bulk_deltaf_kind == 4)
                      {
                          double E_over_T = pdotu/Tdec;
                          delta_f_bulk = -1.0*(1.-sign*f0)*(bulkvisCoefficients[0] - bulkvisCoefficients[1]/E_over_T)*bulkPi;
                      }
                  }

                  double result;
                  //if(1. + delta_f_shear + delta_f_bulk < 0.0) //set results to zero when delta f turns whole expression to negative
                  //   result = 0.0;
                  //else
                     result = prefactor*degen*f0*(1. + delta_f_shear + delta_f_bulk)*pdsigma*tau;

                  dN_ptdptdphidy_tmp += result*delta_eta;
                  if(CALCULATEDED3P == 1) dE_ptdptdphidy_tmp += result*delta_eta*mT;
              } // k
          } // l

          dN_ptdptdphidy_tab[i][j] = dN_ptdptdphidy_tmp;
          if (CALCULATEDED3P == 1) dE_ptdptdphidy_tab[i][j] = dE_ptdptdphidy_tmp;
          if (AMOUNT_OF_OUTPUT>0) print_progressbar((i*phi_tab_length+j)/progress_total);
      }
  //cout << int(100.0*(i+1)/pT_tab_length) << "% completed" << endl;
  }
  if (AMOUNT_OF_OUTPUT>0) print_progressbar(1);

  for (int i=0; i<pT_tab_length; i++)
  for (int j=0; j<phi_tab_length; j++)
  {
    dN_ptdptdphidy->set(i+1,j+1,dN_ptdptdphidy_tab[i][j]);
    if(CALCULATEDED3P == 1) dE_ptdptdphidy->set(i+1,j+1,dE_ptdptdphidy_tab[i][j]);
  }

  delete [] bulkvisCoefficients;

  sw.toc();
  cout << endl << "Finished " << sw.takeTime() << " seconds." << endl;
}

void EmissionFunctionArray::write_dN_ptdptdphidy_toFile()
// Append the dN_xxx results to file.
{
  ofstream of1(dN_ptdptdphidy_filename.c_str(), ios_base::app);
  dN_ptdptdphidy->printTable(of1);
  of1.close();

  if(CALCULATEDED3P == 1)
  {
     ofstream of2(dE_ptdptdphidy_filename.c_str(), ios_base::app);
     dE_ptdptdphidy->printTable(of2);
     of2.close();
  }
}


//***************************************************************************
void EmissionFunctionArray::calculate_flows(int to_order, string flow_differential_filename_in, string flow_integrated_filename_in)
// Calculate flow from order from_order to to_order and store them to files.
{
  /*
  cout << endl
       <<"*************************************************"
       << endl
       << "Function calculate_flows started... " << endl;*/
  Stopwatch sw;
  sw.tic();

  int from_order = 1;

  int number_of_flows = to_order-from_order+1;

  Table vn_diff(3+number_of_flows*3, pT_tab_length); // line format: pT, mT, dN/(pT dpT), flow_1_real, flow_1_imag, flow_1_norm, ...
  Table vn_inte(6, to_order+1); // line format: order# (starting from 0), numerator_real, numerator_imag, flow_real, flow_imag, flow_norm

  double mass = particles[last_particle_idx].mass;

  //---------------------
  // differential flow
  //---------------------
  //cout << "Calculating differential flows... ";

  double normalization[pT_tab_length]; // normalization factor
  for (int i=0; i<pT_tab_length; i++) normalization[i] = 0.0;

  double vn[pT_tab_length][number_of_flows][2]; // diff_flow numerators; 2: 0,1->real,imag
  for (int i=0; i<pT_tab_length; i++)
  for (int t=0; t<number_of_flows; t++)
    {vn[i][t][0]=0; vn[i][t][1]=0;}

  for (int i=0; i<pT_tab_length; i++)
  //for (int i=0; i<1; i++) // for debugging
  {
    double pT = pT_tab->get(1,i+1); // pT_weight = pT_tab->get(2,i+1);
    double mT = sqrt(mass*mass + pT*pT);

    // phi integration
    for(int j=0; j<phi_tab_length; j++)
    //for(int j=0; j<1; j++) // for debugging
    {
      double phi = phi_tab->get(1,j+1), phi_weight = phi_tab->get(2,j+1);
      double dN = dN_ptdptdphidy->get(i+1,j+1);

      normalization[i] += dN*phi_weight;
      for (int order=from_order; order<=to_order; order++)
      {
        vn[i][order-from_order][0] += dN*phi_weight*cos(order*phi);
        vn[i][order-from_order][1] += dN*phi_weight*sin(order*phi);
      }
    }

    normalization[i] = normalization[i] + 1e-30;
    // store values
    vn_diff.set(1, i+1, pT);
    vn_diff.set(2, i+1, mT-mass);
    vn_diff.set(3, i+1, normalization[i]/(2.0*M_PI)); // 2*pi: azimuthal angle averaged
    for (int t=0; t<number_of_flows; t++)
    {
      vn_diff.set(4+t*3, i+1, vn[i][t][0]/normalization[i]);
      vn_diff.set(5+t*3, i+1, vn[i][t][1]/normalization[i]);
      vn_diff.set(6+t*3, i+1, sqrt(vn[i][t][0]*vn[i][t][0]+vn[i][t][1]*vn[i][t][1])/normalization[i]);
    }

  }
  //cout << "done." << endl;


  //---------------------
  // integrated flow
  //---------------------
  //cout << "Calculating integrated flows... ";

  double normalizationi = 0;

  double vni[number_of_flows][2]; // integrated_flow numerators; 2: 0,1->real,imag
  for (int t=0; t<number_of_flows; t++) {vni[t][0]=0; vni[t][1]=0;}

  for (int i=0; i<pT_tab_length; i++)
  //for (int i=0; i<1; i++) // for debugging
  {
    double pT = pT_tab->get(1,i+1), pT_weight = pT_tab->get(2,i+1);

    normalizationi += normalization[i]*pT*pT_weight;

    for (int order=from_order; order<=to_order; order++)
    {
      vni[order-from_order][0] += vn[i][order-from_order][0]*pT*pT_weight;
      vni[order-from_order][1] += vn[i][order-from_order][1]*pT*pT_weight;
    }

  }


  // store values
  // To mimic:
  // WRITE(30,941) " ", N, " ",X(N)," ",Y(N),
  //   &  " ",X(N)/X(0)," ",Y(N)/X(0)," ",
  //   &  sqrt(X(N)*X(N)+Y(N)*Y(N))/X(0)
  vn_inte.set(1, 1, 0);
  vn_inte.set(2, 1, normalizationi);
  vn_inte.set(3, 1, 0);
  vn_inte.set(4, 1, 1);
  vn_inte.set(5, 1, 0);
  vn_inte.set(6, 1, 1);

  for (int t=0; t<number_of_flows; t++)
  {
    vn_inte.set(1, t+2, from_order+t);
    vn_inte.set(2, t+2, vni[from_order+t-1][0]);
    vn_inte.set(3, t+2, vni[from_order+t-1][1]);
    vn_inte.set(4, t+2, vni[from_order+t-1][0]/normalizationi);
    vn_inte.set(5, t+2, vni[from_order+t-1][1]/normalizationi);
    vn_inte.set(6, t+2, sqrt(vni[from_order+t-1][0]*vni[from_order+t-1][0]+vni[from_order+t-1][1]*vni[from_order+t-1][1])/normalizationi);
  }
  //cout << "done." << endl;

  // save to files
  //cout << "Writing to files... ";
  ofstream of1(flow_differential_filename_in.c_str(), ios_base::app);
  vn_diff.printTable(of1);
  of1.close();

  ofstream of2(flow_integrated_filename_in.c_str(), ios_base::app);
  vn_inte.printTable(of2);
  of2.close();
  //cout << "done." << endl;

  sw.toc();
  //cout << "calculate_flows finishes " << sw.takeTime() << " seconds." << endl;

}
//***************************************************************************


//***************************************************************************
void EmissionFunctionArray::calculate_Energyflows(int to_order, string flow_differential_filename_in, string flow_integrated_filename_in)
// Calculate flow from order from_order to to_order and store them to files.
{
  Stopwatch sw;
  sw.tic();

  int from_order = 1;

  int number_of_flows = to_order-from_order+1;

  Table vn_diff(3+number_of_flows*3, pT_tab_length); // line format: pT, mT, dN/(pT dpT), flow_1_real, flow_1_imag, flow_1_norm, ...
  Table vn_inte(6, to_order+1); // line format: order# (starting from 0), numerator_real, numerator_imag, flow_real, flow_imag, flow_norm

  double mass = particles[last_particle_idx].mass;

  //---------------------
  // differential flow
  //---------------------
  //cout << "Calculating differential flows... ";

  double normalization[pT_tab_length]; // normalization factor
  for (int i=0; i<pT_tab_length; i++) normalization[i] = 0.0;

  double vn[pT_tab_length][number_of_flows][2]; // diff_flow numerators; 2: 0,1->real,imag
  for (int i=0; i<pT_tab_length; i++)
  for (int t=0; t<number_of_flows; t++)
    {vn[i][t][0]=0; vn[i][t][1]=0;}

  for (int i=0; i<pT_tab_length; i++)
  //for (int i=0; i<1; i++) // for debugging
  {
    double pT = pT_tab->get(1,i+1); // pT_weight = pT_tab->get(2,i+1);
    double mT = sqrt(mass*mass + pT*pT);

    // phi integration
    for(int j=0; j<phi_tab_length; j++)
    //for(int j=0; j<1; j++) // for debugging
    {
      double phi = phi_tab->get(1,j+1), phi_weight = phi_tab->get(2,j+1);
      double dE = dE_ptdptdphidy->get(i+1,j+1);

      normalization[i] += dE*phi_weight;
      for (int order=from_order; order<=to_order; order++)
      {
        vn[i][order-from_order][0] += dE*phi_weight*cos(order*phi);
        vn[i][order-from_order][1] += dE*phi_weight*sin(order*phi);
      }
    }

    normalization[i] = normalization[i] + 1e-30;
    // store values
    vn_diff.set(1, i+1, pT);
    vn_diff.set(2, i+1, mT-mass);
    vn_diff.set(3, i+1, normalization[i]/(2.0*M_PI)); // 2*pi: azimuthal angle averaged
    for (int t=0; t<number_of_flows; t++)
    {
      vn_diff.set(4+t*3, i+1, vn[i][t][0]/normalization[i]);
      vn_diff.set(5+t*3, i+1, vn[i][t][1]/normalization[i]);
      vn_diff.set(6+t*3, i+1, sqrt(vn[i][t][0]*vn[i][t][0]+vn[i][t][1]*vn[i][t][1])/normalization[i]);
    }

  }
  //cout << "done." << endl;


  //---------------------
  // integrated flow
  //---------------------
  //cout << "Calculating integrated flows... ";

  double normalizationi = 0;

  double vni[number_of_flows][2]; // integrated_flow numerators; 2: 0,1->real,imag
  for (int t=0; t<number_of_flows; t++) {vni[t][0]=0; vni[t][1]=0;}

  for (int i=0; i<pT_tab_length; i++)
  //for (int i=0; i<1; i++) // for debugging
  {
    double pT = pT_tab->get(1,i+1), pT_weight = pT_tab->get(2,i+1);

    normalizationi += normalization[i]*pT*pT_weight;

    for (int order=from_order; order<=to_order; order++)
    {
      vni[order-from_order][0] += vn[i][order-from_order][0]*pT*pT_weight;
      vni[order-from_order][1] += vn[i][order-from_order][1]*pT*pT_weight;
    }

  }

  vn_inte.set(1, 1, 0);
  vn_inte.set(2, 1, normalizationi);
  vn_inte.set(3, 1, 0);
  vn_inte.set(4, 1, 1);
  vn_inte.set(5, 1, 0);
  vn_inte.set(6, 1, 1);

  for (int t=0; t<number_of_flows; t++)
  {
    vn_inte.set(1, t+2, from_order+t);
    vn_inte.set(2, t+2, vni[from_order+t-1][0]);
    vn_inte.set(3, t+2, vni[from_order+t-1][1]);
    vn_inte.set(4, t+2, vni[from_order+t-1][0]/normalizationi);
    vn_inte.set(5, t+2, vni[from_order+t-1][1]/normalizationi);
    vn_inte.set(6, t+2, sqrt(vni[from_order+t-1][0]*vni[from_order+t-1][0]+vni[from_order+t-1][1]*vni[from_order+t-1][1])/normalizationi);
  }

  // save to files
  ofstream of1(flow_differential_filename_in.c_str(), ios_base::app);
  vn_diff.printTable(of1);
  of1.close();

  ofstream of2(flow_integrated_filename_in.c_str(), ios_base::app);
  vn_inte.printTable(of2);
  of2.close();

  sw.toc();
  //cout << "calculate_flows finishes " << sw.takeTime() << " seconds." << endl;

}
//***************************************************************************




//*********************************************************************************************
void EmissionFunctionArray::calculate_dN_ptdptdphidy_and_flows_4all(int to_order)
// Calculate dNArrays and flows for all particles given in chosen_particle file.
{
    if (USE_HISTORIC_FORMAT == 1)
    {
          cout << endl
               <<"**********************************************************"
               << endl
               << "Function calculate_dN_ptdptdphidy_and_flows_4all(old) started... " << endl;
          Stopwatch sw;
          sw.tic();

          remove(dN_ptdptdphidy_filename.c_str());
          remove(flow_differential_filename_old.c_str());
          remove(flow_integrated_filename_old.c_str());
          if(CALCULATEDED3P == 1)
          {
             remove(dE_ptdptdphidy_filename.c_str());
             remove(energyflow_differential_filename_old.c_str());
             remove(energyflow_integrated_filename_old.c_str());
          }

          particle_info* particle = NULL;

          for (int n=0; n<Nparticles; n++)
          {
            particle = &particles[n];
            cout << "Index: " << n << ", Name: " << particle->name << ", Monte-carlo index: " << particle->monval;

            // first, dN_xxx arrays:
            if (chosen_particles_01_table[n]==0)
            {
              cout << " -- Skipped." << endl;
              dN_ptdptdphidy->setAll(0.0);
              if(CALCULATEDED3P == 1) dE_ptdptdphidy->setAll(0.0);
              last_particle_idx = n; // fake a "calculation"
            }
            else
            {
              cout << " -- Processing... " << endl;
              calculate_dN_ptdptdphidy(n);
            }
            write_dN_ptdptdphidy_toFile();

            // next flows:
            ofstream of1(flow_differential_filename_old.c_str(), ios_base::app);
            of1 << "# Output for particle: " << particle->name << endl;
            of1 << "#                 " << particle->monval << endl;
            of1.close();

            ofstream of2(flow_integrated_filename_old.c_str(), ios_base::app);
            of2 << "# For: " << particle->name << endl;
            of2.close();
            calculate_flows(to_order, flow_differential_filename_old, flow_integrated_filename_old);
            if(CALCULATEDED3P == 1)
            {
               ofstream of1(energyflow_differential_filename_old.c_str(), ios_base::app);
               of1 << "# Output for particle: " << particle->name << endl;
               of1 << "#                 " << particle->monval << endl;
               of1.close();
               ofstream of2(energyflow_integrated_filename_old.c_str(), ios_base::app);
               of2 << "# For: " << particle->name << endl;
               of2.close();
               calculate_Energyflows(to_order, energyflow_differential_filename_old, energyflow_integrated_filename_old);
            }
          }


          sw.toc();
          cout << "calculate_dN_ptdptdphidy_and_flows_4all finishes " << sw.takeTime() << " seconds." << endl;
    }
    else
    {
            cout << endl
            << "****************************************************************"
            << endl
            << "Function calculate_dN_ptdptdphidy_and_flows_4all(new) started... " << endl;
            Stopwatch sw;
            sw.tic();

            // prepare a huge array to store calculated dN_ptdptdphidy
            Table* dNs[Nparticles];
            for (int n=0; n<Nparticles; n++) dNs[n]=NULL;
            
            Table* dEs[Nparticles];
            for (int n=0; n<Nparticles; n++) dEs[n]=NULL;

            // loop over chosen particles
            particle_info* particle = NULL;
            for (int m=0; m<number_of_chosen_particles; m++)
            {
                int particle_idx = chosen_particles_sampling_table[m];
                particle = &particles[particle_idx];
                int monval = particle->monval;
                cout << "Index: " << m << ", Name: " << particle->name << ", Monte-carlo index: " << monval << endl;
                // Calculate dN / (ptdpt dphi dy)
                if (m>0 && particles_are_the_same(particle_idx, chosen_particles_sampling_table[m-1]))
                {
                   cout << " -- Using dN_ptdptdphidy from previous calculation... " << endl;
                   last_particle_idx = particle_idx; // fake a calculation
                }
                else
                {
                    cout << " -- Calculating dN_ptdptdphidy... " << endl;
                    calculate_dN_ptdptdphidy(particle_idx);
                }

                // Store calculated table
                dNs[particle_idx] = new Table(*dN_ptdptdphidy);
                if(CALCULATEDED3P == 1) dEs[particle_idx] = new Table(*dE_ptdptdphidy);

                char buffer_diff[500], buffer_inte[500];
                sprintf(buffer_diff, flow_differential_filename.c_str(), monval);
                remove(buffer_diff);
                sprintf(buffer_inte, flow_integrated_filename.c_str(), monval);
                remove(buffer_inte);
                calculate_flows(to_order, buffer_diff, buffer_inte);
                
                if(CALCULATEDED3P == 1)
                {
                   sprintf(buffer_diff, energyflow_differential_filename.c_str(), monval);
                   remove(buffer_diff);
                   sprintf(buffer_inte, energyflow_integrated_filename.c_str(), monval);
                   remove(buffer_inte);
                   calculate_Energyflows(to_order, buffer_diff, buffer_inte);
                }
            }

            // write out dN / (ptdpt dphi dy) matrices
            remove(dN_ptdptdphidy_filename.c_str());
            ofstream of(dN_ptdptdphidy_filename.c_str(), ios_base::app);
            Table zero(dN_ptdptdphidy->getNumberOfCols(), dN_ptdptdphidy->getNumberOfRows(), 0);
            for (int n=0; n<Nparticles; n++)
            {
                if (dNs[n]==NULL)
                {
                    zero.printTable(of);
                }
                else
                {
                    dNs[n]->printTable(of);
                    delete dNs[n];
                }
            }
            of.close();

            if(CALCULATEDED3P == 1)
            {
               remove(dE_ptdptdphidy_filename.c_str());
               ofstream of_E(dE_ptdptdphidy_filename.c_str(), ios_base::app);
               Table zero(dE_ptdptdphidy->getNumberOfCols(), dE_ptdptdphidy->getNumberOfRows(), 0);
               for (int n=0; n<Nparticles; n++)
               {
                   if (dEs[n]==NULL)
                   {
                       zero.printTable(of_E);
                   }
                   else
                   {
                       dEs[n]->printTable(of_E);
                       delete dEs[n];
                   }
               }
               of_E.close();
            }

            sw.toc();
            cout << " -- Calculate_dN_ptdptdphidy_and_flows_4all finishes " << sw.takeTime() << " seconds." << endl;

    }

}





//**********************************************************************************************
bool EmissionFunctionArray::particles_are_the_same(int idx1, int idx2)
{
    if (particles[idx1].sign!=particles[idx2].sign) return false;
    if (particles[idx1].baryon!=particles[idx2].baryon && INCLUDE_MUB) return false;
    if (abs(particles[idx1].mass-particles[idx2].mass) / (particles[idx2].mass+1e-30) > PARTICLE_DIFF_TOLERANCE) return false;
    for (long l=0; l<FO_length; l++)
    {
        double chem1 = FOsurf_ptr[l].particle_mu[idx1], chem2 = FOsurf_ptr[l].particle_mu[idx2];
        if (abs(chem1-chem2)/(chem2+1e-30) > PARTICLE_DIFF_TOLERANCE)
        {
          return false;
        }

    }

    return true;
}

void EmissionFunctionArray::getbulkvisCoefficients(double Tdec, double* bulkvisCoefficients)
{
   double Tdec_fm = Tdec/hbarC;  // [1/fm]
   double Tdec_fm_power[11];    // cache the polynomial power of Tdec_fm
   Tdec_fm_power[1] = Tdec_fm;
   for(int ipower = 2; ipower < 11; ipower++)
       Tdec_fm_power[ipower] = Tdec_fm_power[ipower-1]*Tdec_fm;
   if(bulk_deltaf_kind == 0)       // 14 moment expansion
   {
        // load from file
        bulkvisCoefficients[0] = bulkdf_coeff->interp(1, 2, Tdec_fm, 5)/pow(hbarC, 3);  //B0 [fm^3/GeV^3]
        bulkvisCoefficients[1] = bulkdf_coeff->interp(1, 3, Tdec_fm, 5)/pow(hbarC, 2);  // D0 [fm^3/GeV^2]
        bulkvisCoefficients[2] = bulkdf_coeff->interp(1, 4, Tdec_fm, 5)/pow(hbarC, 3);  // E0 [fm^3/GeV^3]
        // parameterization for mu = 0
        //bulkvisCoefficients[0] = exp(-15.04512474*Tdec_fm + 11.76194266)/pow(hbarC, 3); //B0[fm^3/GeV^3]
        //bulkvisCoefficients[1] = exp( -12.45699277*Tdec_fm + 11.4949293)/hbarC/hbarC;  // D0 [fm^3/GeV^2]
        //bulkvisCoefficients[2] = -exp(-14.45087586*Tdec_fm + 11.62716548)/pow(hbarC, 3);  // E0 [fm^3/GeV^3]
   }
   else if(bulk_deltaf_kind == 1)  // relaxation type
   {
       // parameterization from JF
       // A Polynomial fit to each coefficient -- X is the temperature in fm^-1
       // Both fits are reliable between T=100 -- 180 MeV , do not trust it beyond
       bulkvisCoefficients[0] = (  642096.624265727 
                                 - 8163329.49562861*Tdec_fm_power[1] 
                                 + 47162768.4292073*Tdec_fm_power[2] 
                                 - 162590040.002683*Tdec_fm_power[3] 
                                 + 369637951.096896*Tdec_fm_power[4] 
                                 - 578181331.809836*Tdec_fm_power[5] 
                                 + 629434830.225675*Tdec_fm_power[6] 
                                 - 470493661.096657*Tdec_fm_power[7] 
                                 + 230936465.421*Tdec_fm_power[8] 
                                 - 67175218.4629078*Tdec_fm_power[9] 
                                 + 8789472.32652964*Tdec_fm_power[10]);

       bulkvisCoefficients[1] = (  1.18171174036192 
                                 - 17.6740645873717*Tdec_fm_power[1]
                                 + 136.298469057177*Tdec_fm_power[2] 
                                 - 635.999435106846*Tdec_fm_power[3] 
                                 + 1918.77100633321*Tdec_fm_power[4] 
                                 - 3836.32258307711*Tdec_fm_power[5] 
                                 + 5136.35746882372*Tdec_fm_power[6] 
                                 - 4566.22991441914*Tdec_fm_power[7] 
                                 + 2593.45375240886*Tdec_fm_power[8] 
                                 - 853.908199724349*Tdec_fm_power[9]
                                 + 124.260460450113*Tdec_fm_power[10]);
   }
   else if (bulk_deltaf_kind == 2)
   {
       // A Polynomial fit to each coefficient -- Tfm is the temperature in fm^-1
       // Both fits are reliable between T=100 -- 180 MeV , do not trust it beyond
       bulkvisCoefficients[0] = (  
               21091365.1182649 - 290482229.281782*Tdec_fm_power[1] 
             + 1800423055.01882*Tdec_fm_power[2] - 6608608560.99887*Tdec_fm_power[3] 
             + 15900800422.7138*Tdec_fm_power[4] - 26194517161.8205*Tdec_fm_power[5] 
             + 29912485360.2916*Tdec_fm_power[6] - 23375101221.2855*Tdec_fm_power[7] 
             + 11960898238.0134*Tdec_fm_power[8] - 3618358144.18576*Tdec_fm_power[9] 
             + 491369134.205902*Tdec_fm_power[10]);

       bulkvisCoefficients[1] = (  
               4007863.29316896 - 55199395.3534188*Tdec_fm_power[1] 
             + 342115196.396492*Tdec_fm_power[2] - 1255681487.77798*Tdec_fm_power[3] 
             + 3021026280.08401*Tdec_fm_power[4] - 4976331606.85766*Tdec_fm_power[5] 
             + 5682163732.74188*Tdec_fm_power[6] - 4439937810.57449*Tdec_fm_power[7] 
             + 2271692965.05568*Tdec_fm_power[8] - 687164038.128814*Tdec_fm_power[9] 
             + 93308348.3137008*Tdec_fm_power[10]);
   }
   else if (bulk_deltaf_kind == 3)
   {
       bulkvisCoefficients[0] = (
               160421664.93603 - 2212807124.97991*Tdec_fm_power[1] 
             + 13707913981.1425*Tdec_fm_power[2] - 50204536518.1767*Tdec_fm_power[3] 
             + 120354649094.362*Tdec_fm_power[4] - 197298426823.223*Tdec_fm_power[5] 
             + 223953760788.288*Tdec_fm_power[6] - 173790947240.829*Tdec_fm_power[7] 
             + 88231322888.0423*Tdec_fm_power[8] - 26461154892.6963*Tdec_fm_power[9] 
             + 3559805050.19592*Tdec_fm_power[10]);
       bulkvisCoefficients[1] = (
               33369186.2536556 - 460293490.420478*Tdec_fm_power[1] 
             + 2851449676.09981*Tdec_fm_power[2] - 10443297927.601*Tdec_fm_power[3] 
             + 25035517099.7809*Tdec_fm_power[4] - 41040777943.4963*Tdec_fm_power[5] 
             + 46585225878.8723*Tdec_fm_power[6] - 36150531001.3718*Tdec_fm_power[7] 
             + 18353035766.9323*Tdec_fm_power[8] - 5504165325.05431*Tdec_fm_power[9] 
             + 740468257.784873*Tdec_fm_power[10]);
   }
   else if (bulk_deltaf_kind == 4)
   {
       bulkvisCoefficients[0] = (  
               1167272041.90731 - 16378866444.6842*Tdec_fm_power[1] 
             + 103037615761.617*Tdec_fm_power[2] - 382670727905.111*Tdec_fm_power[3] 
             + 929111866739.436*Tdec_fm_power[4] - 1540948583116.54*Tdec_fm_power[5] 
             + 1767975890298.1*Tdec_fm_power[6] - 1385606389545*Tdec_fm_power[7] 
             + 709922576963.213*Tdec_fm_power[8] - 214726945096.326*Tdec_fm_power[9] 
             + 29116298091.9219*Tdec_fm_power[10]);
       bulkvisCoefficients[1] = (
               5103633637.7213 - 71612903872.8163*Tdec_fm_power[1] 
             + 450509014334.964*Tdec_fm_power[2] - 1673143669281.46*Tdec_fm_power[3] 
             + 4062340452589.89*Tdec_fm_power[4] - 6737468792456.4*Tdec_fm_power[5] 
             + 7730102407679.65*Tdec_fm_power[6] - 6058276038129.83*Tdec_fm_power[7] 
             + 3103990764357.81*Tdec_fm_power[8] - 938850005883.612*Tdec_fm_power[9] 
             + 127305171097.249*Tdec_fm_power[10]);
   }
   return;
}
