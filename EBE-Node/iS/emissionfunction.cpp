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

#define AMOUNT_OF_OUTPUT 0 // smaller value means less outputs
#define F0_IS_NOT_SMALL 0 // set to 0 to agree with Azspectra; set to 1 for reality
#define USE_HISTORIC_FORMAT 1 // 0: use new way of outputting
#define GROUPING_PARTICLES 1 // set to 1 to perform calculations for similar particles together
#define PARTICLE_DIFF_TOLERANCE 0.01 // particles with mass and chemical potential (for each FZ-cell) difference less than this value will be considered to be identical (b/c Cooper-Frye)
#define INCLUDE_DELTAF 1 // include delta f correction to particle distribution function in Cooper-Frye Formula

using namespace std;


// Class EmissionFunctionArray ------------------------------------------
EmissionFunctionArray::EmissionFunctionArray(double particle_y_in, Table* chosen_particles_in, Table* pT_tab_in, Table* phi_tab_in, Table* eta_tab_in, particle_info* particles_in, int Nparticles_in, FO_surf* FOsurf_ptr_in, long FO_length_in)
{
  particle_y = particle_y_in;
  pT_tab = pT_tab_in; pT_tab_length = pT_tab->getNumberOfRows();
  phi_tab = phi_tab_in; phi_tab_length = phi_tab->getNumberOfRows();
  eta_tab = eta_tab_in; eta_tab_length = eta_tab->getNumberOfRows();

  dN_ptdptdphidy = new Table(pT_tab_length, phi_tab_length);
  dN_ptdptdphidy_filename = "results/dN_ptdptdphidy.dat";

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
  // first copy the chosen_particles table, but now using indecies instead of monval
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
  if (GROUPING_PARTICLES) // sort particles according to their mass; bubble-sorting
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
  last_particle_idx = -1;

}


EmissionFunctionArray::~EmissionFunctionArray()
{
  delete dN_ptdptdphidy;
  delete[] chosen_particles_01_table;
  delete[] chosen_particles_sampling_table;
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

  double prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;

  FO_surf* surf = &FOsurf_ptr[0];

  // for intermedia results
  double dN_ptdptdphidy_tab[pT_tab_length][phi_tab_length];
  for (int i=0; i<pT_tab_length; i++)
  for (int j=0; j<phi_tab_length; j++)
    dN_ptdptdphidy_tab[i][j] = 0.0;

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
  //for (int i=0; i<1; i++) // for debug
  {
      double pT = pT_tab->get(1,i+1); // pT_weight = pT_tab->get(2,i+1); // unused
      double mT = sqrt(mass*mass + pT*pT);

      for (int j=0; j<phi_tab_length; j++)
      //for (int j=0; j<1; j++) // for debug
      {
          // double phi = phi_tab->get(1,j+1), phi_weight = phi_tab->get(2,j+1); // unused

          // old way
          //double cos_phi = cos(phi);
          //double sin_phi = sin(phi);
          //double px = pT*cos_phi;
          //double py = pT*sin_phi;
          // new way
          double px = pT*trig_phi_table[j][0];
          double py = pT*trig_phi_table[j][1];

          double dN_ptdptdphidy_tmp = 0.0;

          for (long l=0; l<FO_length; l++)
          {
              surf = &FOsurf_ptr[l];

              double Tdec = surf->Tdec;
              double Pdec = surf->Pdec;
              double Edec = surf->Edec;
              double deltaf_prefactor = 1.0/(2.0*Tdec*Tdec*(Edec+Pdec))*INCLUDE_DELTAF;

              double mu = surf->particle_mu[last_particle_idx];

              double tau = surf->tau;
              double vx = surf->vx;
              double vy = surf->vy;

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

              double v2 = vx*vx + vy*vy;
              double gammaT = 1.0/sqrt(1.0 - v2);


              for (int k=0; k<eta_tab_length; k++)
              {
                  //double eta_s = eta_tab->get(1,k+1); // unused

                  // old way
                  //double delta_eta = eta_tab->get(2,k+1);
                  // new way
                  double delta_eta = delta_eta_tab[k];

                  // old way
                  //double cosh_y_eta = cosh(y-eta_s);
                  //double sinh_y_eta = sinh(y-eta_s);

                  // old way
                  //double pt = mT*cosh_y_eta;
                  //double pz = mT*sinh_y_eta;

                  // new way
                  double pt = mT*hypertrig_etas_table[k][0];
                  double pz = mT*hypertrig_etas_table[k][1];

                  double expon = (gammaT*(pt*1 - px*vx - py*vy) - mu) / Tdec;
                  double f0 = 1./(exp(expon)+sign);       //thermal equilibrium distributions

                  // Must adjust this to be correct for the p*del \tau term.  The plus sign is
                  // due to the fact that the DA# variables are for the covariant surface integration
                  double pdsigma = pt*da0 + px*da1 + py*da2;

                  //viscous corrections
                  double Wfactor = pt*pt*pi00 - 2.0*pt*px*pi01 - 2.0*pt*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33;
                  double deltaf = (1 - F0_IS_NOT_SMALL*sign*f0)*Wfactor*deltaf_prefactor;
                  double result;
                  if(1+deltaf < 0.0) //set results to zero when delta f turns whole expression to negative
                     result = 0.0;
                  else
                     result = prefactor*degen*f0*(1+deltaf)*pdsigma*tau;

                  dN_ptdptdphidy_tmp += result*delta_eta;
              } // k
          } // l

          dN_ptdptdphidy_tab[i][j] = dN_ptdptdphidy_tmp;
          if (AMOUNT_OF_OUTPUT>0) print_progressbar((i*phi_tab_length+j)/progress_total);
      }
  //cout << int(100.0*(i+1)/pT_tab_length) << "% completed" << endl;
  }
  if (AMOUNT_OF_OUTPUT>0) print_progressbar(1);

  for (int i=0; i<pT_tab_length; i++)
  for (int j=0; j<phi_tab_length; j++)
    dN_ptdptdphidy->set(i+1,j+1,dN_ptdptdphidy_tab[i][j]);

  sw.toc();
  cout << endl << "Finished " << sw.takeTime() << " seconds." << endl;
}

void EmissionFunctionArray::write_dN_ptdptdphidy_toFile()
// Append the dN_xxx results to file.
{
  ofstream of1(dN_ptdptdphidy_filename.c_str(), ios_base::app);
  dN_ptdptdphidy->printTable(of1);
  of1.close();
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






//*********************************************************************************************
void EmissionFunctionArray::calculate_dN_ptdptdphidy_and_flows_4all(int to_order)
// Calculate dNArrays and flows for all particles given in chosen_particle file.
{
    if (USE_HISTORIC_FORMAT)
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

                char buffer_diff[500], buffer_inte[500];
                sprintf(buffer_diff, flow_differential_filename.c_str(), monval);
                remove(buffer_diff);
                sprintf(buffer_inte, flow_integrated_filename.c_str(), monval);
                remove(buffer_inte);
                calculate_flows(to_order, buffer_diff, buffer_inte);
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

            sw.toc();
            cout << " -- Calculate_dN_ptdptdphidy_and_flows_4all finishes " << sw.takeTime() << " seconds." << endl;

    }

}





//**********************************************************************************************
bool EmissionFunctionArray::particles_are_the_same(int idx1, int idx2)
{
    if (particles[idx1].sign!=particles[idx2].sign) return false;
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
