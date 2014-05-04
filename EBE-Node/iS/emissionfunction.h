#ifndef EMISSIONFUNCTION_H
#define EMISSIONFUNCTION_H

#include<string>
#include<vector>
#include "Table.h"
#include "main.h"
using namespace std;

class EmissionFunctionArray
{
private:
  double particle_y;
  Table *pT_tab, *phi_tab, *eta_tab;
  int pT_tab_length, phi_tab_length, eta_tab_length;
  long FO_length;
  Table *dN_ptdptdphidy;
  Table *dE_ptdptdphidy;
  int number_of_chosen_particles;
  int *chosen_particles_01_table; // has length Nparticle, 0 means miss, 1 means include
  int *chosen_particles_sampling_table; // store particle index; the sampling process follows the order specified by this table
  int Nparticles;
  particle_info* particles;
  FO_surf* FOsurf_ptr;
  int last_particle_idx; // store the last particle index being used by calculate_dN_ptdptdphidy function
  bool particles_are_the_same(int, int);

  //array for bulk delta f coefficients
  Table *bulkdf_coeff;

public:
  EmissionFunctionArray(double particle_y_in, Table* chosen_particle, Table* pT_tab_in, Table* phi_tab_in, Table* eta_tab_in, particle_info* particles_in, int Nparticles, FO_surf* FOsurf_ptr_in, long FO_length_in);
  ~EmissionFunctionArray();

  void calculate_dN_ptdptdphidy(int);
  void write_dN_ptdptdphidy_toFile();
  string dN_ptdptdphidy_filename; // where to save
  string dE_ptdptdphidy_filename; // where to save

  void calculate_flows(int to_order, string, string);
  void calculate_Energyflows(int to_order, string, string);
  string flow_differential_filename_old, flow_integrated_filename_old;
  string flow_differential_filename, flow_integrated_filename;
  string energyflow_differential_filename_old, energyflow_integrated_filename_old;
  string energyflow_differential_filename, energyflow_integrated_filename;

  void calculate_dN_ptdptdphidy_and_flows_4all(int to_order=9);
  void getbulkvisCoefficients(double Tdec, double* bulkvisCoefficients);

};

#endif
