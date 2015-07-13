#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Regge96.h"
#include "profile_3d.h"


using namespace std;

profile_3d::profile_3d(
    vector<participant_info> part_list_in, 
    vector<participant_info> binary_list_in, 
    int grid_nx_in, int grid_ny_in, int grid_neta_in, 
    double grid_dx_in, double grid_dy_in, double grid_deta_in, 
    double ecm_in, int random_flag_in)
{
    ecm = ecm_in;
    y_beam = atanh(sqrt(1. - 1./pow(ecm/2., 2)));
    double sig = hadronxsec::totalXsection(ecm, 0);   // in [mb]
    double sigel = hadronxsec::elasticXsection(sig, ecm, 0, 0);  // in [mb]
    sigma_inelastic = (sig - sigel)*0.1;   // in [fm^2]

    cout << "ecm = " << ecm << " A GeV, "
         << "y_beam = " << y_beam << ", "
         << "sigmaNN = " << sigma_inelastic*10. << " mb. "<< endl;

    for(int i = 0; i < part_list_in.size(); i++)
    {
       participant_info temp_part;
       temp_part.x = part_list_in[i].x;
       temp_part.y = part_list_in[i].y;
       temp_part.id = part_list_in[i].id;
       participant_list.push_back(temp_part);
    }

    for(int i = 0; i < binary_list_in.size(); i++)
    {
       participant_info temp_part;
       temp_part.x = binary_list_in[i].x;
       temp_part.y = binary_list_in[i].y;
       temp_part.id = -1;
       binary_collision_list.push_back(temp_part);
    }

    grid_nx = grid_nx_in;
    grid_ny = grid_ny_in;
    grid_neta = grid_neta_in;

    rho_part = new double** [grid_neta];
    rho_binary = new double** [grid_neta];
    for(int i = 0; i < grid_neta; i++)
    {
        rho_part[i] = new double* [grid_nx];
        rho_binary[i] = new double* [grid_nx];
        for(int j = 0; j < grid_nx; j++)
        {
            rho_part[i][j] = new double [grid_ny];
            rho_binary[i][j] = new double [grid_ny];
            for(int k = 0; k < grid_ny; k++)
            {
                rho_part[i][j][k] = 0.0;
                rho_binary[i][j][k] = 0.0;
            }
        }
    }

    grid_dx = grid_dx_in;
    grid_dy = grid_dy_in;
    grid_deta = grid_deta_in;

    x_array = new double [grid_nx];
    for(int i = 0; i < grid_nx; i++)
        x_array[i] = (i - (grid_nx - 1)/2.)*grid_dx;
    
    y_array = new double [grid_ny];
    for(int i = 0; i < grid_ny; i++)
        y_array[i] = (i - (grid_ny - 1)/2.)*grid_dy;

    eta_array = new double [grid_neta];
    for(int i = 0; i < grid_neta; i++)
        eta_array[i] = (i - (grid_neta - 1)/2.)*grid_deta;

    random_flag = random_flag_in;

    gsl_rng_env_setup();
    gsl_random_type = gsl_rng_default;
    gsl_random_rng = gsl_rng_alloc (gsl_random_type);
    long seed = time (NULL);
    cout << "random seed: " << seed << endl;
    gsl_rng_set (gsl_random_rng, seed);
    
    set_eta_distribution();
    set_variables();
}

profile_3d::~profile_3d()
{
    participant_list.clear();
    binary_collision_list.clear();

    eta_distribution_array_eta.clear();
    eta_distribution_array_projectile.clear();
    eta_distribution_array_target.clear();

    for(int i = 0; i < grid_neta; i++)
    {
        for(int j = 0; j < grid_nx; j++)
        {
            delete[] rho_part[i][j];
            delete[] rho_binary[i][j];
        }
        delete[] rho_part[i];
        delete[] rho_binary[i];
    }
    delete[] rho_part;
    delete[] rho_binary;

    delete[] x_array;
    delete[] y_array;
    delete[] eta_array;
}

void profile_3d::set_eta_distribution()
{
    double eta_peak = 2.5;
    double sigma_eta_out = 0.5;

    int length = 1000;
    double deta = 2*y_beam/(length - 1);

    eta_distribution_proj_max = 0.0;
    eta_distribution_targ_max = 0.0;
    for(int i = 0; i < length; i++)
    {
        double eta_local = - y_beam + i*deta;
        double exp_factor = 1.0;
        if(fabs(eta_local) > eta_peak)
            exp_factor = exp( - pow((fabs(eta_local) - eta_peak), 2)
                                /(2.*sigma_eta_out*sigma_eta_out));
        double eta_dis_left = (1. - eta_local/y_beam)*exp_factor;
        double eta_dis_right = (1. + eta_local/y_beam)*exp_factor;
        if (eta_dis_left > eta_distribution_proj_max)
            eta_distribution_proj_max = eta_dis_left;
        if (eta_dis_right > eta_distribution_targ_max)
            eta_distribution_targ_max = eta_dis_right;
        eta_distribution_array_eta.push_back(eta_local);
        eta_distribution_array_projectile.push_back(eta_dis_left);
        eta_distribution_array_target.push_back(eta_dis_right);
    }
}

double profile_3d::sample_eta_distribution_from_array(int participant_id)
{
    double eta_sample;
    double probability, prob_max;
    if (participant_id == 1)
        prob_max = eta_distribution_proj_max;
    else
        prob_max = eta_distribution_targ_max;

    do
    {
        eta_sample = (
            y_beam*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));

        int eta_idx = binarySearch(&eta_distribution_array_eta, eta_sample, 
                                   true);
        double eta_frac = ((eta_sample - eta_distribution_array_eta[eta_idx])
                           /(eta_distribution_array_eta[eta_idx+1] 
                             - eta_distribution_array_eta[eta_idx]));
        
        if (participant_id == 1)
        {
            probability = (
                eta_distribution_array_projectile[eta_idx]*(1. - eta_frac) 
                + eta_distribution_array_projectile[eta_idx+1]*eta_frac);
        }
        else
        {
            probability = (
                eta_distribution_array_target[eta_idx]*(1. - eta_frac) 
                + eta_distribution_array_target[eta_idx+1]*eta_frac);
        }
    }while(prob_max*gsl_rng_uniform(gsl_random_rng) > probability);

    return(eta_sample);
}


void profile_3d::set_variables()
{
    double eta_0 = 2.0;
    double sigma_x_0 = sqrt(sigma_inelastic/(8*M_PI));
    double sigma_y_0 = sqrt(sigma_inelastic/(8*M_PI));
    double sigma_eta_0 = 0.5;

    double dsigma_x = 0.3;
    double dsigma_y = 0.3;
    double dsigma_eta = 0.3;

    for(int i = 0; i < participant_list.size(); i++)
    {
        if(random_flag == 0)
        {
            if(participant_list[i].id == 1)
                participant_list[i].eta = eta_0;
            else
                participant_list[i].eta = -eta_0;
            participant_list[i].sigma_x = sigma_x_0;
            participant_list[i].sigma_y = sigma_y_0;
            participant_list[i].sigma_eta = sigma_eta_0;
        }
        else if(random_flag == 1)
        {
            //participant_list[i].eta = (
            //    eta_0*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));
            participant_list[i].eta = 
                sample_eta_distribution_from_array(participant_list[i].id);
            participant_list[i].sigma_x = sigma_x_0;
            participant_list[i].sigma_y = sigma_y_0;
            participant_list[i].sigma_eta = sigma_eta_0;
        }
        else if (random_flag == 2)
        {
            participant_list[i].eta = 
                sample_eta_distribution_from_array(participant_list[i].id);
            participant_list[i].sigma_eta = (sigma_eta_0 
                + dsigma_eta*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));

            double sigma_trans = (sigma_x_0 
                + dsigma_x*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));
            participant_list[i].sigma_x = sigma_trans;
            participant_list[i].sigma_y = sigma_trans;
        }
        else if (random_flag == 3)
        {
            participant_list[i].eta = 
                sample_eta_distribution_from_array(participant_list[i].id);
            participant_list[i].sigma_x = (sigma_x_0 
                + dsigma_x*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));
            participant_list[i].sigma_y = (sigma_y_0 
                + dsigma_y*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));
            participant_list[i].sigma_eta = (sigma_eta_0 
                + dsigma_eta*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));
        }
    }

    for(int i = 0; i < binary_collision_list.size(); i++)
    {
        if(random_flag == 0)
        {
            binary_collision_list[i].eta = eta_0;
            binary_collision_list[i].sigma_x = sigma_x_0;
            binary_collision_list[i].sigma_y = sigma_y_0;
            binary_collision_list[i].sigma_eta = sigma_eta_0;
        }
        else if(random_flag == 1)
        {
            binary_collision_list[i].eta = 
                sample_eta_distribution_from_array(participant_list[i].id);
            binary_collision_list[i].sigma_x = sigma_x_0;
            binary_collision_list[i].sigma_y = sigma_y_0;
            binary_collision_list[i].sigma_eta = sigma_eta_0;
        }
        else if (random_flag == 2)
        {
            binary_collision_list[i].eta = 
                sample_eta_distribution_from_array(participant_list[i].id);
            binary_collision_list[i].sigma_eta = (sigma_eta_0 
                + dsigma_eta*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));

            double sigma_trans = (sigma_x_0 
                + dsigma_x*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));
            binary_collision_list[i].sigma_x = sigma_trans;
            binary_collision_list[i].sigma_y = sigma_trans;
        }
        else if (random_flag == 3)
        {
            binary_collision_list[i].eta = 
                sample_eta_distribution_from_array(participant_list[i].id);
            binary_collision_list[i].sigma_x = (sigma_x_0 
                + dsigma_x*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));
            binary_collision_list[i].sigma_y = (sigma_y_0 
                + dsigma_y*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));
            binary_collision_list[i].sigma_eta = (sigma_eta_0 
                + dsigma_eta*(1. - 2.*gsl_rng_uniform(gsl_random_rng)));
        }
    }
}

void profile_3d::generate_3d_profile()
{
    for(int i = 0; i < participant_list.size(); i++)
    {
        double part_x = participant_list[i].x;
        double part_y = participant_list[i].y;
        double part_eta = participant_list[i].eta;
        int idx_x0 = (int)(part_x/grid_dx + (grid_nx - 1)/2);
        int idx_y0 = (int)(part_y/grid_dy + (grid_ny - 1)/2);
        int idx_eta0 = (int)(part_eta/grid_deta + (grid_neta - 1)/2);

        double part_sigma_x = participant_list[i].sigma_x;
        int idx_x_range = (int)(6*part_sigma_x/grid_dx);
        int idx_x_left = max(idx_x0 - idx_x_range, 0);
        int idx_x_right = min(idx_x0 + idx_x_range, grid_nx);
        double part_sigma_y = participant_list[i].sigma_y;
        int idx_y_range = (int)(6*part_sigma_y/grid_dy);
        int idx_y_left = max(idx_y0 - idx_y_range, 0);
        int idx_y_right = min(idx_y0 + idx_y_range, grid_ny);
        double part_sigma_eta = participant_list[i].sigma_eta;
        int idx_eta_range = (int)(6*part_sigma_eta/grid_deta);
        int idx_eta_left = max(idx_eta0 - idx_eta_range, 0);
        int idx_eta_right = min(idx_eta0 + idx_eta_range, grid_neta);

        for(int j = idx_eta_left; j < idx_eta_right; j++)
        {
            double dis_eta = (
                (eta_array[j] - part_eta)*(eta_array[j] - part_eta)
                /(2.*part_sigma_eta*part_sigma_eta));

            double norm_eta = 1./sqrt(2.*M_PI*part_sigma_eta*part_sigma_eta);
            for(int k = idx_x_left; k < idx_x_right; k++)
            {
                double dis_x = ((x_array[k] - part_x)*(x_array[k] - part_x)
                    /(2.*part_sigma_x*part_sigma_x));
                double norm_x = 1./sqrt(2.*M_PI*part_sigma_x*part_sigma_x);
                for(int l = idx_y_left; l < idx_y_right; l++)
                {
                    double dis_y = ((y_array[l] - part_y)*(y_array[l] - part_y)
                        /(2.*part_sigma_y*part_sigma_y));
                    double norm_y = 1./sqrt(2.*M_PI*part_sigma_y*part_sigma_y);
                    double result = (
                        exp(- dis_eta - dis_x - dis_y)*norm_eta*norm_x*norm_y);
                    rho_part[j][k][l] += result;
                }
            }
        }
    }
}

void profile_3d::generate_3d_profile_binary()
{
    for(int i = 0; i < binary_collision_list.size(); i++)
    {
        double part_x = binary_collision_list[i].x;
        double part_y = binary_collision_list[i].y;
        double part_eta = binary_collision_list[i].eta;
        int idx_x0 = (int)(part_x/grid_dx + (grid_nx - 1)/2);
        int idx_y0 = (int)(part_y/grid_dy + (grid_ny - 1)/2);
        int idx_eta0 = (int)(part_eta/grid_deta + (grid_neta - 1)/2);

        double part_sigma_x = binary_collision_list[i].sigma_x;
        int idx_x_range = (int)(6*part_sigma_x/grid_dx);
        int idx_x_left = max(idx_x0 - idx_x_range, 0);
        int idx_x_right = min(idx_x0 + idx_x_range, grid_nx);
        double part_sigma_y = binary_collision_list[i].sigma_y;
        int idx_y_range = (int)(6*part_sigma_y/grid_dy);
        int idx_y_left = max(idx_y0 - idx_y_range, 0);
        int idx_y_right = min(idx_y0 + idx_y_range, grid_ny);
        double part_sigma_eta = binary_collision_list[i].sigma_eta;
        int idx_eta_range = (int)(6*part_sigma_eta/grid_deta);
        int idx_eta_left = max(idx_eta0 - idx_eta_range, 0);
        int idx_eta_right = min(idx_eta0 + idx_eta_range, grid_neta);

        for(int j = idx_eta_left; j < idx_eta_right; j++)
        {
            double dis_eta = (
                (eta_array[j] - part_eta)*(eta_array[j] - part_eta)
                /(2.*part_sigma_eta*part_sigma_eta));

            double norm_eta = 1./sqrt(2.*M_PI*part_sigma_eta*part_sigma_eta);
            for(int k = idx_x_left; k < idx_x_right; k++)
            {
                double dis_x = ((x_array[k] - part_x)*(x_array[k] - part_x)
                    /(2.*part_sigma_x*part_sigma_x));
                double norm_x = 1./sqrt(2.*M_PI*part_sigma_x*part_sigma_x);
                for(int l = idx_y_left; l < idx_y_right; l++)
                {
                    double dis_y = ((y_array[l] - part_y)*(y_array[l] - part_y)
                        /(2.*part_sigma_y*part_sigma_y));
                    double norm_y = 1./sqrt(2.*M_PI*part_sigma_y*part_sigma_y);
                    double result = (
                        exp(- dis_eta - dis_x - dis_y)*norm_eta*norm_x*norm_y);
                    rho_binary[j][k][l] += result;
                }
            }
        }
    }
}

void profile_3d::output_3d_entropy_profile(string filename, double alpha)
{
    ofstream of(filename.c_str());
    for(int i = 0; i < grid_neta; i++)
    {
        for(int j = 0; j < grid_nx; j++)
        {
            for(int k = 0; k < grid_ny; k++)
            {
                double sd_local = (rho_part[i][j][k]*(1. - alpha)/2. 
                                   + rho_binary[i][j][k]*alpha);
                of << scientific << setw(18) << setprecision(8)
                   << sd_local << "   ";
            }
            of << endl;
        }
    }
    of.close();
}

void profile_3d::output_3d_rhob_profile(string filename)
{
    ofstream of(filename.c_str());
    for(int i = 0; i < grid_neta; i++)
    {
        for(int j = 0; j < grid_nx; j++)
        {
            for(int k = 0; k < grid_ny; k++)
            {
                double rhob_local = rho_part[i][j][k];
                of << scientific << setw(18) << setprecision(8)
                   << rhob_local << "   ";
            }
            of << endl;
        }
    }
    of.close();
}

long profile_3d::binarySearch(vector<double>* A, double value, bool skip_out_of_range)
// Return the index of the largest number less than value in the list A
// using binary search. Index starts with 0.
// If skip_out_of_range is set to true, then it will return -1 for those
// samples that are out of the table range.
{
   int length = A->size();
   int idx_i, idx_f, idx;
   idx_i = 0;
   idx_f = length-1;
   if(value > (*A)[idx_f])
   {
      if (skip_out_of_range) return -1;
      cout << "binarySearch: desired value is too large, exceeding the end of the table." << endl;
      exit(-1);
   }
   if(value < (*A)[idx_i])
   {
      if (skip_out_of_range) return -1;
      cout << "binarySearch: desired value is too small, exceeding the beginning of table." << endl;
      exit(-1);
   }
   idx = (int) floor((idx_f+idx_i)/2.);
   while((idx_f-idx_i) > 1)
   {
     if((*A)[idx] < value)
        idx_i = idx;
     else
        idx_f = idx;
     idx = (int) floor((idx_f+idx_i)/2.);
   }
   return(idx_i);
}
