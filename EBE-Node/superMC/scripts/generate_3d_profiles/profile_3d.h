#ifndef profile_3d_H
#define profile_3d_H

#include <vector>
#include <string>

#include <gsl/gsl_rng.h>

using namespace std;

struct participant_info
{
    int id;
    double x, y, eta;
    double sigma_x, sigma_y, sigma_eta;
};

class profile_3d
{
    private:
        vector<participant_info> participant_list;
        vector<participant_info> binary_collision_list;

        int grid_nx, grid_ny, grid_neta;
        double grid_dx, grid_dy, grid_deta;
        double *x_array, *y_array, *eta_array;
        double ***rho_part, ***rho_binary;

        double eta_distribution_proj_max, eta_distribution_targ_max;
        vector<double> eta_distribution_array_eta;
        vector<double> eta_distribution_array_projectile;
        vector<double> eta_distribution_array_target;

        double ecm, y_beam, sigma_inelastic;

        int random_flag;
    
        const gsl_rng_type * gsl_random_type;
        gsl_rng * gsl_random_rng;

    public:
        profile_3d(vector<participant_info> part_list_in, 
                   vector<participant_info> binary_list_in, 
                   int grid_nx_in, int grid_ny_in, int grid_neta_in, 
                   double grid_dx_in, double grid_dy_in, double grid_deta_in, 
                   double ecm_in, int random_flag_in);
        ~profile_3d();

        void set_variables();
        void set_eta_distribution();
        void generate_3d_profile();
        void generate_3d_profile_binary();
        void output_3d_entropy_profile(string filename, double alpha);
        void output_3d_rhob_profile(string filename);
        double sample_eta_distribution_from_array(int participant_id);

        long binarySearch(vector<double>* A, double value, bool skip_out_of_range);

};

#endif
