#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include "tensor_trans.h"
using namespace std;

//lorentz matrix to first boost along longitudinal direction and then transverse
//direction
void boost_matrix(double** lambda_munu, double vx, double vy, double vz)
{
      double gamma_perp, gamma_L, gamma;
      double beta, beta_1, beta_2, beta_3;
      double eps = 1e-100;
      
      gamma_perp = 1/sqrt(1 - vx*vx - vy*vy + eps);
      gamma_L = 1/sqrt(1 - vz*vz + eps);
      gamma = gamma_perp * gamma_L;
      /*if(isnan(gamma))
      {
         cout<<"gamma is nan!" << endl;
         cout<<"v:" << vx << "  " << vy << "   " << vz << endl;
         exit(1);
      }*/
      beta_1 = vx/gamma_L;
      beta_2 = vy/gamma_L;
      beta_3 = vz;
      beta = sqrt(beta_1*beta_1 + beta_2*beta_2 + beta_3*beta_3 + eps);
      
      lambda_munu[0][0] = gamma;
      lambda_munu[0][1] = -gamma*beta_1;
      lambda_munu[0][2] = -gamma*beta_2;
      lambda_munu[0][3] = -gamma*beta_3;
      lambda_munu[1][0] = lambda_munu[0][1];
      lambda_munu[1][1] = 1+(gamma-1)*beta_1*beta_1/beta/beta;
      lambda_munu[1][2] = (gamma-1)*beta_1*beta_2/beta/beta;
      lambda_munu[1][3] = (gamma-1)*beta_1*beta_3/beta/beta;
      lambda_munu[2][0] = lambda_munu[0][2];
      lambda_munu[2][1] = lambda_munu[1][2];
      lambda_munu[2][2] = 1+(gamma-1)*beta_2*beta_2/beta/beta;
      lambda_munu[2][3] = (gamma-1)*beta_2*beta_3/beta/beta;
      lambda_munu[3][0] = lambda_munu[0][3];
      lambda_munu[3][1] = lambda_munu[1][3];
      lambda_munu[3][2] = lambda_munu[2][3];
      lambda_munu[3][3] = 1+(gamma-1)*beta_3*beta_3/beta/beta;
}

void getTransverseflow_u_mu_low(double* flow_u_mu_low, double vx, double vy)
{
      double gamma;
      double eps = 1e-100;
      
      gamma = 1./sqrt(1. - vx*vx - vy*vy + eps);
      
      flow_u_mu_low[0] = gamma;
      flow_u_mu_low[1] = - gamma*vx;
      flow_u_mu_low[2] = - gamma*vy;
      flow_u_mu_low[3] = 0.0;
}

void boost_vec_trans(double* p, double* p_prime, double** lambda)
{
      for(int i=0;i<4;i++) p_prime[i] = 0;
      for(int i=0;i<4;i++) 
         for(int j=0;j<4;j++)
             p_prime[i] += lambda[i][j]*p[j];
}

void boost_Tensor2_trans(double** M, double** M_prime, double** lambda)
{
      for(int i=0;i<4;i++) 
         for(int j=0;j<4;j++)
            M_prime[i][j] = 0;
      for(int i=0; i<4; i++)
         for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
               for(int l=0;l<4;l++)
                  M_prime[i][j] +=  lambda[i][k]*lambda[j][l]*M[k][l];
}

void RotationMatrix(double** RotationM, double phi, 
                    double theta, double psi)
{
      double cos_phi = cos(phi);
      double sin_phi = sin(phi);
      double cos_theta = cos(theta);
      double sin_theta = sin(theta);
      double cos_psi = cos(psi);
      double sin_psi = sin(psi);
      for(int i=0;i<4;i++) RotationM[0][i] = RotationM[i][0] = 0.0;
      RotationM[0][0] = 1.0;
      RotationM[1][1] = cos_phi*cos_psi - sin_phi*cos_theta*sin_psi;
      RotationM[1][2] = -cos_phi*sin_psi-sin_phi*cos_theta*cos_psi;
      RotationM[1][3] = sin_theta*sin_phi;
      RotationM[2][1] = sin_phi*cos_psi + cos_phi*cos_theta*sin_psi;
      RotationM[2][2] = -sin_phi*sin_psi + cos_phi*cos_theta*cos_psi;
      RotationM[2][3] = -sin_theta*cos_phi;
      RotationM[3][1] = sin_theta*sin_psi;
      RotationM[3][2] = sin_theta*cos_psi;
      RotationM[3][3] = cos_theta;
}

void Rotation_vec_trans(double* p, double* p_prime, double** RotationM)
{
      for(int i=0;i<4;i++) p_prime[i] = 0;
      for(int i=0; i<4; i++)
         for(int j=0; j<4; j++)
              p_prime[i] +=  RotationM[i][j]*p[j];
}

void Rotation_Tensor2_trans(double** M, double** M_prime, double** RotationM)
{
      for(int i=0;i<4;i++) 
         for(int j=0;j<4;j++)
            M_prime[i][j] = 0;
      for(int i=0; i<4; i++)
         for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
               for(int l=0;l<4;l++)
                  M_prime[i][j] +=  RotationM[i][k]*RotationM[j][l]*M[k][l];
}

void Rotation_Matrix_R_z_i(double* R_z_i, double* vec)
{
      double cos_theta = vec[3]/vec[0];
      double sin_theta = sqrt(1 - cos_theta*cos_theta);
      double demon_inverse_temp = 1./(vec[0]*sin_theta);
      double cos_phi = vec[1]*demon_inverse_temp;
      double sin_phi = vec[2]*demon_inverse_temp;
      R_z_i[0] = 0.0e0;
      R_z_i[1] = sin_theta*cos_phi;
      R_z_i[2] = sin_theta*sin_phi;
      R_z_i[3] = cos_theta;
}

double Rotation_Tensor_zz(double** M, double* R_z_i)
{
      double pi_zz = 2*R_z_i[1]*R_z_i[2]*M[1][2] + 2*R_z_i[1]*R_z_i[3]*M[1][3] 
                     + 2*R_z_i[2]*R_z_i[3]*M[2][3] + R_z_i[1]*R_z_i[1]*M[1][1] 
                     + R_z_i[2]*R_z_i[2]*M[2][2] + R_z_i[3]*R_z_i[3]*M[3][3];
      return(pi_zz);
}
