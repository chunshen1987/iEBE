#ifndef TENSOR_TRANS_H
#define TNESOR_TRANS_H

void boost_matrix(double** , double , double , double);
void getTransverseflow_u_mu_low(double* flow_u_mu_low, double vx, double vy);
void boost_vec_trans(double* , double* , double** );
void boost_Tensor2_trans(double** , double** , double** );

void RotationMatrix(double** , double , double , double );
void Rotation_vec_trans(double* , double* , double** );
void Rotation_Tensor2_trans(double** , double** , double** );
void Rotation_Matrix_R_z_i(double* R_z_i, double* vec);
double Rotation_Tensor_zz(double** M, double* R_z_i);

#endif
