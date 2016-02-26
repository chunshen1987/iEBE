#!/usr/bin/env python3

import sys
from numpy import *
from scipy import interpolate

input_filename = sys.argv[1]

input_data = loadtxt(input_filename, skiprows=1)
n_row, n_col = input_data.shape
n_dim = int(sqrt(n_row))
dx = input_data[1, 2] - input_data[0, 2]
print("input grid size: %d, dx = %g fm" % (n_dim, dx))

input_e_field = zeros([n_dim, n_dim])
input_utau_field = zeros([n_dim, n_dim])
input_ux_field = zeros([n_dim, n_dim])
input_uy_field = zeros([n_dim, n_dim])

idx = 0
for i in range(n_dim):
    for j in range(n_dim):
        input_e_field[i][j] = input_data[idx, 3]
        input_utau_field[i][j] = input_data[idx, 4]
        input_ux_field[i][j] = input_data[idx, 5]
        input_uy_field[i][j] = input_data[idx, 6]
        idx += 1


input_x_array = linspace(input_data[0, 1], input_data[-1, 1], n_dim)
input_y_array = linspace(input_data[0, 2], input_data[n_dim-1, 2], n_dim)
input_x_grid, input_y_grid = meshgrid(input_x_array, input_y_array)

print("from input file:")
print("-"*80)
print("total dE/dy: %g" % (sum(input_e_field)*dx*dx) )
r_grid = sqrt(input_x_grid**2. + input_y_grid**2.)
phi_grid = arctan2(input_y_grid, input_x_grid)
ecc_1_real = (
    sum(r_grid**3.*cos(phi_grid)*input_e_field)/sum(r_grid**3.*input_e_field))
ecc_1_imag = (
    sum(r_grid**3.*sin(phi_grid)*input_e_field)/sum(r_grid**3.*input_e_field))
print("ecc_1: %g" % (sqrt(ecc_1_real**2. + ecc_1_imag**2.)))
for iorder in range(2, 6):
    ecc_n_real = (
        sum(r_grid**(iorder)*cos(iorder*phi_grid)*input_e_field)
        /sum(r_grid**(iorder)*input_e_field))
    ecc_n_imag = (
        sum(r_grid**(iorder)*sin(iorder*phi_grid)*input_e_field)
        /sum(r_grid**(iorder)*input_e_field))
    print("ecc_%d: %g" % (iorder, sqrt(ecc_n_real**2. + ecc_n_imag**2.)))
print("-"*80)

interp_e_field = interpolate.RegularGridInterpolator(
                     (input_x_array, input_y_array), input_e_field,
                     bounds_error=False, fill_value=1e-18)
interp_utau_field = interpolate.RegularGridInterpolator(
                     (input_x_array, input_y_array), input_utau_field,
                     bounds_error=False, fill_value=1e-18)
interp_ux_field = interpolate.RegularGridInterpolator(
                     (input_x_array, input_y_array), input_ux_field,
                     bounds_error=False, fill_value=1e-18)
interp_uy_field = interpolate.RegularGridInterpolator(
                     (input_x_array, input_y_array), input_uy_field,
                     bounds_error=False, fill_value=1e-18)

dim = 341
output_x_array = linspace(-17., 17., dim)
output_y_array = linspace(-17., 17., dim)
dx = output_x_array[1] - output_x_array[0]
print("output grid size: %d, dx = %g fm" % (dim, dx))
output_x_grid, output_y_grid = meshgrid(output_x_array, output_y_array)
output_grid = []
for i in range(dim):
    for j in range(dim):
        output_grid.append([output_x_array[i], output_y_array[j]])
output_e_field = interp_e_field(array(output_grid))
output_utau_field = interp_utau_field(array(output_grid))
output_ux_field = interp_ux_field(array(output_grid))
output_uy_field = interp_uy_field(array(output_grid))

output_e_field_mat = output_e_field.reshape(dim, dim)

print("from output file:")
print("-"*80)
print("total dE/dy: %g" % (sum(output_e_field_mat)*dx*dx) )
r_grid = sqrt(output_x_grid**2. + output_y_grid**2.)
phi_grid = arctan2(output_y_grid, output_x_grid)
ecc_1_real = (
    sum(r_grid**3.*cos(phi_grid)*output_e_field_mat)
    /sum(r_grid**3.*output_e_field_mat))
ecc_1_imag = (
    sum(r_grid**3.*sin(phi_grid)*output_e_field_mat)
    /sum(r_grid**3.*output_e_field_mat))
print("ecc_1: %g" % (sqrt(ecc_1_real**2. + ecc_1_imag**2.)))
for iorder in range(2, 6):
    ecc_n_real = (
        sum(r_grid**(iorder)*cos(iorder*phi_grid)*output_e_field_mat)
        /sum(r_grid**(iorder)*output_e_field_mat))
    ecc_n_imag = (
        sum(r_grid**(iorder)*sin(iorder*phi_grid)*output_e_field_mat)
        /sum(r_grid**(iorder)*output_e_field_mat))
    print("ecc_%d: %g" % (iorder, sqrt(ecc_n_real**2. + ecc_n_imag**2.)))
print("-"*80)

output_data = []
idx = 0
for i in range(dim):
    for j in range(dim):
        output_data.append(
            [output_x_array[i], output_y_array[j], output_e_field[idx], 
             output_ux_field[idx], output_uy_field[idx]])
        idx += 1
output_filename = input_filename.split('.dat')[0] + "_regulated.dat"
savetxt(output_filename, array(output_data), fmt = '%.10e', delimiter = '  ')
