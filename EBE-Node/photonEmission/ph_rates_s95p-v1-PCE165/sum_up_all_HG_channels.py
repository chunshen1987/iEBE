#!/usr/bin/env python

from numpy import *

# full list of the mesonic 2 to 2 channels
#channel_list = ['K_Kstar_to_pion_gamma', 
#                'pion_K_to_Kstar_gamma', 
#                'pion_Kstar_to_K_gamma',
#                'pion_pion_to_rho_gamma', 
#                'pion_rho_to_omega_to_pion_gamma', 
#                'pion_rho_to_pion_gamma',
#                'rho_K_to_K_gamma',
#                'rho_to_pion_pion_gamma']

# the list of the mesonic 2 to 2 channels when pipi bremsstrahlung is included
channel_list = ['K_Kstar_to_pion_gamma', 
                'pion_K_to_Kstar_gamma', 
                'pion_Kstar_to_K_gamma',
                'pion_rho_to_omega_to_pion_gamma', 
                'pion_rho_to_pion_gamma',
                'rho_K_to_K_gamma']

rate_type_list = ['eqrate', 'viscous', 'bulkvis']

sum_data = 0.0*loadtxt('rate_%s_%s.dat' % (channel_list[0], rate_type_list[0]))
for itype in range(len(rate_type_list)):
    for ich in range(len(channel_list)):
        filename = 'rate_%s_%s.dat' % (channel_list[ich], rate_type_list[itype])
        data = loadtxt(filename)
        sum_data += data
    savetxt('rate_HG_2to2_meson_total_%s.dat' % rate_type_list[itype], 
            sum_data, fmt='%.10e', delimiter = '   ')
