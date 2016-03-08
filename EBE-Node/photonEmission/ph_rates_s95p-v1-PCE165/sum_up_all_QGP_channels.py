#!/usr/bin/env python

from numpy import *

channel_list = ['QGP_2to2_hard', 'QGP_2to2_soft' ]

rate_type_list = ['eqrate', 'viscous', 'bulkvis']

sum_data = 0.0*loadtxt('rate_%s_%s.dat' % (channel_list[0], rate_type_list[0]))
for itype in range(len(rate_type_list)):
    for ich in range(len(channel_list)):
        filename = 'rate_%s_%s.dat' % (channel_list[ich], rate_type_list[itype])
        data = loadtxt(filename)
        sum_data += data
    savetxt('rate_QGP_2to2_total_%s.dat' % rate_type_list[itype], 
            sum_data, fmt='%.10e', delimiter = '   ')
