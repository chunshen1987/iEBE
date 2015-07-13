#! /usr/bin/env python

import sys, shutil
from numpy import *
from os import path, makedirs
from glob import glob
import subprocess

class color:
    """
    define colors in the terminal
    """
    purple = '\033[95m'
    cyan = '\033[96m'
    darkcyan = '\033[36m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    bold = '\033[1m'
    underline = '\033[4m'
    end = '\033[0m'

grid_nx = 261
grid_ny = 261
grid_nrap = 101
grid_dx = 0.1
grid_dy = 0.1
grid_drap = 0.1

rand_flag = 1

# the width of the Gaussian in the transverse plane
sigma_perp = 0.5; delta_sigma_perp = 0.3

# peak position in the longitudinal direction
eta_0 = 2.0; fluct_eta_0 = 1.0
# the width of the Gaussian in the longitudinal direction
sigma_beam_in = 1.0; delta_sigma_beam_in = 0.5
sigma_beam_out = 0.5; delta_sigma_beam_out = 0.3

# wounded nucleon/binary collision mixing ratio
alpha = 0.0


def get_density(
    eta, x, y, participant_trans_list, participant_eta_list, 
    binary_list, alpha):

    eps = 1e-8

    distance_trans = (
        (   (x - participant_trans_list[:, 0])**2. 
          + (y - participant_trans_list[:, 1])**2.)
        /(2.*participant_trans_list[:, 2]**2.)
    )

    idx = distance_trans < 25.
    dis_cut = distance_trans[idx]
    sigma_trans = participant_trans_list[idx, 2]
    sigma_eta = participant_eta_list[idx, 1:3]
    eta_0_cut = participant_eta_list[idx, 0]

    idx_left = eta_0_cut > eta
    idx_right = eta_0_cut <= eta
    dis_eta_left = (
        (eta - eta_0_cut[idx_left])**2./(2.*sigma_eta[idx_left, 0]**2.))
    dis_eta_right = (
        (eta - eta_0_cut[idx_right])**2./(2.*sigma_eta[idx_right, 1]**2.))

    rho_part = (
        sum(exp(-dis_cut[idx_left])/(2*pi*sigma_trans[idx_left]**2.)
            *exp(-dis_eta_left)/(sqrt(pi*sigma_eta[idx_left, 0]**2./2.) 
                                 + sqrt(pi*sigma_eta[idx_left, 1]**2./2.))) 
        + sum(exp(-dis_cut[idx_right])/(2.*pi*sigma_trans[idx_right]**2.)
              *exp(-dis_eta_right)/(sqrt(pi*sigma_eta[idx_right, 0]**2./2.)
                                   + sqrt(pi*sigma_eta[idx_right, 1]**2./2.)))
    )

    rho_binary = 0.0
    if abs(alpha) > eps: 
        for ibin in range(len(binary_list)):
            binary_x = binary_list[ibin, 0]
            binary_y = binary_list[ibin, 1]
            rho_binary += (
               (exp( - (eta - eta_0)**2./sigma_beam**2.) + 
                exp( - (eta + eta_0)**2./sigma_beam**2.))*0.5
               *exp( - ((x - participant_x)**2. + (y - participant_y)**2.)
                       /sigma_perp**2.)
            )
        rho_binary = prefactor_beam*prefactor_perp*rho_binary

    rho = rho_part*(1. - alpha)/2. + rho_binary*alpha
    return(rho)


def generate_3d_profile(data_path):
    random.seed()
    event_list = glob(path.join(data_path, 'ParticipantTable_event*.dat'))

    for iev in range(1, len(event_list)+1):

        participant_list = loadtxt(
            path.join(data_path, "ParticipantTable_event_%d.dat" % iev))

        participant_trans_list = zeros([len(participant_list[:, 0]), 3])
        participant_eta_list = zeros([len(participant_list[:, 0]), 3])
        participant_trans_list[:, 0:2] = participant_list[:, 0:2]
        for ipart in range(len(participant_list)):
            if rand_flag == 0:
                participant_trans_list[ipart, 2] = sigma_perp
            else:
                participant_trans_list[ipart, 2] = (
                    random.uniform(sigma_perp - delta_sigma_perp, 
                                   sigma_perp + delta_sigma_perp))
            if participant_list[ipart, 2] == 1:
                if rand_flag == 0:
                    participant_eta_list[ipart, 0] = eta_0
                    participant_eta_list[ipart, 1] = sigma_beam_in    # left
                    participant_eta_list[ipart, 2] = sigma_beam_out   # right
                else:
                    participant_eta_list[ipart, 0] = (
                        random.normal(eta_0, fluct_eta_0))
                    participant_eta_list[ipart, 1] = (               #left
                        random.uniform(sigma_beam_in - delta_sigma_beam_in, 
                                       sigma_beam_in + delta_sigma_beam_in))
                    participant_eta_list[ipart, 2] = (               #right
                        random.uniform(sigma_beam_out - delta_sigma_beam_out,
                                       sigma_beam_out + delta_sigma_beam_out)) 
            else:
                if rand_flag == 0:
                    participant_eta_list[ipart, 0] = -eta_0
                    participant_eta_list[ipart, 1] = sigma_beam_out   # left
                    participant_eta_list[ipart, 2] = sigma_beam_in    # right
                else:
                    participant_eta_list[ipart, 0] = (
                        random.normal(-eta_0, fluct_eta_0))
                    participant_eta_list[ipart, 1] = (               #left
                        random.uniform(sigma_beam_out - delta_sigma_beam_out, 
                                       sigma_beam_out + delta_sigma_beam_out))
                    participant_eta_list[ipart, 2] = (               #right
                        random.uniform(sigma_beam_in - delta_sigma_beam_in,
                                       sigma_beam_in + delta_sigma_beam_in)) 

        binary_list = loadtxt(
            path.join(data_path, "BinaryCollisionTable_event_%d.dat" % iev))

        entropy_density = zeros([grid_nrap, grid_nx, grid_ny])

        grid_eta = linspace( -(grid_nrap - 1.)/2.*grid_drap, 
                              (grid_nrap - 1.)/2.*grid_drap, grid_nrap)
        grid_x = linspace( -(grid_nx - 1.)/2.*grid_dx, 
                            (grid_nx - 1.)/2.*grid_dx, grid_nx)
        grid_y = linspace( -(grid_ny - 1.)/2.*grid_dy, 
                            (grid_ny - 1.)/2.*grid_dy, grid_ny)

        for ieta in range(len(grid_eta)):
            eta_local = grid_eta[ieta]
            print eta_local
            for ix in range(len(grid_x)):
                x_local = grid_x[ix]
                for iy in range(len(grid_y)):
                    y_local = grid_y[iy]

                    entropy_density[ieta, ix, iy] = get_density(
                        eta_local, x_local, y_local, 
                        participant_trans_list, participant_eta_list, 
                        binary_list, alpha)
         
        with file('sd_event_%d_block_3d.dat' % iev, 'w') as outfile:
            for slice_2d in entropy_density:
                savetxt(outfile, slice_2d)

def print_help_message():
    print "Usage : "
    print(color.bold
          + "./generateAvgprofile.py -ecm ecm "
          + "-cen cen_bounds"
          + "[-model model -collision_system collsys -cut_type cut_type]"
          + color.end)
    print "Usage of generateAvgprofile.py command line arguments: "
    print(color.bold + "-cen" + color.end
          + "   centrality bounds(%): "
          + color.purple + "20-30" + color.end)
    print(color.bold + "-ecm" + color.end
          + "   collision energy (GeV): "
          + color.purple + "7.7, 11.5, 19.6, 27, 39, 62.4, 200, 2760, 5500"
          + color.end)
    print(color.bold + "-cut_type" + color.end
          + "   centrality cut type: "
          + color.purple + color.bold + "total_entropy[default]" + color.end
          + color.purple + ", Npart" + color.end)
    print(color.bold + "-model" + color.end + " initial condition model: "
          + color.purple + color.bold + " MCGlb[default]" + color.end
          + color.purple + ", MCKLN" + color.end)
    print(color.bold + "-collision_system" + color.end
          + " type of collision system: "
          + color.purple + color.bold + " Pb+Pb[default]" + color.end
          + color.purple + ", Au+Au, Cu+Au, U+U, p+Pb, p+Au, d+Au, He+Au"
          + color.end)


if __name__ == "__main__":
    data_path = path.abspath(str(sys.argv[1]))

    generate_3d_profile(data_path)
