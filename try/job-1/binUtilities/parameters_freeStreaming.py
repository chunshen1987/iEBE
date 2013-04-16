#! /usr/bin/env python

# shared
time_interval = 0.1
starting_time = 0.62

# used only by hydro streaming
ed_dec = 0.18
ed_smallness = 0.1
hydro_endding_time = 22.0
hydro_extended_time = -1 # use fireball evolution

# used only by particle streaming
particle_endding_time = hydro_endding_time + 20
coordinate_xy_boundary = 13.0
coordinate_z_boundary = hydro_endding_time # speed of light assumption
