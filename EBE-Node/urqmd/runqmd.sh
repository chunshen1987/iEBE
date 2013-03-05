#!/bin/sh

export ftn09=uqmd.burner
export ftn10=OSCAR.input
export ftn13=particle_list.dat
export ftn14=particle_list.f14
export ftn15=particle_list.f15
export ftn16=particle_list.f16
export ftn19=particle_list.f19
export ftn20=particle_list.f20

time ./urqmd.$(uname)
