#!/usr/bin/env python

"""
    This module free streams particles produced mainly by the iSS code.
"""

from numpy import matrix, cos, sin, sqrt, floor # matrix is used as a string to array convertor
from formatter import getReorderingFunction

default_time_interval = 0.02
default_starting_time = 0.6
default_endding_time = 30
default_extended_time = 20

default_max_time_steps = None
default_coordinate_xy_boundary = 13.0
default_coordinate_z_boundary = None

smallness = 0.1

#################################################################
def freeStreaming2d(inStream, in_format_dict, outStream,
                  time_interval=default_time_interval,
                  starting_time=default_starting_time,
                  endding_time=default_endding_time,
                  max_time_steps=default_max_time_steps,
                  coordinate_xy_boundary=default_coordinate_xy_boundary,
                  t_sig="tau", x_sig="FZ_x", y_sig="FZ_y",
                  E_sig="E", pT_sig="pT", phi_sig="phi"):
    """ This function reads in the sample data stream for particles
    from inStream, then free streams the particles and write the trajectory
    of particles into ourStream.
    -- inStream, outStream: where to read and write (string streams).
    -- in_format_dict: dictionaries that gives the format of data streams
       (see assignmentFormat.assignmentExprStream2IndexDict).
    -- time_interval: by which time interval should particles be transported.
    -- starting_time, endding_time: if not None, only particle trajectory
       within time range [starting_time, endding_time] are outputted.
    -- max_time_steps: if not None, then each particles are transported by
       a maximum of this number of time steps.
    -- coordinate_xy_boundary: if not None, then only particle trajectory
       whose |x| or |y| values are outputted.
    -- xxx_sig: keys to be used with in_format_dict to identify which columns
       give the desired data.
    Besides writing data to outStream, this function also returns the
    format dictionary for output data stream.
    """
    # first read-in particle info
    extract_info_func = getReorderingFunction(in_format_dict, {t_sig:0, x_sig:1, y_sig:2, E_sig:3, pT_sig:4, phi_sig:5}, indexStartsWith1=False)
    particles = []; smallest_time_idx = 0;
    for aParticle in inStream:
        aParticle = matrix(aParticle).tolist()[0] # convert to numerical data

        # read particle info
        t, x, y, pT, E, phi = extract_info_func(aParticle)
        v = pT/E
        vx = v*cos(phi);
        vy = v*sin(phi);
        t_idx = floor((t-starting_time)/time_interval)
        prop_number = 0

        # add particle info to the global list
        particles.append([t_idx, prop_number, t,x,y,vx,vy]) # the 2nd element is used to store how many times this particle has been propergated
        if t_idx<smallest_time_idx: smallest_time_idx = t_idx

    # then simulate the free-streaming
    time_idx = smallest_time_idx
    while True:
        time = starting_time + time_idx*time_interval
        if time_idx>=0:
            # write to output stream

            # generate a "ghost particle" at the center
            outStream.write("%d %e  %e  %e  %e  %e\n" % (time_idx, time, 0, 0, 0, 0))

            # check time
            if time>=starting_time:
                # write other particles
                for particle_idx in range(len(particles)-1,-1,-1): # inverse order b/c some elements will be deleted during the loop
                    t_idx, prop_number, t,x,y,vx,vy = particles[particle_idx]

                    # check number of propergation on this particle
                    if max_time_steps!=None:
                        if prop_number>max_time_steps:
                            particles.pop(particle_idx)
                            continue

                    # check (x,y)
                    if coordinate_xy_boundary!=None:
                        if abs(x)>coordinate_xy_boundary or abs(y)>coordinate_xy_boundary:
                            particles.pop(particle_idx)
                            continue

                    # sync time
                    if t_idx>time_idx: continue

                    # output particle to stream
                    outStream.write("%d %e  %e  %e  %e  %e\n" % (t_idx, t,x,y,vx,vy))

        # propergate particles
        for particle_idx in range(len(particles)):
            t_idx, prop_number, t,x,y,vx,vy = particles[particle_idx]
            if t_idx>time_idx: continue # only those particles generated before this time will be propergated
            # free-streaming to next time
            t_idx += 1
            prop_number += 1
            t += time_interval
            x += time_interval*vx
            y += time_interval*vy
            particles[particle_idx] = [t_idx, prop_number, t,x,y,vx,vy]

        time_idx += 1
        print("Processing time step: " + str(time_idx))
        if not particles: break # empty particle list
        if time>endding_time: break
    print("Done.")
    return {"t_idx":1, "t":2, "x":3, "y":4, "v_x":5, "v_y":6}



#################################################################
def freeStreaming3d(inStream, in_format_dict, outStream,
                  time_interval=default_time_interval,
                  starting_time=default_starting_time,
                  endding_time=default_endding_time,
                  max_time_steps=default_max_time_steps,
                  coordinate_xy_boundary=default_coordinate_xy_boundary,
                  coordinate_z_boundary=default_coordinate_z_boundary,
                  t_sig="t", x_sig="FZ_x", y_sig="FZ_y", z_sig="z",
                  E_sig="E", pT_sig="pT", phi_sig="phi", pz_sig="p_z"):
    """ This function reads in the sample data stream for particles
    from inStream, then free streams the particles and write the trajectory
    of particles into ourStream.
    -- inStream, outStream: where to read and write (string streams).
    -- in_format_dict: dictionaries that gives the format of data streams
       (see assignmentFormat.assignmentExprStream2IndexDict).
    -- time_interval: by which time interval should particles be transported.
    -- starting_time, endding_time: if not None, only particle trajectory
       within time range [starting_time, endding_time] are outputted.
    -- max_time_steps: if not None, then each particles are transported by
       a maximum of this number of time steps.
    -- coordinate_xy_boundary: if not None, then only particle trajectory
       whose |x| or |y| values are smaller than this value are outputted.
    -- coordinate_z_boundary: if not None, then only particle trajectory
       whose |z| value is smaller than this value is outputted.
    -- xxx_sig: keys to be used with in_format_dict to identify which columns
       give the desired data.
    Besides writing data to outStream, this function also returns the
    format dictionary for output data stream.
    This is the 3d version.
    """
    # first read-in particle info
    extract_info_func = getReorderingFunction(in_format_dict, {t_sig:0, x_sig:1, y_sig:2, z_sig:3, E_sig:4, pT_sig:5, phi_sig:6, pz_sig:7}, indexStartsWith1=False)
    particles = []; smallest_time_idx = 0;
    for aParticle in inStream:
        aParticle = matrix(aParticle).tolist()[0] # convert to numerical data

        # read particle info
        t, x, y, z, E, pT, phi, pz = extract_info_func(aParticle)
        p = sqrt(pz*pz+pT*pT)
        v = p/E
        vT = pT/p*v
        vz = pz/p*v
        vx = vT*cos(phi);
        vy = vT*sin(phi);
        t_idx = floor((t-starting_time)/time_interval)
        prop_number = 0

        # add particle info to the global list
        particles.append([t_idx, prop_number, t,x,y,z,vx,vy,vz]) # the 2nd element is used to store how many times this particle has been propergated
        if t_idx<smallest_time_idx: smallest_time_idx = t_idx

    # then simulate the free-streaming
    time_idx = smallest_time_idx
    while True:
        time = starting_time + time_idx*time_interval
        if time_idx>=0:
            # write to output stream

            # generate a "ghost particle" at the center
            outStream.write("%d %e %e %e %e %e %e %e\n" % (time_idx, time, 0, 0, 0, 0, 0, 0))

            # check time
            if time>=starting_time:
                # write other particles
                for particle_idx in range(len(particles)-1,-1,-1): # inverse order b/c some elements will be deleted during the loop
                    t_idx, prop_number, t,x,y,z,vx,vy,vz = particles[particle_idx]

                    # check number of propergation on this particle
                    if max_time_steps!=None:
                        if prop_number>max_time_steps:
                            particles.pop(particle_idx)
                            continue

                    # check (x,y)
                    if coordinate_xy_boundary!=None:
                        if abs(x)>coordinate_xy_boundary or abs(y)>coordinate_xy_boundary:
                            particles.pop(particle_idx)
                            continue

                    # check z
                    if coordinate_z_boundary!=None:
                        if abs(z)>coordinate_z_boundary:
                            particles.pop(particle_idx)
                            continue

                    # sync time
                    if t_idx>time_idx: continue

                    # output particle to stream
                    outStream.write("%d %e %e %e %e %e %e %e\n" % (t_idx, t,x,y,z,vx,vy,vz))

        # propergate particles
        for particle_idx in range(len(particles)):
            t_idx, prop_number, t,x,y,z,vx,vy,vz = particles[particle_idx]
            if t_idx>time_idx: continue # only those particles generated before this time will be propergated
            # free-streaming to next time
            t_idx += 1
            prop_number += 1
            t += time_interval
            x += time_interval*vx
            y += time_interval*vy
            z += time_interval*vz
            #if abs(z)>abs(t):
            #    print(t_idx, prop_number, t,x,y,z,vx,vy,vz)
            particles[particle_idx] = [t_idx, prop_number, t,x,y,z,vx,vy,vz]

        time_idx += 1
        print("Processing time step: " + str(time_idx))
        if not particles: break # empty particle list
        if time>endding_time: break
    print("Done.")
    return {"t_idx":1, "t":2, "x":3, "y":4, "z":5, "v_x":6, "v_y":7, "v_z":8}



#################################################################
# OBSOLETED --- DO NOT USE
def streamingHydroZ_old(inStream, in_format_dict, outStream,
                  time_interval=default_time_interval,
                  starting_time=default_starting_time,
                  endding_time=default_endding_time,
                  extended_evolusion=True,
                  t_sig="tau", x_sig="x", y_sig="y",
                  ed_sig="e",
                  use_compact_output=True, ed_dec=0.18, ed_smallness=smallness,
                  slice1_stream=None, slice2_stream=None,
                  ):
    """ This function reads in the sample data stream for hydro evolution
    in the plane z=0, then extends it to z-direction using free streaming
    in z direction. The evolution is between starting_time and endding_time.
    -- inStream, outStream: where to read and write (string streams).
    -- in_format_dict: dictionaries that gives the format of data streams
       (see assignmentFormat.assignmentExprStream2IndexDict).
    -- time_interval: by which time interval should particles be transported.
    -- starting_time, endding_time: if not None, only particle trajectory
       within time range [starting_time, endding_time] are outputted.
    -- xxx_sig: keys to be used with in_format_dict to identify which columns
       give the desired data.
    -- use_compact_output: when set to True, only those cells in the "leading"
       slices or those have ed close to ed_dec (within smallness) will be
       outputted. The two slices will be outputted into slice1_stream and
       slice2_stream separately, and the rest of freeze-out surface will be
       outputted into the outStream.
    Besides writing data to outStream, this function also returns the
    format dictionary for output data stream.
    Note that a more efficient way would be to store those readed cells in
    memory then write them out without re-do the checking etc. (FUTURE)
    """

    # convert in stream to desired data
    extract_info_func = getReorderingFunction(in_format_dict, {t_sig:0, x_sig:1, y_sig:2, ed_sig:3}, indexStartsWith1=False)

    # start to build the evolution
    max_time_index = int(floor((endding_time-starting_time)/time_interval))
    for time_idx in range(max_time_index):
        time = starting_time + time_idx*time_interval
        inStream.seek(0) # reset the pointer to the 1st character
        for aLine in inStream:
            aLine = matrix(aLine).tolist()[0] # convert to numerical data
            t, x, y, ed = extract_info_func(aLine) # get info from hydro file
            if t>time: break # assume that the data are stored in forward time order, data beyond this line should be ignored by causality
            z = starting_time + (time-t) # z in the forward direction
            # output to stream
            if use_compact_output: # skip most of the output
                if abs(t-starting_time)<1e-10:
                    # leading slices
                    slice1_stream.write("%e  %e  %e  %e  %e\n" % (time,x,y,z,ed))
                    slice2_stream.write("%e  %e  %e  %e  %e\n" % (time,x,y,-z,ed))
                    continue
                if abs(ed-ed_dec)<ed_smallness:
                    # cell close to freeze-out either
                    outStream.write("%e  %e  %e  %e  %e\n" % (time,x,y,z,ed))
                    outStream.write("%e  %e  %e  %e  %e\n" % (time,x,y,-z,ed))
                    continue
            else:
                # forward direction
                outStream.write("%e  %e  %e  %e  %e\n" % (time,x,y,z,ed))
                # backward direction
                if abs(z)>1e-10: outStream.write("%e  %e  %e  %e  %e\n" % (time,x,y,-z,ed))
                continue

        print("Processing time step: " + str(time_idx))
    print("Done.")
    return {"t":1, "x":2, "y":3, "z":4, "e":5}



#################################################################
def streamingHydroZ(inStream, in_format_dict, outStream,
                  time_interval=default_time_interval,
                  starting_time=default_starting_time,
                  endding_time=default_endding_time,
                  extended_time=-1,
                  t_sig="tau", x_sig="x", y_sig="y",
                  ed_sig="e", ed_dec=0.18,
                  ed_smallness=smallness,
                  ):
    """ This function reads in the sample data stream for hydro evolution
    in the plane z=0, then extends it to z-direction using free streaming
    in z direction. The evolution is between starting_time and endding_time.
    -- inStream, outStream: where to read and write (string streams).
    -- in_format_dict: dictionaries that gives the format of data streams
       (see assignmentFormat.assignmentExprStream2IndexDict).
    -- time_interval: by which time interval should particles be transported.
    -- starting_time, endding_time: if not None, only particle trajectory
       within time range [starting_time, endding_time] are outputted.
    -- extended_time: whether to extend time evolution. The inner structure
       will be will as the "leading" slices during this extended evolution.
       If set to -1, the evolution time of the fireball will be used.
    -- xxx_sig: keys to be used with in_format_dict to identify which columns
       give the desired data.
    Only those cells in the "leading" slices or those have ed close to
    ed_dec (within smallness) will be outputted.
    Besides writing data to outStream, this function also returns the
    format dictionary for output data stream.
    This version stores those readed cells in memory then write them out
    without re-do the checking etc.
    """

    # convert in stream to desired data
    extract_info_func = getReorderingFunction(in_format_dict, {t_sig:0, x_sig:1, y_sig:2, ed_sig:3}, indexStartsWith1=False)

    # first read-in all relavant cells to memory, also store the lastest time
    stored_cells = [] # store those already registered cells (only in + direction)
    flag_first_line = True
    print("Read hydro movie and extract freeze-out surface info...")
    counter = 1
    for aLine in inStream:
        aLine = matrix(aLine).tolist()[0] # convert to numerical data
        t, x, y, ed = extract_info_func(aLine) # get info from hydro file
        if flag_first_line:
            earliest_time=t
            latest_time=t
            flag_first_line=False
        if t>latest_time: latest_time=t
        # output to stream
        if t<=earliest_time:
            # leading slices
            stored_cells.append([t,x,y,ed])
            continue
        if abs(ed-ed_dec)<ed_smallness:
            # cell close to freeze-out either
            stored_cells.append([t,x,y,ed])
            continue
        if (counter % 100000)==0: print("Line %d reached..." % counter)
        counter += 1
    print("Total number of cells readed: " + str(len(stored_cells)) )

    # start to build the evolution
    max_time_index = int(floor((endding_time-starting_time)/time_interval))
    for time_idx in range(max_time_index):
        time = starting_time + time_idx*time_interval
        for tau,x,y,ed in stored_cells:
            if tau>time: break # assume that the data are stored in forward time order, data beyond this line should be ignored by causality
            z = sqrt(time*time - tau*tau) # z in the forward direction; time: where z should be if it has speed of light; (t-earliest_time): delay
            outStream.write("%e  %e  %e  %e  %e\n" % (time,x,y,z,ed))
            outStream.write("%e  %e  %e  %e  %e\n" % (time,x,y,-z,ed))
        print("Processing time step: " + str(time_idx))

    # Extend the evolusion until all the fireball is out-of the view point (z>endding_time) the leading-slices are replaced by the inner structure of the fireball at the corresponding time
    shift_time = time-earliest_time # all the time in the inStream will be shifted by this value; this is the "current time" before the "post process"
    if extended_time<1: extended_time = latest_time # those slices with time larger than this value will be ignored
    fixed_z = sqrt(time*time-earliest_time*earliest_time) # all the slices will be fixed at this z value
    inStream.seek(0) # reset cell index in inStream
    last_time = -1 # will trigger the variable flag_time_change to set to true, which will make the code to scan cells in the stored array and output them into file. The flag_time_change variable is usually triggered when the time read from inStream is changed.
    search_stored_cell_from_idx = 0 # the stored array is searched starting from this index
    search_stored_cell_until_idx = len(stored_cells) # last index
    post_process_time_idx = 0 # keep track of steps
    for aLine in inStream:
        aLine = matrix(aLine).tolist()[0] # convert to numerical data
        tau, x, y, ed = extract_info_func(aLine) # get info from hydro file
        current_time = tau + shift_time # the current time (time is evolved while reading the inStream)
        if current_time > shift_time + extended_time: break # ignore the rest
        if current_time > last_time: flag_time_change = True
        if flag_time_change:
            flag_time_change = False
            last_time = current_time
            # scan and output cells from the stored array
            for idx in range(search_stored_cell_from_idx, search_stored_cell_until_idx):
                s_t,s_x,s_y,s_ed = stored_cells[idx] # "s" for stored
                s_z = sqrt(current_time*current_time - s_t*s_t) # the would-be z location: slice_t is the actual elapsed time, (s_t-earliest_time) is the decay due to the fact this cell is not the earliest cell
                if s_z >= fixed_z:
                    search_stored_cell_from_idx = idx+1 # store the last skipped index + 1: the next search starts here
                    continue # such cells are out-of-boundary
                outStream.write("%e  %e  %e  %e  %e\n" % (current_time,s_x,s_y,s_z,s_ed))
                outStream.write("%e  %e  %e  %e  %e\n" % (current_time,s_x,s_y,-s_z,s_ed))
            print("Processing post process time step: " + str(post_process_time_idx))
            post_process_time_idx += 1
        # write out leading slices
        outStream.write("%e  %e  %e  %e  %e\n" % (current_time,x,y,fixed_z,ed))
        outStream.write("%e  %e  %e  %e  %e\n" % (current_time,x,y,-fixed_z,ed))


    print("Done.")
    return {"t":1, "x":2, "y":3, "z":4, "e":5}



#################################################################
def freeStream2dAlongZ(inStream, outStream, starting_z, velocity_z, evolution_time,
                       time_interval = default_time_interval,
                       ):
    """
    This function takes particle from inStream in the format "x y", place them
    at constant z=starting_z and give them together with speed velocity_z.
    This function then calls freeStreaming3d function to generate the actual
    evolution. The evolution is done in terms of time_interval until
    evolution_time.
    This function returns the format of the output.
    """
    format_dict = {"t":0, "FZ_x":1, "FZ_y":2, "z":3, "E":4, "pT":5, "phi":6, "p_z":7}
    def dataStream():
        for coord in inStream:
            x,y = matrix(coord).tolist()[0] # convert to numerical data
            yield "0 %e %e %e 1 0 0 %e\n" % (x,y,starting_z,velocity_z)
    return freeStreaming3d(dataStream(),format_dict,outStream,
                    time_interval=time_interval,
                    starting_time=0,
                    endding_time=evolution_time,
                    max_time_steps=None,
                    coordinate_xy_boundary=None,
                    coordinate_z_boundary=None)
