#!/usr/bin/env python

"""
    Thie module collects high level functions build on top of EBER and
    momCombinationR. It should contain functions that mostly deal with user
    interface than any actual analysis processing work.

    Also, all functions in this module should use as much as possible the
    default values set in other modules.
"""

from sys import argv

import EBER
import momCombinationR

amount_of_feedback = 99 # number of printed lines on the terminal

def collectMomDataFromSubfolder(folder_to_collect, particle_name_list = EBER.default_particle_names,
    allow_thermal_particles = True, allow_RD_particles = True,
    allow_ecc = True, form_combinations = True,
    use_historical_thermal_filename = False,

    # set up for initial eccentricity related files
    sn_file_path = EBER.default_sn_file_path,
    en_file_path = EBER.default_en_file_path,
    rp_sn_file_path = EBER.default_rp_sn_file_path, # read exotic ecc files from ebe folder
    save_rp_sn_as = EBER.default_save_rp_sn_as, # saved exotic ecc files
    save_sn_combination_as = momCombinationR.default_save_sn_combination_as, # saved ecc combinations
    rp_en_file_path = EBER.default_rp_en_file_path,
    save_rp_en_as = EBER.default_save_rp_en_as,
    save_en_combination_as = momCombinationR.default_save_en_combination_as,

    # set up for flow related files
    vn_file_path_abr = EBER.default_vn_file_path_abr,
    diff_vn_file_path_abr = EBER.default_diff_vn_file_path_abr,
    historical_thermal_vn_file_path_abr = EBER.default_historical_thermal_vn_file_path_abr,
    historical_thermal_diff_vn_file_path_abr = EBER.default_historical_thermal_diff_vn_file_path_abr,
    save_11P5N_as = EBER.default_save_as_abr,
    save_flow_combination_as = momCombinationR.default_save_flow_combination_as,

    echo_mode = amount_of_feedback):
    """
        Collect all relavant quantities from the folder "folder_to_collect".

        -- folder_to_collect: From which folder to collect data.

        -- particle_name_list: For what particles will the flow be collected and
            processed.

        -- allow_thermal_particles, allow_RD_particles: Whether to include
            thermal or RD particles in the process. Affects both the collecting
            and forming combination processes.

        -- allow_ecc: Whether to process for exotic eccentricities. Affects both
            the collecting and forming combination processes.

        -- form_combinations: Whether to form combinations from the collected
            moments. Affects both particles and eccentricites.

        -- use_historical_thermal_filename: When set to true,
            EBER.default_historical_thermal_(diff_)vn_file_path_abr will be used
            instead of default_(diff_)vn_file_path_abr for thermal particles.

        -- echo_mode: The larger this value is, the more output to the screen.

        -- others: For load/save files.
    """
    if echo_mode>8: print("\n===== collectMomDataFromSubfolder =====\n")

    particle_name_list = eval(str(particle_name_list))

    allow_thermal_particles = eval(str(allow_thermal_particles))
    allow_RD_particles = eval(str(allow_RD_particles))
    allow_ecc = eval(str(allow_ecc))
    form_combinations = eval(str(form_combinations))
    use_historical_thermal_filename = eval(str(use_historical_thermal_filename))

    # separate thermal and RD particle lists
    # will call EBER.collect11P5NColMomDataFromFolder twice for thermal and RD particles separately.
    thermal_particle_name_list = filter(EBER.particleNameIsThermal, particle_name_list)
    RD_particle_name_list = filter(lambda x: not EBER.particleNameIsThermal(x), particle_name_list)

    #--------------------------------------------------------------------------
    # Start collecting moments: 11P5N
    if (not use_historical_thermal_filename) and allow_thermal_particles and allow_RD_particles:
        # collect both thermal and RD particles using one call
        if echo_mode>9:
            print('\n---------------------------------------------------------------------\n')
            print('\n        Collecting thermal and RD particles                          \n')
            print('\n---------------------------------------------------------------------\n')
        EBER.collect11P5NColMomDataFromFolder(folder_to_collect, particle_name_list,
            save_11P5N_as, vn_file_path_abr, diff_vn_file_path_abr, sn_file_path, en_file_path,
            echo_mode = echo_mode)
    else:
        # separte the collection of thermal and RD particles
        # Thermal:
        if allow_thermal_particles:
            if use_historical_thermal_filename:
                vn_file_path_abr_to_use = historical_thermal_vn_file_path_abr
                diff_vn_file_path_abr_to_use = historical_thermal_diff_vn_file_path_abr
            else:
                vn_file_path_abr_to_use = vn_file_path_abr
                diff_vn_file_path_abr_to_use = diff_vn_file_path_abr
            if echo_mode>9:
                print('\n---------------------------------------------------------------------\n')
                print('\n        Collecting thermal particles                                 \n')
                print('\n---------------------------------------------------------------------\n')
            EBER.collect11P5NColMomDataFromFolder(folder_to_collect, thermal_particle_name_list,
                save_11P5N_as, vn_file_path_abr_to_use, diff_vn_file_path_abr_to_use,
                sn_file_path, en_file_path,
                echo_mode = echo_mode)
        # RD:
        if allow_RD_particles:
            if echo_mode>9:
                print('\n---------------------------------------------------------------------\n')
                print('\n        Collecting RD particles                                      \n')
                print('\n---------------------------------------------------------------------\n')
            EBER.collect11P5NColMomDataFromFolder(folder_to_collect, RD_particle_name_list,
                save_11P5N_as, vn_file_path_abr, diff_vn_file_path_abr,
                sn_file_path, en_file_path,
                echo_mode = echo_mode)

    #--------------------------------------------------------------------------
    # Then form combinations of flows
    if form_combinations:
        if allow_thermal_particles and allow_RD_particles:
            form_combination_using_particle_names = particle_name_list
        elif allow_thermal_particles:
            form_combination_using_particle_names = thermal_particle_name_list
        elif allow_RD_particles:
            form_combination_using_particle_names = RD_particle_name_list
        else:
            form_combination_using_particle_names = []

        if echo_mode>9:
            print('\n---------------------------------------------------------------------\n')
            print('\n        Forming combinations of flows (thermal)                      \n')
            print('\n---------------------------------------------------------------------\n')
        for particle_name in form_combination_using_particle_names:
            if echo_mode>10: print('-- For particle %s' % particle_name)
            momCombinationR.produceVnCombinationFrom11P5NColMomFiles(
                folder_to_collect, save_11P5N_as % particle_name,
                save_flow_combination_as % particle_name,
                echo_mode = echo_mode)

    #--------------------------------------------------------------------------
    # Then exotic eccentricity data
    if allow_ecc:
        # first collect
        if echo_mode>9:
            print('\n---------------------------------------------------------------------\n')
            print('\n        Collecting exotic eccentricity data                          \n')
            print('\n---------------------------------------------------------------------\n')
        EBER.collectExoticEccFromFolder(folder_to_collect,
            rp_sn_file_path, save_rp_sn_as,
            echo_mode = echo_mode)
        EBER.collectExoticEccFromFolder(folder_to_collect,
            rp_en_file_path, save_rp_en_as,
            echo_mode = echo_mode)

        # then form combinations
        if form_combinations:
            if echo_mode>9:
                print('\n---------------------------------------------------------------------\n')
                print('\n        Forming combinations of cplx eccentricities                  \n')
                print('\n---------------------------------------------------------------------\n')
            momCombinationR.produceEccCombinationFromExoticEccFiles(
                folder_to_collect, save_rp_sn_as, save_sn_combination_as,
                echo_mode = echo_mode)
            momCombinationR.produceEccCombinationFromExoticEccFiles(
                folder_to_collect, save_rp_en_as, save_en_combination_as,
                echo_mode = echo_mode)

    if echo_mode>8: print("collectMomDataFromSubfolder Done.\n")


#------------------------------------------------------------------------------
if __name__ == "__main__":
    print("----------------------------------------------------")
    print("Welcome! -- Zhi Qiu")
    print("----------------------------------------------------")
    if len(argv) == 1:
        print("Use one of the following:\n")
        print("collectMomDataFromSubfolder folder_path")
    else:
        print("Executing: "+argv[1]+"('"+"','".join(argv[2:])+"')\n")
        exec(argv[1]+"('"+"','".join(argv[2:])+"')")
