controlParameterList = {
    'simulation_type'       :   'hydroEM', # 'hybrid' or 'hydro'
}

superMCParameters = {
    'which_mc_model'                :   5,
    'sub_model'                     :   1,
    'Npmin'                         :   223,
    'Npmax'                         :   417,
    'bmin'                          :   0.0,
    'bmax'                          :   7.7,
    'ecm'                           :   2760,
    'finalFactor'                   :   56.763,
    'alpha'                         :   0.118,
    'lambda'                        :   0.288,
}

hydroParameters = {
    'vis'       :   0.08,
    'T0'        :   0.6, # tau_0
    'Edec'      :   0.18,
    'IhydroJetoutput' :   1,  # switch for output hydro evolution history into hdf5 file
}

iSSParameters = {
    'calculate_vn'                  :   1,
    'MC_sampling'                   :   0,
    'number_of_repeated_sampling'   :   10,
    'y_LB'                          :   -2.5,
    'y_RB'                          :   2.5,
}

photonEmissionParameters = {
    'dx'          :   0.5,
    'dy'          :   0.5,
    'dTau'        :   0.1,
    'calHGIdFlag' :   0,
}
