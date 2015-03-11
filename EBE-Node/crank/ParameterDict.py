controlParameterList = {
    'simulation_type'       :   'hydroEM_with_decaycocktail', # 'hybrid', 'hydro', 'hydroEM', 'hydroEM_with_decaycocktail', 'hydroEM_preEquilibrium'
    'niceness'              :   0,       # range from 0 to 19 for process priority, 0 for the highest priority
}

centralityParameters = {
    'centrality': '20-30%',  # centrality bin
    'cut_type': 'total_entropy',
    # centrality cut variable: total_entropy or Npart
}

superMCParameters = {
    'which_mc_model'                :   5,
    'sub_model'                     :   1,
    'Aproj'                         :   208,
    'Atarg'                         :   208,
    'ecm'                           :   2760,
    'finalFactor'                   :   56.763,
    'alpha'                         :   0.118,
    'lambda'                        :   0.218,
    'operation'                     :   1,
}

# only effective when simulation_type == hydroEM_preEquilibrium
preEquilibriumParameters = {
    'event_mode'            :    1,  
    'taumin'                :    0.6,
    'taumax'                :    0.6,
    'dtau'                  :    0.2,
}

hydroParameters = {
    'vis'       :   0.08,
    'T0'        :   0.6,      # tau_0
    'Edec'      :   0.18,
    'IhydroJetoutput' :   1,  # switch for output hydro evolution history into hdf5 file
    'InitialURead'    :   0,  # set it to be 1 when simulation_type == hydroEM_preEquilibrium
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
    'T_sw_high'   :   0.170,
    'T_sw_low'    :   0.150,
    'calHGIdFlag' :   0,
}
