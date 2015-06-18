controlParameterList = {
    'simulation_type'       :   'hydroEM_with_decaycocktail', # 'hybrid', 'hydro', 'hydroEM', 'hydroEM_with_decaycocktail', 'hydroEM_preEquilibrium'
    'niceness'              :   0,       # range from 0 to 19 for process priority, 0 for the highest priority
}

centralityParameters = {
    'centrality': '0-10%',  # centrality bin
    'cut_type': 'total_entropy',
    # centrality cut variable: total_entropy or Npart
}

superMCParameters = {
    'model_name'                    :   'MCGlb',    # MCGlb or MCKLN
    'Aproj'                         :   208,
    'Atarg'                         :   208,
    'ecm'                           :   2760,
    'finalFactor'                   :   56.763,
    'alpha'                         :   0.118,      # WN/BC mixing ratio in MCGlb
    'lambda'                        :   0.218,      # saturation scale parameter in MCKLN
    'operation'                     :   1,
    'include_NN_correlation'        :   0,
    'cc_fluctuation_model'          :   6,
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
    'visflag'   :   0,        # flag to use temperature dependent eta/s(T)
    'T0'        :   0.6,      # tau_0
    'Edec'      :   0.18,
    'IhydroJetoutput' :   1,  # switch for output hydro evolution history into hdf5 file
    'InitialURead'    :   0,  # set it to be 1 when simulation_type == hydroEM_preEquilibrium
    'Initialpitensor' :   0,  # initialization of pi tensor
}

iSSParameters = {
    'calculate_vn'                  :   1,
    'MC_sampling'                   :   0,
    'number_of_repeated_sampling'   :   10,
    'y_LB'                          :   -2.5,
    'y_RB'                          :   2.5,
    'sample_y_minus_eta_s_range'    :   2.0,
}

photonEmissionParameters = {
    'dx'          :   0.1,
    'dy'          :   0.1,
    'dTau'        :   0.05,
    'tau_start'   :   0.6,
    'tau_end'     :   20.0,
    'T_dec'       :   0.12,
    'T_cuthigh'   :   0.8,
    'T_cutlow'    :   0.1,
    'T_sw_high'   :   0.180,
    'T_sw_low'    :   0.1795,
    'calHGIdFlag' :   0,
    'differential_flag'   :  2,
    'enable_polyakov_suppression'   :    0,
}
