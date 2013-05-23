controlParameterList = {
    'simulation_type'       :   'hybrid', # 'hybrid' or 'hydro'
}

superMCParameters = {
    'which_mc_model'                :   5,
    'sub_model'                     :   1,
    'Npmin'                         :   0,
    'Npmax'                         :   1000,
    'bmin'                          :   13,
    'bmax'                          :   20,
    'ecm'                           :   2760,
    'finalFactor'                   :   56.763,
    'alpha'                         :   0.118,
    'lambda'                        :   0.288,
}

hydroParameters = {
    'vis'       :   0.08,
    'T0'        :   0.6, # tau_0
    'Edec'      :   0.18,
}

iSSParameters = {
    'number_of_repeated_sampling'   :   10,
    'y_LB'                          :   -2.5,
    'y_RB'                          :   2.5,
}
