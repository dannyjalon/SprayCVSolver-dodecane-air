def set_simulation_parameters(gamma, Qf, cpg, Ea):
    """
    Defines the dictionary of simulation parameters.
    """
    params = {}

    params['gamma'] = gamma
    params['Qf'] = Qf
    params['phi'] = 0.956
    params['s'] = 15.0 / params['phi']
    params['Qmix'] = params['Qf'] / (1 + params['s'])
    params['latentv'] = 260628
    params['B'] = 5.0e6
    params['Ea'] = Ea
    params['rd0'] = 5.0e-6
    params['T1'] = 298
    params['P1'] = 100000
    params['Wf'] = 0.170
    params['rhod'] = 750
    params['Tb'] = 489
    params['Wmix'] = 0.0304
    params['cpg'] = cpg
    params['R'] = cpg * (gamma - 1) / gamma
    params['cpd'] = 2176
    params['k'] = 0.08
    params['mu'] = 8e-4
    params['Le'] = 1.0
    params['Pr'] = cpg * params['mu'] / params['k']
    params['Wo'] = 0.029
    params['R_univ'] = 8.314

    return params