import numpy as np


def post_shock_conditions(params):
    """
    Computes post-shock conditions based on simulation parameters.
    """
    gamma = params['gamma']
    T1 = params['T1']
    P1 = params['P1']
    R = params['R']
    rd0 = params['rd0']
    rhod = params['rhod']
    Qmix = params['Qmix']

    # Pre-shock conditions
    rho0 = P1 / (R * T1)
    c0 = np.sqrt(gamma * R * T1)

    # Detonation velocity calculation
    alpha1 = 1 / params['s']
    Qmix_s = (gamma + 1) * alpha1 / (2 * params['cpg'] * T1) * (params['Qf'] - params['latentv'] + params['cpd'] * T1)
    M0 = np.sqrt((Qmix_s + 1) / (1 + alpha1)) + np.sqrt(Qmix_s / (1 + alpha1))
    Dcj = M0 * c0

    # Von-Neumann conditions
    shock = {}
    shock['Mvn'] = np.sqrt((1 + 0.5 * (gamma - 1) * M0 ** 2) / (gamma * M0 ** 2 - 0.5 * (gamma - 1)))
    shock['Pvn'] = P1 * (1 + (2 * gamma / (gamma + 1)) * (M0 ** 2 - 1))
    shock['Tvn'] = T1 * (1 + 2 * (gamma - 1) / (gamma + 1) ** 2 * (gamma * (M0 ** 2 - 1) + 1 - 1 / M0 ** 2))
    shock['rhovn'] = rho0 * (M0 ** 2 * (gamma + 1)) / (2 + (gamma - 1) * M0 ** 2)
    u0 = Dcj
    shock['uvn'] = u0 * rho0 / shock['rhovn']
    shock['Dcj'] = Dcj

    #print(f"\nD_CJ: {Dcj:.6f}\nM_CJ: {M0:.6f}\nM_vn: {shock['Mvn']:.6f}")
    #print(f"T_vn: {shock['Tvn']:.6f}\nP_vn: {shock['Pvn']:.6f}")
    #print(f"rho_vn: {shock['rhovn']:.6f}\nu_vn: {shock['uvn']:.6f}")

    # Post-shock species fractions and droplet number density
    shock['YF0'] = 0
    shock['YO0'] = 1
    alpha0 = params['phi'] / (params['s'] * shock['rhovn'] / rho0)
    shock['nd0'] = alpha0 / (4 / 3 * np.pi * rd0 ** 3 * rhod / shock['rhovn'])

    return shock