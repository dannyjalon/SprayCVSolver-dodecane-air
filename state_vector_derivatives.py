import numpy as np
from dragcoefficient import drag_coefficient


def get_vapor_pressure(Td):
    """
    Compute vapor pressure using an empirical correlation from NIST (bar).
    """
    A = 4.10549
    B = 1625.928
    C = -92.839
    # returns pressure in Pascals
    return (10 ** (A - B / (Td + C))) * 100000


def state_vector_derivatives(x, y, params, shock):
    """
    Evaluates the derivatives of the state vector.
    State vector y is ordered as: [ug, ud, Tg, Td, rhog, rd, YF, YO]
    """
    # Unpack state vector
    ug, ud, Tg, Td, rhog, rd, YF, YO = y

    # Unpack parameters
    gamma = params['gamma']
    Lg = params['latentv']
    cpg = params['cpg']
    cpd = params['cpd']
    k = params['k']
    mu = params['mu']
    R = params['R']
    Ea = params['Ea']
    Qf = params['Qf']
    Wf = params['Wf']
    Wmix = params['Wmix']
    Le = params['Le']
    B = params['B']
    s = params['s']
    Ru = params['R_univ']
    rd0 = params['rd0']
    rhod = params['rhod']
    nd0 = shock['nd0']

    # Compute gas pressure and droplet number density
    Pg = rhog * R * Tg
    nd = shock['nd0'] * shock['Dcj'] / ud

    # --- Reaction rate ---
    # Ensure mass fractions are non-negative before calculating omega
    YF_safe = max(YF, 0)
    YO_safe = max(YO, 0)
    omega = B * rhog * YF_safe * YO_safe * np.exp(-Ea / (Ru * Tg))

    # --- Check if droplet has evaporated ---
    if (rd / rd0 > 1e-2) and (nd / nd0 > 1e-4):
        Pfs = get_vapor_pressure(Td)  # Vapor pressure (Pa)
        Xfs = Pfs / Pg  # Molar fraction on droplet surface
        YFs = Xfs * Wf / (Xfs * Wf + (1 - Xfs) * Wmix)  # Vapor fuel mass fraction on surface

        # Avoid log of non-positive number
        log_arg = (1 - YF) / (1 - YFs)
        if log_arg <= 0:
            lamda = 0
        else:
            lamda = (1 / Le) * np.log(log_arg)

        # Vaporization rate
        md = 4 * np.pi * rd * k / cpg * lamda

        # Heat transfer rate
        qd = 4 * np.pi * k * rd * ((Tg - Td) / (np.exp(lamda) - 1) - Lg / cpg) * lamda

        # Drag force
        fx = 6 * np.pi * mu * rd * (ug - ud)

        # Liquid-phase evolution derivatives
        Dt = k / (rhog * cpg)  # Gas thermal diffusivity
        drddx = -((rhog * Dt) / (rhod * rd * ud)) * lamda
        duddx = (9 / 2 * (mu) / (rhod * rd ** 2 * ud)) * (ug - ud)
        dTddx = (3 * (rhog * cpg * Dt) / (rhod * cpd * rd ** 2 * ud)) * qd / (4 * np.pi * k * rd)
    else:
        # Droplet evaporated
        lamda = 0.0
        md = 0.0
        qd = 0.0
        fx = 0.0
        drddx = 0.0
        duddx = 0.0
        dTddx = 0.0

    # Gas-vapor evolution
    if (rd / rd0 > 1e-2) and (nd / nd0 > 1e-4):
        E = (nd * md * ud - nd * fx) * (ud - ug) - nd * (md * (Lg - cpg * Tg) + qd) + Qf * omega
    else:
        E = Qf * omega

    # Gas phase derivatives
    denominator = (rhog * ug - gamma * Pg / ug)
    if abs(denominator) < 1e-9:  # Avoid division by zero, especially at stagnation points
        denominator = 1e-9

    dugdx = (1 / denominator) * ((1 - gamma) * E / ug + nd * md * (ud - ug) - nd * fx)
    dTgdx = gamma / (cpg * rhog * ug) * (E - Pg * dugdx) - Tg / (ug * rhog) * nd * md
    drhogdx = nd * md / ug - rhog / ug * dugdx
    dYFdx = (nd * md - omega) / (rhog * ug) - YF / rhog * drhogdx - YF / ug * dugdx
    dYOdx = -s * omega / (rhog * ug) - YO / rhog * drhogdx - YO / ug * dugdx

    return np.array([dugdx, duddx, dTgdx, dTddx, drhogdx, drddx, dYFdx, dYOdx])