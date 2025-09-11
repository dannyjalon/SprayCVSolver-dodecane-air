import numpy as np
from scipy.optimize import least_squares


def cal1(DCJ_target, TvN_target, TCJ_target, alpha0=0.06667, T0=298, Lv=260628, cp_d=2176):
    """
    Solves for thermodynamic properties [gamma, Qf, cp_g] based on target values.
    This version receives parameters as arguments instead of user input.
    """
    x0 = np.array([1.30, 4.4e7, 1000])  # Initial guesses: [gamma, Qf, cp_g]

    def residuals(z, DCJ_target, alpha0, T0, TCJ, TvN, Lv, cp_d):
        gamma, Qf, cp_g = z
        Qnet = Qf - (Lv - cp_d * T0)
        q = alpha0 * (gamma + 1) * Qnet / (cp_g * T0)

        sqrt_arg = q ** 2 + 2 * q
        if sqrt_arg < 0:
            sqrt_arg = 0

        MCJ = 1 / np.sqrt(1 + alpha0) * np.sqrt(1 + q + np.sqrt(sqrt_arg))
        a0 = np.sqrt((gamma - 1) * cp_g * T0)
        DCJ_calc = MCJ * a0
        TCJ2 = T0 * ((1 + (alpha0 + 1) * gamma * MCJ ** 2) / ((1 + gamma) * (alpha0 + 1) * MCJ)) ** 2
        TvN3 = T0 * (1 + (2 * (gamma - 1) / (gamma + 1) ** 2) * (gamma * (MCJ ** 2 - 1) + 1 - 1 / MCJ ** 2))
        return [DCJ_calc - DCJ_target, TCJ2 - TCJ, TvN3 - TvN]

    lb = [1.1, 1.0e5, 500]
    ub = [1.6, 9.0e8, 2000]
    result = least_squares(residuals, x0, args=(DCJ_target, alpha0, T0, TCJ_target, TvN_target, Lv, cp_d), bounds=(lb,
                                                                                                                   ub))

    gamma, Qf, cp_g = result.x

    Qnet = Qf - (Lv - cp_d * T0)
    q = alpha0 * (gamma + 1) * Qnet / (cp_g * T0)
    MCJ_calc = 1 / np.sqrt(1 + alpha0) * np.sqrt(1 + q + np.sqrt(q ** 2 + 2 * q))
    a0 = np.sqrt((gamma - 1) * cp_g * T0)
    DCJ_calc = MCJ_calc * a0
    TCJ_calc = T0 * ((1 + (alpha0 + 1) * gamma * MCJ_calc ** 2) / ((1 + gamma) * (alpha0 + 1) * MCJ_calc)) ** 2
    TvN_calc = T0 * (1 + (2 * (gamma - 1) / (gamma + 1) ** 2) * (gamma * (MCJ_calc ** 2 - 1) + 1 - 1 / MCJ_calc ** 2))

    out1 = {'gamma': gamma, 'Qf': Qf, 'cp_g': cp_g}
    calc_results = {'DCJ_calc': DCJ_calc, 'TvN_calc': TvN_calc, 'TCJ_calc': TCJ_calc}

    return out1, calc_results
