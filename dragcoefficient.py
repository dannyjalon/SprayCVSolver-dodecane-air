import numpy as np


def drag_coefficient(rhog, ug, ud, mu, rd, c):
    """
    Calculates the drag coefficient for a spherical particle, considering
    compressibility and rarefaction effects.
    Based on "Compressibility and Rarefaction Effects on Drag of a Spherical Particle" by Loth.

    Args:
        rhog (float): Gas density
        ug (float): Gas velocity
        ud (float): Droplet velocity
        mu (float): Dynamic viscosity
        rd (float): Droplet radius
        c (float): Speed of sound

    Returns:
        float: Drag coefficient (Cd)
    """
    # Reynolds Number
    Re = 2 * rd * abs(ud - ug) * rhog / mu

    # Relative Mach number
    Ma = abs(ud - ug) / c

    # Empirical correlation for Gm
    if Ma < 0.89:
        Gm = 1 - 1.525 * Ma ** 4
    else:
        Gm = 1e-4 * (2 + 8 * np.tanh(12.77 * (Ma - 2.02)))

    # Drag ratio for Cm
    if Ma < 1.45:
        Cm = 1 / 3 * (5 + 2 * np.tanh(3 * np.log(Ma + 0.1)))
    else:
        Cm = 2.044 + 0.2 * np.exp(-1.8 * np.log(Ma / 1.5) ** 2)

    # Empirical correlation for Hm
    Hm = 1 - 0.258 * Cm / (1 + 514 * Gm)

    # Drag coefficient calculation
    if Re < 0.1:
        Cd = 24 / Re
    elif Re < 45:
        Cd = (24 / Re) * (1 + 0.15 * Re ** 0.687)
    else:
        Cd = Hm * (24 / Re) * (1 + 0.15 * Re ** 0.687) + 0.42 * Cm / (1 + 42500 * Gm * Re ** (-1.16))

    return Cd