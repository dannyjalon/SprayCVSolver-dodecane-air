import numpy as np

def set_initial_conditions(shock, params):
    """
    Sets the initial conditions vector for the ODE solver.

    Args:
        shock (dict): Dictionary with post-shock conditions.
        params (dict): Dictionary with simulation parameters.

    Returns:
        np.array: Initial state vector [ug, ud, Tg, Td, rhog, rd, YF, YO]
    """
    ug0 = shock['uvn']
    ud0 = shock['Dcj']
    Tg0 = shock['Tvn']
    Td0 = params['T1']
    rhog0 = shock['rhovn']
    rd0 = params['rd0']
    y0 = np.array([ug0, ud0, Tg0, Td0, rhog0, rd0, shock['YF0'], shock['YO0']])
    return y0