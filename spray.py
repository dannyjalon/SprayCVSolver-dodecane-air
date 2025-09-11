import numpy as np
from scipy.integrate import solve_ivp
import time

# IMPORTANT: This script now relies on a NON-INTERACTIVE version of cal1.py
from cal1 import cal1
from set_simulation_parameters import set_simulation_parameters
from post_shock_conditions import post_shock_conditions
from set_initial_conditions import set_initial_conditions
from state_vector_derivatives import state_vector_derivatives


def secant(B_guess, const_params):

    params = const_params.copy()
    params['B'] = B_guess

    shock = post_shock_conditions(params)
    y0 = set_initial_conditions(shock, params)

    t_span = [0, 0.1]
    sol = solve_ivp(
        fun=lambda x, y: state_vector_derivatives(x, y, params, shock),
        t_span=t_span,
        y0=y0,
        method='Radau',
        rtol=1e-6,
        atol=1e-6,
    )

    if sol.status!=0:
        return None

    # Calculate induction length based on maximum omega (heat release rate)
    omega_vals = (params['B'] * sol.y[4] * np.maximum(sol.y[6], 0) * np.maximum(sol.y[7], 0) *
                  np.exp(-params['Ea'] / (params['R_univ'] * sol.y[2])))

    if not np.any(omega_vals):
        return None

    idx_max_omega = np.argmax(omega_vals)
    li_pred = sol.t[idx_max_omega]
    return li_pred


def SprayModel(model_input_params):

    print("\n---Solving for B ---")
    li_target = model_input_params['Lind_target'] / 1000

    # Step 2: Set up constant parameters from the input dictionary
    # Use the non-interactive cal1 to get thermo properties
    out1, _ = cal1(
        DCJ_target=model_input_params['DCJ_target'],
        TvN_target=model_input_params['TvN_target'],
        TCJ_target=model_input_params['TCJ_target']
    )
    Ea = model_input_params['Ea']
    const_params = set_simulation_parameters(out1['gamma'], out1['Qf'], out1['cp_g'], Ea)

    # Step 3: Run the Secant Method Solver
    B0, B1 = 1.0e5, 1.0e7
    tolerance = 1e-3
    max_iterations = 25

    print(f"\nSolving for B to match li_target = {li_target:.6f} mm...")

    li0 = secant(B0, const_params)
    if li0 is None:
        print("\033[91mERROR: Simulation failed for initial guess B0. Aborting.\033[0m")
        return None
    li1 = secant(B1, const_params)
    if li1 is None:
        print("\033[91mERROR: Simulation failed for initial guess B1. Aborting.\033[0m")
        return None

    err0, err1 = li0 - li_target, li1 - li_target
    print(f"Initial B0={B0:.2e} -> li={li0:.6f} m (err={err0:.2e})")
    print(f"Initial B1={B1:.2e} -> li={li1:.6f} m (err={err1:.2e})")

    final_B = None
    for i in range(max_iterations):
        if abs(err1) < tolerance:
            print(f"\n\032[92mConvergence achieved after {i + 1} iterations!\033[0m")
            final_B = B1
            break

        denominator = err1 - err0
        if abs(denominator) < 1e-15:
            print("\033[91mERROR: Denominator is zero. Cannot continue.\033[0m")
            break

        B_next = B1 - err1 * (B1 - B0) / (denominator * 1000000)
        B0, B1 = B1, B_next
        err0 = err1

        li_next = secant(B_next, const_params)
        if li_next is None:
            print(f"\033[91mERROR: Simulation failed for B = {B_next:.4e}. Stopping.\033[0m")
            break

        err1 = li_next - li_target
        print(f"Iter {i + 1}: B={B_next:.4e} -> li={li_next:.6f} m (err={err1:.2e})")

    if final_B is None:
        print(f"\n\033[91mSolver did not converge within {max_iterations} iterations.\033[0m")
        final_B = B1

    print(f"\n--- Running Final Simulation with Optimal B = {final_B:.3e} ---")

    # Step 4: Run a final, detailed simulation
    params_final = const_params.copy()
    params_final['B'] = final_B
    shock_final = post_shock_conditions(params_final)
    y0_final = set_initial_conditions(shock_final, params_final)
    sol_final = solve_ivp(
        fun=lambda x, y: state_vector_derivatives(x, y, params_final, shock_final),
        t_span=[0, 0.1],
        y0=y0_final,
        method='Radau', rtol=1e-6, atol=1e-6, dense_output=True
    )

    if sol_final.status!=0:
        print((f"\n\033[97;41m {'ERROR: Final simulation failed even with the found B value!'} \033[0m\n"))
        return None

    # Step 5: Assemble and return the final predictions dictionary
    li_final = secant(final_B, const_params)
    predictions = {
        'DCJ_pred': f"{shock_final['Dcj']:.2f}",
        'TvN_pred': f"{shock_final['Tvn']:.2f}",
        'TCJ_pred': f"{sol_final.y[2, -1]:.2f}",
        'li_pred': f"{li_final:.6f}"
    }

    return predictions

