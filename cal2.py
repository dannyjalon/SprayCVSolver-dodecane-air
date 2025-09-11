import numpy as np


def cal2():
    """
    Calculates Activation Energy (Ea) based on induction times at two temperatures.
    """
    print("\n--- Activation Energy Calculation ---")
    tau1 = float(input("Enter tau1 (1.01U) [s] (default: 6.470e-5): ") or "6.470e-5")
    tau2 = float(input("Enter tau2 (0.99U) [s] (default: 6.549e-5): ") or "6.549e-5")
    T1 = float(input("Enter T1 [K] (default: 1723): ") or "1723")
    T2 = float(input("Enter T2 [K] (default: 1670): ") or "1670")

    R_univ = 8.314  # Universal gas constant

    Ea = R_univ * (np.log(tau2) - np.log(tau1)) / (1 / T2 - 1 / T1)

    out2 = {'Ea': Ea}

    print(f"\nActivation energy Ea (J/mol) = {Ea:.6f}")

    return out2