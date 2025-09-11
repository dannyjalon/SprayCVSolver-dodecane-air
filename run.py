import os
import numpy as np
from spray import SprayModel
from soot_foil_image_tool import measure_image

# Spray detonation modeling
input_params = {
    'DCJ_target': 1805.0,   # Target CJ velocity [m/s]
    'TvN_target': 1689.0,   # Target von Neumann temperature [K]
    'TCJ_target': 2875.0,   # Target CJ temperature [K]
    'Ea': 5546,             # Activation energy [J/mol]
    'Lind_target': 20.096,   # Target induction length [mm]
    'cell_size_target': 17  # Target Cell Size [cm]
}

print("Starting simulation with the following parameters:")
print(input_params)

try:
    # The function will run the simulation and return the results
    results = SprayModel(input_params)
    if not results:
        print(f"\n\033[97;41m ERROR: Simulation failed! Aborting ... \033[0m")
        exit()

    print((f"\n\033[97;42m {'Final simulation successful! Displaying results...'} \033[0m\n"))
    print(results)

    # Soot foil analysis
    # Inputs

    image_path = '/Users/dnyylwnzqy/Desktop/My Files/MSc/Thesis/my files/CV/more images/spray1.png'
    dimension = 'h'
    size_cm = 50
    step = 10

    print("\n--- Soot Foil Image Analysis ---")
    print("\nAnalyzing image... (this may take a moment)")
    try:
        annotated_image, stats = measure_image(image_path, dimension, size_cm, step=step)
        # The 'l_abs_values_converted' key corresponds to the transverse cell width measurements
        if dimension == 'h':
            mean_cell_size = np.mean(stats['l_euclidean_values_converted'])
            print("Image analysis complete.")
            print(f"  > Mean Cell Size: {mean_cell_size:.4f} cm")
        if dimension == 'w':
            mean_cell_size = np.mean(stats['w_euclidean_values_converted'])
            print("Image analysis complete.")
            print(f"  > Mean Cell Size: {mean_cell_size:.4f} cm")

    except Exception as e:
        print(f"\nAn error occurred: {e}")

    # Calculate Loss
    loss_DCJ = (((float(results['DCJ_pred']) - float(input_params['DCJ_target'])) / float(input_params['DCJ_target']))
                ** 2)
    loss_TCJ = (((float(results['TCJ_pred']) - float(input_params['TCJ_target'])) / float(input_params['TCJ_target']))
                ** 2)
    loss_TvN = (((float(results['TvN_pred']) - float(input_params['TvN_target'])) / float(input_params['TvN_target']))
                ** 2)
    loss_Lind = (((float(results['li_pred']) - float(input_params['Lind_target'])) / float(input_params['Lind_target']))
                ** 2)
    loss_cell = (((float(mean_cell_size) - float(input_params['cell_size_target'])) / float(input_params['cell_size_target']))
                ** 2)
    loss_combined = (loss_DCJ + loss_TCJ + loss_TvN + loss_Lind + loss_cell)

    print(f"\n\033[97;41m Combined Loss: {loss_combined:.6f} \033[0m")

except Exception as e:
    print(f"\nAn error occurred: {e}")