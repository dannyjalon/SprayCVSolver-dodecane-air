# Spray Detonation Modeling and Analysis Suite

1. Overview

This project provides a comprehensive suite for the one-dimensional numerical simulation of two-phase (liquid spray/gas) detonations. It implements a ZND model that accounts for droplet evaporation and simplified one-step chemical kinetics.

The suite recieves the following target value parameters, either from experiments or numerical simulations, and runs a calibration procedure to best match these values.   
* Detonation velocity, Dcj, in m/s.
* von Neumann temperature, TvN, in K.
* Chapman-Jouguet temperature, Tcj, in K.
* Induction length, Lind, in mm.
* Activation Energy, Ea, in J/mol.

The suite also receives an image of the soot foil (experiments or numerical) and measures the mean cell height or width.  

2. Features
* 1D Two-Phase ZND Model: Simulates the detonation structure, including gas and droplet phases.
* Soot Foil Image Analysis: Integrates a powerful computer vision tool (soot_foil_image_tool.py) to measure detonation cell sizes from experimental or numerical images.
* Loss Function Calculation: The main driver script (run.py) calculates a final loss value, quantifying the model's accuracy against target values (experimental or numerical) for detonation velocity, von Neumann temperature, CJ temperature and cell size (length or width).
* Modular Structure: The code is separated into a core simulation library (spray.py) and a user-facing driver script (run.py).

3. Prerequisites
This project is written in Python 3. The following libraries are required:

NumPy

SciPy

OpenCV for Python (opencv-python)

Matplotlib

tqdm

You can install all dependencies at once by creating a requirements.txt file with the contents above and running:

pip install -r requirements.txt

4. How to Use
The entire analysis is orchestrated through the main driver script, run.py.

Step-by-step instructions:

* Navigate to the Directory: Open a terminal and navigate to the root directory of the project.

* Execute the Script: Run the following command:

python run.py

* Provide target parameters into "input_params" inside 'run.py', for example:

  input_params = {
    'DCJ_target': 1805.0,   # Target CJ velocity [m/s]
  
    'TvN_target': 1689.0,   # Target von Neumann temperature [K]
  
    'TCJ_target': 2875.0,   # Target CJ temperature [K]
  
    'Ea': 5546,             # Activation energy [J/mol]
  
    'Lind_target': 20.096,   # Target induction length [mm]
  
    'cell_size_target': 17  # Target Cell Size [cm]
  
}

* Provide Image Analysis Input: 

Image Path: Provide the full, absolute path to your .png or .jpg soot foil image.

Dimension: Enter h or w depending on whether the known physical size corresponds to the image's height or width.

Physical Size: Enter the known physical size in centimeters (e.g., 50.0).

Step: Deafult is 10. Used for the artifact filtration function (parameter optimization). Default value is 10. Decreasing the step size may provide better results but will increase the running time, and vice versa.

* Monitor the Solver: The script will try to run the simulation with the provided target values. If it fails, an error will be raised and the loss function will manually set to an extremely high value. After the simulation is complete, the script will analyze the image to determine the cell size.

View Results: the solver prints two key outputs:

A dictionary of the final predicted model values (DCJ_pred, TvN_pred, etc.).

A final combined loss value, displayed on a red background, which quantifies the model's error against the experimental\numerical data.

5. File Structure
run.py: Main execution script. Handles user input, calls the image analysis and simulation modules, and calculates the final loss.

spray.py: The core simulation library. Contains the SprayModel function which orchestrates the B-factor solver.

soot_foil_image_tool.py: The computer vision module for analyzing soot foils.

state_vector_derivatives.py: Defines the system of ordinary differential equations (ODEs) for the ZND model.

cal1.py / cal2.py: Helper scripts to determine thermodynamic properties and activation energy.

post_shock_conditions.py, set_initial_conditions.py, etc.: Other helper modules for the simulation setup.

