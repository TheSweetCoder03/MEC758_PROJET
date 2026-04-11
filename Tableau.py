import numpy as np
from scipy.interpolate import RegularGridInterpolator

# Figure 1 : beta1 = 0
sc_axis_f1 = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
alpha2_f1 = np.array([40, 50, 60, 65, 70, 75, 80])

# Les données sont entrées comme : [Ligne = Alpha, Colonne = S/C]
data_f1 = np.array([
    [0.0475, 0.0325, 0.0250, 0.0210, 0.0190, 0.0180, 0.0175], # Angle 40
    [0.0500, 0.0350, 0.0280, 0.0240, 0.0220, 0.0215, 0.0230], # Angle 50
    [0.0513, 0.0380, 0.0320, 0.0290, 0.0275, 0.0290, 0.0340], # Angle 60
    [0.0538, 0.0410, 0.0370, 0.0350, 0.0355, 0.0385, 0.0440], # Angle 65
    [0.0538, 0.0465, 0.0440, 0.0445, 0.0470, 0.0520, 0.0585], # Angle 70
    [0.0563, 0.0520, 0.0510, 0.0530, 0.0570, 0.0635, 0.0715], # Angle 75
    [0.0638, 0.0600, 0.0590, 0.0620, 0.0675, 0.0760, 0.0845]  # Angle 80
])

# Figure 2 : beta1 = alpha2
sc_axis_f2 = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
alpha2_f2 = np.array([40, 50, 60, 65, 70]) 

data_f2 = np.array([
    [0.114, 0.082, 0.071, 0.065, 0.066, 0.071, 0.078], # Angle 40
    [0.118, 0.089, 0.078, 0.072, 0.076, 0.084, 0.093], # Angle 50
    [0.138, 0.111, 0.100, 0.100, 0.106, 0.117, 0.131], # Angle 60
    [0.150, 0.125, 0.118, 0.123, 0.136, 0.153, 0.174], # Angle 65
    [0.162, 0.138, 0.135, 0.147, 0.166, 0.187, 0.198], # Angle 70
])

# Figure 14 : Trailing Edge Energy Coefficient
to_axis = np.array([0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40])
type_f14 = np.array([0, 1]) 

data_f14 = np.array([
    [0.000, 0.005, 0.017, 0.031, 0.049, 0.071, 0.092, 0.116, 0.140], # Nozzle
    [0.000, 0.002, 0.008, 0.017, 0.029, 0.041, 0.053, 0.065, 0.076]  # Impulse
])

# Moteurs d'interpolation
interp_f1 = RegularGridInterpolator((sc_axis_f1, alpha2_f1), data_f1.T, bounds_error=False, fill_value=None)
interp_f2 = RegularGridInterpolator((sc_axis_f2, alpha2_f2), data_f2.T, bounds_error=False, fill_value=None)

interp_f14 = RegularGridInterpolator((type_f14, to_axis), data_f14, method='linear', bounds_error=False, fill_value=None)

def extraire_donnee_graphique(x, y, figure=1):
    """
    x : s/c ou t/o
    y : alpha2 (pour fig 1 et 2)
    """
    if figure == 1:
        return interp_f1([x, y])[0]
    
    elif figure == 2:
        return interp_f2([x, y])[0]
    
    elif figure == 140: # Figure 14 - Axial Entry Nozzle
        return interp_f14([0, x])[0]
    
    elif figure == 141: # Figure 14 - Impulse Blading
        return interp_f14([1, x])[0]