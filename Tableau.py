import numpy as np
from scipy.interpolate import RegularGridInterpolator

# Figure 1 : beta1 = 0 
sc_axis_f1 = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
alpha2_f1 = np.array([40, 50, 60, 65, 70, 75, 80])
data_f1 = np.array([
    [0.044, 0.046, 0.048, 0.051, 0.053, 0.057, 0.063],
    [0.032, 0.034, 0.038, 0.043, 0.047, 0.054, 0.060],
    [0.024, 0.026, 0.032, 0.038, 0.044, 0.053, 0.057],
    [0.020, 0.022, 0.029, 0.035, 0.043, 0.056, 0.062],
    [0.018, 0.021, 0.028, 0.034, 0.046, 0.062, 0.070],
    [0.017, 0.021, 0.030, 0.040, 0.053, 0.071, 0.081],
    [0.017, 0.022, 0.035, 0.048, 0.063, 0.082, 0.095]
])

# Figure 2 : beta1 = alpha2
sc_axis_f2 = np.array([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
alpha2_f2 = np.array([40, 50, 55, 60, 65, 70])
data_f2 = np.array([
    [0.125, 0.130, 0.135, 0.140, 0.155, 0.160],
    [0.092, 0.098, 0.103, 0.110, 0.130, 0.140],
    [0.075, 0.082, 0.088, 0.098, 0.118, 0.138],
    [0.068, 0.075, 0.083, 0.095, 0.120, 0.145],
    [0.068, 0.078, 0.088, 0.102, 0.130, 0.160],
    [0.072, 0.085, 0.098, 0.115, 0.145, 0.180],
    [0.080, 0.095, 0.110, 0.125, 0.160, 0.200]
])

# Figure 14 (Pertes bord de fuite)
# Axe X : Ratio t/o (épaisseur bord fuite / ouverture col) 
# Axe Y : Type (0 = Axial Entry Nozzle (beta1=0), 1 = Impulse Blading (beta1=alpha2))
to_axis = np.array([0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40])
type_f14 = np.array([0, 1]) 
data_f14 = np.array([
    [0.000, 0.003, 0.012, 0.028, 0.050, 0.076, 0.106, 0.125, 0.142], # Nozzle
    [0.000, 0.002, 0.007, 0.016, 0.028, 0.043, 0.060, 0.072, 0.082]  # Impulse
])

# Moteurs d'interpolation
interp_f1 = RegularGridInterpolator((sc_axis_f1, alpha2_f1), data_f1, bounds_error=False, fill_value=None)
interp_f2 = RegularGridInterpolator((sc_axis_f2, alpha2_f2), data_f2, bounds_error=False, fill_value=None)
interp_f14 = RegularGridInterpolator((type_f14, to_axis), data_f14, method='linear', bounds_error=False, fill_value=None)

def extraire_donnee_graphique(x, y, figure=1):
    point = [[x, y]]
    
    if figure == 1:
        return interp_f1(point).item()
    
    elif figure == 2:
        return interp_f2(point).item()
    
    elif figure == 140: # Figure 14 - Axial Entry Nozzle
        return interp_f14([[0, x]]).item()
    
    elif figure == 141: # Figure 14 - Impulse Blading
        return interp_f14([[1, x]]).item()

