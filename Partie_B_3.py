import numpy as np
def deg2rad(angle):
    return angle * np.pi / 180.0

def etape_3(donnees_hpt):
    print(f"\nÉTAPES 3: Paramètres des aubes")
    fs = 0.7
    fr = 1.3
    zweifels = 0.75
    zweifelr = 0.9
    h1 = rt1 - rr1
    h2 = rt2 - rr2
    h3 = rt3 - rr3
    rm2 = (rt2 + rr2) / 2
    rm3 = (rt3 + rr3) / 2
    a1 = 1
    a2 = 1
    ar1 = 1
    ar2 = 1
    b1 = 1
    b2 = 1
    b3 = 1

    # Calcul hauteur moyenne ailette
    hs = (h1 + h2) / 2
    hr = (h2 + h3) / 2

    # Calcul de la corde axiale
    cs = hs / fs
    cr = hr / fr
    ys = (a1 + a2)  / 2
    yr = (b1 + b2) / 2
    cas = cs * np.cos(ys)
    car = cr * np.cos(yr)

    # Calcul du pas
    pas_s = (zweifels * cas) / (2 * (np.tan(ar1)+np.tan(ar2)) * (np.cos(ar2))**2)
    pas_r = (zweifelr * car) / (2 * abs(np.tan(b2)+np.tan(b3)) * (np.cos(b3))**2)

    # Calcul du nombre d'ailette 
    ns = 2 * np.pi * rm2 / pas_s
    nr = 2 * np.pi * rm3 / pas_r






