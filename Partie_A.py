import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

cte_amb= {
    "pa": 87560, # Pression atmosphérique à 1219m altitude en Pa
    "ta": 37.8, # Température à 1219m altitude en degres celcius
    "vinf": 0, # Vitesse d'entrée en m/s
    "ve": 0, # Vitesse de sortie en m/s
    "altitude": 1219, # Altitude en mètres
}

cte_comp= {
    "yc": 1.4, 
    "cpc": 1005, # Capacité calorifique en J/kg.K
    "rpbp": 4.25, # Rapport de pression compresseur basse pression
    "ncbp": 0.86, # Rendement du compresseur basse pression
    "rphp": 2.65, # Rapport de pression compresseur haute pression
    "nchp": 0.84, # Rendement du compresseur haute pression
    "mpt": 5.443, # Débit d'air massique total en kg/s
    "pap": 0.10, # Pourcentage d'air utilisé pour le refroidissement des turbines
    "Ur": 300, # Vitesse de roation à la racine en m/s
    "rr": 0.5, # Rapport de rayon des ailettes du compresseur"
    "haller": 0.72, # Nombre de Haller
    "rpmlpc": 32500, # Vitesse de rotation du compresseur basse pression en rpm (valeur posée)
    "rpmhpc": 37500, # Vitesse de rotation du compresseur haute pression en rpm (valeur posée)
    "alpha1": 0 # Angle d'entrée de l'air dans le compresseur basse pression en degrés (valeur posée)
}

cte_cc= {
    "ycc": 1.32, 
    "cpcc": 1172, # Capacité calorifique en J/kg.K
    "f": 0.02, # Rapport de carburant à l'air
    "qr": 38.984e6, # Pouvoir calorifique du carburant en J/kg
    "ncc": 0.99, # Rendement de la chambre de combustion
    "deltapcc": 0.02, # Perte de pression dans la chambre de combustion
}

cte_turb= {
    "yt": 1.31, 
    "cpt": 1214, # Capacité calorifique en J/kg.K
    "ntbp": 0.90, # Rendement de la turbine basse pression
    "nthp": 0.88, # Rendement de la turbine haute pression
    "ntp": 0.93, # Rendement de la turbine de puissance
    "deltapit": 0.02, # Perte de pression dans les canaux inter-turbines
    "deltapet": 0.02, # Perte de pression dans le canal d'échappement
}

station = {}

def Station_1 ():
    # Station 1 : Entrée de l'air dans le compresseur basse pression
    p1 = cte_amb["pa"]
    t1 = cte_amb["ta"] + 273.15 # Convertir en Kelvin
    va = cte_amb["vinf"]
    cp = cte_comp["cpc"]
    y = cte_comp["yc"]
    s1 = 0 # Entropie de référence
    v1 = va

    to1 = t1 + (v1**2) / (2 * cp)
    po1 = p1* (to1/t1)**(y/(y-1))

    # Enregistrer les résultats de la station 1
    station[1] = {"Po": po1, "To": to1, "s": s1}

    return po1, to1, s1

def Station_2 (po1, to1, s1):
    # Station 2 : Entrée de l'air dans le compresseur haute pression
    rpbp = cte_comp["rpbp"]
    ncbp = cte_comp["ncbp"]
    y = cte_comp["yc"]
    mp = cte_comp["mpt"]
    cp = cte_comp["cpc"]
    haller = cte_comp["haller"]
    Ur = cte_comp["Ur"]
    rr = cte_comp["rr"] # Ratio des rayons des ailettes du compresseur
    rpmlpc = cte_comp["rpmlpc"] * (2 * np.pi / 60) #rpm du compresseur basse pression converti en rad/s
    alpha1 = cte_comp["alpha1"] # Angle d'entrée de l'air dans le compresseur basse pression en degrés
    altitude = cte_amb["altitude"]
    r = cp * ((y-1)/y)
    rho = po1 / (r * to1)

    po2 = po1 * rpbp
    to2s = to1 * (rpbp)**((y-1)/y)
    to2 = to1 + (to2s - to1) / ncbp

    wlpc = mp * cp * (to2 - to1)

    s2 = s1 + cp * np.log(to2/to1) - r * np.log(po2/po1)

    #Calcul du nombre d'étages
    Rr = Ur / rpmlpc
    Rt = Rr / rr
    rm = (Rt + Rr) / 2
    A = np.pi * (Rt**2 - Rr**2)
    v1 = mp / (rho * A)
    va = v1 * np.cos(alpha1)
    Um = rm * rpmlpc
    Vru2 = fsolve(lambda Vru2s : 0.72 - (np.sqrt(va**2 + Vru2s**2)/(np.sqrt(va**2 + (Um - Vru2s)**2))), x0 = [Um/2])[0]
    Vru1 = Um - Vru2
    deltaT0max = (Um * (Vru1 - Vru2)) / cp
    n_etages = int((to2 - to1) / deltaT0max) + 1

    station[2] = {"Po": po2, "To": to2, "w": wlpc, "s": s2, "Nb_etage": n_etages}

    return po2, to2, wlpc, s2

def Station_3 (po2, to2, s2):
    # Station 3 : Entrée de l'air dans la chambre à combustion
    rphp = cte_comp["rphp"]
    nchp = cte_comp["nchp"]
    y = cte_comp["yc"]
    mp = cte_comp["mpt"]
    cp = cte_comp["cpc"]
    rpmhpc = cte_comp["rpmhpc"] * (2 * np.pi / 60) #Conversion du rpm du HPC en rad/s
    Ur = cte_comp["Ur"]
    rr = cte_comp["rr"] # Ratio des rayons des ailettes du compresseur
    alpha1 = cte_comp["alpha1"]
    haller = cte_comp["haller"]

    po3 = po2 * rphp
    to3s = to2 * (rphp)**((y-1)/y)
    to3 = to2 + (to3s - to2) / nchp

    whpc = mp * cp * (to3 - to2)

    r = cp * ((y-1)/y)

    s3 = s2 + cp * np.log(to3/to2) - r * np.log(po3/po2)

    # Calcul du nombre d'étages

    alpha1hpc = alpha1 # On suppose que l'angle d'entrée dans le HPC est le même que dans le LPC
    Ma1 = 0.5 # Hypothèse que le nombre de Mach dans le compresseur haute pression est de 0.5
    to1hpc = to2 # On ramène la température à l'entrée du HPC, donc station 2, à la station 1 du HPC
    t1hpc = to1hpc / (1 + (Ma1**2 * y * r) / (2 * cp)) # Température statique à l'entrée du HPC
    v1hpc = Ma1 * np.sqrt(y * r * t1hpc)
    va1hpc = v1hpc * np.cos(alpha1hpc)

    Rr = Ur / rpmhpc
    Rt = Rr / rr
    rm = (Rt + Rr) / 2
    Um = rm * rpmhpc
    Vru2 = fsolve(lambda Vru2s : 0.72 - (np.sqrt(va1hpc**2 + Vru2s**2)/(np.sqrt(va1hpc**2 + (Um - Vru2s)**2))), x0 = [Um/2])[0]
    Vru1 = Um - Vru2
    deltaT0max = (Um * (Vru1 - Vru2)) / cp
    n_etages = int((to3 - to1hpc) / deltaT0max) + 1

    station[3] = {"Po": po3, "To": to3, "w": whpc, "s": s3, "Nb_etage": n_etages}

    return po3, to3, whpc, s3

def Station_4 (po3, to3, s3):
    # Station 4 : Entrée de l'air dans la turbine haute pression
    f = cte_cc["f"]
    qr = cte_cc["qr"]
    ncc = cte_cc["ncc"]
    deltapcc = cte_cc["deltapcc"]
    cp4 = cte_cc["cpcc"]
    cp3 = cte_comp["cpc"]
    y = cte_cc["ycc"]
    mp = cte_comp["mpt"]
    pap = cte_comp["pap"]
    cpt = cte_turb["cpt"]
    cpc = cte_comp["cpc"]

    to4 = (f*qr*ncc+cp3*to3)/(cp4*(1+f))
    po4 = po3 * (1 - deltapcc)

    r = cp4 * ((y-1)/y)

    s4 = s3 + cp4 * np.log(to4/to3) - r * np.log(po4/po3)

    #Calcul de la température To4 une fois que l'aire de refroidissement est mixée
    mpt = mp * (1 - pap) * (1 + f) + mp * pap
    to4u = ((mp * (1 - pap) * (1 + f)) * to4 * cpt + mp * pap * to3 * cpc)/ (mpt * cpt)

    station[4] = {"Po": po4, "To": to4u, "s": s4}

    return po4, to4u, s4, mpt

def Station_5 (po4, to4u, whpc, s4, mpt):
    # Station 5 : Entrée de l'air dans la turbine basse pression
    nthp = cte_turb["nthp"]
    y = cte_turb["yt"]
    cpt = cte_turb["cpt"]
    
    to5 = to4u-(whpc/(mpt*cpt))
    to5s = to4u - (to4u-to5)/nthp
    po5 = po4 * (to5s/to4u)**(y/(y-1))

    whpt = mpt * cpt * (to4u - to5)

    r = cpt * ((y-1)/y)

    s5 = s4 + cpt * np.log(to5/to4u) - r * np.log(po5/po4)

    station[5] = {"Po": po5, "To": to5, "w": whpt, "s": s5, "mpt": mpt}

    return po5, to5, whpt, s5, mpt

def Station_6 (po5, to5, wlpc, s5, mpt):
    # Station 6 : Entrée de l'air dans la power turbine
    ntbp = cte_turb["ntbp"]
    y = cte_turb["yt"]
    cpt = cte_turb["cpt"]
    perte_it = cte_turb["deltapit"]

    po5p = po5 * (1 - perte_it)

    to6 = to5-(wlpc/(mpt*cpt))
    to6s = to5 - (to5-to6)/ntbp
    po6 = po5p * (to6s/to5)**(y/(y-1))

    wlpt = mpt * cpt * (to5 - to6)

    r = cpt * ((y-1)/y)

    s6 = s5 + cpt * np.log(to6/to5) - r * np.log(po6/po5p)

    station[6] = {"Po": po6, "To": to6, "w": wlpt, "s": s6}

    return po6, to6, wlpt, s6

def Station_7 (po6, to6, s6, mpt):
    # Station 7 : Sortie de la turbine de puissance
    ntp = cte_turb["ntp"]
    y = cte_turb["yt"]
    cpt = cte_turb["cpt"]
    perte_et = cte_turb["deltapet"]
    pa = cte_amb["pa"]
    f = cte_cc["f"]
    mp = cte_comp["mpt"]
    pap = cte_comp["pap"]
    perte_it = cte_turb["deltapit"]

    po6p = po6 * (1 - perte_it)

    p7 = pa * (1 + perte_et)

    t7s = to6 * (p7/po6p)**((y-1)/y)
    t7 = to6 - (to6-t7s)*ntp

    wpt = mpt * cpt * (to6 - t7)
    hp = wpt / 745.7

    r = cpt * ((y-1)/y)

    s7 = s6 + cpt * np.log(t7/to6) - r * np.log(p7/po6p)

    mpcc = mp * (1 - pap)

    sfc = (mpcc * f * 3600) / (wpt / 1000)

    station[7] = {"Po": p7, "To": t7, "w": wpt, "hp": hp, "s": s7, "sfc": sfc}

    return p7, t7, wpt, hp, s7, sfc

def print_station():
    print(f"{'Station':<10} | {'Po (Pa)':<10} | {'To (K)':<10} | {'W (Watts)':<10}")
    print("-" * 34)
    for i in range(1, 8):
        if i in station:
            if 'w' in station[i]:
                print(f"{i:<10} | {station[i]['Po']:<10.2f} | {station[i]['To']:<10.2f} | {station[i]['w']:<10.2f}")
            else:
                print(f"{i:<10} | {station[i]['Po']:<10.2f} | {station[i]['To']:<10.2f} | {'N/A':<10}")
    print(f"\nSFC: {station[7]['sfc']:.6f} kg/kW.h")
    print(f"Puissance de la turbine de puissance: {station[7]['hp']:.2f} HP")
    print(f"Le nombre d'étages LPC : {station[2]['Nb_etage']:.2f} étages")
    print(f"Le nombre d'étages HPC : {station[3]['Nb_etage']:.2f} étages")

def plot_cycle():

    s_vals = []
    t_vals = []

    for i in range(1,8):
        if i in station:
            
            s_actuel = station[i]['s']
            t_actuel = station[i]['To']

            s_vals.append(s_actuel)
            t_vals.append(t_actuel)
    
    plt.figure(1, figsize=(10, 6))
    plt.plot(s_vals[0:3], t_vals[0:3], marker='o', linestyle='-', color='red', linewidth=2, label='Avant la C.C.')
    plt.plot(s_vals[3:], t_vals[3:], marker='o', linestyle='-', color='blue', linewidth=2, label='Après la C.C.')

    xmin, xmax = plt.xlim()

    for i, (s, t) in enumerate(zip(s_vals, t_vals), start=1):
        plt.annotate(f'Station {i}', (s, t), textcoords="offset points", xytext=(8,-5), ha='left')
        plt.plot([xmin-10,s], [t, t], linestyle='--', color='gray', linewidth=1, alpha=0.5)
        plt.text(xmin-10, t, f'To{i}', ha='right', va='center', fontsize=9, color='gray')
    
    plt.xlim(xmin, xmax)
    plt.title('Diagramme T-s')
    plt.xlabel('Entropie (s)')
    plt.ylabel('Température de stagnation (To) [K]')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.show()


def calcul():
    po1, to1, s1 = Station_1()
    po2, to2, wlpc, s2 = Station_2(po1, to1, s1)
    po3, to3, whpc, s3 = Station_3(po2, to2, s2)
    po4, to4u, s4, mpt = Station_4(po3, to3, s3)
    po5, to5, whpt, s5, mpt = Station_5(po4, to4u, whpc, s4, mpt)
    po6, to6, wlpt, s6 = Station_6(po5, to5, wlpc, s5, mpt)
    p7, t7, wpt, hp, s7, sfc = Station_7(po6, to6, s6, mpt)

def export_donnees_hpt():
    mpt_calc = station[5]['mpt']
    donnees_hpt = {
        'T04': station[4]['To'],                  
        'P04': station[4]['Po'],                  
        'T05': station[5]['To'],                  
        'P05': station[5]['Po'],                  
        'm_dot': mpt_calc,                
        'W_hpt': station[5]['w'],               
        'cp': cte_turb["cpt"],       
        'gamma': cte_turb["yt"],     
        'eta_iso': cte_turb["nthp"],  
        'dh0_hpt': station[5]['w'] / mpt_calc, # Travail spécifique (J/kg)
        'rpmhpc': cte_comp["rpmhpc"]
    }

    return donnees_hpt