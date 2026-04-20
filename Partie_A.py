import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

cte_amb = {
    "pa": 87560, # Pression atmosphérique à 1219m altitude en Pa
    "ta": 37.8, # Température à 1219m altitude en degres celcius
    "vinf": 0, # Vitesse d'entrée en m/s
    "ve": 0, # Vitesse de sortie en m/s
    "altitude": 1219, # Altitude en mètres
}

cte_comp = {
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
    "rpmlpc": 32500, # Vitesse de rotation en rpm (LPC)
    "rpmhpc": 37500, # Vitesse de rotation en rpm (HPC)
    "alpha1": 0 # Angle d'entrée en degrés
}

cte_cc = {
    "ycc": 1.32, 
    "cpcc": 1172, # Capacité calorifique en J/kg.K
    "f": 0.02, # Rapport de carburant à l'air
    "qr": 38.984e6, # Pouvoir calorifique du carburant en J/kg
    "ncc": 0.99, # Rendement de la chambre de combustion
    "deltapcc": 0.02, # Perte de pression dans la chambre de combustion
}

cte_turb = {
    "yt": 1.31, 
    "cpt": 1214, # Capacité calorifique en J/kg.K
    "ntbp": 0.90, # Rendement turbine basse pression
    "nthp": 0.88, # Rendement turbine haute pression
    "ntp": 0.93, # Rendement turbine de puissance
    "deltapit": 0.02, # Perte de pression dans les canaux inter-turbines
    "deltapet": 0.02, # Perte de pression dans le canal d'échappement
}

station = {}

def Station_1():
    p1 = cte_amb["pa"]
    t1 = cte_amb["ta"] + 273.15 
    va = cte_amb["vinf"]
    cp = cte_comp["cpc"]
    y = cte_comp["yc"]
    s1 = 0 
    v1 = va

    to1 = t1 + (v1**2) / (2 * cp)
    po1 = p1* (to1/t1)**(y/(y-1))

    station[1] = {"Po": po1, "To": to1, "s": s1, "cp": cp}
    return po1, to1, s1

def Station_2(po1, to1, s1):
    rpbp = cte_comp["rpbp"]
    ncbp = cte_comp["ncbp"]
    y = cte_comp["yc"]
    mp = cte_comp["mpt"]
    cp = cte_comp["cpc"]
    Ur = cte_comp["Ur"]
    rr = cte_comp["rr"] 
    rpmlpc = cte_comp["rpmlpc"] * (2 * np.pi / 60)
    alpha1 = cte_comp["alpha1"]
    r = cp * ((y-1)/y)
    rho = po1 / (r * to1)

    po2 = po1 * rpbp
    to2s = to1 * (rpbp)**((y-1)/y)
    to2 = to1 + (to2s - to1) / ncbp

    wlpc = mp * cp * (to2 - to1)
    s2 = s1 + cp * np.log(to2/to1) - r * np.log(po2/po1)

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

    station[2] = {"Po": po2, "To": to2, "w": wlpc, "s": s2, "Nb_etage": n_etages, "cp": cp}
    return po2, to2, wlpc, s2

def Station_3(po2, to2, s2):
    rphp = cte_comp["rphp"]
    nchp = cte_comp["nchp"]
    y = cte_comp["yc"]
    mp = cte_comp["mpt"]
    cp = cte_comp["cpc"]
    rpmhpc = cte_comp["rpmhpc"] * (2 * np.pi / 60)
    Ur = cte_comp["Ur"]
    rr = cte_comp["rr"]
    alpha1 = cte_comp["alpha1"]

    po3 = po2 * rphp
    to3s = to2 * (rphp)**((y-1)/y)
    to3 = to2 + (to3s - to2) / nchp

    whpc = mp * cp * (to3 - to2)
    r = cp * ((y-1)/y)
    s3 = s2 + cp * np.log(to3/to2) - r * np.log(po3/po2)

    alpha1hpc = alpha1
    Ma1 = 0.5 
    to1hpc = to2
    t1hpc = to1hpc / (1 + (Ma1**2 * y * r) / (2 * cp)) 
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

    station[3] = {"Po": po3, "To": to3, "w": whpc, "s": s3, "Nb_etage": n_etages, "cp": cp}
    return po3, to3, whpc, s3

def Station_4(po3, to3, s3):
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

    mpt = mp * (1 - pap) * (1 + f) + mp * pap
    to4u = ((mp * (1 - pap) * (1 + f)) * to4 * cp4 + mp * pap * to3 * cpc)/ (mpt * cpt)

    station[4] = {"Po": po4, "To": to4u, "s": s4, "cp": cpt}
    return po4, to4u, s4, mpt

def Station_5(po4, to4u, whpc, s4, mpt):
    nthp = cte_turb["nthp"]
    y = cte_turb["yt"]
    cpt = cte_turb["cpt"]
    
    to5 = to4u-(whpc/(mpt*cpt))
    to5s = to4u - (to4u-to5)/nthp
    po5 = po4 * (to5s/to4u)**(y/(y-1))

    whpt = whpc
    r = cpt * ((y-1)/y)
    s5 = s4 + cpt * np.log(to5/to4u) - r * np.log(po5/po4)

    station[5] = {"Po": po5, "To": to5, "w": whpt, "s": s5, "mpt": mpt, "cp": cpt}
    
    # Perte inter-turbine 1
    perte_it = cte_turb["deltapit"]
    po5p = po5 * (1 - perte_it)
    s5p = s5 - r * np.log(po5p/po5) # Isenthalpique (T0 = cte, S augmente)
    station["5p"] = {"Po": po5p, "To": to5, "s": s5p, "cp": cpt}

    return po5, to5, whpt, s5, mpt

def Station_6(po5, to5, wlpc, s5, mpt):
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

    station[6] = {"Po": po6, "To": to6, "w": wlpt, "s": s6, "cp": cpt}

    # Perte inter-turbine 2 (selon votre algorithme station 7)
    po6p = po6 * (1 - perte_it)
    s6p = s6 - r * np.log(po6p/po6)
    station["6p"] = {"Po": po6p, "To": to6, "s": s6p, "cp": cpt}

    return po6, to6, wlpt, s6

def Station_7(po6, to6, s6, mpt):
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

    station[7] = {"Po": p7, "To": t7, "w": wpt, "hp": hp, "s": s7, "sfc": sfc, "cp": cpt}
    
    # Chute d'échappement à Patm
    p_atm = pa
    s_atm = s7 - r * np.log(p_atm/p7)
    station["atm"] = {"Po": p_atm, "To": t7, "s": s_atm, "cp": cpt}

    return p7, t7, wpt, hp, s7, sfc

def plot_cycle_enhanced():
    s_vals = []
    t_vals = []
    for i in range(1,8):
        if i in station:
            s_vals.append(station[i]['s'])
            t_vals.append(station[i]['To'])
            
    plt.figure(figsize=(14, 9))

    # Lignes (Compression/Combustion)
    plt.plot(s_vals[0:3], t_vals[0:3], marker='o', linestyle='-', color='blue', linewidth=2.5, label='Compression')
    plt.plot(s_vals[2:4], t_vals[2:4], marker='o', linestyle='-', color='red', linewidth=2.5, label='Chambre de Combustion')
    
    # Lignes (Turbines + Pertes intermédiaires)
    s_turb = [station[4]['s'], station[5]['s'], station['5p']['s'], station[6]['s'], station['6p']['s'], station[7]['s'], station['atm']['s']]
    t_turb = [station[4]['To'], station[5]['To'], station['5p']['To'], station[6]['To'], station['6p']['To'], station[7]['To'], station['atm']['To']]
    plt.plot(s_turb, t_turb, marker='o', linestyle='-', color='darkgreen', linewidth=2.5, label='Détente, Pertes IT & Échappement')
    
    # Ajustement des limites
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    xmax_plot = xmax + 200
    ymax_plot = ymax + 100
    plt.xlim(xmin - 50, xmax_plot)
    plt.ylim(ymin - 50, ymax_plot)
    
    s_range = np.linspace(min(s_vals)-100, max(s_vals)+200, 100)
    
    top_offset = 15 # Compteur pour décaler les textes au plafond
    
    def plot_isobar(s_ref, t_ref, p_val, cp, color='gray', linewidth=1.0, alpha=0.5, is_red=False):
        nonlocal top_offset
        t_iso = t_ref * np.exp((s_range - s_ref) / cp)
        plt.plot(s_range, t_iso, linestyle=':', color=color, linewidth=linewidth, alpha=alpha)
        
        s_end = s_range[-1]
        t_end = t_iso[-1]
        
        # Stagger les textes qui sortent par le haut
        if t_end > ymax_plot:
            t_end = ymax_plot - top_offset
            s_end = s_ref + cp * np.log(t_end / t_ref)
            top_offset += 25 
            
        # MODIFICATION : Décaler le texte rouge vers le BAS (axe Y)
        y_shift = -20 if is_red else 0 # Descend de 20 unités sur l'axe de la température
        
        # Affichage du texte
        plt.text(s_end, t_end + y_shift, f"{p_val/1000:.1f} kPa", 
                 color=color, fontsize=9, va='center', ha='right',
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))

    # Tracé des isobares et marqueurs
    stations_to_plot = [1, 2, 3, 4, 5, '5p', 6, '6p', 7, 'atm']
    for k in stations_to_plot:
        st = station[k]
        
        if 'p' in str(k) or str(k) == 'atm':
            plt.plot(st['s'], st['To'], 'rx', markersize=8) 
            plot_isobar(st['s'], st['To'], st['Po'], st['cp'], color='red', linewidth=1.5, alpha=0.7, is_red=True)
        else:
            plt.plot(st['s'], st['To'], 'ko', markersize=5)
            plot_isobar(st['s'], st['To'], st['Po'], st['cp'], color='gray', is_red=False)

        # Placement intelligent des étiquettes de stations
        ha_val = 'left'
        va_val = 'bottom'
        label = str(k)
        
        if label == '1':
            xytext = (10, -20)
            va_val = 'top'
            label_text = f"St. {label}\nTo={st['To']:.1f} K\nPo={st['Po']/1000:.1f} kPa"
        elif label in ['2', '3']:
            xytext = (-10, 10)
            ha_val = 'right'
            label_text = f"St. {label}\nTo={st['To']:.1f} K\nPo={st['Po']/1000:.1f} kPa"
        elif label == '4':
            xytext = (10, 10)
            label_text = f"St. {label}\nTo={st['To']:.1f} K\nPo={st['Po']/1000:.1f} kPa"
        elif label in ['5', '6', '7']:
            xytext = (-15, 5)
            ha_val = 'right'
            label_text = f"St. {label}\nTo={st['To']:.1f} K\nPo={st['Po']/1000:.1f} kPa"
        elif label in ['5p', '6p']:
            xytext = (15, -15)
            ha_val = 'left'
            va_val = 'top'
            label_text = f"Pertes IT\nPo={st['Po']/1000:.1f} kPa"
        elif label == 'atm':
            xytext = (15, -15)
            ha_val = 'left'
            va_val = 'top'
            label_text = f"Échappement\nPo={st['Po']/1000:.1f} kPa"

        plt.annotate(label_text, (st['s'], st['To']), textcoords="offset points", xytext=xytext, 
                     ha=ha_val, va=va_val, fontsize=8, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
        
        plt.axhline(y=st['To'], color='gray', linestyle='--', linewidth=0.5, alpha=0.4)

    plt.title('Diagramme T-s avec Pressions de Stagnation et Pertes', fontsize=16)
    plt.xlabel('Entropie spécifique (s) [J/kg.K]', fontsize=14)
    plt.ylabel('Température de stagnation (To) [K]', fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend(loc='upper left', fontsize=12)
    plt.tight_layout()
    plt.show()

def verifier_etranglement(nom_turbine, po_in, to_in, po_out, m_dot, cp, y):
    """
    Vérifie si la turbine est étranglée (choked) et calcule la section au col (A*).
    """
    R = cp * (y - 1) / y
    
    pr_crit = ((y + 1) / 2) ** (y / (y - 1))
    
    pr_reel = po_in / po_out
    
    # La turbine est étranglée si le rapport de pression réel dépasse le critique
    est_etrangle = pr_reel >= pr_crit
    
    # Calcul des propriétés au col (Mach = 1) pour déterminer la section critique A*
    T_star = to_in * (2 / (y + 1))
    P_star = po_in * (2 / (y + 1)) ** (y / (y - 1))
    
    rho_star = P_star / (R * T_star)
    V_star = np.sqrt(y * R * T_star) # Vitesse du son locale
    
    A_col = m_dot / (rho_star * V_star) # Section au col en m^2
    
    print(f"{nom_turbine:<30} | PR Réel: {pr_reel:<5.2f} | PR Crit: {pr_crit:<5.2f} | État: {'ÉTRANGLÉ' if est_etrangle else 'NON ÉTRANGLÉ':<12} | Section col A*: {A_col*10000:.2f} cm²")
    
    return est_etrangle, A_col

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

# Exécution du script mis à jour
def calcul():
    po1, to1, s1 = Station_1()
    po2, to2, wlpc, s2 = Station_2(po1, to1, s1)
    po3, to3, whpc, s3 = Station_3(po2, to2, s2)
    po4, to4u, s4, mpt = Station_4(po3, to3, s3)
    po5, to5, whpt, s5, mpt = Station_5(po4, to4u, whpc, s4, mpt)
    po6, to6, wlpt, s6 = Station_6(po5, to5, wlpc, s5, mpt)
    p7, t7, wpt, hp, s7, sfc = Station_7(po6, to6, s6, mpt)
    print("ANALYSE D'ÉTRANGLEMENT DES TURBINES (CHOKING)")
    # Turbine Haute Pression (HPT) : Station 4 vers Station 5
    verifier_etranglement("Turbine Haute Pression (HPT)", po4, to4u, po5, mpt, cte_turb["cpt"], cte_turb["yt"])
    # Turbine Basse Pression (LPT) : Station 5p vers Station 6
    verifier_etranglement("Turbine Basse Pression (LPT)", station["5p"]["Po"], station["5p"]["To"], po6, mpt, cte_turb["cpt"], cte_turb["yt"])
    # Turbine de Puissance (PT) : Station 6p vers Station 7
    verifier_etranglement("Turbine de Puissance (PT)", station["6p"]["Po"], station["6p"]["To"], p7, mpt, cte_turb["cpt"], cte_turb["yt"])
    print("="*80 + "\n")

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