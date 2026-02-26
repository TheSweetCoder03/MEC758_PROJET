import numpy as np
#import matplotlib.pyplot as plt

cte_amb= {
    "pa": 87260, # Pression atmosphérique à 1219m altitude en Pa
    "ta": 38.7, # Température à 1219m altitude en degres celcius
    "va": 0, # Vitesse d'entrée en m/s
    "ve": 0, # Vitesse de sortie en m/s
}

cte_comp= {
    "yc": 1.4, 
    "cpc": 1005, # Capacité calorifique en J/kg.K
    "rpbp": 4.25, # Rapport de pression compresseur basse pression
    "ncbp": 0.86, # Rendement du compresseur basse pression
    "rphp": 2.65, # Rapport de pression compresseur haute pression
    "nchp": 0.84, # Rendement du compresseur haute pression
    "mpt": 5.443, # Débit d'air massique en kg/s
    "pap": 0.10, # Pourcentage d'air utilisé pour le refroidissement des turbines
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
    v1 = cte_amb["va"]
    cp = cte_comp["cpc"]
    y = cte_comp["yc"]
    
    to1 = t1 + (v1**2) / (2 * cp)
    po1 = p1* (to1/t1)**(y/(y-1))

    station[1] = {"Po": po1, "To": to1}

    return po1, to1

def Station_2 (po1, to1):
    # Station 2 : Sortie du compresseur basse pression
    rpbp = cte_comp["rpbp"]
    ncbp = cte_comp["ncbp"]
    y = cte_comp["yc"]
    mp = cte_comp["mpt"]
    cp = cte_comp["cpc"]


    po2 = po1 * rpbp
    to2s = to1 * (rpbp)**((y-1)/y)
    to2 = to1 + (to2s - to1) / ncbp

    wlpc = mp * cp * (to2 - to1)

    station[2] = {"Po": po2, "To": to2, "wlpc": wlpc}

    return po2, to2, wlpc

def Station_3 (po2, to2):
    # Station 3 : Sortie du compresseur haute pression
    rphp = cte_comp["rphp"]
    nchp = cte_comp["nchp"]
    y = cte_comp["yc"]
    mp = cte_comp["mpt"]
    cp = cte_comp["cpc"]

    po3 = po2 * rphp
    to3s = to2 * (rphp)**((y-1)/y)
    to3 = to2 + (to3s - to2) / nchp

    whpc = mp * cp * (to3 - to2)

    station[3] = {"Po": po3, "To": to3, "whpc": whpc}

    return po3, to3, whpc

def Station_4 (po3, to3):
    # Station 4 : Sortie de la chambre de combustion
    f = cte_cc["f"]
    qr = cte_cc["qr"]
    ncc = cte_cc["ncc"]
    deltapcc = cte_cc["deltapcc"]
    cp4 = cte_cc["cpcc"]
    cp3 = cte_comp["cpc"]
    mpt = cte_comp["mpt"]
    pap = cte_comp["pap"]

    mpr = mpt * (1 - pap) # Masse d'air réel passant par la CC

    to4 = (f*qr*ncc+cp3*to3)//(cp4*(1+f))
    po4 = po3 * (1 - deltapcc)
    sfc = f / mpr

    station[4] = {"Po": po4, "To": to4, "sfc": sfc}

    return po4, to4, mpr, sfc

def Station_5 (po4, to4,whpc, mpr):
    # Station 5 : Sortie de la turbine haute pression
    nthp = cte_turb["nthp"]
    y = cte_turb["yt"]
    cpt = cte_turb["cpt"]
    cpc = cte_comp["cpc"]

    to5 = to4-(whpc/(mpr*cpt))
    to5s = to4 - (to4-to5)/nthp
    po5 = po4 * (to5s/to4)**(y/(y-1))

    whpt = mpr * cpt * (to5 - to4)

    station[5] = {"Po": po5, "To": to5, "whpt": whpt}

    return po5, to5, whpt, mpr

def Station_6 (po5, to5, mpr, wlpc):
    # Station 6 : Sortie de la turbine basse pression
    ntbp = cte_turb["ntbp"]
    y = cte_turb["yt"]
    cpt = cte_turb["cpt"]
    perte = cte_turb["deltapit"]

    po5p = po5 * (1 - perte)

    to6 = to5-(wlpc/(mpr*cpt))
    to6s = to5 - (to5-to6)/ntbp
    po6 = po5p * (to6s/to5)**(y/(y-1))

    wlpt = mpr * cpt * (to6 - to5)

    station[6] = {"Po": po6, "To": to6, "wlpt": wlpt}

    return po6, to6, wlpt

def Station_7 (po6, to6):
    # Station 7 : Sortie de la turbine de puissance
    ntp = cte_turb["ntp"]
    y = cte_turb["yt"]
    cpt = cte_turb["cpt"]
    perte_it = cte_turb["deltapet"]
    perte_et = cte_turb["deltapet"]
    p7 = cte_amb["pa"]

    po6p = po6 * (1 - perte_it)
    p7p = p7 * (1 + perte_et)

    t7s = to6 * (p7p/po6p)**((y-1)/y)
    t7 = to6 - (to6-t7s)*ntp

    wpt = cpt * (to6 - t7)
    hp = wpt / 745.7

    station[7] = {"Po": p7, "To": t7, "wpt": wpt, "hp": hp}

    return p7, t7, wpt, hp

def main():
    print(f"{'Station':<10} | {'P (kPa)':<10} | {'T (K)':<10}")
    print("-" * 34)
    for i in range(1, 8):
        if i in station:
            print(f"{i:<10} | {station[i]['Po']:<10.2f} | {station[i]['To']:<10.2f}")
    print(f"\nSFC: {station[4]['sfc']:.6f} kg/Ns")
    print(f"Puissance de la turbine de puissance: {station[7]['hp']:.2f} HP")

if __name__ == "__main__":
    po1, to1 = Station_1()
    po2, to2, wlpc = Station_2(po1, to1)
    po3, to3, whpc = Station_3(po2, to2)
    po4, to4, mpr, sfc = Station_4(po3, to3)
    po5, to5, whpt, mpr = Station_5(po4, to4, whpc, mpr)
    po6, to6, wlpt = Station_6(po5, to5, mpr, wlpc)
    p7, t7, wpt, hp = Station_7(po6, to6)
    main()