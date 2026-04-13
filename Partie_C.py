# Partie_C.py
import numpy as np
import Partie_A
import Partie_B

# Résultats hors-conception (rempli par partie_c())
donnees_hc = {}

def partie_c():
    print("\nÉTAPE C.1 : Triangles de vitesses hors-conception (RPM -20%)")

    # Réduction de 20% du RPM
    rpm_design = Partie_A.export_donnees_hpt()['rpmhpc']
    rpm_hc     = rpm_design * 0.80
    omega_hc   = (rpm_hc * np.pi) / 30.0

    print(f"RPM design       : {rpm_design:.0f} RPM")
    print(f"RPM hors-concept.: {rpm_hc:.0f} RPM (-20%)")

    # --- Géométrie figée (point de conception) ---
    r_m2 = Partie_B.donnees['r_m'][2]
    r_m3 = Partie_B.donnees['r_m'][3]

    # --- Vitesses de l'aubage hors-conception ---
    U2_hc = omega_hc * r_m2
    U3_hc = omega_hc * r_m3

    # --- Vitesses absolues inchangées (stator fixe, même débit, même thermo) ---
    Va2 = Partie_B.donnees['Va'][2]
    Va3 = Partie_B.donnees['Va'][3]
    Vu2 = Partie_B.donnees['Vu'][2]
    Vu3 = Partie_B.donnees['Vu'][3]

    # --- Nouveaux triangles relatifs ---
    Vru2_hc = Vu2 - U2_hc
    Vru3_hc = Vu3 - U3_hc
    Vr2_hc  = np.sqrt(Va2**2 + Vru2_hc**2)
    Vr3_hc  = np.sqrt(Va3**2 + Vru3_hc**2)

    alpha_rel2_hc = np.arctan(Vru2_hc / Va2)
    alpha_rel3_hc = np.arctan(Vru3_hc / Va3)

    # --- Incidence sur le rotor (entrée) ---
    alpha_rel2_design = Partie_B.donnees['alpha_relatif'][2]
    incidence = np.degrees(alpha_rel2_hc) - np.degrees(alpha_rel2_design)

    print(f"\nU2 design        : {Partie_B.donnees['U'][2]:.2f} m/s")
    print(f"U2 hors-concept. : {U2_hc:.2f} m/s")
    print(f"U3 design        : {Partie_B.donnees['U'][3]:.2f} m/s")
    print(f"U3 hors-concept. : {U3_hc:.2f} m/s")
    print(f"\nAngle relatif α_rel2 design      : {np.degrees(alpha_rel2_design):.2f}°")
    print(f"Angle relatif α_rel2 hors-concept: {np.degrees(alpha_rel2_hc):.2f}°")
    print(f"Incidence sur le rotor (i)       : {incidence:.2f}°")

    # --- Sauvegarde ---
    donnees_hc['rpm']          = rpm_hc
    donnees_hc['omega']        = omega_hc
    donnees_hc['U']            = {2: U2_hc, 3: U3_hc}
    donnees_hc['Va']           = {2: Va2,   3: Va3}
    donnees_hc['Vu']           = {2: Vu2,   3: Vu3}
    donnees_hc['Vru']          = {2: Vru2_hc, 3: Vru3_hc}
    donnees_hc['Vr']           = {2: Vr2_hc,  3: Vr3_hc}
    donnees_hc['alpha_relatif']= {2: alpha_rel2_hc, 3: alpha_rel3_hc}
    donnees_hc['incidence']    = incidence
