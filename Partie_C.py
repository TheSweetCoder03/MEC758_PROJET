# Partie_C.py
import numpy as np
import Partie_A
import Partie_B
import matplotlib.pyplot as plt

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
    r_m1 = Partie_B.donnees['r_m'][1]
    r_m2 = Partie_B.donnees['r_m'][2]
    r_m3 = Partie_B.donnees['r_m'][3]

    # --- Vitesses de l'aubage hors-conception ---
    U1_hc = omega_hc * r_m1
    U2_hc = omega_hc * r_m2
    U3_hc = omega_hc * r_m3

    # --- Vitesses absolues inchangées (stator fixe, même débit, même thermo) ---
    Va1 = Partie_B.donnees['Va'][1]
    Va2 = Partie_B.donnees['Va'][2]
    Va3 = Partie_B.donnees['Va'][3]
    Vu1 = Partie_B.donnees['Vu'][1]
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
    donnees_hc['U']            = {1: U1_hc, 2: U2_hc, 3: U3_hc}
    donnees_hc['Va']           = {1: Va1,   2: Va2,   3: Va3}
    donnees_hc['Vu']           = {1: Vu1,   2: Vu2,   3: Vu3}
    donnees_hc['Vru']          = {2: Vru2_hc, 3: Vru3_hc}
    donnees_hc['Vr']           = {2: Vr2_hc,  3: Vr3_hc}
    donnees_hc['alpha_relatif']= {2: alpha_rel2_hc, 3: alpha_rel3_hc}
    donnees_hc['incidence']    = incidence

def tracer_triangles_vitesses():
    # 1. Extraction des variables du dictionnaire 'vitesses' (qui doit être global ici)
    U1 = donnees_hc['U'][1]
    U2 = donnees_hc['U'][2]
    U3 = donnees_hc['U'][3]
    
    Va1 = donnees_hc['Va'][1]
    Va2 = donnees_hc['Va'][2]
    Va3 = donnees_hc['Va'][3]

    Vu1 = donnees_hc['Vu'][1]
    Vu2 = donnees_hc['Vu'][2]
    Vu3 = donnees_hc['Vu'][3]
    
    # 2. Création de la figure avec 3 sous-graphiques alignés
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    
    # Fonction interne pour tracer un vecteur proprement
    def tracer_vecteur(ax, x_depart, y_depart, dx, dy, couleur, label, decalage_y=0):
        ax.quiver(x_depart, y_depart, dx, dy, angles='xy', scale_units='xy', scale=1, color=couleur, width=0.012)
        
        # Placement du texte aux 3/4 du vecteur (0.75) et retrait du fond blanc
        ax.text(x_depart + dx * 0.75, y_depart + dy * 0.75 + decalage_y, label, color=couleur, fontsize=12, fontweight='bold', ha='center', va='center')

    # --- STATION 1 : Entrée Stator ---
    # Vecteurs
    tracer_vecteur(axs[0], 0, Va1, 0, -Va1, 'black', 'Va1') # Ajout de Va1
    tracer_vecteur(axs[0], 0, Va1, Vu1, -Va1, 'blue', 'V1')
    tracer_vecteur(axs[0], 0, 0, U1, 0, 'green', 'U1 (Ref)', decalage_y=-10)
    axs[0].set_title("Station 1 (Entrée Stator)")
    
    # --- STATION 2 : Sortie Stator / Entrée Rotor ---
    # Vecteurs
    tracer_vecteur(axs[1], 0, Va2, 0, -Va2, 'black', 'Va2') # Ajout de Va2
    tracer_vecteur(axs[1], 0, Va2, Vu2, -Va2, 'blue', 'V2', decalage_y=15)            
    tracer_vecteur(axs[1], 0, Va2, Vu2 - U2, -Va2, 'red', 'Vr2', decalage_y=15)        
    tracer_vecteur(axs[1], Vu2 - U2, 0, U2, 0, 'green', 'U2', decalage_y=-10)         
    axs[1].set_title("Station 2 (Entrée Rotor)")
    
    # --- STATION 3 : Sortie Rotor ---
    # Vecteurs
    tracer_vecteur(axs[2], 0, Va3, 0, -Va3, 'black', 'Va3') # Ajout de Va3
    tracer_vecteur(axs[2], 0, Va3, Vu3, -Va3, 'blue', 'V3', decalage_y=15)            
    tracer_vecteur(axs[2], 0, Va3, Vu3 - U3, -Va3, 'red', 'Vr3', decalage_y=15)        
    tracer_vecteur(axs[2], Vu3 - U3, 0, U3, 0, 'green', 'U3', decalage_y=-10)         
    axs[2].set_title("Station 3 (Sortie Rotor)")

    # 3. Mise en forme et uniformisation des axes
    all_x = [0, U1, U2, U3, Vu1, Vu2, Vu3, Vu2-U2, Vu3-U3]
    all_y = [0, Va1, Va2, Va3]
    
    x_min, x_max = min(all_x) - 50, max(all_x) + 50
    y_min, y_max = -30, max(all_y) + 30

    for ax in axs:
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel("Vitesse Tangentielle (Vu, Vru, U) [m/s]")
        ax.set_ylabel("Vitesse Axiale (Va) [m/s]")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.axhline(0, color='black', linewidth=1.2) # Base du triangle
        ax.axvline(0, color='black', linewidth=1.2) # Axe Y principal

    plt.tight_layout()
    plt.show()

def pertes_incidence():

    # -----------------
    # Données
    # -----------------
    alpha_1 = Partie_B.donnees['alpha_relatif'][2]
    alpha_1_deg = np.degrees(alpha_1) 
    alpha_1_des = donnees_hc['alpha_relatif'][2]
    alpha_1_des_deg = np.degrees(alpha_1_des)
    incidence_deg = donnees_hc['incidence']
    d_c = 0.032 # Ratio - Leading eadge Diameter / Chord
    s_c = 0.56 # Ratio - Pitch / Chord
    beta_1 = Partie_B.donnees['alpha'][1] 
    beta_2 = Partie_B.donnees['alpha_relatif'][2]
    Yp_des = Partie_B.donnes['coefficient_perte'][1] # Coefficient de pertes de profil de AMDC
    Ys_des = Partie_B.donnes['coefficient_perte'][2] # Coefficient de pertes secondaire de AMDC

    # ---------------------------------
    # Calcul pour les pertes de profil
    # ----------------------------------
    d_s = d_c / s_c # Ratio - Leading edge diameter / Pitch
    ratio_cos_beta = np.cos(beta_1) / np.cos(beta_2) # Calcul du ratio pour simplifier les formules suivantes
    x_p = (d_s**-1.6)*(ratio_cos_beta**-2)*(alpha_1_deg - alpha_1_des_deg)

    if x_p > 0 :
        d_phi = 0.778 * (10**-5) * x_p + 0.56 * (10**-7) * (x_p**2) + 0.4 * (10**-10) * (x_p**3) + 2.054 * (10**-19) * (x_p**6)
    else :
        d_phi = -5.1734 * (10**-6) * x_p + 7.6902 * (10**-9) * (x_p**2 )
    
    Yp = Yp_des + d_phi

    # ---------------------------------
    # Calcul pour les pertes secondaires
    # ----------------------------------
    x_et = (incidence_deg / (180 - (beta_1 + beta_2))) * (ratio_cos_beta**-1.5) * (d_c**-0.3)

    if x_et > 0: 
        ratio_perte_sec = np.exp(0.9 * x_et) + 13 * (x_et**2) + 400 * (x_et**4)
    else:
        ratio_perte_sec = np.exp(0.9 * x_et)

    Ys = Ys_des * ratio_perte_sec