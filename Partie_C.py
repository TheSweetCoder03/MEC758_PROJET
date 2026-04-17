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
    r_m2 = Partie_B.donnees['r_m'][2]
    r_m3 = Partie_B.donnees['r_m'][3]

    # --- Vitesses de l'aubage hors-conception ---
    U2_hc = omega_hc * r_m2
    U3_hc = omega_hc * r_m3

    a1 = Partie_B.donnees['alpha'][1]

    # --- Vitesses absolues inchangées (stator fixe, même débit, même thermo) ---
    Va1 = Partie_B.donnees['Va'][1]
    Va2 = Partie_B.donnees['Va'][2]
    Va3 = Partie_B.donnees['Va'][3]
    Vu2 = Partie_B.donnees['Vu'][1]
    Vu3 = Partie_B.donnees['Vu'][2]

    # --- Nouveaux triangles relatifs ---
    Vru2_hc = Vu2 - U2_hc
    Vru3_hc = U3_hc - Vu3
    Vr2_hc  = np.sqrt(Va2**2 + Vru2_hc**2)
    Vr3_hc  = np.sqrt(Va3**2 + Vru3_hc**2)

    alpha_rel2_hc = np.arctan(Vru2_hc / Va2)
    alpha_rel3_hc = np.arctan(Vru3_hc / Va3)

    # --- Incidence sur le rotor (entrée) ---
    alpha_rel2_design = Partie_B.donnees['alpha_relatif'][2]
    incidence = np.degrees(alpha_rel2_hc) - np.degrees(alpha_rel2_design)

    print(f"\nU2 design        : {Partie_B.donnees['U'][1]:.2f} m/s")
    print(f"U2 hors-concept. : {U2_hc:.2f} m/s")
    print(f"U3 design        : {Partie_B.donnees['U'][2]:.2f} m/s")
    print(f"U3 hors-concept. : {U3_hc:.2f} m/s")
    print(f"\nAngle de métal beta2 de design      : {np.degrees(alpha_rel2_design):.2f}°")
    print(f"Angle relatif α_rel2 hors-concept: {np.degrees(alpha_rel2_hc):.2f}°")
    print(f"Incidence sur le rotor (i)       : {incidence:.2f}°")

    # --- Sauvegarde ---
    donnees_hc['alpha']        = {1: a1}
    donnees_hc['rpm']          = rpm_hc
    donnees_hc['omega']        = omega_hc
    donnees_hc['U']            = {2: U2_hc, 3: U3_hc}
    donnees_hc['Va']           = {1: Va1,   2: Va2,   3: Va3}
    donnees_hc['Vu']           = {1: Vu2,   2: Vu3}
    donnees_hc['Vru']          = {2: Vru2_hc, 3: Vru3_hc}
    donnees_hc['Vr']           = {2: Vr2_hc,  3: Vr3_hc}
    donnees_hc['alpha_relatif']= {2: alpha_rel2_hc, 3: alpha_rel3_hc}
    donnees_hc['incidence']    = incidence

def tracer_triangles_vitesses():
    a1 = donnees_hc['alpha'][1]
    # 1. Extraction des variables du dictionnaire 'vitesses' (qui doit être global ici)
    U2 = donnees_hc['U'][2]
    U3 = donnees_hc['U'][3]
    
    Va1 = donnees_hc['Va'][1]
    Va2 = donnees_hc['Va'][2]
    Va3 = donnees_hc['Va'][3]

    Vu2 = donnees_hc['Vu'][1]
    Vu3 = donnees_hc['Vu'][2]
    
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
    tracer_vecteur(axs[0], 0, Va1, Va1*np.tan(a1), -Va1, 'blue', 'V1')
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
    all_x = [0, U2, U3, Vu2, Vu3, Vu2-U2, Vu3-U3]
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
    #Calcul des pertes d'incidence avec la méthode de corrélation de Moustapha

    # -----------------
    # Données
    # -----------------
    alpha_1 = donnees_hc['alpha_relatif'][2]
    alpha_1_deg = np.degrees(alpha_1) 
    alpha_1_des = Partie_B.donnees['alpha_relatif'][2]
    alpha_1_des_deg = np.degrees(alpha_1_des)
    incidence_deg = donnees_hc['incidence']
    d_c = 0.053 # Ratio - Leading eadge Diameter / Chord
    s_c = 0.75 # Ratio - Pitch / Chord
    beta_1 = Partie_B.donnees['alpha_relatif'][2]
    beta_1_deg = np.degrees(beta_1) 
    beta_2 = Partie_B.donnees['alpha_relatif'][3]
    beta_2_deg = np.degrees(beta_2) 
    Yp_des = Partie_B.donnees['coefficient_perte_rotor'][1] # Coefficient de pertes de profil de AMDC on design du rotor
    Ys_des = Partie_B.donnees['coefficient_perte_rotor'][2] # Coefficient de pertes secondaire de AMDC on design du rotor
    Ytet_r = Partie_B.donnees['coefficient_perte_rotor'][3] # Coefficient de pertes traling edge de AMDC on design du rotor
    Ytc_r = Partie_B.donnees['coefficient_perte_rotor'][4] # Coefficient de pertes tip clearance de AMDC on design du rotor

    # ---------------------------------
    # Calcul pour les pertes de profil
    # ----------------------------------
    d_s = d_c / s_c # Ratio - Leading edge diameter / Pitch
    ratio_cos_beta = np.cos(beta_1) / np.cos(beta_2) # Calcul du ratio pour simplifier les formules suivantes
    x_p = (d_s**-1.6)*(ratio_cos_beta**-2)*(alpha_1_deg - alpha_1_des_deg)

    if x_p > 0 :
        d_phi = 0.778e-5 * x_p + 0.56e-7 * (x_p**2) + 0.4e-10 * (x_p**3) + 2.054e-19 * (x_p**6)
    else :
        d_phi = -5.1734e-6 * x_p + 7.6902e-9 * (x_p**2 )
    
    Yp_r = Yp_des + d_phi

    # ---------------------------------
    # Calcul pour les pertes secondaires
    # ----------------------------------
    x_et = (incidence_deg / (180 - (beta_1_deg + beta_2_deg))) * (ratio_cos_beta**-1.5) * (d_c**-0.3)

    if x_et > 0: 
        ratio_perte_sec = np.exp(0.9 * x_et) + 13 * (x_et**2) + 400 * (x_et**4)
    else:
        ratio_perte_sec = np.exp(0.9 * x_et)

    Ys_r = Ys_des * ratio_perte_sec

    Y_rotor = Yp_r + Ys_r + Ytet_r + Ytc_r

    donnees_hc['Y_rotor']= Y_rotor

    print(f"\nLes pertes de profil d'incidence sont de : {Yp_r:.2f}")
    print(f"Le delta_phi pour profil est de : {d_phi:.4f}")
    print(f"Les pertes secondaires d'incidence sont de : {Ys_r:.2f}")
    print(f"Le ratio de perte pour secondaire est de : {ratio_perte_sec:.2f}")
    print(f"Les pertes totales au rotor (Y_R) sont : {Y_rotor:.4f}")

def rendement_incidence(donnees_hpt):
    #Calcul du rendement de l'étage avec la réduction du RPM (et donc les pertes d'incidences)

    # ------------
    # Données 
    # ------------
    T2 = Partie_B.donnees['T'][2]
    T3 = Partie_B.donnees['T'][3]
    R_gaz = Partie_B.donnees['R_gaz']
    gamma = donnees_hpt['gamma']
    U2_hc = donnees_hc['U'][2]
    U3_hc = donnees_hc['U'][3]
    Vu2 = donnees_hc['Vu'][1]
    Vu3 = donnees_hc['Vu'][2]
    Va2 = donnees_hc['Va'][2]
    Vr3 = donnees_hc['Vr'][3]
    T01 = Partie_B.donnees['To'][1]
    cp = donnees_hpt['cp']
    Y_rotor = donnees_hc['Y_rotor']
    Y_stator = Partie_B.donnees['Y_perte_stator']


    # ------------
    # Calculs
    # ------------ 
    V2 = np.sqrt((Va2**2) + (Vu2**2))
    M2 = V2 / np.sqrt(gamma * R_gaz * T2)
    lambda_stator = Y_stator / (1 + 0.5 * gamma * M2**2)

    Mr3 = Vr3 / np.sqrt(gamma * R_gaz * T3)
    lambda_rotor = Y_rotor / (1 + 0.5 * gamma * Mr3**2)

    #On calcul le nouveau T03 pour prendre en compte la diminution du travail d'Euler
    w_hc = U2_hc * Vu2 - U3_hc * Vu3
    T03_hc = T01 - w_hc / cp


    rendement = 1 / (1 + (lambda_stator * V2**2 + lambda_rotor * Vr3**2) / (2 * cp * (T01 - T03_hc)))

    print(f"\nLe rendement de l'étage est de : {rendement:.4f}")

    

