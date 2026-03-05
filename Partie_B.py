import numpy as np

# --- Constantes et Contraintes du Projet ---
contraintes = {
    # Paramètres de l'étage
    'M1': 0.14,                    # Mach entrée 
    'M3': 0.37,                    # Mach sortie 
    'alpha_1': -10.0,              # Angle absolu entrée stator 
    'alpha_3': 22.0,               # Angle absolu sortie rotor 
    'reaction': 0.63,              # Degré de réaction 
    'AN2_min': 1.6129e7,           # Surface fois vitesse au carré (converti en m^2 RPM^2) 
    'AN2_max': 3.2258e7,           
    'U_emplanture_min': 335.28,    # Vitesse min aubage à l'emplanture (converti en m/s) 
    'U_emplanture_max': 365.76,    
    
    # Paramètres de l'aube fixe (Stator)
    'stator_h_c': 0.70,            # Facteur de forme (h/c)
    'stator_zweifel': 0.75,        # Coefficient de Zweifel (moyen)
    'stator_te_thick': 0.001016,   # Épaisseur au bord de fuite (m)

     # Paramètres de l'aube mobile (Rotor)
    'rotor_h_c': 1.30,             # Facteur de forme (h/c)
    'rotor_zweifel': 0.90,         # Coefficient de Zweifel (moyen)
    'rotor_te_thick': 0.000508,    # Épaisseur au bord de fuite (m)
    'vie_heures': 300              # Durée de vie 
}

geom = {}      # Pour stocker A, r_root, r_tip, r_moyen
vitesses = {}  # Pour stocker U, Va, Vw, etc.

def deg2rad(angle):
    return angle * np.pi / 180.0

def etape_1(donnees_hpt, racine_constante=True, Va2_guess=150.0):
    print(f"\nÉTAPES 1.a, 1.b, 1.c : Géométrie et Triangles")
    print(f"Stratégie de veine : {'Racine Constante' if racine_constante else 'Bout (Tip) Constant'}")
    
    # Récupération des données
    gamma = donnees_hpt['gamma']
    R_gaz = donnees_hpt['cp'] * (gamma - 1) / gamma
    m_dot = donnees_hpt['m_dot']
    dh0 = donnees_hpt['dh0_hpt']
    cp = donnees_hpt['cp']
    
    T01 = donnees_hpt['T01']
    P01 = donnees_hpt['P01']
    T03 = donnees_hpt['T03']
    P03 = donnees_hpt['P03']

    # 1. STATION 3 (Sortie Rotor - Base de la géométrie)
    M3 = contraintes['M3']
    T3 = T03 / (1 + ((gamma - 1) / 2) * M3**2)
    P3 = P03 / ((1 + ((gamma - 1) / 2) * M3**2)**(gamma / (gamma - 1)))
    rho3 = P3 / (R_gaz * T3)

    V3 = M3 * np.sqrt(gamma * R_gaz * T3)
    alpha3_rad = deg2rad(contraintes['alpha_3'])
    Va3 = V3 * np.cos(alpha3_rad)
    Vw3 = V3 * np.sin(alpha3_rad)

    A3 = m_dot / (rho3 * Va3)

    # Vitesse de rotation (N) basée sur Station 3
    AN2_cible = (contraintes['AN2_min'] + contraintes['AN2_max']) / 2.0
    N_rpm = np.sqrt(AN2_cible / A3)
    omega = (N_rpm * np.pi) / 30.0
    U_root_cible = (contraintes['U_emplanture_min'] + contraintes['U_emplanture_max']) / 2.0

    r_root3 = U_root_cible / omega
    r_tip3 = np.sqrt((A3 / np.pi) + r_root3**2)
    r_m3 = (r_root3 + r_tip3) / 2.0
    U3 = omega * r_m3

    # 2. STATION 1 (Entrée Stator)
    M1 = contraintes['M1']
    T1 = T01 / (1 + ((gamma - 1) / 2) * M1**2)
    P1 = P01 / ((1 + ((gamma - 1) / 2) * M1**2)**(gamma / (gamma - 1)))
    rho1 = P1 / (R_gaz * T1)

    V1 = M1 * np.sqrt(gamma * R_gaz * T1)
    alpha1_rad = deg2rad(contraintes['alpha_1'])
    Va1 = V1 * np.cos(alpha1_rad)

    A1 = m_dot / (rho1 * Va1)

    if racine_constante:
        r_root1 = r_root3
        r_tip1 = np.sqrt((A1 / np.pi) + r_root1**2)
    else:
        r_tip1 = r_tip3
        r_root1 = np.sqrt(r_tip1**2 - (A1 / np.pi))
        
    r_m1 = (r_root1 + r_tip1) / 2.0
    U1 = omega * r_m1

    # 3. STATION 2 (Solveur Itératif pour Va2 et géométrie)
    Va2 = Va2_guess
    r_m2 = r_m3 # Hypothèse de départ
    
    for _ in range(20): # Boucle de convergence
        U2 = omega * r_m2
        
        # Équation d'Euler ajustée pour U2 != U3
        Vw2 = (dh0 + U3 * Vw3) / U2
        V2 = np.sqrt(Va2**2 + Vw2**2)
        
        # Thermodynamique Station 2
        T2 = T01 - (V2**2) / (2 * cp)
        P2 = P01 * (T2 / T01)**(gamma / (gamma - 1)) # Approx isentropique pour la géométrie
        rho2 = P2 / (R_gaz * T2)
        
        A2 = m_dot / (rho2 * Va2)
        
        r_m2_old = r_m2
        if racine_constante:
            r_root2 = r_root3
            r_tip2 = np.sqrt((A2 / np.pi) + r_root2**2)
        else:
            r_tip2 = r_tip3
            r_root2 = np.sqrt(r_tip2**2 - (A2 / np.pi))
            
        r_m2 = (r_root2 + r_tip2) / 2.0
        
        # Condition de sortie : la géométrie ne bouge plus
        if abs(r_m2 - r_m2_old) < 1e-6:
            break

    # 4. Finalisation des Triangles et Pertes
    alpha2 = np.arctan(Vw2 / Va2)

    Ww2 = Vw2 - U2
    beta2 = np.arctan(Ww2 / Va2)
    W2 = np.sqrt(Va2**2 + Ww2**2)

    Ww3 = Vw3 - U3
    beta3 = np.arctan(Ww3 / Va3)
    W3 = np.sqrt(Va3**2 + Ww3**2)

    # Calcul des pertes (1.c)
    eta_hpt = donnees_hpt['eta_hpt']
    dh0_is = dh0 / eta_hpt
    perte_totale = dh0_is - dh0
    zeta_s = (0.40 * perte_totale) / (0.5 * V2**2)
    zeta_r = (0.60 * perte_totale) / (0.5 * W3**2)

    # Calcul du degré de réaction RÉEL final
    Reaction = (cp * (T2 - T3)) / dh0

    # Sauvegarde structurée dans les dictionnaires
    geom['A'] = {1: A1, 2: A2, 3: A3}
    geom['r_root'] = {1: r_root1, 2: r_root2, 3: r_root3}
    geom['r_tip'] = {1: r_tip1, 2: r_tip2, 3: r_tip3}
    geom['r_m'] = {1: r_m1, 2: r_m2, 3: r_m3}
    geom['h'] = {1: r_tip1-r_root1, 2: r_tip2-r_root2, 3: r_tip3-r_root3}

    vitesses['N_rpm'] = N_rpm
    vitesses['omega'] = omega
    vitesses['U'] = {1: U1, 2: U2, 3: U3}
    vitesses['Va'] = {1: Va1, 2: Va2, 3: Va3}
    vitesses['Vw'] = {2: Vw2, 3: Vw3}
    vitesses['alpha'] = {1: alpha1_rad, 2: alpha2, 3: alpha3_rad}
    vitesses['beta'] = {2: beta2, 3: beta3}

    # Affichage des résultats
    print(f"Régime : {N_rpm:.0f} RPM")
    print(f"Rayons moyens [m] : r_m1={r_m1:.4f} | r_m2={r_m2:.4f} | r_m3={r_m3:.4f}")
    print(f"Vitesses U [m/s]  : U1={U1:.2f} | U2={U2:.2f} | U3={U3:.2f}")
    print(f"Vitesses Va [m/s] : Va1={Va1:.2f} | Va2={Va2:.2f} | Va3={Va3:.2f}")
    print(f"Angles [deg]      : Alpha2={np.degrees(alpha2):.2f}° | Beta2={np.degrees(beta2):.2f}° | Beta3={np.degrees(beta3):.2f}°")
    print(f"Pertes            : Zeta_S={zeta_s:.4f} | Zeta_R={zeta_r:.4f}")
    print(f"Degré de réaction : {Reaction:.3f}")


def etape_2():
    print("\n--- Étape 2 : Triangles des vitesses (Vortex Libre) ---")
    # Constantes du vortex libre : r * Vw = Cte, Va = Cte
    r_m = geom['r_moyen']
    K2 = r_m * vitesses['Vw2_m']
    K3 = r_m * vitesses['Vw3_m']
    
    rayons = {'Root': geom['r_root'], 'Tip': geom['r_tip']}
    
    for nom, r in rayons.items():
        U = vitesses['omega'] * r
        Vw2 = K2 / r
        Vw3 = K3 / r
        
        # En conservant Va2 et Va3 constants radialement (hypothèse de base)
        Va2 = vitesses['V2'] * np.cos(vitesses['alpha2']) # = Va2_guess
        Va3 = vitesses['Va3']
        
        alpha2 = np.arctan(Vw2 / Va2)
        beta2 = np.arctan((Vw2 - U) / Va2)
        beta3 = np.arctan((Vw3 - U) / Va3)
        
        # Degré de réaction : R = (Va / 2U) * (tan(beta3) - tan(beta2)) # Approx
        R = 1 - (Va2 / (2*U)) * (np.tan(alpha2) + (Va3/Va2)*np.tan(contraintes['alpha_3'])) 
        
        print(f"Position {nom} (r={r:.4f}m) : U={U:.2f} m/s | Alpha2={np.degrees(alpha2):.2f}° | Beta2={np.degrees(beta2):.2f}° | Reaction={R:.2f}")


def etape_3():
    print("\n--- Étape 3 : Paramètres des aubes (Zweifel) ---")
    # --- Rotor ---
    # Zweifel Rotor: Z_R = 2 * (s/Cx) * cos^2(beta3) * (tan(beta2) + tan(beta3))
    beta2_m = vitesses['beta2']
    beta3_m = abs(vitesses['beta3']) # Prendre la valeur absolue pour la géométrie
    Z_R = contraintes['rotor_zweifel']
    
    s_cx_rotor = Z_R / (2 * (np.cos(beta3_m)**2) * (np.tan(beta2_m) + np.tan(beta3_m)))
    
    # Calcul de Cx via le facteur de forme h/c. 
    Cx_rotor = geom['h'] / contraintes['rotor_h_c'] 
    pas_rotor = s_cx_rotor * Cx_rotor
    N_rotor = 2 * np.pi * geom['r_moyen'] / pas_rotor
    
    print(f"Rotor  -> s/Cx: {s_cx_rotor:.3f}, Corde axiale: {Cx_rotor*100:.2f} cm, Nombre d'aubes: {int(np.ceil(N_rotor))} (arrondi)")

    # --- Stator ---
    alpha1_m = deg2rad(contraintes['alpha_1'])
    alpha2_m = vitesses['alpha2']
    Z_S = contraintes['stator_zweifel']
    
    s_cx_stator = Z_S / (2 * (np.cos(alpha2_m)**2) * (np.tan(alpha2_m) - np.tan(alpha1_m)))
    Cx_stator = geom['h'] / contraintes['stator_h_c']
    pas_stator = s_cx_stator * Cx_stator
    N_stator = 2 * np.pi * geom['r_moyen'] / pas_stator
    
    print(f"Stator -> s/Cx: {s_cx_stator:.3f}, Corde axiale: {Cx_stator*100:.2f} cm, Nombre d'aubes: {int(np.ceil(N_stator))} (arrondi)")


def calcul(donnees_hpt):
    # Adaptation des variables du dictionnaire de la Partie A pour la Partie B
    donnees_hpt['dh0_hpt'] = donnees_hpt['W_hpt'] / donnees_hpt['m_dot'] # Travail spécifique (J/kg)
    donnees_hpt['eta_hpt'] = donnees_hpt['eta_iso']
    
    etape_1(donnees_hpt, racine_constante=True, Va2_guess=150.0)
    etape_2()
    etape_3()