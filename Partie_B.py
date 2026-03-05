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

# --- Données de cycle fictives (À Remplacer par vos résultats de la Partie A) ---
donnees_cycle = {
    'gamma': 1.31,                 # Air chaud HPT 
    'cp': 1214.0,                  # Air chaud HPT
    'T03': 1300.0,                 # Température de stagnation sortie HPT (K) - ESTIMATION
    'P03': 400000.0,               # Pression de stagnation sortie HPT (Pa) - ESTIMATION
    'm_dot': 5.443 * (1 + 0.02 - 0.10), # Débit massique (ajusté fuel et refroidissement) 
    'dh0_hpt': 350000.0,           # Travail spécifique requis par HPT (J/kg) - ESTIMATION PARTIE A
    'eta_hpt': 0.88                # Efficacité HPT 
}

geom = {}      # Pour stocker A, r_root, r_tip, r_moyen
vitesses = {}  # Pour stocker U, Va, Vw, etc.

def deg2rad(angle):
    return angle * np.pi / 180.0

def etape_1a(donnees_cycle):
    print("--- Étape 1.a : Établissement de la vitesse de rotation et géométrie ---")
    gamma = donnees_cycle['gamma']
    R = donnees_cycle['cp'] * (gamma - 1) / gamma
    
    # 1. Calcul des conditions statiques à la sortie
    T03 = donnees_cycle['T03']
    P03 = donnees_cycle['P03']
    M3 = contraintes['M3']
    
    T3 = T03 / (1 + ((gamma - 1) / 2) * M3**2)
    P3 = P03 / ((1 + ((gamma - 1) / 2) * M3**2)**(gamma / (gamma - 1)))
    rho3 = P3 / (R * T3)
    
    # Vitesse absolue V3 et sa composante axiale Va3
    V3 = M3 * np.sqrt(gamma * R * T3)
    alpha3_rad = deg2rad(contraintes['alpha_3'])
    Va3 = V3 * np.cos(alpha3_rad) 
    
    # 2. Calcul de la surface A3
    A3 = donnees_cycle['m_dot'] / (rho3 * Va3)
    geom['A3'] = A3
    geom['A2'] = A3 # Contrainte du cahier des charges : A2 = A3 [cite: 85]
    
    # 3. Choix du RPM (N) respectant le critère AN^2
    AN2_cible = (contraintes['AN2_min'] + contraintes['AN2_max']) / 2.0
    N_rpm = np.sqrt(AN2_cible / A3)
    omega = (N_rpm * np.pi) / 30.0
    
    # 4. Choix de U_root
    U_root_cible = (contraintes['U_emplanture_min'] + contraintes['U_emplanture_max']) / 2.0
    
    # 5. Détermination des rayons
    r_root = U_root_cible / omega
    r_tip = np.sqrt((A3 / np.pi) + r_root**2)
    r_moyen = (r_root + r_tip) / 2.0
    
    # 6. Vitesse d'entraînement au rayon moyen
    U_moyen = omega * r_moyen
    geom['h'] = r_tip - r_root # Hauteur de l'aube
    
    # Sauvegarde
    vitesses['N_rpm'] = N_rpm
    vitesses['omega'] = omega
    vitesses['U_moyen'] = U_moyen
    vitesses['Va3'] = Va3
    geom['r_root'] = r_root
    geom['r_tip'] = r_tip
    geom['r_moyen'] = r_moyen
    
    print(f"RPM: {N_rpm:.0f}, r_moyen: {r_moyen:.4f} m, U_moyen: {U_moyen:.2f} m/s")


def etape_1b_1c(donnees_cycle, Va2_guess):
    print("\n--- Étape 1.b et 1.c : Triangles des vitesses et Pertes (Rayon Moyen) ---")
    U_m = vitesses['U_moyen']
    Va3 = vitesses['Va3']
    Va2 = Va2_guess
    dh0 = donnees_cycle['dh0_hpt']
    
    alpha3 = deg2rad(contraintes['alpha_3'])
    
    # Vitesses sortie (Station 3)
    Vw3 = Va3 * np.tan(alpha3) # Convention
    Ww3 = Vw3 - U_m
    beta3 = np.arctan(Ww3 / Va3)
    
    # Euler pour trouver Vw2 : dh0 = U * (Vw2 - Vw3) (convention où les Vw s'ajoutent si sens opposés)
    Vw2 = (dh0 / U_m) + Vw3
    alpha2 = np.arctan(Vw2 / Va2)
    
    Ww2 = Vw2 - U_m
    beta2 = np.arctan(Ww2 / Va2)
    
    # Enregistrement Station 2 (Moyen)
    vitesses['V2'] = np.sqrt(Va2**2 + Vw2**2)
    vitesses['W2'] = np.sqrt(Va2**2 + Ww2**2)
    vitesses['W3'] = np.sqrt(Va3**2 + Ww3**2)
    vitesses['alpha2'] = alpha2
    vitesses['beta2'] = beta2
    vitesses['beta3'] = beta3
    vitesses['Vw2_m'] = Vw2 # Pour le vortex libre
    vitesses['Vw3_m'] = Vw3
    
    # Pertes (1.c) basées sur l'efficacité
    eta = donnees_cycle['eta_hpt']
    dh0_is = dh0 / eta
    perte_totale = dh0_is - dh0
    perte_stator = 0.40 * perte_totale # Hypothèse courante
    perte_rotor = 0.60 * perte_totale
    
    zeta_s = perte_stator / (0.5 * vitesses['V2']**2)
    zeta_r = perte_rotor / (0.5 * vitesses['W3']**2)
    
    print(f"Alpha 2: {np.degrees(alpha2):.2f}°, Beta 2: {np.degrees(beta2):.2f}°, Beta 3: {np.degrees(beta3):.2f}°")
    print(f"Zeta Stator: {zeta_s:.4f}, Zeta Rotor: {zeta_r:.4f}")


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
    
    # Calcul de Cx via le facteur de forme h/c. Approximation : c approx Cx / cos(gamma_stagger)
    # Pour simplifier en 1ere approx, disons h/Cx approx h/c (c axiale et vraie corde sont liées, on utilise h/Cx direct ici pour la démo)
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


def calcul():
    etape_1a(donnees_cycle)
    etape_1b_1c(donnees_cycle, Va2_guess=150.0) # Va2_guess est l'hypothèse demandée à l'étape 1c
    etape_2()
    etape_3()