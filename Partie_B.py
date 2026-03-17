import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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
vitesses = {}  # Pour stocker U, Va, Vu, etc.

def deg2rad(angle):
    return angle * np.pi / 180.0

def etape_1(donnees_hpt, racine_constante, Va2_guess=150.0, tolerance = 1e-6):
    print(f"\nÉTAPES 1.a, 1.b, 1.c : Géométrie et Triangles")
    print(f"Stratégie de veine : {'Racine Constante' if racine_constante else 'Bout (Tip) Constant'}")
    
    # Récupération des données
    gamma = donnees_hpt['gamma']
    R_gaz = donnees_hpt['cp'] * (gamma - 1) / gamma
    m_dot = donnees_hpt['m_dot']
    dh0 = donnees_hpt['dh0_hpt']
    cp = donnees_hpt['cp']
    
    T01 = donnees_hpt['T04']
    P01 = donnees_hpt['P04']
    T03 = donnees_hpt['T05']
    P03 = donnees_hpt['P05']

    # STATION 3 (Sortie Rotor)
    M3 = contraintes['M3']
    T3 = T03 / (1 + ((gamma - 1) / 2) * M3**2)
    P3 = P03 * ((T3/T03)**(gamma / (gamma - 1)))
    rho3 = P3 / (R_gaz * T3)

    V3 = M3 * np.sqrt(gamma * R_gaz * T3)
    Va3 = V3 * np.cos(deg2rad(contraintes['alpha_3']))
    Vu3 = V3 * np.sin(deg2rad(contraintes['alpha_3']))

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
    
    def equation_rm2(rm2_hypothese):
        # fsolve passe un tableau (array), on extrait la valeur
        r_in = rm2_hypothese[0] 
        
        U2_temp = omega * r_in
        
        # Équation d'Euler ajustée
        Vu2_temp = (dh0 + U3 * Vu3) / U2_temp
        V2_temp = np.sqrt(Va2**2 + Vu2_temp**2)
        
        # Thermodynamique Station 2
        T2_temp = T01 - (V2_temp**2) / (2 * cp)
        P2_temp = P01 * (T2_temp / T01)**(gamma / (gamma - 1))
        rho2_temp = P2_temp / (R_gaz * T2_temp)
        
        A2_temp = m_dot / (rho2_temp * Va2)
        
        if racine_constante:
            r_root2_temp = r_root3
            r_tip2_temp = np.sqrt((A2_temp / np.pi) + r_root2_temp**2)
        else:
            r_tip2_temp = r_tip3
            r_root2_temp = np.sqrt(r_tip2_temp**2 - (A2_temp / np.pi))
            
        r_out = (r_root2_temp + r_tip2_temp) / 2.0
        
        # Le solveur cherche la valeur où cette différence est exactement 0
        return r_out - r_in

    # Résolution : on donne la fonction et une valeur de départ (x0 = r_m3)
    rm2_solution = fsolve(equation_rm2, x0=[r_m2], xtol=tolerance)
    r_m2 = rm2_solution[0]
    
    # Recalcul final des variables
    U2 = omega * r_m2
    Vu2 = (dh0 + U3 * Vu3) / U2
    V2 = np.sqrt(Va2**2 + Vu2**2)
    T2 = T01 - (V2**2) / (2 * cp)
    P2 = P01 * (T2 / T01)**(gamma / (gamma - 1))
    rho2 = P2 / (R_gaz * T2)
    A2 = m_dot / (rho2 * Va2)
    
    if racine_constante:
        r_root2 = r_root3
        r_tip2 = np.sqrt((A2 / np.pi) + r_root2**2)
    else:
        r_tip2 = r_tip3
        r_root2 = np.sqrt(r_tip2**2 - (A2 / np.pi))

    # 4. Finalisation des Triangles et Pertes
    alpha2 = np.arctan(Vu2 / Va2)

    Ww2 = Vu2 - U2
    beta2 = np.arctan(Ww2 / Va2)
    W2 = np.sqrt(Va2**2 + Ww2**2)

    Ww3 = Vu3 - U3
    beta3 = np.arctan(Ww3 / Va3)
    W3 = np.sqrt(Va3**2 + Ww3**2)

    # Calcul des pertes (1.c)
    eta_hpt = donnees_hpt['eta_iso']
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
    vitesses['Vu'] = {2: Vu2, 3: Vu3}
    vitesses['alpha'] = {1: contraintes['alpha_1'], 2: alpha2, 3: contraintes['alpha_3']}
    vitesses['beta'] = {2: beta2, 3: beta3}

    # Affichage des résultats
    print(f"Régime : {N_rpm:.0f} RPM")
    print(f"Rayons moyens [m] : r_m1={r_m1:.4f} | r_m2={r_m2:.4f} | r_m3={r_m3:.4f}")
    print(f"Vitesses U [m/s]  : U1={U1:.2f} | U2={U2:.2f} | U3={U3:.2f}")
    print(f"Vitesses Va [m/s] : Va1={Va1:.2f} | Va2={Va2:.2f} | Va3={Va3:.2f}")
    print(f"Angles [deg]      : Alpha2={np.degrees(alpha2):.2f}° | Beta2={np.degrees(beta2):.2f}° | Beta3={np.degrees(beta3):.2f}°")
    print(f"Pertes            : Zeta_S={zeta_s:.4f} | Zeta_R={zeta_r:.4f}")
    print(f"Degré de réaction : {Reaction:.3f}")

def plot_geometrie_turbine(geom_dict):
    # Extraction des données du dictionnaire
    stations = [1, 2, 3]
    # On définit des positions axiales arbitraires pour la visualisation
    x = [0, 1, 2] 
    
    r_root = [geom_dict['r_root'][s] for s in stations]
    r_tip = [geom_dict['r_tip'][s] for s in stations]
    r_m = [geom_dict['r_m'][s] for s in stations]

    plt.figure(figsize=(10, 6))

    # Tracer les limites de la veine
    plt.plot(x, r_tip, 'k-', linewidth=2, label='Bout (Tip)')
    plt.plot(x, r_root, 'k-', linewidth=2, label='Emplanture (Root)')
    
    # Tracer le rayon moyen
    plt.plot(x, r_m, 'r--', alpha=0.7, label='Rayon moyen (Mean)')

    # Remplissage de la zone de flux (la veine)
    plt.fill_between(x, r_root, r_tip, color='skyblue', alpha=0.2, label='Veine fluide')

    # Ajout des points aux stations
    plt.scatter(x, r_tip, color='black')
    plt.scatter(x, r_root, color='black')
    plt.scatter(x, r_m, color='red', s=20)

    # Décoration du graphique
    plt.title("Profil méridien de l'étage de turbine", fontsize=14)
    plt.xlabel("Stations (Positions axiales relatives)", fontsize=12)
    plt.ylabel("Rayon [m]", fontsize=12)
    
    # Configuration des axes
    plt.xticks(x, ['Station 1\n(Entrée Stator)', 'Station 2\n(Interface)', 'Station 3\n(Sortie Rotor)'])
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='best')
    
    # Ajuster les limites pour mieux voir
    margin = (max(r_tip) - min(r_root)) * 0.2
    plt.ylim(min(r_root) - margin, max(r_tip) + margin)
    
    plt.tight_layout()
    plt.show()


def calcul(donnees_hpt):
    
    etape_1(donnees_hpt, racine_constante=False, Va2_guess=150.0)
    plot_geometrie_turbine(geom)