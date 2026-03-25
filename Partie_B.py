import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

contraintes = {
    # Paramètres de l'étage
    'M1': 0.14,                    # Mach entrée 
    'M3': 0.37,                    # Mach sortie 
    'alpha_1': -10.0,              # Angle absolu entrée stator 
    'alpha_3': 22.0,               # Angle absolu sortie rotor 
    'reaction': 0.63,              # Degré de réaction 
    'AN2_min': 1.6129e7,           # Surface fois vitesse au carré (converti en m^2 RPM) 
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
pertes = {'stator': 0.4,'rotor' : 0.6}

def deg2rad(angle):
    return angle * np.pi / 180.0

def etape_1(donnees_hpt, racine_constante, tolerance=1e-6):
    print(f"\nÉTAPES 1.a, 1.b, 1.c : Géométrie et Triangles")
    print(f"Stratégie de veine : {'Racine Constante' if racine_constante else 'Bout (Tip) Constant'}")
    
    # Récupération des données
    gamma = donnees_hpt['gamma']
    R_gaz = donnees_hpt['cp'] * (gamma - 1) / gamma
    m_dot = donnees_hpt['m_dot']
    dh0 = donnees_hpt['dh0_hpt']
    cp = donnees_hpt['cp']
    rpm = donnees_hpt["rpmhpc"]
    T01 = donnees_hpt['T04']
    P01 = donnees_hpt['P04']
    T03 = donnees_hpt['T05']
    P03 = donnees_hpt['P05']
    
    # Calcul données universel
    U_root_cible = (contraintes['U_emplanture_min'] + contraintes['U_emplanture_max']) / 2.0
    omega = (rpm * np.pi) / 30.0

    # Station 3 (Sortie Rotor) 
    M3 = contraintes['M3']
    T3 = T03 / (1 + ((gamma - 1) / 2) * M3**2)
    P3 = P03 * ((T3/T03)**(gamma / (gamma - 1)))
    rho3 = P3 / (R_gaz * T3)

    V3 = M3 * np.sqrt(gamma * R_gaz * T3)
    Va3 = V3 * np.cos(deg2rad(contraintes['alpha_3']))
    Vu3 = V3 * np.sin(deg2rad(contraintes['alpha_3']))

    A3 = m_dot / (rho3 * Va3)
    r_root3 = U_root_cible / omega
    r_tip3 = np.sqrt((A3 / np.pi) + r_root3**2)
    r_m3 = (r_root3 + r_tip3) / 2.0
    U3 = omega * r_m3

    # Station 1 (Entrée Stator) 
    M1 = contraintes['M1']
    T1 = T01 / (1 + ((gamma - 1) / 2) * M1**2)
    P1 = P01 / ((1 + ((gamma - 1) / 2) * M1**2)**(gamma / (gamma - 1)))
    rho1 = P1 / (R_gaz * T1)

    V1 = M1 * np.sqrt(gamma * R_gaz * T1)
    alpha1_rad = deg2rad(contraintes['alpha_1'])
    Va1 = V1 * np.cos(alpha1_rad)
    Vu1 = V1 * np.sin(deg2rad(contraintes['alpha_1']))

    A1 = m_dot / (rho1 * Va1)
    if racine_constante:
        r_root1 = r_root3
        r_tip1 = np.sqrt((A1 / np.pi) + r_root1**2)
    else:
        r_tip1 = r_tip3
        r_root1 = np.sqrt(r_tip1**2 - (A1 / np.pi))  
    r_m1 = (r_root1 + r_tip1) / 2.0
    U1 = omega * r_m1

    # Calcul des pertes totales
    eta_hpt = donnees_hpt['eta_iso']
    dh0_is = dh0 / eta_hpt
    perte_totale = dh0_is - dh0

    # Thermodynamique de la station 2 (fixée par la réaction)
    T2 = T3 + (contraintes['reaction'] * dh0) / cp
    V2 = np.sqrt(2 * cp * (T01 - T2))
    
    # Pertes au stator (Hypothèse 40%)
    lambda_N = (pertes['stator'] * perte_totale) / (0.5 * V2**2)
    
    r_m2_guess = r_m3
    
    def equation_rm2(rm2_hypothese):
        r_in = rm2_hypothese[0] 
        U2_temp = omega * r_in
        
        Vu2_temp = (dh0 + U3 * Vu3) / U2_temp
        Va2_temp = np.sqrt(V2**2 - Vu2_temp**2)
        
        M2_temp = V2 / np.sqrt(gamma * R_gaz * T2)
        
        YN_temp = lambda_N * (1 + 0.5 * gamma * M2_temp**2)
        
        PR2_temp = (1 + ((gamma - 1) / 2) * M2_temp**2)**(gamma / (gamma - 1))
        
        # 4. Calcul de P2_temp en isolant P2 de l'équation de Y_N
        P2_temp = P01 / (PR2_temp + YN_temp * (PR2_temp - 1))
        
        # 5. Densité réelle
        rho2_temp = P2_temp / (R_gaz * T2)
        
        A2_temp = m_dot / (rho2_temp * Va2_temp)
        if racine_constante:
            r_root2_temp = r_root3
            r_tip2_temp = np.sqrt((A2_temp / np.pi) + r_root2_temp**2)
        else:
            r_tip2_temp = r_tip3
            r_root2_temp = np.sqrt(r_tip2_temp**2 - (A2_temp / np.pi))
            
        r_out = (r_root2_temp + r_tip2_temp) / 2.0
        
        return r_out - r_in

    rm2_solution = fsolve(equation_rm2, x0=[r_m2_guess], xtol=tolerance)
    r_m2 = rm2_solution[0]
    
    # Calcul des variables avec la solution trouvée
    U2 = omega * r_m2
    Vu2 = (dh0 + U3 * Vu3) / U2
    Va2 = np.sqrt(V2**2 - Vu2**2)
    M2 = V2 / np.sqrt(gamma * R_gaz * T2)
    Y_N = lambda_N * (1 + 0.5 * gamma * M2**2)
    PR2 = (1 + ((gamma - 1) / 2) * M2**2)**(gamma / (gamma - 1))
    P2 = P01 / (PR2 + Y_N * (PR2 - 1))
    rho2 = P2 / (R_gaz * T2)
    A2 = m_dot / (rho2 * Va2)
    if racine_constante:
        r_root2 = r_root3
        r_tip2 = np.sqrt((A2 / np.pi) + r_root2**2)
    else:
        r_tip2 = r_tip3
        r_root2 = np.sqrt(r_tip2**2 - (A2 / np.pi))

    # Calcul pour les triangles
    Vru1 = Vu1 - U1
    beta1 = np.arctan(Vru1 / Va1)
    Vr1 = np.sqrt(Va1**2 + Vru1**2)

    alpha2 = np.arctan(Vu2 / Va2)

    Vru2 = Vu2 - U2
    beta2 = np.arctan(Vru2 / Va2)
    Vr2 = np.sqrt(Va2**2 + Vru2**2)

    Vru3 = Vu3 - U3
    beta3 = np.arctan(Vru3 / Va3)
    Vr3 = np.sqrt(Va3**2 + Vru3**2)

    U_moyen = (U2 + U3) / 2.0 
    psi = (2 * dh0) / (U_moyen**2)

    # Calcul du coefficient de perte ciblé pour le rotor (hypothèse 60%)
    lambda_R = (pertes['rotor'] * perte_totale) / (0.5 * Vr3**2)

    # --- 5. Sauvegarde structurée dans les dictionnaires ---
    geom['A'] = {1: A1, 2: A2, 3: A3}
    geom['r_root'] = {1: r_root1, 2: r_root2, 3: r_root3}
    geom['r_tip'] = {1: r_tip1, 2: r_tip2, 3: r_tip3}
    geom['r_m'] = {1: r_m1, 2: r_m2, 3: r_m3}
    geom['h'] = {1: r_tip1-r_root1, 2: r_tip2-r_root2, 3: r_tip3-r_root3}

    vitesses['Rpm'] = rpm
    vitesses['omega'] = omega
    vitesses['U'] = {1: U1, 2: U2, 3: U3}
    vitesses['Va'] = {1: Va1, 2: Va2, 3: Va3}
    vitesses['Vu'] = {1: Vu1, 2: Vu2, 3: Vu3}
    vitesses['Vr'] = {1: Vr1, 2: Vr2, 3: Vr3}
    
    # Uniformisation de tous les angles en degrés
    vitesses['alpha'] = {1: contraintes['alpha_1'], 2: np.degrees(alpha2), 3: contraintes['alpha_3']}
    vitesses['beta'] = {1: np.degrees(beta1), 2: np.degrees(beta2), 3: np.degrees(beta3)}

    # --- 6. Affichage des résultats ---
    print(f"Régime : {rpm:.0f} RPM")
    print(f"Rayons moyens [m] : r_m1={r_m1:.4f} | r_m2={r_m2:.4f} | r_m3={r_m3:.4f}")
    print(f"Vitesses U [m/s]  : U1={U1:.2f} | U2={U2:.2f} | U3={U3:.2f}")
    print(f"Vitesses Va [m/s] : Va1={Va1:.2f} | Va2={Va2:.2f} | Va3={Va3:.2f}")
    print(f"Angles [deg]      : Alpha2={np.degrees(alpha2):.2f}° | Beta2={np.degrees(beta2):.2f}° | Beta3={np.degrees(beta3):.2f}°")
    print(f"Pertes            : Zeta_S={lambda_N:.4f} | Zeta_R={lambda_R:.4f}")
    print(f"Degré de réaction : {contraintes['reaction']:.3f}")
    print(f"Coefficient de chargement (Psi) : {psi:.2f}")

def plot_geometrie_turbine():
    # Extraction des données du dictionnaire
    stations = [1, 2, 3]
    # On définit des positions axiales arbitraires pour la visualisation
    x = [0, 1, 2] 
    
    r_root = [geom['r_root'][s] for s in stations]
    r_tip = [geom['r_tip'][s] for s in stations]
    r_m = [geom['r_m'][s] for s in stations]

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

def tracer_limites_rpm():

    A3 = geom['A'][3]
    rpm = vitesses['Rpm']

    # Calcul des RPM correspondants aux contraintes AN^2 max et min
    N_min = np.sqrt(contraintes['AN2_min'] / A3)
    N_max = np.sqrt(contraintes['AN2_max'] / A3)
    
    fig, ax = plt.subplots(figsize=(10, 3))
    
    # Tracer la zone acceptable (ligne verte)
    ax.plot([N_min, N_max], [0, 0], color='lightgreen', linewidth=10, solid_capstyle='round', label='Plage acceptable (Critère AN²)')
    
    # Tracer le RPM imposé
    couleur_point = 'blue' if (N_min <= rpm <= N_max) else 'red'
    ax.plot(rpm, 0, marker='o', color=couleur_point, markersize=12, label='RPM Imposé')
    
    # Annotations
    ax.text(N_min, 0.05, f'Min: {N_min:.0f} RPM', ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax.text(N_max, 0.05, f'Max: {N_max:.0f} RPM', ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax.text(rpm, -0.05, f'Actuel: {rpm:.0f}', ha='center', va='top', color=couleur_point, fontsize=10, fontweight='bold')
    
    # Mise en forme du graphique
    ax.set_yticks([]) # Cacher l'axe Y
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Ajuster l'axe X pour bien voir les limites
    marge = (N_max - N_min) * 0.2
    ax.set_xlim(min(N_min, rpm) - marge, max(N_max, rpm) + marge)
    
    ax.set_title("Vérification du Régime (RPM) vs Contraintes Structurelles (AN²)")
    ax.set_xlabel("Vitesse de rotation N [RPM]")
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    plt.show()

def tracer_triangles_vitesses():
    # 1. Extraction des variables du dictionnaire 'vitesses' (qui doit être global ici)
    U1 = vitesses['U'][1]
    U2 = vitesses['U'][2]
    U3 = vitesses['U'][3]
    
    Va1 = vitesses['Va'][1]
    Va2 = vitesses['Va'][2]
    Va3 = vitesses['Va'][3]

    Vu1 = vitesses['Vu'][1]
    Vu2 = vitesses['Vu'][2]
    Vu3 = vitesses['Vu'][3]
    
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