import numpy as np
from scipy.optimize import fsolve, minimize_scalar
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Tableau

contraintes = {
    # Paramètres de l'étage
    'M1': 0.14,                    # Mach entrée 
    'M3': 0.37,                    # Mach sortie 
    'alpha_1': -10.0,              # Angle absolu entrée stator 
    'alpha_3': -22.0,              # Angle absolu sortie rotor 
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

donnees = {}       # Pour stocker des variables
resultats_radiaux = {}

def deg2rad(angle):
    return angle * np.pi / 180.0

def etape_1(donnees_hpt, tolerance=1e-6):
    print(f"\nÉTAPES 1: Géométrie et Triangles")
    
    # Constantes de Sutherland (Air/Gaz)
    mu0 = 1.716e-5
    T0_suth = 273.15
    S = 110.4
    
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

    Vru3 = U3 - Vu3
    alpha_rel3 = np.arctan(Vru3 / Va3)
    Vr3 = np.sqrt(Va3**2 + Vru3**2)
    Mr3 = Vr3 / np.sqrt(gamma * R_gaz * T3)
    P0r3 = P3 * (1 + ((gamma -1)/2) * Mr3**2)**(gamma / (gamma -1))

    # Station 1 (Entrée Stator) 
    M1 = contraintes['M1']
    T1 = T01 / (1 + ((gamma - 1) / 2) * M1**2)
    P1 = P01 / ((1 + ((gamma - 1) / 2) * M1**2)**(gamma / (gamma - 1)))
    rho1 = P1 / (R_gaz * T1)

    V1 = M1 * np.sqrt(gamma * R_gaz * T1)
    Va1 = V1 * np.cos(deg2rad(contraintes['alpha_1']))
    Vu1 = V1 * np.sin(deg2rad(contraintes['alpha_1']))

    A1 = m_dot / (rho1 * Va1)
    r_tip1 = r_tip3
    r_root1 = np.sqrt(r_tip1**2 - (A1 / np.pi))  
    r_m1 = (r_root1 + r_tip1) / 2.0

    # Rendement visé
    eta_hpt = donnees_hpt['eta_iso']
    #Calcul des pertes totales de la HPT
    dh0s = dh0 / eta_hpt
    pertes_totale = dh0s - dh0

    # Thermodynamique de la station 2 (fixée par la réaction)
    dh = dh0 - (V1**2 - V3**2) / 2
    T2 = T3 + (contraintes['reaction'] * dh) / cp
    V2 = np.sqrt(2 * cp * (T01 - T2))
    T02 = T01

    M2 = V2 / np.sqrt(gamma * R_gaz * T2)

    A2 = A3
    
    def tracer_diagnostic_va2(v_min=50, v_max=350):
        print(f"Génération du graphique de diagnostic...")
        v_test = np.linspace(v_min, v_max, 100)
        y_rendement, y_pertes_N, y_pertes_R = [], [], []

        for v in v_test:
            # Calcul simplifié identique au solveur
            rho2_tmp = m_dot / (v * A2)
            P2_tmp = rho2_tmp * R_gaz * T2
            P02_tmp = P2_tmp * (1 + (gamma - 1) * M2**2 / 2)**(gamma / (gamma - 1))
            yn = (P01 - P02_tmp) / (P02_tmp - P2_tmp)
            
            vu2_tmp = np.sqrt(max(0, V2**2 - v**2))
            u2_tmp = (dh0 + U3 * Vu3) / (vu2_tmp + 1e-6)
            vru2_tmp = vu2_tmp - u2_tmp
            vr2_tmp = np.sqrt(v**2 + vru2_tmp**2)
            mr2_tmp = vr2_tmp / np.sqrt(gamma * R_gaz * T2)
            p0r2_tmp = P2_tmp * (1 + (gamma-1)/2 * mr2_tmp**2)**(gamma/(gamma-1))
            yr = (p0r2_tmp - P0r3) / (P0r3 - P3)
            
            ln = yn / (1 + 0.5 * gamma * M2**2)
            lr = yr / (1 + 0.5 * gamma * Mr3**2)
            rend = 1 / (1 + (ln * V2**2 + lr * Vr3**2) / (2 * cp * (T01 - T03)))
            
            y_rendement.append(rend)
            y_pertes_N.append(yn)
            y_pertes_R.append(yr)

        plt.figure(figsize=(10, 6))
        plt.subplot(2,1,1)
        plt.plot(v_test, y_rendement, label='Rendement (eta)')
        plt.axhline(eta_hpt, color='r', linestyle='--', label='Cible 0.88')
        plt.ylabel('Rendement'); plt.legend(); plt.grid(True)
        
        plt.subplot(2,1,2)
        plt.plot(v_test, y_pertes_N, label='Y_Stator')
        plt.plot(v_test, y_pertes_R, label='Y_Rotor')
        plt.ylabel('Pertes Y'); plt.xlabel('Va2 [m/s]'); plt.legend(); plt.grid(True)
        plt.axhline(1, color='r', linestyle='--', label='Cible 0.88')
        plt.ylim(0, 2)
        plt.show()

    tracer_diagnostic_va2(v_min=25, v_max=350)


    def objectif_Va2(Va2):
        rho2 = m_dot / (Va2 * A2)
        P2 = rho2 * R_gaz * T2
        P02 = P2 * (1 + (gamma - 1) * M2**2 / 2)**(gamma / (gamma - 1))
        Y_N = (P01 - P02) / (P02 - P2)

        Vu2 = np.sqrt(max(0, V2**2 - Va2**2)) 
        U2 = (dh0 + U3 * Vu3) / (Vu2 + 1e-6)
        Vru2 = Vu2 - U2
        Vr2 = np.sqrt(Va2**2 + Vru2**2)
        Mr2 = Vr2 / np.sqrt(gamma * R_gaz * T2)
        P0r2 = P2 * (1 + ((gamma - 1)/2) * Mr2**2)**(gamma / (gamma - 1))
        Y_R = (P0r2 - P0r3) / (P0r3 - P3)

        if Y_N < 0.001 or Y_R < 0.001:
            return 1e4 + abs(Y_N)*1000 + abs(Y_R)*1000 

        lambda_N = Y_N / (1 + 0.5 * gamma * M2**2)
        lambda_R = Y_R / (1 + 0.5 * gamma * Mr3**2)
        rendement = 1 / (1 + (lambda_N * V2**2 + lambda_R * Vr3**2) / (2 * cp * (T01 - T03)))

        return (rendement - eta_hpt)**2

    res = minimize_scalar(objectif_Va2, bounds=(100, 250), method='bounded')
    
    # Vérification de la convergence
    if not res.success:
        raise ValueError(f"[ERREUR] minimize_scalar n'a pas convergé : {res.message}")

    Va2 = res.x
    
    # --- Reprise des calculs de vérification avec la valeur validée ---
    rho2 = m_dot / (Va2 * A2)
    P2 = rho2 * R_gaz * T2
    P02 = P2 * (1 + (gamma - 1) * M2**2 /2)**(gamma / (gamma -1))
    Y_N = (P01 - P02) / (P02 - P2)

    Vu2 = np.sqrt(V2**2 - Va2**2)
    U2 = (dh0 + U3 * Vu3) / Vu2
    Vru2 = Vu2 - U2
    Vr2 = np.sqrt(Va2**2 + Vru2**2)
    Mr2 = Vr2 / np.sqrt(gamma * R_gaz * T2)
    P0r2 = P2 * (1 + ((gamma -1)/2) * Mr2**2)**(gamma / (gamma -1))
    Y_R = (P0r2 - P0r3) / (P0r3 - P3)

    lambda_N = Y_N / (1 + 0.5 * gamma * M2**2)
    lambda_R = Y_R / (1 + 0.5 * gamma * Mr3**2)
    rendement = 1 / (1 + (lambda_N * V2**2 + lambda_R * Vr3**2) / (2 * cp * (T01 - T03)))

    # Calcul dimensions ailette à station 2
    r_m2 = U2 / omega
    h2 = A2 / (2 * np.pi * r_m2)
    r_root2 = r_m2 - (h2 / 2)
    r_tip2 = r_m2 + (h2 / 2)

    # CALCUL DES ANGLES (ABSOLUS ET RELATIFS)

    # Station 2
    alpha2 = np.arctan(Vu2 / Va2)
    alpha_rel2 = np.arctan(Vru2 / Va2)

    coeff_perte_temp = (2 / U2**2) * (U2 * (U2 - Va2 * np.tan(alpha_rel2)) - U3 * (U3 - Va3 * np.tan(alpha_rel3)))
    coeff_ecoulement = (Va2 + Va3) / (2 * U3)

    # Viscosité stator/rotor
    Visc_s = mu0 * (T2 / T0_suth)**1.5 * (T0_suth + S) / (T2 + S)
    Visc_r = mu0 * (T3 / T0_suth)**1.5 * (T0_suth + S) / (T3 + S)

    # --- Sauvegarde structurée ---
    donnees['P'] = {1: P1, 2: P2, 3: P3}
    donnees['M'] = {1: M1, 2: M2, 3: M3}
    donnees['Mr'] = {2: Mr2, 3: Mr3}
    donnees['To'] = {1: T01, 2: T02, 3: T03}
    donnees['T'] = {1: T1, 2: T2, 3: T3}
    donnees['p'] = {1: rho1, 2: rho2, 3: rho3}
    donnees['alpha'] = {1: deg2rad(contraintes['alpha_1']), 2: alpha2, 3: deg2rad(contraintes['alpha_3'])}
    donnees['alpha_relatif'] = { 2: alpha_rel2, 3: alpha_rel3}

    donnees['A'] = {1: A1, 2: A2, 3: A3}
    donnees['r_root'] = {1: r_root1, 2: r_root2, 3: r_root3}
    donnees['r_tip'] = {1: r_tip1, 2: r_tip2, 3: r_tip3}
    donnees['r_m'] = {1: r_m1, 2: r_m2, 3: r_m3}
    donnees['h'] = {1: r_tip1-r_root1, 2: r_tip2-r_root2, 3: r_tip3-r_root3}

    donnees['Rpm'] = rpm
    donnees['omega'] = omega
    donnees['U'] = {1: U2, 2: U3}
    donnees['Va'] = {1: Va1, 2: Va2, 3: Va3}
    donnees['Vu'] = {1: Vu2, 2: Vu3, 3: Vu1}
    donnees['Vr'] = {2: Vr2, 3: Vr3}

    donnees['Y'] = {1: Y_N, 2: Y_R}
    donnees['nst'] = rendement
    donnees['Visc'] = {'stator': Visc_s, 'rotor': Visc_r}
    donnees['R_gaz'] = R_gaz

    print(f"Régime : {rpm:.0f} RPM")
    print(f"Rayons à l'emplanture: Rr1={r_root1:.3f} | Rr2={r_root2:.3f} | Rr3={r_root3:.3f}")
    print(f"Rayons à l'extréminté: Rt1={r_tip1:.3f} | Rt2={r_tip2:.3f} | Rt3={r_tip3:.3f}")
    print(f"Rayons moyen: Rm1={r_m1:.3f} | Rm2={r_m2:.3f} | Rm3={r_m3:.3f}")
    print(f"Vitesses Va [m/s] : Va1={Va1:.2f} | Va2={Va2:.2f} | Va3={Va3:.2f}")
    print(f"Vitesses V [m/s] : V1={V1:.2f} | V2={V2:.2f} | V3={V3:.2f}")
    print(f"Angles Abs.[deg]  : Alpha1={contraintes['alpha_1']:.2f}° | Alpha2={np.degrees(alpha2):.2f}° | Alpha3={contraintes['alpha_3']:.2f}°")
    print(f"Angles Rel.[deg]  :  Alpha_rel2={np.degrees(alpha_rel2):.2f}° | Alpha_rel3={np.degrees(alpha_rel3):.2f}°")
    print(f"Coefficient de pertes : Stator (Y_N) = {Y_N:.5f}, Rotor (Y_R) = {Y_R:.5f}")
    print(f"Coefficient de charge : {coeff_perte_temp:.3f}")
    print(f"Coefficient d'écoulement : {coeff_ecoulement:.3f}")
    print(f"Rendement de l'étage: {rendement:.4f}")
    print(f"Mr3 : {Mr3:.3f}")

def plot_geometrie_turbine():
    # Extraction des données du dictionnaire
    stations = [1, 2, 3]
    # On définit des positions axiales arbitraires pour la visualisation
    x = [0, 1, 2] 
    
    r_root = [donnees['r_root'][s] for s in stations]
    r_tip = [donnees['r_tip'][s] for s in stations]
    r_m = [donnees['r_m'][s] for s in stations]

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

    A3 = donnees['A'][3]
    rpm = donnees['Rpm']

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
    a1 = donnees['alpha'][1]
    
    U2 = donnees['U'][1]
    U3 = donnees['U'][2]
    
    Va1 = donnees['Va'][1]
    Va2 = donnees['Va'][2]
    Va3 = donnees['Va'][3]

    Vu2 = donnees['Vu'][1]
    Vu3 = donnees['Vu'][2]
    
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

def etape_2(donnees_hpt):
    print(f"\nÉTAPES 2.a, b, c : Répartition des triangles de vitesse")

    # Valeurs lues depuis les dictionnaires globaux 'donnees' et 'contraintes'
    Vu_1 = donnees['Vu'][3]
    Vu_2 = donnees['Vu'][1]
    Vu_3 = donnees['Vu'][2]
    rm_1 = donnees['r_m'][1]
    rm_2 = donnees['r_m'][2]
    rm_3 = donnees['r_m'][3]
    Pm_1 = donnees['P'][1]
    Pm_2 = donnees['P'][2]
    Pm_3 = donnees['P'][3]
    pm_1 = donnees['p'][1]
    pm_2 = donnees['p'][2]
    pm_3 = donnees['p'][3]
    alpham_1 = donnees['alpha'][1]
    alpham_2 = donnees['alpha'][2]
    alpham_3 = donnees['alpha'][3]
    rr_1 = donnees['r_root'][1]
    rr_2 = donnees['r_root'][2]
    rr_3 = donnees['r_root'][3]
    rt_1 = donnees['r_tip'][1]
    rt_2 = donnees['r_tip'][2]
    rt_3 = donnees['r_tip'][3]
    f_t_m = contraintes['reaction']
    Va_1 = donnees['Va'][1]
    Va_2 = donnees['Va'][2]
    Va_3 = donnees['Va'][3]
    To_1 = donnees['To'][1]
    To_2 = donnees['To'][2]
    To_3 = donnees['To'][3]
    omega = donnees['omega']
    
    cp = donnees_hpt['cp']
    gamma = donnees_hpt['gamma']
    R_gaz = cp * (gamma - 1) / gamma

    # Initialisation
    donnees['M_hub'] = {}

    # Hypothèse de Free Vortex
    n = -1

    # Trouver les constantes à chaque station
    stat_1 = [Vu_1, rm_1, Pm_1, alpham_1, pm_1, rr_1, rt_1, Va_1, To_1]
    stat_2 = [Vu_2, rm_2, Pm_2, alpham_2, pm_2, rr_2, rt_2, Va_2, To_2]
    stat_3 = [Vu_3, rm_3, Pm_3, alpham_3, pm_3, rr_3, rt_3, Va_3, To_3]
    
    for i in [stat_1, stat_2, stat_3]:
        k1 = i[0] * i[1]
        k2 = (i[1]**2) * (1 - f_t_m)
        k3 = np.tan(i[3]) * i[1]
        i.extend([k1, k2, k3])

    # Calculer variables le long du rayon
    for idx, j in enumerate([stat_1, stat_2, stat_3], start=1):
        r = np.linspace(j[5], j[6], 100)
        f_t = 1 - j[10] / r**2
        if idx == 1:
            Vu = np.full_like(r, j[0])
            P = np.full_like(r, j[2])
        else:
            Vu = j[9] / r
            P = (((j[9]**2)*j[4])/2) * (1/j[1]**2 - 1/r**2) + j[2]
        V = np.sqrt(j[7]**2 + Vu**2)
        T = j[8] - V**2 / (2*cp)
        if idx == 2 or idx == 3:
            u = r * omega
            alpha = np.degrees(np.arctan((Vu - u) / j[7]))
        else:
            alpha = np.full_like(r, np.degrees(j[3]))
        a = np.sqrt(gamma * R_gaz * T)
        M = V / a
        
        donnees['M_hub'][idx] = M[0]

        # Sauvegarde pour les graphiques
        resultats_radiaux[idx] = {
            'r': r,
            'f_t': f_t,
            'Vu': Vu,
            'P': P,
            'alpha': alpha,
            'T': T
        }

def etape_2_graphique():
    for idx, res in resultats_radiaux.items():
        r = res['r']
        f_t = res['f_t']
        Vu = res['Vu']
        P = res['P']
        alpha = res['alpha']
        T = res['T']

        fig, axes = plt.subplots(2, 3, figsize=(10, 8))
        fig.suptitle(f"Station {idx}", fontsize=14, fontweight='bold')

        axes[0, 0].plot(f_t, r)
        axes[0, 0].set_ylabel("r (m)")
        axes[0, 0].set_xlabel("Λ")
        axes[0, 0].set_title("Λ vs r")
        axes[0, 0].tick_params(axis='x', rotation=90)

        axes[0, 1].plot(Vu, r)
        axes[0, 1].set_ylabel("r (m)")
        axes[0, 1].set_xlabel("Vu")
        axes[0, 1].set_title("Vu vs r")
        axes[0, 1].tick_params(axis='x', rotation=90)
        
        # Désactiver le graphique vide
        axes[0, 2].axis('off')

        axes[1, 0].plot(P, r)
        axes[1, 0].set_ylabel("r (m)")
        axes[1, 0].set_xlabel("P")
        axes[1, 0].set_title("P vs r")
        axes[1, 0].tick_params(axis='x', rotation=90)

        axes[1, 1].plot(alpha, r)
        axes[1, 1].set_ylabel("r (m)")
        axes[1, 1].set_xlabel("Vrillage")
        axes[1, 1].set_title("Vrillage vs r")
        axes[1, 1].tick_params(axis='x', rotation=90)

        axes[1, 2].plot(T, r)
        axes[1, 2].set_ylabel("r (m)")
        axes[1, 2].set_xlabel("T")
        axes[1, 2].set_title("Température vs r")
        axes[1, 2].tick_params(axis='x', rotation=90)

        plt.tight_layout()
        plt.show()

def etape_3(donnees_hpt):
   
    fs = contraintes['stator_h_c']
    fr = contraintes['rotor_h_c']
    zweifels = contraintes['stator_zweifel']
    zweifelr = contraintes['rotor_zweifel']
    h1 = donnees['r_tip'][1] - donnees['r_root'][1]
    h2 = donnees['r_tip'][2] - donnees['r_root'][2]
    h3 = donnees['r_tip'][3] - donnees['r_root'][3] 
    rm2 = donnees['r_m'][2]
    rm3 = donnees['r_m'][3]
    alpha_s1 = abs(donnees['alpha'][1])
    alpha_s2 = abs(donnees['alpha'][2])
    alpha_r2 = abs(donnees['alpha_relatif'][2])
    alpha_r3 = abs(donnees['alpha_relatif'][3])
    gamma_s = (alpha_s1 + alpha_s2) / 2
    gamma_r = (alpha_r2 + alpha_r3) / 2

    # Calcul hauteur moyenne ailette
    hms = (h1 + h2) / 2
    hmr = (h2 + h3) / 2

    # Calcul des cordes (c = h / Facteur de forme)
    cs = hms / fs
    cr = hmr / fr
    
    # Calcul de la corde axiale (c_a = c * cos(gamma))
    cas = cs * np.cos(gamma_s)
    car = cr * np.cos(gamma_r)

    # Calcul du pas
    pas_s = (zweifels * cas) / (2 * (np.tan(alpha_s1) + np.tan(alpha_s2)) * (np.cos(alpha_s2))**2)
    pas_r = (zweifelr * car) / (2 * (np.tan(alpha_r2) + np.tan(alpha_r3)) * (np.cos(alpha_r3))**2)

    # Calcul du nombre d'ailettes (Périmètre moyen / pas)
    ns = int(np.ceil(2 * np.pi * rm2 / pas_s))
    nr = int(np.ceil(2 * np.pi * rm3 / pas_r))

    

    # --- Sauvegarde structurée ---
    donnees['hm'] = {'stator': hms, 'rotor': hmr}
    donnees['c'] = {'stator': cs, 'rotor': cr}
    donnees['ca'] = {'stator': cas, 'rotor': car}
    donnees['pas'] = {'stator': pas_s, 'rotor': pas_r}
    
    print(f"\nÉTAPES 3: Paramètres des aubes")
    print(f"Stator | Corde axiale: {cas:.4f} m | Pas: {pas_s:.4f} m | Aubes: {ns:.1f}")
    print(f"Rotor  | Corde axiale: {car:.4f} m | Pas: {pas_r:.4f} m | Aubes: {nr:.1f}")

def etape_4(donnees_hpt):
    # Extraction des variables du dictionnaire 'donnees'
    beta_1 = donnees['alpha'][1] 
    beta_2 = donnees['alpha_relatif'][2]
    beta_3 = donnees['alpha_relatif'][3]

    alpha_1 = donnees['alpha'][1]
    alpha_2 = donnees['alpha'][2]
    alpha_3 = donnees['alpha'][3]

    P1 = donnees['P'][1]
    P2 = donnees['P'][2]
    P3 = donnees['P'][3]

    M1 = donnees['M'][1]
    M2 = donnees['M'][2]
    M3 = donnees['M'][3]

    Mr2 = donnees['Mr'][2]
    Mr3 = donnees['Mr'][3]

    V2 = np.sqrt((donnees['Va'][2])**2 + (donnees['Vu'][2])**2)
    Vr3 = donnees['Vr'][3]

    rho2 = donnees['p'][2]
    rho3 = donnees['p'][3]

    M1_hub = donnees['M_hub'][1]
    M2_hub = donnees['M_hub'][2]
    M3_hub = donnees['M_hub'][3]

    h_s = donnees['hm']['stator']
    h_r = donnees['hm']['rotor']

    c_s = donnees['c']['stator']
    c_r = donnees['c']['rotor']

    ca_s = donnees['ca']['stator']
    ca_r = donnees['ca']['rotor']

    pas_s = donnees['pas']['stator']
    pas_r = donnees['pas']['rotor']

    rrm_s = (donnees['r_root'][1] + donnees['r_root'][2]) / 2 # Rayon à la racine moyen au stator
    rtm_s = (donnees['r_tip'][1] + donnees['r_tip'][2]) / 2   # Rayon à la pointe moyen au stator
    rrm_r = (donnees['r_root'][2] + donnees['r_root'][3]) / 2 # Rayon à la racine moyen au rotor
    rtm_r = (donnees['r_tip'][2] + donnees['r_tip'][3]) / 2   # Rayon à la pointe moyen au rotor

    Visc_s = donnees['Visc']['stator']
    Visc_r = donnees['Visc']['rotor']

    # Nombre de Reynolds stator/rotor
    Re_s = (rho2 * V2 * c_s) / Visc_s
    Re_r = (rho3 * Vr3 * c_r) / Visc_r

    # Paramètres de conception
    nbre_seal = 3 # Nombre de seal au bout ailette
    tmax_c_s = 0.2    # Ratio Épaisseur max ailette sur corde  stator
    tmax_c_r = 0.2    # Ratio Épaisseur max ailette sur corde rotor
    k_s = 0                # Jeu radial ailette stator
    k_r = 0.0005188         # Jeu radial ailette rotor, à ajuster une fois que le jeu nécessaire est trouvé, pour atteindre le même coefficient de perte au rotor

    # Extraction des variables du dictionnaire 'donnees_hpt'
    y = donnees_hpt['gamma']

    # ====================================
    # PERTE DANS LE STATOR (Station 1 à 2)
    # ====================================
    
    # Extractions graphiques pour le stator
    Yp_1_s = Tableau.extraire_donnee_graphique(x=(pas_s/c_s), y=abs(np.degrees(alpha_2)), figure=1)
    Yp_2_s = Tableau.extraire_donnee_graphique(x=(pas_s/c_s), y=abs(np.degrees(alpha_2)), figure=2)

    # Perte du profil (Yp) -
    Yp_AMDC_s = (Yp_1_s + abs(beta_1/alpha_2)*(beta_1/alpha_2) * (Yp_2_s - Yp_1_s)) * ((tmax_c_s) / 0.2)**(beta_1/alpha_2)

    k1_s = 1 - 1.25 * abs(M2 - 0.2) if M2 > 0.2 else 1
    k2_s = abs(M1 / M2)**2
    kp_s = 1 - k2_s * (1 - k1_s)

    dP_hub_s = 0.75 * (M1_hub - 0.4)**1.75 if M1_hub > 0.4 else 0

    dp_shock_s = (rrm_s / rtm_s) * dP_hub_s
    Yshock_s = dp_shock_s * (P1 / P2) * ((1 - (1 + (y - 1) / 2 * M1**2 )**(y / (y - 1))) / (1 - (1 + (y - 1) / 2 * M2**2 )**(y / (y - 1))))
    Yp_moderne_s = 0.914 * ((2 / 3) * Yp_AMDC_s * kp_s + Yshock_s)

    # Perte due aux écoulements secondaires (Ys)
    f_ar_s = (1 - 0.25 * np.sqrt(2 - (h_s / c_s))) / (h_s / c_s) if (h_s / c_s) <= 2 else 1 / (h_s / c_s)

    alpha_m_s = np.arctan((1 / 2) * (np.tan(alpha_1) - np.tan(alpha_2)))
    Cl_sc_s = 2 * (np.tan(abs(alpha_1)) + np.tan(abs(alpha_2))) * np.cos(alpha_m_s)

    Ys_AMDC_s = 0.0334 * f_ar_s * (np.cos(alpha_2) / np.cos(beta_1)) * Cl_sc_s**2 * (np.cos(alpha_2))**2 / (np.cos(alpha_m_s))**3

    k3_s = (1 / (h_s / ca_s))**2
    ks_s = 1 - (k3_s * (1 - kp_s))    
    
    Ys_moderne_s = 1.2 * Ys_AMDC_s * ks_s

    # 2. Extraction des données spécifiques au stator
    ratio_te_s = contraintes['stator_te_thick'] / pas_s
    d_tet_0_s = Tableau.extraire_donnee_graphique(x=ratio_te_s, y=None, figure=140)
    d_tet_alpha_s = Tableau.extraire_donnee_graphique(x=ratio_te_s, y=None, figure=141)

    # Perte annulaire (Ytet)
    d_tet_s = d_tet_0_s + abs(beta_1 / alpha_2) * (beta_1 / alpha_2) * (d_tet_alpha_s - d_tet_0_s) 
    denom_s = 1 - (1 + (y - 1) / 2 * M2**2)**(-y / (y - 1))
    Ytet_s = ((1 - (y - 1) / 2 * M2**2 * (1 / (1 - d_tet_s) - 1))**(-y / (y - 1)) - 1) / denom_s

    # Perte due au jeu radial (Ytc)
    k_2_s = k_s / (nbre_seal)**0.42
    Ytc_s = 0.37 * (c_s / h_s) * (k_2_s / c_s)**0.78 * (Cl_sc_s)**2 * ((np.cos(alpha_2)**2) / (np.cos(alpha_m_s))**3)

    # Perte totale Stator (Ytot_s)
    if Re_s <= 200000:
        f_re_s = (Re_s / 200000)**-0.4
    elif 200000 < Re_s <= 1000000:
        f_re_s = 1
    else:
        f_re_s = (Re_s / 1000000)**-0.2

    Ytot_s = Yp_moderne_s * f_re_s + Ys_moderne_s + Ytet_s + Ytc_s

    print(f"Coeff de pertes stator: Yp={Yp_moderne_s:.3f} | Fre={f_re_s:.3f} | Ys={Ys_moderne_s:.3f} | Ytet={Ytet_s:.3f} | Ytc={Ytc_s:.3f} | Re={Re_s:.3f}")


    # ===================================
    # PERTE DANS LE ROTOR (Station 2 à 3)
    # ===================================
    
    # Extractions graphiques pour le rotor
    Yp_1_r = Tableau.extraire_donnee_graphique(x=(pas_r/c_r), y=abs(np.degrees(beta_3)), figure=1)
    Yp_2_r = Tableau.extraire_donnee_graphique(x=(pas_r/c_r), y=abs(np.degrees(beta_3)), figure=2)

    # Perte du profil (Yp)
    Yp_AMDC_r = (Yp_1_r + abs(beta_2 / beta_3) * (beta_2 / beta_3) * (Yp_2_r - Yp_1_r)) * ((tmax_c_r) / 0.2)**(beta_2 / beta_3)

    k1_r = 1 - 1.25 * abs(Mr3 - 0.2) if Mr3 > 0.2 else 1
    k2_r = abs(Mr2 / Mr3)**2
    kp_r = 1 - k2_r * (1 - k1_r)

    dP_hub_r = 0
    if M2_hub > 0.4:
        dP_hub_r = 0.75 * (M2_hub - 0.4)**1.75

    dp_shock_r = (rrm_r / rtm_r) * dP_hub_r

    Yshock_r = dp_shock_r * (P2 / P3) * ((1 - (1 + (y - 1) / 2 * Mr2**2 )**(y / (y - 1))) / (1 - (1 + (y - 1) / 2 * Mr3**2 )**(y / (y - 1))))

    Yp_moderne_r = 0.914 * ((2 / 3) * Yp_AMDC_r * kp_r + Yshock_r)


    # Perte secondaire (Ys) 
    f_ar_r = (1 - 0.25 * np.sqrt(2 - (h_r / c_r))) / (h_r / c_r) if (h_r / c_r) <= 2 else 1 / (h_r / c_r)

    beta_m_r = np.arctan(0.5 * (np.tan(beta_2) + np.tan(beta_3)))
    Cl_sc_r = 2 * (np.tan(abs(beta_2)) + np.tan(abs(beta_3))) * np.cos(beta_m_r)

    Ys_AMDC_r = 0.0334 * f_ar_r * (np.cos(beta_3) / np.cos(beta_2)) * Cl_sc_r**2 * (np.cos(beta_3))**2 / (np.cos(beta_m_r))**3

    k3_r = (1 / (h_r / ca_r))**2
    ks_r = 1 - (k3_r * (1 - kp_r))    
    
    Ys_moderne_r = 1.2 * Ys_AMDC_r * ks_r

    # 1. Calcul du ratio d'épaisseur pour le rotor
    ratio_te_r = contraintes['rotor_te_thick'] / pas_r

    d_tet_0_r = Tableau.extraire_donnee_graphique(x=ratio_te_r, y=None, figure=140)
    d_tet_alpha_r = Tableau.extraire_donnee_graphique(x=ratio_te_r, y=None, figure=141)

    # Perte annulaire (Ytet)
    d_tet_r = d_tet_0_r + abs(beta_2 / beta_3) * (beta_2 / beta_3) * (d_tet_alpha_r - d_tet_0_r)
    denom_r = 1 - (1 + (y - 1) / 2 * Mr3**2)**(-y / (y - 1))
    Ytet_r = ((1 - (y - 1) / 2 * Mr3**2 * (1 / (1 - d_tet_r) - 1))**(-y / (y - 1)) - 1) / denom_r

    # Perte due au jeu radial (Ytc)
    k_2_r = k_r / (nbre_seal)**0.42
    Ytc_r = 0.37 * (c_r / h_r) * (k_2_r / c_r)**0.78 * (Cl_sc_r)**2 * ((np.cos(beta_3)**2) / (np.cos(beta_m_r))**3)

    # Perte totale Rotor (Ytot_r)
    if Re_r <= 200000:
        f_re_r = (Re_r / 200000)**-0.4
    elif 200000 < Re_r <= 1000000:
        f_re_r = 1
    else:
        f_re_r = (Re_r / 1000000)**-0.2

    Ytot_r = Yp_moderne_r * f_re_r + Ys_moderne_r + Ytet_r + Ytc_r

    print(f"Coeff de pertes rotor: Yp={Yp_moderne_r:.3f} | Fre={f_re_r:.3f} | Ys={Ys_moderne_r:.3f} | Ytet={Ytet_r:.3f} | Ytc={Ytc_r:.3f} | Re={Re_r:.3f}")  

    # ==========================
    # 4.b : Trouver le jeu rotor 
    # ==========================
    
    Ytot_R2 = donnees['Y'][2]
    Ytc_requis = Ytot_R2 - (Yp_moderne_r * f_re_r) - Ys_moderne_r - Ytet_r
    
    print(f"\nÉTAPES 4: Coefficient de perte")
    if Ytc_requis <= 0:
        k_requis = 0
        Jeu_requis = 0
        print(f"ATTENTION : Budget de perte insuffisant ! (Ytc_requis = {Ytc_requis:.4f})")
        print(f"Les pertes de profil et secondaires sont déjà supérieures au total visé.")
    
    else:
        denominateur = (0.37 * (c_r / h_r) * (Cl_sc_r**2) * (np.cos(beta_3)**2 / np.cos(beta_m_r)**3))
        k_requis = c_r * (Ytc_requis / denominateur)**(1 / 0.78)
        Jeu_requis = k_requis * (nbre_seal)**(0.42)
        
    donnees['coefficient_perte_rotor'] = {1: Yp_moderne_r, 2: Ys_moderne_r, 3: Ytet_r, 4: Ytc_r}
    donnees['Y_perte_stator'] = Ytot_s

    print(f"Coefficient de pertes : Stator (Y_N) = {Ytot_s:.5f}, Rotor (Y_R) = {Ytot_r:.5f}")
    print(f"Jeu radial rotor requis : {Jeu_requis * 1000:.4f} mm")

    # ===================================================
    # 4.c : Proportion aube AH/AT pour durée de vie rotor
    # ===================================================

    rho_m = 0.315
    k1 = 15
    k2 = 55.6
    k3 = -5.2
    k4 = 0.6

    # Récupération des variables importantes
    cp = donnees_hpt['cp']
    t = contraintes['vie_heures']
    vr2 = donnees['Vr'][2]

    Tm = (donnees['T'][2] + (vr2**2 / (2 * cp)))*1.8
    
    an_2 = (donnees['A'][2] * 1550) * donnees['Rpm']**2
    
    LM_M = ((k1 + np.log10(t)) / 10**3) * (Tm + 175)
    
    if LM_M <= 44.34:
        Ar_At = 4
        k5_values = 12.6
        sigma_c = rho_m * k5_values * (an_2*10**-10)**1.08
        print("Température de métal inférieur à 1312 °C")

    else:
        # 1. Calcul de la contrainte admissible sigma_c (KSI) 
        sigma_c = fsolve(lambda s: k2 + k3 * np.log(s) + k4 * (np.log(s))**2 - LM_M, x0=0.00001)[0] 

        # 2. Calcul du coefficient géométrique K5 nécessaire 
        k5_req = sigma_c / (rho_m * (an_2 * 10**-10)**1.08)

        # 3. Identification du ratio Ar/At 
        k5_values = np.array([12.60, 11.86, 11.03, 10.68, 9.93, 9.62, 9.40])
        ratio_values = np.array([4, 5, 6, 7, 8, 9, 10])
        Ar_At = np.interp(k5_req, k5_values[::-1], ratio_values[::-1])
    
    print(f"Contrainte admissible : {sigma_c:.2f} KSI")
    print(f"Température du métal de l'aube : {Tm:.2f}")
    print(f"Durée de vie maximale LM_m : {LM_M:.2f}")
    print(f"Ratio Ar/At nécessaire : {Ar_At:.2f}")


