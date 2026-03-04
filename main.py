# Importation de vos deux modules
import Partie_A# Assurez-vous que votre fichier de la partie A s'appelle bien "partie_a.py"
from Partie_B import ConceptionTurbineHPT

def main():
    print("EXÉCUTION DE LA PARTIE A : CYCLE THERMODYNAMIQUE")
    
    # 1. Exécution du cycle thermodynamique
    # Note : N'oubliez pas de fermer la fenêtre du graphique T-s pour continuer !
    Partie_A.main()

    print("LIAISON ET EXÉCUTION DE LA PARTIE B : CONCEPTION HPT")
    
    # 2. Recalcul du débit massique de la HPT (car non stocké dans le dictionnaire de base)
    mp = Partie_A.cte_comp["mpt"]
    perte = Partie_A.cte_comp["pap"]
    f = Partie_A.cte_cc["f"]
    mpt_calc = mp * (1 - perte) * (1 + f)

    # 3. Préparation du dictionnaire de liaison
    donnees_hpt = {
        # Données de la thermodynamique (Partie A)
        'T01': Partie_A.station[4]['To'],                  
        'P01': Partie_A.station[4]['Po'],                  
        'T03': Partie_A.station[5]['To'],                  
        'P03': Partie_A.station[5]['Po'],                  
        'm_dot': mpt_calc,                
        'W_hpt': Partie_A.station[5]['w'],               
        'cp': Partie_A.cte_turb["cpt"],       
        'gamma': Partie_A.cte_turb["yt"],     
        'eta_iso': Partie_A.cte_turb["nthp"], 
        
        # Données imposées par le cahier des charges
        'N_rpm': 15000    # Hypothèse de départ pour respecter le critère AN^2
    }

    print("--- Résumé des conditions aux limites transférées à la HPT ---")
    print(f"T01 (Entrée HPT) : {donnees_hpt['T01']:.2f} K")
    print(f"P01 (Entrée HPT) : {donnees_hpt['P01']:.2f} Pa")
    print(f"Débit massique   : {donnees_hpt['m_dot']:.3f} kg/s\n")

    # 4. Exécution de la conception aérodynamique
    turbine = ConceptionTurbineHPT(donnees_hpt)
    
    # Lancement de la première étape avec une hypothèse de vitesse axiale
    turbine.etape_1_corde_moyenne(V_a2_guess=150.0)

if __name__ == "__main__":
    main()