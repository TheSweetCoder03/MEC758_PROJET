import Partie_A
import Partie_B

def main():
    print("EXÉCUTION DE LA PARTIE A : CYCLE THERMODYNAMIQUE")
    
    Partie_A.calcul()
    Partie_A.print_station()
    Partie_A.plot_cycle()

    print("")
    print("EXÉCUTION DE LA PARTIE B : CONCEPTION HPT")
    # 1. Recalcul du débit massique de la HPT (car non stocké dans le dictionnaire de base)
    mp = Partie_A.cte_comp["mpt"]
    perte = Partie_A.cte_comp["pap"]
    f = Partie_A.cte_cc["f"]
    mpt_calc = mp * (1 - perte) * (1 + f)

    # 2. Préparation du dictionnaire de liaison
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
    }

    # 3. Exécution de la conception aérodynamique
    Partie_B.calcul()

if __name__ == "__main__":
    main()