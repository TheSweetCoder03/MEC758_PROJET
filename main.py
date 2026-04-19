import Partie_A
import Partie_B
import Partie_C

def main():
    print("EXÉCUTION DE LA PARTIE A : CYCLE THERMODYNAMIQUE")
    Partie_A.calcul()
    Partie_A.print_station()
    Partie_A.plot_cycle_enhanced()

    print("\nEXÉCUTION DE LA PARTIE B : CONCEPTION HPT")
    donnees_hpt = Partie_A.export_donnees_hpt()
    Partie_B.etape_1(donnees_hpt)
    #Partie_B.plot_geometrie_turbine()
    #Partie_B.tracer_limites_rpm()
    #Partie_B.tracer_triangles_vitesses()
    Partie_B.etape_2(donnees_hpt)
    Partie_B.etape_2_graphique()
    Partie_B.etape_3(donnees_hpt)
    Partie_B.etape_4(donnees_hpt)

    print("\nEXÉCUTION DE LA PARTIE C : OFF-DESIGN")
    Partie_C.partie_c()
    Partie_C.tracer_triangles_vitesses()
    Partie_C.pertes_incidence()
    Partie_C.rendement_incidence(donnees_hpt)

if __name__ == "__main__":
    main()
