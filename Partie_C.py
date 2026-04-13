# Partie_C.py
import Partie_A
import Partie_B

def partie_c():
    donnees_hpt = Partie_A.export_donnees_hpt()
    
    # Réduction de 20% du RPM
    rpm_design = donnees_hpt['rpmhpc']
    rpm_hors_conception = rpm_design * 0.80
    donnees_hpt_hc = {**donnees_hpt, 'rpmhpc': rpm_hors_conception}
    
    # Recalculer les triangles de vitesses avec le nouveau RPM
    Partie_B.etape_1(donnees_hpt_hc, racine_constante=False)

    #Tracer les triangles de vitesses
    Partie_B.tracer_triangles_vitesses()
