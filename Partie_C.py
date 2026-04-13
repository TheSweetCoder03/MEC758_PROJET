# Partie_C.py
import Partie_A
import Partie_B

def partie_c():
    donnees_hpt = Partie_A.export_donnees_hpt()

    # Récupérer la géométrie fixe
    r_root3_conception = Partie_B.donnees['r_root'][3]

    # Réduit de 20% le RPM de la HPT
    rpm_design = donnees_hpt['rpmhpc']
    donnees_hpt_hc = {**donnees_hpt, 'rpmhpc': rpm_design * 0.80}

    # Recalculer les triangles avec géométrie figée
    Partie_B.etape_1(donnees_hpt_hc, racine_constante=False, r_root3_fixe=r_root3_conception)

    # Tracer les triangles de vitesses
    Partie_B.tracer_triangles_vitesses()
