import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def etape_2():
    
    print(f"\nÉTAPES 2.a, b, c : Répartition des triangles de vitesse")

    #Valeurs calculées précédemment
    Vu_1 = 1
    Vu_2 = 1
    Vu_3 = 1
    rm_1 = 1
    rm_2 = 1
    rm_3 = 1
    Pm_1 = 1
    Pm_2 = 1
    Pm_3 = 1
    pm_1 = 1
    pm_2 = 1
    pm_3 = 1
    alpham_1 = 1
    alpham_2 = 1
    alpham_3 = 1
    rr_1 = 0.5
    rr_2 = 0.5
    rr_3 = 0.5
    rt_1 = 2
    rt_2 = 2
    rt_3 = 2
    f_t_m = 0.63 # facteur de travail moyen
    Va_1 = 1
    Va_2 = 1
    Va_3 = 1
    To_1 = 1
    To_2 = 1
    To_3 = 1
    cp = 1


    #Hypothèse de Free Vortex
    n = -1

    #Trouver les constantes à chaque station
    stat_1 = [Vu_1, rm_1, Pm_1, alpham_1, pm_1, rr_1, rt_1, Va_1, To_1]
    stat_2 = [Vu_2, rm_2, Pm_2, alpham_2, pm_2, rr_2, rt_2, Va_2, To_2]
    stat_3 = [Vu_3, rm_3, Pm_3, alpham_3, pm_3, rr_3, rt_3, Va_3, To_3]
    
    for i in [stat_1, stat_2, stat_3] :
        k1 = i[0] * i[1]
        k2 = (i[1]**2) * (1 - f_t_m)
        k3 = np.tan(i[3]) * i[1]
        i.extend([k1,k2,k3])

    #Calculer variables le long du rayon
    for j in [stat_1, stat_2, stat_3]:
        r = np.linspace(j[5], j[6], 100)
        f_t = 1 - j[10] / r**2
        Vu = j[9] / r
        V = np.sqrt(j[10]**2 + Vu**2)
        T = j[8] - V**2 / (2*cp)
        P = (((j[9]**2)*j[4])/2) * (1/j[1]**2 - 1/r**2) + j[2]
        alpha = np.arctan(j[11] / r)

        #Faire graphiques

        fig, axes = plt.subplots(2,3, figsize=(10,8))

        axes[0, 0].plot(r, f_t)
        axes[0, 0].set_xlabel("r (m)")
        axes[0, 0].set_ylabel("Λ")
        axes[0, 0].set_title("Λ  vs r")

        axes[0, 1].plot(r, Vu)
        axes[0, 1].set_xlabel("r (m)")
        axes[0, 1].set_ylabel("Vu")
        axes[0, 1].set_title("Vu vs r")

        axes[1, 0].plot(r, P)
        axes[1, 0].set_xlabel("r (m)")
        axes[1, 0].set_ylabel("P")
        axes[1, 0].set_title("P vs r")

        axes[1, 1].plot(r, alpha)
        axes[1, 1].set_xlabel("r (m)")
        axes[1, 1].set_ylabel("Vrillage")
        axes[1, 1].set_title("Vrillage vs r")

        axes[1, 2].plot(r, T)
        axes[1, 2].set_xlabel("r (m)")
        axes[1, 2].set_ylabel("T")
        axes[1, 2].set_title("Température vs r")

        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    etape_2()
    