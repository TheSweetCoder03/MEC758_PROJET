# MEC758 - Projet de Conception : Turboshaft à Trois Arbres

Ce projet, réalisé dans le cadre du cours MEC758 - Systèmes de propulsion et turbomachines à l'École de technologie supérieure, porte sur la conception aérodynamique d'un moteur turboshaft pour un hélicoptère en vol stationnaire.



## Conditions de Conception
Le moteur est dimensionné pour les conditions atmosphériques suivantes à une altitude de 4 000 pieds (1 219 m):
* **Pression absolue** : 12,7 PSI (87,56 kPa).
* **Température** : 100 °F (37,8 °C).
* **Vitesse d'entrée/sortie** : Négligée (vol stationnaire).

---

## Spécifications des Composants

### Compresseurs
| Paramètre | Compresseur Basse Pression (LPC) | Compresseur Haute Pression (HPC) |
| :--- | :--- | :--- |
| **Rapport de pression ($r_p$)** | 4,25 | 2,65 |
| **Rendement ($\eta$)** | 0,86 | 0,84 |
| **Propriétés de l'air** | $\gamma=1,4$ ; $C_p=1005$ J/kg·K | $\gamma=1,4$ ; $C_p=1005$ J/kg·K |

* **Débit massique ($\dot{m}$)** : 5,443 kg/s.
* **Refroidissement** : 10% du débit d'entrée est prélevé pour le refroidissement.

### Combustion et Turbines
* **Chambre de combustion** : $f=0,02$ , $\eta=0,99$ , $\Delta P/P=0,02$.
* **Propriétés (Gaz Chauds)** : $\gamma=1,31$ et $C_p=1214$ J/kg·K.
* **Rendements Isentropiques** :
    * Turbine Haute Pression (HPT) : 0,88.
    * Turbine Basse Pression (LPT) : 0,90.
    * Turbine de Puissance (PT) : 0,93.

---

## Objectifs du Projet

### Partie A : Analyse du Cycle
* Calcul des points de cycle (températures et pressions de stagnation).
* Évaluation du travail des compresseurs et turbines.
* Calcul de la puissance moteur (HP) et de la consommation spécifique.

### Partie B : Conception de la Turbine (HPT)
* Établissement des triangles de vitesses et de la vitesse de rotation.
* Calcul de la distribution radiale (vortex libre) et des coefficients de pertes (AMDC modifiée).
* Dimensionnement géométrique (corde axiale, pas, nombre d'ailettes).
* Analyse de la durée de vie visée de 300 heures.

### Partie C : Performance Hors Conception
* Simulation d'une réduction de 20% du RPM à débit massique constant.
* Calcul des pertes d'incidence via les corrélations de Moustapha.
* Évaluation de la dégradation du rendement.

---

## Nomenclature des Stations
* **Station 1** : Entrée LPC.
* **Station 2** : Entrée HPC.
* **Station 3** : Entrée Chambre de Combustion.
* **Station 4** : Entrée HPT.
* **Station 5** : Entrée LPT.
* **Station 6** : Entrée PT.
* **Station 7** : Sortie PT.
