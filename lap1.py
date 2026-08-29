# Résoudre les examples 4.1, 4.2 et 4.3 de Versteeg

import numpy as np
from solveur import  solve, Settings, BC
import matplotlib.pyplot as plt
from cycler import cycler

# Parametre des figures
plt.rcParams.update({"font.size": 22})
plt.rc("legend", fontsize=15)
plt.rcParams["figure.figsize"] = [10.5, 8.5]

# Nombre de points pour le calcul de l'erreur
N_LIST = [5, 10, 20, 50, 100, 200, 500, 1000]
# Nombre de points pour la résolution utilisant la méthode FVM
N_LIST_PLOT = [5, 10, 20]

# Ex 4.1 
THERM_COND1 = 1000 # (W / m K)
AIRE1 = 10e-3 # m²
L1 = 0.5 # m
T_A1 = 100 # °C
T_B1 = 500 # °C

# Ex 4.2
T_A2 = 100 # °C
T_B2 = 200 # °C
L2 = 0.02 # m
AIRE2 = 1 # m²
THERM_COND2 = 0.5 # W / mK
Q2 = 1e6 # W/m^3

# Ex 4.3
T_A3 = 100 # °C
Q_B3 = 0 # Car extremité isolée donc pas de transfert de chaleur
AIRE3 = 1 # m²
L3 = 1 # m
THERM_COND3 = 1 # Choisi pour s'adapter à l'utilisation de la constante n sans changer la fonction définie dans solveur.py
T_inf = 20 # °C 
n3 = 5 # Car n² = 25
QU_3 = n3**2 * T_inf # D'après le tableau suivant l'équation (4.46) dans le livre de Versteeg
QP_3 = -n3**2 # Idem

# Définition des solution analytiques pour chaque exemple

N_ANALYTIQUE = 100

sol_analytique_4_1 = lambda x: 800 * x + 100
sol_analytique_4_2 = lambda x: ((T_B2 - T_A2) / L2 + Q2 / 2 / THERM_COND2 * (L2 - x))* x + T_A2
sol_analytique_4_3 = lambda x: (T_A3 - T_inf) * ((np.cosh(n3 * (L3 - x)))/np.cosh(n3 * L3)) + T_inf

# Exemple 4.1

def ex4_1():
    param = Settings()
    param.therm_cond = THERM_COND1
    param.aire = AIRE1
    param.left_BC = BC("DIRICHLET", T_A1)
    param.right_BC = BC("DIRICHLET", T_B1)
    param.q_u = 0
    param.q_p = 0
    param.length = L1

    plt.figure("ex 4.1")
    
    # Pour changer automatiquelent la représentation de chaque résultat sur la même courbe
    plt.rc(
        "axes",
        prop_cycle=(
            cycler("color", ["#000000", "#ff6347", "#1f77b4", "#2ca02c", "#d62728"])
            + cycler("ls", ["-", " ", " ", " ", " "])
        )
        + cycler("marker", [" ", "v", "8", "p", "D"]),
    )
    
    # Tracé des résultats par VFM et la solution analytique
    
    x = np.linspace(0, L1, N_ANALYTIQUE)
    plt.plot(x, sol_analytique_4_1(x), label="Solution analytique")
    for n in N_LIST_PLOT:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.xlabel("x (m)")
    plt.ylabel("T (°C)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("./figures/ex4_1.png")
    
    # Calcul de l'erreur en fonction du nombre de points

    epsilon2 = []
    h = []
    for n in N_LIST:
        param.n = n
        h.append(1/n)
        X, T = solve(param)
        epsilon2.append(np.sqrt(((T - sol_analytique_4_1(X))**2).mean()))
    plt.figure("err 4.1")
    plt.plot(h, epsilon2)
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.xlabel("h (m)")
    plt.ylabel("ε (°C)")
    plt.tight_layout()
    plt.savefig("./figures/epsilon4_1.png")

# Exemple 4.2

def ex4_2():
    param = Settings()
    param.length = L2
    param.therm_cond = THERM_COND2
    param.aire = AIRE2
    param.q_u = Q2
    param.left_BC = BC("DIRICHLET", T_A2)
    param.right_BC = BC("DIRICHLET", T_B2)
    T = solve(param)

    # Tracé des résultats par VFM et la solution analytique

    plt.figure("ex 4.2")
    x = np.linspace(0, L2, N_ANALYTIQUE)
    plt.plot(x, sol_analytique_4_2(x), label="Solution analytique")
    for n in N_LIST_PLOT:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.xlabel("x (m)")
    plt.ylabel("T (°C)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("./figures/ex4_2.png")

    # Calcul de l'erreur en fonction du nombre de points

    epsilon2 = []
    h = []
    for n in N_LIST:
        param.n = n
        h.append(1/n)
        X, T = solve(param)
        epsilon2.append(np.sqrt(((T - sol_analytique_4_2(X))**2).mean()))
    plt.figure("err 4.2")
    plt.plot(h, epsilon2)
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.xlabel("h (m)")
    plt.ylabel("ε (°C)")
    plt.tight_layout()
    plt.savefig("./figures/epsilon4_2.png")

    # Calcul d'ordre de convergence:
    p = float(np.log(epsilon2[-1]/epsilon2[-2])/np.log(N_LIST[-2]/N_LIST[-1]))
    print(f"Ordre de convergence pour ex 4.2: {p:.3}")

#Exemple 4.3

def ex4_3():
    param = Settings()
    param.length = L3
    param.therm_cond = THERM_COND3
    param.aire = AIRE3
    param.q_u = QU_3
    param.q_p = QP_3
    param.left_BC = BC("DIRICHLET", T_A3)
    param.right_BC = BC("NEUMANN", Q_B3)
    T = solve(param)

    # Tracé des résultats par VFM et la solution analytique

    plt.figure("ex 4.3")
    x = np.linspace(0, L3, N_ANALYTIQUE)
    plt.plot(x, sol_analytique_4_3(x), label="Solution analytique")
    for n in N_LIST_PLOT:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.xlabel("x (m)")
    plt.ylabel("T (°C)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("./figures/ex4_3.png")
    
    # Calcul de l'erreur en fonction du nombre de points
    
    epsilon2 = []
    h = []
    for n in N_LIST:
        param.n = n
        h.append(1/n)
        X, T = solve(param)
        epsilon2.append(np.sqrt(((T - sol_analytique_4_3(X))**2).mean()))
    plt.figure("err 4.3")
    plt.plot(h, epsilon2)
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.xlabel("h (m)")
    plt.ylabel("ε (°C)")
    plt.tight_layout()
    plt.savefig("./figures/epsilon4_3.png")

    # Calcul d'ordre de convergence:
    p = float(np.log(epsilon2[-1]/epsilon2[-2])/np.log(N_LIST[-2]/N_LIST[-1]))
    print(f"Ordre de convergend pour ex 4.3: {p:.3}")


def main():
    ex4_1()
    ex4_2()
    ex4_3()
    plt.show()

if __name__ == '__main__':
    main()


