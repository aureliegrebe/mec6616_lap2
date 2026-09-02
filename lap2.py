# Résoudre les examples 5.1 et 5.2 de Versteeg

from math import exp
import numpy as np
from solveur import solve, Settings, BC
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

# Constantes
PHI_0 = 1.
PHI_L = 0.
L = 1.0 # m
GAMMA = 0.1 # kg/(ms)
RHO_U_I = 0.1 # kg/(m²s)
RHO_U_II = 2.5 # kg/(m²s)

# Définition des solution analytiques pour chaque exemple

N_ANALYTIQUE = 100

sol_analytique_5_1i = lambda x: PHI_0 + (np.exp(RHO_U_I * x / GAMMA) - 1) / \
    (np.exp(RHO_U_I * L / GAMMA) - 1) * (PHI_L - PHI_0)
sol_analytique_5_1ii = lambda x: PHI_0 + (np.exp(RHO_U_II * x / GAMMA) - 1) / \
    (np.exp(RHO_U_II * L / GAMMA) - 1) * (PHI_L - PHI_0)


def ex5_1i():
    param = Settings()
    param.diffusivity = GAMMA
    param.density = 1.0
    param.u = RHO_U_I
    param.left_BC = BC("DIRICHLET", PHI_0)
    param.right_BC = BC("DIRICHLET", PHI_L)
    param.length = L
    param.conv_scheme = "CENTRAL"

    plt.figure("ex 5.1i")
    
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
    
    x = np.linspace(0, L, N_ANALYTIQUE)
    plt.plot(x, sol_analytique_5_1i(x), label="Solution analytique")
    for n in N_LIST_PLOT:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.xlabel("x (m)")
    plt.ylabel("Φ")
    plt.legend()
    plt.tight_layout()
    plt.savefig("./figures/ex5_1i.png")
    
    # Calcul de l'erreur en fonction du nombre de points

    epsilon2 = []
    h = []
    for n in N_LIST:
        param.n = n
        h.append(1/n)
        X, T = solve(param)
        epsilon2.append(np.sqrt(((T - sol_analytique_5_1i(X))**2).mean()))
    plt.figure("err 5.1i")
    plt.plot(h, epsilon2)
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.xlabel("h (m)")
    plt.ylabel("ε")
    plt.tight_layout()
    plt.savefig("./figures/epsilon5_1i.png")

    # Calcul d'ordre de convergence:
    p = float(np.log(epsilon2[-1]/epsilon2[-2])/np.log(N_LIST[-2]/N_LIST[-1]))
    print(f"Ordre de convergence pour ex 5.1: {p:.5f}")

def ex5_1ii():
    param = Settings()
    param.diffusivity = GAMMA
    param.density = 1.0
    param.u = RHO_U_II
    param.left_BC = BC("DIRICHLET", PHI_0)
    param.right_BC = BC("DIRICHLET", PHI_L)
    param.length = L
    param.conv_scheme = "CENTRAL"

    plt.figure("ex 5.1ii")
    
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
    
    x = np.linspace(0, L, N_ANALYTIQUE)
    plt.plot(x, sol_analytique_5_1ii(x), label="Solution analytique")
    for n in N_LIST_PLOT:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.xlabel("x (m)")
    plt.ylabel("Φ")
    plt.legend()
    plt.tight_layout()
    plt.savefig("./figures/ex5_1ii.png")
    
    # Calcul de l'erreur en fonction du nombre de points

    epsilon2 = []
    h = []
    for n in N_LIST:
        param.n = n
        h.append(1/n)
        X, T = solve(param)
        epsilon2.append(np.sqrt(((T - sol_analytique_5_1ii(X))**2).mean()))
    plt.figure("err 5.1ii")
    plt.plot(h, epsilon2)
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.xlabel("h (m)")
    plt.ylabel("ε")
    plt.tight_layout()
    plt.savefig("./figures/epsilon5_1ii.png")

    # Calcul d'ordre de convergence:
    p = float(np.log(epsilon2[-1]/epsilon2[-2])/np.log(N_LIST[-2]/N_LIST[-1]))
    print(f"Ordre de convergence pour ex 5.1: {p:.5f}")

# Exemple 5.2

def ex5_2i():
    param = Settings()
    param.diffusivity = GAMMA
    param.density = 1.0
    param.u = RHO_U_I
    param.left_BC = BC("DIRICHLET", PHI_0)
    param.right_BC = BC("DIRICHLET", PHI_L)
    param.length = L
    param.conv_scheme = "UPWIND"

    plt.figure("ex 5.2i")
    
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
    
    x = np.linspace(0, L, N_ANALYTIQUE)
    plt.plot(x, sol_analytique_5_1i(x), label="Solution analytique")
    for n in N_LIST_PLOT:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.xlabel("x (m)")
    plt.ylabel("Φ")
    plt.legend()
    plt.tight_layout()
    plt.savefig("./figures/ex5_2i.png")
    
    # Calcul de l'erreur en fonction du nombre de points

    epsilon2 = []
    h = []
    for n in N_LIST:
        param.n = n
        h.append(1/n)
        X, T = solve(param)
        epsilon2.append(np.sqrt(((T - sol_analytique_5_1i(X))**2).mean()))
    plt.figure("err 5.2i")
    plt.plot(h, epsilon2)
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.xlabel("h (m)")
    plt.ylabel("ε")
    plt.tight_layout()
    plt.savefig("./figures/epsilon5_2i.png")

    # Calcul d'ordre de convergence:
    p = float(np.log(epsilon2[-1]/epsilon2[-2])/np.log(N_LIST[-2]/N_LIST[-1]))
    print(f"Ordre de convergence pour ex 5.2: {p:.5f}")

def ex5_2ii():
    param = Settings()
    param.diffusivity = GAMMA
    param.density = 1.0
    param.u = RHO_U_II
    param.left_BC = BC("DIRICHLET", PHI_0)
    param.right_BC = BC("DIRICHLET", PHI_L)
    param.length = L
    param.conv_scheme = "UPWIND"

    plt.figure("ex 5.2ii")
    
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
    
    x = np.linspace(0, L, N_ANALYTIQUE)
    plt.plot(x, sol_analytique_5_1ii(x), label="Solution analytique")
    for n in N_LIST_PLOT:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.xlabel("x (m)")
    plt.ylabel("Φ")
    plt.legend()
    plt.tight_layout()
    plt.savefig("./figures/ex5_2ii.png")
    
    # Calcul de l'erreur en fonction du nombre de points

    epsilon2 = []
    h = []
    for n in N_LIST:
        param.n = n
        h.append(1/n)
        X, T = solve(param)
        epsilon2.append(np.sqrt(((T - sol_analytique_5_1ii(X))**2).mean()))
    plt.figure("err 5.2ii")
    plt.plot(h, epsilon2)
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.xlabel("h (m)")
    plt.ylabel("ε")
    plt.tight_layout()
    plt.savefig("./figures/epsilon5_2ii.png")

    # Calcul d'ordre de convergence:
    p = float(np.log(epsilon2[-1]/epsilon2[-2])/np.log(N_LIST[-2]/N_LIST[-1]))
    print(f"Ordre de convergence pour ex 5.2: {p:.5f}")



def main():
    ex5_1i()
    ex5_1ii()
    ex5_2i()
    ex5_2ii()
    plt.show()

if __name__ == '__main__':
    main()


