# Résoudre les examples 4.1, 4.2 et 4.3 de Versteeg

import numpy as np
from solveur import  solve, Settings, BC
import matplotlib.pyplot as plt
from cycler import cycler

# figure parameters
plt.rcParams.update({"font.size": 22})
plt.rc("legend", fontsize=15)
plt.rcParams["figure.figsize"] = [15, 8.5]

N_LIST = [5, 10, 20, 50, 100, 200, 500, 1000]
N_LIST_PLOT = [5, 10, 20]

# Ex 4.1 
THERM_COND1 = 1000 # (W / m K)
AIRE1 = 10e-3 # m²
L1 = 0.5 # m
T_A1 = 100
T_B1 = 500

# Ex 4.2
T_A2 = 100
T_B2 = 200
L2 = 0.02
THERM_COND2 = 0.5 # W / mK
Q2 = 1e6 # W/m^3

N_ANALYTIQUE = 100

sol_analytique_4_1 = lambda x: 800 * x + 100
sol_analytique_4_2 = lambda x: ((T_B2 - T_A2) / L2 + Q2 / 2 / THERM_COND2 * (L2 - x))* x + T_A2

def ex4_1():
    param = Settings()
    param.therm_cond = THERM_COND1
    param.aire = AIRE1
    param.left_BC = BC("DIRICHLET", T_A1)
    param.right_BC = BC("DIRICHLET", T_B1)

    plt.figure("ex 4.1")
    plt.rc(
        "axes",
        prop_cycle=(
            cycler("color", ["#000000", "#ff6347", "#1f77b4", "#2ca02c", "#d62728"])
            + cycler("ls", ["-", (0, (0, 1)), (0, (0, 1)),(0, (0, 1)),(0, (0, 1))])
        )
        + cycler("marker", [" ", "v", "8", "p", "D"]),
    )
    x = np.linspace(0, L1, N_ANALYTIQUE)
    plt.plot(x, sol_analytique_4_1(x), label="Solution analytique")
    for n in N_LIST_PLOT:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.xlabel("x (m)")
    plt.ylabel("T (°C)")
    plt.legend()

    epsilon2 = []
    dx = []
    for n in N_LIST:
        param.n = n
        dx.append(L1 / n)
        X, T = solve(param)
        epsilon2.append(np.sqrt((T - sol_analytique_4_1(X))**2).mean())
    plt.figure("err 4.1")
    plt.plot(dx, epsilon2)
    plt.xlabel("dx (m)")
    plt.ylabel("ΔT (°C)")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()


def ex4_2():
    param = Settings()
    param.length = L2
    param.therm_cond = THERM_COND2
    param.aire = 1
    param.q_u = Q2
    param.left_BC = BC("DIRICHLET", T_A2)
    param.right_BC = BC("DIRICHLET", T_B2)
    T = solve(param)

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

    epsilon2 = []
    dx = []
    for n in N_LIST:
        param.n = n
        dx.append(L1 / n)
        X, T = solve(param)
        epsilon2.append(np.sqrt((T - sol_analytique_4_2(X))**2).mean())
    plt.figure("err 4.2")
    plt.plot(dx, epsilon2)
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.xlabel("dx (m)")
    plt.ylabel("ΔT (°C)")

def ex4_3():
    param = Settings()
    param.length = 0.02
    param.therm_cond = 0.5
    param.aire = 1
    param.q_u = 1e6
    param.left_BC = BC("DIRICHLET", 100)
    param.right_BC = BC("NEUMANN", 0)
    T = solve(param)

    plt.figure("ex 4.3")
    for n in N_LIST:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.legend()

def main():
    ex4_1()
    ex4_2()
    plt.show()

if __name__ == '__main__':
    main()


