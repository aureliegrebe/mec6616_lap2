# Résoudre les examples 4.1, 4.2 et 4.3 de Versteeg

import numpy as np
from solveur import solve, Settings, BC
import matplotlib.pyplot as plt

N_LIST = [5, 10, 20, 50, 100]

def ex4_1():
    param = Settings()
    param.right_BC = BC("DIRICHLET", 500)
    T = solve(param)

    plt.figure("ex 4.1")
    for n in N_LIST:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.legend()

def ex4_2():
    param = Settings()
    param.length = 0.02
    param.therm_cond = 0.5
    param.aire = 1
    param.q_u = 1e6
    param.left_BC = BC("DIRICHLET", 100)
    param.right_BC = BC("DIRICHLET", 200)
    T = solve(param)

    plt.figure("ex 4.2")
    for n in N_LIST:
        param.n = n
        X, T = solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.legend()

def ex4_3():
    param = Settings()
    param.length = 1
    param.therm_cond = 1
    param.aire = 1
    param.q_u = 25*20 # n^2 * T_inf
    param.q_p = -25 # -n^2
    param.left_BC = BC("DIRICHLET", 100)
    param.right_BC = BC("NEUMANN", 0)
    T = solve(param)

    plt.figure("ex 4.3")
    for n in N_LIST:
        param.n = n
        X, T , A= solve(param)
        plt.plot(X, T, label=f"n = {n}")
    plt.legend()
    plt.xlabel('x (en m)')
    plt.ylabel('T (en °C)')    
    plt.title('Température le long de la barre')

def main():
    ex4_3() # à changer en fonction de l'exemple
    plt.show()

if __name__ == '__main__':
    main()


