# Example 4.1 de Versteeg
import numpy as np

THERM_COND = 1000 # (W / m K)
AIRE = 10e-3 # m²
N = 5
LENGTH = 0.5 # m
Q_U = 0 # W/m³
Q_P = 0 # W/m³

# Initial temperature
T_0 = 100 # C

WEST_BC_TYPE = "DIRICHLET" # DIRIECHLET ou NEUMANN
WEST_BC_VAL = 100 # C ou C/m
EAST_BC_TYPE = "DIRICHLET" # 
EAST_BC_VAL = 500 # C ou C/m


def main():
    # initialise les arrays
    # T = np.ones(N) * T_0
    S = np.zeros(N)
    A = np.zeros((N,N))

    dx = LENGTH / N

    # cellules intérieurs
    for i in range (1,N-1):
        a_w = a_e = THERM_COND * AIRE / dx
        s_u = Q_U
        s_p = Q_P
        A[i, i-1] = -a_w # a_w
        A[i, i] = a_w + a_e - s_p # a_p
        A[i, i+1] = -a_e # a_e

        S[i] = s_u

    # Conditions aux frontières
    # west boundary
    a_e = THERM_COND * AIRE / dx
    if WEST_BC_TYPE == "DIRICHLET":
        s_u = Q_U + 2 * THERM_COND * AIRE * WEST_BC_VAL / dx
        s_p = Q_P - 2 * THERM_COND * AIRE / dx
    elif WEST_BC_TYPE == "NEUMANN":
        s_u = Q_U + AIRE * WEST_BC_VAL
        s_p = Q_P
    else:
        raise(TypeError("Invalid boundary type"))

    A[0, 0] = a_e - s_p
    A[0, 1] = -a_e
    S[0] = s_u

    # east boundary
    a_w = THERM_COND * AIRE / dx
    if EAST_BC_TYPE == "DIRICHLET":
        s_u = Q_U + 2 * THERM_COND * AIRE * EAST_BC_VAL / dx
        s_p = Q_P - 2 * THERM_COND * AIRE / dx
    elif EAST_BC_TYPE == "NEUMANN":
        s_u = Q_U + AIRE * EAST_BC_VAL
        s_p = Q_P
    else:
        raise(TypeError("Invalid boundary type"))

    A[-1, -1] = a_e - s_p
    A[-1, -2] = -a_w
    S[-1] = s_u

    # Resolution
    T = np.linalg.solve(A, S)

    # Post-processing
    print(A)
    print(S)
    print(T)



if __name__ == '__main__':
    main()

