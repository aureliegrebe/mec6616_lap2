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

def init_arrays(n=N, l=LENGTH):
    S = np.zeros(n)
    A = np.zeros((n,n))

    dx = float(l) / n

    return S, A, dx

def set_inner_cells(S, A, dx, n=N, k=THERM_COND, aire=AIRE, q_u=Q_U, q_p=Q_P):
    for i in range (1,N-1):
        a_w = a_e = THERM_COND * AIRE / dx
        s_u = Q_U
        s_p = Q_P
        A[i, i-1] = -a_w # a_w
        A[i, i] = a_w + a_e - s_p # a_p
        A[i, i+1] = -a_e # a_e

        S[i] = s_u

def set_BC(S, A, dx, type, val, left, k=THERM_COND, aire=AIRE, q_u=Q_U,
           q_p=Q_P):
    a_in = k * aire / dx
    if type == "DIRICHLET":
        s_u = q_u + 2 * k * aire * val / dx
        s_p = q_p - 2 * k * aire / dx
    elif type == "NEUMANN":
        s_u = q_u + aire * val
        s_p = q_p
    else:
        raise(TypeError("Invalid boundary type"))

    if left:
        A[0, 0] = a_in - s_p
        A[0, 1] = -a_in
        S[0] = s_u
    else:
        A[-1, -1] = a_in - s_p
        A[-1, -2] = -a_in
        S[-1] = s_u

def main():
    # initialise les arrays
    S, A, dx = init_arrays(N, LENGTH)

    # cellules intérieurs
    set_inner_cells(S, A, dx)

    # Conditions aux frontières
    set_BC(S, A, dx, WEST_BC_TYPE, WEST_BC_VAL, True)
    set_BC(S, A, dx, EAST_BC_TYPE, EAST_BC_VAL, False)
    # # west boundary
    # a_e = THERM_COND * AIRE / dx
    # if WEST_BC_TYPE == "DIRICHLET":
    #     s_u = Q_U + 2 * THERM_COND * AIRE * WEST_BC_VAL / dx
    #     s_p = Q_P - 2 * THERM_COND * AIRE / dx
    # elif WEST_BC_TYPE == "NEUMANN":
    #     s_u = Q_U + AIRE * WEST_BC_VAL
    #     s_p = Q_P
    # else:
    #     raise(TypeError("Invalid boundary type"))

    # A[0, 0] = a_e - s_p
    # A[0, 1] = -a_e
    # S[0] = s_u

    # # east boundary
    # a_w = THERM_COND * AIRE / dx
    # if EAST_BC_TYPE == "DIRICHLET":
    #     s_u = Q_U + 2 * THERM_COND * AIRE * EAST_BC_VAL / dx
    #     s_p = Q_P - 2 * THERM_COND * AIRE / dx
    # elif EAST_BC_TYPE == "NEUMANN":
    #     s_u = Q_U + AIRE * EAST_BC_VAL
    #     s_p = Q_P
    # else:
    #     raise(TypeError("Invalid boundary type"))

    # A[-1, -1] = a_e - s_p
    # A[-1, -2] = -a_w
    # S[-1] = s_u

    # Resolution
    T = np.linalg.solve(A, S)

    # Post-processing
    print(A)
    print(S)
    print(T)


if __name__ == '__main__':
    main()

