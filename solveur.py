# Example 4.1 de Versteeg
import numpy as np
from dataclasses import dataclass, field


THERM_COND = 1000 # (W / m K)
AIRE = 1 # m²
N = 5
LENGTH = 1 # m
Q_U = 0 # W/m³
Q_P = 0 # W/m³

WEST_BC_TYPE = "DIRICHLET" # DIRIECHLET ou NEUMANN
WEST_BC_VAL = 100 # C ou C/m
EAST_BC_TYPE = "NEUMANN" # 
EAST_BC_VAL = 0 # C ou C/m

@dataclass
class BC():
    type: str = WEST_BC_TYPE
    val: float = WEST_BC_VAL

@dataclass
class Settings():
    therm_cond: float = THERM_COND
    aire: float = AIRE
    n: int = 5
    length: float = LENGTH
    q_u: float = Q_U
    q_p: float = Q_P
    left_BC: BC = field(default_factory=BC)
    right_BC: BC = field(default_factory=BC)
    conv_scheme: str = "UPWIND" # CENTRAL ou UPWIND
    rho_u: float = 0

def init_arrays(param: Settings):
    S = np.zeros(param.n)
    A = np.zeros((param.n,param.n))

    dx = param.length / param.n
    X = np.array([dx / 2 + i * dx for i in range(param.n)])

    return S, A, dx, X

def set_diff_coeffs(S, A, dx, param: Settings):
    for i in range (1,param.n-1):
        a_w = a_e = param.therm_cond * param.aire / dx
        s_u = param.q_u * param.aire * dx
        s_p = param.q_p * dx
        A[i, i-1] += -a_w # a_w
        A[i, i] += a_w + a_e - s_p # a_p
        A[i, i+1] += -a_e # a_e

        S[i] += s_u

def set_conv_coeffs_central(S, A, dx, param: Settings):
    # TODO: implement central scheme
    pass

def set_conv_coeffs_upwind(S, A, dx, param: Settings):
    for i in range (1,param.n-1):
        f_w = param.rho_u
        f_e = param.rho_u
        a_w = max(f_w, 0)
        a_e = max(0, -f_e)
        a_p = a_w + a_e + f_e - f_w
        A[i, i-1] += -a_w # a_w
        A[i, i] += a_p
        A[i, i+1] += -a_e # a_e


def set_inner_cells(S, A, dx, param: Settings):
    set_diff_coeffs(S, A, dx, param)
    if param.conv_scheme == "CENTRAL":
        set_conv_coeffs_central(S, A, dx, param)
    elif param.conv_scheme == "UPWIND":
        set_conv_coeffs_upwind(S, A, dx, param)
    else:
        raise(TypeError("Invalid convective scheme!!!"))

def set_BC(S, A, dx, param: Settings, left=True):
    if left:
        type = param.left_BC.type
        val = param.left_BC.val
    else:
        type = param.right_BC.type
        val = param.right_BC.val
    a_in = param.therm_cond * param.aire / dx
    if type == "DIRICHLET":
        s_u = param.q_u * param.aire * dx + 2 * param.therm_cond * param.aire * val / dx
        s_p = param.q_p * dx - 2 * param.therm_cond * param.aire / dx
    elif type == "NEUMANN":
        s_u = param.q_u * param.aire * dx + param.aire * val
        s_p = param.q_p * dx
    else:
        raise(TypeError("Invalid boundary type!!!"))

    if left:
        A[0, 0] = a_in - s_p
        A[0, 1] = -a_in
        S[0] = s_u
    else:
        A[-1, -1] = a_in - s_p
        A[-1, -2] = -a_in
        S[-1] = s_u

def solve(settings: Settings):
    # initialise les arrays
    S, A, dx, X = init_arrays(settings)

    # cellules intérieurs
    set_inner_cells(S, A, dx, settings)

    # CL gauche
    set_BC(S, A, dx, settings, left=True)
    set_BC(S, A, dx, settings, left=False)

    T = np.linalg.solve(A, S)

    return X, T

def main():
    param = Settings()
    param.right_BC = BC(EAST_BC_TYPE, EAST_BC_VAL)
    T = solve(param)

    # Post-processing
    print(T)


if __name__ == '__main__':
    main()



