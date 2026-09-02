import numpy as np
from dataclasses import dataclass, field


GAMMA = 0.1 # kg/m.s
N = 5
LENGTH = 1 # m
U = 0.1 # m/s
rho = 1 # kg/m3

WEST_BC_TYPE = "DIRICHLET" # DIRICHLET ou NEUMANN
WEST_BC_VAL = 1
EAST_BC_TYPE = "DIRICHLET" # DIRICHLET ou NEUMANN
EAST_BC_VAL = 0 

@dataclass
class BC():
    type: str = WEST_BC_TYPE
    val: float = WEST_BC_VAL

@dataclass
class Settings():
    gamma: float = GAMMA
    n: int = 5
    length: float = LENGTH
    u: float = U
    density: float = rho
    left_BC: BC = field(default_factory=BC)
    right_BC: BC = field(default_factory=BC)

def init_arrays(param: Settings):
    S = np.zeros(param.n)
    A = np.zeros((param.n,param.n))

    dx = param.length / param.n
    X = np.array([dx / 2 + i * dx for i in range(param.n)])

    return S, A, dx, X

def set_inner_cells(S, A, dx, param: Settings):
    for i in range (1,param.n-1):
        a_w = param.gamma / dx + param.density*param.u / 2
        a_e = param.gamma / dx - param.density*param.u / 2
        A[i, i-1] = -a_w # a_w
        A[i, i] = a_w + a_e # a_p
        A[i, i+1] = -a_e # a_e

def set_BC(S, A, dx, param: Settings, left=True):
    if left:
        type = param.left_BC.type
        val = param.left_BC.val
    else:
        type = param.right_BC.type
        val = param.right_BC.val
    F = param.density * param.u # En 1D, F_e = F_w = F_A = F_B
    D = param.gamma / dx
    if type == "DIRICHLET":
        if left:
            s_p = -(2 * D + F)
            s_u = (2 * D + F) * val
        else:
            s_p = -(2 * D - F)
            s_u = (2 * D - F) * val
    elif type == "NEUMANN":
        s_p = 0
        if left:
            s_u = (2 * D + F) * val * dx / 2
        else:
            s_u = (2 * D - F) * val * dx / 2
    else:
        raise(TypeError("Invalid boundary type"))

    if left:
        a_w = 0
        a_e = D - F / 2
        A[0, 0] = a_e + a_w - s_p # a_p
        A[0, 1] = - a_e
        S[0] = s_u
    else:
        a_w = D + F / 2
        a_e = 0
        A[-1, -1] = a_e + a_w - s_p # a_p
        A[-1, -2] = - a_w
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



