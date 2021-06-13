import numpy as np
import qutip as qu
from utils.p_operators_qutrit import *


def channel_E_0(rho_L, prob, n_qubit):
    Xq = X[n_qubit]
    rhof = (prob['1a'] * Id * rho_L * Id +
            prob['1b'] * Xq * rho_L  * Xq +
            prob['1c'] * Xq * Xa * rho_L * Xq * Xa +
            prob['1d'] * Xa * rho_L * Xa
            )
    return rhof


def channel_E_1(rho_L, prob, n_qubit):
    Xq = X[n_qubit]
    rhof = (prob['2a'] * Id * rho_L * Id +
            prob['2b'] * Xa * rho_L * Xa
            )
    return rhof


def new_channel(rho_L, prob, n_qubit):
    dimHa = 2
    dimHq = 3

    proj_2 = qu.tensor([qu.qeye(dimHq)] * n_qubit
                        + [proj(2, 2, dimH = 3)]
                        + [qu.qeye(dimHq)] * (L - n_qubit - 1)
                        + [qu.qeye(dimHa)])
    P01 =  qu.tensor([qu.qeye(dimHq)] * n_qubit
                      + [qu.qeye(dimHq) - proj(2, 2, dimH = 3)]
                      + [qu.qeye(dimHq)] * (L - n_qubit - 1)
                      + [qu.qeye(dimHa)])
    Xq = X[n_qubit]
#     print(proj_2)
#     print()
#     print(rho_L)
    rhof =  (prob['1a'] * Id * P01 * rho_L * P01.dag() * Id +
            prob['1b'] * Xq * P01 * rho_L * P01.dag()  * Xq +
            prob['1c'] * Xq * Xa * P01 * rho_L * P01.dag() * Xq * Xa +
            prob['1d'] * Xa * P01 * rho_L * P01.dag() * Xa +
            prob['2a'] * Xa * proj_2 *  rho_L * proj_2 * Xa +
            prob['2b'] * Xa * proj_2 *  rho_L * proj_2 * Xa
            )
    return rhof