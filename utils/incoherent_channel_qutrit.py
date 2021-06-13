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
    Xq = X[n_qubit]
    rhof = (prob['1a'] * Id * rho_L * Id +
            prob['1b'] * Xq * rho_L  * Xq +
            prob['1c'] * Xq * Xa * rho_L * Xq * Xa +
            prob['1d'] * Xa * rho_L * Xa +
            prob['2a'] * Id * rho_L * Id +
            prob['2b'] * Xa * rho_L * Xa
            )
    return rhof