import numpy as np
import qutip as qu
from itertools import product

L = 7


def DepolQubit(i, p):  # Qubit Kraus-operators for depolarizing channel
    dimHq = 3
    index = [
            np.sqrt(1 - 3 * p / 4) * qu.qeye(2),
            np.sqrt(p / 4) * qu.sigmax(),
            np.sqrt(p / 4) * qu.sigmay(),
            np.sqrt(p / 4) * qu.sigmaz()
            ]
    return qu.tensor([qu.qeye(dimHq)] * 7 + [index[i]])


def DepolQutrit(qutrit_n, i, p):
    # Qutrit Kraus-operators for depolarizing channel
    dimHq = 3
    dimHa = 2
    dimTot = [[3], [3]]
    index = [
            np.sqrt(3 - 8 * p / 3) / np.sqrt(3) *
            qu.Qobj([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dims=dimTot),
            np.sqrt(p / 3) / np.sqrt(2) *
            qu.Qobj([[1, 0, 0], [0, 0, 0], [0, 0, -1]], dims=dimTot),
            np.sqrt(p / 3) *
            qu.Qobj([[0, 0, 0], [0, 0, 1], [0, 0, 0]], dims=dimTot),
            np.sqrt(p / 3) *
            qu.Qobj([[0, 0, 1], [0, 0, 0], [0, 0, 0]], dims=dimTot),
            np.sqrt(p / 3) *
            qu.Qobj([[0, 0, 0], [0, 0, 0], [0, 1, 0]], dims=dimTot),
            np.sqrt(p / 3) / np.sqrt(6) *
            qu.Qobj([[1, 0, 0], [0, -2, 0], [0, 0, 1]], dims=dimTot),
            np.sqrt(p / 3) *
            qu.Qobj([[0, 1, 0], [0, 0, 0], [0, 0, 0]], dims=dimTot),
            np.sqrt(p / 3) *
            qu.Qobj([[0, 0, 0], [0, 0, 0], [1, 0, 0]], dims=dimTot),
            np.sqrt(p / 3) *
            qu.Qobj([[0, 0, 0], [1, 0, 0], [0, 0, 0]], dims=dimTot)
            ]
    temp = ([qu.qeye(dimHq)] * qutrit_n +
            [index[i]] +
            [qu.qeye(dimHq)] * (L - qutrit_n - 1) +
            [qu.qeye(dimHa)]
            )

    return qu.tensor(temp)


def UnitaryQubitQutritQNDDepol(p, state_rho, qutrit_n):
    # depolarizing channel on ancilla qubit & system qutrit
    Depol = 0
    for i, j in product(range(4), range(9)):
        Depol = (Depol +
                 DepolQubit(i, p) *
                 DepolQutrit(qutrit_n, j, p) *
                 state_rho *
                 DepolQutrit(qutrit_n, j, p).dag() *
                 DepolQubit(i, p).dag()
                 )

    return Depol / np.trace(Depol.data.toarray())


def CorrelatedOverRotQubit(qutrit_n, alpha):
    dimHq = 3
    dimHa = 2
    X_qutrit = qu.Qobj([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    ket22bra = qu.Qobj([[0, 0, 0], [0, 0, 0], [0, 0, 1]])
    Id = [qu.qeye(dimHq)] * L + [qu.qeye(dimHa)]

    XX_operators_1 = ([qu.qeye(dimHq)] * qutrit_n +
                      [X_qutrit] +
                      [qu.qeye(dimHq)] * (L - qutrit_n - 1) +
                      [qu.sigmax()]
                      )
    XX_operators_2 = ([qu.qeye(dimHq)] * qutrit_n +
                      [ket22bra] +
                      [qu.qeye(dimHq)] * (L - qutrit_n - 1) +
                      [qu.qeye(dimHa)]
                      )

    corr = (np.cos(alpha / 2) * qu.tensor(Id) +
            1j * np.sin(alpha / 2) * (qu.tensor(XX_operators_1) +
            qu.tensor(XX_operators_2))
            )

    return corr


def CorrelatedOverRotQubitAll(alpha):
    return [CorrelatedOverRotQubit(qutrit_n, alpha) for qutrit_n in range(L)]


def SingleOverRotQubit(qutrit_n, theta):
    dimHq = 3
    dimHa = 2
    X_qutrit = qu.Qobj([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    ket22bra = qu.Qobj([[0, 0, 0], [0, 0, 0], [0, 0, 1]])

    R1q = (np.cos(theta / 2) * (qu.qeye(dimHq) - ket22bra) -
           1j * np.sin(theta / 2) * X_qutrit +
           ket22bra
           )
    R1a = (np.cos(theta / 2) * qu.qeye(dimHa) -
           1j * np.sin(theta / 2) * qu.sigmax()
           )

    OverRotSingle = ([qu.qeye(dimHq)] * qutrit_n +
                     [R1q] +
                     [qu.qeye(dimHq)] * (L - qutrit_n - 1) +
                     [R1a]
                     )
    corr = qu.tensor(OverRotSingle)

    return corr


def SingleOverRotQutritAll(theta):
    return [SingleOverRotQutrit(qutrit_n, theta) for qutrit_n in range(L)]


def SingleOverRotQutrit(qutrit_n, theta):
    dimHq = 3
    dimHa = 2
    X_qutrit = qu.Qobj([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    ket22bra = qu.Qobj([[0, 0, 0], [0, 0, 0], [0, 0, 1]])

    R1q = (np.cos(theta / 2) * (qu.qeye(dimHq) - ket22bra) -
           1j * np.sin(theta / 2) * X_qutrit +
           ket22bra
           )
    R1a = qu.qeye(dimHa)

    OverRotSingle = ([qu.qeye(dimHq)] * qutrit_n +
                     [R1q] +
                     [qu.qeye(dimHq)] * (L - qutrit_n - 1) +
                     [R1a]
                     )
    corr = qu.tensor(OverRotSingle)

    return corr


def SingleOverRotQubitAll(theta):
    return [SingleOverRotQubit(qutrit_n, theta) for qutrit_n in range(L)]



if __name__ == "__main__":
    print("c")
    p_corr = 0.5

    vacuum = qu.tensor([qu.basis(3, 0)] * L + [qu.basis(2, 0)])
    X_o = CorrelatedOverRotQubitAll(p_corr)
    print(len(X_o))
#     state_rho = vacuum*vacuum.dag()
#     for data_q in range(7):
#         state_rho = UnitaryQubitQutritQNDDepol(p, state_rho, data_q)
#
