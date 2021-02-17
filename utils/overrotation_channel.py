import numpy as np
import qutip as qu


L = 7

def DepolQubit(i,p):  # Qubit Kraus-operators for depolarizing channel
    dimHq = 3
    index = [
            np.sqrt(1 - 3 * p / 4) * qu.qeye(2),
            np.sqrt(p / 4) * qu.sigmax(),
            np.sqrt(p / 4) * qu.sigmay(),
            np.sqrt(p / 4) * qu.sigmaz()
            ]
    return qu.tensor([qu.qeye(dimHq)] * 7 + [index[i]])

def DepolQutrit(qutrit_n, i, p):  # Qutrit Kraus-operators for depolarizing channel
    dimHq = 3
    dimHa = 2
    index = [
            np.sqrt(3 - 8 * p / 3) / np.sqrt(3) * qu.Qobj([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dims=[[3], [3]]),
            np.sqrt(p / 3) / np.sqrt(2) * qu.Qobj([[1, 0, 0], [0, 0, 0], [0, 0, -1]], dims=[[3], [3]]),
            np.sqrt(p / 3) * qu.Qobj([[0, 0, 0], [0, 0, 1], [0, 0, 0]], dims=[[3], [3]]),
            np.sqrt(p / 3) * qu.Qobj([[0, 0, 1], [0, 0, 0], [0, 0, 0]], dims=[[3], [3]]),
            np.sqrt(p / 3) * qu.Qobj([[0, 0, 0], [0, 0, 0], [0, 1, 0]], dims=[[3], [3]]),
            np.sqrt(p / 3) / np.sqrt(6) * qu.Qobj([[1, 0, 0], [0, -2, 0], [0, 0, 1]], dims=[[3], [3]]),
            np.sqrt(p / 3) * qu.Qobj([[0, 1, 0], [0, 0, 0], [0, 0, 0]], dims=[[3], [3]]),
            np.sqrt(p / 3) * qu.Qobj([[0, 0, 0], [0, 0, 0], [1, 0, 0]], dims=[[3], [3]]),
            np.sqrt(p / 3) * qu.Qobj([[0, 0, 0], [1, 0, 0], [0, 0, 0]], dims=[[3], [3]])
            ]
    temp = [qu.qeye(dimHq)] * qutrit_n + [index[i]] + [qu.qeye(dimHq)] * (L - qutrit_n - 1) + [qu.qeye(dimHa)]

    return qu.tensor(temp)
    
    
def UnitaryQubitQutritQNDDepol(p, state_rho, qutrit_n): # depolarizing channel on ancilla qubit & system qutrit
    Depol = 0
    for i in range(4):
        for j in range(9):
            Depol = (Depol + 
                DepolQubit(i,p) * DepolQutrit(qutrit_n, j,p) * state_rho 
                * DepolQutrit(qutrit_n, j,p).dag() * DepolQubit(i,p).dag() )

    return Depol / np.trace(Depol.data.toarray())

def CorrelatedOverRotQubit(qutrit_n, p_corr): 
    dimHq = 3
    dimHa = 2
    X_qutrit = qu.Qobj([[0, 1 ,0], [1, 0, 0], [0, 0, 0]])
    ket22bra = qu.Qobj([[0, 0 ,0], [0, 0 ,0], [0, 0, 1]])
    Id =  [qu.qeye(dimHq)] * L + [qu.qeye(dimHa)]

    XX_operators_1 = [qu.qeye(dimHq)] * qutrit_n + [X_qutrit] + [qu.qeye(dimHq)] * (L - qutrit_n - 1) + [qu.sigmax()]
    XX_operators_2 = [qu.qeye(dimHq)] * qutrit_n + [ket22bra] + [qu.qeye(dimHq)] * (L - qutrit_n - 1) + [qu.qeye(dimHa)]

    corr = np.cos(p_corr * np.pi / 2) * qu.tensor(Id) + \
            1j * np.sin(p_corr * np.pi / 2) * (qu.tensor(XX_operators_1) + qu.tensor(XX_operators_2))

    return corr
    
def CorrelatedOverRotQubitAll(p_corr):
    return [CorrelatedOverRotQubit(qutrit_n, p_corr) for qutrit_n in range(L)]

if __name__ == "__main__":
    print("c")
    p_corr = 0.5

    vacuum = qu.tensor([qu.basis(3,0)] * L + [qu.basis(2,0)])
    X_o = CorrelatedOverRotQubitAll(p_corr)
    print(len(X_o))
#     state_rho = vacuum*vacuum.dag()
#     for data_q in range(7):
#         state_rho = UnitaryQubitQutritQNDDepol(p, state_rho, data_q)
#         
