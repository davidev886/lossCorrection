"""

`choi_repr` contains several functions for computing and managing the 
 the choi operator of a qubit-qutrit system

"""

import numpy as np
import qutip as qu

from itertools import product


def proj(ket, bra, dimH = 2):
    if isinstance(ket, str):
        states_ket = [int(_) for _ in ket]
        ket_s = qu.basis([dimH] * len(states_ket), states_ket)        
    elif isinstance(ket, int):
        states_ket = ket
        ket_s = qu.basis(dimH, states_ket)        
    if isinstance(bra, str):
        states_bra = [int(_) for _ in bra]
        bra_s = qu.basis([dimH] * len(states_bra), states_bra).dag()
    elif isinstance(bra, int):
        states_bra = bra
        bra_s = qu.basis(dimH, states_bra).dag()
           
    return ket_s * bra_s


def apply_choi(rho, choiState):
    dimTot = stateTry.shape[0]
    if rho.type == "ket":
        rho = (rho * rho.dag()).unit()


    finalState = sum([rho[i,j] * choiState[dimTot*i : dimTot*i + dimTot, dimTot*j: dimTot*j + dimTot] 
                   for i,j in product(range(dimTot), range(dimTot))])    
    
    finalState = qu.Qobj(inpt=finalState, dims=rho.dims)
    return finalState


if __name__ == "__main__":
    dimHa = 2
    dimHq = 3
    dimTot = dimHa * dimHq
    GammaState = sum([qu.tensor(qu.basis(dimHa, ja), qu.basis(dimHq, jq), qu.basis(dimHa, ja), qu.basis(dimHq, jq)) 
                            for ja,jq in product(range(dimHa), range(dimHq))])
              
    rhoGamma = GammaState * GammaState.dag()        

    A0 = proj(1, 1, dimHq) + proj(0, 0, dimHq)
    A1 = proj(2, 2, dimHq)

    U = qu.tensor(qu.qeye(dimHa), A0) + qu.tensor(qu.sigmax(), A1)

    choiState = qu.tensor(qu.qeye(dimHa), qu.qeye(dimHq), U) * rhoGamma * qu.tensor(qu.qeye(dimHa), qu.qeye(dimHq), U.dag())

    stateQubit = (np.random.random() * qu.basis(dimHq, 0) + np.random.random() * qu.basis(dimHq, 1)).unit()

    stateTry_ket = qu.tensor(qu.basis(dimHa, 0), stateQubit)
    stateTry = stateTry_ket * stateTry_ket.dag()

    finalState = apply_choi(stateTry, choiState)

    print(finalState)

    print(apply_choi(stateTry_ket, choiState))

