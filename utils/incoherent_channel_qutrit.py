import numpy as np
import qutip as qu
from utils.p_operators_qutrit import *

def inc_channel(prob, n_qubit, noloss=True):
    Xq = X[n_qubit]
    if noloss:
        channel = (prob['1a'] * qu.spre(Id) * qu.spost(Id) +
                   prob['1b'] * qu.spre(Xq) * qu.spost(Xq) +
                   prob['1c'] * qu.spre(Xq * Xa) * qu.spost(Xq * Xa) +
                   prob['1d'] * qu.spre(Xa) * qu.spost(Xa)
                   )
        return channel
    else:
        channel_0 = (prob['1a'] * qu.spre(Id) * qu.spost(Id) +
                   prob['1b'] * qu.spre(Xq) * qu.spost(Xq) +
                   prob['1c'] * qu.spre(Xq * Xa) * qu.spost(Xq * Xa) +
                   prob['1d'] * qu.spre(Xa) * qu.spost(Xa)
                   )

        channel_1 = (prob['2a'] * qu.spre(Id) * qu.spost(Id) +
                   prob['2b'] * qu.spre(Xa) * qu.spost(Xa)
                   )

        return channel_0, channel_1



if __name__ == "__main__":
    print("c")
    p = 1.00
    qutrit_n = 6
    vacuum = qu.tensor([qu.basis(3,0)] * L + [qu.basis(2,0)])

    state_rho = vacuum*vacuum.dag()
    for data_q in range(7):
        state_rho = UnitaryQubitQutritQNDDepol(p, state_rho, data_q)
