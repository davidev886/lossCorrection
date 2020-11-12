import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu

from utils.qnd_error_gen import pick_qnd_error

from utils.p_operators_qutrit import *

from random import randint
seme = randint(0,100)

seme=5
np.random.seed(seme)

print("\nseed:", seme)

a = np.random.random()  + np.random.random() * 1j
b = np.random.random()  + np.random.random() * 1j

psiL = (a * ZeroL + b * OneL).unit()

#w_0 = psiL.ptrace(kept_qubits)
#w = (qu.tensor([qu.fock_dm(2,0)] * len(losses) + [w_0])).unit()
phi = 0.1 * np.pi

R_loss_list = Rloss(phi)

for data_q in [0]:
    #apply Rloss with an angle phi
    print(R_loss_list[data_q] * psiL)