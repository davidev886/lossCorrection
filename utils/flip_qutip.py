
import numpy as np
from qutip import *
from numpy import linalg as LA


def bit_flip_channel(p):
    X = sigmax()
    I = qeye(2)
    return (1-p) * sprepost(I, I)  + p * sprepost(X, X) 
    
print(to_choi(bit_flip_channel(0.3)))

chi_matrix = to_chi(to_choi(bit_flip_channel(0.3)))

print(chi_matrix.shape)
