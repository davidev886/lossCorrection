from  qutip import *
import numpy as np
import scipy.linalg
from itertools import product
import matplotlib.pyplot as plt



x = sigmax()
y = sigmay()
z = sigmaz()

p = 1
Eps1 = ((1-3*p/4) * spre(qeye(2)) * spost(qeye(2).dag())  + 
        0.25 * p* (spre(x) * spost(x.dag()) + spre(y) * spost(y.dag()) + spre(z) * spost(z.dag())) )
Eps2 = Qobj(np.eye(4)/2, dims = Eps1.dims, type = Eps1.type)

state_s = (basis(2,0)+0.45 * basis(2,1)).unit()
print(Eps1)
print(Eps2)
print("to choi")
print(to_choi(Eps1))
print(Eps1(state_s))


choi = Qobj(np.eye(4)/2, dims = Eps1.dims, type = Eps1.type)

print("this", choi_to_super(choi))
