import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu

from itertools import combinations

from utils.p_operators import *
from utils.qnd_detection import check_correctable_state

import os
if not os.path.exists("data"):
    os.makedirs("data")

p_qnd = 0.0
p_loss = 0.4
p_stab = 0.4
#       1 means lost, 0 means not lost
random_losses = np.random.binomial(1, p_loss, 7)
#       1 means qnd error
qnd_errors = np.random.binomial(1, p_qnd, 7)

stab_errors = np.random.binomial(1, p_stab, 3)


print(random_losses)
print(qnd_errors)
print(stab_errors)