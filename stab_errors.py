import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu

from itertools import combinations

from utils.p_operators import *
from utils.corrections import check_correctable_state

import os
if not os.path.exists("data"):
    os.makedirs("data")

p_qnd = 0.2
p_loss = 0.2
p_stab = 0.4
#       1 means lost, 0 means not lost
random_losses = np.random.binomial(1, p_loss, 7)
#       1 means qnd error
qnd_errors = np.random.binomial(1, p_qnd, 7)

stab_errors_binary = np.random.binomial(1, p_stab, 3)
faulty_stab_qubits = [el for j,el in enumerate(stab_qubits) if stab_errors_binary[j]]

position_loss = np.where(random_losses)[0].tolist()
position_qnd = np.where(qnd_errors)[0].tolist()

guessed_loss = [(loss + qnd_err) % 2 for loss, qnd_err in zip(random_losses, qnd_errors)]
num_fresh_qubits = sum(guessed_loss)
position_guessed_loss = np.where(guessed_loss)[0].tolist() 
qnderror_hit_loss = any([(_ in position_loss) for _ in position_qnd])    

print(f"{'random_losses':40}", random_losses)
print(f"{'qnd_errors':40}", qnd_errors)
print(f"{'stab_errors_binary':40}", stab_errors_binary)
print(f"{'position_guessed_loss':40}", position_guessed_loss)
print(f"{'faulty_stab_qubits':40}", faulty_stab_qubits)

if qnderror_hit_loss:
    non_correctable = True
    correctable = False
    print("qnderror_hit_loss")
    print(f"{'correctable':40}",  correctable)    
else:
    correctable = not (any([any([(loss in stab) for loss in position_guessed_loss]) for stab in faulty_stab_qubits]))
    print(f"{'correctable':40}",  correctable)
#    print(f"{'not_correctable':40}", not correctable)