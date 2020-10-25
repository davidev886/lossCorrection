import matplotlib.pyplot as plt
    
import numpy as np
#from qutip import *
import qutip as qu

from itertools import combinations, product

from utils.p_operators import *
from utils.qnd_detection import check_correctable_state

import os
if not os.path.exists("data"):
    os.makedirs("data")


n = 7
allowed_states_single = []
allowed_states_double = []
for i in range(2**7):
    state = [int(_) for _ in reversed(f"{bin(i)[2:].zfill(7)}")]
    if sum(state) == 1:
        allowed_states_single.append(state)
    elif sum(state) == 2:
        allowed_states_double.append(state)    
   
index_single_errors = 0
for random_losses, qnd_errors in product(allowed_states_single, allowed_states_single):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if non_correctable_event:
        index_single_errors += 1
        print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_single_errors,)

print()

non_correctable_3_events = [[0, 1, 4], [0, 2, 5], [0, 3, 6], [1, 2, 6], [2, 3, 4], [4, 5, 6], [1, 3, 5]]
index_same_hit = 0
index_log_op = 0
for random_losses, qnd_errors in product(allowed_states_single, allowed_states_double):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if non_correctable_event:
        index_same_hit += 1
        print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_same_hit,)


print()       
for random_losses, qnd_errors in product(allowed_states_single, allowed_states_double):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if to_check_event:
        if sorted(position_loss + position_qnd) in non_correctable_3_events:
            index_log_op += 1
            print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_log_op,)

            
print(index_same_hit)                
print(index_log_op)



print()
print()
index_same_hit = 0
index_log_op = 0
for random_losses, qnd_errors in product(allowed_states_double, allowed_states_single):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if non_correctable_event:
        index_same_hit += 1
        print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_same_hit,)


print()       
for random_losses, qnd_errors in product(allowed_states_double, allowed_states_single):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if to_check_event:
        if sorted(position_loss + position_qnd) in non_correctable_3_events:
            index_log_op += 1
            print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_log_op,)

            
print(index_same_hit)                
print(index_log_op)

