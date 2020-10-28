import matplotlib.pyplot as plt
    
import numpy as np
#from qutip import *
import qutip as qu

from itertools import combinations, product

from utils.p_operators import *
from utils.binary_conf import *
from utils.corrections import *

import os
if not os.path.exists("data"):
    os.makedirs("data")

binary_conf = binary_configurations()    
#print(list(binary_conf.generate_configuration_loss_qnderror(0, 1)))
count_correctable_dic = {}
for num_l in range(8):
    for num_q in range(5-num_l):
        count_correctable = 0
        print(num_l, num_q)
        loss_error_conf = binary_conf.generate_configuration_loss_qnderror(num_l, num_q)
        for random_losses, qnd_errors in loss_error_conf:
            [correctable, non_correctable] = check_correctable_state_analytics(random_losses, qnd_errors)
            print(random_losses, qnd_errors, [correctable, non_correctable])
            if correctable:
                count_correctable += 1
        count_correctable_dic[(num_l, num_q)] = count_correctable
print(count_correctable_dic)


count_correctable_array = []
for (num_l, num_q) in count_correctable_dic:
    count_correctable_array.append([num_l, num_q, count_correctable_dic[(num_l, num_q)]])
    print(f"{count_correctable_dic[(num_l, num_q)]} * p^{num_l} * (1-p)^{7-num_l} * q^{num_q} * (1-q)^{7-num_q} +",)


np.savetxt("count_correctable_array.dat", count_correctable_array, fmt="%10d\t%10d\t%10d")

exit()
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
index_same_hit_pq2 = 0
index_log_op_pq2 = 0
for random_losses, qnd_errors in product(allowed_states_single, allowed_states_double):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if non_correctable_event:
        index_same_hit_pq2 += 1
        print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_same_hit_pq2,)


print()       
for random_losses, qnd_errors in product(allowed_states_single, allowed_states_double):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if to_check_event:
        if sorted(position_loss + position_qnd) in non_correctable_3_events:
            index_log_op_pq2 += 1
            print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_log_op_pq2,)

print(index_same_hit_pq2)
print(index_log_op_pq2)
print()
print()
index_same_hit_p2q = 0
index_log_op_p2q = 0
for random_losses, qnd_errors in product(allowed_states_double, allowed_states_single):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if non_correctable_event:
        index_same_hit_p2q += 1
        print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_same_hit_p2q,)


print()       
for random_losses, qnd_errors in product(allowed_states_double, allowed_states_single):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if to_check_event:
        if sorted(position_loss + position_qnd) in non_correctable_3_events:
            index_log_op_p2q += 1
            print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_log_op_p2q,)

            
print(index_same_hit_p2q)                
print(index_log_op_p2q)


print(f"{index_single_errors}* p * q + {index_same_hit_p2q+index_log_op_p2q}* p**2 * q + {index_same_hit_pq2 +index_log_op_pq2}* p * q**2")



exit()
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
index_same_hit_pq2 = 0
index_log_op_pq2 = 0
for random_losses, qnd_errors in product(allowed_states_single, allowed_states_double):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if non_correctable_event:
        index_same_hit_pq2 += 1
        print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_same_hit_pq2,)


print()       
for random_losses, qnd_errors in product(allowed_states_single, allowed_states_double):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if to_check_event:
        if sorted(position_loss + position_qnd) in non_correctable_3_events:
            index_log_op_pq2 += 1
            print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_log_op_pq2,)

print(index_same_hit_pq2)
print(index_log_op_pq2)
print()
print()
index_same_hit_p2q = 0
index_log_op_p2q = 0
for random_losses, qnd_errors in product(allowed_states_double, allowed_states_single):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if non_correctable_event:
        index_same_hit_p2q += 1
        print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_same_hit_p2q,)


print()       
for random_losses, qnd_errors in product(allowed_states_double, allowed_states_single):
    [correctable_event, non_correctable_event, to_check_event] = check_correctable_state(random_losses, qnd_errors)
    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    if to_check_event:
        if sorted(position_loss + position_qnd) in non_correctable_3_events:
            index_log_op_p2q += 1
            print(random_losses, qnd_errors, sorted(position_loss), sorted(position_qnd), sorted(position_loss + position_qnd), check_correctable_state(random_losses, qnd_errors), index_log_op_p2q,)

            
print(index_same_hit_p2q)                
print(index_log_op_p2q)


print(f"{index_single_errors}* p * q + {index_same_hit_p2q+index_log_op_p2q}* p**2 * q + {index_same_hit_pq2 +index_log_op_pq2}* p * q**2")