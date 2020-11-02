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

#

binary_conf = binary_configurations_loss_qnd_stab()    

# loss_error_stab_conf = binary_conf.generate_configuration_loss_qnderror(0, 1, 2)
# 
# for random_losses, qnd_errors, stab_errors_binary  in loss_error_stab_conf:
#     print(random_losses, qnd_errors, stab_errors_binary )


count_correctable_dic = {}
for num_l in range(8):
    for num_q in range(5-num_l):
        for num_s in range(4):
            count_correctable = 0
            loss_error_stab_conf = binary_conf.generate_configuration_loss_qnderror(num_l, num_q, num_s)

            for random_losses, qnd_errors, stab_errors_binary  in loss_error_stab_conf:

                [correctable, non_correctable] = check_stabilizers_measurement(random_losses, qnd_errors, stab_errors_binary)            
                print(random_losses, qnd_errors, stab_errors_binary, "correctable=", correctable)
                if correctable:
                    count_correctable += 1
            count_correctable_dic[(num_l, num_q, num_s)] = count_correctable
            print(num_l, num_q, num_s, count_correctable)            
print(count_correctable_dic)


count_correctable_array = []
for (num_l, num_q, num_s) in count_correctable_dic:
    count_correctable_array.append([num_l, num_q, num_s, count_correctable_dic[(num_l, num_q, num_s)]])
    print(f"{count_correctable_dic[(num_l, num_q, num_s)]} * p^{num_l} * (1-p)^{7-num_l} * q^{num_q} * (1-q)^{7-num_q} * s^{num_s} * (1-s)^{3-num_s} +",)


np.savetxt("count_correctable_stab_measurement_array.dat", count_correctable_array, fmt="%10d\t%10d\t%10d\t%10d")

exit()