import numpy as np
import qutip as qu
from itertools import product
from p_operators_qutrit import *



def deterministic_correct(losses, stab_qubit_list, stab_outcome):
    if len(losses) == 1:
        print(losses[0])
    elif len(losses) == 2:
        
        
    for j_outcome, outcome in enumerate(stab_outcome):
        if outcome == -1:
            affected_stabilizer = stab_qubit_list[j_outcome]
            causing_qubits = list(set(affected_stabilizer) & set(losses)) 
            print("affected_stabilizer", affected_stabilizer)
            print("causing_qubits", causing_qubits)
            

stab_qubits = [[0,1,2,3], [1,2,4,5], [2,3,5,6]]
    

losses = [0,1, 6]


for x, y, z in product([1,-1],[1,-1],[1,-1]):
    print(x, y, z)
    deterministic_correct(losses, stab_qubits, [x,y,z])
    print()