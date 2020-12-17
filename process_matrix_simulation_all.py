import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu

from utils.qnd_error_gen import pick_qnd_error
from utils.binary_conf import binary_configurations, binary_raw_configurations
from utils.p_operators_qutrit import *

import datetime
now = datetime.datetime.now()
final_data_name = now.strftime("%Y%m%d%H%M")

import argparse
#python process_matrix_simulation_all.py --phi_tilde  --epsilon_choi
parser = argparse.ArgumentParser(description = "Simulate qubit losses with QND measurement qubit+7qutrit system")
parser.add_argument('--phi_tilde',  type=float, default=0.01, help = "Rotation angle")
parser.add_argument('--epsilon_choi',  type=float, default=0.0, help = "epsilon_choi")
parser.add_argument('--logical_state',  type=int, default=0, help = "logical state integer corresponding to: 0, 1, +, -, +i, -i")
parser.add_argument('--chi_threshold',  type=float, default=1e-3, help = "threshold for discarding Kraus operators in the chi matrix")
args = parser.parse_args()

phi_tilde = args.phi_tilde
epsilon_choi = args.epsilon_choi
jLog = args.logical_state
chi_threshold = args.chi_threshold


from random import randint
seme = randint(0,100)
#seme = 44
np.random.seed(seme)

print("seed:", seme)

choi_ideal = np.loadtxt("choiFinal_ideal.dat")

choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')

import os
folder_name = f'chi_{chi_threshold:.01e}_eps_{epsilon_choi:1.3f}'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

choi = np.real((1 - epsilon_choi) * choi_ideal + 6 * epsilon_choi * choi_experiment)

T_matrix = give_transformation_matrix()
chi_matrix = get_chi_from_choi(choi, T_matrix)
np.set_printoptions(precision=3) 


final_p_loss = []
index_confs = 0
rotation_ops = Rloss_all(phi_tilde * np.pi)
result_correction = []

LogicalStates = [ZeroL, OneL, (ZeroL + OneL)/np.sqrt(2), (ZeroL - OneL)/np.sqrt(2), (ZeroL + 1j* OneL)/np.sqrt(2), (ZeroL - 1j*OneL)/np.sqrt(2)]

LogicalStates_str = ["0", "1", "+", "-", "+i", "-i"]


file_data_name = os.path.join(folder_name, final_data_name + f"_state_{LogicalStates_str[jLog]}_phi_{phi_tilde}_eps_{epsilon_choi}.dat")    


print(f"logical state |{LogicalStates_str[jLog]}_L>")

#for outcomes_ancilla in  binary_raw_configurations(n=7).configurations:

for num_loss, loss_confs in binary_configurations().configurations.items():  
     final_prob = []
     result_correction = []
     num_losses = []

     for outcomes_ancilla in loss_confs:
        index_confs += 1
         
        num_loss = sum(outcomes_ancilla)


        projectors_ancilla = 1 - 2*np.array(outcomes_ancilla)

        phi = phi_tilde * np.pi 

        loss_pattern = []

        prob_total_event = 1.0
        prob_correction_logical_state = []
        psiL = LogicalStates[jLog]

    

        list_qubits = list(range(L))

    #    print("XL", qu.expect(XL, rho_L))
    #    print("ZL", qu.expect(ZL, rho_L))
    #    print("YL", qu.expect(1j * XL * ZL, rho_L))   
        null_state = False    
        rho_L = psiL * psiL.dag()

        for data_q in list_qubits:
            #apply Rloss with an angle phi
            rho_L = rotation_ops[data_q] * rho_L * rotation_ops[data_q].dag()
            #apply the QND detection unit
            rho_L = apply_qnd_process_unit(chi_matrix, rho_L, data_q, chi_threshold)

            if projectors_ancilla[data_q] == +1:
                prob_outcome = (rho_L * Pp_ancilla).tr()
                if abs(prob_outcome.imag) > 1e-5: print("warning: im prob_outcome = {prob_outcome}")
                if prob_outcome == 0:
                    #the state cannot be projected in the +1 eigenstate of the ancilla
                    null_state = True
                else:    
                    rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / abs(prob_outcome)
            elif  projectors_ancilla[data_q] == -1:
                prob_outcome = (rho_L * Pm_ancilla).tr()
                if abs(prob_outcome.imag) > 1e-5: print("warning: im prob_outcome = {prob_outcome}")
                if prob_outcome == 0:
                    null_state = True
                else:               
                    rho_L = Pm_ancilla * rho_L * Pm_ancilla.dag() / abs(prob_outcome)
                    rho_L = Xa  * rho_L * Xa.dag() #reinitializing ancilla
            loss_pattern.append(outcomes_ancilla[data_q])
            prob_total_event *= prob_outcome
            print(data_q, projectors_ancilla[data_q], outcomes_ancilla[data_q], prob_outcome)


        losses = np.where(loss_pattern)[0].tolist()
        kept_qubits = list(set(range(L)) - set(losses))

        if sum(outcomes_ancilla) >= 7 or null_state:
            print(prob_total_event)
            correction_successful = 0.0
            prob_correction_logical_state.append(correction_successful)
        else:
            w_0 = rho_L.ptrace(kept_qubits)

            rho_L = qu.tensor([qu.maximally_mixed_dm(3)] * len(losses) + [w_0] + [qu.fock_dm(2,0)])

            permutation_order_q = {}
            for j, el in enumerate(losses + kept_qubits):
                permutation_order_q[el] = j
    #        print("permutation_order_q", permutation_order_q)

            stab_qubits_new_order = []
            for stab in stab_qubits:
                stab_qubits_new_order.append([permutation_order_q[q] for q in stab])
    #            print(stab, [permutation_order_q[q] for q in stab])

            Sx = [X[j1] * X[j2] * X[j3] * X[j4] for j1,j2,j3,j4 in stab_qubits_new_order]
            Sz = [Z[j1] * Z[j2] * Z[j3] * Z[j4] for j1,j2,j3,j4 in stab_qubits_new_order]

            Px = [(Id + el) / 2 for el in Sx]
            Pz = [(Id + el) / 2 for el in Sz]

            Pmx = [(Id - el) / 2 for el in Sx]
            Pmz = [(Id - el) / 2 for el in Sz]

            state_after_measure = rho_L

            stabZ_eigenvalues = []
            for meas in range(3):
                #TO DO ADD a check on the imaginary part of prob_plus that should be VERY small
                prob_plus =  (Pz[meas] * state_after_measure).tr()
                if abs(prob_plus.imag) > 1e-5: print("warning: im prob_plus = {prob_plus}")
                prob_plus = abs(prob_plus)
    
                if prob_plus > 1: prob_plus = 1
                if prob_plus < 0: prob_plus = 0    

                result = 2 * np.random.binomial(1, prob_plus) - 1
                stabZ_eigenvalues.append(result)
                print(f"RESULT Z {meas} {prob_plus} {result}")
                if result == + 1:
                    state_after_measure =  Pz[meas] * state_after_measure * Pz[meas] / prob_plus
                else:
                    state_after_measure =  Pmz[meas] * state_after_measure * Pmz[meas] / (1 - prob_plus)

            stabX_eigenvalues = []
            for meas in range(3):
                prob_plus =  (Px[meas] * state_after_measure).tr()
                if abs(prob_plus.imag) > 1e-5: print("warning: im prob_plus = {prob_plus}")
                prob_plus = abs(prob_plus)
    
                if prob_plus > 1: prob_plus = 1
                if prob_plus < 0: prob_plus = 0    
                result = 2 * np.random.binomial(1, prob_plus) - 1
                stabX_eigenvalues.append(result)
                print(f"RESULT X {meas} {prob_plus} {result}")
                if result == +1:
                    state_after_measure =  Px[meas] * state_after_measure * Px[meas] / prob_plus
                else:
                    state_after_measure = Pmx[meas] * state_after_measure * Pmx[meas] / (1 - prob_plus)


            print("XL", qu.expect(XL, state_after_measure))
            print("ZL", qu.expect(ZL, state_after_measure))
            print("YL", qu.expect(1j * XL * ZL, state_after_measure))

            if jLog in (0,1):
                correction_successful = (1 +  abs(qu.expect(ZL, state_after_measure))) / 2
            elif jLog in (2,3):
                correction_successful = (1 +  abs(qu.expect(XL, state_after_measure))) / 2
            elif jLog in (4,5):
                correction_successful = (1 +  abs(qu.expect(1j * XL * ZL, state_after_measure))) / 2

        print("prob_correction_logical_state:", correction_successful)


        conf_loss = int("".join(str(_) for _ in outcomes_ancilla)) 

        final_p_loss.append([phi_tilde, conf_loss, correction_successful, num_loss, prob_total_event])
        np.savetxt(file_data_name, final_p_loss, fmt= '%1.3f\t' + '%07d\t' + '%.10e\t' +'%d\t' + '%1.10f\t')
        
