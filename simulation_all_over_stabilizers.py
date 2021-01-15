import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu
np.set_printoptions(precision=4, suppress=True)

from utils.qnd_error_gen import pick_qnd_error
from utils.binary_conf import binary_configurations, binary_raw_configurations
from utils.p_operators_qutrit import *

import datetime
now = datetime.datetime.now()
final_data_name = now.strftime("%Y%m%d%H%M")

from sys import getsizeof

import time



import argparse
#python process_matrix_simulation_all.py --phi_tilde  --epsilon_choi
parser = argparse.ArgumentParser(description = "Simulate qubit losses with QND measurement qubit+7qutrit system")
parser.add_argument('--phi_tilde',  type=float, default=0.0, help = "Rotation angle")
parser.add_argument('--epsilon_choi',  type=float, default=0.0, help = "epsilon_choi")
parser.add_argument('--logical_state',  type=int, default=0, help = "logical state integer corresponding to: 0, 1, +, -, +i, -i")
parser.add_argument('--chi_threshold',  type=float, default=0.0, help = "threshold for discarding Kraus operators in the chi matrix")
args = parser.parse_args()

phi_tilde = args.phi_tilde
epsilon_choi = args.epsilon_choi
jLog = args.logical_state
chi_threshold = args.chi_threshold


from random import randint
seme = randint(0,100)
#seme = 44
np.random.seed(seme)

choi_ideal = np.loadtxt("choiFinal_ideal.dat")

choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')

import os
folder_name = f'chi_{chi_threshold:.01e}_eps_{epsilon_choi:1.3f}_over_stab'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

choi = np.real((1 - epsilon_choi) * choi_ideal + 6 * epsilon_choi * choi_experiment)

T_matrix = give_transformation_matrix()
chi_matrix = get_chi_from_choi(choi, T_matrix)
chi_matrix[np.abs(chi_matrix) < chi_threshold] = 0

final_p_loss = []
index_confs = 0
rotation_ops = Rloss_all(phi_tilde * np.pi)
result_correction = []

LogicalStates = [ZeroL, OneL, (ZeroL + OneL)/np.sqrt(2), (ZeroL - OneL)/np.sqrt(2), (ZeroL + 1j* OneL)/np.sqrt(2), (ZeroL - 1j*OneL)/np.sqrt(2)]

LogicalStates_str = ["0", "1", "+", "-", "+i", "-i"]


file_data_name = os.path.join(folder_name, final_data_name + f"_state_{LogicalStates_str[jLog]}_phi_{phi_tilde:1.2f}_eps_{epsilon_choi}.dat")    
file_data_name_meas_report = os.path.join(folder_name, "stab_" + final_data_name + f"_state_{LogicalStates_str[jLog]}_phi_{phi_tilde:1.2f}_eps_{epsilon_choi}.dat")

print(f"logical state |{LogicalStates_str[jLog]}_L>")

#for outcomes_ancilla in  binary_raw_configurations(n=7).configurations:

for num_loss, loss_confs in binary_configurations().configurations.items():  
     final_prob = []
     result_correction = []
     num_losses = []

     for outcomes_ancilla in loss_confs:
        print(outcomes_ancilla)
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
            rho_L.tidyup(atol = 1e-6)
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
            rho_L = qu.tensor([qu.fock_dm(3,0)] * len(losses) + [w_0] + [qu.fock_dm(2,0)])

            permutation_order_q = {}
            for j, el in enumerate(losses + kept_qubits):
                permutation_order_q[el] = j
    #        print("permutation_order_q", permutation_order_q)

            stab_qubits_new_order = []
            for stab in stab_qubits:
                stab_qubits_new_order.append([permutation_order_q[q] for q in stab])


            Sx = [X[j1] * X[j2] * X[j3] * X[j4] for j1,j2,j3,j4 in stab_qubits_new_order]
            Sz = [Z[j1] * Z[j2] * Z[j3] * Z[j4] for j1,j2,j3,j4 in stab_qubits_new_order]
            
            PPx = [[(Id + el) / 2, (Id - el) / 2] for el in Sx]
            PPz = [[(Id + el) / 2, (Id - el) / 2] for el in Sz]            
            
            
            average_value_each_stab_meas = []
            correction_each_measurement = []
            report_on_stab_measurements = []
            for meas_binary_X,meas_binary_Z in product(range(8), range(8)): 
                state_after_measure = qu.Qobj(rho_L[:], dims = rho_L.dims)
                configuration_str_X = bin(meas_binary_X)[2:].zfill(3)
                configuration_int_X = [int(_) for _ in configuration_str_X]
                configuration_str_Z = bin(meas_binary_Z)[2:].zfill(3)
                configuration_int_Z = [int(_) for _ in configuration_str_Z]                
                probability_each_measurement = []
                for stab_num, outcome_stab in enumerate(configuration_int_X):
                    prob = (PPx[stab_num][outcome_stab] * state_after_measure).tr() 
                    if np.abs(prob) > 0:
                        state_after_measure = PPx[stab_num][outcome_stab] * state_after_measure * PPx[stab_num][outcome_stab].dag() / prob
                        probability_each_measurement.append(prob)                        
                    else:
                        probability_each_measurement.append(0)

                for stab_num, outcome_stab in enumerate(configuration_int_Z):
                    prob = (PPz[stab_num][outcome_stab] * state_after_measure).tr() 
                    if np.abs(prob) > 0:
                        state_after_measure = PPz[stab_num][outcome_stab] * state_after_measure * PPz[stab_num][outcome_stab].dag() / prob
                        probability_each_measurement.append(prob)                        
                    else:
                        probability_each_measurement.append(0)
                print(configuration_int_X, configuration_int_Z, np.prod(probability_each_measurement), qu.expect(XL, state_after_measure), qu.expect(ZL, state_after_measure), qu.expect(1j * XL * ZL, state_after_measure))                        

                if jLog in (0,1):
                    correction_successful = (1 +  abs(qu.expect(ZL, state_after_measure))) / 2
                elif jLog in (2,3):
                    correction_successful = (1 +  abs(qu.expect(XL, state_after_measure))) / 2
                elif jLog in (4,5):
                    correction_successful = (1 +  abs(qu.expect(1j * XL * ZL, state_after_measure))) / 2
    
                average_value_each_stab_meas.append(np.prod(probability_each_measurement) * correction_successful)
                conf_loss = int("".join(str(_) for _ in outcomes_ancilla))
                conf_stab_meas = int("".join(str(_) for _ in configuration_int_X + configuration_int_Z))
                
                report_on_stab_measurements.append([phi_tilde, 
                                                    conf_loss, 
                                                    conf_stab_meas] + 
                                                    probability_each_measurement + 
                                                    [np.real(np.prod(probability_each_measurement)),
                                                    np.real(correction_successful),
                                                    np.real(np.prod(probability_each_measurement) *  correction_successful)
                                                    ])
            with open(file_data_name_meas_report, 'a') as file_rep_stab:
                np.savetxt(file_rep_stab, report_on_stab_measurements, fmt= '%1.3f\t' + '%07d\t'+ '%06d\t'  + '%1.8f\t' * (len(probability_each_measurement) + 3))
            print("prob_of_succ_correction", np.sum(average_value_each_stab_meas))
            conf_loss = int("".join(str(_) for _ in outcomes_ancilla)) 

            final_p_loss.append([phi_tilde, conf_loss, np.real(np.sum(average_value_each_stab_meas)), num_loss, np.real(prob_total_event)])
            np.savetxt(file_data_name, final_p_loss, fmt= '%1.3f\t' + '%07d\t' + '%.10e\t' +'%d\t' + '%1.12f\t')

np.savetxt(file_data_name, final_p_loss, fmt= '%1.3f\t' + '%07d\t' + '%.10e\t' +'%d\t' + '%1.12f\t')
