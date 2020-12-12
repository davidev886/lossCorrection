import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu
from itertools import product
from utils.qnd_error_gen import pick_qnd_error
from utils.binary_conf import binary_configurations, binary_raw_configurations
from utils.p_operators_qutrit import *
import datetime
now = datetime.datetime.now()
final_data_name = now.strftime("%Y%m%d%H%M")

#python process_matrix_stab_errors_all.py  --phi_tilde 0.05 --epsilon_choi 0.0 --logical_state 2 --chi_threshold 1e-3 


import argparse

parser = argparse.ArgumentParser(description = "Simulate qubit losses with QND measurement qubit+7qutrit system")
parser.add_argument('--phi_tilde',  type=float, default=0.0, help = "Rotation angle")
parser.add_argument('--epsilon_choi',  type=float, default=0.0, help = "Rotation angle")
parser.add_argument('--logical_state',  type=int, default=0, help = "logical state integer corresponding to: 0, 1, +, -, +i, -i")
parser.add_argument('--chi_threshold',  type=float, default=1e-3, help = "threshold for discarding Kraus operators in the chi matrix")
args = parser.parse_args()

phi_tilde = args.phi_tilde
epsilon_choi = args.epsilon_choi
jLog = args.logical_state
chi_threshold = args.chi_threshold

p_err_stab = 0.04
            
LogicalStates = [ZeroL, OneL, (ZeroL + OneL)/np.sqrt(2), (ZeroL - OneL)/np.sqrt(2), (ZeroL + 1j* OneL)/np.sqrt(2), (ZeroL - 1j*OneL)/np.sqrt(2)]

LogicalStates_str = ["0", "1", "+", "-", "+i", "-i"]
            


from random import randint
seme = randint(0,100)

np.random.seed(seme)

print("seed:", seme)

choi_ideal = np.loadtxt("choiFinal_ideal.dat")

choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')



import os
folder_name = f'chi_{chi_threshold:.01e}_eps_{epsilon_choi:1.3f}_stab_error_{p_err_stab:1.3f}'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

file_data_name = os.path.join(folder_name, final_data_name + f"_stab_{p_err_stab}_state_{LogicalStates_str[jLog]}_phi_{phi_tilde}_eps_{epsilon_choi}.dat")    

print(f"logical state |{LogicalStates_str[jLog]}_L>")

choi = np.abs((1 - epsilon_choi) * choi_ideal + 6 * epsilon_choi * choi_experiment)

T_matrix = give_transformation_matrix()
chi_matrix = get_chi_from_choi(choi, T_matrix)
    

stab_errors_list = binary_raw_configurations(n=3).configurations


final_p_loss = []
index_confs = 0
rotation_ops = Rloss_all(phi_tilde * np.pi)





for num_loss, loss_confs in binary_configurations().configurations.items():  

    for outcomes_ancilla in loss_confs:
        index_confs += 1
    
        projectors_ancilla = 1 - 2*np.array(outcomes_ancilla)

        phi = phi_tilde * np.pi 

        loss_pattern = []
    
        print("outcomes_ancilla", outcomes_ancilla)
        prob_single_loss = []
        prob_total_event = 1.0

        prob_correction_logical_state = []
        psiL = LogicalStates[jLog]

    

        list_qubits = list(range(L))

    #    print("XL", qu.expect(XL, rho_L))
    #    print("ZL", qu.expect(ZL, rho_L))
    #    print("YL", qu.expect(1j * XL * ZL, rho_L))   
        null_state = False    
        rho_L = psiL * psiL.dag()
        print("outcomes_ancilla", num_loss, outcomes_ancilla)                            
        for data_q in range(L):

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
            prob_single_loss.append(prob_outcome)
            prob_total_event *= prob_outcome
            print(index_confs, data_q, projectors_ancilla[data_q], outcomes_ancilla[data_q], prob_outcome)
        
        losses = np.where(loss_pattern)[0].tolist()
        kept_qubits = list(set(range(L)) - set(losses))
        print(kept_qubits)
        if sum(outcomes_ancilla) >= 5 or null_state:
            # loop over the 2^(3+3) stabilizers errors
            for stab_errors_binary0_x, stab_errors_binary0_z in product(stab_errors_list, stab_errors_list):
                stabX_errors = np.array(stab_errors_binary0_x)
                stabZ_errors = np.array(stab_errors_binary0_z)

                prob_stab_event_X = np.prod(p_err_stab**stabX_errors * (1 - p_err_stab)**(1 - stabX_errors ))
                prob_stab_event_Z = np.prod(p_err_stab**stabZ_errors * (1 - p_err_stab)**(1 - stabZ_errors ))
                prob_stab_event =   prob_stab_event_X * prob_stab_event_Z   
                        
                correction_successful = 0.0

                print("qu.expect(ZL, state_after_measure)", f"{qu.expect(ZL, psiL):1.4}", f"{qu.expect(ZL, state_after_measure):1.4}" )
                print("qu.expect(XL, state_after_measure)", f"{qu.expect(XL, psiL):1.4}", f"{qu.expect(XL, state_after_measure):1.4}")
                conf_loss = int("".join(str(_) for _ in outcomes_ancilla))
                conf_stab_error = int("".join(str(_) for _ in stab_errors_binary0_x + stab_errors_binary0_z))
                print("correction_successful:", correction_successful, f"{conf_loss:07d}", f"{conf_stab_error:06d}")
                final_p_loss.append([phi_tilde, conf_loss, conf_stab_error, correction_successful, num_loss, prob_total_event, prob_stab_event] + prob_single_loss)
                np.savetxt(file_data_name, final_p_loss, fmt= '%1.3f\t' + '%07d\t' + '%06d\t' + '%.10e\t' +'%d\t' + '%1.10f\t' * 2 + '%1.10f\t' * len(prob_single_loss))
        
        elif sum(outcomes_ancilla) <= 4:

            w_0 = rho_L.ptrace(kept_qubits)#.unit()
            print("w_0.tr()", w_0.tr())
            rho_L = qu.tensor([qu.fock_dm(3,0)] * len(losses) + [w_0] + [qu.fock_dm(2,0)])

            permutation_order_q = {}
            for j, el in enumerate(losses + kept_qubits):
                permutation_order_q[el] = j

    
            stab_qubits_new_order = []
            for stab in stab_qubits:
                stab_qubits_new_order.append([permutation_order_q[q] for q in stab])
    #                print(stab, [permutation_order_q[q] for q in stab])

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
                
            # loop over the 2^(3+3) stabilizers errors
            for stab_errors_binary0_x, stab_errors_binary0_z in product(stab_errors_list, stab_errors_list):
                stabX_errors = np.array(stab_errors_binary0_x)
                stabZ_errors = np.array(stab_errors_binary0_z)
                #check if a loss happens on a faulty stabilizer
                #if so, the state is not correctable (as measuring a stabilizer on a lost qubit is impossible)
                faulty_stabX_qubits = [el for j,el in enumerate(stab_qubits) if stabX_errors[j]]
                faulty_stabZ_qubits = [el for j,el in enumerate(stab_qubits) if stabZ_errors[j]]
                print("stab_errors", stabX_errors, stabZ_errors)
                print(faulty_stabX_qubits, faulty_stabZ_qubits)

                loss_on_faulty_stabX = any([any([(loss in stab) for loss in losses]) for stab in faulty_stabX_qubits])
                loss_on_faulty_stabZ = any([any([(loss in stab) for loss in losses]) for stab in faulty_stabZ_qubits])
                prob_stab_event_X = np.prod(p_err_stab**stabX_errors * (1 - p_err_stab)**(1 - stabX_errors ))
                prob_stab_event_Z = np.prod(p_err_stab**stabZ_errors * (1 - p_err_stab)**(1 - stabZ_errors ))
                prob_stab_event =   prob_stab_event_X * prob_stab_event_Z           
                print("loss_on_faulty_stabX =", loss_on_faulty_stabX, "\tloss_on_faulty_stabZ =", loss_on_faulty_stabZ)

                if loss_on_faulty_stabX or loss_on_faulty_stabZ: 
                    correction_successful = 0.0 
                else: #now check is the logical operators have the same expectation values of the original state

                    if jLog in (0,1):
                        correction_successful = (1 +  abs(qu.expect(ZL, state_after_measure))) / 2
                    elif jLog in (2,3):
                        correction_successful = (1 +  abs(qu.expect(XL, state_after_measure))) / 2
                    elif jLog in (4,5):
                        correction_successful = (1 +  abs(qu.expect(1j * XL * ZL, state_after_measure))) / 2

#               print("qu.expect(ZL, state_after_measure)", f"{qu.expect(ZL, psiL):1.4}", f"{qu.expect(ZL, state_after_measure):1.4}" )
#               print("qu.expect(XL, state_after_measure)", f"{qu.expect(XL, psiL):1.4}", f"{qu.expect(XL, state_after_measure):1.4}")
            
                conf_loss = int("".join(str(_) for _ in outcomes_ancilla)) 
                conf_stab_error = int("".join(str(_) for _ in stab_errors_binary0_x + stab_errors_binary0_z))
                print("correction_successful:", correction_successful, f"{conf_loss:07d}", f"{conf_stab_error:06d}")
                final_p_loss.append([phi_tilde, conf_loss, conf_stab_error, correction_successful, num_loss, prob_total_event, prob_stab_event] + prob_single_loss)
            np.savetxt(file_data_name, final_p_loss, fmt= '%1.3f\t' + '%07d\t' + '%06d\t' + '%.10e\t' +'%d\t' + '%1.10f\t' * 2 + '%1.10f\t' * len(prob_single_loss))
    
np.savetxt(file_data_name, final_p_loss, fmt= '%1.3f\t' + '%07d\t' + '%06d\t' + '%.10e\t' +'%d\t' + '%1.10f\t' * 2 + '%1.10f\t' * len(prob_single_loss))
