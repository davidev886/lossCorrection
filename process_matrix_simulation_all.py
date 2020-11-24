import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu

from utils.qnd_error_gen import pick_qnd_error
from utils.binary_conf import binary_configurations
import datetime
now = datetime.datetime.now()
final_data_name = now.strftime("%Y%m%d%H%M")

import argparse
#python process_matrix_simulation_all.py --phi_tilde  --epsilon_choi
parser = argparse.ArgumentParser(description = "Simulate qubit losses with QND measurement qubit+7qutrit system")
parser.add_argument('--phi_tilde',  type=float, default=0.0, help = "Rotation angle")
parser.add_argument('--epsilon_choi',  type=float, default=0.0, help = "Rotation angle")
args = parser.parse_args()

from utils.p_operators_qutrit import *

from random import randint
seme = randint(0,100)


np.random.seed(seme)

print("seed:", seme)

choi_ideal = np.loadtxt("choiFinal_ideal.dat")

choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')

epsilon_choi = args.epsilon_choi
choi = (1 - epsilon_choi) * choi_ideal + epsilon_choi * choi_experiment

T_matrix = give_transformation_matrix()
chi_matrix = get_chi_from_choi(choi, T_matrix)
    
phi_tilde = args.phi_tilde

final_p_loss = []
index_confs = 0
rotation_ops = Rloss_all(phi_tilde * np.pi)
for num_loss, loss_confs in binary_configurations().configurations.items():  

        result_correction = []
        num_losses = []
        for outcomes_ancilla in loss_confs:
            index_confs += 1
        
            projectors_ancilla = 1 - 2*np.array(outcomes_ancilla)

            phi = phi_tilde * np.pi 

            loss_pattern = []

            a_0 = np.random.random()  + np.random.random() * 1j
            b_0 = np.random.random()  + np.random.random() * 1j

            a = a_0 / np.sqrt(abs(a_0)**2 + abs(b_0)**2)
            b = b_0 / np.sqrt(abs(a_0)**2 + abs(b_0)**2)

            psiL = a * ZeroL + b * OneL
        
        
            print("outcomes_ancilla", outcomes_ancilla)
            prob_single_loss = []
            prob_total_event = 1.0
            
            
            for data_q in range(L):


                rho_L = psiL * psiL.dag()
                #apply Rloss with an angle phi
                rho_L = rotation_ops[data_q] * rho_L * rotation_ops[data_q].dag()
                rho_L = apply_qnd_process_unit(chi_matrix, rho_L, data_q, T_matrix)
 
                if projectors_ancilla[data_q] == +1:
                    prob_outcome = (rho_L * Pp_ancilla).tr()
                    if abs(prob_outcome.imag) > 1e-5: print("warning: im prob_outcome = {prob_outcome}")
                    prob_outcome = abs(prob_outcome)
                    rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / prob_outcome
                elif  projectors_ancilla[data_q] == -1:
                    prob_outcome = (rho_L * Pm_ancilla).tr()
                    if abs(prob_outcome.imag) > 1e-5: print("warning: im prob_outcome = {prob_outcome}")
                    prob_outcome = abs(prob_outcome)                    
                    rho_L = Pm_ancilla * rho_L * Pm_ancilla.dag() / prob_outcome                 

                loss_pattern.append(outcomes_ancilla[data_q])
                prob_single_loss.append(prob_outcome)
                prob_total_event *= prob_outcome
                print(index_confs, data_q, projectors_ancilla[data_q], outcomes_ancilla[data_q], prob_outcome)
            
            losses = np.where(loss_pattern)[0].tolist()
            kept_qubits = list(set(range(L)) - set(losses))
            print(kept_qubits)
            if sum(outcomes_ancilla) == 7:
                correction_successful = False
                print("correction_successful:", correction_successful)
                print("qu.expect(ZL, state_after_measure)", stabZ_eigenvalues, f"{qu.expect(ZL, psiL):1.4}", f"{qu.expect(ZL, state_after_measure):1.4}" )
                print("qu.expect(XL, state_after_measure)", stabX_eigenvalues,f"{qu.expect(XL, psiL):1.4}", f"{qu.expect(XL, state_after_measure):1.4}")
                conf_loss = int("".join(str(_) for _ in outcomes_ancilla)) 
                final_p_loss.append([phi_tilde, conf_loss, int(correction_successful), num_loss, prob_total_event] + prob_single_loss)
                result_correction.append([int(correction_successful), num_loss, prob_total_event])
                np.savetxt(final_data_name + f"_loss_process_matrix_exact.dat", final_p_loss, fmt= '%1.3f\t' + '%07d\t' + '%d\t' * 2 + '%1.10f\t' + '%1.10f\t' * len(prob_single_loss))
                break

            w_0 = rho_L.ptrace(kept_qubits)
            rho_L = qu.tensor([qu.fock_dm(3,0)] * len(losses) + [w_0] + [qu.fock_dm(2,0)])

            permutation_order_q = {}
            for j, el in enumerate(losses + kept_qubits):
                permutation_order_q[el] = j
    #        print("permutation_order_q", permutation_order_q)            
        
            stab_qubits_new_order = []
            for stab in stab_qubits:
                stab_qubits_new_order.append([permutation_order_q[q] for q in stab])
                print(stab, [permutation_order_q[q] for q in stab])

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
    #            print(f"RESULT Z {meas} {prob_plus} {result}")    
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
    #            print(f"RESULT X {meas} {prob_plus} {result}")      
                if result == +1:
                    state_after_measure =  Px[meas] * state_after_measure * Px[meas] / prob_plus
                else:
                    state_after_measure = Pmx[meas] * state_after_measure * Pmx[meas] / (1 - prob_plus)

            expected_Z = np.abs(qu.expect(ZL, psiL))
            expected_X = np.abs(qu.expect(XL, psiL))

            measured_Z = np.abs(qu.expect(ZL, state_after_measure))
            measured_X = np.abs(qu.expect(XL, state_after_measure))

            correction_successful = ((np.abs(expected_Z - measured_Z) < 1e-7)
                                        and 
                                        (np.abs(expected_X - measured_X) < 1e-7))

            print("correction_successful:", correction_successful)
            print("qu.expect(ZL, state_after_measure)", stabZ_eigenvalues, f"{qu.expect(ZL, psiL):1.4}", f"{qu.expect(ZL, state_after_measure):1.4}" )
            print("qu.expect(XL, state_after_measure)", stabX_eigenvalues,f"{qu.expect(XL, psiL):1.4}", f"{qu.expect(XL, state_after_measure):1.4}")
            conf_loss = int("".join(str(_) for _ in outcomes_ancilla)) 
            final_p_loss.append([phi_tilde, conf_loss, int(correction_successful), num_loss, prob_total_event] + prob_single_loss)
            result_correction.append([int(correction_successful), num_loss, prob_total_event])
            np.savetxt(final_data_name + f"_loss_procmat_exact_phi_{phi_tilde}_eps_{epsilon_choi}.dat", final_p_loss, fmt= '%1.3f\t' + '%07d\t' + '%d\t' * 2 + '%1.10f\t' + '%1.10f\t' * len(prob_single_loss))

        result_correction = np.array(result_correction)


np.savetxt(final_data_name + f"_loss_procmat_exact_phi_{phi_tilde}_eps_{epsilon_choi}.dat", final_p_loss, fmt= '%1.3f\t' + '%07d\t' + '%d\t' * 2 + '%1.10f\t' + '%1.10f\t' * len(prob_single_loss))
