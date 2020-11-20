import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu

from utils.qnd_error_gen import pick_qnd_error
from utils.binary_conf import binary_configurations
import datetime
now = datetime.datetime.now()
final_data_name = now.strftime("%Y%m%d%H%M")


from utils.p_operators_qutrit import *

from random import randint
seme = randint(0,100)

seme=35
np.random.seed(seme)

print("\nseed:", seme)

choi_ideal = np.loadtxt("choiFinal_ideal.dat")

choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')

epsilon_choi = 0.0
choi = (1 - epsilon_choi) * choi_ideal + epsilon_choi * choi_experiment

final_p_loss = []


index_confs = 0
for num_loss, loss_confs in binary_configurations().configurations.items():  
    result_correction = []
    num_losses = []

    for outcomes_ancilla in loss_confs:
        index_confs += 1

        projectors_ancilla = 1 - 2*np.array(outcomes_ancilla)

#        phi = phi_tilde * np.pi 

        loss_pattern = []

        a_0 = np.random.random()  + np.random.random() * 1j
        b_0 = np.random.random()  + np.random.random() * 1j

        a = a_0 / np.sqrt(abs(a_0)**2 + abs(b_0)**2)
        b = b_0 / np.sqrt(abs(a_0)**2 + abs(b_0)**2)

        psiL = a * ZeroL + b * OneL
        print("outcomes_ancilla", outcomes_ancilla)
        for data_q in range(L):

            #apply Rloss with an angle phi
#            rho_L = Rloss(psiL, phi, data_q)

            rho_L = apply_qnd_process_unit(choi, psiL, qu_data = data_q)
    
            if projectors_ancilla[data_q] == +1:
                prob_outcome = (rho_L * Pp_ancilla).tr()
                rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / prob_outcome
            elif  projectors_ancilla[data_q] == -1:
                prob_outcome = (rho_L * Pp_ancilla).tr()            
                rho_L = Pm_ancilla * rho_L * Pm_ancilla.dag() / prob_outcome                 

            loss_pattern.append(outcomes_ancilla[data_q])
            print(index_confs, f"prob_outcome={prob_outcome:1.3f}", f"QND_outcome={outcomes_ancilla[data_q]}")
            

        losses = np.where(loss_pattern)[0].tolist()
        kept_qubits = list(set(range(L)) - set(losses))
        
        if sum(losses) == 7:
            correction_successful = False
            result_correction.append(correction_successful + 0)
            continue 

        w_0 = rho_L.ptrace(kept_qubits)
        rho_L = qu.tensor([qu.fock_dm(3,0)] * len(losses) + [w_0] + [qu.fock_dm(2,0)])

        permutation_order_q = {}
        for j, el in enumerate(losses + kept_qubits):
            permutation_order_q[el] = j
#        print("permutation_order_q", permutation_order_q)            
        
        stab_qubits_new_order = []
        for stab in stab_qubits:
            stab_qubits_new_order.append([permutation_order_q[q] for q in stab])
            #print(stab, [permutation_order_q[q] for q in stab])

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
            prob_plus =  abs((Pz[meas] * state_after_measure).tr())
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
            prob_plus =  abs((Px[meas] * state_after_measure).tr())
            if prob_plus > 1: prob_plus = 1
            if prob_plus < 0: prob_plus = 0    
            result = 2 * np.random.binomial(1, prob_plus) - 1
            stabX_eigenvalues.append(result)
            print(f"RESULT X {meas} {prob_plus} {result}")      
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

        final_p_loss.append([int(bin(index_confs - 1)[2:]), int(correction_successful), num_loss])
        print(final_p_loss)
        np.savetxt(final_data_name + f"_loss_process_matrix_exact.dat", final_p_loss, fmt= '%07d\t' + '%d\t' * 2)

np.savetxt(final_data_name + f"_loss_process_matrix_exact.dat", final_p_loss, fmt= '%07d\t' + '%d\t' * 2)
