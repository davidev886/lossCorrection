import matplotlib.pyplot as plt
    
import numpy as np
#from qutip import *
import qutip as qu

from itertools import combinations

from utils.p_operators import *

from random import randint
seme = randint(0,100)

np.random.seed(seme)

print("\nseed:", seme)

num_trials = 100

final_p_error = []

p_qnd = 0.05

for p_error in np.arange(0.0,0.9,0.05):
    trial = 0
    result_correction = []  
    num_losses = []
    num_qnd_errors = []
    while trial < num_trials:

#       1 means lost, 0 means not lost
        random_loss = np.random.binomial(1, p_error, 7)
#       1 means qnd error
        qnd_errors = np.random.binomial(1, p_qnd, 7)

        A_loss = (sum(random_loss) == 0)
        B_loss = (sum(random_loss) == 1)
        C_loss = (2 <= sum(random_loss) <= 4)
        D_loss = (5 <= sum(random_loss) <= 7)


        A_qnd = (sum(qnd_errors) == 0)
        B_qnd = (sum(qnd_errors) == 1)
        C_qnd = (2 <= sum(qnd_errors) <= 4)
        D_qnd = (5 <= sum(qnd_errors) <= 7)


        correctable_events = (A_loss and A_qnd) or (A_loss and B_qnd) or (B_loss and A_qnd)
        non_correctable_events = ( (A_loss and D_qnd)
                                or (B_loss and (B_qnd or  C_qnd or D_qnd))
                                or (C_loss and (B_qnd or  C_qnd or D_qnd))
                                or D_loss )
                                    
        to_check_no_losses     = A_loss and C_qnd
        to_check_no_qnd_errors = C_loss and A_qnd
        
        
        
        if correctable_events:
            trial += 1
            correction_successful = True
            result_correction.append(correction_successful + 0)
            num_losses.append(sum(random_loss))
            num_qnd_errors.append(sum(qnd_errors))
            continue
        elif non_correctable_events:
            trial += 1
            correction_successful = False
            result_correction.append(correction_successful + 0)
            num_losses.append(sum(random_loss))
            num_qnd_errors.append(sum(qnd_errors))            
            continue
        elif to_check_no_losses or to_check_no_qnd_errors:
            trial += 1
            if to_check_no_losses:
                #transform the wrongly detected losses in actual losses and make the correction
                losses = np.where(qnd_errors)[0].tolist()
            if to_check_no_qnd_errors:                
                losses = np.where(random_loss)[0].tolist()                

            kept_qubits = list(set(range(L)) - set(losses))                
                
            a = np.random.random()  + np.random.random() * 1j
            b = np.random.random()  + np.random.random() * 1j

            psiL = (a * ZeroL + b * OneL).unit()

            w_0 = psiL.ptrace(kept_qubits)
            w = (qu.tensor([qu.fock_dm(2,0)] * len(losses) + [w_0])).unit()

            permutation_order_q = {}
            print(kept_qubits)
            for j, el in enumerate(losses + kept_qubits):
                permutation_order_q[el] = j
        
#            print("permutation_order_q", permutation_order_q)

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

            stabZ_eigenvalues = []

            state_after_measure = w
            for meas in range(3):
                prob_plus =  (Pz[meas] * state_after_measure).tr()
                print(prob_plus)
                if prob_plus > 1: prob_plus = 1
                if prob_plus < 0: prob_plus = 0    

                result = 2 * np.random.binomial(1, prob_plus) - 1
                stabZ_eigenvalues.append(result)
            #    print(f"RESULT Z {meas} {prob_plus} {result}")    
                if result == + 1:
                    state_after_measure =  Pz[meas] * state_after_measure * Pz[meas] / prob_plus
                else:
                    state_after_measure =  Pmz[meas] * state_after_measure * Pmz[meas] / (1 - prob_plus)

            stabX_eigenvalues = []
            for meas in range(3):
                prob_plus =  (Px[meas] * state_after_measure).tr()
                if prob_plus > 1: prob_plus = 1
                if prob_plus < 0: prob_plus = 0    
                result = 2 * np.random.binomial(1, prob_plus) - 1
                stabX_eigenvalues.append(result)
            #    print(f"RESULT X {meas} {prob_plus} {result}")      
                if result == +1:
                    state_after_measure =  Px[meas] * state_after_measure * Px[meas] / prob_plus
                else:
                    state_after_measure = Pmx[meas] * state_after_measure * Pmx[meas] / (1 - prob_plus)
    ########################################################################
                 
            correction_qubits = [0, 4, 6]
        
    #         for j_stab, syndrome in enumerate(stabZ_eigenvalues):
    #             if syndrome == -1:
    #                 qubit_correction = permutation_order_q[correction_qubits[j_stab]]
    #                 #print("Z syndrome corrected by X[", correction_qubits[j_stab], "]")
    #                 #print("on the new state: Z syn ecorrected by X[", qubit_correction, "]")
    #                 op_correction = X[qubit_correction]
    #                 state_after_measure = op_correction * state_after_measure * op_correction   



    #         for j_stab, syndrome in enumerate(stabX_eigenvalues):
    #             if syndrome == -1:
    #                 qubit_correction = permutation_order_q[correction_qubits[j_stab]]
    #                 print("X syndrome corrected by Z[", correction_qubits[j_stab], "]")
    #                 print("on the new state: X syn ecorrected by Z[", qubit_correction, "]")
    #                 op_correction = Z[qubit_correction]
    #                 state_after_measure = op_correction * state_after_measure * op_correction        

#            print("before correction")                               
#            print("Sz", [(state_after_measure * Sz[j]).tr() for j in range(3)])
#            print("Sx", [(state_after_measure * Sx[j]).tr() for j in range(3)])    



            if 0:
                for corr in [6]:
                    op_correction = X[permutation_order_q[corr]]
                    state_after_measure = op_correction * state_after_measure * op_correction        

                for corr in [1,3,5]:
                    op_correction = Z[permutation_order_q[corr]]
                    state_after_measure = op_correction * state_after_measure * op_correction   
            
                print("after correction")                               
                print("Sz", [(state_after_measure * Sz[j]).tr() for j in range(3)])
                print("Sx", [(state_after_measure * Sx[j]).tr() for j in range(3)])    
    
                                                    
            state_after_measure = state_after_measure.unit()
        
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

            result_correction.append(correction_successful + 0)

            done_trials = trial
            num_losses.append(len(losses))
            num_qnd_errors.append(sum(qnd_errors))        

    final_p_error.append([p_error, p_qnd, np.mean(result_correction), np.std(result_correction), np.mean(num_losses), np.mean(num_qnd_errors)])

np.savetxt(f"final_qnd_faulty_{num_trials}_pqnd_{p_qnd:1.3f}.dat", final_p_error, fmt='%1.3f\t%1.3f\t' + '%1.6f\t' * 4)
final = np.array(final_p_error)

import matplotlib.pyplot as plt

x_data = final[:,0]
y_data = final[:,2]
y_error = final[:,3] / np.sqrt(num_trials)
plt.errorbar(x_data, y_data,  yerr=y_error, fmt='o-')


x_data = np.linspace(0,1,100)

y_data = 1 - 7*x_data**3 + 21*x_data**5 - 21*x_data**6 + 6*x_data**7
plt.plot(x_data, y_data, '-')

y_data = 1 - x_data
plt.plot(x_data, y_data, '-')


plt.xlabel("p")
plt.ylabel("p(success)")
plt.savefig(f"final_qnd_faulty_{num_trials}_pqnd_{p_qnd:1.3f}_seed_{seme}.pdf")

plt.show()

exit()