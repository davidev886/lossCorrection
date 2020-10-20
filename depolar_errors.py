import matplotlib.pyplot as plt
    
import numpy as np
#from qutip import *
import qutip as qu

from utils_qnd import pick_qnd_error

from p_operators import *

from random import randint
seme = randint(0,100)

seme=47
np.random.seed(seme)

print("\nseed:", seme)

num_trials = 1

final_p_loss = []

p_qnd = 0.0
p_loss = 0.3
for p_loss in np.arange(0.00,0.9,0.05):
    trial = 0
    result_correction = []  
    num_losses = []
    num_qnd_errors = []
    while trial < num_trials:
 #   if 1:



        a = np.random.random()  + np.random.random() * 1j
        b = np.random.random()  + np.random.random() * 1j

        psiL = (a * ZeroL + b * OneL).unit()

        #1 means lost, 0 means not lost

        random_loss = []
        qnd_errors = []
        qubit_detected = []        
        for qubit in range(L):
            qnd_err_str = pick_qnd_error(p_qnd)
            err1_str, err2_str = qnd_err_str
            qnd_error = (err2_str in ("X", "Y")) + 0        
            
            qnd_detection = (not (err2_str in ("X", "Y"))) + 0
            loss_bool = np.random.binomial(1, p_loss)
            random_loss.append(loss_bool)
            qnd_errors.append(qnd_error)
            qubit_detected.append(qnd_detection)


#         print(f"{'random_loss':30s}", random_loss)
#         print(f"{'qnd_errors':30s}", qnd_errors)  
#         print(f"{'qubit_detected':30s}", qubit_detected)               


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
                
#         print(f'{"correctable_events":30s}', correctable_events)
#         print(f'{"non_correctable_events":30s}', non_correctable_events)
#         print(f'{"to_check_no_loss":30s}', to_check_no_losses)        
#         print(f'{"to_check_no_qnd":30s}', to_check_no_qnd_errors)                
#         exit()
    for _ in range(5):
        if correctable_events:
            trial += 1
            correction_successful = True
            result_correction.append(correction_successful + 0)
            continue
        elif non_correctable_events:
            trial += 1
            correction_successful = False
            result_correction.append(correction_successful + 0)
            continue
        elif to_check_no_qnd_errors or to_check_no_losses:

            trial += 1
                    

#         if 1 < len(losses) < 7:
#             w_0 = psiL.ptrace(kept_qubits)
#             w = (qu.tensor([qu.fock_dm(2,0)] * len(losses) + [w_0])).unit()
#         elif len(losses) == 7:
#             w = (qu.tensor([qu.fock_dm(2,0)] * len(losses))).unit()        
#         else:
#             w = (qu.ket2dm(psiL)).unit()

        w = psiL.ptrace(kept_qubits)
         
        #1 means qnd error
        qnd_errors = np.random.binomial(1, p_qnd, 7)
    


        detected_state =[ (qnd_error + loss) % 2 for qnd_error,loss in zip(qnd_errors, random_loss)]
            
                    

        detected_losses = np.where(detected_state)[0].tolist()

        # detection error: no loss detected
        detected_losses = []

        guessed_kept_qubits = list(set(range(L)) - set(detected_losses))
        
        print(f"{'random_loss':25}", random_loss)
        print(f"{'qnd_errors':25}", qnd_errors)
        print(f"{'detected_state':25}", detected_state)
        print(f"{'detected_losses':25}", detected_losses)        
        print(f"{'guessed_kept_qubits':25}", guessed_kept_qubits)                


#         if (len(detected_losses) == 0):
#             continue
#         elif (len(detected_losses) > 4):
#             trial += 1
#             correction_successful = False
#             result_correction.append(correction_successful + 0)
#             continue
#         else:
#             trial += 1
#             print("trial", trial, p_loss, p_qnd)

 
        pemutation_order_q = {}
        print(kept_qubits)
        for j, el in enumerate(detected_losses + guessed_kept_qubits):
            pemutation_order_q[el] = j
        
        print("pemutation_order_q", pemutation_order_q)

        stab_qubits_new_order = []
        for stab in stab_qubits:
            stab_qubits_new_order.append([pemutation_order_q[q] for q in stab])
            print(stab, [pemutation_order_q[q] for q in stab])
            

        Sx = [X[j1] * X[j2] * X[j3] * X[j4] for j1,j2,j3,j4 in stab_qubits_new_order]
        Sz = [Z[j1] * Z[j2] * Z[j3] * Z[j4] for j1,j2,j3,j4 in stab_qubits_new_order]

        Id = qu.tensor([qu.qeye(2)] * (L - len(losses)))
        
        temp = [[qu.qeye(2)] * j + [qu.sigmax()] + [qu.qeye(2)] * (L  - len(losses) - j - 1) for j in range(L - len(losses))]
        X = [qu.tensor(temp[j]) for j in range(L - len(losses))]

        temp = [[qu.qeye(2)] * j + [qu.sigmay()] + [qu.qeye(2)] * (L - len(losses) - j - 1) for j in range(L - len(losses))]
        Y = [qu.tensor(temp[j]) for j in range(L - len(losses))]

        temp = [[qu.qeye(2)] * j + [qu.sigmaz()] + [qu.qeye(2)] * (L  - len(losses) - j - 1) for j in range(L - len(losses))]
        Z = [qu.tensor(temp[j]) for j in range(L - len(losses))]
        
        stab_qubits = [[0,1,2,3], [1,2,4], [2,3,5]]

        Sx = [None, None, None]
        Sz = [None, None, None]
                
        for j_el, st_el in enumerate(stab_qubits):
            if len(st_el) == 4:
                j1,j2,j3,j4 = st_el
                Sx[j_el] = X[j1] * X[j2] * X[j3] * X[j4]
                Sz[j_el] = Z[j1] * Z[j2] * Z[j3] * Z[j4]
            elif len(st_el) == 3:
                j1,j2,j3 = st_el
                Sx[j_el] = X[j1] * X[j2] * X[j3]
                Sz[j_el] = Z[j1] * Z[j2] * Z[j3]                
            elif len(st_el) == 2:
                j1,j2 = st_el
                Sx[j_el] = X[j1] * X[j2]
                Sz[j_el] = Z[j1] * Z[j2]                
            elif len(st_el) == 1:
                j1 = st_el
                Sx[j_el] = X[j1]
                Sz[j_el] = Z[j1]                
            elif len(st_el) == 0:
                Sx[j_el] = Id    
                Sz[j_el] = Id


        ZL = Z[0] * Z[1] * Z[2] * Z[3] * Z[4] * Z[5]

        XL = X[0] * X[1] * X[2] * X[3] * X[4] * X[5]
            
        Px = [(Id + el) / 2 for el in Sx]
        Pz = [(Id + el) / 2 for el in Sz]

        Pmx = [(Id - el) / 2 for el in Sx]
        Pmz = [(Id - el) / 2 for el in Sz]                

        stabZ_eigenvalues = []

        state_after_measure = w
        for meas in range(3):
            prob_plus =  (Pz[meas] * state_after_measure).tr()
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
            if prob_plus > 1: prob_plus = 1
            if prob_plus < 0: prob_plus = 0    
            result = 2 * np.random.binomial(1, prob_plus) - 1
            stabX_eigenvalues.append(result)
            print(f"RESULT X {meas} {prob_plus} {result}")      
            if result == +1:
                state_after_measure =  Px[meas] * state_after_measure * Px[meas] / prob_plus
            else:
                state_after_measure = Pmx[meas] * state_after_measure * Pmx[meas] / (1 - prob_plus)
########################################################################
                 
        correction_qubits = [0, 4, 6]
        
#         for j_stab, syndrome in enumerate(stabZ_eigenvalues):
#             if syndrome == -1:
#                 qubit_correction = pemutation_order_q[correction_qubits[j_stab]]
#                 #print("Z syndrome corrected by X[", correction_qubits[j_stab], "]")
#                 #print("on the new state: Z syn ecorrected by X[", qubit_correction, "]")
#                 op_correction = X[qubit_correction]
#                 state_after_measure = op_correction * state_after_measure * op_correction   



#         for j_stab, syndrome in enumerate(stabX_eigenvalues):
#             if syndrome == -1:
#                 qubit_correction = pemutation_order_q[correction_qubits[j_stab]]
#                 print("X syndrome corrected by Z[", correction_qubits[j_stab], "]")
#                 print("on the new state: X syn ecorrected by Z[", qubit_correction, "]")
#                 op_correction = Z[qubit_correction]
#                 state_after_measure = op_correction * state_after_measure * op_correction        

        print("before correction")                               
        print("Sz", [(state_after_measure * Sz[j]).tr() for j in range(3)])
        print("Sx", [(state_after_measure * Sx[j]).tr() for j in range(3)])    



        if 0:
            for corr in [6]:
                op_correction = X[pemutation_order_q[corr]]
                state_after_measure = op_correction * state_after_measure * op_correction        

            for corr in [1,3,5]:
                op_correction = Z[pemutation_order_q[corr]]
                state_after_measure = op_correction * state_after_measure * op_correction   
            
            print("after correction")                               
            print("Sz", [(state_after_measure * Sz[j]).tr() for j in range(3)])
            print("Sx", [(state_after_measure * Sx[j]).tr() for j in range(3)])    
    
                                                    
        state_after_measure = state_after_measure.unit()
        
        
        measured_Z = np.abs(qu.expect(ZL, state_after_measure))
        measured_X = np.abs(qu.expect(XL, state_after_measure))

        correction_successful = ((np.abs(expected_Z - measured_Z) < 1e-7)
                                    and 
                                    (np.abs(expected_X - measured_X) < 1e-7))

        print("correction_successful:", correction_successful)
        print("qu.expect(ZL, state_after_measure)", stabZ_eigenvalues, f"{expected_Z:1.4}", f"{qu.expect(ZL, state_after_measure):1.4}" )
        print("qu.expect(XL, state_after_measure)", stabX_eigenvalues,f"{expected_X:1.4}", f"{qu.expect(XL, state_after_measure):1.4}")

        result_correction.append(correction_successful + 0)

        done_trials = trial
        num_losses.append(len(losses))
        num_qnd_errors.append(sum(qnd_errors))        

    final_p_loss.append([p_loss, p_qnd, np.mean(result_correction), np.std(result_correction), np.mean(num_losses), np.mean(num_qnd_errors)])

exit()



np.savetxt(f"final_qnd_faulty_{num_trials}_pqnd_{p_qnd:1.3f}.dat", final_p_loss, fmt='%1.3f\t%1.3f\t' + '%1.6f\t' * 4)
final = np.array(final_p_loss)

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
plt.savefig(f"final_qnd_faulty_{num_trials}_pqnd_{p_qnd:1.3f}.pdf")

plt.show()

exit()