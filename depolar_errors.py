import matplotlib.pyplot as plt
    
import numpy as np
#from qutip import *
import qutip as qu

from utils.utils_qnd import pick_qnd_error

from utils.p_operators import *

from random import randint
seme = randint(0,100)

seme=5
np.random.seed(seme)

print("\nseed:", seme)

num_trials = 1

final_p_loss = []

p_qnd = 0.25
p_loss = 0.5
#for p_loss in np.arange(0.00,0.9,0.05):
if 1:
    trial = 0
    result_correction = []  
    num_losses = []
    num_qnd_errors = []
    errors_on_qubits = []
    while trial < num_trials:
        random_loss = []
        qnd_errors = []
        for qubit in range(L):
            qnd_err_str = pick_qnd_error(p_qnd)
            print(qnd_err_str)
            err1_str, err2_str = qnd_err_str
            errors_on_qubits.append(err1_str)
            qnd_error = (err2_str in ("X", "Y")) + 0
            loss_bool = np.random.binomial(1, p_loss)
            random_loss.append(loss_bool)
            qnd_errors.append(qnd_error)

        random_loss=                    [1, 0, 0, 0, 0, 1, 0]
        qnd_errors=                     [0, 0, 0, 0, 0, 0, 0]
        errors_on_qubits=               ['X', 'I', 'I', 'X', 'I', 'I', 'I']
        
        print(f"{'random_loss':30s}", random_loss)
        print(f"{'qnd_errors':30s}", qnd_errors)
        print(f"{'errors_on_qubits':30s}", errors_on_qubits)

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
                
        print(f'{"correctable_events":30s}', correctable_events)
        print(f'{"non_correctable_events":30s}', non_correctable_events)
        print(f'{"to_check_no_loss":30s}', to_check_no_losses)        
        print(f'{"to_check_no_qnd":30s}', to_check_no_qnd_errors)                


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
        elif to_check_no_losses or to_check_no_qnd_errors:
            trial += 1
            if to_check_no_losses:
                #transform the wrongly detected losses in actual losses and make the correction
                losses = np.where(qnd_errors)[0].tolist()
            if to_check_no_qnd_errors:                
                losses = np.where(random_loss)[0].tolist()
                # remove from the list of errors, the ones that occur on lost qubits:
                # this becasue lost qubits are replaced by fresh new qubits and this process 
                # happens with no errors                 
                for loss in losses:
                    errors_on_qubits[loss] = 'I'


                
                

            kept_qubits = list(set(range(L)) - set(losses))                
                
            a = np.random.random()  + np.random.random() * 1j
            b = np.random.random()  + np.random.random() * 1j

            psiL = (a * ZeroL + b * OneL).unit()

            w_0 = psiL.ptrace(kept_qubits)
            w = (qu.tensor([qu.fock_dm(2,0)] * len(losses) + [w_0])).unit()

            permutation_order_q = {}
            print(losses , kept_qubits)
            for j, el in enumerate(losses + kept_qubits):
                permutation_order_q[el] = j
            print("permutation_order_q", permutation_order_q)            
            #Apply to the state w the errors in the list errors_on_qubits
            errors_p = "XYZ"            
            for qubit_n, err_str in enumerate(errors_on_qubits):
                if err_str != 'I':
                    data_qubit = permutation_order_q[qubit_n]
                    errs_q = [X[data_qubit], Y[data_qubit], Z[data_qubit]]
                    errs_q_t = [f"X[{data_qubit}]",f"Y[{data_qubit}]",f"Z[{data_qubit}]"]                    
                    op1_index = errors_p.index(err_str)
                    w = errs_q[op1_index] * w * errs_q[op1_index]
                    print(errs_q_t[op1_index])

            w = w.unit()
            
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
            exit()
    exit()
    final_p_error.append([p_error, p_qnd, np.mean(result_correction), np.std(result_correction), np.mean(num_losses), np.mean(num_qnd_errors)])

np.savetxt(f"final_errors_{num_trials}_pqnd_{p_qnd:1.3f}_seed_{seme}.dat", final_p_error, fmt='%1.3f\t%1.3f\t' + '%1.6f\t' * 4)
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
plt.savefig(f"final_errors_{num_trials}_pqnd_{p_qnd:1.3f}_seed_{seme}.pdf")

plt.show()

