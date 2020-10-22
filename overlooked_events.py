import matplotlib.pyplot as plt
    
import numpy as np
#from qutip import *
import qutip as qu

from itertools import combinations

from utils.p_operators import *
from utils.qnd_detection import check_correctable_state


from random import randint
seme = randint(0,100)

np.random.seed(seme)

print("\nseed:", seme)

num_trials = 1000

final_p_error = []

import argparse
#python overlooked_events.py --p_qnd 
parser = argparse.ArgumentParser(description = "Simulate qubit losses and qnd error with a phenomenological model")
parser.add_argument('--p_qnd',  type=float, default=0.0, help = "probability of a false negative")
args = parser.parse_args()

p_qnd = args.p_qnd

for p_error in np.arange(0.0,0.9,0.05):
    trial = 0
    result_correction = []  
    num_losses = []
    num_qnd_errors = []
    while trial < num_trials:

#       1 means lost, 0 means not lost
        random_losses = np.random.binomial(1, p_error, 7)
#       1 means qnd error
        qnd_errors = np.random.binomial(1, p_qnd, 7)

        print("random_losses   ", random_losses)
        print("qnd_errors      ", qnd_errors)        

        dict_corr = check_correctable_state(random_losses, qnd_errors)
        
        correctable_event = dict_corr["correctable"]
        non_correctable_event = dict_corr["non_correctable_event"]
        event_to_check = (correctable_event == False) and (non_correctable_event == False)
        
        if correctable_event:
            trial += 1
            correction_successful = True
            result_correction.append(correction_successful + 0)
            num_losses.append(sum(random_losses))
            num_qnd_errors.append(sum(qnd_errors))
            print(f"{'correctable_event':30}", correctable_event)

            
            continue
        elif non_correctable_event:
            trial += 1
            correction_successful = False
            result_correction.append(correction_successful + 0)
            num_losses.append(sum(random_losses))
            num_qnd_errors.append(sum(qnd_errors))            
            print(f"{'NON correctable_event':30}", non_correctable_event)
            continue
        elif event_to_check:
            trial += 1
            #transform the wrongly detected losses in actual losses and make the correction            
            losses_binary = [(loss + qnd_err) for loss, qnd_err in zip(random_losses, qnd_errors)]
            
            losses = np.where(losses_binary)[0].tolist()            


            kept_qubits = list(set(range(L)) - set(losses))      
            print(f"{'TO CHECK':30}")
            print("losses          ", losses)
            print("kept_qubits     ", kept_qubits)
            
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

np.savetxt(f"data/final_qnd_faulty_{num_trials}_pqnd_{p_qnd:1.3f}.dat", final_p_error, fmt='%1.3f\t%1.3f\t' + '%1.6f\t' * 4)
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
plt.savefig(f"data/final_qnd_faulty_{num_trials}_pqnd_{p_qnd:1.3f}_seed_{seme}.pdf")
plt.title("$p_\mathrm{qnd} = " +  f"{p_qnd:1.2f}$")
plt.show()

exit()