import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu

from itertools import combinations

from utils.p_operators import *
from utils.corrections import check_correctable_state, check_stabilizers_measurement

import os
if not os.path.exists("data"):
    os.makedirs("data")



from random import randint
seme = randint(0,100)

np.random.seed(seme)

print("\nseed:", seme)

num_trials = 5000

final_p_error = []

import argparse
#python overlooked_events.py --p_qnd 
parser = argparse.ArgumentParser(description = "Simulate qubit losses and qnd error with a phenomenological model")
parser.add_argument('--p_qnd',  type=float, default=0.0, help = "probability of a false negative")
parser.add_argument('--p_stab',  type=float, default=0.0, help = "probability of a measurement error")
args = parser.parse_args()

p_qnd = args.p_qnd

p_stab = args.p_stab

#for p_loss in np.arange(0,0.2+0.005, 0.005): #np.arange(0.0, 0.2 + 1e-5, 1e-4): #
for p_loss in np.arange(0.0,0.95,0.05):
    trial = 0
    result_correction = []  
    num_losses = []
    num_qnd_errors = []
    while trial < num_trials:

#       1 means lost, 0 means not lost
        random_losses = np.random.binomial(1, p_loss, 7)
        
#       1 means qnd error
        qnd_errors = np.random.binomial(1, p_qnd, 7)

#       1 means measurement error
        stab_errors_binary = np.random.binomial(1, p_stab, 3)


        print(f"{'random_losses':40}", random_losses)
        print(f"{'qnd_errors':40}", qnd_errors)
        print(f"{'stab_errors_binary':40}", stab_errors_binary)

        [correctable_event, non_correctable_event] = check_stabilizers_measurement(random_losses, qnd_errors, stab_errors_binary)
        
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


    final_p_error.append([p_loss, p_qnd, p_stab, np.mean(result_correction), np.std(result_correction), np.mean(num_losses), np.mean(num_qnd_errors)])

np.savetxt(f"data/final_qnd_faulty_{num_trials}_pqnd_{p_qnd:1.3f}_pstab_{p_stab:1.3f}.dat", final_p_error, fmt='%1.5f\t%1.5f\t%1.5f\t' + '%1.8f\t' * 4)
final = np.array(final_p_error)

import matplotlib.pyplot as plt

x_data = final[:,0]
y_data = final[:,3]
y_error = final[:,4] / np.sqrt(num_trials)
plt.errorbar(x_data, y_data,  yerr=y_error, fmt='o-')


x_data = np.linspace(0,1,100)

y_data = 1 - 7*x_data**3 + 21*x_data**5 - 21*x_data**6 + 6*x_data**7
plt.plot(x_data, y_data, '-')

y_data = 1 - x_data
plt.plot(x_data, y_data, '-')


plt.xlabel("p")
plt.ylabel("p(success)")
plt.title("$p_\mathrm{qnd} = " +  f"{p_qnd:1.4f}$, " + " $p_\mathrm{st} = " +  f"{p_stab:1.4f}$")

plt.savefig(f"data/final_qnd_faulty_{num_trials}_pqnd_{p_qnd:1.3f}_pstab_{p_stab:1.3f}_seed_{seme}.pdf")
plt.show()

exit()
