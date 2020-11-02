import matplotlib.pyplot as plt
    
import numpy as np
#from qutip import *
import qutip as qu

from itertools import combinations

from utils.p_operators import *
from utils.corrections import check_correctable_state, check_correctable_state_analytics
from utils.qnd_error_gen import *

import os
if not os.path.exists("data"):
    os.makedirs("data")
    

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

#p_qnd = args.p_qnd
fp_prob = 0.01
fn_prob = fp_prob / 2
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
        qnd_errors = generate_qnd_error_fpn(random_losses, fp_prob, fn_prob)

        print("random_losses   ", random_losses)
        print("qnd_errors      ", qnd_errors)        

        [correctable_event, non_correctable_event] = check_correctable_state_analytics(random_losses, qnd_errors)
        print(check_correctable_state_analytics(random_losses, qnd_errors))

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

    final_p_error.append([p_loss, fp_prob, fn_prob, np.mean(result_correction), np.std(result_correction), np.mean(num_losses), np.mean(num_qnd_errors)])

np.savetxt(f"data/final_qnd_faulty_{num_trials}_fp_{fp_prob:1.3f}_fn_{fn_prob:1.3f}.dat", final_p_error, fmt='%1.5f\t%1.5f\t%1.5f\t' + '%1.8f\t' * 4)
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
plt.savefig(f"data/final_qnd_faulty_{num_trials}_fp_{fp_prob:1.3f}_fn_{fn_prob:1.3f}.pdf")
plt.title("$p_\mathrm{qnd} = " +  f"{fp_prob:1.3f} / " +  f"{fn_prob:1.3f}$")
plt.show()

exit()
