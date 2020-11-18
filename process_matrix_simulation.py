import matplotlib.pyplot as plt
    
import numpy as np
import qutip as qu

from utils.qnd_error_gen import pick_qnd_error

import datetime
now = datetime.datetime.now()
final_data_name = now.strftime("%Y%m%d%H%M")


from utils.p_operators_qutrit import *

from random import randint
seme = randint(0,100)

#seme=72
np.random.seed(seme)

print("\nseed:", seme)



choi_ideal = np.loadtxt("choiFinal_ideal.dat")

final_p_loss = []

num_trials = 20


for phi_tilde in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.45]:
    phi = phi_tilde * np.pi 
    result_correction = []
    num_losses = []
    for trial in range(num_trials):
        print(f"trial={trial + 1: 5d} out of {num_trials}")    
        loss_pattern = []

        a_0 = np.random.random()  + np.random.random() * 1j
        b_0 = np.random.random()  + np.random.random() * 1j

        a = a_0 / np.sqrt(abs(a_0)**2 + abs(b_0)**2)
        b = b_0 / np.sqrt(abs(a_0)**2 + abs(b_0)**2)

        psiL = a * ZeroL + b * OneL

        for data_q in range(L):

            #apply Rloss with an angle phi
            rho_L = Rloss(psiL, phi, data_q)

            rho_L = apply_qnd_process_unit(choi_ideal, rho_L, qu_data = data_q)
    
        #    print(rho_L.tr())

            p_0 = (rho_L * Pp_ancilla).tr()
            if p_0 >= 1: p_0 = 1
            if p_0 <= 0: p_0 = 0
            QND_outcome = np.random.binomial(1, 1 - p_0)
            loss_pattern.append(QND_outcome)
            print(f"p_0={p_0:1.3f}", f"QND_outcome={QND_outcome}", f"p1={1-p_0:1.4f}",f"p_loss={np.sin(phi/2)**2/2:1.4f}")
            if QND_outcome == 0:
                rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / p_0
            elif QND_outcome == 1:
                rho_L = Pm_ancilla * rho_L * Pm_ancilla.dag() / (1 - p_0)    

        state_after_measure = rho_L

        stabZ_eigenvalues = []
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
        num_losses.append(sum(loss_pattern))
    


    print("result_correction")
    print(result_correction)


    final_p_loss.append([phi, np.mean(result_correction), np.std(result_correction), np.mean(num_losses)])
    np.savetxt(final_data_name + f"_loss_process_matrix_trials_{num_trials}.dat", final_p_loss, fmt='%1.6f')



np.savetxt(final_data_name + f"_loss_process_matrix_trials_{num_trials}.dat", final_p_loss, fmt='%1.6f')
