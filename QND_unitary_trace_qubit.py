#
# Case 2: incoherent model for the overrotations
# with no loss rotations
# run on the first most probable events
#

import qutip as qu
import numpy as np
from itertools import product

import os
from utils.p_operators_qutrit import *
# from utils.p_operators import *
from utils.binary_conf import create_random_event
from utils.parameters import parse_command_line
from utils.incoherent_channel_qutrit import inc_channel
import datetime


import glob

np.set_printoptions(precision=4, suppress=True)

args = parse_command_line()

phi_tilde = args.phi_tilde
print(phi_tilde)
phi = phi_tilde * np.pi
epsilon_choi = args.epsilon_choi
jLog = args.logical_state
chi_threshold = args.chi_threshold
eta = args.p_overrot_2 * np.pi
eps = args.p_overrot_1 * np.pi
folder_name = args.dir_name
num_trials = args.num_trials
VERBOSE = args.verbose
choi_ideal = np.loadtxt("choi_op/choiFinal_ideal.dat")

choi_experiment = np.genfromtxt("choi_op/qubitqutrit_choi_noloss.csv",
                                dtype=complex,
                                delimiter=','
                                )

if not os.path.exists(folder_name):
    os.makedirs(folder_name)

choi = choi_ideal

# T_matrix = give_transformation_matrix()
# chi_matrix = get_chi_from_choi(choi, T_matrix)

final_p_loss = []

# rotation_ops = Rloss_all(phi_tilde * np.pi)
choi = np.real((1 - epsilon_choi) * choi_ideal +
               6 * epsilon_choi * choi_experiment
               )

T_matrix = give_transformation_matrix()
chi_matrix = get_chi_from_choi(choi, T_matrix)
chi_matrix[np.abs(chi_matrix) < chi_threshold] = 0

final_p_loss = []

rotation_ops = Rloss_all(phi_tilde * np.pi)
result_correction = []

LogicalStates_str = ["0", "1", "+", "-", "+i", "-i"]

now = datetime.datetime.now()
final_data_name = (now.strftime("%Y%m%d%H%M") +
                   f"_state_{LogicalStates_str[jLog]}_" +
                   f"phi_{phi_tilde:1.5f}_eps_{epsilon_choi}.dat"
                   )

basic_event_probs = {'1a': (3 + np.cos(2*eps) + 4*np.cos(eps)*np.cos(eta))/8.,  # 1(a)
                     '1b': np.sin(eps)**2/4.,  # 1(b)
                     '1c': (3 + np.cos(2*eps) - 4*np.cos(eps)*np.cos(eta))/8.,  # 1(c)
                     '1d': np.sin(eps)**2/4.,  # 1(d)
                     '2a': np.cos(eps/2.)**2,  # 2(a)
                     '2b': np.sin(eps/2.)**2,  # 2(b)
                     }
# try with these probability inverted as from Mathematica notebook
# original definition
#  If the ancilla is found in 1 (loss detected) we
# (a) leave the qutrit state (ρ = |2⟩ ⟨2|) as it is with probability np.cos(eps/2.)**2
# (b) flip the ancilla to the no-loss detection state with probability np.sin(eps/2.)**2  (false negative).
#                     '2a': np.cos(eps/2.)**2,  # 2(a)
#                     '2b': np.sin(eps/2.)**2,  # 2(b)


prob_loss = np.sin(phi / 2)**2 / 2

all_events = product([0, 1], repeat=L)

file_data_name = os.path.join(folder_name,
                              final_data_name
                              )

# get the incoherent channels once
channels_0 = []
channels_1 = []
# for data_q in range(L):
#     print(data_q)
#     ch_0, ch_1 = inc_channel(basic_event_probs,
#                              data_q,
#                              noloss=False
#                              )
#     channels_0.append(ch_0)
#     channels_1.append(ch_1)

print(f"logical state |{LogicalStates_str[jLog]}_L>")
print(basic_event_probs)
index_conf = 0
cumulative_probability = 0
for outcomes_ancilla in all_events:

    psiL = LogicalStates[jLog]

    null_state = False
    rho_L = psiL * psiL.dag()
    do_nothing = []
    replace_qubits = []

    probs_outcome = []
    probs_incoherent_process = []
    print(outcomes_ancilla)
    ancilla_after_channel = []
    for data_q in range(L):
#        print(data_q)
        # apply Rloss with an angle phi
        rho_L = rotation_ops[data_q] * rho_L * rotation_ops[data_q].dag()
        # apply the QND detection unit
        rho_L = apply_qnd_process_unit(chi_matrix,
                             rho_L,
                             data_q,
                             chi_threshold
                             )
        # rho_L.tidyup(atol = 1e-8)

        # apply the effective incoherent noise model

        Xq = X[data_q]
#         rho_Lp = rho_L #Pp_ancilla * rho_L * Pp_ancilla.dag()
#         rho_Lm = rho_L #Pm_ancilla * rho_L * Pm_ancilla.dag()
#         channel_0 = (basic_event_probs['1a'] * Id * rho_Lp * Id +
#                      basic_event_probs['1b'] * Xq * rho_Lp  * Xq +
#                      basic_event_probs['1c'] * Xq * Xa * rho_Lp * Xq * Xa +
#                      basic_event_probs['1d'] * Xa * rho_Lp * Xa
#                      )
#
#         channel_1 = (basic_event_probs['2a'] * Id * rho_Lm * Id +
#                      basic_event_probs['2b'] * Xa * rho_Lm * Xa
#                      )
#        print("prob 0 on channel 0", data_q, (channel_0 * Pp_ancilla).tr())
#        print("prob 1 on channel 1", data_q, (channel_1 * Pp_ancilla).tr())


        if outcomes_ancilla[data_q] == 0:  # ancilla in 0 state before the channel
            prob_outcome = (rho_L * Pp_ancilla).tr()
            probs_outcome.append(prob_outcome)
            if abs(prob_outcome.imag) > 1e-5:
                print("warning: in prob_outcome = {prob_outcome}")
            if prob_outcome == 0:
                # the state cannot be projected
                # in the +1 eigenstate of the ancilla
                null_state = True
                print("check null state")
                break  # exit()

            rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / prob_outcome
            print("prob_outcome 0" , prob_outcome)
           # rho_Lp = Pp_ancilla * rho_L * Pp_ancilla.dag()
            rho_L = (basic_event_probs['1a'] * Id * rho_L * Id +
                     basic_event_probs['1b'] * Xq * rho_L  * Xq +
                     basic_event_probs['1c'] * Xq * Xa * rho_L * Xq * Xa +
                     basic_event_probs['1d'] * Xa * rho_L * Xa
                     )
            print("prob 0 on channel 0 rho", data_q, (rho_L * Pp_ancilla).tr())
            print("prob 1 on channel 0 rho", data_q, (rho_L * Pm_ancilla).tr())
            ancilla_final = np.random.binomial(1, 1 - (rho_L * Pp_ancilla).tr())
            ancilla_after_channel.append(ancilla_final)
#             if ancilla_final == 0:
#                 do_nothing.append(data_q)
#             elif ancilla_final == 1:
#                 replace_qubits.append(data_q)

        elif outcomes_ancilla[data_q] == 1:  # ancilla in 1 state
            prob_outcome = (rho_L * Pm_ancilla).tr()
#            rho_Lm = Pm_ancilla * rho_L * Pm_ancilla.dag()
            probs_outcome.append(prob_outcome)
            if abs(prob_outcome.imag) > 1e-5:
                print("warning: in prob_outcome = {prob_outcome}")
            if prob_outcome == 0:
                # the state cannot be projected
                # in the +1 eigenstate of the ancilla
                null_state = True
                print("check null state outcome 1", prob_outcome)
                break  # exit()
            # print("before 0 Pm_ancilla rho_L.tr()" , rho_L.tr(), prob_outcome)
            rho_L = Pm_ancilla * rho_L * Pm_ancilla.dag() / prob_outcome
            print("prob_outcome 1" , prob_outcome)
            replace_qubits.append(data_q)
            rho_L = (basic_event_probs['2a'] * Id * rho_L * Id +
                     basic_event_probs['2b'] * Xa * rho_L * Xa
                     )
            print("prob 0 on channel 1 rho", data_q, (rho_L * Pp_ancilla).tr())
            print("prob 1 on channel 1 rho", data_q, (rho_L * Pm_ancilla).tr())
            ancilla_final = np.random.binomial(1, 1 - (rho_L * Pp_ancilla).tr())
            ancilla_after_channel.append(ancilla_final)
#             if ancilla_final == 0:
#                 do_nothing.append(data_q)
#             elif ancilla_final == 1:
#                 replace_qubits.append(data_q)
            rho_L = Xa * rho_L * Xa.dag()  # reinitializing ancilla


        # renormalize if the trace of rho_L is bigger than 1
        # because of accumulated errorr
        traccia = rho_L.tr()
        if traccia > 1:
            rho_L = rho_L / traccia

    print("ancilla_final", ancilla_after_channel)
    ancilla_after_channel_arr = np.array(ancilla_after_channel)
    replace_qubits = np.where(ancilla_after_channel_arr == 1)[0].tolist()
    do_nothing = np.where(ancilla_after_channel_arr == 0)[0].tolist()
    print("replace_qubits", replace_qubits)
    print("do_nothing", do_nothing)
    prob_total_event = np.prod(probs_outcome)
    cumulative_probability += prob_total_event
    print("probs_outcome", np.array(probs_outcome))

    if sum(outcomes_ancilla) >= 7 or null_state or len(do_nothing) == 0:
        correction_successful = 0.0

        conf_loss = int("".join(str(_) for _ in outcomes_ancilla))
        conf_ancilla = int("".join(str(_) for _ in ancilla_after_channel))
        res = ([phi_tilde,
                conf_loss,
                correction_successful,
                np.real(prob_total_event),
                conf_ancilla
                ])
        print(index_conf, outcomes_ancilla,
              "null_state"
              )
    else:
        print(index_conf, outcomes_ancilla,
              do_nothing,
              replace_qubits,
              f"{prob_total_event:.4}",
              f"{1-cumulative_probability:.4e}"
              )
        w_0 = rho_L.ptrace(do_nothing)
        rho_L = qu.tensor([qu.fock_dm(dimQ, 0)] * len(replace_qubits)
                           + [w_0]
                           + [qu.fock_dm(dimA, 0)])

        print(replace_qubits,
              do_nothing
              )
        permutation_order_q = {}
        # the order in the for is important because redefine the state as
        # ket(0) detected_losses , ket(2), kept_qubits
        for j, el in enumerate(replace_qubits + do_nothing):
            permutation_order_q[el] = j
            # print("permutation_order_q", permutation_order_q)

        stab_qubits_new_order = []
        for stab in stab_qubits:
            stab_qubits_new_order.append([permutation_order_q[q] for q in stab])

        Sx = [X[j1] * X[j2] * X[j3] * X[j4] for j1, j2, j3, j4 in stab_qubits_new_order]
        Sz = [Z[j1] * Z[j2] * Z[j3] * Z[j4] for j1, j2, j3, j4 in stab_qubits_new_order]

        PPx = [[(Id + el) / 2, (Id - el) / 2] for el in Sx]
        PPz = [[(Id + el) / 2, (Id - el) / 2] for el in Sz]

        average_value_each_stab_meas = []
        correction_each_measurement = []
        index_stab_measurement = 0

        cumulative_probability_stabilizers = 0.0

        for meas_binary_X, meas_binary_Z in product(range(8), range(8)):
            # if abs(1 - cumulative_probability_stabilizers) < 1e-4:
                # exit if the cumulative probability
                # during the stabilizer measurement
                # is already close to 1
                # break
            state_after_measure = qu.Qobj(rho_L[:], dims=rho_L.dims)
            configuration_str_X = bin(meas_binary_X)[2:].zfill(3)
            configuration_int_X = [int(_) for _ in configuration_str_X]
            configuration_str_Z = bin(meas_binary_Z)[2:].zfill(3)
            configuration_int_Z = [int(_) for _ in configuration_str_Z]
            probability_each_measurement = []
            for stab_num, outcome_stab in enumerate(configuration_int_X):
                prob = (PPx[stab_num][outcome_stab] *
                        state_after_measure).tr()
                if np.abs(prob) > 0:
                    state_after_measure = (PPx[stab_num][outcome_stab] *
                                            state_after_measure *
                                          PPx[stab_num][outcome_stab].dag() / prob)
                    probability_each_measurement.append(np.real(prob))
                else:
                    probability_each_measurement.append(0)

            for stab_num, outcome_stab in enumerate(configuration_int_Z):
                prob = (PPz[stab_num][outcome_stab] * state_after_measure).tr()
                if np.abs(prob) > 0:
                    state_after_measure = (PPz[stab_num][outcome_stab] *
                                           state_after_measure *
                                           PPz[stab_num][outcome_stab].dag() / prob)
                    probability_each_measurement.append(np.real(prob))
                else:
                    probability_each_measurement.append(0)

            # place where we can apply corrections but we don't

            if VERBOSE:
                print(f"{index_stab_measurement: 4d}",
                      configuration_int_X,
                      configuration_int_Z,
                      f"{np.prod(probability_each_measurement):1.4f}",
                      f"{qu.expect(XL, state_after_measure):+1.4f}",
                      f"{qu.expect(ZL, state_after_measure):+1.4f}",
                      f"{qu.expect(1j * XL * ZL, state_after_measure):+1.4f}"
                      )

            prob_stabilizers = np.prod(probability_each_measurement)
            cumulative_probability_stabilizers += prob_stabilizers

            if jLog in (0, 1):
                correction_successful = (1 + abs(qu.expect(ZL, state_after_measure))) / 2
            elif jLog in (2, 3):
                correction_successful = (1 + abs(qu.expect(XL, state_after_measure))) / 2
            elif jLog in (4, 5):
                correction_successful = (1 + abs(qu.expect(1j * XL * ZL, state_after_measure))) / 2

            average_value_each_stab_meas.append(prob_stabilizers *
                                                correction_successful
                                                )
            # conf_stab_meas = int("".join(str(_) for _ in configuration_int_X + configuration_int_Z))
            index_stab_measurement += 1

        print("cumulative_probability_stabilizers: ",
              f"{cumulative_probability_stabilizers:.4f}")
        print("1-cumulative_probability_stabilizers: ",
              f"{1-cumulative_probability_stabilizers:.4e}")
        print("prob_of_succ_correction", np.sum(average_value_each_stab_meas))
        conf_loss = int("".join(str(_) for _ in outcomes_ancilla))
        conf_ancilla = int("".join(str(_) for _ in ancilla_after_channel))
        res = ([phi_tilde,
                conf_loss,
                np.real(np.sum(average_value_each_stab_meas)),
                np.real(prob_total_event),
                conf_ancilla
                ])

    index_conf += 1
    final_p_loss.append(res)
    np.savetxt(file_data_name, final_p_loss, fmt='%1.5f\t' +
                                                 '%07d\t' +
                                                 '%.18e\t' +
                                                 '%.18e\t' +
                                                 '%07d\t')


# A = np.array(final_p_loss)
# prob_success = A[:, 2]
# prob_event = A[:, 3]
# print("sum(prob_event)", sum(prob_event))
# p_success = np.sum(prob_success * prob_event)
# print("p_success   ", f"{p_success:.4f}")
# print("log_err_rate", f"{1-p_success:.4e}")
#
