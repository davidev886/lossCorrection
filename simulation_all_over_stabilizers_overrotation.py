#
# Case 1: Coherent errors (single and two qubit overrotation)
# with no stabilizers measurement errors
#

import numpy as np
import qutip as qu
import os
from itertools import product
from utils.binary_conf import binary_configurations
from utils.p_operators_qutrit import *
from utils.overrotation_channel import (CorrelatedOverRotQubitAll,
                                        SingleOverRotQubitAll,
                                        SingleOverRotQutritAll
                                        )
import datetime
from utils.parameters import parse_command_line

import argparse

np.set_printoptions(precision=4, suppress=True)

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

p_overrot_2 = args.p_overrot_2
p_overrot_1 = args.p_overrot_1
eta = p_overrot_2 * np.pi
eps = p_overrot_1 * np.pi

folder_name = args.dir_name
num_trials = args.num_trials

choi_ideal = np.loadtxt("choi_op/choiFinal_ideal.dat")
choi_experiment = np.genfromtxt("choi_op/qubitqutrit_choi_noloss.csv",
                                dtype=complex,
                                delimiter=',')

if not os.path.exists(folder_name):
    os.makedirs(folder_name)

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

OverRotationOperators = CorrelatedOverRotQubitAll(p_overrot_2 * np.pi)
#SingleOverRotations = SingleOverRotQubitAll(p_overrot_1 * np.pi)

SingleOverRotations = SingleOverRotQutritAll(p_overrot_1 * np.pi)

now = datetime.datetime.now()
final_data_name = (now.strftime("%Y%m%d%H%M") +
                   f"_state_{LogicalStates_str[jLog]}_" +
                   f"phi_{phi_tilde:1.5f}_eps_{epsilon_choi}.dat"
                   )
file_data_name = os.path.join(folder_name,
                              final_data_name
                              )

print(f"logical state |{LogicalStates_str[jLog]}_L>")

index_confs = 0
final_prob = []
result_correction = []
num_losses = []

all_loss_events = binary_configurations().configurations_list

for outcomes_ancilla in all_loss_events:
    index_confs += 1

    prob_total_event = []
    prob_correction_logical_state = []
    psiL = LogicalStates[jLog]

    list_qubits = list(range(L))

    null_state = False
    rho_L = psiL * psiL.dag()

    for data_q in list_qubits:
        # apply Rloss with an angle phi
        rho_L = rotation_ops[data_q] * rho_L * rotation_ops[data_q].dag()
        # apply the QND detection unit
        rho_L = apply_qnd_process_unit(chi_matrix,
                                       rho_L,
                                       data_q,
                                       chi_threshold
                                       )

        # apply the overrotations
        if p_overrot_2 or p_overrot_1:
            rho_L = (OverRotationOperators[data_q] *
                     rho_L *
                     OverRotationOperators[data_q].dag()
                     )
            rho_L = (SingleOverRotations[data_q] *
                     rho_L *
                     SingleOverRotations[data_q].dag()
                     )

        if outcomes_ancilla[data_q] == 0:
            prob_outcome = (rho_L * Pp_ancilla).tr()
            if abs(prob_outcome.imag) > 1e-5:
                print("warning: im prob_outcome = {prob_outcome}")
            if prob_outcome == 0:
                # the state cannot be projected in the
                # +1 eigenstate of the ancilla
                null_state = True
                print("check null state outcome 0", prob_outcome)
                probs_outcome.append(prob_outcome)
                break
            else:
                prob_total_event.append(prob_outcome)
                rho_L = (Pp_ancilla *
                         rho_L *
                         Pp_ancilla.dag() / prob_outcome
                         )
        elif outcomes_ancilla[data_q] == 1:
            prob_outcome = (rho_L * Pm_ancilla).tr()
            if abs(prob_outcome.imag) > 1e-5:
                print("warning: im prob_outcome = {prob_outcome}")
            if prob_outcome == 0:
                null_state = True
                print("check null state outcome 1", prob_outcome)
                probs_outcome.append(prob_outcome)
                break
            else:
                prob_total_event.append(prob_outcome)
                rho_L = (Pm_ancilla *
                         rho_L *
                         Pm_ancilla.dag() / prob_outcome
                         )
                # reset ancilla
                rho_L = Xa * rho_L * Xa.dag()

        print(data_q,
              outcomes_ancilla[data_q],
              prob_outcome)
    prob_loss_event = np.abs(np.prod(prob_total_event))

    losses = np.where(outcomes_ancilla)[0].tolist()
    kept_qubits = list(set(range(L)) - set(losses))

    if sum(outcomes_ancilla) >= 7 or null_state:
        print(prob_loss_event)
        correction_successful = 0.0
        prob_correction_logical_state.append(correction_successful)
        conf_loss = int("".join(str(_) for _ in outcomes_ancilla))
        final_p_loss.append([phi_tilde,
                             conf_loss,
                             correction_successful,
                             prob_loss_event
                             ])
    else:
        w_0 = rho_L.ptrace(kept_qubits)
        rho_L = (qu.tensor([qu.fock_dm(3, 0)] * len(losses) +
                 [w_0] +
                 [qu.fock_dm(2, 0)])
                 )

        permutation_order_q = {}
        for j, el in enumerate(losses + kept_qubits):
            permutation_order_q[el] = j
        # print("permutation_order_q", permutation_order_q)

        stab_qubits_new_order = []
        for stab in stab_qubits:
            stab_qubits_new_order.append([permutation_order_q[q]
                                         for q in stab]
                                         )

        Sx = [X[j1] * X[j2] * X[j3] * X[j4]
              for j1, j2, j3, j4 in stab_qubits_new_order
              ]
        Sz = [Z[j1] * Z[j2] * Z[j3] * Z[j4]
              for j1, j2, j3, j4 in stab_qubits_new_order
              ]

        PPx = [[(Id + el) / 2, (Id - el) / 2] for el in Sx]
        PPz = [[(Id + el) / 2, (Id - el) / 2] for el in Sz]

        index_stab_measurement = 0

        cumulative_probability_stabilizers = 0.0

        for meas_binary_X, meas_binary_Z in product(range(8), range(8)):
            # if abs(1 - cumulative_probability_stabilizers) < 1e-4:
                # exit if the cumulative probability
                # during the stabilizer measurement
                # is already close to 1
            #    break
            state_after_measure = qu.Qobj(rho_L[:], dims=rho_L.dims)
            conf_str_X = bin(meas_binary_X)[2:].zfill(3)
            conf_int_X = [int(_) for _ in conf_str_X]
            conf_str_Z = bin(meas_binary_Z)[2:].zfill(3)
            conf_int_Z = [int(_) for _ in conf_str_Z]
            probability_each_measurement = []
            for stab_num, outcome_stab in enumerate(conf_int_X):
                prob = (PPx[stab_num][outcome_stab] *
                        state_after_measure).tr()
                if np.abs(prob) > 0:
                    state_after_measure = (PPx[stab_num][outcome_stab] *
                                           state_after_measure *
                                           PPx[stab_num][outcome_stab].dag() /
                                           prob)
                    probability_each_measurement.append(np.real(prob))
                else:
                    probability_each_measurement.append(0)

            for stab_num, outcome_stab in enumerate(conf_int_Z):
                prob = (PPz[stab_num][outcome_stab] * state_after_measure).tr()
                if np.abs(prob) > 0:
                    state_after_measure = (PPz[stab_num][outcome_stab] *
                                           state_after_measure *
                                           PPz[stab_num][outcome_stab].dag() /
                                           prob)
                    probability_each_measurement.append(np.real(prob))
                else:
                    probability_each_measurement.append(0)
            prob_stabilizers = np.prod(probability_each_measurement)
            cumulative_probability_stabilizers += prob_stabilizers
            exp_x = np.real(qu.expect(XL, state_after_measure))
            exp_z = np.real(qu.expect(ZL, state_after_measure))
            exp_y = np.real(qu.expect(1j * XL * ZL, state_after_measure))

            if jLog in (0, 1):
                correction_successful = (1 + abs(exp_z)) / 2
            elif jLog in (2, 3):
                correction_successful = (1 + abs(exp_x)) / 2
            elif jLog in (4, 5):
                correction_successful = (1 + abs(exp_y)) / 2

            if VERBOSE:
                print(conf_int_X,
                      conf_int_Z,
                      f"{prob_stabilizers:1.8f}",
                      f"{correction_successful:1.8f}",
                      f"{exp_z:+1.8}",
                      f"{exp_x:+1.8}",
                      f"{exp_y:+1.8}"
                      )

            prob_correction_logical_state.append(prob_stabilizers *
                                                 correction_successful
                                                 )
            conf_stab_meas = int("".join(str(_)
                                 for _ in conf_int_X + conf_int_Z)
                                 )

        prob_correction = np.real(np.sum(prob_correction_logical_state))
        conf_loss = int("".join(str(_) for _ in outcomes_ancilla))

        final_p_loss.append([phi_tilde,
                             conf_loss,
                             prob_correction,
                             prob_loss_event
                             ])

        np.savetxt(file_data_name,
                   final_p_loss,
                   fmt='%1.5f\t' +
                       '%07d\t' +
                       '%.18e\t' +
                       '%.18e\t')

np.savetxt(file_data_name,
           final_p_loss,
           fmt='%1.5f\t%07d\t%.18e\t%.18e\t')