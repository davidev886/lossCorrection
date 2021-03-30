import qutip as qu
import numpy as np

import os
from utils.p_operators_qutrit import *
from utils.binary_conf import create_random_event
import datetime
import argparse

np.set_printoptions(precision=4, suppress=True)

# python process_matrix_simulation_all.py --phi_tilde  --epsilon_choi
parser = argparse.ArgumentParser(description="Simulate qubit losses"
                                 "with QND measurement qubit+7qutrit system")
parser.add_argument('--phi_tilde',
                    type=float,
                    default=0.05,
                    help="Rotation angle"
                    )
parser.add_argument('--epsilon_choi',
                    type=float,
                    default=0.0,
                    help="epsilon_choi"
                    )
parser.add_argument('--logical_state',
                    type=int,
                    default=0,
                    help="logical state corresponding"
                         "to: 0, 1, +, -, +i, -i"
                    )
parser.add_argument('--chi_threshold',
                    type=float,
                    default=0.0,
                    help="threshold for discarding Kraus"
                         "operators in the chi matrix"
                    )
parser.add_argument('--dir_name',
                    type=str,
                    default="./",
                    help="directory for saving data"
                    )
parser.add_argument('--p_overrot_2',
                    type=float,
                    default=0.136,
                    help="over rotation MS gate"
                    )
parser.add_argument('--p_overrot_1',
                    type=float,
                    default=0.010,
                    help="over rotation single-qubit gates"
                    )
parser.add_argument('--num_trials',
                    type=int,
                    default=10000,
                    help="Number of Monte Carlo samples"
                    )
parser.add_argument('--p_err_stab',
                    type=float,
                    default=0.05,
                    help="error in stab measurements"
                    )

args = parser.parse_args()
phi_tilde = args.phi_tilde
phi = phi_tilde * np.pi
epsilon_choi = args.epsilon_choi
jLog = args.logical_state
chi_threshold = args.chi_threshold
eta = args.p_overrot_2 * np.pi
eps = args.p_overrot_1 * np.pi
folder_name = args.dir_name
num_trials = args.num_trials
p_err_stab = args.p_err_stab

choi_ideal = np.loadtxt("choiFinal_ideal.dat")

choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv",
                                dtype=complex,
                                delimiter=','
                                )

if not os.path.exists(folder_name):
    os.makedirs(folder_name)

choi = choi_ideal

T_matrix = give_transformation_matrix()
chi_matrix = get_chi_from_choi(choi, T_matrix)

final_p_loss = []

rotation_ops = Rloss_all(phi_tilde * np.pi)

LogicalStates_str = ["0", "1", "+", "-", "+i", "-i"]

basis_events = ([[0, _] for _ in range(4)] +
                [[1, _] for _ in range(2)]
                )

basic_event_str = {'0': (0, 0),
                   '1': (0, 1),
                   '2': (0, 2),
                   '3': (0, 3),
                   '4': (1, 0),
                   '5': (1, 1)
                   }

basic_event_probs = {'0': (1 - eps**2 / 2 - eta**2 / 4),
                     '1': eta**2 / 4,
                     '2': eps**2 / 4,
                     '3': eps**2 / 4,
                     '4': (1 - eps**2 / 4),
                     '5': eps**2 / 4
                     }

prob_loss = np.sin(phi / 2)**2 / 2

# all_events = product(basis_events, repeat=L)

# trial_list = [randrange(6**7) for _ in range(num_trials)]
# print(trial_list)

now = datetime.datetime.now()
final_data_name = (now.strftime("%Y%m%d%H%M") +
                   f"_state_{LogicalStates_str[jLog]}_" +
                   f"phi_{phi_tilde:1.2f}_eps_{epsilon_choi}.dat"
                   )
file_data_name = os.path.join(folder_name,
                              final_data_name
                              )

print(f"logical state |{LogicalStates_str[jLog]}_L>")

index_conf = 0
cumulative_probability = 0

done_events = []
while len(done_events) < num_trials:
    event, event_str = create_random_event(prob_loss, basic_event_probs)
    if event_str in done_events:
        continue
    done_events.append(event_str)
    print(index_conf, event_str)
    index_conf += 1
    if index_conf == num_trials: break

    outcomes_ancilla = [el[0] for el in event]
    sub_case_ancilla = [el[1] for el in event]

    print(event)
    print(outcomes_ancilla)

    prob_correction_logical_state = []
    psiL = LogicalStates[jLog]

    null_state = False
    rho_L = psiL * psiL.dag()
    do_nothing = []
    replace_qubits = []
    false_negative = []
    probs_outcome = []
    probs_incoherent_process = []
    for data_q in range(L):
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

        if outcomes_ancilla[data_q] == 0:  # ancilla in 0 state
            prob_outcome = (rho_L * Pp_ancilla).tr()
            if abs(prob_outcome.imag) > 1e-5:
                print("warning: im prob_outcome = {prob_outcome}")
            if prob_outcome == 0:
                # the state cannot be projected
                # in the +1 eigenstate of the ancilla
                null_state = True
                print("check null state")
                exit()
            if sub_case_ancilla[data_q] == 0:  # 1 - eps**2/2 - eta**2/4
                rho_L = (Pp_ancilla * rho_L * Pp_ancilla.dag() /
                         abs(prob_outcome))
                do_nothing.append(data_q)
                probs_incoherent_process.append(1 - eps**2 / 2 - eta**2 / 4)
            elif sub_case_ancilla[data_q] == 1:  # eta**2 / 4
                rho_L = (Pm_ancilla * rho_L * Pm_ancilla.dag() /
                         (1 - abs(prob_outcome)))
                rho_L = X[data_q] * rho_L * X[data_q].dag()
                rho_L = Xa * rho_L * Xa.dag()  # reinitializing ancilla
                replace_qubits.append(data_q)
                probs_incoherent_process.append(eta**2 / 4)
            elif sub_case_ancilla[data_q] == 2:  # epsilon**2/4
                rho_L = (Pm_ancilla * rho_L * Pm_ancilla.dag() /
                         (1 - abs(prob_outcome)))
                rho_L = Xa * rho_L * Xa.dag()  # reinitializing ancilla
                replace_qubits.append(data_q)
                probs_incoherent_process.append(eps**2 / 4)
            elif sub_case_ancilla[data_q] == 3:  # epsilon**2/4
                rho_L = (Pp_ancilla * rho_L * Pp_ancilla.dag()
                            / abs(prob_outcome))
                rho_L = X[data_q] * rho_L * X[data_q].dag()
                do_nothing.append(data_q)
                probs_incoherent_process.append(eps**2 / 4)
        elif outcomes_ancilla[data_q] == 1:  # ancilla in 1 state
            prob_outcome = (rho_L * Pm_ancilla).tr()
            if abs(prob_outcome.imag) > 1e-5:
                print("warning: im prob_outcome = {prob_outcome}")
            if prob_outcome == 0:
                null_state = True
                print("check null state")
                exit()
            if sub_case_ancilla[data_q] == 0:  # 1 - eps**2 / 4
                rho_L = (Pm_ancilla * rho_L * Pm_ancilla.dag()
                        / abs(prob_outcome))
                rho_L = Xa * rho_L * Xa.dag()  # reinitializing ancilla
                replace_qubits.append(data_q)
                probs_incoherent_process.append(1 - eps**2 / 4)
            elif sub_case_ancilla[data_q] == 1:  # eps**2 / 4 false negative
                rho_L = (Pp_ancilla * rho_L * Pp_ancilla.dag()
                        / (1-abs(prob_outcome)))
                rho_L = Xa * rho_L * Xa.dag()  # reinitializing ancilla
                do_nothing.append(data_q)
                probs_incoherent_process.append(eps**2 / 4)

        probs_outcome.append(prob_outcome)

    prob_total_event = np.prod(probs_outcome) * np.prod(probs_incoherent_process)
    cumulative_probability += prob_total_event
    print("probs_outcome", np.array(probs_outcome))
    print("probs_incoherent_process", np.array(probs_incoherent_process))
    print(index_conf, outcomes_ancilla, sub_case_ancilla,
          do_nothing,
          replace_qubits,
          false_negative,
          # np.array(probs_outcome),
          # np.array(probs_incoherent_process),
          # np.prod(probs_outcome),
          # f"{np.prod(probs_incoherent_process):4}",
          f"{np.prod(probs_outcome)*np.prod(probs_incoherent_process):.4}",
          f"{cumulative_probability:.4}"
          )
    losses = replace_qubits

    if sum(outcomes_ancilla) >= 7 or null_state or len(do_nothing) == 0:
        print(prob_total_event)
        correction_successful = 0.0
        prob_correction_logical_state.append(correction_successful)
    else:
        w_0 = rho_L.ptrace(do_nothing)
        rho_L = qu.tensor([qu.fock_dm(3,0)] * len(replace_qubits)
                           + [qu.fock_dm(3,2)] * len(false_negative)
                           + [w_0]
                           + [qu.fock_dm(2,0)])

        print(replace_qubits,
              false_negative,
              do_nothing
              )
        permutation_order_q = {}
        # the order in the for is important because redefine the state as
        # ket(0) detected_losses , ket(2) false negative, kept_qubits
        for j, el in enumerate(replace_qubits + false_negative + do_nothing):
            permutation_order_q[el] = j
            # print("permutation_order_q", permutation_order_q)

        stab_qubits_new_order = []
        for stab in stab_qubits:
            stab_qubits_new_order.append([permutation_order_q[q] for q in stab])

        Sx = [X[j1] * X[j2] * X[j3] * X[j4] for j1,j2,j3,j4 in stab_qubits_new_order]
        Sz = [Z[j1] * Z[j2] * Z[j3] * Z[j4] for j1,j2,j3,j4 in stab_qubits_new_order]

        PPx = [[(Id + el) / 2, (Id - el) / 2] for el in Sx]
        PPz = [[(Id + el) / 2, (Id - el) / 2] for el in Sz]

        average_value_each_stab_meas = []
        correction_each_measurement = []
        index_stab_measurement = 0

        cumulative_probability_stabilizers = 0.0
        index_stab = 0
        for meas_binary_X, meas_binary_Z in product(range(8), range(8)):
            if abs(1 - cumulative_probability_stabilizers) < 1e-4:
                # exit if the cumulative probability
                # during the stabilizer measurement
                # is already close to 1
                break
            state_after_measure = qu.Qobj(rho_L[:], dims = rho_L.dims)
            conf_str_X = bin(meas_binary_X)[2:].zfill(3)
            conf_int_X = [int(_) for _ in conf_str_X]
            conf_str_Z = bin(meas_binary_Z)[2:].zfill(3)
            conf_int_Z = [int(_) for _ in conf_str_Z]
            probability_each_measurement = []
            for stab_num, outcome_stab in enumerate(conf_int_X):
                prob = (PPx[stab_num][outcome_stab] * state_after_measure).tr()
                if np.abs(prob) > 0:
                    state_after_measure = (PPx[stab_num][outcome_stab] *
                                            state_after_measure *
                                          PPx[stab_num][outcome_stab].dag() / prob)
                    probability_each_measurement.append(np.real(prob))
                else:
                    probability_each_measurement.append(0)

            for stab_num, outcome_stab in enumerate(conf_int_Z):
                prob = (PPz[stab_num][outcome_stab] * state_after_measure).tr()
                if np.abs(prob) > 0:
                    state_after_measure = (PPz[stab_num][outcome_stab] *
                                           state_after_measure *
                                           PPz[stab_num][outcome_stab].dag() / prob)
                    probability_each_measurement.append(np.real(prob))
                else:
                    probability_each_measurement.append(0)

            # place where we can apply corrections but we don't

            print(f"{index_stab_measurement: 4d}",
                  conf_int_X,
                  conf_int_Z,
                  f"{np.prod(probability_each_measurement):1.4f}",
                  f"{qu.expect(XL, state_after_measure):+1.4f}",
                  f"{ qu.expect(ZL, state_after_measure):+1.4f}",
                  f"{qu.expect(1j * XL * ZL, state_after_measure):+1.4f}"
                  )
            prob_stabilizers = np.prod(probability_each_measurement)
            cumulative_probability_stabilizers += prob_stabilizers

            successful_correction_total = []
            for error_in_stab_X, error_in_stab_Z in product(range(8), range(8)):

                conf_str_error_X = bin(error_in_stab_X)[2:].zfill(3)
                conf_int_error_X = [int(_) for _ in conf_str_error_X]
                conf_str_error_Z = bin(error_in_stab_Z)[2:].zfill(3)
                conf_int_error_Z = [int(_) for _ in conf_str_error_Z]

                stabX_errors = np.array(conf_int_error_X)
                stabZ_errors = np.array(conf_int_error_Z)
                # check if a loss happens on a faulty stabilizer
                # if so, the state is not correctable
                # (as measuring a stabilizer on a lost qubit is impossible)
                faulty_stabX_qubits = [el for j,el in enumerate(stab_qubits)
                                       if stabX_errors[j]
                                       ]
                faulty_stabZ_qubits = [el for j,el in enumerate(stab_qubits)
                                       if stabZ_errors[j]
                                       ]
                # print("stab_errors", stabX_errors, stabZ_errors)
                # print(faulty_stabX_qubits, faulty_stabZ_qubits)

                loss_on_faulty_stabX = any([any([(loss in stab)
                                            for loss in losses])
                                            for stab in faulty_stabX_qubits]
                                            )
                loss_on_faulty_stabZ = any([any([(loss in stab)
                                            for loss in losses])
                                            for stab in faulty_stabZ_qubits]
                                            )
                prob_stab_event_X = np.prod(p_err_stab**stabX_errors *
                                            (1 - p_err_stab)**(1 - stabX_errors)
                                            )
                prob_stab_event_Z = np.prod(p_err_stab**stabZ_errors *
                                            (1 - p_err_stab)**(1 - stabZ_errors)
                                            )
                prob_stab_error_event = prob_stab_event_X * prob_stab_event_Z


                if loss_on_faulty_stabX or loss_on_faulty_stabZ:
                    # print(faulty_stabX_qubits, [([(loss in stab) for loss in losses]) for stab in faulty_stabX_qubits])
                    # print(faulty_stabZ_qubits, [([(loss in stab) for loss in losses]) for stab in faulty_stabZ_qubits])
                    correction_successful = 0.0
                else:
                    if jLog in (0, 1):
                        correction_successful = (1 + abs(qu.expect(ZL, state_after_measure))) / 2
                    elif jLog in (2, 3):
                        correction_successful = (1 + abs(qu.expect(XL, state_after_measure))) / 2
                    elif jLog in (4, 5):
                        correction_successful = (1 + abs(qu.expect(1j * XL * ZL, state_after_measure))) / 2
                conf_stab_errors = int("".join(str(_) for _ in conf_int_error_X + conf_int_error_Z))
                conf_stab_meas = int("".join(str(_) for _ in conf_int_X + conf_int_Z))
                # print("------------------>",
                #      f"{prob_stabilizers:1.5f}",
                #      f"{prob_stab_error_event:1.5f}",
                #      f"{correction_successful:1.5f}",
                #      f"{conf_stab_errors:06d}" )

                average_value_each_stab_meas.append(prob_stabilizers *
                                                    prob_stab_error_event *
                                                    correction_successful
                                                    )
            index_stab += 1

        print("prob_of_succ_correction", np.sum(average_value_each_stab_meas))
        conf_loss = int("".join(str(_) for _ in outcomes_ancilla))
        conf_case_ancilla = int("".join(str(_) for _ in sub_case_ancilla))
        res = ([phi_tilde,
                conf_loss,
                conf_case_ancilla,
                np.real(np.sum(average_value_each_stab_meas)),
                np.real(prob_total_event)
                ])

        final_p_loss.append(res)
        np.savetxt(file_data_name, final_p_loss, fmt='%1.3f\t' +
                                                     '%07d\t' +
                                                     '%07d\t' +
                                                     '%.10e\t' +
                                                     '%1.14f\t')

np.savetxt(file_data_name, final_p_loss, fmt='%1.3f\t' +
                                             '%07d\t' +
                                             '%07d\t' +
                                             '%.10e\t' +
                                             '%1.14f\t')
