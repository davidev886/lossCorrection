import datetime
now = datetime.datetime.now()
final_data_name = now.strftime("%Y%m%d%H%M")


import time

import argparse
#python process_matrix_simulation_all.py --phi_tilde  --epsilon_choi
parser = argparse.ArgumentParser(description = "Simulate qubit losses with QND measurement qubit+7qutrit system")
parser.add_argument('--phi_tilde',  type=float, default=0.1, help = "Rotation angle")
parser.add_argument('--epsilon_choi',  type=float, default=0.0, help = "epsilon_choi")
parser.add_argument('--logical_state',  type=int, default=0, help = "logical state integer corresponding to: 0, 1, +, -, +i, -i")
parser.add_argument('--chi_threshold',  type=float, default=0.0, help = "threshold for discarding Kraus operators in the chi matrix")
parser.add_argument('--p_err_stab',  type=float, default=0.05, help = "error in stab measurements")
args = parser.parse_args()

phi_tilde = args.phi_tilde
epsilon_choi = args.epsilon_choi
jLog = args.logical_state
chi_threshold = args.chi_threshold
p_err_stab = args.p_err_stab

import os
folder_name = f'chi_{chi_threshold:.01e}_eps_{epsilon_choi:1.3f}_p_stab_{p_err_stab:1.3f}'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

