import argparse

# python process_matrix_simulation_all.py --phi_tilde  --epsilon_choi
def parse_command_line():
    parser = argparse.ArgumentParser(description="Simulate qubit losses"
                                     "with QND measurement qubit+7qutrit system")
    parser.add_argument('--phi_tilde',
                        type=float,
                        default=0.00,
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
                        default=4000,
                        help="Number of Monte Carlo samples"
                        )

    parser.add_argument('--verbose',
                        action="store_true",
                        help="verbose"
                        )
    args = parser.parse_args()
    return parser.parse_args()