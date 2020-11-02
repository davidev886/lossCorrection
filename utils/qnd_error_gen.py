import numpy as np
from itertools import product


def generate_qnd_error_fpn(random_losses, fp_prob, fn_prob, num_qubits = 7):
    """
    generate a pattern for qnd error detection with false positive and false negative
    """

    false_positive = np.random.binomial(1, fp_prob, num_qubits)
    false_negative = np.random.binomial(1, fn_prob, num_qubits)

    qnd_error = [0] * num_qubits
    for j,x in enumerate(random_losses):
        if x == 0:
            if false_positive[j] == 1:
                #make a qnd error by detecting a loss even if it has not happened (false positive)
                qnd_error[j] = 1
        elif x == 1:
            if false_negative[j] == 1:
                #make a qnd error by discarding this loss even if it has happened (false negative)
                qnd_error[j] = 1
    return qnd_error
        


def pick_qnd_error(p_qnd):
    """
    generate a depolarizing channel with probability p_qnd
    returns a string with "Eq Ea" with a Pauli operator on the data qubit and the ancilla
    representing the error
    """
    random_toss_qnd = np.random.binomial(1, p_qnd, 1)
    if not (sum(random_toss_qnd)):
        #no qnd_error
        return "II"
    else:
        depol_errors = [e1 + e2 for e1, e2 in product("IXYZ", "IXYZ")]
        depol_errors.remove("II")    
        rng = np.random.default_rng()
        return rng.choice(depol_errors)