import numpy as np
from itertools import product

def pick_qnd_error(p_qnd):

    random_toss_qnd = np.random.binomial(1, p_qnd, 1)
    if not (sum(random_toss_qnd)):
        #no qnd_error
        return "II"
    else:
        depol_errors = [e1 + e2 for e1, e2 in product("IXYZ", "IXYZ")]
        depol_errors.remove("II")    
        rng = np.random.default_rng()
        return rng.choice(depol_errors)