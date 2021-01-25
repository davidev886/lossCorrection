import qutip as qu
import numpy as np
import scipy.linalg
from scipy.linalg import eigh, svd

from itertools import product
import matplotlib.pyplot as plt


def proj(ket, bra, dimH = 2):
    """
    Define basis projector operators|ket >< bra|
    """
    if isinstance(ket, str):
        states_ket = [int(_) for _ in ket]
        ket_s = qu.basis([dimH] * len(states_ket), states_ket)        
    elif isinstance(ket, int):
        states_ket = ket
        ket_s = qu.basis(dimH, states_ket)        
    if isinstance(bra, str):
        states_bra = [int(_) for _ in bra]
        bra_s = qu.basis([dimH] * len(states_bra), states_bra).dag()
    elif isinstance(bra, int):
        states_bra = bra
        bra_s = qu.basis(dimH, states_bra).dag()
           
    return ket_s * bra_s



if __name__ == "__main__":
    np.set_printoptions(precision = 4, suppress = True, linewidth=100000)

    choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')

#Define the ideal channel with A0 and A1 operators:
#    A0 = |1><1|+|0><0|
    A0 = proj(1, 1, dimH = 3) + proj(0, 0, dimH = 3)

#    A1 = |2><2|
    A1 = proj(2, 2, dimH = 3)

    
#   U_ideal = 1_a  A0 + X_a  A1
    U_ideal = qu.tensor(qu.qeye(2), A0) + qu.tensor(qu.sigmax(), A1)
#   U_ideal = 
#     [[1. 0. 0. 0. 0. 0.]
#      [0. 1. 0. 0. 0. 0.]
#      [0. 0. 0. 0. 0. 1.]
#      [0. 0. 0. 1. 0. 0.]
#      [0. 0. 0. 0. 1. 0.]
#      [0. 0. 1. 0. 0. 0.]]
    
# Define Omega state for computing Choi matrix 
# omega = [00,01,02,10,11,12] * [00,01,02,10,11,12] ancilla first
    omega = sum([qu.tensor(qu.basis(2, i), qu.basis(3, j), qu.basis(2, i), qu.basis(3, j)) 
                    for i, j in product(range(2), range(3))])
    dimHtot= 6
    omega_state = omega * omega.dag()
          
#    QND_ideal =  qu.tensor(qu.tensor(qu.qeye(2),qu.qeye(3)), U_ideal) * omega_state * qu.tensor(qu.tensor(qu.qeye(2),qu.qeye(3)), U_ideal.dag())        
    

    p_dep = 1


    choi_depolarizing = (p_dep * qu.tensor(qu.qeye(2), qu.qeye(3), qu.qeye(2) , qu.qeye(3))/dimHtot 
                            + (1 - p_dep) * omega_state)

    vals, vecs = eigh(choi_depolarizing.data.todense())

    kraus_matrix = [np.sqrt(vals[j]) * vecs[:, j].reshape((dimHtot, dimHtot)).T for j in range(len(vecs))]
#    print(sum([np.dot(np.conjugate(Op.T), Op) for Op in kraus_matrix]))    

    system_length = len(choi_depolarizing.dims[0])//2
    dim_kraus = [choi_depolarizing.dims[0][-system_length:], choi_depolarizing.dims[0][-system_length:]]   
    kraus_ops = [qu.Qobj(inpt=op, dims=dim_kraus) for op in kraus_matrix]

    channel_dep = sum([qu.spre(op) * qu.spost(op.dag()) for op in kraus_ops])
    print(channel_dep)
    print()

    a = 1/np.sqrt(2)
    b = 1/np.sqrt(2)
#    a = 1
#    b = 0
    state_qutrit =qu.tensor(qu.basis(2,0),  (a * qu.basis(3,0) + b * qu.basis(3,1))).unit()
    rho = state_qutrit * state_qutrit.dag()
    print(channel_dep(state_qutrit))
#    state_0 = (qu.tensor([state_qutrit, qu.basis(2,0)])).unit()


#    print(    qu.process_fidelity(choi_ideal, choi_full_mixed))
    
    exit()


    distances_list = []
    list_prob = np.linspace(0,1) #p of depolarizing noise
    
    for p_dep in list_prob:
        # noisy channel
        model_process = (1 - p_dep) * choi_ideal + p_dep * choi_full_mixed

        row, cols = model_process.shape
        u1 = model_process.full().reshape((row * cols, 1))
        u2 = choi_experiment.reshape((row * cols, 1))
        distance_2 = scipy.linalg.norm(u1-u2)
        distances_list.append([p_dep, distance_2])

    distances_array = np.array(distances_list)
    arg_min = np.argmin(distances_array[:,1])
    print(distances_array[arg_min]) 

    # pdepol_min, min_distance
    # = [0.08163 0.31092]

    plt.plot(distances_array[:,0],distances_array[:,1], '-')
    plt.xlabel("depolarizing error p")
    plt.ylabel("cost function")    
#    plt.savefig("distances_p_depol.pdf")
    plt.show()
