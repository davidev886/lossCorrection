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

def depolarization_channel(p_dep):
# Define Omega state for computing Choi matrix 
# omega = [00,01,02,10,11,12] * [00,01,02,10,11,12] ancilla first
    omega = sum([qu.tensor(qu.basis(2, i), qu.basis(3, j), 
                 qu.basis(2, i), qu.basis(3, j)) 
                for i, j in product(range(2), range(3))])

    omega_state = omega * omega.dag()

    dimHtot = 2 * 3
    
    choi_depolarizing = (p_dep * qu.tensor(qu.qeye(2), qu.qeye(3), qu.qeye(2) , qu.qeye(3))/dimHtot 
                            + (1 - p_dep) * omega_state)

    vals, vecs = eigh(choi_depolarizing.data.todense())

    kraus_matrix = [np.sqrt(vals[j]) * vecs[:, j].reshape((dimHtot, dimHtot)).T for j in range(len(vecs))]


    system_length = len(choi_depolarizing.dims[0])//2
    dim_kraus = [choi_depolarizing.dims[0][-system_length:], choi_depolarizing.dims[0][-system_length:]]   
    kraus_ops = [qu.Qobj(inpt=op, dims=dim_kraus) for op in kraus_matrix]

    channel_dep = sum([qu.spre(op) * qu.spost(op.dag()) for op in kraus_ops])
    return channel_dep
    

if __name__ == "__main__":
    np.set_printoptions(precision = 4, suppress = True, linewidth=100000)

    choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')

#Define the ideal channel with A0 and A1 operators:
#    A0 = |1><1|+|0><0|
    A0 = proj(1, 1, dimH = 3) + proj(0, 0, dimH = 3)

#    A1 = |2><2|
    A1 = proj(2, 2, dimH = 3)

    id_a = qu.qeye(2)
    X_a = qu.sigmax()
    id_all = qu.tensor(qu.qeye(2), qu.qeye(3))

    QND_ideal = qu.tensor(id_a, A0) + qu.tensor(X_a, A1)

#   QND_ideal = 
#     [[1. 0. 0. 0. 0. 0.]
#      [0. 1. 0. 0. 0. 0.]
#      [0. 0. 0. 0. 0. 1.]
#      [0. 0. 0. 1. 0. 0.]
#      [0. 0. 0. 0. 1. 0.]
#      [0. 0. 1. 0. 0. 0.]]

    choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')

    distances_list = []
    list_prob = np.linspace(0,1) #p of depolarizing noise
    
    for p_dep in list_prob:
        channel_dep = depolarization_channel(p_dep)
# construct the choi matrix of the ideal QND followed by the depolarizing channel       
# choi_ideal_and_depol = sum |i_1 j_1, i_2 j_2> <i_1 j_1, i_2 j_2| x Dep(QND(|i_1 j_1, i_2 j_2> <i_1 j_1, i_2 j_2|))
# i_1, i_2 qubit basis  j_1, j_2 qutrit basis

        choi_parts_ideal_depol = []        
        for i1, j1 in product(range(2), range(3)):
            for i2, j2 in product(range(2), range(3)):    
                state_ket = qu.tensor(qu.basis(2, i1), qu.basis(3, j1))
                state_bra = qu.tensor(qu.basis(2, i2), qu.basis(3, j2)).dag()           
                Qij = QND_ideal * state_ket * state_bra * QND_ideal.dag()
                tenso_ij_depol_Qij = qu.tensor(state_ket * state_bra, channel_dep(Qij))
                choi_parts_ideal_depol.append(tenso_ij_depol_Qij)
        choi_ideal_and_depol = sum(choi_parts_ideal_depol) / 6.

        row, cols = choi_ideal_and_depol.shape
        u1 = choi_ideal_and_depol.full().reshape((row * cols, 1))
        u2 = choi_experiment.reshape((row * cols, 1))
        distance_2 = scipy.linalg.norm(u1-u2)
        distances_list.append([p_dep, distance_2])

    distances_array = np.array(distances_list)
    arg_min = np.argmin(distances_array[:,1])
    print(distances_array[arg_min]) 


    plt.plot(distances_array[:,0],distances_array[:,1], '-')
    plt.xlabel("depolarizing error p")
    plt.ylabel("cost function")    
    plt.title(f"depolarizing error p min = {distances_array[arg_min][0]:.3}")

    from matplotlib.ticker import AutoMinorLocator
    ax = plt.gca()
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    plt.plot([distances_array[arg_min][0]], [distances_array[arg_min][1]], 'o')
#    plt.savefig("distances_p_depol.pdf")
    plt.show()
