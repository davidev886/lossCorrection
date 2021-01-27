import qutip as qu
import numpy as np
import scipy.linalg
from scipy.linalg import eigh, svd
from p_operators_qutrit import get_chi_from_choi, give_transformation_matrix, normalize_operators
from itertools import product
import matplotlib.pyplot as plt

_lambdas_GM = [[[ 1, 0, 0], [ 0, 1, 0], [0, 0, 1]],
             [[ 0 , 1 , 0], [1 , 0 , 0], [0 , 0 , 0 ]],
             [[ 0 , -1j , 0], [1j , 0 , 0], [0 , 0 , 0 ] ],
             [[ 1 , 0 , 0], [0 , -1 , 0], [0 , 0 , 0 ] ],
             [[ 0 , 0 , 1], [0 , 0 , 0], [1 , 0 , 0 ] ],
             [[ 0 , 0 , -1j], [0 , 0 , 0], [1j , 0 , 0 ] ],
             [[ 0 , 0 , 0], [0 , 0 , 1], [0 , 1 , 0 ] ],
             [[ 0 , 0 , 0], [0 , 0 , -1j], [0 , 1j , 0 ] ],
             [[ 1 , 0 , 0], [0 , 1 , 0], [0 , 0 , -2 ]]]

_sigmas_P =   [[[1,  0],  [0, 1]],
             [[0,  1],  [1, 0]],
             [[0, -1j], [1j,0]],
             [[1,0], [0,-1]]]




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

def experimental_qnd_channel(chi_matrix):
    on_basis_lambda = normalize_operators(_lambdas_GM)
    on_basis_Pauli = normalize_operators(_sigmas_P)

    rows, cols = chi_matrix.shape
    final_state_list = []

    for alpha, beta in product(range(rows), range(cols)):
        if chi_matrix[alpha, beta]:
            a_GM = alpha % 9
            a_Pauli = (alpha - a_GM) // 9

            OP_temp = [on_basis_Pauli[a_Pauli], on_basis_lambda[a_GM]]
            OP_1 = qu.tensor(OP_temp)

            b_GM = beta % 9
            b_Pauli = (beta - b_GM) // 9

            OP_temp = [on_basis_Pauli[b_Pauli], on_basis_lambda[b_GM]]
            OP_2 = qu.tensor(OP_temp)

            action = chi_matrix[alpha, beta] * qu.spre(OP_1) *  qu.spost(OP_2.dag())
            final_state_list.append(action)

    QND_exp_channel = sum(final_state_list)
    return QND_exp_channel

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

    omega = sum([qu.tensor(qu.basis(2, i), qu.basis(3, j),
                 qu.basis(2, i), qu.basis(3, j))
                for i, j in product(range(2), range(3))])

    omega_state = omega * omega.dag()

    dimHtot = 2 * 3

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
    fidelities_list = []
    list_prob = np.linspace(0,1) #p of depolarizing noise
    
    for p_dep in list_prob:
        channel_dep = depolarization_channel(p_dep)
#construct the choi matrix of the ideal QND followed by the depolarization channel
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

        choi_depolarizing = (p_dep * qu.tensor(qu.qeye(2), qu.qeye(3), qu.qeye(2) , qu.qeye(3))/dimHtot
                            + (1 - p_dep) * omega_state)
        choi_experiment_op = qu.Qobj(inpt=choi_experiment, dims=choi_depolarizing.dims, type=choi_depolarizing.type)
        print(qu.process_fidelity(choi_depolarizing, choi_experiment_op))
        fidelities_list.append([p_dep, qu.process_fidelity(choi_depolarizing, choi_experiment_op)])

    distances_array = np.array(distances_list)
    fidelities_array =  np.array(fidelities_list)
    arg_min = np.argmin(distances_array[:,1])
    print(distances_array[arg_min]) 
    print("-----", distances_array[arg_min][0])
    choi_experiment = 6 * np.real(np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=','))
    T_matrix = give_transformation_matrix()
    chi_matrix = np.real(get_chi_from_choi(choi_experiment, T_matrix))

    exp_QND = experimental_qnd_channel(chi_matrix)
    channel_dep_close = depolarization_channel(distances_array[arg_min][0])


    p_dep_min = distances_array[arg_min][0]

    plt.plot(distances_array[:,0],distances_array[:,1], '-', label="Cost function $C(p)$")
    plt.plot(distances_array[:,0],fidelities_array[:,1], '-', label="Process fidelity")
    plt.xlabel("depolarizing error p")
#    plt.ylabel("cost function")
    plt.title(f"depolarizing error p min = {distances_array[arg_min][0]:.3}")
    ax = plt.gca()
    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
#    ax.tick_params(axis='x', which='minor', bottom=False)
    plt.plot([distances_array[arg_min][0]], [distances_array[arg_min][1]], 'o')
    plt.legend()
    plt.savefig("distances_p_depol.pdf")
    plt.show()
