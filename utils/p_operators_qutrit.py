import qutip as qu
import numpy as np
from numpy import linalg as LA
from itertools import product
import matplotlib.pyplot as plt

_lambdas_GM = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
               [[0, 1, 0], [1, 0, 0], [0, 0, 0]],
               [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]],
               [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
               [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
               [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]],
               [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
               [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]],
               [[1, 0, 0], [0, 1, 0], [0, 0, -2]]
               ]

_sigmas_P = [[[1, 0], [0, 1]],
             [[0, 1], [1, 0]],
             [[0, -1j], [1j, 0]],
             [[1, 0], [0, -1]]
             ]


# ancilla always the last qubit
L = 7

Id = qu.tensor([qu.qeye(3)] * L + [qu.qeye(2)])

x_qutrit = qu.Qobj([[0, 1, 0], [1, 0, 0], [0, 0, 1]])

y_qutrit = qu.Qobj([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]])

z_qutrit = qu.Qobj([[1, 0, 0], [0, -1, 0], [0, 0, 0]])

temp = [[qu.qeye(3)] * j +
        [x_qutrit] +
        [qu.qeye(3)] * (L - j - 1) +
        [qu.qeye(2)]
        for j in range(L)
        ]
X = [qu.tensor(temp[j]) for j in range(L)]

temp = [[qu.qeye(3)] * j +
        [y_qutrit] +
        [qu.qeye(3)] * (L - j - 1) +
        [qu.qeye(2)]
        for j in range(L)
        ]
Y = [qu.tensor(temp[j]) for j in range(L)]

temp = [[qu.qeye(3)] * j +
        [z_qutrit] +
        [qu.qeye(3)] * (L - j - 1) +
        [qu.qeye(2)]
        for j in range(L)
        ]
Z = [qu.tensor(temp[j]) for j in range(L)]

# ancilla operators
temp = [qu.qeye(3)] * L + [qu.sigmax()]
Xa = qu.tensor(temp)

temp = [qu.qeye(3)] * L + [qu.sigmaz()]
Za = qu.tensor(temp)

Pp_ancilla = (Id + Za) / 2
Pm_ancilla = (Id - Za) / 2


if L == 7:
    stab_qubits = [[0, 1, 2, 3],
                   [1, 2, 4, 5],
                   [2, 3, 5, 6]
                   ]

    Sx = [X[j1] * X[j2] * X[j3] * X[j4]
          for j1, j2, j3, j4 in stab_qubits
          ]
    Sz = [Z[j1] * Z[j2] * Z[j3] * Z[j4]
          for j1, j2, j3, j4 in stab_qubits
          ]

    ZL = Z[0] * Z[1] * Z[2] * Z[3] * Z[4] * Z[5] * Z[6]

    XL = X[0] * X[1] * X[2] * X[3] * X[4] * X[5] * X[6]

    Px = [(Id + el) / 2 for el in Sx]
    Pz = [(Id + el) / 2 for el in Sz]

    Pmx = [(Id - el) / 2 for el in Sx]
    Pmz = [(Id - el) / 2 for el in Sz]

    vacuum = qu.tensor([qu.basis(3, 0)] * L + [qu.basis(2, 0)])

    # logical states
    ZeroL = (Px[0] * Px[1] * Px[2] * vacuum).unit()
    OneL = (XL * ZeroL).unit()

    LogicalStates = [ZeroL,
                     OneL,
                     (ZeroL + OneL) / np.sqrt(2),
                     (ZeroL - OneL) / np.sqrt(2),
                     (ZeroL + 1j * OneL) / np.sqrt(2),
                     (ZeroL - 1j * OneL) / np.sqrt(2)
                     ]


def proj(ket, bra, dimH = 2):
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


def Rloss(initial_state, phi, qu_data):
    if initial_state.type == "ket":
        initial_state = initial_state * initial_state.dag()

    dimHq = 3  # Hilbert space data qutrit

    dimHa = 2  # Hilbert space ancilla qubit

    rloss = (proj(1, 1, dimHq) +
             np.cos(phi/2) * (proj(0, 0, dimHq) +
             proj(2 ,2, dimHq)) +
             np.sin(phi/2) * (proj(0, 2, dimHq) - proj(2, 0, dimHq))
             )

    rot_temp = [qu.qeye(dimHq)] * qu_data + [rloss] + [qu.qeye(dimHq)] * (L - qu_data - 1) + [qu.qeye(dimHa)]

    rot = qu.tensor(rot_temp)

    return rot * initial_state * rot.dag()


def Rloss_1(initial_state, phi, qu_data):
    if initial_state.type == "ket":
        initial_state = initial_state * initial_state.dag()

    dimHq = 3 # Hilbert space data qutrit

    dimHa = 2 # Hilbert space ancilla qubit

    rloss = (proj(0, 0, dimHq) + np.cos(phi/2) * (proj(1, 1, dimHq) + proj(2 ,2, dimHq))
                + np.sin(phi/2) * (proj(1, 2, dimHq) - proj(2, 1, dimHq)))

    rot_temp = [qu.qeye(dimHq)] * qu_data + [rloss] + [qu.qeye(dimHq)] * (L - qu_data - 1) + [qu.qeye(dimHa)]

    rot = qu.tensor(rot_temp)

    return rot * initial_state * rot.dag()


def Rloss_all_from_0(phi):
    # return  a list with 7 Rloss gates one for each data qutrit
    dimHq = 3 # Hilbert space data qutrit

    dimHa = 2 # Hilbert space ancilla qubit

    rloss = (proj(1, 1, dimHq) + np.cos(phi/2) * (proj(0, 0, dimHq) + proj(2 ,2, dimHq))
                + np.sin(phi/2) * (proj(0, 2, dimHq) - proj(2, 0, dimHq)) )

    temp = [[qu.qeye(dimHq)] * j + [rloss] + [qu.qeye(dimHq)] * (L - j - 1) + [qu.qeye(dimHa)] for j in range(L)]
    R_loss_list = [qu.tensor(temp[j]) for j in range(L)]

    return R_loss_list


def Rloss_all_from_1(phi):
    # return  a list with 7 Rloss gates one for each data qutrit
    dimHq = 3 # Hilbert space data qutrit

    dimHa = 2 # Hilbert space ancilla qubit

    rloss = (proj(0, 0, dimHq) + np.cos(phi/2) * (proj(1, 1, dimHq) + proj(2 ,2, dimHq))
                + np.sin(phi/2) * (proj(1, 2, dimHq) - proj(2, 1, dimHq)) )

    temp = [[qu.qeye(dimHq)] * j + [rloss] + [qu.qeye(dimHq)] * (L - j - 1) + [qu.qeye(dimHa)] for j in range(L)]
    R_loss_list = [qu.tensor(temp[j]) for j in range(L)]

    return R_loss_list


def Rloss_all(phi):
    # return  a list with 7 Rloss gates one for each data qutrit
    dimHq = 3 # Hilbert space data qutrit

    dimHa = 2 # Hilbert space ancilla qubit

    rloss = (proj(1, 1, dimHq) + np.cos(phi/2) * (proj(0, 0, dimHq) + proj(2 ,2, dimHq))
                + np.sin(phi/2) * (proj(0, 2, dimHq) - proj(2, 0, dimHq)) )

    temp = [[qu.qeye(dimHq)] * j + [rloss] + [qu.qeye(dimHq)] * (L - j - 1) + [qu.qeye(dimHa)] for j in range(L)]
    R_loss_list = [qu.tensor(temp[j]) for j in range(L)]

    return R_loss_list


def normalize_operators(matrices_list):
    return [qu.Qobj(np.array(el) / LA.norm(el)) for el in matrices_list]


def give_transformation_matrix():
    basis_elements_list = []
    # the order of these for loops is important to define the T_matrix because the convention is
    # sigma_j x lambda_k in the choi matrix we get from experiments
    for j in range(len(_sigmas_P)):
        for i in range(len(_lambdas_GM)):
            _lambda = _lambdas_GM[i]
            _sigma = _sigmas_P[j]
            # we use the column vectorization convention. that's why we transpose the basis_element
            basis_operator = (qu.tensor(qu.Qobj(_sigma), qu.Qobj(_lambda))).full()

            basis_element = np.transpose(basis_operator)
            rows, cols = basis_element.shape
            dimH  = rows * cols

            vector_basis = basis_element.reshape(1 , dimH)[0]

            # the paper 1111.6950 defines the change of basis matrix not in the standard way:
            # In linear algebra, the standard definition for the matrix P(e->v) for going to base e to basis v is
            # P(e->v)[x]_v = [x]_e where [x]_v is a vector x written in the basis v
            # and the matrix P(e->v) is given by putting the components of the basis v written in the basis e by columns
            # By looking at appendix B formula 2.2 they define T(sigma -> w) as
            # T(sigma -> w) |A>>_sigma = |A>>_w.
            # so we need to put the compoments of the basis w written in the basis sigma by row
            normalized_vector = np.conjugate(vector_basis / LA.norm(vector_basis))

            basis_elements_list.append(normalized_vector.tolist())

    T_matrix = np.array(basis_elements_list)
    return T_matrix

def orthonormal_basis_operator(T_matrix):
    on_basis = []
    for el in T_matrix:
        dimH = int(np.sqrt(len(el)))
        # we need to transpose as we were working with col-vec convention
        on_basis.append(el.reshape((dimH, dimH)).T)
    return on_basis


def get_chi_from_choi(choi, T_matrix):
    return np.dot(T_matrix, np.dot(choi, T_matrix.conj().T))

def apply_qnd_process_unit(chi_matrix, state_total, qu_data, chi_threshold):

    if state_total.type == "ket":
        state_total = state_total * state_total.dag()

    on_basis_lambda = normalize_operators(_lambdas_GM)
    on_basis_Pauli = normalize_operators(_sigmas_P)

    rows, cols = chi_matrix.shape
    final_state_list = []

    for alpha, beta in product(range(rows), range(cols)):
        if chi_matrix[alpha, beta]:
            a_GM = alpha % 9
            a_Pauli = (alpha - a_GM) // 9

            OP_temp = [qu.qeye(3)] * qu_data + [on_basis_lambda[a_GM]] + [qu.qeye(3)] * (L - qu_data - 1) + [on_basis_Pauli[a_Pauli]]
            OP_1 = qu.tensor(OP_temp)

            b_GM = beta % 9
            b_Pauli = (beta - b_GM) // 9

            OP_temp = [qu.qeye(3)] * qu_data + [on_basis_lambda[b_GM]] + [qu.qeye(3)] * (L - qu_data - 1) + [on_basis_Pauli[b_Pauli]]
            OP_2 = qu.tensor(OP_temp)

            #partial_state = chi_matrix[alpha, beta] * OP_1 * state_total * OP_2.dag()
            #final_state_list.append(partial_state)

            if alpha > beta:
                action = chi_matrix[alpha, beta] * OP_1 * state_total * OP_2.dag()
                final_state_list.append(action + action.dag())
            elif alpha == beta:
                final_state_list.append(chi_matrix[alpha, beta] * OP_1 * state_total * OP_2.dag())
    final_state = sum(final_state_list)
    return final_state


def find_false_positive():

    phi_tilde = 0
    rotation_ops_0 = Rloss_all_from_0(phi_tilde * np.pi)

    rotation_ops_1 = Rloss_all_from_1(phi_tilde * np.pi)

    on_basis_lambda = normalize_operators(_lambdas_GM)
    on_basis_Pauli = normalize_operators(_sigmas_P)
    choi_ideal = np.loadtxt("choiFinal_ideal.dat")
    choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')


    T_matrix = give_transformation_matrix()
    chi_matrix = get_chi_from_choi(6 * choi_experiment, T_matrix) #.round(15)
    a = 1/np.sqrt(2)
    b = 1/np.sqrt(2)
    a = 1
    b = 0
    c = 0

    if 1:
#    for a, b in [[1,0], [0,1], [1/np.sqrt(2), 1/np.sqrt(2)]]:
        state_qutrit = (a * qu.basis(3,0) + b * qu.basis(3,1) + c * qu.basis(3,2)).unit()

        state_0 = (qu.tensor([state_qutrit, qu.basis(2,0)])).unit()
        rho_L = state_0 * state_0.dag()

        rho_L = rotation_ops_0[0] * rho_L * rotation_ops_0[0].dag()
        print(rho_L)
        exit()
        rho_L = apply_qnd_process_unit(chi_matrix, rho_L, 0, 0)

        print(f"prob ancilla 0 {(rho_L * Pp_ancilla).tr()}")
        print(f"prob ancilla 1 {(rho_L * Pm_ancilla).tr()}")
        p_10 = (rho_L * Pp_ancilla).tr()
        p_11 = (rho_L * Pm_ancilla).tr()
        print("--------")
        rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / p_10
        print(rho_L)
        w_0 = rho_L.ptrace([0])
        print(w_0)
#         rho_L = apply_qnd_process_unit(chi_matrix, rho_L, 0, 0)
#         print(f"prob ancilla 0 {(rho_L * Pm_ancilla).tr()}")
#         print(f"prob ancilla 1 {(rho_L * Pp_ancilla).tr()}")
#         p_20 = (rho_L * Pp_ancilla).tr()
#         p_21 = (rho_L * Pm_ancilla).tr()
#         print(p_10*p_20)
#         print(p_11*p_21)
#         print()
#         print(p_10*p_20 + p_11*p_21)
    exit()
    prob_outcome_1 = (rho_L * Pp_ancilla).tr()
    print(prob_outcome_1)
    rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / abs(prob_outcome_1)
    print(rho_L)
    print("rotation from 1")
    rho_L = rotation_ops_1[0] * rho_L * rotation_ops_1[0].dag()
    rho_L = apply_qnd_process_unit(chi_matrix, rho_L, 0, 0)
    print("apply second qnd")
    print(rho_L)
    prob_outcome_2 = (rho_L * Pp_ancilla).tr()
    rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / abs(prob_outcome_2)
    print("project on 0 second time")
    print(Pp_ancilla)
    print(rho_L)


    print(state_0 * state_0.dag())


if __name__ == "__main__":


    np.set_printoptions(precision = 5, suppress = True,  linewidth=100000)
#     find_false_positive()
#     exit()

    phi_tilde = 1 / 2
    rotation_ops_0 = Rloss_all_from_0(phi_tilde * np.pi)

    rotation_ops_1 = Rloss_all_from_1(phi_tilde * np.pi)

    on_basis_lambda = normalize_operators(_lambdas_GM)
    on_basis_Pauli = normalize_operators(_sigmas_P)
    choi_ideal = np.loadtxt("choiFinal_ideal.dat")
    choi_experiment = 6*  np.real(np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=','))

    print(choi_experiment.round(3))
    exit()
    T_matrix = give_transformation_matrix()
    chi_matrix = get_chi_from_choi(choi_experiment.round(4), T_matrix).round(15)
    a = 1/np.sqrt(2)
    b = 1/np.sqrt(2)
#    a = 1
#    b = 0
#    state_qutrit = (a * qu.basis(3,0) + b * qu.basis(3,1)).unit()

#    state_0 = (qu.tensor([state_qutrit, qu.basis(2,0)])).unit()
#    state_0 = (qu.tensor([ZeroL, qu.basis(2,0)])).unit()
    state_0 = OneL
    rho_L = state_0 * state_0.dag()
    print(rho_L)
    rho_L = rotation_ops_0[0] * rho_L * rotation_ops_0[0].dag()

    rho_L = apply_qnd_process_unit(chi_matrix, rho_L, 0, 0)

    print("----")
    print(rho_L)
    print("----")

    prob_outcome_1 = (rho_L * Pp_ancilla).tr()
    print(prob_outcome_1)
    rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / abs(prob_outcome_1)
    print(rho_L)
    print("rotation from 1")
    rho_L = rotation_ops_1[0] * rho_L * rotation_ops_1[0].dag()
    rho_L = apply_qnd_process_unit(chi_matrix, rho_L, 0, 0)
    print("apply second qnd")
    print(rho_L)
    prob_outcome_2 = (rho_L * Pp_ancilla).tr()
    rho_L = Pp_ancilla * rho_L * Pp_ancilla.dag() / abs(prob_outcome_2)
    print("project on 0 second time")
    print(rho_L)
    print(rho_L == OneL * OneL.dag())
