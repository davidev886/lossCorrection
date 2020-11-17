import qutip as qu
import numpy as np
from numpy import linalg as LA
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
             
    
#ancilla always the last qubit
L = 1
ancilla_num = 1

Id = qu.tensor([qu.qeye(3)] * L + [qu.qeye(2)])

x_qutrit = qu.Qobj([[0,1,0] , [1,0,0], [0,0,1]])

y_qutrit = qu.Qobj([[0,-1j,0] , [1j,0,0], [0,0,1]])

z_qutrit = qu.Qobj([[1,0,0] , [0,-1,0], [0,0,1]])

temp = [[qu.qeye(3)] * j + [x_qutrit] + [qu.qeye(3)] * (L - j - 1) + [qu.qeye(2)] for j in range(L)]
X = [qu.tensor(temp[j]) for j in range(L)]

temp = [[qu.qeye(3)] * j + [y_qutrit] + [qu.qeye(3)] * (L - j - 1) + [qu.qeye(2)] for j in range(L)]
Y = [qu.tensor(temp[j]) for j in range(L)]

temp = [[qu.qeye(3)] * j + [z_qutrit] + [qu.qeye(3)] * (L - j - 1) + [qu.qeye(2)] for j in range(L)]
Z = [qu.tensor(temp[j]) for j in range(L)]

if L == 7:
    stab_qubits = [[0,1,2,3], [1,2,4,5], [2,3,5,6]]

    Sx = [X[j1] * X[j2] * X[j3] * X[j4] for j1,j2,j3,j4 in stab_qubits]
    Sz = [Z[j1] * Z[j2] * Z[j3] * Z[j4] for j1,j2,j3,j4 in stab_qubits]

    ZL = Z[0] * Z[1] * Z[2] * Z[3] * Z[4] * Z[5] * Z[6]

    XL = X[0] * X[1] * X[2] * X[3] * X[4] * X[5] * X[6]

    Px = [(Id + el) / 2 for el in Sx]
    Pz = [(Id + el) / 2 for el in Sz]

    Pmx = [(Id - el) / 2 for el in Sx]
    Pmz = [(Id - el) / 2 for el in Sz]

    vacuum = qu.tensor([qu.basis(3,0)] * L + [qu.basis(2,0)])

    #logical states
    ZeroL = (Px[0] * Px[1] * Px[2] * vacuum).unit()
    OneL = (XL * ZeroL).unit()


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

def Rloss(phi):
    # return  a list with 7 Rloss gates one for each data qutrit
    dimHq = 3 # Hilbert space data qutrit
    
    dimHa = 2 # Hilbert space ancilla qubit
    
    rloss = (proj(1, 1, dimHq) + np.cos(phi/2) * (proj(0, 0, dimHq) + proj(2 ,2, dimHq)) 
                + np.sin(phi/2) * (proj(0, 2, dimHq) - proj(2, 0, dimHq)) )
    
    temp = [[qu.qeye(dimHq)] * j + [rloss] + [qu.qeye(dimHq)] * (L - j - 1) + [qu.qeye(dimHa)] for j in range(L)]
    R_loss_list = [qu.tensor(temp[j]) for j in range(L)]    
    
    return R_loss_list

def normalize_operators(matrices):
    return [qu.Qobj(np.array(el) / LA.norm(el)) for el in matrices]
             
def give_transformation_matrix():
    basis_elements_list = []
    #the order of these for loops is important to define the T_matrix
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

def apply_qnd_process_unit(choi, state_total, qu_data):

    if state_total.type == "ket":
        state_total = state_total * state_total.dag()

    T_matrix = give_transformation_matrix()
    np.savetxt("T_matrix.dat", T_matrix, fmt="%1.4f")

    chi_matrix = get_chi_from_choi(choi, T_matrix)
    np.savetxt("chi_matrix.dat", chi_matrix.real, fmt="%1.5f")
    op_label = [[str(_) for _ in range(36)]] * 2 
    fig, ax = qu.qpt_plot_combined(chi_matrix / 6,op_label)
    plt.show()
    on_basis_lambda = normalize_operators(_lambdas_GM)
    on_basis_Pauli = normalize_operators(_sigmas_P)

    rows, cols = chi_matrix.shape
    final_state_list = []

    for alpha, beta in product(range(rows), range(cols)):   
        a_GM = alpha % 9
        a_Pauli = (alpha - a_GM) // 9
    
        OP_temp = [qu.qeye(3)] * qu_data + [on_basis_lambda[a_GM]] + [qu.qeye(3)] * (L - qu_data - 1) + [on_basis_Pauli[a_Pauli]]
        OP_1 = qu.tensor(OP_temp)

        b_GM = beta % 9
        b_Pauli = (beta - b_GM) // 9

        OP_temp = [qu.qeye(3)] * qu_data + [on_basis_lambda[b_GM]] + [qu.qeye(3)] * (L - qu_data - 1) + [on_basis_Pauli[b_Pauli]]
        OP_2 = qu.tensor(OP_temp)

        final_state_list.append(chi_matrix[alpha, beta] * OP_1 * state_total * OP_2.dag())
        
    final_state = sum(final_state_list)
    return final_state        
        
        
if __name__ == "__main__":
#    on_basis_lambda = normalize_operators(_lambdas_GM)
#    on_basis_Pauli = normalize_operators(_sigmas_P)
    choi_ideal = np.loadtxt("choiFinal_ideal.dat")
    print("qui")
    a = np.random.random()  + np.random.random() * 1j
    b = np.random.random()  + np.random.random() * 1j
    state_qutrit = (a * qu.basis(3,0) + b * qu.basis(3,1)).unit()
    state_total = (qu.tensor([state_qutrit] + [qu.basis(2,0)])).unit()

    final = apply_qnd_process_unit(choi_ideal, state_total, qu_data = 0)
    np.savetxt("final.dat", final.full())
    np.savetxt("initial.dat", state_total * state_total.dag().full())
    