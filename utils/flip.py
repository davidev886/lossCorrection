import qutip as qu
import numpy as np
from numpy import linalg as LA
from itertools import product
_sigmas_P =   [[[1,  0],  [0, 1]], 
             [[0,  1],  [1, 0]], 
             [[0, -1j], [1j,0]], 
             [[1,0], [0,-1]]]
             


def bit_flip_channel(p):
    X = sigmax()
    I = qeye(2)
    return (1-p) * sprepost(I, I)  + p * sprepost(X, X) 
    

             
np.set_printoptions(precision=4,suppress=True)

def choiFlip(p):
    Lambda = [[1 - p, 0, 0, 1 - p], [0, p, p, 0], [0, p, p, 0], [1 - p, 0, 0, 1 - p]]
    return np.array(Lambda)
  
def normalize_operators(matrices):
    return [np.array(el) / LA.norm(el) for el in matrices]

  
def give_transformation_matrix():
    basis_elements_list = []
    for _sigma in _sigmas_P:
            # we use the column vectorization. that's why we transpose the basis_operator
            basis_operator = qu.Qobj(_sigma).full()
            basis_element = np.transpose(basis_operator)

            rows, cols = basis_element.shape
            dimH  = rows * cols 
            
            vector_basis = basis_element.reshape(1 , dimH)[0]
        
            # the paper defines the change of basis matrix not in the tandard way:
            # In linear algebra, the standard definition for the matrix P(e->v) for going to base e to basis v is
            # P(e->v)[x]_v = [x]_e where [x]_v is a vector x written in the basis v
            # and the matrix P(e->v) is given by putting by columns the components of the basis v written in the basis e 
            # By looking at appendix B formula 2.2 they define T(sigma -> w) as
            #  T(sigma -> w) |A>>_sigma = |A>>_w.
            #so we need to put the compoments of the basis w written in the basis sigma by row
            normalized_vector = np.conjugate(vector_basis) / LA.norm(vector_basis)
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
    
def apply_qnd_process_unit(choi, state_total):
    T_matrix = give_transformation_matrix()
    chi_matrix = get_chi_from_choi(choi, T_matrix)
    on_basis_Pauli = normalize_operators(_sigmas_P)
    rows, cols = chi_matrix.shape
    final_state_list = []
    for alpha, beta in product(range(rows), range(cols)):   
    
        OP_1 = qu.Qobj(on_basis_Pauli[alpha])
        OP_2 = qu.Qobj(on_basis_Pauli[beta]).dag()
        final_state_list.append(chi_matrix[alpha, beta] * OP_1 * state_total * OP_2)
        
    final_state = sum(final_state_list)
    return final_state


state_t = qu.basis(2,0) * qu.basis(2,0).dag()

choi = choiFlip(0.3)

final_state = apply_qnd_process_unit(choi, state_t)
print(final_state)


#function of qutip
print(to_choi(bit_flip_channel(0.3)))

chi_matrix = to_chi(to_choi(bit_flip_channel(0.3)))

print(chi_matrix.shape)


