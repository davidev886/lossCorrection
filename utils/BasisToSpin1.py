import numpy as np
from qutip import *
import itertools

def BasisToSpin1(process):

    pauli_ops = [qeye(2).data.toarray(), # Qubit basis: Commonly known Pauli Basis
                 sigmax().data.toarray(),
                 sigmay().data.toarray(),
                 sigmaz().data.toarray()]

    gell_mann_ops = [qeye(3).data.toarray(), # Qutrit basis: https://en.wikipedia.org/wiki/Gell-Mann_matrices
                          get_qudit_from_qubit_operation(sigmax(), 3, (0, 1), True).data.toarray(),
                          get_qudit_from_qubit_operation(sigmay(), 3, (0, 1), True).data.toarray(),
                          get_qudit_from_qubit_operation(sigmaz(), 3, (0, 1), True).data.toarray(),
                          get_qudit_from_qubit_operation(sigmax(), 3, (0, 2), True).data.toarray(),
                          get_qudit_from_qubit_operation(sigmay(), 3, (0, 2), True).data.toarray(),
                          get_qudit_from_qubit_operation(sigmax(), 3, (1, 2), True).data.toarray(),
                          get_qudit_from_qubit_operation(sigmay(), 3, (1, 2), True).data.toarray(),
                          (get_qudit_from_qubit_operation(sigmaz(), 3, (0, 1), True).data.toarray()+
                          2*get_qudit_from_qubit_operation(sigmaz(), 3, (1, 2), True).data.toarray()) / np.sqrt(3)]


    TOp = []
    for i in pauli_ops:
        for j in gell_mann_ops:
            TOp.append(np.kron(i, j))

    dim = 2*3 # Dimension for qubit-qutrit
    print(TOp)
    exit()
    TOp = np.reshape(TOp, (dim**2, dim**2))
    TOp = np.conjugate(TOp)
    return dim**(-1)*(np.dot(TOp, np.dot(process, np.conjugate(np.transpose(TOp)))))



def get_qudit_from_qubit_operation(operation, dim = 4, transition = (0,2), is_hamiltonian = True):
    """
    Return a higher dimensional operation for a given qubit operation.

    Parameters
    ----------
    dim : int
        Dimension of the returned operation (dim=2 is qubit).

    operation : Qobj
        Qubit operation that is expanded to higher dimension.

    transition : tuple
        Tuple of states that defines on which transition of the qudit the operation acts.

    is_hamiltonian : bool
        Indicates if the input operation is a Hamiltonian (gets exponentiated later) or or time evolution operator.

    Returns
    -------
    qudit_operation : Qobj
        Higher dimensional operation object.
    """

    if is_hamiltonian:
        qudit_operation = np.zeros((dim, dim), dtype=complex)
    else:
        qudit_operation = np.identity(dim, dtype=complex)

    mask = np.zeros((dim, dim))
    mask[tuple(np.array(list(itertools.product(transition, repeat=2))).T)] = 1

    np.place(qudit_operation, mask, operation)
    qudit_operation = Qobj(qudit_operation)
    return qudit_operation


A = np.loadtxt("choiFinal_ideal.dat")/6
for _ in BasisToSpin1(A):
    print(_)

print( get_qudit_from_qubit_operation(sigmay(), 3, (1, 2), True).data.toarray(),)