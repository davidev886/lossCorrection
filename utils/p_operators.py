import qutip as qu

#ancilla always the last qubit
L = 7

Id = qu.tensor([qu.qeye(2)] * L)

temp = [[qu.qeye(2)] * j + [qu.sigmax()] + [qu.qeye(2)] * (L - j -1) for j in range(L)]
X = [qu.tensor(temp[j]) for j in range(L)]

temp = [[qu.qeye(2)] * j + [qu.sigmay()] + [qu.qeye(2)] * (L - j -1) for j in range(L)]
Y = [qu.tensor(temp[j]) for j in range(L)]

temp = [[qu.qeye(2)] * j + [qu.sigmaz()] + [qu.qeye(2)] * (L - j -1) for j in range(L)]

Z = [qu.tensor(temp[j]) for j in range(L)]

stab_qubits = [[0,1,2,3], [1,2,4,5], [2,3,5,6]]

Sx = [X[j1] * X[j2] * X[j3] * X[j4] for j1,j2,j3,j4 in stab_qubits]
Sz = [Z[j1] * Z[j2] * Z[j3] * Z[j4] for j1,j2,j3,j4 in stab_qubits]

ZL = Z[0] * Z[1] * Z[2] * Z[3] * Z[4] * Z[5] * Z[6]

XL = X[0] * X[1] * X[2] * X[3] * X[4] * X[5] * X[6]

Px = [(Id + el) / 2 for el in Sx]
Pz = [(Id + el) / 2 for el in Sz]

Pmx = [(Id - el) / 2 for el in Sx]
Pmz = [(Id - el) / 2 for el in Sz]

vacuum = qu.tensor([qu.basis(2,0)] * L)

#logical states
ZeroL = (Px[0] * Px[1] * Px[2] * vacuum).unit()
OneL = (XL * ZeroL).unit()


# def get_qnd_error_op(qnd_err_str, data_qubit, ancilla_qubit = 7):
#     errors_p = "IXYZ"
#     err1_str, err2_str = qnd_err_str
#     errs_q = [Id, X[data_qubit], Y[data_qubit], Z[data_qubit]]
#     errs_a = [Id, X[ancilla_qubit], Y[ancilla_qubit], Z[ancilla_qubit]]
#     errs_q_t = ["Id", f"X[{data_qubit}]",f"Y[{data_qubit}]",f"Z[{data_qubit}]"]
#     errs_a_t = ["Id", f"X[{ancilla_qubit}]",f"Y[{ancilla_qubit}]",f"Z[{ancilla_qubit}]"]
# 
#     op1_index = errors_p.index(err1_str)
#     op2_index = errors_p.index(err2_str)
#     
#     return errs_q[op1_index] * errs_a[op2_index]
