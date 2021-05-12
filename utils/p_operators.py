import qutip as qu
import numpy as np

# ancilla always the last qubit
L = 7
dimQ = 2

Id = qu.tensor([qu.qeye(dimQ)] * L + [qu.qeye(2)])

temp = [[qu.qeye(dimQ)] * j +
        [qu.sigmax()] +
        [qu.qeye(dimQ)] * (L - j - 1) +
        [qu.qeye(2)]
        for j in range(L)
        ]
X = [qu.tensor(temp[j]) for j in range(L)]

temp = [[qu.qeye(dimQ)] * j +
        [qu.sigmay()] +
        [qu.qeye(dimQ)] * (L - j - 1) +
        [qu.qeye(2)]
        for j in range(L)
        ]
Y = [qu.tensor(temp[j]) for j in range(L)]

temp = [[qu.qeye(dimQ)] * j +
        [qu.sigmaz()] +
        [qu.qeye(dimQ)] * (L - j - 1) +
        [qu.qeye(2)]
        for j in range(L)
        ]
Z = [qu.tensor(temp[j]) for j in range(L)]

# ancilla operators
temp = [qu.qeye(dimQ)] * L + [qu.sigmax()]
Xa = qu.tensor(temp)

temp = [qu.qeye(dimQ)] * L + [qu.sigmaz()]
Za = qu.tensor(temp)

Pp_ancilla = (Id + Za) / 2
Pm_ancilla = (Id - Za) / 2

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

vacuum = qu.tensor([qu.basis(2, 0)] * L + [qu.basis(2, 0)])

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

# ancilla operators
temp = [qu.qeye(dimQ)] * L + [qu.sigmax()]
Xa = qu.tensor(temp)

temp = [qu.qeye(dimQ)] * L + [qu.sigmaz()]
Za = qu.tensor(temp)

Pp_ancilla = (Id + Za) / 2
Pm_ancilla = (Id - Za) / 2
