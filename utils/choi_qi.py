import numpy as np
import qutip as qu
import sys
from itertools import product
import matplotlib.pyplot as plt
from matplotlib import colors


PLOT = True
SAVEFILE = False


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
    
    
dimHa = 2
dimHq = 3

GammaState = sum([qu.tensor(qu.basis(dimHa, ja), qu.basis(dimHq, jq), qu.basis(dimHa, ja), qu.basis(dimHq, jq)) 
                        for ja,jq in product(range(dimHa), range(dimHq))])/np.sqrt(dimHa * dimHq)
          
rhoGamma = GammaState * GammaState.dag()

for phi_tilde in [0.0, 1/2, 3/4]:
    phi =  phi_tilde * np.pi

    rot =  [ [np.cos(phi/2), 0, np.sin(phi/2)], 
              [0, 1, 0],  
              [-np.sin(phi/2), 0, np.cos(phi/2)]
            ]
    Rloss =  qu.tensor(qu.qeye(dimHa), qu.Qobj(rot))
    
    msGate = qu.tensor(qu.sigmax(), qu.Qobj([[0, 1, 0], [1, 0, 0], [0, 0, 0]])) + qu.tensor(qu.qeye(dimHa), proj(2, 2, 3)) 
    RxA = qu.tensor(qu.sigmax(), qu.qeye(dimHq)) # bit flip ancilla
    temp_Xq = proj(2, 2, dimHq) +  proj(1, 0, dimHq) +  proj(0, 1, dimHq)    # bit flip qutrit
    RxQ = qu.tensor(qu.qeye(dimHa), temp_Xq)
    U = RxA * RxQ * msGate * Rloss
    
    choiState = qu.tensor(qu.qeye(dimHa), qu.qeye(dimHq), U) * rhoGamma * qu.tensor(qu.qeye(dimHa), qu.qeye(dimHq), U.dag())

    A = choiState.full().real
    if PLOT:
        base_S = [(x,y) for x,y in product(range(2), range(3))]
        label_axis =[ "".join([str(_) for _ in bR])+ "," + "".join([str(_) for _ in bS]) for (bR, bS) in product(base_S, base_S) ]
    
        fig, ax = plt.subplots()

        img = ax.matshow(A, cmap='RdGy' ) #, interpolation='nearest', cmap=cmap, norm=norm)#, extent=[36,0,36,0])
#        img.set_clim(vmin=-1/np.sqrt(6), vmax=1/np.sqrt(6))
        img.set_clim(vmin=-1/6, vmax=1/6)        
        fig.colorbar(img)

        plt.title(f"$\phi={phi_tilde:1.2f}\pi$" + "  $p_\mathrm{loss}"+ f"={np.sin(phi_tilde*np.pi/2)**2:1.2f}$")
        ax.set_yticks(range(len(label_axis)))#, minor = True)
        ax.set_yticklabels(label_axis, fontsize = 8)#), minor = True)
        ax.set_xticks(range(len(label_axis)), [])
        ax.set_xticklabels([])

        ax=plt.gca()
        ax.set_xticks([x-0.5 for x in range(1,36)],minor=True )
        ax.set_yticks([y-0.5 for y in range(1,36)],minor=True)

        plt.grid(which="minor",ls="-",lw=0.2, color='w')
        ax.tick_params(which='minor', length=0, color='w')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('top')

        for j in range(6):
            ax.axhline(y=6*j - 0.5,color='lightgrey', linewidth=1, linestyle='-')
            ax.axvline(x=6*j - 0.5,color='lightgrey', linewidth=1, linestyle='-')
        plt.show()
#        plt.savefig(f"choiFinal_ideal_{phi_tilde}.pdf", bbox_inches='tight')

    if SAVEFILE:
        np.savetxt(f"choiFinal_ideal_{phi_tilde:1.2f}.dat", choiState, delimiter=',')   

