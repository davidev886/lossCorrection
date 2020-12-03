"""

`choi_repr` contains several functions for computing and managing the 
 the choi operator of a qubit-qutrit system

"""

import numpy as np
import qutip as qu
import sys
from itertools import product
from p_operators_qutrit import proj



def apply_choi(rho, choiState):
    dimTot = stateTry.shape[0]
    if rho.type == "ket":
        rho = (rho * rho.dag()).unit()


    finalState = sum([rho[i,j] * choiState[dimTot*i : dimTot*i + dimTot, dimTot*j: dimTot*j + dimTot] 
                   for i,j in product(range(dimTot), range(dimTot))])    
    
    finalState = qu.Qobj(inpt=finalState, dims=rho.dims)
    return finalState

# for plotting experimental choi operators
if 0: #__name__ == "__main__":
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from itertools import product
    base_S = [(x,y) for x,y in product(range(2), range(3))]
    label_axis =[ "".join([str(_) for _ in bR])+ "," + "".join([str(_) for _ in bS]) for (bR, bS) in product(base_S, base_S) ]

    choi_experiment = np.genfromtxt("qubitqutrit_choi_noloss.csv", dtype=complex, delimiter=',')

    np.set_printoptions(precision=4, suppress=True, threshold=sys.maxsize)

    ratio = np.abs(choi_experiment.imag)/np.abs(choi_experiment.real)
    #compute and print ratio imaginary / real part 
    for i,j in product(range(ratio.shape[0]), range(ratio.shape[1])) :
        print(f"{i: 3d},{j: 3d}, {abs(choi_experiment.real[i,j]):1.2e},  {abs(choi_experiment.imag[i,j]):1.2e}, {ratio[i,j]:6.2f}")


    for index, A in enumerate([choi_experiment.real, choi_experiment.imag, np.abs(choi_experiment.imag)/np.abs(choi_experiment.real) ]):
        str_img = ["real", "imag", "ratio"][index]
        fig, ax = plt.subplots()

        img = ax.matshow(A, cmap='RdGy' ) #, interpolation='nearest', cmap=cmap, norm=norm)#, extent=[36,0,36,0])
        if index < 2:
            img.set_clim(vmin=-0.15, vmax=0.15)
        fig.colorbar(img)

        plt.title(f"choi experiment no loss {str_img}")
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
        #ax.yaxis.grid(True, which='major')
    #        plt.show()



        plt.savefig(f"choi_exp_noloss_{str_img}.pdf", bbox_inches='tight')
    

# for plotting ideal choi operators
if __name__ == "__main__":
    dimHa = 2
    dimHq = 3
    dimTot = dimHa * dimHq
    GammaState = sum([qu.tensor(qu.basis(dimHa, ja), qu.basis(dimHq, jq), qu.basis(dimHa, ja), qu.basis(dimHq, jq)) 
                            for ja,jq in product(range(dimHa), range(dimHq))])
              
    rhoGamma = GammaState * GammaState.dag()        
    phi_tilde = 1 / 2 ## loss 50%
#    phi_tilde = 3 / 4. ## loss 85%  


    for phi_tilde in [0.0, 1/2, 3/4]:
        phi =  phi_tilde * np.pi
        A0 = proj(1, 1, dimHq) + np.cos(phi/2) * proj(0, 0, dimHq) + np.sin(phi/2) * proj(0, 2, dimHq)
        A1 =  - np.sin(phi/2) * proj(2, 0, dimHq) + np.cos(phi/2) * proj(2, 2, dimHq)

#        U = qu.tensor(qu.qeye(dimHa), A0) + qu.tensor(qu.sigmax(), A1)
        rot =  [ [np.cos(phi/2), 0, np.sin(phi/2)], 
                  [0, 1, 0],  
                  [-np.sin(phi/2), 0, np.cos(phi/2)]
                ]
        Rloss =  qu.tensor(qu.qeye(dimHa), qu.Qobj(rot))
        
        msGate = qu.tensor(qu.qeye(dimHa), proj(2, 2, 3)) +  qu.tensor(qu.sigmax(), qu.Qobj([[0, 1, 0], [1, 0, 0], [0, 0, 0]]))
        RxA = qu.tensor(qu.sigmax(), qu.qeye(dimHq)) # bit flip ancilla
        RxQ = qu.tensor(qu.qeye(dimHa), qu.Qobj([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) + qu.Qobj([[0, 0, 0], [0, 0, 0], [0, 0, 1]]))  # bit flip qutrit       
        print(RxA.dims)
        print(RxQ.dims)        
        print(msGate.dims)                
        print(Rloss.dims)
        print(RxQ)
        U = RxA * RxQ * msGate * Rloss
        
        choiState = qu.tensor(qu.qeye(dimHa), qu.qeye(dimHq), U) * rhoGamma * qu.tensor(qu.qeye(dimHa), qu.qeye(dimHq), U.dag())

        import matplotlib.pyplot as plt
        from matplotlib import colors
        from itertools import product

        base_S = [(x,y) for x,y in product(range(2), range(3))]
        label_axis =[ "".join([str(_) for _ in bR])+ "," + "".join([str(_) for _ in bS]) for (bR, bS) in product(base_S, base_S) ]
        CC = choiState.full() / 6


        for A in [choiState.full().real / 6]: #, choiState.full().imag / 6] :   
            fig, ax = plt.subplots()

            img = ax.matshow(A, cmap='RdGy' ) #, interpolation='nearest', cmap=cmap, norm=norm)#, extent=[36,0,36,0])
            img.set_clim(vmin=-0.15, vmax=0.15)
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
            #ax.yaxis.grid(True, which='major')
            plt.show()
    
    
#         np.savetxt(f"choiFinal_ideal_{phi_tilde:1.2f}.dat", CC, delimiter=',')   
    #    plt.savefig(f"choiFinal_ideal_{phi_tilde}.pdf", bbox_inches='tight')
