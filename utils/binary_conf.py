import numpy as np

from itertools import product
    
class binary_configurations(object):
    def __init__(self, n = 7):
        self.n = n
        configurations = {}
        for i in range(2**n):
            configuration_str = bin(i)[2:].zfill(n)
            configuration_int = [int(_) for _ in configuration_str]
            num_particles = sum(configuration_int)
            if num_particles in configurations:
                configurations[num_particles].append(configuration_int)
            else:
                configurations[num_particles]= [configuration_int]
        self.configurations = configurations
    
    def generate_configuration_loss_qnderror(self, num_l, num_q):
        configuration_loss = []
        configuration_qnd = []
        configuration_loss = self.configurations[num_l]
        configuration_qnd = self.configurations[num_q]   
        
        return product(configuration_loss, configuration_qnd)     


class binary_configurations_loss_qnd_stab(object):
    def __init__(self, n = 7, stab = 3 ):
        self.n = n
        configurations = {}
        for i in range(2**n):
            configuration_str = bin(i)[2:].zfill(n)
            configuration_int = [int(_) for _ in configuration_str]
            num_particles = sum(configuration_int)
            if num_particles in configurations:
                configurations[num_particles].append(configuration_int)
            else:
                configurations[num_particles]= [configuration_int]
        self.configurations = configurations

        configurations_stab = {}
        for i in range(2**stab):
            configuration_str = bin(i)[2:].zfill(stab)
            configuration_int = [int(_) for _ in configuration_str]
            num_particles = sum(configuration_int)
            if num_particles in configurations_stab:
                configurations_stab[num_particles].append(configuration_int)
            else:
                configurations_stab[num_particles]= [configuration_int]
        self.configurations_stab = configurations_stab
        
        
    
    def generate_configuration_loss_qnderror(self, num_l, num_q, num_stab):
        configuration_loss = []
        configuration_qnd = []
        configuration_loss = self.configurations[num_l]
        configuration_qnd = self.configurations[num_q]   
        configurations_stab = self.configurations_stab[num_stab]   
    
        return product(configuration_loss, configuration_qnd, configurations_stab)        
        

def check_correctable_state_analytics(random_losses, qnd_errors):

    non_correctable_3_events = [[0, 1, 4], [0, 2, 5], [0, 3, 6], [1, 2, 6], [2, 3, 4], [4, 5, 6], [1, 3, 5]]
    correctable_4_events = [[0, 1, 2, 3], [1, 2, 4, 5], [2, 3, 5, 6], [0, 3, 4, 5], [0, 1, 5, 6], [1, 3, 4, 6], [0, 2, 4, 6]]

    num_loss = sum(random_losses)
    num_qnd = sum(qnd_errors)  

    guessed_loss = [(loss + qnd_err) % 2 for loss, qnd_err in zip(random_losses, qnd_errors)]
    num_fresh_qubits = sum(guessed_loss)

    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()
    
    #find if a qnd_error hits a loss    
    qnderror_hit_loss = any([(_ in position_loss) for _ in position_qnd])
    
    if qnderror_hit_loss:
        non_correctable = True
        correctable = False
    else:
        if  num_fresh_qubits in (0,1,2):
            non_correctable = False
            correctable = True
        elif num_fresh_qubits == 3:
            if sorted(position_loss + position_qnd) in non_correctable_3_events:
                non_correctable = True
                correctable = False
            else:
                non_correctable = False
                correctable = True

        elif num_fresh_qubits == 4:
            if sorted(position_loss + position_qnd) in correctable_4_events:
                non_correctable = False
                correctable = True
            else:
                non_correctable = True
                correctable = False

        elif num_fresh_qubits in (5,6,7):
            non_correctable = True
            correctable = False

    return [correctable, non_correctable]
    
    
    
    
def check_correctable_state(random_losses, qnd_errors):
    """
    Check whether a loss + qnd error event is correctable.
    The cases that are correctable with no doubts are the ones where
    a. no loss happens and 0,1,2 qnd error happen
    or
    b. no qnd errors happens and 0,1,2 losses happen
    or
    c. no qnd errors happens on the losses and the number of replaced qubits is 0,1,2
    
     
    The cases that are NOT correctable with no doubts are the ones where
    d. one qnd error happens on the position of a loss
    or
    e. five, six or seven fresh qubits are introduced for replacing the actual losses or the guessed ones

    In all the other cases, one should check
    """
    
    num_loss = sum(random_losses)
    num_qnd = sum(qnd_errors)  

    guessed_loss = [(loss + qnd_err) % 2 for loss, qnd_err in zip(random_losses, qnd_errors)]
    num_fresh_qubits = sum(guessed_loss)

    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()
    
    #find if a qnd_error hits a loss    
    qnderror_hit_loss = any([(_ in position_loss) for _ in position_qnd])
    
    if qnderror_hit_loss:
        non_correctable = True
        correctable = False
        to_check = False
    else:
        if   num_fresh_qubits in (0,1,2):
            non_correctable = False
            correctable = True
            to_check = False
        elif num_fresh_qubits in (3,4):
            non_correctable = False
            correctable = False
            to_check = True

        elif num_fresh_qubits in (5,6,7):
            non_correctable = True
            correctable = False
            to_check = False

    return [correctable, non_correctable, to_check]