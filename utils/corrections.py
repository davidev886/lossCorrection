import numpy as np
from utils.p_operators import *

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
    
def check_stabilizers_measurement(random_losses, qnd_errors, stab_errors_binary):
    [correctable_0, non_correctable_0] = check_correctable_state_analytics(random_losses, qnd_errors)
    if non_correctable_0:
        return [correctable_0, non_correctable_0]
    else:
        #Stabilizer that are affected by a measurement error
        faulty_stab_qubits = [el for j,el in enumerate(stab_qubits) if stab_errors_binary[j]]

        position_loss = np.where(random_losses)[0].tolist()
        position_qnd = np.where(qnd_errors)[0].tolist()

        guessed_loss = [(loss + qnd_err) % 2 for loss, qnd_err in zip(random_losses, qnd_errors)]
        position_guessed_loss = np.where(guessed_loss)[0].tolist() 


        correctable = not (any([any([(loss in stab) for loss in position_guessed_loss]) for stab in faulty_stab_qubits]))
        print(f"{'correctable':40}",  correctable)
        return [correctable, not correctable]
        
                
    
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