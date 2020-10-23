import numpy as np



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