import numpy as np


def check_correctable_state(random_losses, qnd_errors):
    """
    Check whether a loss + qnd error event is correctable.
    The cases that are correctable with no doubts are the ones where
    a. no loss happens and 0,1,2 qnd error happen
    or
    b. no qnd errors happens and 0,1,2 losses happen

    The cases that are NOT correctable with no doubts are the ones where
    c. one qnd error happens on the position of a loss
    or
    d. five, six or seven fresh qubits are introduced for replacing the actual losses or the guessed ones

    In all the other cases, one should check
    """
    
    num_loss = sum(random_losses)
    num_qnd = sum(qnd_errors)  
       
    case_a = (num_loss == 0) and (num_qnd in (0,1,2))
    case_b = (num_loss in (0,1,2)) and (num_qnd == 0)
    
    correctable_event = case_a or case_b


    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    case_c = any([(_ in position_loss) for _ in position_qnd])

    guessed_loss = [(loss + qnd_err) % 2 for loss, qnd_err in zip(random_losses, qnd_errors)]
    num_fresh_qubits = sum(guessed_loss)
    if not case_c: print("num_fresh_qubits", num_fresh_qubits)
    case_d = num_fresh_qubits in (5,6,7)

    non_correctable_event = case_c or case_d

    return {"correctable": correctable_event, "non_correctable_event": non_correctable_event }