import numpy as np

from itertools import product


class binary_raw_configurations(object):
    def __init__(self, n=6):
        self.n = n
        configurations = []
        for i in range(2**n):
            configuration_str = bin(i)[2:].zfill(n)
            configuration_int = [int(_) for _ in configuration_str]
            configurations.append(configuration_int)
        self.configurations = configurations


class binary_configurations(object):
    def __init__(self, n=7):
        self.n = n
        configurations = {}
        for i in range(2**n):
            configuration_str = bin(i)[2:].zfill(n)
            configuration_int = [int(_) for _ in configuration_str]
            num_particles = sum(configuration_int)
            if num_particles in configurations:
                configurations[num_particles].append(configuration_int)
            else:
                configurations[num_particles] = [configuration_int]
        self.configurations = configurations

        conf_temp = []
        for j in range(n + 1):
            conf_temp.extend(configurations[j])
        self.configurations_list = conf_temp

    def generate_configuration_loss_qnderror(self, num_l, num_q):
        configuration_loss = []
        configuration_qnd = []
        configuration_loss = self.configurations[num_l]
        configuration_qnd = self.configurations[num_q]

        return product(configuration_loss, configuration_qnd)


class binary_configurations_loss_qnd_stab(object):
    def __init__(self, n=7, stab=3):
        self.n = n
        configurations = {}
        for i in range(2**n):
            configuration_str = bin(i)[2:].zfill(n)
            configuration_int = [int(_) for _ in configuration_str]
            num_particles = sum(configuration_int)
            if num_particles in configurations:
                configurations[num_particles].append(configuration_int)
            else:
                configurations[num_particles] = [configuration_int]
        self.configurations = configurations

        configurations_stab = {}
        for i in range(2**stab):
            configuration_str = bin(i)[2:].zfill(stab)
            configuration_int = [int(_) for _ in configuration_str]
            num_particles = sum(configuration_int)
            if num_particles in configurations_stab:
                configurations_stab[num_particles].append(configuration_int)
            else:
                configurations_stab[num_particles] = [configuration_int]
        self.configurations_stab = configurations_stab

    def generate_configuration_loss_qnderror(self, num_l, num_q, num_stab):
        configuration_loss = []
        configuration_qnd = []
        configuration_loss = self.configurations[num_l]
        configuration_qnd = self.configurations[num_q]
        configurations_stab = self.configurations_stab[num_stab]

        return product(configuration_loss,
                       configuration_qnd,
                       configurations_stab)


class binary_configurations_false_pn(object):
    def __init__(self, n=7, stab=3):
        self.n = n
        configurations = {}
        for i in range(2**n):
            configuration_str = bin(i)[2:].zfill(n)
            configuration_int = [int(_) for _ in configuration_str]
            num_particles = sum(configuration_int)
            if num_particles in configurations:
                configurations[num_particles].append(configuration_int)
            else:
                configurations[num_particles] = [configuration_int]
        self.configurations = configurations

        configurations_stab = {}
        for i in range(2**stab):
            configuration_str = bin(i)[2:].zfill(stab)
            configuration_int = [int(_) for _ in configuration_str]
            num_particles = sum(configuration_int)
            if num_particles in configurations_stab:
                configurations_stab[num_particles].append(configuration_int)
            else:
                configurations_stab[num_particles] = [configuration_int]
        self.configurations_stab = configurations_stab

    def generate_configuration_loss_qnderror(self, num_l, num_q, num_stab):
        configuration_loss = []
        configuration_qnd = []
        configuration_loss = self.configurations[num_l]
        configuration_qnd = self.configurations[num_q]
        configurations_stab = self.configurations_stab[num_stab]

        return product(configuration_loss,
                       configuration_qnd,
                       configurations_stab)


def check_correctable_state_analytics(random_losses, qnd_errors):

    non_correctable_3_events = [[0, 1, 4], [0, 2, 5], [0, 3, 6],
                                [1, 2, 6], [2, 3, 4], [4, 5, 6], [1, 3, 5]]
    correctable_4_events = [[0, 1, 2, 3], [1, 2, 4, 5],
                            [2, 3, 5, 6], [0, 3, 4, 5],
                            [0, 1, 5, 6], [1, 3, 4, 6],
                            [0, 2, 4, 6]
                            ]

    guessed_loss = [(loss + qnd_err) % 2
                    for loss, qnd_err in zip(random_losses, qnd_errors)]
    num_fresh_qubits = sum(guessed_loss)

    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    # find if a qnd_error hits a loss
    qnderror_hit_loss = any([(_ in position_loss) for _ in position_qnd])

    if qnderror_hit_loss:
        non_correctable = True
        correctable = False
    else:
        if num_fresh_qubits in (0, 1, 2):
            non_correctable = False
            correctable = True
        elif num_fresh_qubits == 3:
            if sorted(position_loss + position_qnd) \
               in non_correctable_3_events:
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

        elif num_fresh_qubits in (5, 6, 7):
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
    c. no qnd errors happens on the losses and the number
       of replaced qubits is 0,1,2

    The cases that are NOT correctable with no doubts are the ones where
    d. one qnd error happens on the position of a loss
    or
    e. five, six or seven fresh qubits are introduced for replacing
       the actual losses or the guessed ones

    In all the other cases, one should check
    """

    guessed_loss = [(loss + qnd_err) % 2 for
                    loss, qnd_err in zip(random_losses, qnd_errors)]
    num_fresh_qubits = sum(guessed_loss)

    position_loss = np.where(random_losses)[0].tolist()
    position_qnd = np.where(qnd_errors)[0].tolist()

    # find if a qnd_error hits a loss
    qnderror_hit_loss = any([(_ in position_loss) for _ in position_qnd])

    if qnderror_hit_loss:
        non_correctable = True
        correctable = False
        to_check = False
    else:
        if num_fresh_qubits in (0, 1, 2):
            non_correctable = False
            correctable = True
            to_check = False
        elif num_fresh_qubits in (3, 4):
            non_correctable = False
            correctable = False
            to_check = True

        elif num_fresh_qubits in (5, 6, 7):
            non_correctable = True
            correctable = False
            to_check = False

    return [correctable, non_correctable, to_check]


def get_ancilla_outcomes_false_negatives(L):
    all_binary_ordered_confs = []
    for num_loss, loss_confs in binary_configurations().configurations.items():
        all_binary_ordered_confs.extend(loss_confs)

    all_configurations = []
    index_outcomes_ancilla = 0
    for outcomes_ancilla in all_binary_ordered_confs:

        false_negative_confs = []
        for false_neg_events_all in all_binary_ordered_confs:
            # print("outcomes_ancilla", outcomes_ancilla, false_neg_events_all)
            false_neg_events = [0] * L
            for _ in range(L):
                if outcomes_ancilla[_] == 0:
                    false_neg_events[_] = false_neg_events_all[_]
                else:
                    false_neg_events[_] = 0

            losses = np.where(outcomes_ancilla)[0].tolist()
            false_negative_qubits = np.where(false_neg_events)[0].tolist()
            kept_qubits = [_ for _ in range(L)
                           if (_ not in false_negative_qubits)
                           and (_ not in losses)]
            if false_neg_events not in false_negative_confs \
               and \
               len(kept_qubits):
                false_negative_confs.append(false_neg_events)

        all_configurations.append([outcomes_ancilla, false_negative_confs])
        index_outcomes_ancilla += 1
    return all_configurations


def get_ancilla_outcomes_false_positives(L):
    all_binary_ordered_confs = []
    for num_loss, loss_confs in binary_configurations().configurations.items():
        all_binary_ordered_confs.extend(loss_confs)

    all_configurations = []
    index_outcomes_ancilla = 0
    for outcomes_ancilla in all_binary_ordered_confs:

        false_positive_confs = []
        for false_pos_events_all in all_binary_ordered_confs:

            false_pos_events = [0] * L
            for _ in range(L):
                if outcomes_ancilla[_] == 1:
                    false_pos_events[_] = false_pos_events_all[_]
                else:
                    false_pos_events[_] = 0

            losses = np.where(outcomes_ancilla)[0].tolist()
            false_positive_qubits = np.where(false_pos_events)[0].tolist()
            kept_qubits = [_ for _ in range(L)
                           if (_ not in false_positive_qubits)
                           and (_ not in losses)]

            if false_pos_events not in false_positive_confs \
                and sum(false_pos_events) < 7:
               false_positive_confs.append(false_pos_events)

        all_configurations.append([outcomes_ancilla, false_positive_confs])
        index_outcomes_ancilla += 1
    return all_configurations

def create_random_event(prob_loss, basic_event_probs):
    id_event = []
    event = []
    for j_qubit in range(7):
        ancilla = np.random.binomial(n=1, p=prob_loss, size=1)

        ancilla_0 = np.random.multinomial(n=1, size=1,
                                          pvals=[basic_event_probs['0'],
                                                 basic_event_probs['1'],
                                                 basic_event_probs['2'],
                                                 basic_event_probs['3']]
                                          )[0]

        ancilla_1 = np.random.multinomial(n=1, size=1,
                                          pvals=[basic_event_probs['4'],
                                                 basic_event_probs['5']]
                                          )[0]
        if ancilla:
            event.append([ancilla[0], np.where(ancilla_1)[0][0]])
            id_event.append(str(4 + np.where(ancilla_1)[0][0]))
        else:
            event.append([ancilla[0], np.where(ancilla_0)[0][0]])
            id_event.append(str(np.where(ancilla_0)[0][0]))

    event_str = "".join(id_event)

    return event, event_str



if __name__ == "__main__":
    cc = 0
    a = get_ancilla_outcomes_false_negatives(7)
    for ifd, [x, y] in enumerate(a):
        print(ifd, x, "||")
        for el in y:
                print(el)
        cc += len(y)
    print(cc)
