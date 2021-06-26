import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D

marker_list = list(Line2D.markers.keys())


LogicalStates_str = ["0", "1", "+", "-", "+i", "-i"]

files = []

directory = './'
list_dirs = []
for filename in os.listdir(directory):
    if os.path.isdir(filename) and filename.startswith("chi"):
        list_dirs.append(filename)

figura = plt.figure()
ax = figura.add_subplot(111)
#ax.yaxis.tick_right()
ax.yaxis.set_ticks_position('both')
x = np.arange(0, 0.5, 0.01)
plt.plot(x, 1 - x, '-')
x_data = np.arange(0, 0.5, 0.01)
y_data = 1 - 7*x_data**3 + 21*x_data**5 - 21*x_data**6 + 6*x_data**7
# plt.plot(x, y_data, '-')

print(list_dirs)

final_prob_list = []
for folder in sorted(list_dirs):

    split_folder_name = folder.split("_")

    err_stab = 0.0 #float(split_folder_name[-1])
    p1 = float(split_folder_name[-1])
    p2 = float(split_folder_name[-3])
    epsilon = float(split_folder_name[-5])
    number_of_files = len([filename for filename in os.listdir(folder) if filename.startswith("2")])
    # if number_of_files == 1: continue
    prob_success_state = {}
    for filename in sorted(os.listdir(folder)):

        if (filename.startswith("2") and
            filename.endswith(".dat") and
            (("state_0_" in filename) or ("state_+_" in filename))
            ):
            epsilon = float(filename[:-4].split("_")[-1])
            log_state = filename[:-4].split("_")[-5]

            A = np.loadtxt(os.path.join(folder, filename))
            if len(A.shape) == 1:
                A = np.array([A])

            phi_tilde = A[0, 0]
            phi_tilde_str = f"{phi_tilde:1.6f}"
            prob_channel = A[:, 3]
            prob_ancilla = A[:, 4]
            prob_success = A[:, 5]

            key = log_state + "_" + phi_tilde_str
            final_prob = sum(prob_channel
                             * prob_ancilla
                             )
            log_rate_success = sum(prob_channel
                                   * prob_ancilla
                                   * prob_success
                                   )
            print(filename,
                  f"{final_prob:.8f}",
                  f"{1-final_prob:.4e}",
                  f"{log_rate_success:.8f}"
                  f"{1-log_rate_success:.4e}"
                  )

            if log_state in prob_success_state:
                prob_success_state[key].append(final_prob)
            else:
                prob_success_state[key] = [final_prob]

    phis = []
    final_prob = {}

exit()
    for key, val in prob_success_state.items():
        log_operator = key.split("_")[0]
        phi_tilde = float(key.split("_")[1])
        log_state = key.split("_")[0]
        final_prob_list.append(
            {'log_operator': log_operator,
             'phi_tilde': phi_tilde,
             'epsilon': epsilon,
             'p2': p2,
             'p1': p1,
             'err_stab': err_stab,
             'prob_success_state': prob_success_state[key][0]
             })


results_df = pd.DataFrame(final_prob_list)
print(results_df)

exit()
results_df.sort_values(by=['log_operator', 'phi_tilde'], inplace=True)

B = results_df.groupby(['p2',
                        'p1',
                        'err_stab',
                        'phi_tilde',
                        'epsilon'
                        ], as_index=False).mean()

C = B.groupby(['p2', 'p1', 'epsilon', 'err_stab'])

marker_idx = 2
out_files = []
for k, gr in C:
    p2, p1, epsilon, err_stab = k
    print(p2, p1, epsilon)
    sorted_data = gr.sort_values(by=['phi_tilde'])
    print(sorted_data)
    x = np.sin(np.pi * sorted_data['phi_tilde']/2)**2 / 2
    y = sorted_data['prob_success_state']
    label_plot = ('p_2 =' + f'{p2:1.5f}' +
                  ' p_1 =' + f'{p1:1.5f}'
#                  ' err_stab =' + f'{err_stab:1.2f}'
                  )
    plt.plot(x, y,
             linestyle='-',
             marker=marker_list[marker_idx],
             label=label_plot
             )
    marker_idx += 1
    name_csv = (f"eps_{epsilon:1.5f}_p2_{p2:1.5f}_" +
                f"p1_{p1:1.5f}_errorstab_{err_stab:1.5f}")
    out_files.append(name_csv + '.csv')
    print(name_csv)
    sorted_data.to_csv(name_csv + '.csv', index=False)
    print(sorted_data.values)
#    print(type(sorted_data.values))

plt.title("Incoherent over-rotations")
plt.xlabel("$p_{loss}$")
plt.ylabel("$p_{success}$")
plt.legend(fontsize=10)
plt.savefig("incoherent_errors_inverted.pdf",
            bbox_inches="tight"
            )
# plt.savefig("../notes/false_pos_neg_success_probability_choi.pdf", bbox_inches="tight")

for namef in out_files:
    print(f'"{namef:s}",')

dfs = [f"{namef:s}" for namef in out_files]
frames = [pd.read_csv(x) for x in dfs]

result = pd.concat(frames)

result.to_csv("raw_data_merged_case_2_zero_loss_qutrit.csv", index=False)

for namef in out_files:
    os.remove(namef)


# folder=("/Users/vodola/Dropbox/ColorCodeLattice/"
#         "losses_steane_code/choi_operators/notes/case_1/")
#
#
# plt.savefig(folder + "coherent_errors.pdf",
#             bbox_inches="tight"
#             )


