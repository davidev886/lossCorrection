#!/bin/bash
export OMP_NUM_THREADS=1

eps=0.05
chi=0.00

mkdir $(printf "chi_%.01e_eps_%1.3f_erasure" $chi $eps )
wait

for state in 0 2 4
do
for c in 0.05  0.1  0.15  0.2  0.25  0.3  0.35  0.4  0.45  0.5  0.55  0.6  0.65  0.7  0.75  0.8  0.85  0.9  0.95
do
        python process_matrix_erasure_channel.py --logical_state ${state}  --phi_tilde ${c}  --epsilon_choi ${eps} --chi_threshold ${chi} 2>&1 >& 0_log_${c}_${eps}_${state}.log &
done
done
