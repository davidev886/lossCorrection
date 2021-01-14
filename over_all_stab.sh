#!/bin/bash
export OMP_NUM_THREADS=1

eps=1.0
chi=1.0e-3
mkdir -p  $(printf "chi_%.01e_eps_%1.3f_over_stab" $chi $eps )
wait

for state in  1 3 5
do
for c in 0.05  0.1  0.15  0.2  0.25  0.3  0.35  0.4  0.45  0.5  0.55  0.6  0.65  0.7  0.75  0.8  0.85  0.9  0.95
do
        python simulation_all_over_stabilizers.py --logical_state ${state}  --phi_tilde ${c}  --epsilon_choi ${eps} --chi_threshold ${chi} 2>&1 >& 0_log_${c}_${eps}_${state}.log &
done
done


