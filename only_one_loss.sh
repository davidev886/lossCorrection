#!/bin/bash
export OMP_NUM_THREADS=1

eps=1.0
chi=1e-4

folder=$(printf "chi_%.01e_eps_%1.3f_one_loss" $chi $eps)
mkdir -p  ${folder}
sleep 1

for state in  0 2 4
do
for c in 0.0 0.1  0.15  0.2  0.25  0.3  0.35  0.4  0.45  0.5  0.55  0.6  0.65  0.7  0.75  0.8  0.85  0.9  0.95
do
     python simulation_false_negative.py \
        --logical_state ${state}  --phi_tilde ${c}  \
        --epsilon_choi ${eps} --chi_threshold ${chi} \
         --dir_name ${folder} \
        2>&1 >& 0_log_${c}_${eps}_${state}_${falseneg}.log &
done
done