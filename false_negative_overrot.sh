#!/bin/bash
export OMP_NUM_THREADS=1

eps=0.0
chi=0.0
poverrot=0.136
falseneg=0.05
folder=$(printf "chi_%.01e_eps_%1.3f_overrot_p_%1.3f_fn_%1.3f" $chi $eps $poverrot $falseneg)
mkdir -p  ${folder}
wait

for state in  2
do
for c in 0.0 # 0.1  0.15  0.2  0.25  0.3  0.35  0.4  0.45  0.5  0.55  0.6  0.65  0.7  0.75  0.8  0.85  0.9  0.95
do
     python simulation_false_negative.py \
        --logical_state ${state}  --phi_tilde ${c}  \
        --epsilon_choi ${eps} --chi_threshold ${chi} \
        --p_overrot ${poverrot} \
        --falseneg ${falseneg} \
         --dir_name ${folder}
#        2>&1 >& 0_log_${c}_${eps}_${state}.log &
done
done
