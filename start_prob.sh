#!/bin/bash
export OMP_NUM_THREADS=1
    
eps=0.0
chi=0.0
folder=$(printf "zzz_chi_%.01e_eps_%1.3f" $chi $eps )
mkdir -p  $folder
wait

for state in 3
do
for c in 0.1
do
        python simulation_all_over_stabilizers.py \
         --logical_state ${state} \
         --phi_tilde ${c} \
         --epsilon_choi ${eps} \
         --chi_threshold ${chi} \
         --dir_name $folder
           #2>&1 >& 0_log_${c}_${eps}_${state}.log & 
done
done



#0.05  0.1  0.15  0.2  0.25  0.3  0.35  0.4  0.45  0.5  0.55  0.6  0.65  0.7  0.75  0.8  0.85  0.9  0.95