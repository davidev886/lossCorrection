#!/bin/bash
export OMP_NUM_THREADS=1

dt=$(date '+%Y%m%d%H%M') 

eps=0.0
chi=0.0
overrot_2=0.136
overrot_1=0.01
folder="$dt"_$(printf "mc_chi_%.01e_eps_%1.3f_p2_%1.3f_p1_%1.3f" $chi $eps ${overrot_2} ${overrot_1})
echo $folder

mkdir -p  $folder
wait

#done states 3
for state in 0 #2
do
for c in $(seq -f "%1.2f" 0.05 0.05 0.95) 
do
        python monte_carlo_simulation_false_positive_negative \ 
         --logical_state ${state} \
         --phi_tilde ${c} \
         --epsilon_choi ${eps} \
         --chi_threshold ${chi} \
         --dir_name $folder \
         --num_trials 10000 \
         --p_overrot_2 ${overrot_2} \
         --p_overrot_1 ${overrot_1} \
           2>&1 >& 0_log_${c}_${eps}_${state}_${overrot_2}_${overrot_1}.log & 
done
done

