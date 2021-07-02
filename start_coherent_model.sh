#!/bin/bash
 export OMP_NUM_THREADS=1

eps=0
chi=0

for overrot1 in $(seq -f "%1.5f" 0.000 0.00100 0.0141421)
do
for overrot2 in $(seq -f "%1.5f" 0.000 0.0100 0.16) 0.136
do
    folder=$(printf "case_1_zero_loss_final/chi_%.01e_eps_%1.3f_p2_%1.5f_p1_%1.5f" $chi $eps ${overrot2} ${overrot1})
    echo ${overrot1}, ${overrot2}, ${folder}
    mkdir -p  $folder
    i=0
    c=0.00
    for state in 4 2 0
    do
    echo $overrot1 $overrot2  $state
    python  simulation_all_over_stabilizers_overrotation.py  \
        --logical_state ${state}  \
        --phi_tilde ${c}  \
        --epsilon_choi ${eps} \
        --chi_threshold ${chi} \
        --p_overrot_1 ${overrot1} \
        --p_overrot_2 ${overrot2} \
        --dir_name $folder \
        --num_trials 0 \
        2>&1 >& 0_log_${c}_${eps}_${state}_${overrot2}_${overrot1}.log &

         pids[${i}]=$!
         i=$((i+1))
    done

done #for overrot2

    for pid in ${pids[*]}; do
       wait $pid
    done

done
