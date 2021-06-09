#!/bin/bash
export OMP_NUM_THREADS=1

eps=0
chi=0

#for overrot1 in 0.00 0.005773502691896258 0.00707107 0.01 0.0141421  # 0.03162277660  #0.01 # 0.01 #0.03162277660
#do
#for overrot2 in 0.00 0.034 0.017  0.0680 0.0785196 0.0961665 0.136 # 0.4300697618 #0 #0.136 0.4300697618
#do

for overrot1 in 0.10000 #$(seq -f "%1.5f" 0.000 0.00400 0.0141421)
do
for overrot2 in 0.20000 #$(seq -f "%1.5f" 0.000 0.0400 0.16) 0.136
do

    folder=$(printf "case_2_zero_loss_commented/chi_%.01e_eps_%1.3f_p2_%1.5f_p1_%1.5f" $chi $eps ${overrot2} ${overrot1})
    echo ${overrot1}, ${overrot2}, ${folder}
    mkdir -p  $folder
    i=0
    c=0.2
    for state in 2
    do
    echo $overrot1 $overrot2  $state
    python qec_incoherent_channel.py \
        --logical_state ${state}  \
        --phi_tilde ${c}  \
        --epsilon_choi ${eps} \
        --chi_threshold ${chi} \
        --p_overrot_1 ${overrot1} \
        --p_overrot_2 ${overrot2} \
        --dir_name $folder \
        --num_trials 0 \
        --verbose
#        2>&1 >& 0_log_${c}_${eps}_${state}_${overrot2}_${overrot1}.log &

         pids[${i}]=$!
         i=$((i+1))
    done

done #for overrot2


done
