#!/bin/bash
export OMP_NUM_THREADS=1

eps=0
chi=0

echo 
for overrot1 in 0.00 0.01 0.02       # $(seq -f "%1.4f" 0.02 0.02 0.100) #0.00 0.01 # 0.03162277660  #0.01 # 0.01 #0.03162277660
do
for overrot2 in 0.00 0.0961665 0.136 # 0.0961665 #0.136 # 0.4300697618 #0 #0.136 0.4300697618
do

folder=$(printf "case_2_new_3/chi_%.01e_eps_%1.5f_p2_%1.5f_p1_%1.5f" $chi $eps ${overrot2} ${overrot1})
mkdir -p  $folder
i=0
    for state in 0 2
    do
    for c in $(seq -f "%1.4f" 0.050 0.050 0.450) #$(seq -f "%1.4f" 0.500 0.050 0.950) #$(seq -f "%1.4f" 0.0050 0.0050 0.0450)
    do
  #   echo ${overrot1}, ${overrot2}, $c, $state, ${folder}
     python incoherent_channel_0_1.py \
        --logical_state ${state}  \
        --phi_tilde ${c}  \
        --epsilon_choi ${eps} \
        --chi_threshold ${chi} \
        --num_trials 4000 \
        --p_overrot_1 ${overrot1} \
        --p_overrot_2 ${overrot2} \
        --dir_name $folder \
        2>&1 >& 0_log_${c}_${eps}_${state}_${overrot2}_${overrot1}.log &

         pids[${i}]=$!

     echo ${pids[$i]}, ${overrot1}, ${overrot2}, $c, $state, ${folder}
         i=$((i+1))
    done # for c
    done # for state

for pid in ${pids[*]}; do
    wait $pid
done


done # for overrot2
done # for overrot1


echo
#
#for overrot1 in 0.01 # 0.03162277660  #0.01 # 0.01 #0.03162277660
#do
#for overrot2 in 0.0961665 #0.136 # 0.4300697618 #0 #0.136 0.4300697618
#do
#
#folder=$(printf "case_2_new_2/chi_%.01e_eps_%1.3f_p2_%1.3f_p1_%1.3f" $chi $eps ${overrot2} ${overrot1})
#echo ${overrot1}, ${overrot2}, ${folder}
#mkdir -p  $folder
#i=0
#for state in 0 2
#do
#for c in $(seq -f "%1.4f" 0.500 0.050 0.950)
#do
#     echo $c $state
#     python incoherent_channel_better_sampling.py \
#        --logical_state ${state}  \
#        --phi_tilde ${c}  \
#        --epsilon_choi ${eps} \
#        --chi_threshold ${chi} \
#        --num_trials 4000 \
#        --p_overrot_1 ${overrot1} \
#        --p_overrot_2 ${overrot2} \
#        --dir_name $folder \
#        2>&1 >& 0_log_${c}_${eps}_${state}_${overrot2}_${overrot1}.log &
#
#         pids[${i}]=$!
#         i=$((i+1))
#done
#done
#
#for pid in ${pids[*]}; do
#    wait $pid
#done
#
#
#done
#done
#
#
