#!/bin/bash
export OMP_NUM_THREADS=1

eps=0
chi=0

for overrot1 in 0.01 # 0.03162277660  #0.01 # 0.01 #0.03162277660
do
for overrot2 in 0.0961665 #0.136 # 0.4300697618 #0 #0.136 0.4300697618
do

folder=$(printf "case_1_new_2/chi_%.01e_eps_%1.3f_p2_%1.3f_p1_%1.3f" $chi $eps ${overrot2} ${overrot1})
echo ${overrot1}, ${overrot2}, ${folder}
mkdir -p  $folder
i=0
for state in  2 #2
do
for c in  $(seq -f "%1.4f" 0.050 0.050 0.50)  
#$(seq -f "%1.4f" 0.050 0.050 0.950) #£$(seq -f "%1.4f" 0.005 0.005 0.045) #  $(seq -f "%1.4f" 0.050 0.050 0.950)
do
    echo $c $state
     python simulation_all_over_stabilizers_overrotation.py \
        --logical_state ${state}  \
        --phi_tilde ${c}  \
        --epsilon_choi ${eps} \
        --chi_threshold ${chi} \
        --p_overrot_1 ${overrot1} \
        --p_overrot_2 ${overrot2} \
        --dir_name $folder \
        2>&1 >& 0_log_${c}_${eps}_${state}_${overrot2}_${overrot1}.log &
    
         pids[${i}]=$!
         i=$((i+1))
done
done

for pid in ${pids[*]}; do
    wait $pid
done

done
done




for overrot1 in 0.01 # 0.03162277660  #0.01 # 0.01 #0.03162277660
do
for overrot2 in 0.0961665 #0.136 # 0.4300697618 #0 #0.136 0.4300697618
do

folder=$(printf "case_1_new_2/chi_%.01e_eps_%1.3f_p2_%1.3f_p1_%1.3f" $chi $eps ${overrot2} ${overrot1})
echo ${overrot1}, ${overrot2}, ${folder}
mkdir -p  $folder
i=0
for state in  2 #2
do
for c in $(seq -f "%1.4f" 0.550 0.050 0.950) 
#$(seq -f "%1.4f" 0.050 0.050 0.950) #£$(seq -f "%1.4f" 0.005 0.005 0.045) #  $(seq -f "%1.4f" 0.050 0.050 0.950)
do
    echo $c $state
     python simulation_all_over_stabilizers_overrotation.py \
        --logical_state ${state}  \
        --phi_tilde ${c}  \
        --epsilon_choi ${eps} \
        --chi_threshold ${chi} \
        --p_overrot_1 ${overrot1} \
        --p_overrot_2 ${overrot2} \
        --dir_name $folder \
        2>&1 >& 0_log_${c}_${eps}_${state}_${overrot2}_${overrot1}.log &

         pids[${i}]=$!
         i=$((i+1))
done
done

for pid in ${pids[*]}; do
    wait $pid
done

done
done


