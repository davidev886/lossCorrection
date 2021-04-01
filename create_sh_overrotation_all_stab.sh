#!/bin/bash

cat << EOM
#!/bin/bash
export OMP_NUM_THREADS=1

eps=0
chi=0

EOM


for overrot1 in 0 0.01 0.03162277660
do
for overrot2 in 0 0.136 0.4300697618
do
 
cat << EndOfMessage
echo ------------------------------------

overrot1=$overrot1
overrot2=$overrot2

EndOfMessage


cat << "EndOfMessage"

echo $overrot1, $overrot2
folder=$(printf "chi_%.01e_eps_%1.3f_p2_%1.3f_p1_%1.3f" $chi $eps ${overrot2} ${overrot1})
echo $folder
mkdir -p  $folder

for state in 0 2 4
do
for c in $(seq -f "%1.2f" 0.05 0.05 0.95)
do
     python simulation_all_over_stabilizers_overrotation.py \
        --logical_state ${state}  \
        --phi_tilde ${c}  \
        --epsilon_choi ${eps}
        --chi_threshold ${chi} \
        --p_overrot_1 ${overrot1} \
        --p_overrot_2 ${overrot2} \
        --dir_name $folder \
        2>&1 >& 0_log_${c}_${eps}_${state}_${overrot2}_${overrot1}.log &

done
done
wait

EndOfMessage

done
done