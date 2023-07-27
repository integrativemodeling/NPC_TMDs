export cluster=0
export state=0

export run_dir=/home/ignacia/Research/yeast_NPC/modeling_2023/mod_TMDs_symmetric/

cp $run_dir/analys/A_models_clust${cluster}_${state}.txt $run_dir/analys/scoresA.txt
cp $run_dir/analys/B_models_clust${cluster}_${state}.txt $run_dir/analys/scoresB.txt
cp $run_dir/density.txt density.txt
cp $run_dir/symm_groups.txt .

nohup python /home/ignacia/SOFTW/imp-sampcon-2023/pyext/src/exhaust.py \
       -n Memb -p ../ -ra A_models_clust${cluster}_${state}.rmf3 -rb B_models_clust${cluster}_${state}.rmf3 -d density.txt \
       -m cpu_omp -c 8 -g 3.0 --align --ambiguity symm_groups.txt > clustering.log &

