#!/bin/bash

vmin=0
vg_avg=0.05
c_avg=1
cg_avg=10
r_avg=1
rg_avg=1000
repeats=10
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.1
f="2D_Array/first_run"

#sbatch run_RC_SET.slurm 1 1 "$vmin" 2 "$vstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "no_disorder_1_1" "$f" "1:00:00"

sleep 3

#sbatch run_RC_SET.slurm 1 1 "$vmin" 2 "$vstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 1  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_1_1" "$f" "1:00:00"   

sleep 3

#sbatch run_RC_SET.slurm 1 10 "$vmin" 20 "$vstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "no_disorder_1_10" "$f" "5:00:00"   

sleep 3


#sbatch run_RC_SET.slurm 1 10 "$vmin" 20 "$vstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 1  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_1_10" "$f" "5:00:00"  

sleep 3

#sbatch run_RC_SET.slurm 10 10 "$vmin" 20 "$vstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "no_disorder_10_10" "$f" "10:00:00"    

sleep 3

#sbatch run_RC_SET.slurm 10 10 "$vmin" 20 "$vstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 1  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_10_10" "$f" "10:00:00"  