#!/bin/bash
vmin=0
vg_avg=0.05
vg_std=0
c_avg=1
c_std=0
cg_avg=10
cg_std=0
r_avg=1
rg_avg=1000
rg_std=0
repeats=10
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.1
f="1D_Array/2D_array_trial"

qsub -V -S /bin/bash -cwd -N no_disorder_10_10 -q jdubi.q run_RC_SET.sh 10 10 "$vmin" 20 "$vstep" "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 0  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "no_disorder_10_10" "$f"   

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder_10_10 -q jdubi.q run_RC_SET.sh 10 10 "$vmin" 20 "$vstep" "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 1  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_10_10" "$f"