#!/bin/bash

N=100
vmin=0
tstep=1000
vg_avg=0
c_avg=1
cg_avg=1
r_avg=1
rg_avg=1000
repeats=10
n_avg=0
n_std=0
q_avg=0
q_std=0
dt=0.01
vmax=100
vstep=1
f="1D_Array/one_dimension_big_array"

qsub -V -S /bin/bash -cwd -N r_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "no_disorder_N_$N" "$f"

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 1  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "r_disorder_N_$N" "$f"   

sleep 3

qsub -V -S /bin/bash -cwd -N c_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 0.1 "$cg_avg" 0 "$r_avg" 0 "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "c_disorder_N_$N" "$f" 

sleep 3


qsub -V -S /bin/bash -cwd -N rg_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 100 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "rg_disorder_N_$N" "$f"  

sleep 3

qsub -V -S /bin/bash -cwd -N cg_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 10 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "cg_disorder_N_$N" "$f"   

sleep 3

qsub -V -S /bin/bash -cwd -N vg_disprder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 10 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "vg_disorder_N_$N" "$f" 

sleep 3

qsub -V -S /bin/bash -cwd -N all_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 10 "$c_avg" 0 "$cg_avg" 1 "$r_avg" 1  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "all_disorder_N_$N" "$f"  