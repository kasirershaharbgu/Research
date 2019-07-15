#!/bin/bash

arrayN = ( 2 10 100 )
vmin = 0
vmax = 6
vstep = 0.1
tstep = 100
vg_avg = 0
c_avg = 1
cg_avg = 1
r_avg = 1
rg_avg = 1000
repeats = 10
n_avg = 0
n_std = 0
q_avg = 0
q_std = 0
dt = 0.01
folder = one_dimension_array_first_run

for N in "${arrayN[@]}"
do

qsub -V -S /bin/bash -cwd -N no_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" no_disorder_N_"$N" "$folder"  

sleep 10

qsub -V -S /bin/bash -cwd -N one_dimension_array_r_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 10  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "r_disorder_N_$N" "$folder"   

sleep 10

qsub -V -S /bin/bash -cwd -N one_dimension_array_c_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 10 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "c_disorder_N_$N" "$folder" 

sleep 10


qsub -V -S /bin/bash -cwd -N one_dimension_array_rg_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 100 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "rg_disorder_N_$N" "$folder"  

sleep 10

qsub -V -S /bin/bash -cwd -N one_dimension_array_cg_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 0 "$c_avg" 0 "$cg_avg" 10 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "cg_disorder_N_$N" "$folder"   

sleep 10

qsub -V -S /bin/bash -cwd -N one_dimension_array_vg_disprder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 10 "$c_avg" 0 "$cg_avg" 0 "$r_avg" 0  "$rg_avg" 0 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "vg_disorder_N_$N" "$folder" 

sleep 10

qsub -V -S /bin/bash -cwd -N one_dimension_array_all_disorder -q jdubi.q run_RC_SET.sh 1 "$N" "$vmin" "$vmax" "$vstep" "$tstep" "$vg_avg" 10 "$c_avg" 10 "$cg_avg" 10 "$r_avg" 10  "$rg_avg" 100 "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$dt" "all_disorder_N_$N" "$folder"  

sleep 10

done