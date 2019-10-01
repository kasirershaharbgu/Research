#!/bin/bash
vmin=0
vg_avg=0
vg_std=0
c_avg=1
c_std=0
cg_avg=10
cg_std=0
r_avg=10
rg_avg=1000
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
f="2D_arrays_different_dist_longer_time"

qsub -V -S /bin/bash -cwd -N no_disorder_1_1 -q jdubi.q run_RC_SET.sh 1 1 "$vmin" 2 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 0  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "no_disorder_big_cg_1_1" "for_weizmann" "uniform"

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder_two_points_1_1 -q jdubi.q run_RC_SET.sh 1 1 "$vmin" 2 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 9  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_big_cg_1_1" "for_weizmann" "two_points"

sleep 3

qsub -V -S /bin/bash -cwd -N no_disorder_small_cg_1_1 -q jdubi.q run_RC_SET.sh 1 1 "$vmin" 2 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" 1 "$cg_std" "$r_avg" 0  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "no_disorder_small_cg_1_1" "for_weizmann" "uniform"

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder_small_cg_1_1 -q jdubi.q run_RC_SET.sh 1 1 "$vmin" 2 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" 1 "$cg_std" "$r_avg" 9  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_small_cg_1_1" "for_weizmann" "two_points"

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder_two_points_10_1 -q jdubi.q run_RC_SET.sh 10 1 "$vmin" 2 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 9  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_two_points_10_1" "$f" "two_points"

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder_uniform_10_1 -q jdubi.q run_RC_SET.sh 10 1 "$vmin" 2 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 9  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_uniform_10_1" "$f" "uniform"

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder_uniform_1_10 -q jdubi.q run_RC_SET.sh 1 10 "$vmin" 20 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 9  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_uniform_1_10" "$f" "uniform"  

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder_two_points_1_10 -q jdubi.q run_RC_SET.sh 1 10 "$vmin" 20 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 9  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_two_points_1_10" "$f" "two_points"

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder_uniform_10_10 -q jdubi.q run_RC_SET.sh 10 10 "$vmin" 20 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 9  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_uniform_10_10" "$f" "uniform"

sleep 3

qsub -V -S /bin/bash -cwd -N r_disorder_two_points_10_10 -q jdubi.q run_RC_SET.sh 10 10 "$vmin" 20 0.01 "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" 9  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "r_disorder_two_points_10_10" "$f" "two_points"