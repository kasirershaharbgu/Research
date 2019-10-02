#!/bin/bash
vmin=0
vr=0
vg_avg=0
vg_std=0
c_avg=1
c_std=0
cg_avg=10
cg_std=0
r_avg=10
r_std=9
rg_avg=1000
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
f="2D_arrays_different_dist_longer_time"


M=1
N=10
dist="uniform"
name="uniform_r_disorder_1_10"
vmax=20
qsub -V -S /bin/bash -cwd -N "$name" -q jdubi.q run_RC_SET.sh "$M" "$N" "$vmin" "$vmax" "$vstep" "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" "$r_std"  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$name" "$f" "$dist" "$vr"

sleep 3

M=10
N=1
dist="uniform"
name="uniform_r_disorder_10_1"
vmax=2
qsub -V -S /bin/bash -cwd -N "$name" -q jdubi.q run_RC_SET.sh "$M" "$N" "$vmin" "$vmax" "$vstep" "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" "$r_std"  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$name" "$f" "$dist" "$vr"

sleep 3

M=10
N=10
dist="uniform"
name="uniform_r_disorder_10_10"
vmax=20
qsub -V -S /bin/bash -cwd -N "$name" -q jdubi.q run_RC_SET.sh "$M" "$N" "$vmin" "$vmax" "$vstep" "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" "$r_std"  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$name" "$f" "$dist" "$vr"

sleep 3

M=1
N=10
dist="two_points"
name="two_points_r_disorder_1_10"
vmax=20
qsub -V -S /bin/bash -cwd -N "$name" -q jdubi.q run_RC_SET.sh "$M" "$N" "$vmin" "$vmax" "$vstep" "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" "$r_std"  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$name" "$f" "$dist" "$vr"

sleep 3

M=10
N=1
dist="two_points"
name="two_points_r_disorder_10_1"
vmax=2
qsub -V -S /bin/bash -cwd -N "$name" -q jdubi.q run_RC_SET.sh "$M" "$N" "$vmin" "$vmax" "$vstep" "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" "$r_std"  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$name" "$f" "$dist" "$vr"

sleep 3

M=10
N=10
dist="two_points"
name="two_points_r_disorder_10_10"
vmax=20
qsub -V -S /bin/bash -cwd -N "$name" -q jdubi.q run_RC_SET.sh "$M" "$N" "$vmin" "$vmax" "$vstep" "$vg_avg" "$vg_std" "$c_avg" "$c_std" "$cg_avg" "$cg_std" "$r_avg" "$r_std"  "$rg_avg" "$rg_std" "$repeats" "$n_avg" "$n_std" "$q_avg" "$q_std" "$name" "$f" "$dist" "$vr"