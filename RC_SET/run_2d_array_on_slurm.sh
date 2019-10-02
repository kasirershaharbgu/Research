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
rg_avg=1000
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
r_std=9
f="2D_Array/slurm_run_different_disorder"

name="uniform_r_disorder_1_10"
rows=1
columns=10
v_max=20
t="10:00:00"
dist="uniform"

sbatch -J "name" -t "$t" --output "$f/$name.out" --ntasks-per-node "$repeats" --export=M="$rows",N="$columns",vmin="$vmin",vmax="$v_max",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",cstd="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist" run_RC_SET_on_slurm.sh

sleep 3

name="uniform_r_disorder_10_1"
rows=10
columns=1
v_max=2
t="4:00:00"
dist="uniform"

sbatch -J "name" -t "$t" --output "$f/$name.out" --ntasks-per-node "$repeats" --export=M="$rows",N="$columns",vmin="$vmin",vmax="$v_max",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",cstd="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist" run_RC_SET_on_slurm.sh

sleep 3

name="uniform_r_disorder_10_10"
rows=10
columns=10
v_max=20
t="24:00:00"
dist="uniform"

sbatch -J "name" -t "$t" --output "$f/$name.out" --ntasks-per-node "$repeats" --export=M="$rows",N="$columns",vmin="$vmin",vmax="$v_max",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",cstd="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist" run_RC_SET_on_slurm.sh

sleep 3

name="two_points_r_disorder_1_10"
rows=1
columns=10
v_max=20
t="10:00:00"
dist="two_points"

sbatch -J "name" -t "$t" --output "$f/$name.out" --ntasks-per-node "$repeats" --export=M="$rows",N="$columns",vmin="$vmin",vmax="$v_max",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",cstd="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist" run_RC_SET_on_slurm.sh

sleep 3

name="two_points_r_disorder_10_1"
rows=10
columns=1
v_max=2
t="4:00:00"
dist="two_points"

sbatch -J "name" -t "$t" --output "$f/$name.out" --ntasks-per-node "$repeats" --export=M="$rows",N="$columns",vmin="$vmin",vmax="$v_max",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",cstd="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist" run_RC_SET_on_slurm.sh

sleep 3

name="two_points_r_disorder_10_10"
rows=10
columns=10
v_max=20
t="24:00:00"
dist="two_points"

sbatch -J "name" -t "$t" --output "$f/$name.out" --ntasks-per-node "$repeats" --export=M="$rows",N="$columns",vmin="$vmin",vmax="$v_max",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",cstd="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist" run_RC_SET_on_slurm.sh

sleep 3


