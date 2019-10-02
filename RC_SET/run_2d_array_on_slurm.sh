#!/bin/bash

vmin=0
vr=0
vg_avg=0
vg_std=1
c_avg=1
c_std=1
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
f="2D_Array/slurm_run"

sbatch -J "no_disorder_10_1" -t "1:00:00" --output "no_disorder_10_1.out" --export=M=10,N=1,vmin="$vmin",vmax=2,vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",cstd="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std=0,rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="no_disorder_1_1",out="$f",vr="$vr$ run_RC_SET_on_slurm.sh

sleep 3