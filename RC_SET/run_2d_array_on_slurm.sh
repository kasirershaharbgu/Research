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

sbatch -J "no_disorder_10_1" -t "1:00:00" --output "no_disorder_10_1.out" --export=M=10,N=1,vmin="$vmin",vmax=2,vstep="$vstep",vg_avg="$vg_avg",vg_std=0,c_avg="$c_avg",cstd=0,cg_avg="$cg_avg",cg_std=0,r_avg="$r_avg",r_std=0,rg_avg="$rg_avg",rg_std=0,repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="no_disorder_1_1",out="$f",time="1:00:00" run_RC_SET_on_slurm.sh

sleep 3

sbatch -J "r_disorder_1_1" -t "1:00:00" --output "r_disorder_1_1.out" --export=M=1,N=1,vmin="$vmin",vmax=2,vstep="$vstep",vg_avg="$vg_avg",vg_std=0,c_avg="$c_avg",cstd=0,cg_avg="$cg_avg",cg_std=0,r_avg="$r_avg",r_std=1,rg_avg="$rg_avg",rg_std=0,repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="r_disorder_1_1",out="$f",time="1:00:00" run_RC_SET_on_slurm.sh  

sleep 3

sbatch -J "no_disorder_1_10" -t "5:00:00" --output "no_disorder_1_10.out" --export=M=1,N=10,vmin="$vmin",vmax=2,vstep="$vstep",vg_avg="$vg_avg",vg_std=0,c_avg="$c_avg",cstd=0,cg_avg="$cg_avg",cg_std=0,r_avg="$r_avg",r_std=0,rg_avg="$rg_avg",rg_std=0,repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="no_disorder_1_10",out="$f",time="5:00:00" run_RC_SET_on_slurm.sh  

sleep 3


sbatch -J "r_disorder_1_10" -t "5:00:00" --output "r_disorder_1_10.out" --export=M=1,N=10,vmin="$vmin",vmax=2,vstep="$vstep",vg_avg="$vg_avg",vg_std=0,c_avg="$c_avg",cstd=0,cg_avg="$cg_avg",cg_std=0,r_avg="$r_avg",r_std=1,rg_avg="$rg_avg",rg_std=0,repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="r_disorder_1_10",out="$f",time="5:00:00" run_RC_SET_on_slurm.sh  

sleep 3

sbatch -J "no_disorder_10_10" -t "5:00:00" --output "no_disorder_10_10.out" --export=M=10,N=10,vmin="$vmin",vmax=2,vstep="$vstep",vg_avg="$vg_avg",vg_std=0,c_avg="$c_avg",cstd=0,cg_avg="$cg_avg",cg_std=0,r_avg="$r_avg",r_std=0,rg_avg="$rg_avg",rg_std=0,repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="no_disorder_10_10",out="$f",time="10:00:00" run_RC_SET_on_slurm.sh   

sleep 3

sbatch -J "r_disorder_10_10" -t "10:00:00" --output "r_disorder_10_10.out" --export=M=10,N=10,vmin="$vmin",vmax=2,vstep="$vstep",vg_avg="$vg_avg",vg_std=0,c_avg="$c_avg",cstd=0,cg_avg="$cg_avg",cg_std=0,r_avg="$r_avg",r_std=1,rg_avg="$rg_avg",rg_std=0,repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="r_disorder_10_10",out="$f",time="10:00:00" run_RC_SET_on_slurm.sh 