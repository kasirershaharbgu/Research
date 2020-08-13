#!/bin/bash

M=50 
N=50
T=0.01
R=1 
alpha=0 
vmin=0 
vmax=3 
vstep=0.01 
vr=0 
vu=0 
vd=0 
dist=exp 
out=thermal_arrays_slurm 
navg=10
vsym=--symmetric-v
t="10:00:00"

for nstd in 9 9.5 9.9
	for run in 1 2 3
		name="thermal_array_$M_$N_nstd_$nstd_run_$run"
		sbatch -J="$name" --nodes=1 --ntasks-per-node=1 --mem=1GB --time="$t" --partition=dept --mail-user=kasirer@post.bgu.ac.il --mail-type=end --output="$out/$name.out"  --export=M="$M",N="$N",T="$T",R="$R",alpha="$alpha",vmin="$vmin",vmax="$vmax",vstep="$vstep",vr="$vr",vu="$vu",vd="$vd",dist="$dist",out="$out",name="$name",navg="$navg",nstd="$nstd",vsym="$vsym" run_thermal_array_slurm.sh
	done
done
