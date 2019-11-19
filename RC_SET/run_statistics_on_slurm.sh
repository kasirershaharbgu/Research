#!/bin/bash

vmin=0
vr=0
vg_avg=0
vg_std=0
c_avg=1
c_std=0
cg_std=0
r_avg=1
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
rg_avg=1000
cg_avg=10
t="03:00:00"
dist="exp"
r_std=0
graph=""
currentmap=""
resume=""
full=""
vSym="--symmetric-v"
resume=""
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
T=0
vmax=2
f="3X3_array_statistics"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
rows=3
columns=3


for c_std in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
	for r_std in 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5
	do
		for run in 1 2 3 4 5 6 7 8 9 10
		do
		name="c_std_${c_std}_r_std_${r_std}_run_${run}"
		sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=5000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym" run_RC_SET_on_slurm.sh
		done
	done
done
