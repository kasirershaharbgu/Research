#!/bin/bash

vmin=0
vr=0
vg_avg=0
vg_std=0
c_avg=1
c_std=0
cg_std=0
r_avg=10
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
rg_avg=1000
cg_avg=10
t="10:00:00"
dist="exp"
r_std=0
graph=""
currentmap=""
resume="--resume"
vSym="--symmetric-v"
resume=""
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
T=0
vmax=2
f="3X3_array_statistics_r_avg_10"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
rows=3
columns=3
c_std=0.4
for r_std in 9
do
	for run in 1
	do
	name="c_std_${c_std}_r_std_${r_std}_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=5000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym" run_RC_SET_on_slurm.sh
	done
done

