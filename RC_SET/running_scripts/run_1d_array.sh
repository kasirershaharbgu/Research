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
vstep=0.001
rg_avg=1000
cg_avg=10
dist="exp"
r_std=9
graph=""
currentmap=""
full="--full"
resume=""
vSym="--symmetric-v"
leaping=""
efermi=""
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
T=0
vmax=2
f="1d_array_small_step_bgu"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
rows=1

for run in 1 2 3
do

	for columns in 1 2 3 4 5
	do
	
		efermi=""
		cg_std=0
		name="array_1_${columns}_r_disorder_run_${run}"
		qsub -V -S /bin/bash -cwd -N "$name" -q jdubi.q -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" run_RC_SET.sh
	
		sleep 1
	done	
done	