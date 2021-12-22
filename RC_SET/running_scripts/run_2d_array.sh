#!/bin/bash

flags="-pe shared 14 -V -S /bin/bash -q jdubi.q -cwd"
vmin=0
vr=0
vg_avg=0
vg_std=0
c_avg=1
r_avg=10
repeats=50
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.001
rg_avg=500
rg_std=0
cg_std=0
cg_avg=10
dist="exp"
graph=""
currentmap=""
full=""
resume=""
vSym="--symmetric-v"
leaping=""
efermi=""
dbg=""
double_time=""
double_loop=""
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
temperature_gradient=0
vmax=2
efermi=""
input=""
gap=0
T=0.001
f="different_size_arrays_finite_temperature"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
custom_rh="\"\""
custom_rv="\"\""
ustom_ch="\"\""
custom_cv="\"\""
c_std=0.9
r_std=9
double_loop=""
double_time=""
for length in 6 7 8
do
	for run in 1 2 3 4 5
	do
		rows=$length
		columns=$length
		name="array_${rows}_${columns}_run_${run}"
		qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",gap="$gap",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",double_time="$double_time",double_loop="$double_loop",dbg="$dbg" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" running_scripts/run_RC_SET.sh
		sleep 1
	done
done	

