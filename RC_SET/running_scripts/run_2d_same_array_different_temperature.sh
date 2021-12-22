#!/bin/bash

flags="-pe shared 14 -V -S /bin/bash -q jdubi.q -cwd -q jdubi.q@sge183"
vmin=0
vr=0
vg_avg=0
vg_std=0
c_avg=2
c_std=0.1
r_avg=10
r_std=9
repeats=50
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.001
rg_avg=500
rg_std=0
cg_avg=10
cg_std=0
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
vmax=1.2
efermi=""
gap=0
f="same_array_different_temperature"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
custom_rh="\"\""
custom_rv="\"\""
ustom_ch="\"\""
custom_cv="\"\""
c_avg=2
c_std=0
r_avg=10
r_std=0
rows=10
columns=10
double_loop=""
double_time=""
input="same_array_different_temperature/runningParameters_array_10_10_disorder_c_std_0.1_r_std_9_run_1.txt"
for T in 0.0001 0.0002 0.0003 0.0004 0.0006 0.0007 0.0008 0.0009
do 
	name="array_10_10_T_${T}"
	qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",gap="$gap",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",double_time="$double_time",double_loop="$double_loop",dbg="$dbg" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" running_scripts/run_RC_SET.sh
	sleep 1
done	
