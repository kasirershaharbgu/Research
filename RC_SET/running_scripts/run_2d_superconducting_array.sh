#!/bin/bash

flags="-pe shared 14 -V -S /bin/bash -q jdubi.q -cwd"
vmin=0
vr=0
vg_avg=0
vg_std=0
c_avg=1
c_std=0.9
r_avg=10
r_std=9
repeats=25
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.001
rg_avg=1000
rg_std=0
cg_avg=10
cg_std=0
dist="exp"
graph=""
currentmap="--current-map"
full="--full"
resume="--resume"
vSym="--symmetric-v"
leaping=""
efermi=""
dbg=""
double_time=""
double_loop=""
input="super_conducting_array_different_gap/runningParameters_sc_array_5_5_T_0.001_cg_10_run_2_gap_0.1.txt"
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
temperature_gradient=0
vmax=5
efermi=""
T=0.001
superconducting="--superconducting"
f="super_conducting_array_different_gap"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
rows=5
columns=5
cg_avg=10
for gap in 0.01 0.05
do
	name="sc_array_${rows}_${columns}_T_${T}_gap_${gap}"
	qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",superconducting="$superconducting",gap=$gap,leaping="$leaping",efermi="$efermi",double_time="$double_time",double_loop="$double_loop",dbg="$dbg" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" running_scripts/run_RC_SET.sh
	sleep 1
done
