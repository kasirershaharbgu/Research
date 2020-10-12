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
vstep=0.01
rg_avg=1000
rg_std=0
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
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
temperature_gradient=0
vmax=5
efermi=""
superconducting=""
gap=0.1
T=0
f="bgu_3_3_custom_arrays"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
c_avg=1
c_std=0
r_avg=10
r_std=0

double_loop=""
double_time=""
input=""

rows=3
columns=3
cg_avg=5



custom_ch="\"[[1,1,1,1],[1,1,1,1],[1,1,1,1]]\""
custom_cv="\"[[1,1,1],[1,1,1]\""
custom_rh="\"[[10,10,10,10],[10,10,10,10],[10,10,10,10]]\""
custom_rv="\"[[10,10,10],[10,10,10]\""
name="array_${rows}_${columns}_const_r_const_c_cg_${cg_avg}"
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",gap="$gap",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",superconducting="$superconducting",double_time="$double_time",double_loop="$double_loop",dbg="$dbg" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" running_scripts/run_RC_SET.sh
sleep 1

cg_avg=3.45
custom_ch="\"[[1,0.1,0.1,1],[1,1,1,1],[1,0.1,0.1,1]]\""
custom_cv="\"[[1,0.01,1],[1,0.01,1]\""
custom_rh="\"[[10,10,10,10],[10,10,10,10],[10,10,10,10]]\""
custom_rv="\"[[10,10,10],[10,10,10]\""
name="array_${rows}_${columns}_const_r_custom_c_cg_${cg_avg}"
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",gap="$gap",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",superconducting="$superconducting",double_time="$double_time",double_loop="$double_loop",dbg="$dbg" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" running_scripts/run_RC_SET.sh
sleep 1

cg_avg=5
custom_ch="\"[[1,1,1,1],[1,1,1,1],[1,1,1,1]]\""
custom_cv="\"[[1,1,1],[1,1,1]\""
custom_rh="\"[[1,10,10,1],[1,1,1,1],[1,10,10,1]]\""
custom_rv="\"[[10,100,10],[10,100,10]\""
name="array_${rows}_${columns}_const_r_custom_c_cg_${cg_avg}"
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",gap="$gap",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",superconducting="$superconducting",double_time="$double_time",double_loop="$double_loop",dbg="$dbg" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" running_scripts/run_RC_SET.sh
sleep 1






