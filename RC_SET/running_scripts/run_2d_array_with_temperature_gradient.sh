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
vstep=0.01
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
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
up_electrode="\"\""
down_electrode="\"\""
left_electrode="\"\""
right_electrode="\"\""
vd=0
vu=0
temperature_gradient=0.01
temperature_gradient_step=0.001
vmax=1.2
efermi=""
gap=0
calc_it=""
f="linear_response_temperature_gradient"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
c_avg=2
c_std=0
r_avg=10
r_std=0
rows=10
columns=10
double_loop=""
double_time=""
T=0.001
input="linear_response_temperature_gradient/runningParameters_array_10_10_T_0.001.txt"
name="array_10_10_T_${T}"
qsub "$flags" -N "${name}" -o "$f/${name}.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",temperature_gradient_step="$temperature_gradient_step",file_name="${name}",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1	
