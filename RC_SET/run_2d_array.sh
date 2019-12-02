#!/bin/bash

vmin=0
vr=0
vg_avg=0
vg_std=0
c_avg=1
c_std=0
r_avg=10
r_std=0
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
rg_avg=1000
cg_avg=10
cg_std=0
f="big_array"
full="" 
graph=""
resume=""
currentmap="--current-map"
vmax=2
vSym="--symmetric-v"
T=0.01
custom_rh="\"\""
custom_rv="\"\""
dist="exp"
leaping="--tau-leaping"

if [ ! -d "$f" ]; then
  mkdir "$f"
fi

rows=10
columns=10

c_std=0.5
r_std=9
name="array_10_10_tau_leaping"
qsub -V -S /bin/bash -cwd -N "$name" -q jdubi.q -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" run_RC_SET.sh