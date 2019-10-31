#!/bin/bash

vmin=0
vr=0
vg_avg=0
c_avg=1
r_avg=10
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=1
rg_avg=1000
cg_avg=10
f="big_array"
full="" 
graph=""
currentmap="--current-map"
vmax=100
custom_rh="\"\""
custom_rv="\"\""
dist="exp"

if [ ! -d "$f" ]; then
  mkdir "$f"
fi

rows=20
columns=20

c_std=0.5
cg_std=1
r_std=5
vg_std=0
name="array_20_20"
qsub -V -S /bin/bash -cwd -N "$name" -q jdubi.q -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" run_RC_SET.sh