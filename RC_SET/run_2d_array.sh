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
f="2D_square_arrays_with_current_map_bgu"
full="" 
graph=""
currentmap="--current-map"
rows=5
columns=5
vmax=5
custom_rh="\"\""
custom_rv="\"\""
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
dist="two_points"
name="no_disorder_5_5"
r_std=0

qsub -S /bin/bash -cwd -N "$name" -q jdubi.q -o "$f/$name.out" -V M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET.sh

dist="uniform"
name="uniform_disorder_5_5"
r_std=9
qsub -S /bin/bash -cwd -N "$name" -q jdubi.q -o "$f/$name.out" -V M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET.sh

dist="two_points"
name="two_points_disorder_5_5"
r_std=9
qsub -S /bin/bash -cwd -N "$name" -q jdubi.q -o "$f/$name.out" -V M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET.sh

dist="exp"
name="exp_disorder_5_5"
r_std=9
qsub -S /bin/bash -cwd -N "$name" -q jdubi.q -o "$f/$name.out" -V M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET.sh

N=2
M=2
custom_rh="\"[[10,10,10],[1,1,1]]\""
custom_rv="\"[[10,10]]\""
dist="uniform"
name="one_nig_one_small_2_2"
qsub -S /bin/bash -cwd -N "$name" -q jdubi.q -o "$f/$name.out" -V M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET.sh

custom_rh="\"[[1,10,10],[10,1,1]]\""
custom_rv="\"[[1,10]]\""
dist="uniform"
name="zigzag_2_2"
qsub -S /bin/bash -cwd -N "$name" -q jdubi.q -o "$f/$name.out" -V M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET.sh

