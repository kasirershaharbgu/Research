#!/bin/bash

flags="-pe shared 14 -V -S /bin/bash -q jdubi.q -cwd"
vmin=0
vr=0
vu=0.1
vd=-0.1
vg_avg=0
c_avg=1
cg_std=0
r_avg=10
rg_std=0
repeats=50
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
rg_avg=500
cg_avg=10
vg_std=0
dist="exp"
graph=""
currentmap="--current-map"
full="--full"
resume=""
vSym="--symmetric-v"
leaping=""
efermi=""
dbg=""
calc_it=""
custom_rh="\"[[1,5,8,12,15,19],[19,15,12,8,5,1]]\""
custom_rv="\"[[1,5,10,15,19],[5,19,10,15,1],[19,15,1,10,5]]\""
custom_ch="\"[[0.1,0.5,0.8,1.2,1.5,1.9],[1.9,1.5,1.2,0.8,0.5.0.1]]\""
custom_cv="\"[[0.1,1,1.5,1.9,0.5],[0.1,1,1.5,1.9,0.5],[0.1,1,1.5,1.9,0.5]\""
rows=2
columns=5
up_electrode="\"[0,0,1,0,0]\""
down_electrode="\"[0,0,1,0,0]\""
T=0.001
vmax=3
efermi=""
f="bgu_2d_small_array_with_perp_current"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
c_std=2
r_std=20
for cg_avg in 5 10 30
do	
	name="array_${rows}_${columns}_cg_${cg_avg}"
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" run_RC_SET_with_perp.sh
sleep 1
done


Tmin=0.001
Tmax=0.1
Tstep=0.001
vl=0.5
vr=-0.5
vu=0.1
vd=-0.1
f="bgu_2d_small_array_with_perp_current_it"
calc_it="--calc-it"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
for cg_avg in 5 10 30
do	
	name="array_${rows}_${columns}_cg_${cg_avg}"
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vl",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" run_RC_SET_with_perp.sh
sleep 1
done
