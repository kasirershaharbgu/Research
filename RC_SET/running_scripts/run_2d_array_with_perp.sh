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
rg_avg=1000
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
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
rows=2
columns=5
up_electrode="\"[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]\""
down_electrode="\"[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]\""
T=0.001
vmax=10
efermi=""
f="2d_long_array_bgu_with_perp_point_contact"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
c_std=0
r_std=0
input_base="2d_long_array_bgu_with_perp_point_contact/runningParameters_array_5_15_r_std_9_run_"
for run in 1 2 3 4
do      
        name="array_5_15_r_std_9_run_${run}"
        input="${input_base}${run}.txt"
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vl",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1
done
name="array_5_15_no_disorder"
input=""
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vl",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1

Tmin=0.001
Tmax=0.1
Tstep=0.001
vl=0.5
vr=-0.5
vu=0.1
vd=-0.1

f="2d_long_array_bgu_with_perp_point_contact_it"
calc_it="--calc-it"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
for run in 1 2 3 4
do      	
	name="array_5_15_r_std_9_run_${run}_it"
        input="${input_base}${run}.txt"
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vl",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1
done
name="array_5_15_no_disorder_it"
input=""
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vl",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1
