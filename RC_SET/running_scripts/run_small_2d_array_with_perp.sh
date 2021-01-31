#!/bin/bash

flags="-pe shared 14 -V -S /bin/bash -q jdubi.q -cwd"
vmin=0
vr=0
vu=0.1
vd=-0.1
vg_avg=0
c_avg=1
c_std=0
cg_std=0
r_avg=10
r_std=0
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
resume="--resume"
vSym="--symmetric-v"
leaping=""
efermi=""
dbg=""
calc_it=""
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
rows=3
columns=7
up_electrode="\"[0,0,0,1,0,0,0]\""
down_electrode="\"[0,0,0,1,0,0,0]\""
T=0.001
vmax=2
efermi=""
vr=0
vu=0.1
vd=-0.1
input=""

f="2d_small_array_bgu_with_perp"
calc_it=""
if [ ! -d "$f" ]; then
  mkdir "$f"
fi

name="array_3_7_r_c_std_cg_1"
cg_avg=1
r_std=9
c_std=0.9
#qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1

name="array_3_7_r_c_std_cg_5"
cg_avg=5
r_std=9
c_std=0.9
#qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1

name="array_no_disorder_cg_1"
cg_avg=1
r_std=0
c_std=0
#qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1


name="array_no_disorder_cg_5"
cg_avg=5
r_std=0
c_std=0
#qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1

#IT runs

calc_it="--calc-it"
Tmin=0
Tmax=0.1
Tstep=0.001


name="array_3_7_r_c_std_cg_1"
input="${f}/runningParameters_${name}.txt"
vmin=1
qsub "$flags" -N "${name}_it" -o "$f/${name_it}.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$Tmin",file_name="${name}_it",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1

name="array_3_7_r_c_std_cg_5"
input="${f}/runningParameters_${name}.txt"
vmin=1.5
qsub "$flags" -N "${name}_it" -o "$f/${name_it}.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$Tmin",file_name="${name}_it",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1

name="array_no_disorder_cg_1"
input="${f}/runningParameters_${name}.txt"
vmin=0.25
qsub "$flags" -N "${name}_it" -o "$f/${name_it}.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$Tmin",file_name="${name}_it",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1


name="array_no_disorder_cg_5"
input="${f}/runningParameters_${name}.txt"
vmin=1
qsub "$flags" -N "${name}_it" -o "$f/${name_it}.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$Tmin",file_name="${name}_it",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v down_electrode="$down_electrode" -v up_electrode="$up_electrode" running_scripts/run_RC_SET_with_perp.sh
sleep 1

