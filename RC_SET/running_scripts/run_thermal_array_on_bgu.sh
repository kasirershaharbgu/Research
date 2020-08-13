#!/bin/bash

flags="-pe shared 14 -V -S /bin/bash -q jdubi.q -cwd"
M= 
N= 
T=
R= 
alpha= 
vmin= 
vmax= 
vstep= 
vr= 
vu= 
vd= 
dist= 
out= 
name= 
navg= 
nstd= 
vsym=

qsub "$flags" -N "$name" -o "$out/$name.out" -v M="$M",N="$N",T="$T",R="$R",alpha="$alpha",vmin="$vmin",vmax="$vmax",vstep="$vstep",vr="$vr",vu="$vu",vd="$vd",dist="$dist",out="$out",name="$name",navg="$navg",nstd="$nstd",vsym="$vsym" run_thermal_array.sh
sleep 1
