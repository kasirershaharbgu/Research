#!/bin/bash

vmin=-2
vr=0
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
dist="two_points"
r_std=0
graph="\"\""
full="--full"
currentmap="\"\""
vmax=3
out="2_dots_arrays_gillespie"
if [ ! -d "$out" ]; then
  mkdir "$out"
fi

N=1
M=2
vg_avg=1
file_name="increasing_r_vg_1"
custom_rh="\"[[1,10,20]]\""
custom_rv="\"[[]]\""
qsub -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$f/$name.out" -V run_RC_SET.sh

file_name="decreasing_r_vg_1"
custom_rh="\"[[20,10,1]]\""
custom_rv="\"[[]]\""
qsub -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$f/$name.out" -V run_RC_SET.sh

vg_avg=1
vg_std=1
file_name="increasing_r_noisy_vg"
custom_rh="\"[[1,10,20]]\""
custom_rv="\"[[]]\""
qsub -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$f/$name.out" -V run_RC_SET.sh

file_name="decreasing_r_noisy_vg"
custom_rh="\"[[20,10,1]]\""
custom_rv="\"[[]]\""
qsub -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$f/$name.out" -V run_RC_SET.sh