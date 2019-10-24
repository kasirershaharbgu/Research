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
out="2D_square_arrays_with_current_map_bgu"
full="\"\"" 
graph="\"\""
currentmap="--current-map"
M=5
N=5
vmax=5
custom_rh="\"\""
custom_rv="\"\""
if [ ! -d "$out" ]; then
  mkdir "$out"
fi
dist="two_points"
file_name="no_disorder_5_5"
r_std=0
qsub -V -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$out/$file_name.out" run_RC_SET.sh

dist="uniform"
name="uniform_disorder_5_5"
r_std=9
qsub -V -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$out/$file_name.out" run_RC_SET.sh

dist="two_points"
file_name="two_points_disorder_5_5"
r_std=9
qsub -V -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$out/$file_name.out" run_RC_SET.sh

dist="exp"
file_name="exp_disorder_5_5"
r_std=9
qsub -V -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$out/$file_name.out" run_RC_SET.sh

N=2
M=2
custom_rh="\"[[10,10,10],[1,1,1]]\""
custom_rv="\"[[10,10]]\""
dist="uniform"
file_name="one_nig_one_small_2_2"
qsub -V -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$out/$file_name.out" run_RC_SET.sh

custom_rh="\"[[1,10,10],[10,1,1]]\""
custom_rv="\"[[1,10]]\""
dist="uniform"
file_name="zigzag_2_2"
qsub -V -S /bin/bash -cwd -N "$file_name" -q jdubi.q -o "$out/$file_name.out" run_RC_SET.sh

