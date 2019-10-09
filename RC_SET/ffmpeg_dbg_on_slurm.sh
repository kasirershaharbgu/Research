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
repeats=1
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.1
rg_avg=1000
cg_avg=10
f="dbg"
full="" 
currentmap="--current-map"
t="00:10:00"
rows=2
columns=2
vmax=1
custom_rh="\"\""
custom_rv="\"\""
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
dist="two_points"
name="current_mp_dbg"
r_std=0

sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=1000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv" run_RC_SET_on_slurm.sh