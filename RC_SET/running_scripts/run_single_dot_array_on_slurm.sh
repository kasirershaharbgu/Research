#!/bin/bash

vmin=-2
vr=0
vg_avg=0
vg_std=0
c_avg=1
c_std=0.5
cg_std=0
r_avg=10
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
vmax=2
dist="uniform"
rows=1
columns=1
t="00:10:00"
f="single_c_disorder"
rg_avg=1000
full="--full" 
currentmap="\"\""
r_std=9
if [ ! -d "$f" ]; then
  mkdir "$f"
fi


cg_avg=10
custom_rh="\"[[1,10]]\""
custom_rv="\"[[]]\""
name="r_disorder_big_cg"

sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=5000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET_on_slurm.sh

cg_avg=1
name="r_disorder_small_cg"

sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=5000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET_on_slurm.sh

cg_avg=10
custom_rh="\"[[1,1]]\""
custom_rv="\"[[]]\""
name="no_disorder_big_cg"

sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=5000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET_on_slurm.sh

cg_avg=1
name="no_disorder_small_cg"

sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=5000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",graph="$graph" run_RC_SET_on_slurm.sh

