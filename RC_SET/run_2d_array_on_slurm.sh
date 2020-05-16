#!/bin/bash

vmin=0
vr=0
vg_avg=0
c_avg=1
cg_std=0
r_avg=1
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
rg_avg=500
cg_avg=10
t="5-12:00:00"
dist="exp"
graph=""
currentmap="--current-map"
full="--full"
resume=""
vSym=""
leaping=""
efermi=""
dbg=""
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
T=0
vmax=4
efermi=""
f="2d_array_slurm_custom_paths"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
columns=10
rows=10
vg_std=0
cg_std=0
c_std=0
r_std=0

name="array_${rows}_${columns}_r_block"
custom_rh="\"[[1,1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,1,10]]\""
custom_rv="\"[[1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,10],[1,1,1,1,1,1,1,1,1,10]]\""
custom_ch="\"\""
custom_cv="\"\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=10GB --time="$t" --partition=dept --mail-user=kasirer@post.bgu.ac.il --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",inpu"$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh

f="2d_array_slurm_big_c_disorder"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
c_avg=10
cg_std=9.5
r_avg=10
rg_std=9

for run in 1 2 3 4 5
do

	name="array_${rows}_${columns}_c_r_disorder_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=10GB --time="$t" --partition=dept --mail-user=kasirer@post.bgu.ac.il --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",inpu"$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh

done	
