#!/bin/bash

vmin=0
vr=0
vg_avg=0
vg_std=0
c_avg=1
c_std=0
cg_std=0
r_avg=1
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
rg_avg=1000
cg_avg=10
t="20:00:00"
dist="uniform"
r_std=0
graph=""
currentmap="--current-map"
resume=""
full="--full"
vmax=2
f="2_vertical_dots_arrays_full_2"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi

rows=2
columns=1

name="increasing_r_one_big_one_small_c_big_connect_c_small_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[1]]\""
custom_ch="\"[[1,1],[0.5,0.5]]\""
custom_cv="\"[[0.5]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_one_big_one_small_c_big_connect_c_big_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[10]]\""
custom_ch="\"[[1,1],[0.5,0.5]]\""
custom_cv="\"[[1]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_one_big_one_small_c_small_connect_c_big_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[10]]\""
custom_ch="\"[[1,1],[0.5,0.5]]\""
custom_cv="\"[[0.5]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_one_big_one_small_c_small_connect_c_small_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[1]]\""
custom_ch="\"[[1,1],[0.5,0.5]]\""
custom_cv="\"[[0.5]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_const_c_big_connect_c_small_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[1]]\""
custom_ch="\"[[0.5,0.5],[0.5,0.5]]\""
custom_cv="\"[[1]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_const_c_big_connect_c_big_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[10]]\""
custom_ch="\"[[0.5,0.5],[0.5,0.5]]\""
custom_cv="\"[[1]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_const_c_small_connect_c_big_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[10]]\""
custom_ch="\"[[0.5,0.5],[0.5,0.5]]\""
custom_cv="\"[[0.5]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_const_c_small_connect_c_small_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[1]]\""
custom_ch="\"[[0.5,0.5],[0.5,0.5]]\""
custom_cv="\"[[0.5]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_zigzag_c_big_connect_c_small_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[1]]\""
custom_ch="\"[[0.5,1],[1,0.5]]\""
custom_cv="\"[[1]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_zigzag_c_big_connect_c_big_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[10]]\""
custom_ch="\"[[0.5,1],[1,0.5]]\""
custom_cv="\"[[1]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_zigzag_c_small_connect_c_big_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[10]]\""
custom_ch="\"[[0.5,1],[1,0.5]]\""
custom_cv="\"[[0.5]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh

name="increasing_r_zigzag_c_small_connect_c_small_connect_r_2_1"
custom_rh="\"[[1,10],[1,10]]\""
custom_rv="\"[[1]]\""
custom_ch="\"[[0.5,1],[1,0.5]]\""
custom_cv="\"[[0.5]]\""
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=2000MB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume" run_RC_SET_on_slurm.sh