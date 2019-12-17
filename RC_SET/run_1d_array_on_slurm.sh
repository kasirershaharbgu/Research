#!/bin/bash

vmin=0
vr=0
vg_avg=0
vg_std=0
c_avg=1
c_std=0
cg_std=1
r_avg=10
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.05
rg_avg=1000
cg_avg=10
t="1:00:00"
dist="exp"
r_std=0
graph=""
currentmap="--current-map"
full="--full"
resume=""
vSym="--symmetric-v"
resume=""
leaping=""
efermi=""
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
T=0.01
vmax=5
f="1d_array"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
rows=1
columns=10
for run in 1 2 3
do
	efermi=""
	r_std=0
	c_std=0
	name="1_10_shorter_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=3GB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh
	name="1_10_shorter_run_${run}"
	
	
	efermi="--variable-ef"
	r_std=0
	c_std=0
	name="1_10_shorter_var_ef_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=3GB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh
	
	efermi=""
	r_std=5
	c_std=0
	name="1_10_shorter_r_disorder_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=3GB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh
	name="1_10_shorter_run_${run}"
	
	r_std=5
	c_std=0
	efermi="--variable-ef"
	name="1_10_shorter_var_ef_r_disorder_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=3GB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh
	
	efermi=""
	c_std=1
	r_std=0
	name="1_10_shorter_c_disorder_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=3GB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh
	name="1_10_shorter_run_${run}"
	
	c_std=1
	r_std=0
	efermi="--variable-ef"
	name="1_10_shorter_var_ef_c_disorder_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=3GB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh
	
	efermi=""
	c_std=1
	r_std=5
	name="1_10_shorter_cr_disorder_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=3GB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh
	name="1_10_shorter_run_${run}"
	
	c_std=1
	r_std=5
	efermi="--variable-ef"
	name="1_10_shorter_var_ef_cr_disorder_run_${run}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=3GB --time="$t" --partition=dept --mail-user=skasirer@princeton.edu --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi" run_RC_SET_on_slurm.sh
	
done	
