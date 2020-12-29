#!/bin/bash

vmin=0
vr=0
temperature_gradient=0
vg_avg=0
vg_std=0
c_avg=1
cg_avg=10
cg_std=0
r_avg=10
rg_std=0
repeats=25
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.01
rg_avg=1000
t="2-12:00:00"
dist="exp"
graph=""
currentmap="--current-map"
full="--full"
resume="--resume"
vSym="--symmetric-v"
leaping=""
efermi=""
dbg=""
double_time=""
double_loop=""
superconducting="--superconducting"
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
T=0.001
f="super_conducting_array_different_temperatures/"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
c_avg=1
c_std=2
r_std=20
rows=5
columns=5


run=1
input="super_conducting_array_different_temperatures/runningParameters_sc_array_5_5_T_0.001_cg_10_run_1.txt"
vmax=4
for gap in 0.01
do
name="sc_array_${rows}_${columns}_T_${T}_cg_${cg_avg}_run_${run}_gap_${gap}"
#sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=10GB --time="$t" --partition=dept --mail-user=kasirer@post.bgu.ac.il --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",double_time="$double_time",double_loop="$double_loop",superconducting="$superconducting",gap="$gap",input="$input",dbg="$dbg" running_scripts/run_RC_SET_on_slurm.sh
sleep 1
done

run=2
input="super_conducting_array_different_temperatures/runningParameters_sc_array_5_5_T_0.001_cg_10_run_2.txt"
vmax=5
for gap in 0.01
do
name="sc_array_${rows}_${columns}_T_${T}_cg_${cg_avg}_run_${run}_gap_${gap}"
#sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=10GB --time="$t" --partition=dept --mail-user=kasirer@post.bgu.ac.il --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",double_time="$double_time",double_loop="$double_loop",superconducting="$superconducting",gap="$gap",input="$input",dbg="$dbg" running_scripts/run_RC_SET_on_slurm.sh
sleep 1
done

run=3
input="super_conducting_array_different_temperatures/runningParameters_sc_array_5_5_T_0.001_cg_10_run_3.txt"
vmax=5.5
for gap in 0.01
do
name="sc_array_${rows}_${columns}_T_${T}_cg_${cg_avg}_run_${run}_gap_${gap}"
sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=10GB --time="$t" --partition=dept --mail-user=kasirer@post.bgu.ac.il --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",temperature_gradient="$temperature_gradient",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",custom_rh="$custom_rh",custom_rv="$custom_rv",custom_ch="$custom_ch",custom_cv="$custom_cv",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",double_time="$double_time",double_loop="$double_loop",superconducting="$superconducting",gap="$gap",input="$input",dbg="$dbg" running_scripts/run_RC_SET_on_slurm.sh
sleep 1
done
