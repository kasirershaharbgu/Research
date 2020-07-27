

flags="-pe shared 14 -V -S /bin/bash -q jdubi.q -cwd"
vg_avg=0
c_avg=1
cg_std=0
r_avg=10
rg_std=0
repeats=20
n_avg=0
n_std=0
q_avg=0
q_std=0
vstep=0.001
rg_avg=500
cg_avg=10
dist="exp"
graph=""
currentmap="--current-map"
full="--full"
resume=""
vSym="--symmetric-v"
leaping=""
efermi=""
dbg=""
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
IT="--calc-it"
T=0.001
vmax=0.1
efermi=""
f="2d_long_array_bgu_with_perp_point_contact_IT"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
columns=15
rows=5
vg_std=0
cg_std=0
c_std=1
r_std=9
for run in 2 3 6 7 9 10
do	

input="2d_long_array_bgu_with_perp_point_contact/runningParameters_array_5_15_r_std_9_run_${run}.txt"
name="array_${rows}_${columns}_r_std_${r_std}_run_${run}_horizontal_res"
vmin=2
vr=-1
vu=0
vd=0
up_electrode="\"[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]\""
down_electrode="\"[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]\""
left_electrode="\"[1,1,1,1,1]\""
right_electrode="\"[1,1,1,1,1]\""
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",input=$"input",leaping="$leaping",efermi="$efermi",dbg="$dbg" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v up_electrode="$down_electrode" -v down_electrode="$up_electrode" run_RC_SET_with_perp.sh
sleep 1
name="array_${rows}_${columns}_r_std_${r_std}_run_${run}_vertical_res"
vmin=0
vr=0
vu=1
vd=-1
up_electrode="\"[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]\""
down_electrode="\"[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0]\""
left_electrode="\"[0,0,0,0,0]\""
right_electrode="\"[0,0,0,0,0]\""
qsub "$flags" -N "$name" -o "$f/$name.out" -v M="$rows",N="$columns",vmin="$vmin",vmax="$vmax",vstep="$vstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",vr="$vr",dist="$dist",full="$full",IT="$IT",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",input=$"input",leaping="$leaping",efermi="$efermi",dbg="$dbg" -v custom_rh="$custom_rh" -v custom_rv="$custom_rv" -v custom_ch="$custom_ch" -v custom_cv="$custom_cv" -v up_electrode="$down_electrode" -v down_electrode="$up_electrode" run_RC_SET_with_perp.sh
sleep 1
done
