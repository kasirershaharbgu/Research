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
rg_avg=500
t="3-12:00:00"
dist="exp"
graph=""
currentmap="--current-map"
full="--full"
resume=""
vSym="--symmetric-v"
leaping=""
efermi=""
dbg=""
double_time=""
double_loop=""
superconducting=""
gap=0
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
T=0.001
vmax=3
efermi=""
f="bgu_it"
input="bgu_it/runningParameters_array_10_10_disorder_c_std_0.1_r_std_9_run_1.txt"
if [ ! -d "$f" ]; then
  mkdir "$f"
fi
custom_rh="\"\""
custom_rv="\"\""
custom_ch="\"\""
custom_cv="\"\""
c_avg=1
c_std=0.9
r_std=9
rows=5
columns=5
Tmin=0.001
Tmax=0.1
Tstep=0.001
vr=0
for cg_avg in 1 5 10
do
	for vl in 0.2 0.6 0.8 1
	do
	name="array_${rows}_${columns}_it_v_${vl}"
	sbatch -J="$name" --nodes=1 --ntasks-per-node="$repeats" --mem=10GB --time="$t" --partition=dept --mail-user=kasirer@post.bgu.ac.il --mail-type=end --output="$f/$name.out"  --export=M="$rows",N="$columns",vmin="$vl",vmax="$Tmax",vstep="$Tstep",vd="$vd",vu="$vu",vg_avg="$vg_avg",vg_std="$vg_std",c_avg="$c_avg",c_std="$c_std",cg_avg="$cg_avg",cg_std="$cg_std",r_avg="$r_avg",r_std="$r_std",rg_avg="$rg_avg",rg_std="$rg_std",repeats="$repeats",n_avg="$n_avg",n_std="$n_std",q_avg="$q_avg",q_std="$q_std",T="$T",file_name="$name",out="$f",input="$input",vr="$vr",dist="$dist",full="$full",currentmap="$currentmap",graph="$graph",resume="$resume",vSym="$vSym",leaping="$leaping",efermi="$efermi",dbg="$dbg",calc_it="$calc_it",custom_rh="$custom_rh",custom_rv="$custom_rv" -v custom_ch="$custom_ch",custom_cv="$custom_cv" running_scripts/run_RC_SET_on_slurm.sh
	sleep 1
done	
done
