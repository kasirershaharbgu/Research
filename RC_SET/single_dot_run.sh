#!/bin/bash

arrayVG=( 0 1 10 )
arrayCG=( 0.1 1 10 )
arrayRG=( 0.1 1 10 )
arrayCR=( 1 3 )
arrayRR=( 1 3 )
for VG in "${arrayVG[@]}"
do

for CG in "${arrayCG[@]}"
do

for RG in "${arrayRG[@]}"
do

for CR in "${arrayCR[@]}"
do

for RR in "${arrayRR[@]}"
do

qsub -V -S /bin/bash -cwd -N different_c -q jdubi.q run_RC_SET.sh 1 1 0 5 0.1 100 "$VG" "$CR" 1 "$CG" "$RR" 1 "$RG" 10 0 0 0.01 "CR_$CR_RR_$RR" "single_dot_VG_$VG_CG_$CG_RG_$RG" 

sleep 10

done
done
done
done
done