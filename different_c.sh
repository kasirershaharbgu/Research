#!/bin/bash

arrayC=( 0.0001 0.001 0.01 0.1 )
arrayR=( 1 )

for R in "${arrayR[@]}"
do

for C in "${arrayC[@]}"
do

qsub -V -S /bin/bash -cwd -N different_c -q jdubi.q runToyModel.sh "$R" "$C" "C_$C" "$1"

sleep 10

done
done