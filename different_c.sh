#!/bin/bash

arrayC=( 0.001 0.01 0.1 1 )
arrayR=( 1 )

for R in "${arrayR[@]}"
do

for C in "${arrayC[@]}"
do

qsub -V -S /bin/bash -cwd -N dbg -q jdubi.q runToyModel.sh "$R" "$C"

sleep 10

done
done