#!/bin/bash

qsub -V -S /bin/bash -cwd -N dbg -q jdubi.q runToyModel.sh