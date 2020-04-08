#!/bin/bash

for file in $@
do
	echo "processing $file"
	name=(`echo "$file" | sed "s/.*\///"`)
	directory=(`echo "$file" | sed "s!\(.*\)/.*!\1!"`)
	/c/Users/shahar/AppData/Local/Programs/Python/Python37/python.exe dots_array.py --file-name "$name" --out "$directory" --plot-current-map --full
	/c/Users/shahar/AppData/Local/Programs/Python/Python37/python.exe dots_array.py --file-name "$name" --out "$directory" --plot-binary-current-map --full
done

