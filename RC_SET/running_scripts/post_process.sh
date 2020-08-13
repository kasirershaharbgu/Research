#!/bin/bash

for file in $@
do
	echo "processing $file"
	name=(`echo "$file" | sed "s/.*\///"`)
	directory=(`echo "$file" | sed "s!\(.*\)/.*!\1!"`)
	python3 dots_array.py --file-name "$name" --out "$directory" --plot-current-map --full
	python3 dots_array.py --file-name "$name" --out "$directory" --plot-binary-current-map --full
done

