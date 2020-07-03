#!/bin/sh

echo "
python3 thermalArray.py -M "$M" -N "$N" -T "$T" -R "$R" --alpha "$alpha" --vmin "$vmin" --vmax "$vmax" --vstep "$vstep" --vr "$vr" --vu "$vu" --vd "$vd" --distribution "$dist" -o "$out" --file-name "$name" --n-avg "$navg" --n-std "$nstd" "$vsym"
"

. /etc/local/set-anaconda3
python3 thermalArray.py -M "$M" -N "$N" -T "$T" -R "$R" --alpha "$alpha" --vmin "$vmin" --vmax "$vmax" --vstep "$vstep" --vr "$vr" --vu "$vu" --vd "$vd" --distribution "$dist" -o "$out" --file-name "$name" --n-avg "$navg" --n-std "$nstd" "$vsym"
