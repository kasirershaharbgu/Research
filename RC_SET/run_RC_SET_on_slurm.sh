#!/bin/sh

echo "
python3 dots_array.py -M "$M" -N "$N" --vmin "$vmin" --vmax "$vmax" --vstep "$vstep" --vr "$vr" --vg-avg "$vg_avg" --vg-std "$vg_std" --c-avg "$c_avg" --c-std "$c_std" --cg-avg "$cg_avg" --cg-std "$cg_std" --r-avg "$r_avg" --r-std "$r_std" --rg-avg "$rg_avg" --rg-std "$rg_std" --repeats "$repeats" --n-avg "$n_avg" --n-std "$n_std" --q-avg "$q_avg" --q-std "$q_std" -T "$T" --distribution "$dist" --custom-rh "$custom_rh" --custom-rv "$custom_rv" --custom-ch "$custom_ch" --custom-cv "$custom_cv" --file-name "$file_name" -o "$out" "$full" "$currentmap" "$graph" "$resume" "$vSym" "$leaping" "$efermi"
"

. /etc/local/set-anaconda3
python3 dots_array.py -M "$M" -N "$N" --vmin "$vmin" --vmax "$vmax" --vstep "$vstep" --vr "$vr" --vg-avg "$vg_avg" --vg-std "$vg_std" --c-avg "$c_avg" --c-std "$c_std" --cg-avg "$cg_avg" --cg-std "$cg_std" --r-avg "$r_avg" --r-std "$r_std" --rg-avg "$rg_avg" --rg-std "$rg_std" --repeats "$repeats" --n-avg "$n_avg" --n-std "$n_std" --q-avg "$q_avg" --q-std "$q_std" -T "$T" --distribution "$dist" --custom-rh "$custom_rh" --custom-rv "$custom_rv" --custom-ch "$custom_ch" --custom-cv "$custom_cv" --file-name "$file_name" -o "$out" "$full" "$currentmap" "$graph" "$resume" "$vSym" "$leaping" "$efermi"
