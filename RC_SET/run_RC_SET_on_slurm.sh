#!/bin/sh

#SBATCH -N 1                               # nodes=1
#SBATCH --ntasks-per-node=1                # ppn=1
#SBATCH -p dept                            # partition/queue name
#SBATCH --mem=2000MB                       # memory in MB
#SBATCH --mail-user=skasirer@princeton.edu # Mail  id of the user
#SBATCH --mail-type=end                    # Slurm will send at the completion of your job

. /etc/local/set-anaconda3
python3 dots_array.py -M "$M" -N "$N" --vmin "$vmin" --vmax "$vmax" --vstep "$vstep" --vg-avg "$vg_avg" --vg-std "$v_std" --c-avg "$c_avg" --c-std "$c_std" --cg-avg "$cg_avg" --cg-std "$cg_std" --r-avg "$r_avg" --r-std "$r_std" --rg-avg "$rg_avg" --rg-std "$rg_std" --repeats "$repeats" --n-avg "$n_avg" --n-std "$n_std" --q-avg "$q_avg" --q-std "$q_std" --file-name "$file_name" -o "$out"