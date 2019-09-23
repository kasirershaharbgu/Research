#!/bin/sh

#SBATCH -N 1                               # nodes=1
#SBATCH --ntasks-per-node=1                # ppn=1
#SBATCH -J "${21}"                         # job name
#SBATCH -t "${23}"                         # jobtime hh:mm:ss
#SBATCH -p dept                            # partition/queue name
#SBATCH --mem=2000MB                       # memory in MB
#SBATCH --output="{21}.out"                # file for STDOUT
#SBATCH --mail-user=skasirer@princeton.edu # Mail  id of the user
#SBATCH --mail-type=end                    # Slurm will send at the completion of your job

python3 dots_array.py -M "$1" -N "$2" --vmin "$3" --vmax "$4" --vstep "$5" --vg-avg "$6" --vg-std "$7" --c-avg "$8" --c-std "$9" --cg-avg "${10}" --cg-std "${11}" --r-avg "${12}" --r-std "${13}" --rg-avg "${14}" --rg-std "${15}" --repeats "${16}" --n-avg "${17}" --n-std "${18}" --q-avg "${19}" --q-std "${20}" --file-name "${21}" -o "${22}"