#!/bin/sh
#SBATCH --qos=fast
#SBATCH --partition=common,dedicated
#SBATCH -n 5
#SBATCH -N 1
#SBATCH --mem=24G
#SBATCH -o resfinder.out.txt -e resfinder.err.txt


module load blast+/2.10.0 kma/1.3.4
module load resfinder/4.0.1 

#path to the assembly file (Fasta)
assembly=$1

#path to the output folder
out_path=$2

run_resfinder.py -ifa "$assembly" -o "$out_path"  -l 0.6 -t 0.9 --acquired

