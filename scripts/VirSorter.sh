#!/bin/sh
#SBATCH --mem=24G
#SBATCH --qos=normal
#SBATCH --partition=common
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -o /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/VirSorter.out.txt -e /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/VirSorter.err.txt

module load R
module load hmmer/3.3 blast+/2.10.0 mcl/14-137 muscle/3.8.1551 diamond/0.9.30 MetaGeneAnnotator/2008-8-19
module load VirSorter/1.0.6


#out_fold
out_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/

#ID of the studied kids
sample=$1

wrapper_phage_contigs_sorter_iPlant.pl -f "$out_dir"/assemblage/"$sample".fa --db 2 --wdir "$out_dir"/annotations/VirSorter/"$sample"/ -ncpu 32
