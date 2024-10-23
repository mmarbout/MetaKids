#!/bin/sh
#SBATCH --partition=depgg
#SBATCH --account=depgg
#SBATCH --mem=256G
#SBATCH -n 96
#SBATCH -N 1
#SBATCH -o /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/Anvio_init.out.txt -e /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/Anvio_init.err.txt

# Martial Marbouty

#path to the assembly file
assembly=$1

#ID of the studied kids
sample=$2

#folder containing the needed databases
db=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/database/

#out_fold
out_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/

module load R
module load hmmer prodigal muscle blast+ mcl muscle diamond
module load VirSorter/2.2.3 
module load diamond/2.0.4 hmmer/3.3 prodigal/2.6.3 prodigal-gv/
module load CheckV/ 
module load infernal
module add blast+ SQLite/3.31.0
module add R/3.6.2 prodigal mcl/14-137 muscle/3.8.31 hmmer/3.2.1 diamond/2.0.4 samtools/1.10 FastTree/2.1.11 centrifuge/1.0.4-beta SQLite/3.31.0 tRNAscan-SE/2.0.7 blast+/2.10.0 megahit/1.2.9 SPAdes/3.15.0 bowtie2/2.3.5.1 bwa/0.7.17 trimal/1.4.1 IQ-TREE/2.0.6
module load anvio/
module add MMseqs2/10-6d92c
module add bedtools/2.29.2 
module load gcc
module load perl/5.30.1

mkdir -p "$out_fold"/Anvio/

anvi-gen-contigs-database -L 0 -T 96 --project-name "$sample" -f "$assembly" -o "$out_dir"/Anvio/"$sample".db

anvi-run-hmms -T 96 -c "$out_dir"/Anvio/"$sample".db

anvi-run-ncbi-cogs -T 96 --cog-version COG20 --cog-data-dir "$db"/COG_2020 -c "$out_dir"/Anvio/"$sample".db

anvi-run-pfams -T 96 --pfam-data-dir "$db"/Pfam_v32 -c "$out_dir"/Anvio/"$sample".db

anvi-export-gene-calls --gene-caller prodigal -c "$out_dir"/Anvio/"$sample".db -o "$out_dir"/Anvio/"$sample"-gene-calls.txt









