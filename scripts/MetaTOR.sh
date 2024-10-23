#!/bin/sh
#SBATCH --qos=normal
#SBATCH --partition=common
#SBATCH --mem=56G
#SBATCH -n 48
#SBATCH -N 1
#SBATCH -o /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/MetaTOR.out.txt -e /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/MetaTOR.err.txt

module load hmmer/3.3 pplacer/v1.1-alpha19 prodigal/2.6.3 
module load CheckM/1.1.3 
module add samtools/1.10 
module add bedtools/2.29.2 
module add bowtie2/2.3.5.1 
module load gcc
module add Python/3.8.3 

project=$1
assembly=$2
reads_for=$3
reads_rev=$4

export LOUVAIN_PATH=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/gen-louvain/

mkdir -p "$project"/temp/

metator pipeline -F -N --threads 48 \
	--assembly "$assembly" \
	--forward "$reads_for" \
	--reverse "$reads_rev"  \
	--tmpdir "$project"/temp/ \
	--outdir "$project"

rm -r "$project"/temp/
rm "$project"/*.bam
