#!/bin/sh
#SBATCH --mem=12G
#SBATCH --qos=normal
#SBATCH --partition=common
#SBATCH -n 32
#SBATCH -o /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/dedup_reads.out.txt -e /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/dedup_reads.err.txt

module load samtools pigz bzip2
module load graalvm/ce-java11-20.0.0
module load bbmap


################# input ###################

out_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/
master_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/scripts/

#sample ID
sample=$1


################ code #####################@

echo "deduplicating reads"

for i in $(ls "$out_dir"/FastQ/"$sample"/ | grep "MK" | sed 's/_R/ /' | awk '{print $1}' | sort -u)
do

clumpify.sh dedupe=t optical=f \
	in1="$out_dir"/FastQ/"$sample"/"$i"_R1.fq.gz \
	in2="$out_dir"/FastQ/"$sample"/"$i"_R2.fq.gz \
	out1="$out_dir"/FastQ/"$sample"/"$i"_dedup_R1.fq.gz \
	out2="$out_dir"/FastQ/"$sample"/"$i"_dedup_R2.fq.gz

	#rm "$out_dir"/FastQ/"$sample"/"$i"_R1.fq.gz
	#rm "$out_dir"/FastQ/"$sample"/"$i"_R2.fq.gz

done


