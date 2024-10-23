#!/bin/sh
#SBATCH --mem=12G
#SBATCH --qos=normal
#SBATCH --partition=common
#SBATCH -n 32
#SBATCH -o /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/reads_cleaning_SG.out.txt -e /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/reads_cleaning_SG.err.txt

module load graalvm/ce-java19-22.3.1
module load Trimmomatic/0.39
module load fastqc/0.11.9
module load fastx_toolkit/0.0.14

################# input ###################

out_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/
master_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/scripts/

#sample ID
sample=$1

################ code #####################@

echo "cleaning reads"

for i in $(ls "$out_dir"/FastQ/"$sample"/ | grep "MK" | grep "Shot" | sed 's/_R/ /' | awk '{print $1}' | sort -u)
do

	Trimmomatic PE -threads 32 "$out_dir"/FastQ/"$sample"/"$i"_dedup_R1.fq.gz "$out_dir"/FastQ/"$sample"/"$i"_dedup_R2.fq.gz \
		"$out_dir"/FastQ/SG_cleaned/"$i"_R1.fq.gz \
		"$out_dir"/FastQ/SG_cleaned/"$i"_unpaired_R1.fq.gz \
		"$out_dir"/FastQ/SG_cleaned/"$i"_R2.fq.gz \
		"$out_dir"/FastQ/SG_cleaned/"$i"_unpaired_R2.fq.gz \
		ILLUMINACLIP:/pasteur/zeus/projets/p02/rsg_fast/mmarbout/Fasta/divers/TruSeq3_PE_adapt.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36


echo "checking quality"

mkdir -p "$out_dir"/rapport/

fastqc -t 32 --nogroup -o "$out_dir"/FastQ/rapport/raw/ "$out_dir"/FastQ/"$sample"/"$i"_R1.fq.gz
fastqc -t 32 --nogroup -o "$out_dir"/FastQ/rapport/raw/ "$out_dir"/FastQ/"$sample"/"$i"_R2.fq.gz

fastqc -t 32 --nogroup -o "$out_dir"/FastQ/rapport/SG/ "$out_dir"/FastQ/SG_cleaned/"$i"_R1.fq.gz
fastqc -t 32 --nogroup -o "$out_dir"/FastQ/rapport/SG/ "$out_dir"/FastQ/SG_cleaned/"$i"_R2.fq.gz

rm "$out_dir"/FastQ/SG_cleaned/"$i"_unpaired_R1.fq.gz
rm "$out_dir"/FastQ/SG_cleaned/"$i"_unpaired_R2.fq.gz

done


sbatch "$master_dir"/Megahit.sh \
	$(ls "$out_dir"/FastQ/SG_cleaned/ | grep "MK" | grep "$sample" | grep "Shot" | sed 's/_R/ /' | awk '{print $1}' | sort -u | awk '{print "'$out_dir'""/FastQ/SG_cleaned/"$1"_R1.fq.gz"}' | paste -s | sed 's/\t/,/g')\
	$(ls "$out_dir"/FastQ/SG_cleaned/ | grep "MK" | grep "$sample" | grep "Shot" | sed 's/_R/ /' | awk '{print $1}' | sort -u | awk '{print "'$out_dir'""/FastQ/SG_cleaned/"$1"_R2.fq.gz"}' | paste -s | sed 's/\t/,/g') \
	"$sample"

rm "$out_dir"/FastQ/"$sample"/"$i"_dedup_R1.fq.gz
rm "$out_dir"/FastQ/"$sample"/"$i"_dedup_R2.fq.gz
