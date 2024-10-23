#!/bin/sh
#SBATCH --mem=12G
#SBATCH --qos=normal
#SBATCH --partition=common
#SBATCH -n 32
#SBATCH -o /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/reads_cleaning_HiC.out.txt -e /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/reads_cleaning_HiC.err.txt

module load graalvm/ce-java19-22.3.1
module load Trimmomatic/0.39
module load fastqc/0.11.9
module load fastx_toolkit/0.0.14
module add samtools/1.10 
module add bedtools/2.29.2 
module add bowtie2/2.3.5.1 
module load gcc
module add Python/3.8.3


################# input ###################

#path to the folder where the cleaned reads will be stored
out_fold=$1

#path to the input forward reads
reads_for=$2

#path to the input reverse reads
reads_rev=$3

#name of your project: i.e. it will be the name of your cleaned reads --> project_cleaned_R1.fq.gz & project_cleaned_R2.fq.gz
project=$4


out_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/

################ code #####################@

echo "cleaning reads"

Trimmomatic PE -threads 32 "$reads_for" "$reads_rev" \
	"$out_fold"/"$project"_cleaned_R1.fq.gz \
	"$out_fold"/"$project"_unpaired_R1.fq.gz \
	"$out_fold"/"$project"_cleaned_R2.fq.gz \
	"$out_fold"/"$project"_unpaired_R2.fq.gz \
	ILLUMINACLIP:/pasteur/zeus/projets/p02/rsg_fast/mmarbout/Fasta/divers/TruSeq3_PE_adapt.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36


echo "checking quality"

mkdir -p "$out_dir"/FastQ/rapport/

fastqc -t 32 --nogroup -o "$out_dir"/FastQ/rapport/raw/ "$reads_for"
fastqc -t 32 --nogroup -o "$out_dir"/FastQ/rapport/raw/ "$reads_rev"

fastqc -t 32 --nogroup -o "$out_dir"/FastQ/rapport/HiC/ "$out_fold"/"$project"_cleaned_R1.fq.gz
fastqc -t 32 --nogroup -o "$out_dir"/FastQ/rapport/HiC/ "$out_fold"/"$project"_cleaned_R2.fq.gz

rm "$out_fold"/"$project"_unpaired_R1.fq.gz
rm "$out_fold"/"$project"_unpaired_R2.fq.gz

hicstuff cutsite --mode all -t 32 -e DpnII,HinfI \
	-1 "$out_fold"/"$project"_cleaned_R1.fq.gz \
	-2 "$out_fold"/"$project"_cleaned_R2.fq.gz \
	-p "$out_dir"/FastQ/HiC_digested/"$project"

rm "$out_fold"/"$project"_cleaned_R1.fq.gz
rm "$out_fold"/"$project"_cleaned_R2.fq.gz

rm "$reads_for"
rm "$reads_rev"
