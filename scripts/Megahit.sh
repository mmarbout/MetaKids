#!/bin/sh
#SBATCH --partition=depgg
#SBATCH -A depgg
#SBATCH --mem=248G
#SBATCH -n 48
#SBATCH -o /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/assembly_MEGAHIT.out.txt -e /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/log/assembly_MEGAHIT.err.txt

module load megahit/1.2.9

# Path to the global output folder and the folder containing the scripts
out_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/
master_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/scripts/

#path to the input forward reads
reads_for=$1

#path to the input reverse reads
reads_rev=$2

#sample ID
sample=$3

rm -r "$out_dir"/assemblage/"$sample"/ 

megahit -t 48 -m 248 -o "$out_dir"/assemblage/"$sample"/ -1 "$reads_for" -2 "$reads_rev"

#NB : you can feed megahit with several reads by adding a comma between them to perform a co-assembly of one kinetic: 
# reads1_R1.fq.gz,reads2_R1.fq.gz,reads3_R1.fq.gz  and reads1_R2.fq.gz,reads2_R2.fq.gz,reads3_R2.fq.gz

cat "$out_dir"/assemblage/"$sample"/final.contigs.fa | sed 's/>/> /' | sed 's/_/ /g' | sed 's/=/ /g' | awk '{if($1==">") print $1"Contig_"$3; else print $0}' > "$out_dir"/assemblage/"$sample"/assembly_raw.fa 

python "$master_dir"/parse_genome_500pb.py "$out_dir"/assemblage/"$sample"/assembly_raw.fa "$out_dir"/assemblage/"$sample".fa

mv "$out_dir"/assemblage/"$sample"/ "$out_dir"/assemblage/raw/

rm -r "$out_dir"/assemblage/raw/"$sample"/intermediate_contigs/

sbatch "$master_dir"/VirSorter.sh "$sample"

sbatch "$master_dir"/ResFinder.sh "$out_dir"/assemblage/"$sample".fa "$out_dir"/annotations/ResFinder/"$sample"/

sbatch "$master_dir"/Anvio.sh "$out_dir"/assemblage/"$sample".fa "$sample"
