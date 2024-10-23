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

#out_dir
out_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/

#Master script diretory
master_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/Master_scripts/

module load R/3.6.2
module load mcl/14-137 muscle/3.8.31 diamond/2.0.4 blast+/2.10.0 MetaGeneAnnotator/2008-8-19
module load SQLite/3.31.0  blast+/2.10.0
module load centrifuge/1.0.4-beta
module load infernal
module load R/3.6.2 prodigal mcl/14-137 muscle/3.8.31 hmmer/3.2.1 diamond/2.0.4 samtools/1.10 FastTree/2.1.11 centrifuge/1.0.4-beta SQLite/3.31.0 tRNAscan-SE/2.0.7 blast+/2.10.0 megahit/1.2.9 SPAdes/3.15.0 bowtie2/2.3.5.1 bwa/0.7.17 trimal/1.4.1 IQ-TREE/2.0.6
module load prodigal-gv/
module load prodigal FastTree/2.1.11
module load anvio/
module add MMseqs2/10-6d92c
module add bedtools/2.29.2 
module load gcc
module load perl/5.30.1

mkdir -p "$out_dir"/Anvio/

####### Anvio contig database ##############

anvi-gen-contigs-database -L 0 -T 96 --project-name "$sample" -f "$assembly" -o "$out_dir"/Anvio/"$sample".db

anvi-run-hmms -T 96 -c "$out_dir"/Anvio/"$sample".db

anvi-run-ncbi-cogs -T 96 --cog-version COG20 --cog-data-dir "$db"/COG_2020 -c "$out_dir"/Anvio/"$sample".db

anvi-run-pfams -T 96 --pfam-data-dir "$db"/Pfam_v32 -c "$out_dir"/Anvio/"$sample".db


####### Anvio BAM profile ##############

bowtie2-build "$assembly" "$out_dir"/index/"$sample"

for lib in $(ls "$out_dir"/FastQ/SG_cleaned/ | grep "$sample" | grep "Shot" | sed 's/_R/ /' | awk '{print $1}' | sort -u)
do
    
    bowtie2 --very-sensitive-local -p 96  -x "$out_dir"/index/"$sample" -1 "$out_dir"/FastQ/SG_cleaned/"$lib"_R1.fq.gz -2 "$out_dir"/FastQ/SG_cleaned/"$lib"_R2.fq.gz -S "$out_dir"/BAM/align.sam
    samtools view --threads 96 -F 0x904 -bS "$out_dir"/BAM/align.sam > "$out_dir"/BAM/align_raw.bam
    anvi-init-bam -T 96 "$out_dir"/BAM/align_raw.bam -o "$out_dir"/BAM/"$lib".bam
    anvi-profile -T 96 -i "$out_dir"/BAM/"$lib".bam -c "$out_dir"/Anvio/"$sample".db --sample-name "$lib" --min-contig-length 2000
    rm  "$out_dir"/BAM/align.sam
    rm "$out_dir"/BAM/align_raw.bam



done


anvi-merge $(ls "$out_dir"/BAM/ | grep "$sample" | grep "PROFILE" | awk -F "." '{print $1}' | sed 's/_T/ /' | sort -k 2,2 -g | \
     awk '{print "'$out_dir'""/BAM/"$1"_T"$2".bam-ANVIO_PROFILE/PROFILE.db"}' | paste -s | sed 's/\t/ /g') \
    -o "$out_dir"/Anvio/"$sample"_cinetiq \
    -c "$out_dir"/Anvio/"$sample".db


####### VirSorter annotation transfert to Anvio ############## a VOIR

anvi-export-table "$out_dir"/Anvio/"$sample".db --table "$out_dir"/Anvio/table/splits_basic_info_"$sample"

anvi-export-gene-calls --gene-caller prodigal -c "$out_dir"/Anvio/"$sample".db -o "$out_dir"/Anvio/table/all_gene_calls_"$sample".txt

"$master_dir"/./virsorter_to_anvio.py \
    --affi-file "$out_dir"/annotations/VirSorter/"$sample"/VIRSorter_affi-contigs.tab \
    --global-file "$out_dir"/annotations/VirSorter/"$sample"/VIRSorter_global-phage-signal.csv \
    --splits-info "$out_dir"/Anvio/table/splits_basic_info_"$sample".txt \
    --db 2 \
    --anvio-gene-calls "$out_dir"/Anvio/table/all_gene_calls_"$sample".txt \
    --hallmark-functions "$master_dir"/db2_hallmark_functions.txt \
    --addl-info "$out_dir"/Anvio/table/virsorter_additional_info_"$sample".txt \
    --output-functions "$out_dir"/Anvio/table/virsorter_functions_"$sample".txt \
    --exclude-cat3 \
    --exclude-prophages



anvi-import-misc-data virsorter_additional_info.txt \
                      -c "$out_dir"/Anvio/"$sample".db \
                      --target-data-table items

anvi-import-functions 
                       -c CONTIGS.db \
                       -i virsorter_annotations.txt






