# MetaKids cohort analysis

this set of scripts allows to analyze metagenomic data from the MetaKids cohort.

## Before starting

### set the different arguments and mandatory directories

```sh
out_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/
master_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/Master_scripts/
sample=XX
```


```sh
mkdir -p "$out_dir"/temp/
mkdir -p "$out_dir"/FastQ/HiC_digested/
mkdir -p "$out_dir"/FastQ/SG_cleaned/
mkdir -p "$out_dir"/FastQ/temp/
mkdir -p "$out_dir"/assemblage/
mkdir -p "$out_dir"/annotations/
mkdir -p "$out_dir"/MetaTOR/metator/"$sample"
mkdir -p "$out_dir"/MetaTOR/metamge/"$sample"
mkdir -p "$out_dir"/anvio/
mkdir -p "$out_dir"/coverage/"$sample"/
```


##########################################################################################################################

## Raw reads treatment

### copy the raw from GAIA storage unit and rename them appropriately

NB: do not forget do that by logging on the sftpcampus server

```sh
bash "$master_dir"/reads_treatment.sh "$out_dir"/FastQ/"$sample"/ "$sample"
```
##########################################################################################################################

## Shotgun reads treatment, assembly and annotation

### clean SG reads using Trimmomatic and start the Assembly and its annotation

```sh
sbatch "$master_dir"/Trimmomatic_SG.sh "$sample"
```

NB: this script will automatically launch the assembly once the reads have been cleaned and then the annotation (Anvio contig db, resfinder, )

### MGEs annotation

launch genomad pipeline (do not forget to activate the conda environment)

```sh
conda activate genomad 
```

```sh
sbatch "$master_dir"/Genomad.sh "$out_dir"/assemblage/"$sample".fa "$out_dir"/annotations/genomad/"$sample"/
```

```sh
conda deactivate
```


### AMR annotation

```sh
sbatch "$master_dir"/ResFinder.sh "$out_dir"/assemblage/"$sample".fa "$out_dir"/annotations/ResFinder/"$sample"/
```

### anvio DB, 

```sh
sbatch "$master_dir"/Anvio_init_MK.sh /pasteur/zeus/projets/p02/rsg_fast/mmarbout/database/ "$out_dir"/assemblage/"$sample".fa "$sample" "$out_dir"/FastQ/ "$out_dir"
```

##########################################################################################################################

## HiC reads process and MetaTOR pipeline for the first binning

### clean HiC reads using Trimmomatic and digest them

NPO - conda activate metator

```sh
conda activate metator
sample=XX
```


```sh
for i in $(ls /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/FastQ/"$sample"/ | grep "MK" | grep "Arima" | sed 's/_R/ /' | awk '{print $1}' | sort -u)
do

	sbatch "$master_dir"/Trimmomatic_HiC.sh "$out_dir"/FastQ/temp/ "$out_dir"/FastQ/"$sample"/"$i"_R1.fq.gz "$out_dir"/FastQ/"$sample"/"$i"_R2.fq.gz "$i" 

done
```

### predigest the HiC reads


```sh
for i in $(ls /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/FastQ/temp/ | grep "$sample" | grep "clean" | grep "Arima" | sed 's/_R/ /' | awk '{print $1}' | sort -u)
do

	sbatch "$master_dir"/reads_digestion.sh "$out_dir"/FastQ/temp/"$i"_R1.fq.gz "$out_dir"/FastQ/temp/"$i"_R2.fq.gz "$out_dir"/FastQ/HiC_digested/"$i" DpnII,HinfI 

done
```





## coverage files


### RPKM data files

count the number of hits per contigs

### reads mapping back and BAM files generation



use the script in perl to process BAM files

```sh
"$master_dir"/jgi_summarize_bam_contig_depths --outputDepth "$out_dir"/coverage/"$sample"/coverage_"$sample".txt --showDepth --minContigLength 500 "$out_dir"/BAM/*.bam
```

move outputfiles and process the data

```sh
mv "$out_dir"/BAM/"$sample"*.depth "$out_dir"/coverage/"$sample"/
```



## Binning of MAGs and MGEs using metator

### MetaTOR

launch metator pipeline (do not forget to activate the conda environment)

```sh
conda activate metator 
```

```sh
sbatch "$master_dir"/MetaTOR.sh "$out_dir"/MetaTOR/metator/"$sample"/ "$out_dir"/assemblage/"$sample".fa $(ls "$out_dir"/FastQ/HiC_digested/ | sed 's/_clean/ /' | grep "$sample" | awk '{print $1}' | sort -u | awk '{print "'$out_dir'""/FastQ/HiC_digested/"$1"_cleaned_R1.fq.gz"}' | paste -s | sed 's/\t/,/g') $(ls "$out_dir"/FastQ/HiC_digested/ | sed 's/_clean/ /' | grep "$sample" | awk '{print $1}' | sort -u | awk '{print "'$out_dir'""/FastQ/HiC_digested/"$1"_cleaned_R2.fq.gz"}' | paste -s | sed 's/\t/,/g')
```

### MetaMGE

first, select the contigs

```sh
cat "$out_dir"/annotations/genomad/"$sample"/"$sample"_summary/"$sample"_plasmid_summary.tsv | sed '1d' | awk '{print $1}' | grep -v "provirus" > "$out_dir"/MetaTOR/metamge/"$sample"/contig_mge.txt
cat "$out_dir"/annotations/genomad/"$sample"/"$sample"_summary/"$sample"_virus_summary.tsv | sed '1d' | awk '{print $1}' | grep -v "provirus" >> "$out_dir"/MetaTOR/metamge/"$sample"/contig_mge.txt
```

launch the pipeline MetaMGE

```sh
sbatch "$master_dir"/MetaMGE.sh \
	"$out_dir"/MetaTOR/metamge/"$sample"/ \
	"$out_dir"/MK/MetaTOR/metator/"$sample"/ \
	"$out_dir"/MetaTOR/metamge/"$sample"/contig_mge.txt \
	"$out_dir"/assemblage/"$sample".fa \
	"$out_dir"/MetaTOR/metator/"$sample"/ \
	alignment 0.8
```

### recovering the Fasta files cleaned for MAGs and MGEs

scrits matrice_generation

### annotation of the retrieved bins (GTDB-tk, CheckM, GeNomad, Pharokka, Iphop and CheckV)

```sh
sbatch "$master_dir"/checkM.sh "$out_dir"/MetaTOR/metator/"$sample"/final_bin/ "$out_dir"/checkM/"$sample"/ "$sample"_MAG
```

### construction of a large contact map encompassing the MGEs and the MAGs




## Contact

### Authors

* martial.marbouty@pasteur.fr

### Research lab

[Spatial Regulation of Genomes](https://research.pasteur.fr/en/team/spatial-regulation-of-genomes/) (Institut Pasteur, Paris)

