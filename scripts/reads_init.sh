#!/bin/sh

################# input ###################

#path to the folder where the reads will be stored
out_fold=$1

#id of the kid
sample=$2


out_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/
master_dir=/pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/scripts/

mkdir -p "$out_fold"/tmp/

################ code #####################

for lib in $(cat /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/data_lib/Banques_MK.csv | awk -F ";" '$1=="'$sample'" {print $4}')
	do
	name=$(cat /pasteur/zeus/projets/p02/rsg_fast/mmarbout/projets/MK/data_lib/Banques_MK.csv | awk -F ";" '$1=="'$sample'" && $4=="'$lib'" {print $6}')
	cp /pasteur/gaia/projets/p02/Rsg_reads/2_METAG/MetaKids/"$lib"_nvq_R1.fq.gz "$out_fold"/"$name"_R1.fq.gz
	cp /pasteur/gaia/projets/p02/Rsg_reads/2_METAG/MetaKids/"$lib"_nvq_R2.fq.gz "$out_fold"/"$name"_R2.fq.gz
done
