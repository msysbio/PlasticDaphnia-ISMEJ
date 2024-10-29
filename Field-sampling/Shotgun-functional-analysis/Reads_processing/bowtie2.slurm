#!/bin/bash

#SBATCH --cluster=genius
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --account=lkulak_daphnia_gp

cd /staging/leuven/stg_00106/PlasticDaphnia-Field-Shotgun/reads/reads/quality-filtered

#run bowtie2
while IFS= read -r id; do
	/staging/leuven/stg_00106/tools/bowtie2-2.5.2-linux-x86_64/bowtie2 -p 36 -x /staging/leuven/stg_00106/PlasticDaphnia-Field-Shotgun/Daphnia_microbiome/contamination_index/index/all_contam -1 "${id}_qc.R1.fastq.gz" -2 "${id}_qc.R2.fastq.gz" -S /staging/leuven/stg_00106/PlasticDaphnia-Field-Shotgun/Daphnia_microbiome/contamination_filtered/"${id}_mapped_and_unmapped.sam"
done < ./Microbiome.txt

#make bowtie2 index
/staging/leuven/stg_00106/tools/bowtie2-2.5.2-linux-x86_64/bowtie2-build --threads 36 /staging/leuven/stg_00106/PlasticDaphnia-Field-Shotgun/Daphnia_microbiome/contamination_index/all_contam.fna index/all_contam
