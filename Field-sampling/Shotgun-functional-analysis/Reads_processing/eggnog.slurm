#!/bin/bash -l

#SBATCH --cluster=genius
#SBATCH --time=24:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --account=lmsysbio

module load Python/3.7.2-foss-2018a

cd /staging/leuven/stg_00106/PlasticDaphnia-Field-Shotgun/Daphnia_microbiome/ORFs

while IFS= read -r id; do	
	python3 /staging/leuven/stg_00106/tools/eggnogg/eggnog-mapper-2.1.6/emapper.py -i /staging/leuven/stg_00106/PlasticDaphnia-Field-Shotgun/Daphnia_microbiome/ORFs/${id}_proteins.fa -o ${id} --output_dir /staging/leuven/stg_00106/PlasticDaphnia-Field-Shotgun/Daphnia_microbiome/Eggnog --cpu 36
done < ./ID


