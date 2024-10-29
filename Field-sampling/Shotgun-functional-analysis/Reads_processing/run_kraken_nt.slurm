#!/bin/bash -l

#SBATCH --cluster=genius
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --account=lmsysbio
#SBATCH --partition=superdome

cd /staging/leuven/stg_00106/PlasticDaphnia-Field-Shotgun/reads/reads/quality-filtered

module load Perl/5.28.1-GCCcore-6.4.0

export PATH="/staging/leuven/stg_00106/tools/kraken2/kraken2-2.1.3:$PATH"

# Directory where Kraken2 reports will be saved
REPORT_DIR="/staging/leuven/stg_00106/PlasticDaphnia-Field-Shotgun/Daphnia_microbiome/KRAKEN/nt_kraken"

# Kraken2 database path
DB_PATH="/staging/leuven/stg_00106/databases/kraken-nt/nt_database"

while IFS= read -r id; do
	kraken2 --db ${DB_PATH} --threads 14 --report "${REPORT_DIR}/${id}.k2report" --quick --output "${REPORT_DIR}/${id}.kraken2.output" --paired "${id}_qc.R1.fastq.gz" "${id}_qc.R2.fastq.gz" 
done < ./Microbiome-druga-tura.txt








