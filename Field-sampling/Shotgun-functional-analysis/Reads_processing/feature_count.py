#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 12:19:42 2024

@author: u0145079
"""

from Bio import SeqIO
import os
import subprocess

orf='/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/ORFs/St.Donatus_proteins.fa'
output='/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/feature_count_output/St.Donatus'
####SAM_BAM FILES:
mapping_files='/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/BAM_counts/St.Donatus'


##Extract .ORFs and convert it to .SAF
#ORF file to .SAF
with open(orf, "r") as ffn, open(os.path.join(output, 'saf_output.saf'), "w") as saf:
    # Write the SAF header
    saf.write("\t".join(["GeneID", "Chr", "Start", "End", "Strand"]) + "\n")
    
    # Parse each FASTA entry in the input file
    for record in SeqIO.parse(ffn, "fasta"):
        # Extract information from the record description
        fields = record.description.split("#")
        gene_id = fields[0].strip()
        start = fields[1].strip()
        end = fields[2].strip()
        strand = fields[3].strip()
        if strand == '1':
            strand = '+'
        elif strand == '-1':
            strand = '-'
        
        # Remove last suffix from Chr
        chr = "_".join(gene_id.split("_")[:-1])
        
        # Write the SAF entry
        saf.write("\t".join([gene_id, chr, start, end, strand]) + "\n")

#Check if featurecounts in the path:
found = False
for dir in os.environ["PATH"].split(os.pathsep):
    fpath = os.path.join(dir, "featureCounts")
    if os.path.exists(fpath) and os.access(fpath, os.X_OK):
        print("featureCounts is in the path")
        found = True
        break

if not found:
    print("featureCounts not in the path")     
        
for bam_file in os.listdir(mapping_files):
    if bam_file.endswith(".bam"):
        full_sam_path = os.path.join(mapping_files, bam_file)
        base_name = os.path.splitext(bam_file)[0]
        fc_output = os.path.join(output, f"{base_name}_fc.out")
        fc_input = os.path.join(output, 'saf_output.saf')
        fc_command = f"featureCounts -a {fc_input} -F SAF -o {fc_output} {full_sam_path}"

        # Run featureCounts
        process = subprocess.Popen(fc_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        # Check for errors
        if process.returncode != 0:
            print(f"featureCounts failed with error code {process.returncode}")
            print(f"Stderr: {stderr.decode()}")
        else:
            print(f"featureCounts succeeded for {bam_file}")
            print(f"Stdout: {stdout.decode()}")

