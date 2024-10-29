#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 18:12:44 2024

@author: u0145079

Calculate abundances for each COG/KO mapping :
    COG - Absolute abundances (for Deseq2)
    KO - Relative normalised abundances (For KEGG patways)

"""


"""
COGs abundances:
"""
import os
import pandas as pd
from collections import defaultdict

# Function to load files into dictionaries
def load_file_to_dict(file_path, skip_rows=0, delimiter="\t"):
    data_dict = {}
    with open(file_path, 'r') as file:
        for _ in range(skip_rows):
            next(file)
        for line in file:
            parts = line.strip().split(delimiter)
            data_dict[parts[0]] = parts[1:]
    return data_dict

# Paths to input files
#Eggnog mapper output file:
eggnog_annotation_path = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/Eggnog/St.Donatus_cleaned.emapper.annotations'
#Path to FeatureCounts outputs:
feature_count_basepath = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/feature_count_output/St.Donatus'
#Path to idxstats outputs:
idx_basepath = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/BAM_counts/St.Donatus'
#Path to the desired output folder:
output_basepath = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/COGs'

# Load Eggnog annotations once as they are common for all samples
eggnog_data = load_file_to_dict(eggnog_annotation_path, skip_rows=5)

# Loop through files in the feature_count folder
for fc_filename in os.listdir(feature_count_basepath):
    if fc_filename.endswith(".sorted_fc.out"):
        sample_name = fc_filename.split('.')[0]  # Extract the sample name (assuming the format is 'GC123456.sorted_fc.out')
        
        # Match the sample name with the corresponding idx file
        idx_filename = sample_name + '_contig_counts.txt'
        if idx_filename in os.listdir(idx_basepath):
            
            # Load data for the specific sample
            feature_count_data = load_file_to_dict(os.path.join(feature_count_basepath, fc_filename), skip_rows=1)
            idx_data = load_file_to_dict(os.path.join(idx_basepath, idx_filename))
            
            # Calculate total reads for the sample
            total_reads_per_sample = sum([int(value[1]) + int(value[2]) for value in idx_data.values()])
            
            # Compute results + normalize the abundances
            results = {}
            for key, value in eggnog_data.items():
                if key in feature_count_data:
                    length = float(feature_count_data[key][4])  # Assuming the 5th position has the length info
                    count = float(feature_count_data[key][5])  # Assuming the 6th position has the read count
                    
                    normalized_count = count / length
                    adjCMR = count / ((length / 1000) * (total_reads_per_sample / 1000000))
                    
                    COG = value[3].split('@')[0]  # Assuming the 4th position has the COG info
                    if COG not in results:
                        results[COG] = {'count': 0}
                    results[COG]['count'] += count
            
            # Write results to output file
            output_file = os.path.join(output_basepath, f"FINAL_{sample_name}_out.tsv")
            with open(output_file, 'w') as outfile:
                outfile.write("COG\tCount\n")
                for key, value in results.items():
                    outfile.write(f"{key}\t{value['count']}\n")


"""
KEGGs abundances:
"""

# Paths
#KEGG output basepath:
output_basepath = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/KEGG-KOs/St.Donatus'

# Loop through files in the feature_count folder
for fc_filename in os.listdir(feature_count_basepath):
    if fc_filename.endswith(".sorted_fc.out"):
        sample_name = fc_filename.split('.')[0]  # Extract the sample name (assuming the format is 'GC123456.sorted_fc.out')
        
        # Match the sample name with the corresponding idx file
        idx_filename = sample_name + '_contig_counts.txt'
        if idx_filename in os.listdir(idx_basepath):
            
            # Load data for the specific sample
            feature_count_data = load_file_to_dict(os.path.join(feature_count_basepath, fc_filename), skip_rows=1)
            idx_data = load_file_to_dict(os.path.join(idx_basepath, idx_filename))
            
            # Calculate total reads for the sample
            total_reads_per_sample = sum([int(value[1]) + int(value[2]) for value in idx_data.values()])
            
            # Compute results
            results = {}
            for key, value in eggnog_data.items():
                if key in feature_count_data:
                    length = float(feature_count_data[key][4])  # Assuming the 5th position has the length info
                    count = float(feature_count_data[key][5])  # Assuming the 6th position has the read count
                    adjCMR = count / ((length / 1000) * (total_reads_per_sample / 1000000))
                    
                    KO = value[10]  # Assuming the 10th position has the KEGG info
                    if KO not in results:
                        results[KO] = {'CMR': 0}
                    results[KO]['CMR'] += adjCMR
          
            # Write results to output file
            output_file = os.path.join(output_basepath, f"FINAL_{sample_name}_out.tsv")
            with open(output_file, 'w') as outfile:
                outfile.write("KEGG\tadjCMR\n")
                for key, value in results.items():
                    outfile.write(f"{key}\t{value['CMR']}\n")






