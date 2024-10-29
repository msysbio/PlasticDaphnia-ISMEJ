#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 17:04:17 2023

@author: u0145079
"""

import os
from Bio import SearchIO
import pandas as pd

folder_path="/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/2_Bacterioplankton/BAM_counts/Blauwe_Poort"
CARD_path="/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/2_Bacterioplankton/Antibiotic_resistence/Blauwe_Poort/Blauwe_Poort-filtered_results.txt"
output_path="/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Publication/Review/Additional-analysis/Field-shotgun/Bacterioplankton/Antibiotic_resistence_genes"

# Initialize a dictionary for ORF counts
CARD_hits= set()
orf_counts = {}
orf_lengths = {}
total_reads = 0
   
qresults = SearchIO.parse(CARD_path, "blast-tab")

for hit in qresults:
    CARD_hits.add(hit.id)

# Process each feature count file
for file in os.listdir(folder_path):
    if file.endswith('_fc.out'):
        fc_file_path = os.path.join(folder_path, file)
        with open(fc_file_path, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if line.startswith("Geneid"):
                    sample_name = line.strip().split('\t')[-1]  # Assumes last column is the sample name
                    continue
                parts = line.strip().split('\t')
                orf_id = parts[0]
                count = int(parts[-1])  # Assumes last column is the count
                length = int(parts[5])
                if orf_id in CARD_hits:
                    orf_counts[orf_id] = orf_counts.get(orf_id, 0) + count
                    orf_lengths[orf_id] = length
                total_reads += count  # Accumulate the total reads

# Calculate adjCMR for each ORF
df = pd.DataFrame(list(orf_counts.items()), columns=['ORF', 'ReadCounts'])
df['Length'] = df['ORF'].map(orf_lengths)
df['adjCMR'] = df['ReadCounts'] / ((df['Length'] / 1000) * (total_reads / 1e6))


# Calculate the sum of adjCMR
sum_adjCMR = df['adjCMR'].sum()
mean_adjCMR = df['adjCMR'].mean()

# Save to a CSV file
base_name = os.path.splitext(os.path.basename(CARD_path))[0]
output = os.path.join(output_path, f"{base_name}_adjCMR_summed.csv")
df[['ORF', 'adjCMR']].to_csv(output, index=False)
print(f"Sum of adjCMR for {base_name}: {sum_adjCMR}")
print(f"Average of adjCMR for {base_name}: {mean_adjCMR}")

print("adjCMR calculations and CSV export completed for all .scores files.")
