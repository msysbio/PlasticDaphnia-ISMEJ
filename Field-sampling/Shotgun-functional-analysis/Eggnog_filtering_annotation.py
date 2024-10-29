#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 17:00:39 2024

@author: u0145079

Extract eggnog annotation for cleaned/reduced fasta set
Fasta after filtering - >100aa limit + non-redundant

"""
import pandas as pd

# File paths
eggnog_annotation = "/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/Eggnog/St.Donatus.emapper.annotations"
fasta_file = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/ORFs/Unique/St.Donatus_proteins_unique_clean.fa'
output_file = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/Eggnog/St.Donatus_cleaned.emapper.annotations'

# Load the eggnog annotation file, skipping header lines that are not needed
annotation_df = pd.read_csv(eggnog_annotation, sep='\t', skiprows=4)

# Read the FASTA file and collect all sequence IDs in a set for quick lookup
fasta_ids = set()
with open(fasta_file, 'r') as fasta:
    for line in fasta:
        if line.startswith('>'):
            # Extract the ID part from the FASTA header
            fasta_id = line.split()[0][1:]  # Remove '>' and split by space, take the first part
            fasta_ids.add(fasta_id)

# Filter the annotations to include only those that match the reduced FASTA list
# Set the '#query' column as the index of the DataFrame
reduced_annotations = annotation_df[annotation_df['#query'].isin(fasta_ids)]

# Save the filtered DataFrame to a new CSV file with the index (which is now your ORF identifiers)
reduced_annotations.to_csv(output_file, sep='\t', index=False)