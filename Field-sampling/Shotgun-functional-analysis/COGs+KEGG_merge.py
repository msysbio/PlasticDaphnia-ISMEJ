#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 18:19:39 2024

@author: u0145079

Merge separate COGs/KOs annotation per sample into one file 

"""

import os
import pandas as pd
from collections import defaultdict


# Base directory containing the subdirectories with 
base_dir = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/3_Daphnia_microbiome/KEGG-KOs'

# Initialize an empty dictionary to hold the data
merged_data = defaultdict(dict)

# Walk through each directory and subdirectory in the base directory
for root, dirs, files in os.walk(base_dir):
    for file in files:
        # Check if the file name matches the pattern
        if file.startswith('FINAL_') and file.endswith('_out.tsv'):
            # Construct the full file path
            file_path = os.path.join(root, file)
            
            # Extract the sample name from the file name
            sample_name = file.split('_')[1]

            # Read the file into a Pandas DataFrame
            df = pd.read_csv(file_path, sep='\t')

            # Drop unnecessary columns, keeping only 'COG' and 'Count'
            df = df[['KEGG', 'adjCMR']] #For KEGG
            #df = df[['COG', 'Count']] #For COGs

            # Rename the 'Count' column to the sample name
            df.rename(columns={'adjCMR': sample_name}, inplace=True) #For KEGG
            #df.rename(columns={'Count': sample_name}, inplace=True) #For COG

            # Merge this DataFrame into the merged_data dictionary
            for index, row in df.iterrows():
                cog = row['KEGG']
                #cog = row['COG']
                count = row[sample_name]
                merged_data[cog][sample_name] = count
            

# Create a DataFrame from the merged_data dictionary
final_df = pd.DataFrame.from_dict(merged_data, orient='index')

# Replace NaN with 0
final_df.fillna(0, inplace=True)

# Reset index name to 'COG'
final_df.index.name = 'KEGG'
#final_df.index.name = 'COG'

# Specify the output directory and file name for the final DataFrame
output_dir = base_dir  # Change this if you want to save the output elsewhere
output_file = os.path.join(output_dir, 'KEGG_all_samples_raw_counts.tsv')
#output_file = os.path.join(output_dir, 'COGs_all_samples_raw_counts.tsv')
# Save the final DataFrame to a new TSV file
#final_df.to_csv(output_file, sep='\t')

"""
Selected samples (Natural pond vs Artifical pond) - adjCMR counts for KEGG mapping
"""

#NP "GC125610", "GC125611", "GC125612", "GC127825", "GC127826", "GC127827", "GC127828", "GC127829", "GC127830", "GC127831", "GC127832", "GC127833", "GC127834", "GC127835", "GC127836", "GC127837", "GC127838", "GC127839", "GC127840", "GC127841", "GC127843", "GC127844", "GC127845", "GC127852", "GC127853", "GC127856", "GC127866", "GC127868"
#AP "GC125618","GC127846","GC127847","GC127848","GC127849","GC127855","GC127864","GC127867","GC127871","GC127874"
selected_samples = ["GC125610", "GC125611", "GC125612", "GC127825", "GC127826", "GC127827", "GC127828", "GC127829", "GC127830", "GC127831", "GC127832", "GC127833", "GC127834", "GC127835", "GC127836", "GC127837", "GC127838", "GC127839", "GC127840", "GC127841", "GC127843", "GC127844", "GC127845", "GC127852", "GC127853", "GC127856", "GC127866", "GC127868"]
selected_df = final_df[selected_samples]

# Add a new column for the sum of adjCMR values across the selected samples
selected_df['Average_adjCMR'] = selected_df.mean(axis=1)

# Remove rows where the sum of adjCMR values is equal to 0
selected_df = selected_df[selected_df['Average_adjCMR'] > 0]

# Save the filtered DataFrame to a new CSV file
output_file = os.path.join(base_dir, 'filtered_COGs_adjCMR_NP.tsv')
selected_df.to_csv(output_file, sep='\t')
