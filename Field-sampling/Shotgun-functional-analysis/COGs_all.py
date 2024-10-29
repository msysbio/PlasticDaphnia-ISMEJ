#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:04:55 2024

@author: u0145079

Merge all COGs: Bacterioplankton + Host
"""

import os
import pandas as pd
from collections import defaultdict


# Base directory containing the subdirectories with 
base_dir = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/14_All_together_datasets/COGs'

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
            df = df[['COG', 'Count']] #For COGs

            # Rename the 'Count' column to the sample name
            df.rename(columns={'Count': sample_name}, inplace=True) #For COG

            # Merge this DataFrame into the merged_data dictionary
            for index, row in df.iterrows():
                cog = row['COG']
                count = row[sample_name]
                merged_data[cog][sample_name] = count
            

# Create a DataFrame from the merged_data dictionary
final_df = pd.DataFrame.from_dict(merged_data, orient='index')

# Replace NaN with 0
final_df.fillna(0, inplace=True)

# Reset index name to 'COG'
#final_df.index.name = 'KEGG'
final_df.index.name = 'COG'

# Specify the output directory and file name for the final DataFrame
output_dir = base_dir  # Change this if you want to save the output elsewhere
output_file = os.path.join(output_dir, 'COGs_all_samples_raw_counts.tsv')
final_df.to_csv(output_file, sep='\t')
