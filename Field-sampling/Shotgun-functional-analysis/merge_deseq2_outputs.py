#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:04:33 2023

@author: u0145079
"""

import os
import pandas as pd
from collections import defaultdict

# Change to the directory containing the files, if needed
os.chdir('/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/2_Bacterioplankton/Eggnog/all')

# Initialize an empty dictionary to hold the data
merged_data = defaultdict(dict)

# List the files in the directory
files = [f for f in os.listdir() if f.startswith('FINAL_') and f.endswith('_out.tsv')]

for file in files:
    # Extract the sample name from the file name
    sample_name = file.split('_')[1]

    # Read the file into a Pandas DataFrame
    df = pd.read_csv(file, sep='\t')

    # Drop unnecessary columns
    df = df[['COG', 'Count']]

    # Rename the 'Count' column to the sample name
    df.rename(columns={'Count': sample_name}, inplace=True)

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
final_df.index.name = 'COG'

# Save the final DataFrame to a new TSV file
final_df.to_csv('all_data.tsv', sep='\t')

# Create and save a DataFrame containing only the selected samples
#selected_samples = ['GC125618','GC127864','GC127853','GC127871','GC127874'] 
#if selected_samples:
    #selected_df = final_df[selected_samples]
    #selected_df.to_csv('water_samples_1.tsv', sep='\t')
