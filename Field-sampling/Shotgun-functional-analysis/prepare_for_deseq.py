from Bio import SeqIO
import os
import subprocess

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

# Paths
eggnog_annotation_path = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/2_Bacterioplankton/Eggnog/WWTP_Olsene.emapper.annotations'
feature_count_basepath = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/2_Bacterioplankton/BAM_counts/WWTP_Olsene'
idx_basepath = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/2_Bacterioplankton/BAM_counts/WWTP_Olsene'
output_basepath = '/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Field/Shotgun/2_Bacterioplankton/Eggnog/WWTP_Olsene'

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
            total_reads_per_sample = sum([int(value[1]) for value in idx_data.values()])
            
            # Compute results
            results = {}
            for key, value in eggnog_data.items():
                if key in feature_count_data:
                    length = float(feature_count_data[key][4])  # Assuming the 5th position has the length info
                    count = float(feature_count_data[key][5])  # Assuming the 6th position has the read count
                    
                    normalized_count = count / length
                    adjCMR = count / ((length / 1000) * (total_reads_per_sample / 1000000))
                    
                    COG = value[3].split('@')[0]  # Assuming the 4th position has the COG info
                    if COG not in results:
                        results[COG] = {'CMR': 0, 'normalized': 0, 'count': 0}
                    results[COG]['CMR'] += adjCMR
                    results[COG]['normalized'] += normalized_count
                    results[COG]['count'] += count
            
            # Write results to output file
            output_file = os.path.join(output_basepath, f"FINAL_{sample_name}_out.tsv")
            with open(output_file, 'w') as outfile:
                outfile.write("COG\tadjCMR\tNormalized\tCount\n")
                for key, value in results.items():
                    outfile.write(f"{key}\t{value['CMR']}\t{value['normalized']}\t{value['count']}\n")
