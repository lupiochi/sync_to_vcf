import os
import pandas as pd
import numpy as np
import sys
import re

"""
This script will convert a SYNC file into a VCF file. The SYNC file should have the following format:
chr2	84	A	30:0:0:0:0:0	30:0:0:0:0:0	30:0:0:0:0:0
chr2	85	T	0:30:0:0:0:0	0:30:0:0:0:0	0:30:0:0:0:0
chr2	86	A	30:0:0:0:0:0	30:0:0:0:0:0	30:0:0:0:0:0
...

The script will generate a VCF file with the following format:
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	POP1	POP2	POP3	pop1_pop3_diff
chr2	977029	.	G	A	30	.	DP=30	AF	1.0	1.0	1.0	0.0
chr2	977031	.	A	A	30	.	DP=30	AF	0.033	0.0	0.0	0.033
chr2	977052	.	T	C	30	.	DP=30	AF	1.0	0.933	1.0	0.0

Remember to adjust the configuration variables below to suit your needs, including the names of the population columns, and if you want to calculate the absolute difference between allele frequencies of some populations.

"""

# Names of the population columns
pop_columns = ['pop1', 'pop2', 'pop3']

# Toggle here if you want to add a column with the absolute difference between the allele frequencies of pop1 and pop3.
add_diff_column = True
# Populations for calculate:
diff_pop1 = 'pop1'
diff_pop2 = 'pop3'

# Toggle here if you want to keep only rows in which the ALT allele is different from the REF allele.
keep_only_non_ref = True

# Toggle here if you want to add the genotype (GT) to the VCF file.
enable_gt = False

### FUNCTIONS

def get_af(row, pop, enable_gt=False):
    ref_value = row[f"{pop}_{row['ref']}"]
    alt_value = 30 - ref_value

    # Calculate AF
    AF = round(alt_value / 30, 3)

    # Format output
    output = f"{AF}"

    if enable_gt:
        # Calculate GT
        if ref_value == 30 and alt_value == 0:
            GT = "0/0"
        elif ref_value == 0 and alt_value == 30:
            GT = "1/1"
        else:
            GT = "0/1"
        
        output = f"{GT}:{AF}"

    return output

def sanitize_vcf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line.strip() + '\n')
            else:
                columns = [col.strip() for col in re.split(r'\t+', line.strip())]
                outfile.write('\t'.join(columns) + '\n')

### MAIN

# arguments
file_path = sys.argv[1]
output_folder_name = sys.argv[2]
regions_file_path = None
if len(sys.argv) == 4:
    regions_file_path = sys.argv[3]
else:
    print("You didn't provide any specific chromosome regions, proceeding to convert the entire SYNC file into VCF.")

# Create output folder if it doesn't exist
output_folder = output_folder_name
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Read input SYNC file
sync_input = pd.read_csv(file_path, sep="\t", header=None)
sync_input.insert(loc=3, column="alt", value='')
sync_input.columns = ["chromosome", "bp_number", "ref", "alt"] + pop_columns

# Filter by regions if provided
if regions_file_path:
    regions_df = pd.read_csv(regions_file_path)
    regions_of_interest = [(row['start'], row['end']) for _, row in regions_df.iterrows()]

    filtered_dfs = []
    for start, end in regions_of_interest:
        filtered_dfs.append(sync_input[(sync_input["bp_number"] >= start) & (sync_input["bp_number"] <= end)])
    sync_input = pd.concat(filtered_dfs)

# Split population columns
for pop in pop_columns:
    sync_input[[f'{pop}_A', f'{pop}_T', f'{pop}_C', f'{pop}_G', f'{pop}_N', f'{pop}_DEL']] = sync_input[pop].str.split(':', expand=True)
sync_input.drop(pop_columns + [f'{pop}_N' for pop in pop_columns] + [f'{pop}_DEL' for pop in pop_columns], axis=1, inplace=True)

# Convert columns to integers
columns_to_int = [f'{pop}_{allele}' for pop in pop_columns for allele in ['A', 'T', 'C', 'G']]
sync_input[columns_to_int] = sync_input[columns_to_int].astype('int32')

# Determine ALT allele
alleles = ['A', 'T', 'C', 'G']
def consensus_alt(row):
    counts = {allele: 0 for allele in alleles}
    for pop in pop_columns:
        for allele in alleles:
            counts[allele] += row[f'{pop}_{allele}']
    return max(counts, key=counts.get)

sync_input['alt'] = sync_input.apply(consensus_alt, axis=1)

# Keep only non-reference alleles if specified
if keep_only_non_ref:
    sync_input = sync_input[sync_input['ref'] != sync_input['alt']]

# Calculate allele frequencies
for pop in pop_columns:
    sync_input[f'{pop}_af'] = sync_input.apply(get_af, pop=pop, axis=1)

# Calculate the absolute difference between specified populations if specified
if add_diff_column:
    sync_input[f'{diff_pop1}_{diff_pop2}_diff'] = (sync_input[f'{diff_pop1}_af'].astype(float) - sync_input[f'{diff_pop2}_af'].astype(float)).abs()

# Save frequencies to CSV
output_csv = os.path.join(output_folder, f"{os.path.splitext(os.path.basename(file_path))[0]}_frequencies.csv")
sync_input.to_csv(output_csv, index=False)

# Prepare VCF format
sync_input.drop([f'{pop}_{allele}' for pop in pop_columns for allele in alleles], axis=1, inplace=True)
sync_input.rename(columns={"chromosome":"#CHROM", "bp_number":"POS", "ref":"REF", "alt":"ALT"}, inplace=True)
sync_input.rename(columns={f'{pop}_af': pop.upper() for pop in pop_columns}, inplace=True)
sync_input.insert(loc=2, column="ID", value='.')
sync_input.insert(loc=5, column="QUAL", value=30)
sync_input.insert(loc=6, column="FILTER", value='.')
sync_input.insert(loc=7, column="INFO", value='DP=30')
sync_input.insert(loc=8, column="FORMAT", value='AF')

# Generate a temporary VCF
temp_vcf_path = f"{os.path.splitext(os.path.basename(file_path))[0]}_temp.vcf"
sync_input.to_csv(temp_vcf_path, index=False, sep="\t")

# Sanitize the VCF file and generate the final output
final_vcf_path = os.path.join(output_folder, f"{os.path.splitext(os.path.basename(file_path))[0]}.vcf")
sanitize_vcf(temp_vcf_path, final_vcf_path)

# Delete the temporary VCF
os.remove(temp_vcf_path)
