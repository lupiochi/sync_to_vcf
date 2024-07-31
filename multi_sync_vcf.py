import os
import sys
import pandas as pd
import subprocess
import re

"""
This script processes several SYNC files to generate VCF and CSV files and concatenates the results.
It is essentially a master script to automate the conversion of multiple SYNC files to VCF format usin the sync_to_vcf.py script.

It includes the following steps:

1. Input Parameters: The script takes three arguments: the SYNC folder, the output folder, and an optional regions folder. The output folder will contain the generated files.

2. VCF Sanitization Function: A function sanitize_vcf is defined to clean up the VCF files by ensuring proper tab-separated formatting and removing any extra whitespace.

3. Processing SYNC Files: The main function process_files_in_folder performs the following tasks:
    a) Checks if the output folder exists and creates it if necessary.
    b) Iterates over all SYNC files in the specified folder.
    c) For each SYNC file, constructs the corresponding chromosome name and forms a command to convert the SYNC file to VCF using an external script (sync_to_vcf.py).
        If a regions folder is provided, it includes the appropriate regions file in the command.
    d) Executes the command and handles any errors during the conversion.
    e) Loads the resulting CSV file for the SYNC file and appends it to a list of dataframes.
    f) Reads the generated VCF file, extracts and stores the header if it's the first VCF file, and appends the data to a temporary concatenated VCF file.

4. Concatenation and Finalization:
    a) After processing all SYNC files, the temporary concatenated VCF file is sanitized and saved as the final concatenated VCF file.
    b) All loaded CSV dataframes are concatenated into a single dataframe and saved as a concatenated CSV file in the output folder.

5. Script Execution: The script checks for the correct number of command-line arguments and initiates the processing function with the provided parameters.

Example Usage:
python multi_sync_vcf.py <sync_folder> <output_folder_name> [regions_folder_name]

The script will output a concatenated VCF file and a concatenated CSV file in the specified output folder.

"""

# Function to sanitize VCF
def sanitize_vcf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line.strip() + '\n')
            else:
                columns = [col.strip() for col in re.split(r'\t+', line.strip())]
                outfile.write('\t'.join(columns) + '\n')

def process_files_in_folder(sync_folder, output_folder, regions_folder=None):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created directory: {output_folder}")
    
    # List to store dataframes for concatenation
    all_dataframes = []
    vcf_header = None
    concatenated_vcf_path = os.path.join(output_folder, "concatenated_output.vcf")
    temp_concatenated_vcf_path = os.path.join(output_folder, "temp_concatenated_output.vcf")
    
    # Ensure a temporary concatenated VCF file is created for sanitization
    if os.path.exists(temp_concatenated_vcf_path):
        os.remove(temp_concatenated_vcf_path)

    # Iterate through all sync files in the sync folder
    for sync_filename in os.listdir(sync_folder):
        if sync_filename.endswith(".sync"):
            sync_file_path = os.path.join(sync_folder, sync_filename)
            chromosome = sync_filename.replace(".sync", "")  # Extract chromosome from filename
            
            command = ["python", "sync_to_vcf.py", sync_file_path, output_folder]
            if regions_folder:
                regions_filename = f"{chromosome}_regions.csv"
                regions_file = os.path.join(regions_folder, regions_filename)
                if os.path.exists(regions_file):
                    command.append(regions_file)
                else:
                    print(f"Regions file not found for {chromosome}, proceeding to convert the entire SYNC file into VCF. This might take a while.")
            else:
                print("You didn't provide a regions folder, proceeding to convert the entire SYNC files into VCF. This might take a while.")

            result = subprocess.run(command, capture_output=True, text=True)

            if result.returncode != 0:
                print(f"Error processing {sync_file_path}")
                continue
            else:
                print(f"Processed {sync_file_path}")
            
            # Load the generated CSV file
            output_csv = os.path.join(output_folder, f"{chromosome}_frequencies.csv")
            if os.path.exists(output_csv):
                out_csv = pd.read_csv(output_csv)
                all_dataframes.append(out_csv)
                print(f"Loaded CSV for {sync_file_path}")
            else:
                print(f"CSV file not found for {sync_file_path}, expected at {output_csv}")
                continue
            
            # Add to concatenated VCF 
            temp_vcf = os.path.join(output_folder, f"{chromosome}.vcf")
            if os.path.exists(temp_vcf):
                with open(temp_vcf, 'r') as infile:
                    lines = infile.readlines()
                    if vcf_header is None:
                        vcf_header = [line for line in lines if line.startswith('#')]
                        with open(temp_concatenated_vcf_path, 'w') as outfile:
                            outfile.writelines(vcf_header)
                    with open(temp_concatenated_vcf_path, 'a') as outfile:
                        for line in lines:
                            if not line.startswith('#'):
                                outfile.write(line)
    
    # Sanitize the concatenated VCF
    if os.path.exists(temp_concatenated_vcf_path):
        sanitize_vcf(temp_concatenated_vcf_path, concatenated_vcf_path)
        os.remove(temp_concatenated_vcf_path)
    
    # Concatenate all dataframes
    if all_dataframes:
        concatenated_df = pd.concat(all_dataframes)
        concatenated_output_path = os.path.join(output_folder, "concatenated_output.csv")
        concatenated_df.to_csv(concatenated_output_path, index=False)
        print("All files processed and concatenated successfully.")

# Main function
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python multi_sync_vcf.py <sync_folder> <output_folder_name> [regions_folder_name]")
        sys.exit(1)

    sync_folder = sys.argv[1]
    output_folder_name = sys.argv[2]
    regions_folder_name = sys.argv[3] if len(sys.argv) == 4 else None

    script_folder = os.path.dirname(os.path.abspath(__file__))
    output_folder = os.path.join(script_folder, output_folder_name)  # output folder for the CSV and VCF files + the concatenated files

    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    print(f"Processing files in folder: {sync_folder}")
    process_files_in_folder(sync_folder, output_folder, regions_folder_name)
