
# SYNC to VCF Conversion

This repo consists of three main scripts designed to process genomic data from [SYNC files](https://sourceforge.net/projects/popoolation2/), generate regions, convert SYNC files to VCF format, and produce consolidated output files.
SYNC files are a popular data format that encodes allele counts for each population sample [Kapun et al. 2021](https://academic.oup.com/mbe/article/38/12/5782/6361628), [Kofler et al. 2011](https://academic.oup.com/bioinformatics/article/27/24/3435/306737). They are integral to the study of population genetics, as they provide detailed allele frequency data from pooled sequencing experiments. Tools like PoPoolation2 rely on SYNC files to enable comprehensive analysis of genetic diversity, evolutionary patterns, and population structure. By converting SYNC files to the widely-used [Variant Call Format (VCF)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format) and generating a summarized output, this brief tool facilitates more accessible and extensive genomic analyses.

### 1. generate_regions.py

This script will output regions based on the given points of interest, taking into account the sliding window size - that is, what is the distance between the start and end of each point in the chromosome.
The text file must have the first two columns as the chromosome and point, although additional columns might be present and won't affect anything.
For example, you may first sort the chromosome locations by the FST values and then generate regions based on the sorted points.

The output will be a CSV file for each chromosome with the start and end of each region.

The script will take into account if consecutive points are within the sliding window size and merge them into a single region.

### Usage:
```bash
python generate_regions.py <input_file> <region_folder_name>
```

### Arguments:
- `<input_file>`: Path to the input text file containing chromosome locations.
- `<region_folder_name>`: Name of the output folder to save the generated region CSV files.

### Example input file format:
```
chromosome    point    fst
chr2          5723000  0.39448213
chr2          5561000  0.34830555
chr2          5725000  0.34593956
```

### Output format:
```
start    end
977000   985000
1067000  1071000
1437000  1439000
```

### 2. sync_to_vcf.py

This script converts a SYNC file into a VCF file. The script processes the SYNC file to generate allele frequencies and other relevant data, and outputs a VCF file.

### Usage:
```bash
python sync_to_vcf.py <sync_file> <output_folder_name> [regions_file]
```

### Arguments:
- `<sync_file>`: Path to the SYNC file.
- `<output_folder_name>`: Name of the output folder to save the VCF and CSV files.
- `[regions_file]`: Optional path to a regions file to filter the SYNC file by specific regions.

### Example SYNC file format:
```
chr2    84    A    30:0:0:0:0:0    30:0:0:0:0:0    30:0:0:0:0:0
chr2    85    T    0:30:0:0:0:0    0:30:0:0:0:0    0:30:0:0:0:0
```

### Output VCF format:
```
#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT    POP1    POP2    POP3    pop1_pop3_diff
chr2      977029 .      G     A      30      .        DP=30    AF        1.0     1.0     1.0     0.0
```

The script can add columns for the absolute difference between allele frequencies of specified populations and filter out rows where the ALT allele is the same as the REF allele.

### 3. multi_sync_vcf.py

This master script automates the conversion of multiple SYNC files to VCF format using the sync_to_vcf.py script, concatenates the results, and produces consolidated output files. SYNC files are derived from pooled sequencing data and contain information about allele frequencies at various positions in the genome, making them crucial in population genetics studies.

### Usage:
```bash
python multi_sync_vcf.py <sync_folder> <output_folder_name> [regions_folder_name]
```

### Arguments:
- `<sync_folder>`: Path to the folder containing SYNC files.
- `<output_folder_name>`: Name of the output folder to save the results.
- `[regions_folder_name]`: Optional path to a folder containing regions files to filter the SYNC files.

### multi_sync_vcf.py will follow the following steps:
1. Checks if the output folder exists and creates it if necessary.
2. Iterates over all SYNC files in the specified folder.
3. For each SYNC file, constructs a command to convert it to VCF using sync_to_vcf.py, optionally including regions files.
4. Executes the command and handles any errors during the conversion.
5. Loads the resulting CSV files and appends them to a list of dataframes.
6. Reads the generated VCF files, extracts the headers, and appends the data to a temporary concatenated VCF file.
7. Sanitizes and saves the final concatenated VCF file.
8. Concatenates all dataframes and saves the final concatenated CSV file.


The script will output a concatenated VCF file and a concatenated CSV file in the specified output folder.

This pipeline ensures efficient processing and conversion of genomic data from SYNC files to VCF and CSV formats, facilitating further analysis and research in population genetics.

### Developed by Luiz F. Piochi | [Website](https://lupiochi.github.io/)
