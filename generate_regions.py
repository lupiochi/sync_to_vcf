import pandas as pd
import sys
import os

"""
This script will generate regions based on the given points (in a text file), taking into account the sliding window size - that is, what is the distance between the start and end of each point in the chromosome.
The text file must have the first two columns as the chromosome and point, although additional columns might be present and won't affect anything.
For example, you may first sort the chromosome locations by the FST values and then generate regions based on the sorted points.

chromosome	point	fst
chr2	5723000	0.39448213
chr2	5561000	0.34830555
chr2	5725000	0.34593956


The output will be a CSV file for each chromosome with the start and end of each region, based on the sliding window size, as shown below:
start	end
977000	985000
1067000	1071000
1437000	1439000

The script will take into account if consecutive points are within the sliding window size and merge them into a single region.

"""

# Select the appropriate sliding window size here.
genome_window_size = 10000

# Function to generate regions based on given points
def generate_regions(points, window_size = genome_window_size):
    if not points:  # Check if points list is empty
        return []
    points = sorted(points)
    regions = []
    start = points[0]
    end = start + window_size

    for i in range(1, len(points)):
        if points[i] <= end + window_size:
            end = max(end, points[i] + window_size)
        else:
            regions.append((start, end))
            start = points[i]
            end = start + window_size

    regions.append((start, end))
    return regions

# Read input
input_file = sys.argv[1]
region_folder_name = sys.argv[2]
region_input = pd.read_csv(input_file, sep="\t")

# Group data by chromosome and generate regions
grouped = region_input.groupby('chromosome')

# Create the output folder with the specified name (it will be created automatically if it doesn't exist)
output_folder = region_folder_name
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for chromosome, group in grouped:
    points = group['point'].tolist()
    regions = generate_regions(points)
    
    # Create a DataFrame for the regions
    regions_df = pd.DataFrame(regions, columns=['start', 'end'])
    
    # Save to a CSV file named by the chromosome inside the "regions" folder
    output_file = os.path.join(output_folder, f"{chromosome}_regions.csv")
    regions_df.to_csv(output_file, index=False)

print(f"Regions files generated successfully. Please find them in the {region_folder_name} folder.")