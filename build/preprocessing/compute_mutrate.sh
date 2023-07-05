#!/bin/bash

###########################
# Example
# ./bgsign.sh <input_dir_to_maf_files> <output_mut_profile_dir> <regions_file> <genome> <cores>
# ./compute_mutrate.sh /workspace/projects/clustering_3d/evaluation/datasets/input/maf/ /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile_bgsign/mut_profile /workspace/projects/clustering_3d/evaluation/datasets/input/mut_profile_bgsign/hg38_wg_regions.tsv hg38 10

# To get region file
# python3 /workspace/projects/intogen_plus/fixdatasets-20230223/intogen-plus/build/datasets/regions/create_wg_regions.py hg38 3 > hg38_wg_regions.tsv
###########################

input_dir="$1"
output_dir="$2"
regions="$3"
genome="$4"
cores="$5"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Check if the count file already exists
if [ ! -f "$output_dir/wg.count.json" ]; then
    # Run the bgsignature count command to generate the output file
    bgsignature count -r "$regions" -s 3 -g "$genome" --cores "$cores" --collapse --exclude-N -o "$output_dir/wg.count.json"
fi

# Loop through each file in the input directory
for file in "$input_dir"/*; do

    # Get the filename without the directory path
    filename=$(basename "$file")
    filename="${filename%%.in.maf}"

    # Get count of mutation per context and normalize by trinucleotide bias (mut profile)
    bgsignature normalize -m "$file" -r "$regions" --normalize "$output_dir/wg.count.json" -s 3 -g "$genome" --collapse --cores 4 -o "$output_dir/${filename}.mutrate.json"

done