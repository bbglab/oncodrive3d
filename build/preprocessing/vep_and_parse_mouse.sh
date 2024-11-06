#!/bin/bash

# Help message
function show_help() {
    printf "Usage: ./vep_and_parse.sh [input_vcf_file] [output_file]\n\n"
    printf "Options:\n"
    printf "  input_vcf_file : Specify the input VCF file.\n"
    printf "  output_file    : Specify the output file (e.g., .maf).\n\n"
    printf "Example:\n"
    printf "  ./vep_and_parse.sh /path/to/input.vcf /path/to/output.vep\n\n"
}

# Check if the help flag is provided
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

# Check if arguments are missing
if [[ -z "$1" || -z "$2" ]]; then
    printf "Error: Missing input_vcf_file or output_file.\n" >&2
    show_help
    exit 1
fi

# Assign arguments to variables for clarity
input_vcf="$1"
output_vep="$2"
temp_output="${output_vep}.temp"
parsed_output="${output_vep}.in.maf.tsv"

# Run VEP
printf "Running VEP...\n"
if ! singularity exec /workspace/datasets/vep/mus_musculus/ensembl-vep_111.0.sif vep --dir /workspace/datasets/vep/ -i "$input_vcf" --offline --cache -o "$output_vep" --species mus_musculus --assembly GRCm39 --fork 8 --symbol --protein --tab --canonical --pick; then
    printf "Error: VEP command failed.\n" >&2
    exit 1
fi

# Remove VEP summary
rm -f "${output_vep}_summary.html"

# Remove comments in the header and create a temporary file
if ! grep -v "^##" "$output_vep" > "$temp_output"; then
    printf "Error: Failed to process VEP output.\n" >&2
    exit 1
fi
printf "Generated %s\n" "$temp_output"

# Parse output
printf "Parsing output...\n"
if ! python3 /workspace/projects/clustering_3d/clustering_3d/build/preprocessing/parse_vep_out.py -i "$temp_output" -o "$parsed_output"; then
    printf "Error: Parsing command failed.\n" >&2
    exit 1
fi
printf "Generated %s\n" "$parsed_output"

# Clean up temporary file
rm -f "$temp_output"