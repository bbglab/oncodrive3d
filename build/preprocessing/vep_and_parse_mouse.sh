#!/bin/bash

# Help message
function show_help() {
    printf "Usage: ./vep_and_parse.sh [input_vcf_file] [output_file]\n\n"
    printf "Options:\n"
    printf "  input_vcf_file : Specify the input VCF file.\n"
    printf "  output_file    : Specify the output file (e.g., .maf).\n\n"
    printf "Example:\n"
    printf "  ./vep_and_parse.sh /path/to/input.vcf /path/to/output.vep.tsv\n\n"
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
input_sorted=${input_vcf}.sorted.tsv
output_vep="$2"

# Sort input for VEP
printf "Sorting input...\n"
sort -k1,1 -k2,2n $input_vcf > $input_sorted

# Run VEP
printf "Running VEP...\n"
if ! singularity exec /workspace/datasets/vep/mus_musculus/ensembl-vep_111.0.sif vep --dir /workspace/datasets/vep/mus_musculus/111_GRCm39/ -i "$input_sorted" --offline --cache -o "$output_vep" --species mus_musculus --assembly GRCm39 --fork 8 --symbol --protein --tab --canonical; then
    printf "Error: VEP command failed.\n" >&2
    exit 1
fi

# Parse output
printf "Parsing output...\n"
if ! python3 /workspace/projects/clustering_3d/clustering_3d/build/preprocessing/parse_vep_out.py -i "$output_vep" -o "${output_vep}.gz"; then
    printf "Error: Parsing command failed.\n" >&2
    exit 1
fi
printf "Generated %s\n" "$parsed_output"

# Remove temp files
printf "Removing temp files...\n"
rm $input_sorted
rm -f "${output_vep}_summary.html"
#rm -f "$output_vep"