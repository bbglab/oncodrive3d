#!/bin/bash

# Help message
function show_help() {
    echo "Usage: ./vep_and_parse.sh [option1] [option2]"
    echo ""
    echo "Options:"
    echo "  option1 : Specify the input vcf file."
    echo "  option2 : Specify the output file (e.g. .maf)."
    echo ""
    echo "Example:"
    echo "/workspace/projects/clustering_3d/clustering_3d/build/preprocessing/vep_and_parse_mouse.sh /input.vcf /output.maf"
    echo ""
}

# Check if the help flag is provided
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

# Check if arguments are missing
if [[ -z "$1" || -z "$2" ]]; then
    echo "Error: Missing option1 or option2."
    show_help
    exit 0
fi

# Run VEP
echo "Running VEP.."
singularity exec /workspace/datasets/vep/mus_musculus/vep102.sif vep --dir /workspace/datasets/vep/ -i $1 --offline --cache -o $2 --species mus_musculus --assembly GRCm38 --fork 8 --symbol --protein --tab --canonical --pick
# Remove VEP summary
rm $2"_summary.html"

# Remove comments in the header
echo $2".temp"
cat $2 | grep -v "^##" > $2".temp"
# mv "$temp_file" "$2"

# Parse output
echo "Parsing output.."
python3 /workspace/projects/clustering_3d/clustering_3d/build/preprocessing/parse_vep_out.py -i $2".temp" -o $2
rm $2