#!/bin/bash

# Help message
function show_help() {
    echo "Usage: get_pae.sh [arg1] [arg2] [arg3] [arg4]"
    echo ""
    echo "The script download all Predicted Aligned Error files of the non-fragmented PDB structures stored in the given input directory"
    echo "Args:"
    echo "  arg1 : Input directory including the PDB structures"
    echo "  arg2 : Output directory where to download the predicted aligned errors"
    echo "  arg3 : AlphaFold 2 version (4 as default)"
    echo "  arg4 : Verbose (False as default)"
    echo ""
    echo "Example:"
    echo "get_pae.sh /workspace/projects/clustering_3d/clustering_3d/datasets_frag/pdb_structures /workspace/projects/clustering_3d/clustering_3d/datasets_frag/pae"
    echo ""
}

input_dir="$1"
output_dir="$2"
af_version="${3:-4}"
verbose="${4:-'False'}"

# Check if the required arguments are provided
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
    show_help
    exit 0
fi

# Create the output directory if it doesn't exist
if [ ! -d "$output_dir" ]; then
    mkdir "$output_dir"
fi

# Iterate through the PDB structures
path="$input_dir/AF-*-model_v$af_version.pdb"
for file in $path; do 

    # Extract the Uniprot ID and AF fragment
    filename=$(basename "$file")
    uniprot_id=$(echo "$filename" | awk -F'[-.]' '{print $2}')
    af_f=$(echo "$filename" | awk -F'[-.]' '{print gensub(/[^0-9]+/, "", "g", $3)}')

    if [ ! "$af_f" -gt 1 ]; then

        file_path="$output_dir/${uniprot_id}-F1-predicted_aligned_error.json"
        download_url="https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-predicted_aligned_error_v${af_version}.json"
        download_attempts=0
        
        # Try downloading PAE until the file exists and it is not truncated
        while [ ! -f "$file_path" ] || ! tail -n 1 "$file_path" | grep -q "}]"; do
            wget -q -O "$file_path" "$download_url"

            ((download_attempts++))
            if [ $download_attempts -ge 2 ]; then
                if [ "$verbose" = "True" ]; then
                    echo "$file_path"
                    echo "Download attempts: $download_attempts"
                fi
                sleep 0.5 
            fi

        done

    fi

done