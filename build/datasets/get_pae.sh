#!/bin/bash

# Help message
function show_help() {
    echo "Usage: ./get_pae.sh [arg1] [arg2]"
    echo ""
    echo "The script download all Predicted Aligned Error files of the non-fragmented PDB structures stored in the given input directory"
    echo "Args:"
    echo "  arg1 : Specify the path of the sequence df (seq_for_mut_prob.csv) storing the Uniprot_IDs"
    echo "  arg1 : Specify the input directory including the PDB structures"
    echo "  arg2 : Specify the output directory where to download the predicted aligned errors"
    echo ""
    echo "Example:"
    echo "./get_pae.sh /workspace/projects/clustering_3d/clustering_3d/datasets_frag/pdb_structures /workspace/projects/clustering_3d/clustering_3d/datasets_frag/pae"
    echo ""
}

input_dir="$1"
output_dir="$2"

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
for file in "$input_dir"/AF-*-model_v4.pdb; do

    filename=$(basename "$file") 

    # Extract the Uniprot ID and AF fragment
    uniprot_id=$(echo "$filename" | awk -F'[-.]' '{print $2}')
    af_f=$(echo "$filename" | awk -F'[-.]' '{print gensub(/[^0-9]+/, "", "g", $3)}')

    if [ ! "$af_f" -gt 1 ]; then

        file_path="$output_dir/${uniprot_id}-F1-predicted_aligned_error.json"
        download_url="https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-predicted_aligned_error_v4.json"
        download_attempts=0
        
        # Try downloading PAE until the file exists and it is not truncated
        while [ ! -f "$file_path" ] || ! tail -n 1 "$file_path" | grep -q "}]"; do
            wget -q -O "$file_path" "$download_url"

            ((download_attempts++))
            if [ $download_attempts -ge 2 ]; then
                echo "$file_path"
                echo "Download attempts: $download_attempts"
                sleep 0.5 
            fi

        done

    fi

done
