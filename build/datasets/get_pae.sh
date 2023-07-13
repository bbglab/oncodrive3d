#!/bin/bash

# Help message
function show_help() {
    echo "Usage: ./get_pae.sh [arg1] [arg2]"
    echo ""
    echo "Args:"
    echo "  arg1 : Specify the path of the sequence df (seq_for_mut_prob.csv) storing the Uniprot_IDs."
    echo "  arg2 : Specify the directory where to save the predicted aligned errors."
    echo ""
    echo "Example:"
    echo "./get_pae.sh /workspace/projects/clustering_3d/clustering_3d/datasets_frag/seq_for_mut_prob.csv /workspace/projects/clustering_3d/clustering_3d/datasets_frag/pae"
    echo ""
}

filename="$1"
output_dir="$2"

# Check if the required arguments are provided
if [ -z "$filename" ] || [ -z "$output_dir" ]; then
    show_help
    exit 0
fi

# Create the output directory if it doesn't exist
if [ ! -d "$output_dir" ]; then
    mkdir "$output_dir"
fi

# Read the seq file and iterate through the (non-fragmented) Uniprot_ID
while IFS=',' read -ra cols
do

    # Check if the third element is equal to 1
    if [[ ${cols[2]} == 1 && ${cols[0]} != "Uniprot_ID" ]]; then

        # Download PAE from AF databasae
        wget -q -O $output_dir/${cols[1]}-F1-predicted_aligned_error.json https://alphafold.ebi.ac.uk/files/AF-${cols[1]}-F1-predicted_aligned_error_v4.json

        # Check that the file is not truncated
        if ! tail -n 1 $output_dir/${cols[1]}-F1-predicted_aligned_error.json | grep -q "}]" ; then
            echo The ${cols[1]} is incomplete.. re-downloading..
            wget -q -O $output_dir/${cols[1]}-F1-predicted_aligned_error.json https://alphafold.ebi.ac.uk/files/AF-${cols[1]}-F1-predicted_aligned_error_v4.json

        fi

    fi

done < "$filename"