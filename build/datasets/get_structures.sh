#!/bin/bash

#################################### FIX THE ERROR #######################

# Help message
function show_help() {
    echo "Usage: ./get_structures.sh [option1] [option2]"
    echo ""
    echo "Options:"
    echo "  option1 : Specify the path where to save PDB structures."
    echo "  option2 : Specify the species (human or mouse)."
    echo ""
    echo "Example:"
    echo "./script.sh ../../datasets_mouse/pdb_structures mouse"
    echo ""
}

# Check if the help flag is provided
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

# Check if option1 is provided
if [ -z "$1" ]; then
    echo "Error: Missing option1. Specify the path where to save PDB structures."
    show_help
    exit 1
fi

# Assign default value if $2 is empty
if [ -z "$2" ]; then
    species=human
else
    species=$2
fi

# Perform actions based on $2
if [ "$species" == "human" ]; then
    proteome=UP000005640_9606_HUMAN_v4
elif [ "$species" == "mouse" ]; then
    proteome=UP000000589_10090_MOUSE_v4
fi

# Check if valid species is provided
if [[ "$species" == "mouse" || "$species" == "human" ]]; then
    echo "Valid species: $species"
    echo "Proteome to download: $proteome"
    echo ""

    AF_URL=https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/$proteome.tar
    rm -rf $1/*.tar
    wget -P $1 $AF_URL

    if [ -n "$(ls -A $1/*.tar 2>/dev/null)" ]; then
        tar xC $1 -f $1/*.tar
        find $1 -name "*pdb.gz" -exec gunzip "{}" \;
        find "$1" \( -name "*.tar" -o -name "*.gz" \) -type f -delete
    else
        echo "Could not download .tar file from $AF_URL"
    fi

    if [ -n "$(ls -A $1/*.pdb 2>/dev/null)" ]; then
        ls $1/*.pdb | sed "s/.*\///" > $1/pdb_lst.txt
        echo "Structures downloaded in directory $1"
        touch $1/.checkpoint
    else
        echo "No .pdb files found in directory $1"
    fi

else
    echo "Invalid species: $species"
fi