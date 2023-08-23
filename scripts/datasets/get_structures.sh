#!/bin/bash


# Help message
function show_help() {
    echo "Script to download all AlphaFold predicted structures (PDB files) for a given organism from AlphaFold DB"
    echo ""
    echo "Usage: get_structures.sh [arg1] [arg2] [arg3]"
    echo ""
    echo "Options:"
    echo "  arg1 : Path where to save PDB structures."
    echo "  arg2 : Species (human (default) or mouse)."
    echo "  arg3 : AlphaFold 2 version (4 as default)"
    echo "  arg4 : Verbose (True (default) or False)"
    echo ""
    echo "Example:"
    echo "get_structures.sh ../../datasets/pdb_structures human True"
    echo ""
}

# Check if the help flag is provided
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

# Check if option 1 is provided
if [ -z "$1" ]; then
    echo "Error: Missing option1. Specify the path where to save PDB structures."                   #### >> Raise error
    show_help
    exit 1
fi

# Assign default value to empty variables
species="${2:-'human'}"
af_version="${3:-4}"
verbose="${4:-'False'}"

# Select proteome
if [ "$species" == "human" ]; then
    proteome=UP000005640_9606_HUMAN_v$af_version
elif [ "$species" == "mouse" ]; then
    proteome=UP000000589_10090_MOUSE_v$af_version
fi

# Download proteome
AF_URL=https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/$proteome.tar
rm -rf $1/*.tar

if [ "$verbose" = "True" ]; then
    echo "[INFO]: Selected species: $species"                                                      #### >> INFO
    echo "[INFO]: Proteome to download: $proteome"                                                 #### >> INFO
    wget -P $1 $AF_URL
    #timeout 2m wget -P $1 $AF_URL || true                                    ################## FOR TESTING  ###########
else
    wget -q -P $1 $AF_URL > /dev/null 2>&1
    #timeout 2m wget -q -P $1 $AF_URL > /dev/null 2>&1 || true                ################## FOR TESTING  ###########
fi

# Unzip files
if [ -n "$(ls -A $1/*.tar 2>/dev/null)" ]; then
    tar xC $1 -f $1/*.tar
    find $1 -name "*pdb.gz" -exec gunzip "{}" \;
    find "$1" \( -name "*.tar" -o -name "*.gz" \) -type f -delete
else
    echo "Could not download .tar file from $AF_URL"                                          #### >> Raise error
fi

# Add list of downloaded structures to txt file
if [ -n "$(ls -A $1/*.pdb 2>/dev/null)" ]; then
    ls $1/*.pdb | sed "s/.*\///" > $1/pdb_lst.txt
    if [ "$verbose" = "True" ]; then
        echo "[INFO]: Structures downloaded in directory $1"                                  #### >> INFO
    fi
    #touch $1/.checkpoint
else
    echo "No .pdb files found in directory $1"                                    #### >> Raise error
fi