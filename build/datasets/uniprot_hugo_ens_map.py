"""
Given a directory including AF predicted structure, it generate a dictionary 
having Uniprot IDs as keys and HUGO symbols as values. If genes_as_keys equal 1,
use HUGO symbols as value and Uniprot IDs as keys. If ensebl_id equal 1, use 
a tuple (Uniprot/HUGO, ENSEMBL ID) as value.

###################################### EXAMPLE USAGE ################################################

python3 uniprot_hugo_ens_map.py -i /workspace/datasets/alphafold_features/AF_homo_sapiens_pred/ \
-o /../../datasets/temp \
-f af_uniprot_to_gene_id.json \
-g 0 -e 0

WARNING: it can overwrite hand-curated dictionary for mapping

#####################################################################################################
"""


import argparse
import json
import os
from utils import uniprot_to_hugo


def main():
    
    ## Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input directory with PDB structures", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output path", type=str, default="/../../datasets/temp")         
    parser.add_argument("-f", "--filename", help="Dictionary filename", type=str, default="af_uniprot_to_gene_id.json")                                  
    parser.add_argument("-g", "--genes_as_keys", help="if 1 use HUGO symbols as keys, if 0 use Uniprot ID as keys", type=int, default=0)   
    parser.add_argument("-e", "--ensembl_id", help="if 1 include Ensembl ID as value", type=int, default=0)
    
    args = parser.parse_args()
    input = args.input
    output = args.output
    filename = args.filename
    genes_as_keys = args.genes_as_keys
    ensembl_id = args.ensembl_id

    print("\nGet dictionary for Uniprot ID to HUGO symbol /+ ENSEMBL ID mapping")
    print("\nInput directory:", input)
    print("Output:", output)
    print("Use HUGO symbol as keys:", genes_as_keys)
    print("Include Ensembl ID as values:", ensembl_id)

    # Get mapping
    dict_mapping = uniprot_to_hugo(hugo_as_keys=genes_as_keys, get_ensembl_id=ensembl_id)

    # Save
    if not os.path.exists(output):
        os.makedirs(output)
    output = f"{output}/{filename}"
    with open(output, "w") as fp:
        json.dump(dict_mapping, fp)
    print("\nFile saved")


if __name__ == "__main__":
    main()
