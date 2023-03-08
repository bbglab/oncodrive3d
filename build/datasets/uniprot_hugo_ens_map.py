"""
Given a directory including AF predicted structure, it generate a dictionary 
having Uniprot IDs as keys and HUGO symbols as values. If genes_as_keys equal 1,
use HUGO symbols as value and Uniprot IDs as keys. If ensebl_id equal 1, use 
a tuple (Uniprot/HUGO, ENSEMBL ID) as value.

###################################### EXAMPLE USAGE ################################################

python3 uniprot_hugo_ens_map.py -i /workspace/datasets/alphafold_features/AF_homo_sapiens_pred/ \
-o /../../datasets/af_uniprot_to_gene_id.json \
-g 0 -e 0

WARNING: it can overwrite hand-curated dictionary for mapping

#####################################################################################################
"""


import argparse
import json
from utils import get_pdb_path_list_from_dir, uniprot_to_hugo


def main():
    
    ## Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input directory with PDB structures", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output path for the Uniprot-HUGO-Ensembl dictionary", type=str, default="/../../datasets/af_uniprot_to_gene_id.json")                                        
    parser.add_argument("-g", "--genes_as_keys", help="if 1 use HUGO symbols as keys, if 0 use Uniprot ID as keys", type=int, default=0)   
    parser.add_argument("-e", "--ensembl_id", help="if 1 include Ensembl ID as value", type=int, default=0)
    
    args = parser.parse_args()
    input = args.input
    output = args.output
    genes_as_keys = args.genes_as_keys
    ensembl_id = args.ensembl_id

    print("\nGet dictionary for Uniprot ID to HUGO symbol /+ ENSEMBL ID mapping")
    print("\nInput directory:", input)
    print("Output:", output)
    print("Use HUGO symbol as keys:", genes_as_keys)
    print("Include Ensembl ID as values:", ensembl_id)

    # Run
    # list_prot_path = get_pdb_path_list_from_dir(input)
    # uniprot_ids = [p.split("AF-")[1].split("-model_v1")[0].split("-F")[0] for p in list_prot_path]
    dict_mapping = uniprot_to_hugo(hugo_as_keys=genes_as_keys, get_ensembl_id=ensembl_id)

    with open(output, "w") as fp:
        json.dump(dict_mapping, fp)
    print("\nFile saved")


if __name__ == "__main__":
    main()
