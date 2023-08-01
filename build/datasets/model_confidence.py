"""
Get per-residue model confidence from all AlphaFold predicted structures contained in a given directory.

###################################### EXAMPLE USAGE ################################################

python3 model_confidence.py -i ../../datasets/pdb_structures/ -o ../../datasets/confidence.csv

#####################################################################################################
"""


import click
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.PDBParser import PDBParser
from utils import get_pdb_path_list_from_dir
from progressbar import progressbar


def get_3to1_protein_id(protein_id):
    """
    Convert a 3 letter protein code into 1 letter.
    """
    return protein_letters_3to1[protein_id.lower().capitalize()]


def get_confidence(chain):
    """
    Get AF model confidence from its predicted structure.
    """
    
    res_ids = []
    confidence_scores = []
    
    # Iterate through the chain
    for pos, res in enumerate(chain):
        res_id = get_3to1_protein_id(res.resname)
        confidence = res["CA"].get_bfactor()
        res_ids.append(res_id)
        confidence_scores.append(confidence)
        
    return pd.DataFrame({"Res" : res_ids, "Confidence" : confidence_scores})


def get_confidence_from_dir(path_dir):
    """
    Get per-residue model confidence from all AlphaFold 
    predicted structures contained in a given directory.
    """
    
    df_list = []
    pdb_path_list = get_pdb_path_list_from_dir(path_dir)
    pdb_path_list = [f"{path_dir}{file}" for file in pdb_path_list]
    
    for file in progressbar(pdb_path_list):
            identifier = file.split("AF-")[1].split("-model")[0].split("-F")
            
            # Get chain
            parser = PDBParser()
            structure = parser.get_structure("ID", file)
            chain = structure[0]["A"]

            # Get confidence
            confidence_df = get_confidence(chain).reset_index().rename(columns={"index": "Pos"})
            confidence_df["Pos"] = confidence_df["Pos"] + 1
            confidence_df["Uniprot_ID"] = identifier[0]
            confidence_df["AF_F"] = identifier[1]            
            df_list.append(confidence_df)
    
    return pd.concat(df_list).reset_index(drop=True)


@click.command(context_settings=dict(help_option_names=['-h', '--help']),
               help='Get per-residue model confidence from all AlphaFold predicted structures contained in a given directory.')
@click.option("-i", "--input", type=click.Path(exists=True), required=True, help="Input directory with PDB structures")
@click.option("-o", "--output", help="Output path", default="../../datasets/confidence.csv")
def main(input, output):
    
    print("\nExtracting model confidence from PDB structures..")
    print("\nInput directory:", input)
    print("Output:", output)
    
    # Get model confidence
    confidence_df = get_confidence_from_dir(input)
    confidence_df.to_csv(output, index=False)

if __name__ == "__main__":
    main()
