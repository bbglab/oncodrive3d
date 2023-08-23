"""
Module to get per-residue model confidence from all AlphaFold 
predicted structures contained in a given directory.
"""


import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.PDBParser import PDBParser
from progressbar import progressbar
import os
from scripts.datasets.utils import get_pdb_path_list_from_dir


def get_3to1_protein_id(protein_id):
    """
    Convert a 3 letter protein code into 1 letter.
    """
    return protein_letters_3to1[protein_id.lower().capitalize()]


def get_confidence_one_chain(chain):
    """
    Get AF model confidence from its predicted structure.
    """
    
    res_ids = []
    confidence_scores = []
    
    # Iterate through the chain
    for res in chain:
        res_id = get_3to1_protein_id(res.resname)
        confidence = res["CA"].get_bfactor()
        res_ids.append(res_id)
        confidence_scores.append(confidence)
        
    return pd.DataFrame({"Res" : res_ids, "Confidence" : confidence_scores})


def get_confidence(input, output, verbose):
    """
    Get per-residue model confidence from all AlphaFold 
    predicted structures contained in a given directory.
    """
    
    if output is None:                                                          ##### TO CHANGE WHEN CACHE PATH TO DATASETS WILL BE FINALIZED
        dir_path = os.path.abspath(os.path.dirname(__file__))
        output = f"{dir_path}/../../datasets/confidence.csv"
        
    if  verbose:
        print("Input directory:", input)
        print("Output:", output)
    
    # Get model confidence

    df_list = []
    pdb_path_list = get_pdb_path_list_from_dir(input)
    
    for file in progressbar(pdb_path_list):
        identifier = file.split("AF-")[1].split("-model")[0].split("-F")
        
        # Get chain
        parser = PDBParser()
        structure = parser.get_structure("ID", file)
        chain = structure[0]["A"]

        # Get confidence
        confidence_df = get_confidence_one_chain(chain).reset_index().rename(columns={"index": "Pos"})
        confidence_df["Pos"] = confidence_df["Pos"] + 1
        confidence_df["Uniprot_ID"] = identifier[0]
        confidence_df["AF_F"] = identifier[1]            
        df_list.append(confidence_df)
        
    confidence_df = pd.concat(df_list).reset_index(drop=True)
    confidence_df.to_csv(output, index=False)