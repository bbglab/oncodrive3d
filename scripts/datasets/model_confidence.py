"""
Module to get per-residue model confidence from all AlphaFold
predicted structures contained in a given directory.
"""


import gzip
import os

import daiquiri
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB.PDBParser import PDBParser
from tqdm import tqdm

from scripts import __logger_name__
from scripts.datasets.utils import get_pdb_path_list_from_dir

logger = daiquiri.getLogger(__logger_name__ + ".build.model_confidence")


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


def get_confidence(input, output_dir, seq_df):
    """
    Get per-residue model confidence from all AlphaFold
    predicted structures contained in a given directory.
    """

    checkpoint = os.path.join(output_dir, '.checkpoint.conf.txt')
    if os.path.exists(checkpoint):
        logger.debug("Confidence extraction performed: Skipping...")

    else:
        output = os.path.join(output_dir, 'confidence.tsv')

        logger.debug(f"Input directory: {input}")
        logger.debug(f"Output: {output}")

        # Get model confidence
        df_list = []
        pdb_path_list = get_pdb_path_list_from_dir(input)

        for file in tqdm(pdb_path_list, total=len(pdb_path_list), desc="Extracting model confidence"):
            try:
                identifier = file.split("AF-")[1].split("-model")[0].split("-F")
            except Exception as e:
                logger.warning(f'Could not extract Uniprot ID from {file}, {e}')

            if identifier[0] in seq_df["Uniprot_ID"].values:

                parser = PDBParser()
                if file.endswith(".gz"):
                    with gzip.open(file, 'rt') as handle:
                        structure = parser.get_structure("ID", handle)
                else:
                    with open(file, 'r') as handle:
                        structure = parser.get_structure("ID", handle)
                chain = structure[0]["A"]

                # Get confidence
                confidence_df = get_confidence_one_chain(chain).reset_index().rename(columns={"index": "Pos"})
                confidence_df["Pos"] = confidence_df["Pos"] + 1
                confidence_df["Uniprot_ID"] = identifier[0]
                confidence_df["AF_F"] = identifier[1]
                df_list.append(confidence_df)

        confidence_df = pd.concat(df_list).reset_index(drop=True)
        confidence_df.to_csv(output, index=False, sep="\t")

        with open(checkpoint, "w") as f:
                    f.write('')

        logger.info("Extraction of model confidence completed!")