import re
import pandas as pd
import numpy as np
import os
import csv
import requests
from difflib import SequenceMatcher
from Bio import SeqIO
from Bio.Seq import Seq



def rounded_up_division(num, den):
    """
    Simply round up the result of the division.
    """
    
    return -(-num // den)


def get_pos_fragments(mut_gene_df):
    """
    Get the corresponding fragment of each position of the protein.
    """
    
    max_f = rounded_up_division(max(mut_gene_df["Pos"]), 1400)
    bins = [n * 1400 for n in range(0, max_f+1)]
    group_names = list(range(1, max_f+1))
    
    return pd.cut(mut_gene_df["Pos"], bins, labels = group_names)


def get_af_id_from_pdb(path_structure):
    """
    Get AlphaFold 2 identifier (UniprotID_F) from path
    """

    return path_structure.split("AF-")[1].split("-model")[0]   


def get_seq_from_pdb(path_structure):
    """
    Get sequense of amino acid residues from PDB structure.
    """
    
    return np.array([record.seq for record in SeqIO.parse(path_structure, 'pdb-seqres')][0])


def get_pdb_path_list_from_dir(path_dir):
    """
    Takes as input a path of a given directory and it 
    outputs a list of paths of the contained PDB files.
    """

    pdb_files = os.listdir(path_dir)
    pdb_path_list = [f"{path_dir}{f}" for f in pdb_files if re.search('.\.pdb$', f) is not None]
    return pdb_path_list


def uniprot_to_hugo(uniprot_ids=None, hugo_as_keys=False, get_ensembl_id=False):
    """
    Returns a dictionary with UniProt IDs as keys and corresponding HUGO symbols as values.
    """
    # Download the latest HGNC gene symbols file
    url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_pub_ensembl_id&col=md_prot_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
    response = requests.get(url)
    csv_data = response.content.decode('utf-8')
    
    # Parse the CSV data into a dictionary of UniProt IDs to HUGO symbols
    reader = csv.reader(csv_data.splitlines(), delimiter='\t')
    next(reader)  # Skip the header row
    dict_output = {}
    for row in reader:
        hugo_symbol, ensembl_id, uniprot_id = row
        if uniprot_id and hugo_symbol:
            
            # Get Uniprot ID, HUGO symbol, Ensembl ID
            if get_ensembl_id:
                if not ensembl_id:
                    ensembl_id = np.nan
                
                if hugo_as_keys:
                    dict_output[hugo_symbol] = uniprot_id, ensembl_id
                else:
                    dict_output[uniprot_id] = hugo_symbol, ensembl_id
                    
            # Get Uniprot ID, HUGO symbol
            else:
                if hugo_as_keys:
                    dict_output[hugo_symbol] = uniprot_id
                else:
                    dict_output[uniprot_id] = hugo_symbol
    
    # Filter the dictionary to only include UniProt IDs that were input
    if uniprot_ids is not None and not hugo_as_keys:
        dict_output = {k: v for k, v in dict_output.items() if k in uniprot_ids}
    
    return dict_output


def translate_dna(dna_seq):
    """
    Translate dna sequence to protein.
    """
    
    dna_seq = Seq(dna_seq)
    return str(dna_seq.translate()) 


def get_seq_similarity(a, b, decimals=3):
    """
    Compute the similarity ratio between sequences a and b.
    """
    
    return round(SequenceMatcher(a=a, b=b).ratio(), decimals)