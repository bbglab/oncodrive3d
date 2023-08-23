"""
Module including a collection of functions that provide 
general-purpose functionalities that can be used across 
different parts of the dataset building process.
"""


import re
import pandas as pd
import numpy as np
import os
import gzip
import io
import time
import requests
from difflib import SequenceMatcher
from Bio import SeqIO
from Bio.Seq import Seq


# General utils

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


# PDB

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
    pdb_path_list = [f"{path_dir}/{f}" for f in pdb_files if re.search('.\.pdb$', f) is not None]
    return pdb_path_list


# Backtranslation

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


# Uniprot ID to Hugo symbols mapping

def get_response_jobid(response):
    """
    Get jobId after submitting ID Mapping job to UniprotKB. 
    """
    
    try:
        data = response.json()
        job_id = data.get("jobId")
    except:
        job_id = None
        
    return job_id


def get_mapping_jobid(uniprot_ids, verbose):
    """
    Submit an ID Mapping job to UniprotKB.
    """
    
    command = f"https://rest.uniprot.org/idmapping/run?from=UniProtKB_AC-ID&to=UniProtKB&ids={','.join(uniprot_ids)}"
    response = requests.post(command)
    job_id = get_response_jobid(response)
    
    i = 60
    while job_id is None:
        time.sleep(1) 
        job_id = get_response_jobid(response)
        if i % 60 == 0:
            if verbose: print(f"Requesting ID mapping job to UniprotKB for IDs.. [waited {i-59}s]")
        i += 1
    
    return job_id


def load_df_from_url(url):
    """
    Load a pandas dataframe from url.
    """

    try:
        response = requests.get(url)
        decompressed_data = gzip.decompress(response.content)
        df = pd.read_csv(io.BytesIO(decompressed_data), sep='\t')
    except:
        df = None
    
    return df


def split_lst_into_chunks(lst, batch_size = 5000):
    """
    Simple split a list into list of list of chunk_size elements.
    """
    
    return [lst[i:i+batch_size] for i in range(0, len(lst), batch_size)]


def uniprot_to_hudo_df(uniprot_ids, verbose):
    """
    Given a list of Uniprot IDs (from any species), request an Id 
    mapping job to UniprotKB to retrieve the corresponding Hugo 
    symbols and additional protein info. Return a pandas dataframe.
    It is recommended to provide batches of IDs up to 5000 elements.
    """
    
    job_id = get_mapping_jobid(uniprot_ids, verbose)
    url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job_id}?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=tsv"
    df = load_df_from_url(url)
    
    i = 60
    while df is None:
        time.sleep(1)
        df = load_df_from_url(url)
        if i % 60 == 0:
            if verbose: print(f"Waiting for UniprotKB mapping job to produce url..")
        i += 1 
        
    return df


def convert_dict_hugo_to_uniprot(dict_uniprot_hugo):
    """
    Convert a Uniprot IDs to Hugo symbol dictionary to a Hugo symbo to 
    Uniprot IDs dictionary, if multiple Hugo symbols are mapped to the 
    same Uniprot ID, add them as multiple keys.
    """
    
    dict_hugo_uniprot = {}

    for uni_id, gene in dict_uniprot_hugo.items():
        if type(gene) == str:
            for g in gene.split(" "):
                dict_hugo_uniprot[g] = uni_id
                
    return dict_hugo_uniprot


def uniprot_to_hugo(uniprot_ids, hugo_as_keys=False, batch_size=5000, verbose=False):
    """
    Given a list of Uniprot IDs (any species.), request an Id mapping 
    job to UniprotKB to retrieve the corresponding Hugo symbols. 
    Return a dictionary of Uniprot IDs to Hugo symbols or vice versa.
    """
    
    # Split uniprot IDs into chunks
    uniprot_ids_lst = split_lst_into_chunks(uniprot_ids, batch_size)
    
    # Get a dataframe including all IDs mapping info
    df_lst = []
    for i, ids in enumerate(uniprot_ids_lst):
        if verbose: print(f"Batch {i+1}/{len(uniprot_ids_lst)} ({len(ids)} IDs)..")
        df = uniprot_to_hudo_df(ids, verbose)
        df_lst.append(df)
    df = pd.concat(df_lst)

    # Get a dictionary for Uniprot ID to Hugo symbols
    dictio = {}
    for i, r in df[["Entry", "Gene Names"]].iterrows():
        uni_id, gene = r 
        dictio[uni_id] = gene
    
    # Convert to a dictionary of Hugo symbols to Uniprot IDs
    if hugo_as_keys:
        dictio = convert_dict_hugo_to_uniprot(dictio)
            
    return dictio