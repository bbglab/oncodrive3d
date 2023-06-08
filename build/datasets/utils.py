import re
import pandas as pd
import numpy as np
import os
import gzip
import io
import csv
import time
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

# def uniprot_to_hugo(uniprot_ids=None, hugo_as_keys=False, get_ensembl_id=False):
#     """
#     Get human Uniprot IDs to Hugo symbols mapping from www.genenames.org.
#     Returns a dictionary with UniProt IDs as keys and corresponding HUGO symbols as values.
#     """
#     # Download the latest HGNC gene symbols file
#     url = "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_pub_ensembl_id&col=md_prot_id&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
#     response = requests.get(url)
#     csv_data = response.content.decode('utf-8')
    
#     # Parse the CSV data into a dictionary of UniProt IDs to HUGO symbols
#     reader = csv.reader(csv_data.splitlines(), delimiter='\t')
#     next(reader)  # Skip the header row
#     dict_output = {}
#     for row in reader:
#         hugo_symbol, ensembl_id, uniprot_id = row
#         if uniprot_id and hugo_symbol:
            
#             # Get Uniprot ID, HUGO symbol, Ensembl ID
#             if get_ensembl_id:
#                 if not ensembl_id:
#                     ensembl_id = np.nan
                
#                 if hugo_as_keys:
#                     dict_output[hugo_symbol] = uniprot_id, ensembl_id
#                 else:
#                     dict_output[uniprot_id] = hugo_symbol, ensembl_id
                    
#             # Get Uniprot ID, HUGO symbol
#             else:   
#                 if hugo_as_keys:
#                     dict_output[hugo_symbol] = uniprot_id
#                 else:
#                     dict_output[uniprot_id] = hugo_symbol
    
#     # Filter the dictionary to only include UniProt IDs that were input
#     if uniprot_ids is not None and not hugo_as_keys:
#         dict_output = {k: v for k, v in dict_output.items() if k in uniprot_ids}
    
#     return dict_output


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


def get_mapping_jobid(uniprot_ids):
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
            print(f"Requesting ID mapping job to UniprotKB for IDs.. [waited {i-59}s]")
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


def uniprot_to_hudo_df(uniprot_ids):
    """
    Given a list of Uniprot IDs (from any species), request an Id 
    mapping job to UniprotKB to retrieve the corresponding Hugo 
    symbols and additional protein info. Return a pandas dataframe.
    It is recommended to provide batches of IDs up to 5000 elements.
    """
    
    job_id = get_mapping_jobid(uniprot_ids)
    url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job_id}?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=tsv"
    df = load_df_from_url(url)
    
    i = 60
    while df is None:
        time.sleep(1)
        df = load_df_from_url(url)
        if i % 60 == 0:
            print(f"Waiting for UniprotKB mapping job to produce url..")
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


def uniprot_to_hugo(uniprot_ids, hugo_as_keys=False, batch_size=5000):
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
        print(f"Batch {i+1}/{len(uniprot_ids_lst)} ({len(ids)} IDs)..")
        df = uniprot_to_hudo_df(ids)
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