"""
Module including a collection of functions that provide 
general-purpose functionalities that can be used across 
different parts of the dataset building process.
"""

import gzip
import hashlib
import io
import os
import time
import tarfile
from difflib import SequenceMatcher

import daiquiri
import numpy as np
import pandas as pd
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from pypdl import Downloader
from zipfile import ZipFile
import subprocess

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".build.utils")


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


def get_species(species):
    """
    Simply change species name to accepted format.
    """
    
    if species == "human" or species.capitalize() == "Homo sapiens":
        species = "Homo sapiens"
    elif species == "mouse" or species.capitalize() == "Mus musculus": 
        species = "Mus musculus"
    else:
        raise RuntimeError(f"Failed to recognize '{species}' as species. Currently accepted ones are 'Homo sapiens' and 'Mus musculus'. Exiting...")

    return species


# Download

def calculate_hash(filepath: str, hash_func=hashlib.sha256) -> str:
    """
    Calculate the hash of a file using the specified hash function.

    Args:
        filepath (str): The path to the file for which to calculate the hash.
        hash_func (hashlib.Hash): The hash function to use (default is hashlib.sha256).

    Returns:
        str: The hexadecimal representation of the calculated hash.
    """
    
    with open(filepath, 'rb') as file:
        hash_obj = hash_func()
        for chunk in iter(lambda: file.read(8192), b''):
            hash_obj.update(chunk)
    return hash_obj.hexdigest()


def download_single_file(url: str, destination: str, threads: int) -> None:
    """
    Downloads a file from a URL and saves it to the specified destination.

    Args:
        url (str): The URL of the file to download.
        destination (str): The local path where the file will be saved.
    """
    
    num_connections = 40 if threads > 40 else threads

    if os.path.exists(destination):
        logger.debug(f"File {destination} already exists: Skipping download...")
    else:
        logger.debug(f'Downloading from {url}')
        logger.debug(f"Downloading to {destination}")
        dl = Downloader()
        dl.start(url, destination, num_connections=num_connections, display=True)
        subprocess.run("clear")
        logger.debug('clear')
        logger.debug('Download complete')


def extract_tar_file(file_path, path):
     
    checkpoint = os.path.join(path, ".checkpoint.txt")

    if os.path.exists(checkpoint):
         logger.debug('Tar already extracted: skipping...')
    else:
        with tarfile.open(file_path, "r") as tar:
            tar.extractall(path)
            logger.debug(f'Extracted { int(len(tar.getnames())/2)} structure.')
            with open(checkpoint, "w") as f:
                f.write('')
                

def extract_zip_file(file_path, path):
    checkpoint = os.path.join(path, ".zip_checkpoint.txt")

    if os.path.exists(checkpoint):
        logger.debug('ZIP file already extracted: skipping...')
    else:
        with ZipFile(file_path, "r") as zip_file:
            zip_file.extractall(path)
            logger.debug(f'Extracted {len(zip_file.namelist())} files.')
            with open(checkpoint, "w") as f:
                f.write('')
                

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
    if path_structure.endswith("gz"):  
        with gzip.open(path_structure, 'rt') as handle:
            return np.array([record.seq for record in SeqIO.parse(handle, 'pdb-seqres')][0])
    else:
        with open(path_structure, 'r') as handle:
            return np.array([record.seq for record in SeqIO.parse(handle, 'pdb-seqres')][0])


def get_pdb_path_list_from_dir(path_dir):
    """
    Takes as input a path of a given directory and it 
    outputs a list of paths of the contained PDB files.
    """

    pdb_files = os.listdir(path_dir)
    pdb_path_list = [f"{os.path.join(path_dir, f)}" for f in pdb_files if ".pdb" in f and not os.path.isdir(os.path.join(path_dir, f))]
    return pdb_path_list


# Backtranslation

def translate_dna(dna_seq):
    """
    Translate dna sequence to protein.
    """

    if pd.isnull(dna_seq):
        return np.nan    
    else:
        dna_seq = Seq(dna_seq)
        return str(dna_seq.translate()) 


def get_seq_similarity(a, b, decimals=3):
    """
    Compute the similarity ratio between sequences a and b.
    """

    if pd.isnull(a) or pd.isnull(b):
        return np.nan
    else:
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


def get_mapping_jobid(uniprot_ids):
    """
    Submit an ID Mapping job to UniprotKB.
    """
    
    command = f"https://rest.uniprot.org/idmapping/run?from=UniProtKB_AC-ID&to=UniProtKB&ids={','.join(uniprot_ids)}"

    response = "INIT"
    while str(response) != "<Response [200]>":
        if response != "INIT":
            time.sleep(10)
        try:
            response = requests.post(command)
        except requests.exceptions.RequestException as e:                          
            response = "ERROR"                                                     
            logger.debug(f"Request failed: {e}")    
            
    job_id = get_response_jobid(response)
    i = 60
    while job_id is None:
        time.sleep(1) 
        job_id = get_response_jobid(response)
        if i % 60 == 0:
            logger.debug(f"Requesting ID mapping job to UniprotKB for IDs...")
        i += 1
    
    return job_id


def load_df_from_url(url):
    """
    Load a pandas dataframe from url.
    """

    try:
        response = requests.get(url, timeout=(160, 160))
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
            logger.debug(f"Waiting for UniprotKB mapping job to produce url...")
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
        logger.debug(f"Batch {i+1}/{len(uniprot_ids_lst)} ({len(ids)} IDs)...")
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