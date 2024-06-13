import logging
import os

import daiquiri
import pandas as pd
import numpy as np
import re
import glob
from progressbar import progressbar
from multiprocessing import Pool
import json

from scripts import __logger_name__
from scripts.datasets.utils import download_single_file, extract_zip_file
from scripts.globals import rm_dir

logger = daiquiri.getLogger(__logger_name__ + ".plotting.stability_change")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)



# ===============================
# Stability change upon mutations
# ===============================


def download_stability_change(path: str,
                              threads: int = 1):
    """
    Downloads stability change upon mutations predicted on AlphaFold 
    structures by RaSP.
    
    Rapid protein stability prediction using deep learning representations
    https://elifesciences.org/articles/82593
    DOI: 10.7554/eLife.82593
    """

    url_website = "https://sid.erda.dk/cgi-sid/ls.py?share_id=fFPJWflLeE"
    filename = "rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.zip"
    download_url = "https://sid.erda.dk/share_redirect/fFPJWflLeE/rasp_preds_alphafold_UP000005640_9606_HUMAN_v2.zip"

    logger.debug(f"Filename: {filename}")
    logger.debug(f"Website url: {url_website}")
    file_path = os.path.join(path, filename)

    try:
        # Download file
        logger.debug(f'Downloading to {file_path}')
        download_single_file(download_url, file_path, threads)
        
        # Extract from zip
        logger.debug(f'Extracting {filename}')
        extract_zip_file(file_path, path)
        if os.path.exists(file_path): 
            logger.debug(f'rm {file_path}')
            os.remove(file_path)                       

        logger.debug('Download stability change: SUCCESS')
        logger.debug(f"Files downloaded in directory {path}")
        
        return file_path.replace(".zip", "")

    except Exception as e:
        logger.error('Download stability change: FAIL')
        logger.error(f"Error while downloading stability change: {e}")
        raise e


def append_ddg_to_dict(ddg_dict, df, frag=False):

    pattern = re.compile(r'([A-Za-z])(\d+)([A-Za-z])')
    
    for _, row in df.iterrows():
        variant, ddg = row.values
        pos, alt = extract_mut(variant, pattern)
        
        if pos not in ddg_dict:
            ddg_dict[pos] = {}
        
        if alt not in ddg_dict[pos] and frag:
            ddg_dict[pos][alt] = []

        if frag:
            ddg_dict[pos][alt].append(ddg)
        else:
            ddg_dict[pos][alt] = ddg

    return ddg_dict


def extract_mut(variant_str, pattern):

    match = pattern.match(variant_str)
    pos = match.group(2)
    alt = match.group(3)

    return pos, alt


def save_json(path_dir, uni_id, dictionary):
    
    with open(os.path.join(path_dir, f"{uni_id}_ddg.json"), "w") as json_file:
        json.dump(dictionary, json_file)


def id_from_ddg_path(path):
    
    return os.path.basename(path).split('-')[1]


def parse_ddg_rasp_worker(args):
    
    file, path_dir, output_path = args
    
    # Get Uniprot_ID
    uni_id = id_from_ddg_path(file)

    # Get paths of all fragments 
    lst_path_prot = glob.glob(os.path.join(path_dir, f"*{uni_id}*"))
    frag = True if len(lst_path_prot) > 1 else False

    # Save a dictionary for each pos with keys as ALT and lst of DDG as values
    ddg_dict = {}
    for path_prot in progressbar(lst_path_prot):
        df = pd.read_csv(path_prot)[["variant", "score_ml"]]
        ddg_dict = append_ddg_to_dict(ddg_dict, df, frag=frag)

    # Iterate through the pos and the ALT and get the mean across frags for each variant
    if frag:
        for pos in ddg_dict:
            for alt in ddg_dict[pos]:
                ddg_dict[pos][alt] = np.mean(ddg_dict[pos][alt])    

    # Save dict
    save_json(output_path, uni_id, ddg_dict)


def parse_ddg_rasp(input_path, output_path, threads=1):
    """
    It iterates through the csv files in <path_dir> and convert each one into 
    a .json dictionary of dictionaries having protein position as keys (str) and 
    ALT amino acid (1-letter) as sub-dictionaries keys whose values are the DDG
    (protein stability change upon mutations) for each variant predicted by RaSP.
    If a the protein is fragmented, the DDG of a variant is computed as average 
    DDG of that variant across the different fragments (fragments are overlapping). 

    Rapid protein stability prediction using deep learning representations
    https://elifesciences.org/articles/82593
    DOI: 10.7554/eLife.82593
    """

    # Get already processed files and available ones for processing
    files_processed = glob.glob(os.path.join(output_path, "*.json"))
    lst_files = [file for file in os.listdir(input_path)
                 if file.endswith(".csv") and os.path.join(output_path, f"{id_from_ddg_path(file)}.json") not in files_processed]
    ## Save dict for each proteins
    logger.debug(f"Input: {input_path}")
    logger.debug(f"Output: {output_path}")
    if len(lst_files) > 0:
        logger.debug(f"Parsing DDG of {len(lst_files)} proteins...")
        
        # TODO: for now it is created a process for each protein, while it would
        #       be better to have chunks of protein processed by the same process
        #       to decrese runtime (at the moment quite slow, 1h40m with 40 cores)
        
        # TODO: also the parsing itself can be optimized
        
        # Create a pool of workers parsing processes
        with Pool(processes=threads) as pool:
            args_list = [(file, input_path, output_path) for file in lst_files]
            # Map the worker function to the arguments list
            pool.map(parse_ddg_rasp_worker, args_list) 
        if len(lst_files) > 50:
            os.system('clear')
            logger.debug(f"clear")
        logger.debug(f"DDG succesfully converted into json files...")
    else:
        logger.debug(f"DDG not found: Skipping...")
        
    # Remove the original folder
    logger.debug(f"Deleting {input_path}")
    rm_dir(input_path)
    logger.info(f"Parsing of DDG completed!")