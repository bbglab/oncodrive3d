import logging
import os

import daiquiri
import subprocess
import pandas as pd
import numpy as np
import re
import glob
#from tqdm import tqdm
from progressbar import progressbar
from multiprocessing import Pool
import json
import click

from scripts import __logger_name__
from scripts.globals import setup_logging_decorator
from scripts.datasets.utils import download_single_file, extract_zip_file

logger = daiquiri.getLogger(__logger_name__ + ".build.annotation")

logging.getLogger('urllib3.connectionpool').setLevel(logging.WARNING)


#=================
# Utils
#=================


def clean_annotations_dir(path: str, loc: str) -> None:
    """
    Clean the annotations directory by removing specific files 
    and subdirectories.

    Args:
        path (str): Path to the directory to be cleaned.
    """

    if loc == "d":

        clean_files = f"rm -rf {os.path.join(path, '*.csv')} {os.path.join(path, '*.json')} {os.path.join(path, '.*.txt')}"
        clean_ddg = ["rm", "-rf", os.path.join(path, "stability_change")]
        clean_log = ["rm", "-rf", os.path.join(path, "log")]

        logger.debug(clean_files)
        subprocess.run(clean_files, shell=True)

        logger.debug(' '.join(clean_ddg))
        subprocess.run(clean_ddg)

        logger.debug(' '.join(clean_log))
        subprocess.run(clean_log)

    elif loc == "r":
        # TODO: implement cleaning function for output
        pass


def clean_annot_dir(path: str, loc: str = 'd') -> None:
    """
    Clean it upon request if it already exists.

    Args:
        path (str): Path to the directory to be created or cleaned.
    """

    if os.listdir(path) != ['log']:
        logger.warning(f"Directory {path} already exists and is not empty.")

        overwrite = "y" if click.get_current_context().params['yes'] else input("Clean existing directory? (y/n): ")
        while overwrite.lower() not in ["y", "yes", "n", "no"]:
            print("Please choose yes or no")
            overwrite = input("Clean existing directory? (y/n): ")

        if overwrite.lower() in ["y", "yes"]:
            clean_annotations_dir(path, loc)
            logger.info(f"Dataset files in {path} have been removed.")
        else:
            logger.warning(f"Dataset files in {path} have not been removed.")
    else:
        pass


#=================
# Stability change
#=================


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
        ## STEP1 --- Download file
        logger.debug(f'Downloading to {file_path}')
        download_single_file(download_url, file_path, threads)
        
        ## STEP2 --- Extract from zip
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
        logger.error(f"Stability change will not be used for annotation but it will not affect the 3D-clustering analysis.")
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
    
    with open(f"{path_dir}/{uni_id}.json", "w") as json_file:
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

    # pattern_af = re.compile(r'-F(\d+)-')

    # Get already processed files and available ones for processing
    files_processed = glob.glob(os.path.join(output_path, "*.json"))
    lst_files = [file for file in os.listdir(input_path)
                 if file.endswith(".csv") and f"{output_path}/{id_from_ddg_path(file)}.json" not in files_processed]
    
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
            args_list = [(file, input_path, output_path, True) for file in fragmented_pdb]
            # Map the worker function to the arguments list
            pool.map(parse_ddg_rasp_worker, args_list) 
        if len(lst_files) > 50:
            subprocess.run("clear")
            logger.debug(f"clear")
        logger.debug(f"DDG succesfully converted into json files...")
    else:
        logger.debug(f"DDG not found: Skipping...")

        
        
    # fragmented_pdb = [
    #     file for file in os.listdir(input_path)
    #     if pattern_af.search(file) and int(pattern_af.search(file).group(1)) == 2 
    #     and f"{output_path}/{id_from_ddg_path(file)}.json" not in uni_id_processed
    # ]

    # ## Save dict for fragmented proteins
    # logger.debug(f"Input: {input_path}")
    # logger.debug(f"Output: {output_path}")
    # if len(fragmented_pdb) > 0:
    #     logger.debug(f"Parsing DDG of {len(fragmented_pdb)} fragments...")
        
    #     # TODO: for now it is created a process for each protein, while it would
    #     #       be better to have chunks of protein processed by the same process
    #     #       to decrese runtime (at the moment quite slow, 1h40m with 40 cores)
        
    #     # TODO: also the parsing itself can be optimized
        
    #     # Create a pool of workers parsing processes
    #     with Pool(processes=threads) as pool:
    #         args_list = [(file, input_path, output_path, True) for file in fragmented_pdb]
    #         # Map the worker function to the arguments list
    #         pool.map(parse_ddg_rasp_worker, args_list) 
    #     if len(fragmented_pdb) > 50:
    #         subprocess.run("clear")
    #         logger.debug(f"clear")
    #     logger.debug(f"Parsing of DDG fragments completed...")
    # else:
    #     logger.debug(f"DDG fragments not found: Skipping...")
    # uni_id_frag = glob.glob(os.path.join(output_path, "*.json"))

    # ## Save dict for non-fragmented proteins
    # uni_id_files_not_frag = [file for file in [f for f in os.listdir(input_path) if f.endswith(".csv")] 
    #                         if os.path.join(output_path, f"{id_from_ddg_path(file)}.json") not in uni_id_frag]
    # logger.debug(f"Parsing DDG of {len(uni_id_files_not_frag)} non-fragmented proteins...")
    
    # # Create a pool of worker processes for non-fragmented proteins
    # with Pool() as pool:
    #     args_list = [(file, input_path, output_path, False) for file in uni_id_files_not_frag]
    #     pool.map(parse_ddg_rasp_worker, args_list)
        
    # if len(uni_id_files_not_frag) > 50:
    #     subprocess.run("clear")
    #     logger.debug(f"clear")
        
    # Remove the original folder
    rm_ddg_dir = ["rm", "-rf", input_path]
    logger.debug(' '.join(rm_ddg_dir))
    subprocess.run(rm_ddg_dir)
    logger.info(f"Parsing of DDG completed!")


# #=================
# # pdb_tool
# #=================


# # Conversion for pdb_tool
# df = pd.read_csv("/workspace/projects/alphafold_features/feature_extraction/pdb_tool/pdb_tool_result_3s_sse_all_prot.csv")
# cols = ['Coil', 'Helix', 'Ladder']
# conditions = [df[cols[0]] == 1,
#               df[cols[1]] == 1,
#               df[cols[2]] == 1]
# df['Structure'] = np.select(conditions, cols, default='OtherLabel')
# df = df.drop(columns=cols)
# df.to_csv("/workspace/projects/clustering_3d/o3d_analysys/datasets/annotation/pdb_tool_annotation.csv", index=False)




# @click.command(context_settings={"help_option_names" : ['-h', '--help']})
# @click.option("-o", "--output_dir", help="Path to dir where to store annotations", type=str)
# @click.option("-c", "--cores", type=click.IntRange(min=1, max=len(os.sched_getaffinity(0)), clamp=False), default=len(os.sched_getaffinity(0)),
#               help="Number of cores to use in the computation")
# @click.option("-v", "--verbose", help="Verbose", is_flag=True)
# @setup_logging_decorator
def get_annotations(output_dir,
                    cores,
                    verbose):
    """
    Main function to build annotations to generate annotated plots.
    """
    
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Cores: {cores}")
    logger.info(f"Verbose: {bool(verbose)}")
    logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    logger.info("")

    # Empty directory
    clean_annot_dir(output_dir, 'd')

    # DDG
    logger.info(f"Downloading stability change...")
    path_dir = f"{output_dir}/stability_change"
    if not os.path.isdir(path_dir):
        os.makedirs(path_dir)
        logger.debug(f'mkdir {path_dir}')
    file_path = download_stability_change(path_dir, cores)
    logger.info(f"Download stability change completed!")
    
    # Parsing DDG
    logger.info(f"Parsing stability change...")
    parse_ddg_rasp(file_path, path_dir, cores)
    logger.info(f"Parsing stability change completed!")
    
    # pdb_tool
    

# if __name__ == "__main__":
#     from utils import download_single_file, extract_zip_file
#     build_annotations()
    
# else:
#     from scripts.datasets.utils import download_single_file, extract_zip_file