import logging
import os

import daiquiri
import subprocess
import pandas as pd
import numpy as np
import re
import glob
from tqdm import tqdm
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
        #clean_log = ["rm", "-rf", os.path.join(path, "log")]

        logger.debug(clean_files)
        subprocess.run(clean_files, shell=True)

        logger.debug(' '.join(clean_ddg))
        subprocess.run(clean_ddg)

        # logger.debug(' '.join(clean_log))
        # subprocess.run(clean_log)

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
            args_list = [(file, input_path, output_path) for file in lst_files]
            # Map the worker function to the arguments list
            pool.map(parse_ddg_rasp_worker, args_list) 
        if len(lst_files) > 50:
            subprocess.run("clear")
            logger.debug(f"clear")
        logger.debug(f"DDG succesfully converted into json files...")
    else:
        logger.debug(f"DDG not found: Skipping...")
        
    # Remove the original folder
    rm_ddg_dir = ["rm", "-rf", input_path]
    logger.debug(' '.join(rm_ddg_dir))
    subprocess.run(rm_ddg_dir)
    logger.info(f"Parsing of DDG completed!")


# #=================
# # pdb_tool
# #=================

def run_pdb_tool(pdb_tool_sif_path, i, o, F="4"):
    """
    Use PDB_Tool to extract features from all pdb files in directory.
    """
    pdb_files = os.listdir(i)
    
    for file in tqdm(pdb_files):
        if re.search('.\.pdb$', file) is not None:
            run = subprocess.run(["singularity", "exec", f"{pdb_tool_sif_path}", "/PDB_Tool/PDB_Tool", "-i", f"{i}/{file}", "-o", f"{o}{file}"])
            
            
def load_pdb_tool_file(path):
    """
    Parse .feature file obtained by PDB_Tool.
    """
    with open(path, "r") as f:
        lines = f.readlines()
        lst_lines = []
        for l in lines:
            lst_lines.append(l.strip().split())
   # return lst_lines
    df = pd.DataFrame(lst_lines[1:], columns = lst_lines[0]).iloc[:,:9]
    df = df.drop("Missing", axis=1)
    df = df.rename(columns={"Num" : "Pos"})

    for c in df.columns:
        if c == "Pos" or c == "ACC" or c == "CNa" or c == "CNb":
            data_type = int
        else:
            data_type = float
        try:
            df[c] = df[c].astype(data_type)
        except:
            pass
    return df


def get_pdb_tool_file_in_dir(path):
    """
    Get the list of PDB_Tool .feature files from directory.
    """
    list_files = os.listdir(path)
    ix = [re.search('.\.feature$', x) is not None for x in list_files]
    return list(np.array(list_files)[np.array(ix)])


def load_all_pdb_tool_files(path):
    """
    Get the features of all proteins in the directory.
    """
    df_list = []
    
    feature_file_list = get_pdb_tool_file_in_dir(path)
    
    for file in tqdm(feature_file_list):
        df = load_pdb_tool_file(f"{path}/{file}")
        identifier = file.split("-")
        print(identifier)
        df["Uniprot_ID"] = identifier[1]
        df["F"] = identifier[2].replace("F", "")
        df_list.append(df)

    return pd.concat(df_list).reset_index(drop=True)


# # Conversion for pdb_tool
# df = pd.read_csv("/workspace/projects/alphafold_features/feature_extraction/pdb_tool/pdb_tool_result_3s_sse_all_prot.csv")
# cols = ['Coil', 'Helix', 'Ladder']
# conditions = [df[cols[0]] == 1,
#               df[cols[1]] == 1,
#               df[cols[2]] == 1]
# df['Structure'] = np.select(conditions, cols, default='OtherLabel')
# df = df.drop(columns=cols)
# df.to_csv("/workspace/projects/clustering_3d/o3d_analysys/datasets/annotation/pdb_tool_annotation.csv", index=False)


#==============
#     MAIN 
#==============     

def get_annotations(path_pdb_structure,
                    path_pdb_tool_sif,
                    output_dir,
                    cores,
                    verbose):
    """
    Main function to build annotations to generate annotated plots.
    """
    
    logger.info(f"Path to PDB structures: {path_pdb_structure}")
    logger.info(f"Output to PDB_Tool SIF: {path_pdb_tool_sif}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Cores: {cores}")
    logger.info(f"Verbose: {bool(verbose)}")
    logger.info(f'Log path: {os.path.join(output_dir, "log")}')
    logger.info("")

    # Empty directory
    clean_annot_dir(output_dir, 'd')

    # Download DDG
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
    
    # TODO: Enable multiprocessing for PDB_Tool 
    
    # Run PDB_Tool
    logger.info(f"Extracting pACC and 2° structure...")
    pdb_tool_output = f"{output_dir}/pdb_tool"
    #path_pdb_structure = /workspace/projects/clustering_3d/clustering_3d/datasets_backup/datasets_distances/pdb_structures
    #pdb_tool_sif_path = /workspace/projects/clustering_3d/clustering_3d/build/containers/pdb_tool.sif
    run_pdb_tool(path_pdb_tool_sif, i=path_pdb_structure, o=pdb_tool_output)
    logger.info(f"pACC and secondary 2° extraction completed!")
    
    # Parse PDB_Tool
    logger.info(f"Parsing pACC and 2° structures...")
    pdb_tool_df = load_all_pdb_tool_files(pdb_tool_output)
    pdb_tool_df.to_csv(f"{output_dir}/pdb_tool_df.csv", index=False)
    logger.debug(f"Deleting {pdb_tool_output}")
    os.rmdir(pdb_tool_output)
    
    # Convert secondary structure to 3 feat only
    
    logger.info(f"pACC and 2° structures parcing completed!")
    
    






##### WIP PDB_Tool MP #######


#     # Run PDB_Tool
#     logger.info(f"Extracting pACC and 2° structure...")
#     pdb_tool_output = f"{output_dir}/pdb_tool"
#     #path_pdb_structure = /workspace/projects/clustering_3d/clustering_3d/datasets_backup/datasets_distances/pdb_structures
#     #pdb_tool_sif_path = /workspace/projects/clustering_3d/clustering_3d/build/containers/pdb_tool.sif
#     run_pdb_tool_mp(input_directory=path_pdb_structure, 
#                     output_directory=pdb_tool_output,
#                     pdb_tool_sif_path=path_pdb_tool_sif,
#                     num_processes=cores)
#     logger.info(f"pACC and secondary 2° extraction completed!")
    
#     # Parse PDB_Tool
#     logger.info(f"Parsing pACC and 2° structures...")
#     pdb_tool_df = load_all_pdb_tool_files(pdb_tool_output)
#     pdb_tool_df.to_csv(f"{output_dir}/pdb_tool_df.csv", index=False)
#     os.rmdir(pdb_tool_output)
#     logger.info(f"pACC and 2° structures parcing completed!")
    

# def run_pdb_tool(lst_path_pdb_files, pdb_tool_sif_path, i, o, F="4"):
#     """
#     Use PDB_Tool to extract features from all pdb files in directory.
#     """
    
#     for file in lst_path_pdb_files:
#         run = subprocess.run(["singularity", "exec", f"{pdb_tool_sif_path}", "/PDB_Tool/PDB_Tool", "-i", f"{i}/{file}", "-o", f"{o}{file}"])
    
    
# def run_pdb_tool_mp(input_directory, output_directory, pdb_tool_sif_path, F="4", num_processes=4):
#     pdb_files = [path for path in os.listdir(input_directory) if path.endswith(".pdb")]

#     # Divide pdb_files into chunks based on num_processes
#     chunks = [pdb_files[i:i + num_processes] for i in range(0, len(pdb_files), num_processes)]

#     with Pool(processes=num_processes) as pool:
#         pool.starmap(run_pdb_tool, [(chunk, pdb_tool_sif_path, input_directory, output_directory, F) for chunk in tqdm(chunks, desc="Processing chunks")])
        
            